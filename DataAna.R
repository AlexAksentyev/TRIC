source("DataPrep.R")
source("EstiStats.R")

library(Hmisc)
library(mosaic); library(plyr); library(dplyr)
library(parallel); library(doParallel)

registerDoParallel(detectCores())



get2016Data <- function(){
  list(
    "3" = cleanData(loadData(from = "2016-6-30 20:20:00", to = "2016-6-30 21:50:00")),
    "1" = cleanData(loadData(from = "2016-6-30 22:50:00", to = "2016-7-1 08:00:00"))
  ) %>% 
    ldply(.id = "T.Spin", .parallel = TRUE) -> Data16 #all data #%>% 
  # filter(FSgl == "T", eCool == "Acc", !is.na(FABS)) -> Data16
  
  base = Data16 %>% filter(FSgl=="F"); ADC2o = median(base$BCT2)
  
  Data16 %>% filter(FSgl == "T", eCool == "Acc", !is.na(FABS)) %>% mutate(BCT2 = BCT2 - ADC2o) %>%
    ddply(
      .(T.Spin, B.Spin, Unit, FABS), 
      function(.sub) {fit.(.sub,plot.it=FALSE) %>% c(Clock = .sub$UTS[1], I0m = .sub$BCT2[60:120]%>%median)}, .parallel = TRUE
    ) %>% mutate(Clock = as.POSIXct(Clock, origin="1970-1-1")) %>% .markOutliers -> slopes16
  
  return(list(Data = Data16, Slopes = slopes16))
}

get2012Data <- function(){
  require(splitstackshape)
  
  .runs = 967:976
  
  ## fetch metadata
  where = "~/Analysis/Sep12/BCTAnalyser/results/"
  readODS::read_ods(paste(where, "RunMetadata(1).ods",sep="/"), sheet=2) %>% 
    filter(Run %in% .runs) %>% mutate(
      Start = paste(day, Start, sep = " ")%>%as.POSIXct(origin="1970-1-1")-30,
      Stop = paste(day, Stop, sep = " ")%>%as.POSIXct(origin="1970-1-1")+30
    ) -> MD
  for(obs in c("eCool.I", "Bunch", "Cell", "Targ", "RF.V", "BB.V", "Beam")) ## could use sapply(, factor) instead of loop
    MD[,obs] <- factor(MD[,obs])
  
  ## load data
  .flds = c("WD", "Month", "Day", "Time", "Year", paste0("BCT",1:8))
  loadData("CosyAll_2012-09-25_15-53-06.dat", "~/Analysis/Sep12/DevReader/Data/", 
           Fields = .flds, from=first(MD$Start), to=last(MD$Stop)
  ) -> Data12
  
  ## find when cycles begin
  filter(Data12, BCT2==0) %>% mutate(H = hour(Clock)) %>% 
    group_by(H) %>% dplyr::summarize(Clock = first(Clock)) -> CycSrt
  
  ## correct metadata
  MD %>% mutate(
    Start = CycSrt$Clock, Stop = c(Start[-1]-1, last(Stop)),
    UTS.Stt = as.numeric(Start),
    UTS.Stp = as.numeric(Stop),
    Size = UTS.Stp-UTS.Stt+1,
    Quality = derivedFactor("bad" = Quality=="bad", "ok" = is.na(Quality))
  ) %>% dplyr::select(-day,-Beam, -Comment) -> MD
  
  MD %>%expandRows("Size") -> MD
  rownames(MD) %>% strsplit(".", fixed=TRUE) %>% 
    laply(function(e) {e <- as.numeric(e); if(length(e) <2) 0 else e[2]}) -> addage
  MD %>% mutate(UTS = UTS.Stt + addage) -> MD
  
  ## add metadata to cycles
  MD %>% dplyr::select(Run, UTS, eCool.I, Bunch, RF.V, BB.V, Cell, Targ, Quality) %>% 
    join(Data12, by="UTS", type="right") %>% filter(!is.na(Run)) -> Data12 ## offset before the first cycle is discarded
  
  .dumbCut <- function(.data, x) mutate(.data, FSgl = derivedFactor("F" = BCT2 < x, .default = "T"))
  
  Data12 %>% .dumbCut(5000) %>% filter(Cell=="In") -> Data12
  
  ## fitting for slopes
  Data12 %>% ddply("Run", function(r) {
    b = filter(r, FSgl=="F", BCT2 != 0)$BCT2 %>% median
    mutate(r, BCT2 = BCT2-b) %>% filter(FSgl=="T") -> .sub
    cat(r$Run[1]); cat("\n")
    
    fit.(.sub, .subset = c(45, 15)) -> .stats; cat(.stats);cat("\n")
    
    data.frame(
      Unit=1,
      Estimate = .stats[1], SE = .stats[2], Chi2 = .stats[3], 
      I0 = .sub$BCT2[60:120]%>%median,
      Quality = .sub$Quality[5],
      eCool.I = .sub$eCool.I[5],
      Targ = .sub$Targ[5],
      Bunch = .sub$Bunch[5],
      Clock = .sub$UTS[1] %>% as.POSIXct(origin="1970-1-1")
    )
  }, .parallel=TRUE) %>% ddply("Targ", .markOutliers) -> slopes12
  
  return(list(Data = Data12, Slopes = slopes12))
}

fit. = function(.sub, plot.it = FALSE, model = FALSE, .legend = c("Chi2"), .subset = NULL, .obs = "Unit", ...){
  .sub %>% mutate(uts = UTS - UTS[1]) -> .sub
  
  if(!is.null(.subset)) {
    .subset = .subset*c(1,-1) + c(0, nrow(.sub)); .subset = .subset[1]:.subset[2]
  } else .subset = 1:nrow(.sub)
  
  m <- lm(log(BCT2) ~ uts, data = .sub, subset = .subset)
  m$model$UTS <- .sub$UTS[.subset]
  
  .stats = .collect.stats(m)

  if(plot.it){
    with(.sub, plot(log(BCT2), xaxt="n", xlab="time", main = paste("#", .sub[,.obs][1], "; Targ. ", .sub$FABS[1]), type="l"));
    ind = seq(1, dim(.sub)[1], length.out = 6); axis(1, at=ind, labels = .sub$Clock[ind])
    abline(m, col="red")
    legend(
      "topright", ncol=2,
      legend = c(
        .legend, 
        formatC(.stats[.legend], ...)
      )
    )
  }  
  
  if(model) return(m) else return(.stats)
  
}

scan4change <- function(.sub, plot.it = FALSE, size = .1, ...){
  require(strucchange); require(xts)
  .sub %>% mutate(uts = UTS - UTS[1]) -> .sub
  
  f = log(BCT2) ~ uts
  
  breakpoints(f, data=.sub, h=size) -> x
  bp = x$breakpoints; tb = .sub$UTS[bp]
  
  if(any(is.na(bp))) bp <- tb <- NULL
  .sub$Break <- rep(c(tb,last(.sub$UTS)), diff(c(0, bp, dim(.sub)[1])))
  
  dlply(.sub, "Break", fit., FALSE, TRUE, ...) -> .sub.mdls
  
  m = lm(f, data=.sub); chi2 = sum(m$residuals^2)/var(m$residuals)/dim(.sub)[1]
  
  if(plot.it){
    with(.sub, plot(log(BCT2) ~ uts, main = paste("#", .sub$Unit[1], "; Targ. ", .sub$FABS[1])))
    abline(m, col="red", lwd=2); abline(v = bp, col="gray", lty=2)
    
    ldply(.sub.mdls, function(x){
          with(x, lines(fitted.values~I(model$UTS-.sub$UTS[1]), col="green", lwd=2));
          return(.collect.stats(x))
      }
    ) ->.stats
    legend(
      "topright", lty=1, col=c("red", rep("green", length(.stats$Chi2))),
      legend = paste("chi2", round(c(chi2, .stats$Chi2), 3))
    )
  }  
  
  # if(sub.$B.Spin == "U")
  #   slp. <- slp.[which.min(slp.$Estimate),]
  # else if(sub.$B.Spin == "D")
  #   slp. <- slp.[which.max(slp.$Estimate),]
  # else
  #   slp. <- slp.[which.max(slp.$Estimate),]
  # 
  return(.stats) 
}

nl.sys2 <- function(x, estr. = WMN){
  y = numeric(2)
  
  y[1] = prod(x) - estr.(d1)
  y[2] = prod(x) - estr.(d2)
  
  return(y*1e27)
}

.QQNorm <- function(X, name ="", plot.it = TRUE){
  
  mu <- WMN(X); n <- dim(X)[1]; #se <- sqrt(sum(X$SE^2))/n
  se = sd(X$Estimate)/sqrt(n)
  chi2 <- with(X, sum(((Estimate - mu)/sd(Estimate))^2))/n
  
  
  if(plot.it){
    qqnorm(X$Estimate, main = paste("Normal Q-Q Plot for", name), xlab = "Normal Quantiles")
    qqline(X$Estimate, col="red"); grid()
    legend(
      "topleft", bty="n", ncol = 2,
      legend = c(
        c(expression(mu), expression(sigma[mu]), bquote(chi[.(n)]^2)),
        round(c(mu, se, chi2),3)
      )
    )
  }
  
  return(c(mean = mu, se = se, chi2 = chi2))
}
