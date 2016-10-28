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
      function(.sub) {fit.(.sub,plot.it=FALSE) %>% c(Clock = .sub$UTS[1])}, .parallel = TRUE
    ) %>% mutate(Clock = as.POSIXct(Clock, origin="1970-1-1")) %>% .markOutliers -> slopes16
  
  return(list(Data = Data16, Slopes = slopes16))
}

get2012Data <- function(){ #not removing offset
  require(forecast); require(splitstackshape)
  
  ## preparing metadata
  where = "~/Analysis/Sep12/BCTAnalyser/results/"
  readODS::read_ods(paste(where, "RunMetadata(1).ods",sep="/"), sheet=2) %>% 
    filter(Run %in% c(935:937, 942:977)) -> MD
  for(obs in c("eCool.I", "Bunch", "Cell", "Targ", "RF.V", "BB.V", "Beam")) ## could use sapply(, factor) instead of loop
    MD[,obs] <- factor(MD[,obs])
  
  MD %>% mutate(
    Start = paste(day, Start, sep = " ")%>%as.POSIXct(origin="1970-1-1"),
    Stop = c(Start[-1]-1, paste(day[nrow(MD)], Stop[nrow(MD)], sep = " ")%>%as.POSIXct(origin="1970-1-1")),
    UTS.Stt = as.numeric(Start),
    UTS.Stp = as.numeric(Stop),
    Size = UTS.Stp-UTS.Stt+1,
    Quality = derivedFactor("bad" = Quality=="bad", "ok" = is.na(Quality))
  ) %>% dplyr::select(-day,-Beam, -Comment) -> MD
  
  .from = MD$Start %>% first; .to = MD$Stop %>% last
  
  ## a cycle starts where current goes to zero
  ## I can load the part of the file .from--.to, find the lines where BCT2 is zero
  ## and mark the times at those lines as Start
  
  MD %>% expandRows("Size") -> MD
  
  rownames(MD) %>% strsplit(".", fixed=TRUE) %>% 
    laply(function(e) {e <- as.numeric(e); if(length(e) <2) 0 else e[2]}) -> addage
  MD %>% mutate(UTS = UTS.Stt + addage) -> MD
  
  #### loading & modifying data 
  .flds = c("WD", "Month", "Day", "Time", "Year", paste0("BCT",1:8))
  loadData("CosyAll_2012-09-25_15-53-06.dat", "~/Analysis/Sep12/DevReader/Data/", Fields = .flds, from=.from, to=.to) -> Data12
  
  MD %>% dplyr::select(Run, UTS, eCool.I, Bunch, RF.V, BB.V, Cell, Targ, Quality) %>% 
    join(Data12, by="UTS", type="right") %>% dlply("Run") -> Data12
  
  names(Data12) %>% as.numeric -> .runs; .runs <- .runs[1:30]
  bri = .runs>=965 | .runs %in% 941:944; 
  sri = !.runs %in% .runs[bri]
  
  .dumbCut <- (function(.data, x) mutate(.data, FSgl = derivedFactor("F" = BCT2 < x, .default = "T")))
  
  list(
    Small = ldply(Data12[which(sri)], .id = "Run") %>% .dumbCut(2700), 
    Big = ldply(Data12[which(bri)], .id = "Run") %>% .dumbCut(5000)
  ) %>% ldply(.id = "Size") %>% ddply(.(Size,Run), .unitize) -> Data12
  
  Data12 %>% ddply("Size", function(u) {filter(u, FSgl=="F")$BCT2 %>% median -> b; mutate(u, BCT2 = BCT2-b)}) -> Data12
  
  ## fitting for slopes
  Data12 %>% filter(FSgl == "T", Cell == "In") %>%
    ddply(.(Run, Unit), function(.sub) {
      fit.(.sub, .subset = c(45, 15)) -> .stats
      
      data.frame(
        Estimate = .stats[1], SE = .stats[2], Chi2 = .stats[3], I0 = .stats[4],
        Quality = .sub$Quality[5],
        eCool.I = .sub$eCool.I[5],
        Targ = .sub$Targ[5],
        Bunch = .sub$Bunch[5],
        Size = .sub$Size[5],
        Clock = .sub$UTS[1] %>% as.POSIXct(origin="1970-1-1")
      )
    }, .parallel = TRUE) %>% .markOutliers -> slopes12
  
  return(list(Data = Data12, Slopes = slopes12))
}

fit. = function(.sub, plot.it = FALSE, model = FALSE, .legend = c("Chi2"), .subset = NULL, .obs = "Unit", ...){
  .sub %>% mutate(uts = UTS - UTS[1]) -> .sub
  
  if(!is.null(.subset)) {
    .subset = .subset*c(1,-1) + c(0, dim(.sub)[1]); .subset = .subset[1]:.subset[2]
  } else .subset = 1:dim(.sub)[1]
  
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
