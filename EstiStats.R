source("DataPrep.R")
source("DataAna.R")

#### experiment-related estimates ####
ResolStD <- function(ErrorStD, IniBeamCurrent) return(ErrorStD/IniBeamCurrent)

# The slope estimate standard deviation (theoretically)
SlopeStD <- function(ResolStD, IntPeriod, StateTime){
  K = StateTime/IntPeriod
  fcr = (2*sqrt(3))/(IntPeriod * sqrt(K*(K^2-1)))
  return(ResolStD*fcr)
}

# inherent slope variation's sd
bStD <- function(p, N, Dt){return(sqrt(p*(1-p)/N)/Dt)}

# The correlation coefficient estimate's sampling distribution's standard deviation
QIMeanEstStD <- function(ResolStD, IntPeriod, StateTime, TotBeamTime, EstCoef){
  C = 2*EstCoef
  return(2*sqrt(3)*C*sqrt(IntPeriod)/(StateTime*sqrt(TotBeamTime))*ResolStD)
}

# Required beam time
BeamTime <- function(ResolStD, IntPeriod, EstCoef, StateTime, QIPrec){
  frac = 4*sqrt(3*IntPeriod)/(StateTime*QIPrec)
  Hmin = (frac*EstCoef*ResolStD)^2
  return(Hmin)
}

StateTime <- function(IntPeriod, EstCoef, ResolStD, QIPrec, BeamTime){
  frac = ResolStD/QIPrec
  return(4*sqrt(3*IntPeriod)*EstCoef*frac/sqrt(BeamTime))
}


###### commonly used statistics ####
library(Hmisc)
MDE = function(x){d = density(x$Estimate[!is.na(x$Estimate)]); d$x[which.max(d$y)]}
MDN = function(x) median(x$Estimate, na.rm = TRUE)
WMN = function(x) wtd.mean(x$Estimate, 1/x$SE^2, na.rm = TRUE)

L2Norm. <- function(X, standardize = TRUE){
  X <- X/sum(X); d = length(X)
  
  sqrt(sum(X^2))-> n
  
  if(standardize)
    (n * sqrt(d) - 1)/(sqrt(d) - 1)
  else
    n
}

.chi2 <- function(x) {m = mean(x); v = var(x); sum((x-m)^2)/v/(length(x)-1)}

.collect.stats <- function(m){
  require(sandwich)
  chi2 = with(m, sum(residuals^2)/var(residuals)/df.residual)
  c(
    Estimate = summary(m)$coefficients[2, 1], 
    SE       = sqrt(vcovHAC(m)[2,2]),
    Chi2     = chi2,
    I0       = exp(summary(m)$coefficients[1, 1]),
    Start    = m$UTS[1]
  )
}

.markOutliers <- function(.data, what = "Estimate", .sigmas = 3){
  require(mosaic)
  r = (.sigmas - qnorm(.75))/(qnorm(.75)-qnorm(.25))
  
  boxplot(.data[,what], plot = FALSE, range = r)$out -> out
  
  .data %>% mutate(
    FOut = derivedFactor(
      "F" = !.data[,what] %in% out,
      "T" = .data[,what] %in% out
    )
  )
}

##### slope statistic ####
library(dplyr); library(plyr); library(splitstackshape)
library(parallel); library(doParallel)
source("DataPrep.R")

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
      ) %>% mutate(Clock = as.POSIXct(Clock, origin="1970-1-1")) %>% .outliers -> slopes16
  
  return(list(Data = Data16, Slopes = slopes16))
}

get2012Data <- function(){
  require(forecast); require(changepoint)

  ## preparing metadata
  where = "~/Analysis/Sep12/BCTAnalyser/results/"
  readODS::read_ods(paste(where, "RunMetadata(1).ods",sep="/"), sheet=2) %>% 
    filter(Run %in% c(935:937, 942:977)) -> MD
  for(obs in c("eCool.I", "Bunch", "Cell", "Targ", "RF.V", "BB.V", "Beam")) ## could use sapply(, factor) instead of loop
    MD[,obs] <- factor(MD[,obs])
  
  MD %>% mutate(
    Start = paste(day, Start, sep = " "),
    Stop = paste(day, Stop, sep = " "),
    UTS.Stt = timeToUTS(Start, format = "%Y-%m-%d %H:%M:%S"),
    UTS.Stp = timeToUTS(Stop, format = "%Y-%m-%d %H:%M:%S"),
    Size = UTS.Stp-UTS.Stt+1,
    Quality = derivedFactor("bad" = Quality=="bad", "ok" = is.na(Quality))
  ) %>% dplyr::select(-day,-Beam, -Comment) -> MD
  
  .from = MD$Start %>% first; .to = MD$Stop %>% last
  
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
  
  # ## finding correct cycle offsets
  # Data12 = Data12 %>% ddply("Size", function(.sub){
  #     .sub %>% filter(FSgl=="F") %>% .markOutliers(what="BCT2") %>% filter(FOut=="F") -> .subb
  #     cpt.mean(.subb$BCT2)@cpts -> brk
  #     
  #     .subb %>% mutate(Group = rep(.subb$Run[brk], diff(c(0,brk)))) -> .subb
  #     mosaic::median(BCT2~Group, data=.subb)->bct2o
  #     bct2o%>%rep(c(length(which(.sub$Run <= names(bct2o)[1])), length(which(.sub$Run > names(bct2o)[1])))) -> bct2o
  #     
  #     .sub %>% mutate(BCT2o = bct2o)
  #   })
  
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

extractARMA <- function(u, out.pct = c(15,5)){
  n = nrow(u)
  if(n >= 270){
    lv <- ceiling(out.pct/100*c(1,-1)*n) + c(1, n)
    us <- slice(u, lv[1]:lv[2])
    us$BCT2 %>% Arima(c(4,0,4), include.drift=TRUE) -> .arima
    .arima$coef -> phima; .arima %>% vcov %>% diag %>% sqrt -> ses; abs(phima/ses) ->tscore
    qt(tscore, n-1, lower.tail=FALSE)->pvals
    
    plot((us$BCT2), type="l", main = paste(u$Run[1], u$Unit[1], sep="-")); lines(fitted(.arima), col="red")
    
    cbind(phima,ses,tscore,pvals) %>% as.data.frame() #%>% cbind(data.frame(CName=names(phima))) %>% print
  }
  
}





#### testing ####
f <- function(u){
  with(u, plot(BCT2~UTS, col=FSgl, main = u$Unit[1]))
  
  ub = filter(u, FSgl == "F")
  us = filter(u, FSgl == "T", eCool == "Acc", !is.na(FABS)) %>% mutate(BCT2 = BCT2 - ub$BCT2 %>% median)
  
  us %>% ddply("FABS", function(.sub) {fit.(.sub,plot.it=FALSE) -> .stats 
                 data.frame(
                   T.Spin = .sub$T.Spin[1],
                   B.Spin = .sub$B.Spin[1],
                   Clock = .sub$UTS[1] %>% as.POSIXct(origin="1970-1-1"),
                   Estimate = .stats[1],
                   SE = .stats[2],
                   Chi2 = .stats[3],
                   I0 = .stats[4]
                 )})
  
}
