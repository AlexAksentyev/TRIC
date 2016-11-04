
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
MDE = function(x){d = density(x$Estimate[!is.na(x$Estimate)]); d$x[which.max(d$y)]}
MDN = function(x) median(x$Estimate, na.rm = TRUE)
WMN = function(x) weighted.mean(x$Estimate, 1/x$SE^2, na.rm = TRUE)

L2Norm. <- function(X, standardize = TRUE){
  X <- X/sum(X); d = length(X)
  
  sqrt(sum(X^2))-> n
  
  if(standardize)
    (n * sqrt(d) - 1)/(sqrt(d) - 1)
  else
    n
}

.chi2 <- function(x) {m = mean(x); v = var(x); sum((x-m)^2)/v/(length(x)-1)}
.se <- function(x) {sd(x$Estimate)/sqrt(nrow(x))}

.collect.stats <- function(m){
  require(sandwich)
  chi2 = with(m, sum(residuals^2)/var(residuals)/df.residual)
  c(
    Estimate = summary(m)$coefficients[2, 1], 
    SE       = sqrt(vcovHAC(m)[2,2]),
    Chi2     = chi2,
    RSE      = summary(m)$sigma,
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
pair.up <- function(.slp, .fields = NULL){
  .f = c("Unit", "Estimate", "SE")
  if(is.null(.fields)) .fields = .f else c(.f, .fields) %>% unique -> .fields
  llply(
    .fields, function(f){
      expand.grid(.slp[[1]][,f], .slp[[2]][,f]) -> x
      names(x) <- paste(f, names(.slp), sep=".")
      return(x)
    }
  ) -> fl
  
  do.call(cbind, fl)
}

dbeta. <- function(.slp, .fields = NULL){
  # .f = c("Unit", "Estimate", "SE")
  # if(is.null(.fields)) .fields = .f else c(.f, .fields) %>% unique -> .fields
  # llply(
  #   .fields, function(f){
  #     expand.grid(.slp[[1]][,f], .slp[[2]][,f]) -> x
  #     names(x) <- paste(f, names(.slp), sep=".")
  #     return(x)
  #   }
  # ) -> fl
  
  # do.call(cbind, fl) -> .data
  pair.up(.slp, .fields) -> .data
  .data %>% mutate(
    Estimate = .data[,3] - .data[,4],
    SE = sqrt(.data[,5]^2 + .data[,6]^2) #ideally this should take account of correlation
  )
}

rbeta. <- function(.slp, .fields = NULL){
  pair.up(.slp, .fields) -> .data
  
  .data %>% mutate(
    Estimate = .data[,3]/.data[,4],
    SE = Estimate * sqrt((.data[,5]/.data[,3])^2 + (.data[,6]/.data[,4])^2)
  ) -> B
  
  B %>% mutate(
    SE = SE * 2/(1+Estimate)^2,
    Estimate = (Estimate - 1)/(Estimate + 1)
  ) -> R
  
  R
}
