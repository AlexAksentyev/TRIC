library(Hmisc)
# library(nlme)
# library(lawstat)
# library(tseries)

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
