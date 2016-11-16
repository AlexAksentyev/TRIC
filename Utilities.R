
# Plotting
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}

LogAxes <- function(X, Y){
  x <- floor(log10(range(X)))
  y <- floor(log10(range(Y)))
  powx <- seq(x[1], x[2]+1)
  powy <- seq(y[1], y[2]+1)
  xticksat <- as.vector(sapply(powx, function(p) (1:10)*10^p))
  yticksat <- as.vector(sapply(powy, function(p) (1:10)*10^p))
  axis(1, 10^powx); axis(2, 10^powy);
  axis(1, xticksat, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1)
  axis(2, yticksat, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1)
  abline(h = yticksat, v = xticksat, col = "gray", lty = 3)
}

present <- function(data.list){
  pres.mat = matrix(
    unlist(data.list),
    byrow = TRUE, nrow=length(data.list),
    dimnames = list(names(data.list), names(data.list[[1]]))
  )
  
  print(pres.mat)
}

PlotWBars <- function(Y, errY = rep(0,length(Y)), X = 1:length(Y), errX = rep(0,length(X)), ...){
  
  extra = list(...)
  Y <- c(Y); errY <- c(errY); X = c(X); errX <- c(errX)
  
  argnames = names(extra)
  if(!"col" %in% argnames) extra$col = "blue"
  if(!"ylab" %in% argnames) extra$ylab = ""
  if(!"xlab" %in% argnames) extra$xlab = ""
  
  maxerrX = max(errX); maxerrY = max(errY)
  
  xmin = min(X)-maxerrX; xmax = max(X)+maxerrX
  ymin = min(Y)-maxerrY; ymax = max(Y)+maxerrY
  
  plot(c(xmin, xmax), c(ymin, ymax), type="n", xlab = extra$xlab, ylab = extra$ylab, main = extra$main, ... = ...)
  
  points(
    Y ~ X, type="o",
    pch= '.', cex=3, lty = 0, col = extra$col, list(...)
  )
  if(max(errX) != 0){
    arrows(
      Y, X-errX,
      Y, X+errX,
      length=0.01, angle=90, code=3, col = extra$col
    )
  }
  if(max(errY) != 0){
    arrows(
      X, Y-errY,
      X, Y+errY,
      length=0.01, angle=90, code=3, col = extra$col
    )
  }
  
}

addGauss <- function(series, color = "darkblue"){
  x = seq(min(series), max(series), length.out = length(series))
  d = dnorm(x,mean(series),sd(series))
  lines(x,d,col=color)
}


# Time handling

deciTots <- function(deci, .format = "m:s"){
  require(chron)
  given = paste(floor(deci), round((deci-floor(deci))*60), sep=":")
  given <- switch(
    .format,
    "h:m" = paste(given, "00", sep=":"),
    paste("00", given, sep=":")
  )
  ts = times(given, format = "h:m:s")
  return(ts)
}

makeClock <- function(Time, summer.time = TRUE){
  require(chron)
  Clock = formatC(Time, width = 6, flag = "0")
  Clock <- chron::times(Clock, format = "hms", out.format = "h:m:s")
  
  if(summer.time) Clock <- Clock + chron::times("02:00:00", format = "h:m:s")
  
  return(Clock)
}

timeToUTS <- function(datetime, ...){
  return(as.numeric(as.POSIXct(datetime, origin="1970-01-01", ...))) 
}
