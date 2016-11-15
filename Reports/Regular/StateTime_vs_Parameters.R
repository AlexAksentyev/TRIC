# rm(list= ls())
# cat("\014")

source("Parameters.R")
source("EstiStats.R")
source("Utilities.R")

StateTimeWrap <- function(.rslsd, QIPrec){
  IntPeriod = Tint
  EstCoef = C
  BeamTime = H*Duty
  
  return(StateTime(IntPeriod, EstCoef, .rslsd, QIPrec, BeamTime))
}

BeamTimeWrap <- function(.rslsd, QIPrec){
  IntPeriod = Tint
  EstCoef = C
  ResolStD = .rslsd
  h = 15*60
  
  return(BeamTime(.rslsd, IntPeriod, EstCoef, h, QIPrec))
}

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

rsd = c(1%o% 10^(-4:-6)); names(rsd) <- as.character(rsd)

Precision <- c(1:10 %o% 10^(-6:-4))

QIprec.trg = 1e-6; trg.i = which(QIprec.trg %in% Precision)

acK = 1 #HAC correction

library(plyr)
ldply(Precision, function(p) BeamTimeWrap(rsd, p)) %>% melt(variable.name="ResSD", value.name="H") %>% 
  cbind(data.frame("Prec" = rep(Precision, times=length(rsd)))) -> Hplot; head(hplot)

ggplot(Hplot, aes(Prec, H/86400, col=ResSD)) + geom_line() +
  scale_x_log10(name=expression(sigma[EA[yy]]~"(a.u.)"), minor_breaks=Precision, labels=fancy_scientific) + 
  scale_y_log10(name="H (days)", minor_breaks=(c(1:10) %o% 10^(-4:4))) + 
  guides(col=guide_legend(title="Resolution")) + 
  theme_minimal() +
  theme(panel.grid.minor = element_line(colour = "gray", linetype = "dashed"), legend.position="top")

## variable slp
vmAe <- function(var.meas,var.slope,h){
  return(
    4/H * C^2 * (12*Tint*var.meas/h^2 + h*var.slope)
  )
}

vb = 1e-15
sdd = c(5e-3, 1e-3, 5e-4, 1e-4, 5e-5); vd = sdd^2; names(vd) <- as.character(sdd)
hbest = (24*Tint* vd/vb)^(1/3)
h = seq(250, 3600, by = 1)

# par(xpd=TRUE, mar=c(5, 5, 1, 3))
# plot(
#   c(0, hmax), c(1e-3, 1e0), type="n", xaxt="n", bty="l",
#   xlab = "state time, seconds",
#   ylab = "mean QoI estimate sd",
#   log = "y"
# )
# 
# for(p in 1:length(vd)) lines(h, sqrt(vmAe(vd[p],vb,h)), col=rainbow(length(vd))[p])

ldply(h, function(h) sqrt(vmAe(vd, vb,h))) %>% melt(variable.name="Resolution",value.name="SE") %>%
  cbind(data.frame("h" = rep(h, times=length(vd)))) -> SEAM; head(SEAM)

ggplot(SEAM, aes(h, SE, col=Resolution)) + geom_line() +
  scale_y_log10(name=expression(sigma[EA[yy]]~"(a.u.)")) +
  scale_x_continuous(name="h (sec)", minor_breaks = seq(200,3600, by=100)) +
  theme_minimal() +
  theme(legend.position="top",panel.grid.minor = element_line(colour = "gray", linetype = "dashed"))

# 
# legend(
#   "topright", inset=c(0,0), title=expression("sd"~delta ~~ "&" ~~ h[min]~"sec"),
#   lty=1, bty="n", col = rainbow(length(vd)), ncol = 2,
#   legend=c(formatC(sdd,1,format="e"),formatC(hmin,0,format="f"))
# )
# legend(
#   "bottomleft", bty="n", legend = bquote(var~beta ~ "=" ~ .(vb))
# )
# 
# axis(1, at = formatC(c(hmin[2:5],hmax),0,format = "f"))


