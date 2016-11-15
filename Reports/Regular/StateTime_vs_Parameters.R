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

ldply(h, function(h) sqrt(vmAe(vd, vb,h))) %>% melt(variable.name="Resolution",value.name="SE") %>%
  cbind(data.frame("h" = rep(h, times=length(vd)))) -> SEAM; head(SEAM)

ggplot(SEAM, aes(h, SE, col=Resolution)) + geom_line() +
  scale_y_log10(name=expression(sigma[EA[yy]]~"(a.u.)")) +
  scale_x_continuous(name="h (sec)", minor_breaks = seq(200,3600, by=100)) +
  theme_minimal() +
  theme(legend.position="top",panel.grid.minor = element_line(colour = "gray", linetype = "dashed"))

daply(SEAM, "Resolution", function(x) min(x$SE)) -> SEBE

data.frame("hbest" = round(hbest/60), "SEbest" = formatC(SEBE,0,format="e"))

## slope variance 

library(MASS)
h = 3600
v = array(NA,h)
for(n in 1:h){
  thick <- d*c(Rho.on , Rho.on * cumprod(1 - rnorm(n-1, mean = 0*rate, sd = 1*rate)))
  lam = Sigma0*nu*thick
  v[n] = var(lam)
}

v<-v[-1]; n = 2:h
v.fit = rlm(v~0+n)
plot(n,v,log="y",
     type="l", 
     xlab = "h (sec)", 
     ylab = expression(sigma[beta]^2~"(h)")
)
lines(v.fit$fitted.values,col="red",lwd = 3)
# abline(h=0, lty=3, col="magenta")
# legend(
#   "topleft", bty="n",
#   legend = c(paste("dissipation rate =", formatC(rate*100,3,format="e"), "*(0 +- 1) %/sec"),
#              paste("v(h) =", formatC(v.fit$coefficients[[1]],2,format="e"), "*h"))
# )

