rm(list = ls(all=TRUE))
source("EstiStats.R")
source("Parameters.R")
source("Utilities.R")

library(reshape2)

lblfnt = 18
thm = theme_bw() + 
  theme(axis.text=element_text(size=lblfnt), axis.title=element_text(size=lblfnt), 
        legend.title=element_text(size=lblfnt), legend.text=element_text(size=lblfnt), legend.position="top")

#### Beam Time ####
Prec = c(1:10%o%10^(-3:-6))
Res = 10^(-4:-6); names(Res) <- Res

ldply(Prec, function(p) c("Precision" = p, BeamTime(Res, Tint, C, 15*60, p)/3600/24)) %>% 
  melt(measure.vars=2:4, variable.name="Resolution", value.name="BeamTime") -> df

ggplot(df, aes(Precision, BeamTime, col=Resolution)) + geom_line() + thm +
  scale_y_log10(name="H (days)", minor_breaks=c(1:10%o%10^(-4:+4)), labels=fancy_scientific) + 
  scale_x_log10(name=expression(sigma[EA[yy]]~"(a.u.)"), minor_breaks=Prec, labels=fancy_scientific) +
  scale_color_discrete(labels=fancy_scientific)

#### State Time ####
vmbe <- function(var.meas,var.slope,h){
  return(
    2*C^2* 2/H * (12*Tint*var.meas/h^2 + h*var.slope)
  )
}

vb = 1e-15
sdd = c(5e-3, 1e-3, 5e-4, 1e-4, 5e-5)
vd = sdd^2; names(vd) <- sdd
hmin = (24*Tint* vd/vb)^(1/3); hmin/60

ldply(seq(100,3600,100), function(h) c("h" = h, sqrt(vmbe(vd, vb, h)))) %>% 
  melt(measure.vars=1+seq(length(sdd)), variable.name="Resolution",value.name="SDQI") -> df

ggplot(df, aes(h, SDQI, col=Resolution)) + geom_line() + thm +
  scale_y_log10(name=expression(sigma[EA[yy]]~"(a.u.)"), breaks=10^(-4:-1), labels=fancy_scientific) +
  scale_x_continuous(name="h (sec)", minor_breaks=seq(100,3600,100)) +
  scale_color_discrete(labels=fancy_scientific)
#### slope variation from rho ####
Rho.ini = Rho.on
rate = 1.157e-6

ldply(seq(3600), function(n){
  Rho <- c(Rho.ini , Rho.ini * cumprod(1 - rnorm(n-1, mean = 1*rate, sd = 1*rate)))
  lam = Sigma0*d*nu*Rho
  c("h" = n, "vb" = var(lam))
}) -> df

ggplot(df, aes(h,vb)) + geom_line() + thm + labs(x="h (sec)", y=expression(sigma[beta]^2~(sec^-2))) +
  scale_y_log10(limits=c(1e-24, 1e-16), labels=fancy_scientific)
