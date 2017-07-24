rm(list = ls(all=TRUE))
source("EstiStats.R")
source("Parameters.R")
source("Utilities.R")

library(reshape2)
library(ggplot2)
library(dplyr); library(plyr)

lblfnt = 16
thm = theme_bw() + 
  theme(axis.text=element_text(size=lblfnt), axis.title=element_text(size=lblfnt), 
        legend.title=element_text(size=lblfnt), legend.text=element_text(size=lblfnt), legend.position="top")

Prec = c(1:10%o%10^(-3:-6))
Res = 10^(-4:-6) 
#### Beam Time ####
h0 <- 15*60
names(Res) <- Res

ldply(Prec, function(p) c("Precision" = p, BeamTime(Res, Tint, C, h0, p)/3600/24)) %>% 
  melt(measure.vars=2:4, variable.name="Resolution", value.name="BeamTime") %>%
  mutate(StateTime = h0, Term="Stat") -> df0

ggplot(df0, aes(Precision, BeamTime, col=Resolution)) + geom_line() + thm +
  scale_y_log10(name="H (days)", minor_breaks=c(1:10%o%10^(-4:+4)), labels=fancy_scientific) +
  scale_x_log10(name=expression(sigma[EA[yy]]~"(a.u.)"), minor_breaks=Prec, labels=fancy_scientific) +
  scale_color_discrete(labels=fancy_scientific) -> p.BT.stat

#### State Time ####
vmbe <- function(var.meas,var.slope,h){
  return(
    2*C^2* 2/H * (12*Tint*var.meas/h^2 + h*var.slope)
  )
}

vb = 1e-15
vd = Res^2; names(vd) <- Res
hmin = (24*Tint* vd/vb)^(1/3); hmin/60
hmin <- data.frame(Time=hmin, Resolution=names(hmin))

ldply(seq(1,650,1), function(h) c("h" = h, sqrt(vmbe(vd, vb, h)))) %>% 
  melt(measure.vars=1+seq(length(Res)), variable.name="Resolution",value.name="SDQI") -> df

ggplot(df, aes(h, SDQI, col=Resolution)) + geom_line() + thm +
  scale_y_log10(name=expression(sigma[EA[]]~"(a.u.)"), breaks=10^(-4:-1), labels=fancy_scientific) +
  scale_x_continuous(name="h (sec)", minor_breaks=seq(100,3600,100)) +
  scale_color_discrete(labels=fancy_scientific) +
  geom_vline(data=hmin, aes(xintercept=Time, col=Resolution), lty=2,show.legend = FALSE) +
  geom_text(data=hmin, aes(x=Time+20, y=5e-5, label=round(Time)),show.legend = FALSE) -> p.ST
#### slope variation from rho ####
Rho.ini = Rho.on
rate = 1.157e-6

ldply(seq(3600), function(n){
  Rho <- c(Rho.ini , Rho.ini * cumprod(1 - rnorm(n-1, mean = 1*rate, sd = 1*rate)))
  lam = Sigma0*d*nu*Rho
  c("h" = n, "vb" = var(lam))
}) -> df

ggplot(df, aes(h,vb)) + geom_line() + thm + labs(x="h (sec)", y=expression(sigma[beta]^2~(sec^-2))) +
  scale_y_log10(limits=c(1e-24, 1e-16), labels=fancy_scientific) -> p.vB.Rho

#### BEAM TIME STATISTICS + SYSTEMATICS ####
vb = 1e-15
vd = Res^2; names(vd) <- Res
hmin = (24*Tint* vd/vb)^(1/3)

ldply(1:length(hmin), function(i) 
  data.frame("StateTime" = hmin[i]/60, "Resolution" = names(hmin[i]), "Precision" = Prec,
             "BeamTime" = BeamTime(as.numeric(names(hmin[i])), Tint, C, hmin[i], Prec)/3600/24)
  ) %>% mutate(Term="Stat+Syst") -> dfd

df <- rbind(df0,dfd)

ggplot(df, aes(Precision, BeamTime, col=Resolution,linetype=Term)) + geom_line() + thm +
  scale_y_log10(name="H (days)", minor_breaks=c(1:10%o%10^(-6:+6)), labels=fancy_scientific) +
  scale_x_log10(name=expression(sigma[EA[]]~"(a.u.)"), minor_breaks=Prec, labels=fancy_scientific) +
  scale_color_discrete(labels=fancy_scientific) -> p.BT.both

print(p.ST)
print(p.vB.Rho)
print(p.BT.both)
