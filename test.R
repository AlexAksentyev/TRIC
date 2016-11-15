source("DataAna.R")
get2016Data() -> Data16
slopes16 <- Data16$Slopes; Data16 <- Data16$Data

source("Parameters.R")
DP = -diff(Pb); SP = sum(Pb)
b1U = filter(slopes16, FABS=="T", B.Spin=="U"); b1D = filter(slopes16, FABS=="T", B.Spin=="D")
lbeta = list("D" = b1D, "U" = b1U)
D = dbeta.(lbeta); d = D$Estimate
R = rbeta.(lbeta); r = R$Estimate; rr = 2*r/Pt/(DP - SP*r); R$Rp <- rr

## comparing the statistics' distributions
df = data.frame(Stat = "D", Estimate = d) %>% rbind(data.frame(Stat = "R", Estimate = rr))
lm(rr[order(rr)]~d[order(d)]) -> m
qqplot(d, rr, xlab = "D-stat", ylab = "R-stat"); abline(m, col="red"); abline(v=0,h=0, lty=4, col="gray")
qqplot(d*m$coefficients[2] + m$coefficients[1], rr, xlab="g(D)-stat", ylab="R-stat"); abline(0,1, col="red"); abline(v=0,h=0, lty=4, col="gray")

#### Ayy(R) -> CS0 ####
lam = (nu*1.1e14*Pt*DP)
mutate(R, Estimate = rr) -> Rp

WMN(Rp) -> AyyR
WMN(D)/lam/AyyR * 1e24; cat(AyyR)
