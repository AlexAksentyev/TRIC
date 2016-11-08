
DP = -diff(Pb); SP = sum(Pb)
b1U = filter(slopes16, FABS=="T", B.Spin=="U"); b1D = filter(slopes16, FABS=="T", B.Spin=="D")
lbeta = list("D" = b1D, "U" = b1U)
D = dbeta.(lbeta); d = D$Estimate
R = rbeta.(lbeta); r = R$Estimate; rr = 2*r/Pt/(DP - SP*r)


## comparing the statistics' distributions
df = data.frame(Stat = "D", Estimate = d) %>% rbind(data.frame(Stat = "R", Estimate = rr))
lm(rr[order(rr)]~d[order(d)]) -> m
qqplot(d, rr, xlab = "D-stat", ylab = "R-stat"); abline(m, col="red"); abline(v=0,h=0, lty=4, col="gray")
qqplot(d*m$coefficients[2] + m$coefficients[1], rr, xlab="g(D)-stat", ylab="R-stat"); abline(0,1, col="red"); abline(v=0,h=0, lty=4, col="gray")


lm(AyyR[order(AyyR)]~D_PPt[order(D_PPt)]) -> m
beta = coef(summary(m))[2,1]; alpha = coef(summary(m))[1,1]
qqplot(D_PPt, AyyR, main = "QQ-plot Ayy estimate from R vs D stats"); abline(m, col="blue")
AyyD <- D_PPt * beta

g = (AyyR-AyyD)/sqrt(var(AyyR) + var(AyyD) + 2*cov(AyyR, AyyD))
