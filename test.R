source("DataAna.R")
get2016Data() -> Data16
slopes16 <- Data16$Slopes; Data16 <- Data16$Data

source("Parameters.R")
DP = -diff(Pb); SP = sum(Pb)
b1U = filter(slopes16, FABS=="T", B.Spin=="U"); b1D = filter(slopes16, FABS=="T", B.Spin=="D")
lbeta = list("D" = b1D, "U" = b1U)
D = dbeta.(lbeta); d = D$Estimate
R = rbeta.(lbeta); r = R$Estimate; rr = 2*r/Pt/(DP - SP*r); R$Rp <- rr

