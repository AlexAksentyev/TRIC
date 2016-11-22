source("DataAna.R")
get2016Data() -> Data16
slopes16 <- Data16$Slopes; Data16 <- Data16$Data

source("Parameters.R")
DP = -diff(Pb); SP = sum(Pb)
b1U = filter(slopes16, FABS=="T", B.Spin=="U"); b1D = filter(slopes16, FABS=="T", B.Spin=="D")
lbeta = list("D" = b1D, "U" = b1U)
D = dbeta.(lbeta); d = D$Estimate
R = rbeta.(lbeta); r = R$Estimate; rr = 2*r/Pt/(DP - SP*r); R$Rp <- rr

##
filter(Data16, Unit==17) -> TRun
ggplot(TRun, aes(Clock, BCT2)) + geom_line() + theme_bw()
TRun %>% filter(eCool=="Acc",FABS=="T") %>% mutate(uts = UTS-UTS[1]) -> TRun
f = log(BCT2) ~ uts
library(lmtest)
lm(f, data = TRun) -> m
lmtest::bptest(m)
lmtest::dwtest(m)[c("statistic", "p.value")]
library(strucchange)
Fstats(f, data=TRun) %>% plot

Data16 %>% ddply(.(B.Spin, FABS, Unit), function(s) Fstats(f, data=s)$breakpoint) -> x; names(x)[4] <- "BP"
ggplot(x, aes(BP, col=B.Spin)) + geom_density() + facet_grid(FABS~.) + 
  theme_bw() + theme(legend.position="top") + labs(x = "Breakpoint (seconds from start)")

library(forecast)
Data16 %>% filter(eCool=="Acc", !is.na(FABS)) %>% 
  dlply(.(B.Spin, FABS, Unit), function(s) auto.arima(log(s$BCT2))%>%coef %>%t %>% as.data.frame() %>%
          cbind(B.Spin = s$B.Spin[1], FABS = s$FABS[1], Unit = s$Unit[1])) %>% 
  rbind.fill() -> AAMC

Data16 %>% ddply(.(B.Spin, FABS, Unit), function(s) dwtest(lm(log(BCT2)~I(UTS-UTS[1]), data=s))$statistic/2) -> x
ggplot(x, aes(DW, col=B.Spin)) + geom_density() + facet_grid(FABS~.) + theme_bw() + labs(x="Correlation at lag 1")

library(dynlm)
dynlm(Y ~ trend(Y), data=mutate(TRun, Y = log(BCT2)))%>%summary
lm(Y ~ uts, data = mutate(TRun, Y = log(BCT2), uts = UTS-UTS[1])) %>% summary
