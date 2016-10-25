## computing the estimate via the cross ratio method

source("EstiStats.R")
source("DataAna.R")
source("Parameters.R")

library(dplyr); library(plyr)

get2016Data() -> Data
slopes <- Data$Slopes; Data <- Data$Data

.closeness.factor <- function(.data) .data %>% mutate(
  Closeness = derivedFactor(Close = abs(as.numeric(Clock.U-Clock.D))<= 60, .default = "Far")
)

dlply(filter(slopes, FABS == "T"), "B.Spin")[c("D","U")] %>% rbeta.("Clock") %>% .closeness.factor -> Ron
dlply(filter(slopes, FABS == "T"), "B.Spin")[c("D","U")] %>% dbeta.("Clock") %>% .closeness.factor -> Don
dlply(filter(slopes, FABS == "F"), "B.Spin")[c("D","U")] %>% rbeta.("Clock") %>% .closeness.factor -> Roff

filter(slopes12, Targ=="Valve", Bunch=="BB") -> slp.vlv
slp.vlv %>%filter(Run>=937) %>% ddply("Size",.markOutliers) -> slp.vlv
slp.vlv%>%filter(FOut=="F")%>%with(boxplot(Estimate~Size))

list(Vlv = slp.vlv, On = filter(slopes, FABS=="T",B.Spin=="N")) %>% rbeta. -> Rx
x = mutate(Rx, SE = .5*SE/Estimate^2, Estimate = -.5*(Estimate+1)/Estimate) %>%
  .markOutliers() %>% filter(FOut=="F")
mx = WMN(x); vx = .se(x)^2

rbind(
  Ron %>% dplyr::select(Estimate, SE, Closeness) %>% mutate(Targ = "On"),
  Roff %>% dplyr::select(Estimate, SE, Closeness) %>% mutate(Targ = "Off")
) -> R

PtP = Pt*mean(abs(Pb))

Ayy.R = mutate(R, SE = sqrt(Estimate^2*vx + mx^2*SE^2 + vx*SE^2)/PtP, Estimate = Estimate/PtP*(1+mx), Stat="R")
a = 2*PtP*nu*d*Rho.on*500e-27
Ayy.D = transmute(Don, Estimate = Estimate/a, SE = SE/a, Targ="On", Stat="D", Closeness = Closeness)

Ayy = rbind(Ayy.D, Ayy.R)

Ayy %>% ggplot(aes(col=Stat, linetype=Targ)) + geom_density(aes(x=Estimate)) + 
  theme_minimal() + ggtitle("Ayy-statistic")

Ayy %>% group_by(Targ, Stat) %>% dplyr::summarise(
  NUM = n(),
  Mean = mean(Estimate),
  WMN = weighted.mean(Estimate, SE^(-2)),
  SD = sd(Estimate),
  SE = SD/sqrt(NUM),
  RSE = SE/WMN,
  Chi2 = (sum(Estimate-WMN)/SD)^2/NUM,
  Chi2Exc = Chi2-1
) %>% write.table(file="./Stats/Ayy.ods", quote=FALSE, row.names=FALSE)
