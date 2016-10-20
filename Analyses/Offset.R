library(dplyr); library(plyr)
library(parallel); library(doParallel)
library(strucchange)

source("./Tests/DataPrep.R")
source("DataAna.R")

########## definitions ##################


##############################################

registerDoParallel(detectCores())

loadData(from = "2016-6-30 02:00:00", to = "2016-6-30 08:00:00") -> base

base %>% ggplot(aes(Clock, BCT2)) + geom_line() + 
  geom_smooth(method="loess") +
  theme_minimal() + ggtitle("Offset")

diff(base$BCT2) -> bct2d

psdcore(base$BCT2%>%diff, ntaper=100, X.freq=1) %>% plot
pacf(bct2d, lag=5000) -> p
p$acf -> pcoef
which(abs(pcoef[40:length(pcoef)]) > qnorm(1-1e-2/2)/sqrt(length(bct2d))) %>% diff %>% table %>% plot


################################
Data %>% filter(Unit==1) -> .sub
.sub %>% filter(FSgl=="T",eCool=="Acc", FABS == "F") -> .sub.off
lm(log(BCT2)~I(UTS-UTS[1]), data=.sub.off) -> m.off
.sub.off %>% with(plot(log(BCT2)~I(UTS-UTS[1]))); abline(m.off, col="red")
m.off$residuals -> res
slice(base, 6000:6248) -> tbase

lm(log(BCT2+tbase$BCT2)~I(UTS-UTS[1]), data=.sub.off) -> tm.off
.sub.off %>% with(plot(log(BCT2+tbase$BCT2)~I(UTS-UTS[1]))); abline(m.off, col="red"); abline(tm.off, col="blue")
