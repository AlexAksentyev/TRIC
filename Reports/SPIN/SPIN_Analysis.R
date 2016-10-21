source("EstiStats.R")
source("DataPrep.R")
source("Parameters.R")

library(dplyr); library(plyr); library(mosaic)
library(parallel); library(doParallel)
library(ggfortify)

thick = 6.92e13

registerDoParallel(detectCores())

get2012Data() -> Data12
Data <- Data12$Data
Data %>% filter(Run %in% c(967:976), Cell=="In") %>% mutate(Run = as.numeric(Run), uts = UTS-UTS[1])-> Data

#### plotting the nine cycles ####
# plot(range(Data$uts), with(Data, c(min(BCT2), max(BCT2))), 
#      type="n", xlab="Time, sec", ylab = "Average Current, ADC"
# ); legend("topleft", bty="n", lty=1, col=c("blue","red"), legend=paste("Target",c("off", "on")))
# Data %>% d_ply("Run", function(r) r%>%with(lines(BCT2~uts, col=c("Chopper" = "blue", "On" = "red")[Targ[1]])))
##same with ggplot
Data %>% ggplot(aes(uts, BCT2, col=Targ)) + 
  scale_color_manual(breaks=c("Chopper","On"), labels=c("Off","On"), values=c("black", "red")) + 
  geom_point() + 
  theme_minimal() + theme(legend.position = "top", legend.title=element_blank()) +
  labs(x="Time, seconds", y="Average current, ADC")

base = Data %>% filter(FSgl=="F") %>% .markOutliers("BCT2") %>% filter(FOut=="F")
base %>% mutate(GRP = derivedFactor("Before" = Run<972, "After" = Run>973, .default="973")) -> base
.summary <- function(.dat) .dat%>%group_by(Run)%>%dplyr::summarise(BCT2o = median(BCT2), SE = sd(BCT2)/sqrt(n()))
bmtot.dat = .summary(base)
bmtot = lm(BCT2o~Run, data=bmtot.dat)
base%>%filter(GRP!="973")%>%ddply("GRP", function(g) lm(BCT2o~Run, data=.summary(g)) %>% .collect.stats) %>% 
  rbind.fill(bmtot %>% .collect.stats %>% t %>% as.data.frame()) 

## chi2 is best when total fit (last row); hence don't separate base into two groups to predict offsets
PredBCT2o = data.frame(
  Run = c(969, 971)%>%rep(each=3), 
  BCT2o = c(
    predict(bmtot, data.frame(Run=c(969)), interval="predict",level=.68), 
    predict(bmtot, data.frame(Run=c(971)), interval="predict",level=.68)
    # predict(bmtot, data.frame(Run=c(973)), interval="predict",level=.68)
  ), 
  FSgl = "F" %>% rep(each=3)
)

PredBCT2o %>% dplyr::select(-FSgl) %>% mutate(Group="Fit") %>% 
  group_by(Run) %>% dplyr::summarise(SE = (max(BCT2o)-min(BCT2o))/2, BCT2o=median(BCT2o), GRP=Group[1]) %>% 
  rbind(mutate(bmtot.dat, GRP="Data")) %>%
  ggplot(aes(Run, BCT2o)) + geom_pointrange(aes(ymin=BCT2o-SE, ymax=BCT2o+SE, col=GRP)) + 
  geom_smooth(method="gam", col="black", lty=3) +
  theme_minimal()  + theme(legend.position = "top", legend.title=element_blank()) + 
  scale_color_manual(values=c("black","red")) + 
  labs(y="Baseline median", x="Cycle")

#### fill the missing offset data ####
PredBCT2o %>% group_by(Run) %>% 
  dplyr::summarise(BCT2=median(BCT2o), FSgl="F") %>% 
  rbind.fill(Data) -> Data

#### fit ####
Data %>% ddply("Run", function(.sub){
  .sub %>% filter(FSgl=="F") %>% .markOutliers("BCT2") %>% filter(FOut=="F") -> b
  .sub %>% mutate(BCT2 = BCT2 - median(b$BCT2)) %>% filter(FSgl=="T") -> .sub
  .sub %>% fit.(.subset = c(45, 15)) -> .stats
  
  data.frame(
    Estimate = .stats[1], SE = .stats[2], Chi2 = .stats[3], I0 = .stats[4],
    Quality = .sub$Quality[5],
    eCool.I = .sub$eCool.I[5],
    Targ = .sub$Targ[5],
    Bunch = .sub$Bunch[5],
    Size = .sub$Size[5],
    Clock = .sub$UTS[1] %>% as.POSIXct(origin="1970-1-1")
  )
}) %>% mutate(Unit=1) %>% ddply("Targ", .markOutliers) -> slopes

#### testing ####
Data %>% filter(Run==969) -> Run969; b = Run969$BCT2[Run969$FSgl=="F"]; Run969%>%mutate(BCT2=BCT2-b)%>%filter(FSgl=="T")->Run969
m = fit.(Run969, model=TRUE, .subset=c(45,15))
ggplot(Run969,aes(Clock, log(BCT2))) + geom_line() + geom_smooth(method="lm", col="red", size=.75) + 
  theme_minimal() + labs(y="Logarithm of average current")

autoplot(m, which=1) + theme_minimal() + ggtitle("")
autoplot(pacf(m$residuals, lag=600, main="", plot=FALSE)) + theme_minimal() + labs(y="Partial ACF")


library(lmtest); library(strucchange)
f = log(BCT2)~uts

harvtest(f, data=Run969)
raintest(f, data=Run969)
sctest(f, from=45, to = nrow(Run969)-15, data=Run969, type="aveF")
sctest(f, from=45, to = nrow(Run969)-15, data=Run969, type="ME")
bptest(f, data=Run969)
dwtest(f, data=Run969)

#### cross section ####
dlply(slopes, "Targ")[c("Chopper","On")] %>% 
  dbeta.(.fields=c("Clock", "FOut", "Run")) %>% 
  mutate(
    Estimate = Estimate/nu/thick*1e27, SE = SE/nu/thick*1e27, 
    TDiff = (as.numeric(Clock.Chopper)-as.numeric(Clock.On))/60,
    Closeness = derivedFactor(
      "Close" = abs(TDiff) <= 1.5*60,
      .default = "Far",
      .method="first"
    ),
    Soundness = derivedFactor(
      "Sound" = FOut.Chopper=="F"&FOut.On=="F",
      .default = "Unsound"
    )
  ) -> cs0mb

cs0mb %>% ggplot(aes(Run.Chopper, Estimate, shape = Soundness, col = Closeness)) + 
  geom_pointrange(aes(ymin=Estimate-SE, ymax=Estimate+SE)) + 
  theme_minimal() + theme(legend.position="top") + 
  scale_color_manual(values=c("black","red"))+
  labs(x="Off-cycle", y="Cross section estimate")

cs0mb%>%group_by(Soundness, Closeness) %>% 
  dplyr::summarise(
    NUM = n(),
    MEAN = mean(Estimate),
    W.MEAN = weighted.mean(Estimate),
    SD = sd(Estimate),
    SE = SD/sqrt(NUM)
  )


#### slopes ####
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
slopes %>% mutate(
  Group = derivedFactor("Used" = format(Clock, format="%H:%M") < "06:00", .default="Unused")#,
  # Estimate = -1/Estimate, SE = SE*Estimate^2
) %>% 
  ggplot(aes(Clock, Estimate, col=Group)) + 
  scale_y_continuous(labels=fancy_scientific)+
  scale_color_manual(values=c("black","red")) +
  geom_pointrange(aes(ymin=Estimate-SE,ymax=Estimate+SE),size=.3) + 
  # geom_smooth(lty=3, col="grey", method="lm", se=FALSE) +
  theme_minimal() +
  theme(legend.position="top", legend.title=element_blank()) + 
  facet_grid(Targ~., space="free_y", scale="free_y", labeller = as_labeller(c("Chopper" = "Off", "On" = "On")),switch="y") + 
  labs(y="Slope estimate")

slopes%>%group_by(Targ) %>% 
  dplyr::summarise(
    NUM = n(),
    MEAN = mean(Estimate),
    W.MEAN = weighted.mean(Estimate),
    SD = sd(Estimate),
    SE = SD/sqrt(NUM)
  )
