rm(list = ls(all=TRUE))
.rename4plot <- function(x){
  mutate(x, 
         `Beam Spin` = factor(B.Spin, levels=c("U","D","N"), labels = c("Up","Down","Null")),
         `Target State` = factor(FABS, levels =c("F","T"), labels=c("Off","On"))
  ) 
}

###############################################################################################
library(parallel); library(doParallel); registerDoParallel(detectCores())
source("DataAna.R")
get2016Data() -> Data16
slopes16 <- Data16$Slopes; Data16 <- Data16$Data

.rename4plot(slopes16) -> slopes16
.rename4plot(Data16) -> Data16

## all cycles
ggplot(Data16, aes(Clock, BCT2, col=`Beam Spin`)) + geom_point() + 
  theme_bw() + theme(legend.position="top") + labs(y = "I (a.u.)")

Data16 <- filter(Data16, eCool=="Acc", !is.na(`Target State`), T.Spin==1)
slopes16 <- filter(slopes16, T.Spin==1)

if(FALSE){
 
  
  ## slopes vs time
  ggplot(slopes16, aes(Clock, Estimate, shape=`Target State`, col=`Beam Spin`)) + geom_point() + 
    facet_grid(`Beam Spin`~`Target State`) + 
    theme_bw() + theme(legend.position="top") + labs(y=expression(hat(beta))) +
    geom_smooth(method="lm", se=FALSE, show.legend=FALSE, size=.4, linetype=3)
  
  ## autocorrelation density curve
  Data16 %>% ddply(.(`Beam Spin`, `Target State`, Unit), function(s) dwtest(lm(log(BCT2)~I(UTS-UTS[1]), data=s))$statistic/2) -> x
  ggplot(x, aes(DW, col=`Beam Spin`)) + geom_density() + facet_grid(`Target State`~.) + 
    theme_bw() + labs(x="Autocorrelation function at lag 1") + theme(legend.position="top")
  
  ## breakpoint density curve
  library(strucchange)
  f = log(BCT2) ~ uts
  Data16 %>% ddply(.(`Beam Spin`, `Target State`, Unit), function(s) data.frame("BP" = Fstats(f, data=mutate(s,uts=UTS-UTS[1]))$breakpoint)) -> x
  ggplot(x, aes(BP, col=`Beam Spin`)) + geom_density() + facet_grid(`Target State`~.) +
    theme_bw() + theme(legend.position="top") + labs(x = "Breakpoint (seconds from start)")
  
  ## linear model fit ####
  filter(Data16, Unit==17) %>% filter(eCool=="Acc",`Target State`=="On") %>% mutate(uts = UTS-UTS[1]) -> TRun
  ggplot(TRun, aes(Clock, BCT2)) + geom_line() + theme_bw() + labs(y="I (a.u.)")
  library(ggfortify)
  lm(f, data = TRun) -> m
  autoplot(m,1) + ggtitle("") + theme_bw() + labs(x = expression(hat(y)~"="~ hat(alpha) + hat(beta)*t), y=expression(y - hat(y)))
  library(lmtest)
  lmtest::bptest(m)
  lmtest::dwtest(m)
  
  ## slopes vs current
  ggplot(slopes16, aes(I0, Estimate, col=`Beam Spin`, shape=`Target State`)) + geom_point() +
    facet_grid(`Beam Spin`~.) +
    theme_bw() + theme(legend.position="top") + labs(x=expression(I[0]~"(a.u.)"), y=expression(hat(beta))) +
    geom_smooth(method="lm",se=FALSE,show.legend=FALSE, aes(linetype=`Target State`), size=.4)
}

## estimates ####
.sumstat <- function(x){
  require(modeest)
  group_by(x, Sound, Close) %>% 
    dplyr::summarise(NUM=n(), WMN = weighted.mean(Estimate, SE^-2), SE_d = sd(Estimate)/sqrt(NUM),
                     Chi2 = sum((Estimate - WMN)^2/SE^2)/(NUM-1),
                     SE_t = sqrt(1/sum(SE^-2))*Chi2,
                     MODE = mlv(Estimate, method="mfv")[["M"]]
    ) %>% print
  group_by(x, Sound) %>% 
    dplyr::summarise(NUM=n(), WMN = weighted.mean(Estimate, SE^-2), SE_d = sd(Estimate)/sqrt(NUM),
                     Chi2 = sum((Estimate - WMN)^2/SE^2)/(NUM-1),
                     SE_t = sqrt(1/sum(SE^-2))*Chi2,
                     MODE = mlv(Estimate, method="mfv")[["M"]]
    ) %>% print
}

## cross section
source("Parameters.R")
thick=1.1e14
filter(slopes16, B.Spin =="N") %>% dlply("`Target State`") %>% dbeta.(c("FOut", "Clock")) %>% 
  mutate(Estimate = Estimate/nu/thick*1e27, SE = SE/nu/thick*1e27, 
         Sound = derivedFactor("No" = FOut.On=="T"|FOut.Off=="T", .default="Yes"),
         Close = derivedFactor("Yes" = abs(as.numeric(Clock.Off)- as.numeric(Clock.On))/60 <= 40, .default="No")
        ) -> cs0mb

ggplot(cs0mb%>%filter(Sound=="Yes"), aes(Estimate,col=Close, alpha=.2)) + geom_density(kernel="rect") + theme_bw() + 
  labs(x=expression(hat(sigma)[0]~"(a.u.)"), y="Rectangular kernel density estimate") + guides(alpha=FALSE)
.sumstat(cs0mb)

filter(cs0mb, Sound=="Yes") %>% daply("Close", function(s) WMN(s)*1e-27) -> cs0est

## asymmetry
a = -nu*cs0est*thick*Pt*diff(Pb[1:2])
(filter(slopes16, FABS=="T") %>% dlply("B.Spin"))[c("D","U")] %>% dbeta.(c("FOut","Clock")) %>%
  mutate(Sound = derivedFactor("No" = FOut.D=="T"|FOut.U=="T", .default="Yes"),
         Close = derivedFactor("Yes" = abs(as.numeric(Clock.D)- as.numeric(Clock.U))/60 <= 40, .default="No")
      ) %>% ddply("Close", function(s) mutate(s, Estimate = Estimate/a[Close[1]], SE = SE/a[Close[1]])) -> Ayy

ggplot(Ayy%>%filter(Sound=="Yes"), aes(Estimate,col=Close, alpha=.2)) + geom_density(kernel="rect") + theme_bw() + 
  labs(x=expression(hat(A)[yy]~"(a.u.)"), y="Rectangular kernel density estimate") + guides(alpha=FALSE)
.sumstat(Ayy)


## testing if current-dependent error fit fits better ####
# library(nlme)
# gls(f, data=TRun, correlation = corAR1(.6, ~uts)) ->mgls
