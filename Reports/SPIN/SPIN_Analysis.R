source("EstiStats.R")
source("DataAna.R")
source("Parameters.R")

library(dplyr); library(plyr)
library(parallel); library(doParallel)
library(ggfortify); library(cowplot)

registerDoParallel(detectCores())

lblfnt = 20

#### getting data ####
get2012Data() -> Data12
slopes12 <- Data12$Slopes; Data12 <- Data12$Data
get2016Data() -> Data16
slopes16 <- Data16$Slopes %>% filter(T.Spin==1); Data16 <- Data16$Data %>% filter(T.Spin==1)

#### plotting the cycles ####
Data12 %>% transmute(
  Year = 2012, Clock, Unit = Run, 
  Target = derivedFactor("Off" = Targ=="Chopper",.default="On"),
  `Beam Spin` = "Null", BCT2
) %>% rbind(
  transmute(
    Data16, Year = 2016, Clock, Unit,
    Target = derivedFactor("Off" = FABS=="F", "On" = FABS=="T", .default = NA),
    `Beam Spin` = derivedFactor("Up" = B.Spin=="U", "Down" = B.Spin=="D", "Null" = B.Spin=="N", .default = NA),
    BCT2
  )
) -> Data

ggplot(filter(Data, Year==2012), aes(Clock, BCT2, col=Target)) + geom_point() + labs(y="I (a.u.)") +
  theme_bw() + 
  theme(legend.position="top", axis.text=element_text(size=lblfnt), axis.title=element_text(size=lblfnt)) + 
  scale_color_manual(values=c("black","red")) -> p12
ggplot(filter(Data, Year==2016), aes(Clock, BCT2, col=`Beam Spin`)) + geom_point() + labs(y="I (a.u.)") +
  theme_bw() + 
  theme(legend.position="top", axis.text=element_text(size=lblfnt), axis.title=element_text(size=lblfnt)) + 
  scale_color_manual(values=c("blue","black","red")) -> p16
plot_grid(p12, p16, ncol=1, labels=c("2012","2016"), label_size=lblfnt)

##################################
## shows that BCT2 = 0 is just stack overflow ##
# Data%>%filter(Run==971) %>% slice(1:60) -> xcheck
# ggplot(xcheck, aes(x=Clock)) + geom_point(aes(y=BCT1, col="BCT1")) + geom_point(aes(y=BCT2/10, col="BCT2/10")) +
#   geom_vline(xintercept=xcheck$UTS[4], col="black") +
#   theme_minimal() + ggtitle("Run 971; first minute") + labs(y="Current, ADC")
#################################

#### testing cycle for statistical properties ####
.sub=c(45,15)
Data12 %>% filter(Run==969) %>% mutate(uts=UTS-UTS[1]) -> Run969
b = filter(Run969, FSgl=="F", BCT2 != 0)$BCT2 %>% median
Run969%>%mutate(BCT2=BCT2-b)%>%filter(FSgl=="T")->Run969
m = fit.(Run969, model=TRUE, .subset=.sub); .sub <- c(0, nrow(Run969))+c(1,-1)*.sub
slice(Run969, 20:nrow(Run969))%>%ggplot(aes(uts, log(BCT2))) + geom_line() + geom_smooth(method="lm", col="red", size=.75) + 
  theme_minimal() + labs(y="ln I", x="Time (seconds)") + 
  geom_vline(xintercept=Run969$uts[.sub], col="gray",linetype=2)

df = data.frame(Res = residuals(m), Fit = fitted(m))
pacf(df$Res, lag=60, plot = FALSE) -> x
xdf = data.frame(Lag = x$lag, pACF = x$acf)

ggplot(df, aes(Fit, Res)) + geom_point() + 
  geom_smooth(method="loess", span=.8, se=FALSE, col="red") + 
  geom_hline(yintercept=0, linetype=2, col="gray") + 
  theme_bw() + 
  theme(legend.position="top", axis.text=element_text(size=lblfnt), axis.title=element_text(size=lblfnt)) + 
  labs(x=expression(hat(alpha) + hat(beta)*t), y=expression(ln~I[t] - group("(",list(hat(alpha) +hat(beta)*t),")") )) -> prf
ggplot(xdf, aes(Lag, pACF)) + geom_bar(stat="identity", width=.3) + 
  geom_hline(yintercept=c(-1,+1)*1.96/sqrt(x$n.used), col="blue", linetype=2) +
  theme_bw() + 
  theme(legend.position="top", axis.text=element_text(size=lblfnt), axis.title=element_text(size=lblfnt)) + 
  labs(y="PACF") -> ppacf
plot_grid(prf, ppacf, ncol=1, labels=c("a)","b)"), label_size=lblfnt)

library(lmtest); library(strucchange); library(compute.es)
f = log(BCT1)~uts
Run969 <- slice(Run969, .sub[1]:.sub[2])

harvtest(f, data=Run969) 
raintest(f, data=Run969) 
sctest(f, data=Run969, type="aveF")
sctest(f, data=Run969, type="ME")
bptest(f, data=Run969) 
dwtest(f, data=Run969)

#### cross section ####
cs0mb. <- function(.slp, thick){
  .slp %>% dbeta.(.fields=c("Clock", "FOut", "Run")) -> db
  names(db) %>% (function(x) gsub(".F", ".Chopper", x)) %>% (function(x) gsub(".T",".On", x)) -> names(db)
  db %>% mutate(
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
    )
}
dlply(slopes12, "Targ")[c("Chopper","On")] %>% cs0mb.(6.92e13) -> cs0mb12
filter(slopes16, B.Spin == "N") %>% mutate(Run = 1) %>% dlply("FABS") %>% cs0mb.(1.1e14) -> cs0mb16

.sumstat <- function(x){
  x %>% group_by(Soundness, Closeness) %>% 
    dplyr::summarise(
      NUM = n(),
      MEAN = mean(Estimate),
      W.MEAN = weighted.mean(Estimate, SE^(-2)),
      SD = sd(Estimate),
      SE = SD/sqrt(NUM)
    ) %>% print
  x %>% group_by(Soundness) %>% 
    dplyr::summarise(
      NUM = n(),
      MEAN = mean(Estimate),
      W.MEAN = weighted.mean(Estimate, SE^(-2)),
      SD = sd(Estimate),
      SE = SD/sqrt(NUM)
    ) %>% print
  x %>% dplyr::summarise(
      NUM = n(),
      MEAN = mean(Estimate),
      W.MEAN = weighted.mean(Estimate, SE^(-2)),
      SD = sd(Estimate),
      SE = SD/sqrt(NUM)
    ) %>% print
}

.sumstat(cs0mb12)
.sumstat(cs0mb16)

#### slopes ####
slopes16 <- filter(slopes16, T.Spin == 1)
slopes12 %>% group_by(Targ) %>% 
  dplyr::summarise(
    NUM = n(),
    MEAN = mean(Estimate),
    W.MEAN = weighted.mean(Estimate),
    SD = sd(Estimate),
    SE = SD/sqrt(NUM)
  )

slopes12<-mutate(slopes12, B.Spin=factor(1, labels = "Null"))
ggplot(slopes12, aes(I0, Estimate)) + geom_pointrange(aes(ymin=Estimate-SE,ymax=Estimate+SE)) + 
  scale_y_continuous(labels=fancy_scientific)+
  facet_grid(B.Spin~Targ, scale="free_x", label=as_labeller(c("Chopper"="Off", "On"="On", "Null" = "Null"))) + 
  theme_bw() + theme(legend.position="top") + 
  labs(x=expression(I[0]~"(a.u.)"), y=expression(hat(beta))) -> p12
ggplot(slopes16, aes(I0, Estimate, col=B.Spin)) + geom_point() + 
  scale_y_continuous(labels=fancy_scientific) +
  scale_color_manual(name="Beam spin", breaks=c("U","D", "N"), labels=c("Up","Down","Null"), values=c("red", "blue","black")) + 
  geom_smooth(method="lm", se=FALSE, show.legend=FALSE, size=.4) +
  facet_grid(B.Spin~FABS, scale="free_x",
             label=as_labeller(c("F" = "Off", "T" = "On", "U" = "Up","D"="Down","N"="Null"))) + 
  theme_bw() + theme(legend.position="none") + 
  labs(y=expression(hat(beta)), x=expression(I[0]~"(a.u.)")) -> p16
plot_grid(p12,p16,ncol=1, labels=c("2012","2016"))


pull <- function(x,y){
  (mean(x$Estimate)-mean(y$Estimate))/sqrt(.se(x)^2+.se(y)^2) %>% abs -> .stat
  pt(.stat, 12, lower.tail = FALSE)
}

pull(filter(slopes16, B.Spin == "D"), filter(slopes16, B.Spin == "N"))

#### Ayy estimation ####
thick = 1.1e14; dP = -diff(Pb); cs0est = (cs0mb16 %>% filter(Soundness=="Sound", Closeness=="Close") %>% WMN)*1e-27
(slopes16 %>% filter(FABS=="T", B.Spin != "N") %>% dlply("B.Spin"))[c("D","U")] %>% dbeta. %>% 
  mutate(Estimate = Estimate/(dP*nu*Pt*thick*cs0est), SE = SE/(dP*nu*Pt*thick*cs0est)) -> Ayy

Ayy%>%WMN
ggplot(Ayy, aes(Estimate)) + geom_histogram(binwidth=.15, fill='white', col='black') + 
  theme_bw() + labs(x=expression(hat(A)[yy])) +
  theme(legend.position="top", axis.text=element_text(size=lblfnt), axis.title=element_text(size=lblfnt))
