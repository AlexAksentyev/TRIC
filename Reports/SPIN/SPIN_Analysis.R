source("EstiStats.R")
source("DataAna.R")
source("Parameters.R")

library(dplyr); library(plyr)
library(parallel); library(doParallel)
library(ggfortify)

registerDoParallel(detectCores())

#### getting data ####
get2012Data() -> Data12
slopes12 <- Data12$Slopes; Data12 <- Data12$Data
get2016Data() -> Data16
slopes16 <- Data16$Slopes; Data16 <- Data16$Data

#### plotting the cycles ####
Data12 %>% ggplot(aes(Clock, BCT2, col=Targ)) + 
  scale_color_manual(name="Target state", breaks=c("Chopper","On"), labels=c("Off","On"), values=c("black", "red")) + 
  geom_point() + 
  theme_minimal() + theme(legend.position = "top", legend.title=element_text()) +
  labs(x="Time (local)", y="I (a.u.)")

Data16 %>% ggplot(aes(Clock, BCT2, col=B.Spin)) + 
  scale_color_manual(name="Beam spin", breaks=c("U","D", "N"), labels=c("Up","Down","Null"), values=c("red", "black","blue")) + 
  geom_point() + 
  theme_minimal() + theme(legend.position = "top", legend.title=element_text()) +
  labs(x="Time (local)", y="I (a.u.)")

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

autoplot(m, which=1) + theme_minimal() + ggtitle("") + labs(y=expression(ln~I[i]~-~(hat(alpha)~+~hat(beta)~t)), x=expression(hat(alpha)~+~hat(beta)~t))
autoplot(pacf(m$residuals, lag=600, main="", plot=FALSE)) + theme_minimal() + labs(y="Partial ACF")


library(lmtest); library(strucchange); library(compute.es)
f = log(BCT1)~uts
Run969 <- slice(Run969, .sub[1]:.sub[2])

.es <- function(ht, type){
  stat = ht$statistic; n1 = ht$parameter; if(length(n1)>1) n2 <- n1[2] else n2 <- n1; n1 <- n1[1]
  match.fun(type)(stat, n1, n2, verbose=FALSE) -> tht
  print(c("Cohen's d" = tht$d, "P-value" = tht$pval.d));cat("\n")
}

harvtest(f, data=Run969) %>% .es("tes")
raintest(f, data=Run969) %>% .es("fes")
Fstats(f, data=Run969) %>% .es("fes")
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
filter(slopes16, B.Spin == "N", T.Spin==1) %>% mutate(Run = 1) %>% dlply("FABS") %>% cs0mb.(1.1e14) -> cs0mb16

ggplot(cs0mb12, aes(Run.Chopper, Estimate, shape = Soundness, col = Closeness)) + 
  geom_pointrange(aes(ymin=Estimate-SE, ymax=Estimate+SE)) + 
  theme_minimal() + theme(legend.position="top") + 
  scale_color_manual(values=c("black","red"))+
  labs(x="Off-cycle", y=expression(hat(sigma)[0]~"(a.u.)"))

ggplot(cs0mb12, aes(Estimate, col=Closeness)) +
  facet_grid(Soundness~., scale="free_y", space="free_y") +
  scale_color_manual(breaks=c("Close","Far"), values=c("black","red")) + 
  geom_density(trim=TRUE, kernel="rect") + 
  # geom_segment(aes(x=Estimate, xend=Estimate, y=0, yend=.00005, col=Closeness)) +
  theme_minimal() + labs(x=expression(hat(sigma)[0]~"(a.u.)"), y="Kernel density") + theme(legend.position="top")

cs0mb12 %>% group_by(Soundness, Closeness) %>% 
  dplyr::summarise(
    NUM = n(),
    MEAN = mean(Estimate),
    W.MEAN = weighted.mean(Estimate, SE^(-2)),
    SD = sd(Estimate),
    SE = SD/sqrt(NUM)
  ) 


ggplot(cs0mb16, aes(Estimate, col=Closeness)) +
  facet_grid(Soundness~., scale="free_y", space="free_y") +
  scale_color_manual(breaks=c("Close","Far"), values=c("black","red")) + 
  geom_density(trim=TRUE, kernel="rect") + 
  # geom_segment(aes(x=Estimate, xend=Estimate, y=0, yend=.00005, col=Closeness)) +
  theme_minimal() + labs(x=expression(hat(sigma)[0]~"(a.u.)"), y="Kernel density") + theme(legend.position="top")

cs0mb16 %>% group_by(Soundness, Closeness) %>% 
  dplyr::summarise(
    NUM = n(),
    MEAN = mean(Estimate),
    W.MEAN = weighted.mean(Estimate, SE^(-2)),
    SD = sd(Estimate),
    SE = SD/sqrt(NUM)
  ) 

rbind(mutate(cs0mb12, Year="2012"), mutate(cs0mb16, Year="2016")) -> cs0mb

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
slopes12 %>% mutate(
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
  labs(y=expression(hat(beta)))

slopes12 %>% group_by(Targ) %>% 
  dplyr::summarise(
    NUM = n(),
    MEAN = mean(Estimate),
    W.MEAN = weighted.mean(Estimate),
    SD = sd(Estimate),
    SE = SD/sqrt(NUM)
  )

ggplot(slopes16, aes(I0, Estimate, col=B.Spin, shape=FABS)) + geom_point() + 
  facet_grid(B.Spin~.) + theme_minimal()

ggplot(slopes12, aes(I0, Estimate, col=Targ)) + geom_point() + geom_text(aes(label=Run, vjust=1.35, hjust=-.25))+
  scale_y_continuous(labels=fancy_scientific)+
  facet_grid(Targ~., scale="free_y", space="free_y", label=as_labeller(c("Chopper"="Off", "On"="On"))) + 
  theme_minimal() + theme(legend.position="top") +
  scale_color_manual(name="Target state", breaks=c("Chopper","On"), labels=c("Off","On"), values=c("black", "red")) +
  labs(x=expression(I[0]~"(a.u.)"), y=expression(hat(beta)))
  
#### Ayy estimation ####
thick = 1.1e14; dP = diff(Pb); cs0est = (cs0mb16 %>% filter(Soundness=="Sound", Closeness=="Close") %>% WMN)*1e-27
slopes16 %>% filter(FABS=="T", B.Spin != "N") %>% dlply("B.Spin") %>% dbeta. %>% 
  mutate(Estimate = Estimate/(dP*nu*Pt*thick*cs0est), SE = -SE/(dP*nu*Pt*thick*cs0est)) -> Ayy

ggplot(Ayy, aes(Estimate)) + geom_density(kernel="gaus") +
  theme_minimal() + labs(x=expression(hat(A)[y,y]~"(a.u.)"), y="Kernel density")
