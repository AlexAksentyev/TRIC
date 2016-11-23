source("DataAna.R")
get2016Data() -> Data16
slopes16 <- Data16$Slopes; Data16 <- Data16$Data

ggplot(Data16, aes(Clock, BCT2, col=B.Spin)) + geom_point() + 
  theme_bw() + theme(legend.position="top") + labs(y = "I (a.u.)") +
  scale_color_discrete(name="Beam Spin", breaks=c("U","D","N"), labels=c("Up","Down","Null"))

Data16 <- filter(Data16, eCool=="Acc", !is.na(FABS), T.Spin==1)
slopes16 <- filter(slopes16, T.Spin==1)
mutate(slopes16, 
       B.Spin = factor(B.Spin, levels=c("U","D","N"), labels = c("Up","Down","Null")),
       FABS = factor(FABS, levels =c("F","T"), labels=c("Off","On"))
) -> slopes16
names(slopes16)[2] <- "Beam Spin"; names(slopes16)[4] <- "Target State"

ggplot(slopes16, aes(Clock, Estimate, shape=`Target State`, col=`Beam Spin`)) + geom_point() + 
  facet_grid(`Beam Spin`~`Target State`) + 
  theme_bw() + theme(legend.position="top") + labs(y=expression(hat(beta))) +
  geom_smooth(method="lm", se=FALSE, show.legend=FALSE)

ggplot(slopes16, aes(I0, Estimate, col=`Beam Spin`, shape=`Target State`)) + geom_point() +
  facet_grid(`Beam Spin`~.) +
  theme_bw() + theme(legend.position="top") + labs(x=expression(I[0]~"(a.u.)"), y=expression(hat(beta))) +
  geom_smooth(method="lm",se=FALSE,show.legend=FALSE, aes(linetype=`Target State`))

source("Parameters.R")
DP = -diff(Pb); SP = sum(Pb)
b1U = filter(slopes16, FABS=="T", B.Spin=="U"); b1D = filter(slopes16, FABS=="T", B.Spin=="D")
lbeta = list("D" = b1D, "U" = b1U)
D = dbeta.(lbeta); d = D$Estimate
R = rbeta.(lbeta); r = R$Estimate; rr = 2*r/Pt/(DP - SP*r); R$Rp <- rr

##
filter(Data16, Unit==17) -> TRun
ggplot(TRun, aes(Clock, BCT2)) + geom_line() + theme_bw() + labs(y="I (a.u.)")
TRun %>% filter(eCool=="Acc",FABS=="T") %>% mutate(uts = UTS-UTS[1]) -> TRun
f = log(BCT2) ~ uts
library(lmtest)
lm(f, data = TRun) -> m
autoplot(m,1) + ggtitle("") + theme_bw() + labs(x = expression(hat(y)~"="~ hat(alpha) + hat(beta)*t), y=expression(y - hat(y)))
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
ggplot(x, aes(DW, col=B.Spin)) + geom_density() + facet_grid(FABS~.) + theme_bw() + labs(x="Autocorrelation with previous error")

library(dynlm)
dynlm(Y ~ trend(Y), data=mutate(TRun, Y = log(BCT2)))%>%summary
lm(Y ~ uts, data = mutate(TRun, Y = log(BCT2), uts = UTS-UTS[1])) %>% summary


## estimates ####
##cross section
thick=1.1e14
filter(slopes16, B.Spin=="Null") %>% dlply("FABS") %>% dbeta.(c("FOut", "Clock")) %>% 
  mutate(Estimate = Estimate/nu/thick*1e27, SE = SE/nu/thick*1e27, 
         Sound = derivedFactor("No" = FOut.On=="T"|FOut.Off=="T", .default="Yes"),
         Close = derivedFactor("Yes" = abs(as.numeric(Clock.Off)- as.numeric(Clock.On))/60 <= 40, .default="No")
        ) -> cs0mb
ggplot(cs0mb%>%filter(Sound=="Yes"), aes(Estimate,fill=Close, alpha=.2)) + geom_density(kernel="rect") + theme_bw() + 
  labs(x=expression(hat(sigma)[0]~"(a.u.)"), y="Rectangular kernel density estimate") + guides(alpha=FALSE)

group_by(cs0mb, Sound, Close) %>% 
  dplyr::summarise(NUM=n(), WMN = weighted.mean(Estimate, SE^-2), SE_d = sd(Estimate)/sqrt(NUM),
                   Chi2 = sum((Estimate - WMN)^2/SE^2)/(NUM-1),
                   SE_t = sqrt(1/sum(SE^-2))*Chi2
  )
group_by(cs0mb, Sound) %>% 
  dplyr::summarise(NUM=n(), WMN = weighted.mean(Estimate, SE^-2), SE_d = sd(Estimate)/sqrt(NUM),
                   Chi2 = sum((Estimate - WMN)^2/SE^2)/(NUM-1),
                   SE_t = sqrt(1/sum(SE^-2))*Chi2
  )
cs0est = (filter(cs0mb, Sound=="Yes",Close=="Yes") %>% WMN)*1e-27
## asymmetry
a = -nu*cs0est*thick*Pt*diff(Pb[1:2])
(filter(slopes16, B.Spin != "Null", FABS=="On") %>% dlply("B.Spin"))[c("Down","Up")] %>% dbeta.(c("FOut","Clock")) %>%
  mutate(Estimate = Estimate/a, SE = SE/a,
         Sound = derivedFactor("No" = FOut.Down=="T"|FOut.Up=="T", .default="Yes"),
         Close = derivedFactor("Yes" = abs(as.numeric(Clock.Down)- as.numeric(Clock.Up))/60 <= 40, .default="No")
      ) -> Ayy
ggplot(Ayy%>%filter(Sound=="Yes"), aes(Estimate,fill=Close, alpha=.2)) + geom_density(kernel="rect") + theme_bw() + 
  labs(x=expression(hat(sigma)[0]~"(a.u.)"), y="Rectangular kernel density estimate") + guides(alpha=FALSE)

group_by(Ayy, Sound, Close) %>% 
  dplyr::summarise(NUM=n(), WMN = weighted.mean(Estimate, SE^-2), SE_d = sd(Estimate)/sqrt(NUM),
                   Chi2 = sum((Estimate - WMN)^2/SE^2)/(NUM-1),
                   SE_t = sqrt(1/sum(SE^-2))*Chi2
  )
group_by(Ayy, Sound) %>% 
  dplyr::summarise(NUM=n(), WMN = weighted.mean(Estimate, SE^-2), SE_d = sd(Estimate)/sqrt(NUM),
                   Chi2 = sum((Estimate - WMN)^2/SE^2)/(NUM-1),
                   SE_t = sqrt(1/sum(SE^-2))*Chi2
  )