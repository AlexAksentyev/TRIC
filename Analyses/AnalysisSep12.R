source("Utilities.R")
source("EstiStats.R")
source("DataPrep.R")
source("DataAna.R")

library(dplyr); library(plyr)
library(splitstackshape)
library(mosaic)
library(beanplot)
library(parallel); library(doParallel)

registerDoParallel(detectCores())

# this should go to testing

get2012Data() -> Data12
# slopes <- Data12$Slopes; Base <- Data12$Data %>% filter(FSgl=="F")
# Data <- Data12$Data %>% filter(FSgl=="T")

## renaming units into consecutive order
# Data %>% daply(.(Run, Unit), function(u) dim(u)[1]) %>%c %>% na.remove() -> uleng
# Data %>% mutate(Unit = rep(1:length(uleng), uleng)) -> Data

if(FALSE){
  #### Analysis of Base ####
  Base %>% dlply("Run") -> Base
  
  list(Small = Base[1:3] %>% ldply(.id="Run"), Big = Base[4:length(Base)] %>% ldply(.id="Run")) %>%
    ldply(function(x) .outliers(x, what="BCT2"), .id="Size") -> Base
  
  Base %>% filter(Size=="Big", FOut=="F") %>% 
    mutate(
      Clock = as.POSIXct(UTS, origin="1970-1-1"), 
      Group = derivedFactor("Good" = format(Clock, format="%H:%M") < "06:00", .default = "Bad")
    ) -> longbase
  
  longbase %>% ggplot(aes(Clock, BCT2)) + 
    geom_point(aes(col=Run)) + 
    geom_smooth(aes(linetype=Group), method="gam") + 
    theme_minimal() + ggtitle("Long Run Base")
  
  ## linear model stats
  longbase %>% ddply("Group", function(g) mutate(g, BCT2 = BCT2 + 2318) %>% fit.) -> .stats; print(.stats)
  
  ## testing assumption
  library(strucchange)
  longbase %>% filter(Run==970) %>% mutate(BCT2 = BCT2 + 2318) -> tlbase
  tlbase %>% fit.  -> .stats; print(.stats)
  tlbase %>% ggplot(aes(Clock,log(BCT2))) + geom_point() + geom_smooth(span=1.5) + 
    annotate("text", x=tlbase$Clock[5], y=7.740, label = paste("Slope =", .stats[1], "+-", .stats[2])) +
    ggtitle("Offset #970")
  
  a = rep(c(+45,0,-45)*1, each=350); b = rnorm(3, mean=.stats[[2]], sd=.stats[[3]]); names(b) <- c("U","N","D")
  x=0:349; f = t(0*exp(b%o%x)+1*1) %>% c
  
  Data %>% filter(Run==tlbase$Run[1]) %>% slice(45:3194) %>% 
    mutate(BCT2 = BCT2 + 2318, BCT2m = rep(f,3)*(BCT2 + rep(a,3)), Clock = as.POSIXct(UTS, origin="1970-1-1")) -> tRun
  
  
  tRun %>% fit. %>% formatC(2,format="e") -> .stats; print(.stats)
  tRun %>% slice(45:(45+nrow(tlbase))) %>% ggplot(aes(Clock, log(BCT2))) + geom_point()+ geom_smooth(span=1.5) + 
    annotate("text", x=tRun$Clock[50], y=9.654, label = paste("Slope =", .stats[1], "+-", .stats[2])) +
    ggtitle("Signal #970")
  
  nls(BCT2 ~ I0*exp(beta*uts) + a*exp(alpha*uts), 
      data=mutate(tRun, uts=UTS-UTS[1]), 
      start=list(I0=tRun$BCT2[5], beta=-1e-4, a=tlbase$BCT2%>%median, alpha=-1e-2)
  ) -> trunnls
  
  
  
  par(mfrow=c(3,1))
  tRun %>% slice(1:1049) %>% with({plot(log(BCT2m)~Clock, type="l",col="red"); lines(log(BCT2)~Clock)})
  Fstats(BCT2~UTS, data=tRun%>%slice(1:4049)) %>% plot()
  Fstats(BCT2m~UTS, data=tRun%>%slice(1:4049)) %>% plot(col="red")
  
  
  #### ARIMA FITTING ####
  Data %>% dlply(.(Size, Run,Unit), extractARMA, .parallel=FALSE) -> arimastr
  
  Base %>% dlply(.(Size, Run,Unit), function(u) auto.arima(u$BCT2)$coef%>%t%>% data.frame(), .parallel=FALSE) %>% 
    ldply(rbind.fill) -> barimastr
}


#### Cross section ####
source("Parameters.R")
thick = 6.92e13

Data <- Data12$Data
Data %>% filter(Run %in% c(967:976), Cell=="In") %>% mutate(Run = as.numeric(Run), uts = UTS-UTS[1])-> Data

plot(range(Data$uts), with(Data, c(min(BCT2), max(BCT2))), 
     type="n", xlab="Time, sec", ylab = "Average Current, ADC"
); legend("topleft", lty=1, col=c("blue","red"), legend=paste("Target",c("off", "on")))
Data %>% d_ply("Run", function(r) r%>%with(lines(BCT2~uts, col=c("Chopper" = "blue", "On" = "red")[Targ[1]])))

base = Data %>% filter(FSgl=="F") %>% .markOutliers("BCT2") %>% filter(FOut=="F")

base %>% filter(Run < 972) -> base1; rlm(BCT2~Run, data=base1)-> m1
base %>% filter(Run > 972) -> base2; rlm(BCT2~Run, data=base2)-> m2

x = rbind(m1$model, m2$model) %>% mutate(Group="Natural")
y = data.frame(
  Run = c(969, 971)%>%rep(each=3), 
  BCT2 = c(
    predict(m1, data.frame(Run=c(969)), interval="predict",level=.68), 
    predict(m1, data.frame(Run=c(971)), interval="predict",level=.68)
    # predict(m2, data.frame(Run=c(973)), interval="predict",level=.68)
  ), 
  FSgl = "F" %>% rep(each=3)
)
y %>% dplyr::select(-FSgl) %>% mutate(Group="Fit") %>% rbind(x) %>% 
  group_by(Run) %>% dplyr::summarise(RMN=median(BCT2), SE = sd(BCT2)/sqrt(n()), GRP=Group[1]) %>% 
  ggplot(aes(Run, RMN, col=GRP)) + geom_pointrange(aes(ymin=RMN-SE, ymax=RMN+SE)) + theme_minimal() + 
  ggtitle("Offset model") + labs(y="Median offset")


y %>% rbind.fill(Data) -> Data

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

dlply(slopes%>%filter(FOut=="F"), "Targ")[c("Chopper","On")] %>% dbeta.(.fields=c("Clock", "Size", "Run")) %>% 
  mutate(
    Estimate = Estimate/nu/thick*1e27, SE = SE/nu/thick*1e27, 
    TDiff = (as.numeric(Clock.Chopper)-as.numeric(Clock.On))/60,
    Group = derivedFactor(
      "<20 mins" = abs(TDiff) <= 20, 
      "<1.5 hours" = abs(TDiff) <= 1.5*60,
      .default = "Far",
      .method="first"
    ),
    Size = derivedFactor(
      "Small" = Size.Chopper=="Small"&Size.On=="Small",
      "Big" = Size.Chopper=="Big"&Size.On=="Big",
      .default = "Mixed"
    )
  ) %>%  
  ddply("Size", function(.sub) .sub %>% mutate(Weight = SE^(-2)/max(SE^(-2)))) %>%
  arrange(abs(TDiff)) -> cs0mb

cs0mb %>% filter(Size != "Mixed") %>% mutate(Goodness = derivedFactor(
  "Good" = (Size=="Big"&Group=="<1.5 hours")|(Size=="Small"&Group=="<20 mins"),
  .default = "Bad"
))%>% 
  ggplot(aes(x=abs(TDiff)/60, y=Weight, col=Goodness, shape=Size)) + 
  geom_point() + 
  facet_grid(Group~., scale="free_x", space="free_x") + 
  labs(x="Time difference in hours") + theme(legend.position="top")

cs0mb %>% filter(Size=="Big", Group=="<1.5 hours") %>% ggplot(aes(x=Clock.Chopper, y=Estimate.Chopper)) + 
  geom_point()

cs0mb %>% filter(Size=="Big", Group=="<1.5 hours") %>% arrange(Clock.Chopper) %>% slice(1:4) %>%
  mutate(TDiff = abs(TDiff)) -> cs0mb.good

cs0mb.good %>% mutate(Close = Group, Group = derivedFactor("Norm" = Run.Chopper<972,"Weird" = Run.Chopper>972)) -> cs0mb.good
cs0mb.good %>%
  group_by(Close, Group) %>% dplyr::summarise(
    NUM = n(),
    MEAN = mean(Estimate),
    W.MEAN = weighted.mean(Estimate),
    SD = sd(Estimate),
    SE = SD/sqrt(NUM)
  ) -> cs0mb.good.summary

write.table(x=cs0mb.good.summary, file='summary_table.ods',sep='\t',row.names=F)

cs0mb.good %>% ggplot(aes(col=Close, shape=Group)) + 
  geom_pointrange(aes(x=Run.Chopper, y=Estimate, ymin=Estimate-SE, ymax=Estimate+SE)) + theme_minimal() + 
  geom_text(labels=)

Data %>% filter(Run %in% c(967,968,969,970,971)) -> Data.good
Data.good %>% d_ply("Run", function(.run){
  .run %>% fit.(TRUE, .legend = c("Estimate","SE","Chi2", "I0"), .subset=c(45,15), .obs="Run")
})

Data.good %>% d_ply("Run", function(.run){efp(log(BCT2)~UTS, data=.run[45:(dim(.run)[1]-15),], type="ME") -> .efp; .efp %>% plot; sctest(.efp)$p.value %>% cat})
