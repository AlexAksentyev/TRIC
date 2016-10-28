source("EstiStats.R")
source("DataAna.R")
source("Parameters.R")

library(dplyr); library(plyr); library(mosaic)
library(parallel); library(doParallel)
library(ggfortify)

thick = 6.92e13

registerDoParallel(detectCores())

#### preparing data for analysis ####
## fetch metadata
.runs = 967:976
where = "~/Analysis/Sep12/BCTAnalyser/results/"
readODS::read_ods(paste(where, "RunMetadata(1).ods",sep="/"), sheet=2) %>% 
  filter(Run %in% .runs) %>% mutate(
    Start = paste(day, Start, sep = " ")%>%as.POSIXct(origin="1970-1-1")-30,
    Stop = paste(day, Stop, sep = " ")%>%as.POSIXct(origin="1970-1-1")+30
  ) -> MD
for(obs in c("eCool.I", "Bunch", "Cell", "Targ", "RF.V", "BB.V", "Beam")) ## could use sapply(, factor) instead of loop
  MD[,obs] <- factor(MD[,obs])

## load data
.flds = c("WD", "Month", "Day", "Time", "Year", paste0("BCT",1:8))
loadData("CosyAll_2012-09-25_15-53-06.dat", "~/Analysis/Sep12/DevReader/Data/", 
         Fields = .flds, from=first(MD$Start), to=last(MD$Stop)
) -> Data12

## find when cycles begin
filter(Data12, BCT2==0) %>% mutate(H = hour(Clock)) %>% 
  group_by(H) %>% dplyr::summarize(Clock = first(Clock)) -> CycSrt

## correct metadata
MD %>% mutate(
  Start = CycSrt$Clock, Stop = c(Start[-1]-1, last(Stop)),
  UTS.Stt = as.numeric(Start),
    UTS.Stp = as.numeric(Stop),
    Size = UTS.Stp-UTS.Stt+1,
    Quality = derivedFactor("bad" = Quality=="bad", "ok" = is.na(Quality))
  ) %>% dplyr::select(-day,-Beam, -Comment) -> MD

MD %>%expandRows("Size") -> MD
rownames(MD) %>% strsplit(".", fixed=TRUE) %>% 
  laply(function(e) {e <- as.numeric(e); if(length(e) <2) 0 else e[2]}) -> addage
MD %>% mutate(UTS = UTS.Stt + addage) -> MD

## add metadata to cycles
MD %>% dplyr::select(Run, UTS, eCool.I, Bunch, RF.V, BB.V, Cell, Targ, Quality) %>% 
  join(Data12, by="UTS", type="right") %>% filter(!is.na(Run)) -> Data12 ## offset before the first cycle is discarded

## classify signal & base, exclude wrong cycle
.dumbCut <- (function(.data, x) mutate(.data, FSgl = derivedFactor("F" = BCT2 < x, .default = "T")))

Data12 %>% .dumbCut(5000) %>% filter(Cell=="In") -> Data

#### plotting the nine cycles ####
Data %>% ggplot(aes(Clock, BCT2, col=Targ)) + 
  scale_color_manual(name="Target state", breaks=c("Chopper","On"), labels=c("Off","On"), values=c("black", "red")) + 
  geom_point() + 
  theme_minimal() + theme(legend.position = "top", legend.title=element_text()) +
  labs(x="Time (local)", y="I (a.u.)")

##################################
## shows that BCT2 = 0 is just stack overflow ##
# Data%>%filter(Run==971) %>% slice(1:60) -> xcheck
# ggplot(xcheck, aes(x=Clock)) + geom_point(aes(y=BCT1, col="BCT1")) + geom_point(aes(y=BCT2/10, col="BCT2/10")) +
#   geom_vline(xintercept=xcheck$UTS[4], col="red") +
#   theme_minimal() + ggtitle("Run 971; first minute") + labs(y="Current, ADC")
#################################

#### baseline analysis ####
base = Data %>% filter(FSgl=="F", BCT2 != 0) %>% .markOutliers("BCT1") %>% filter(FOut=="F")
mosaic::median(BCT1~Run, data=base) -> bmed
data.frame(Run = as.factor(names(bmed)), Med = bmed) -> bmed
## see if there's correlation between beam current and offset !!!!!!!!!!!!!!!!!!!!!##
##!!!!!!!!!!!!!!!!!!!!!! simultaneity => exogeneity !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!##
Data%>%filter(FSgl=="T") %>% daply("Run", function(r) median(r$BCT2[100:140])) -> I0
df = cbind(bmed, I0)

#### fit ####
Data %>% ddply("Run", function(r) {
  b = filter(r, FSgl=="F", BCT2 != 0)$BCT2 %>% median
  mutate(r, BCT2 = BCT2-b) %>% filter(FSgl=="T") -> .sub
  cat(r$Run[1]); cat("\n")
  
  fit.(.sub, .subset = c(45, 15)) -> .stats; cat(.stats);cat("\n")
  
  data.frame(
    Unit=1,
    Estimate = .stats[1], SE = .stats[2], Chi2 = .stats[3], I0 = .stats[4],
    Quality = .sub$Quality[5],
    eCool.I = .sub$eCool.I[5],
    Targ = .sub$Targ[5],
    Bunch = .sub$Bunch[5],
    Clock = .sub$UTS[1] %>% as.POSIXct(origin="1970-1-1")
  )
}, .parallel=TRUE) %>% ddply("Targ", .markOutliers) -> slopes

#### testing ####
.sub=c(45,15)
Data %>% filter(Run==969) %>% mutate(uts=UTS-UTS[1]) -> Run969
b = filter(Run969, FSgl=="F", BCT2 != 0)$BCT2 %>% median
Run969%>%mutate(BCT2=BCT2-b)%>%filter(FSgl=="T")->Run969
m = fit.(Run969, model=TRUE, .subset=.sub); .sub <- c(0, nrow(Run969))+c(1,-1)*.sub
slice(Run969, 20:nrow(Run969))%>%ggplot(aes(uts, log(BCT2))) + geom_line() + geom_smooth(method="lm", col="red", size=.75) + 
  theme_minimal() + labs(y="ln I", x="Time (seconds)") + 
  geom_vline(xintercept=Run969$uts[.sub], col="gray",linetype=2)

autoplot(m, which=1) + theme_minimal() + ggtitle("") + labs(y=expression(ln~I[i]~-~(hat(alpha)~+~hat(beta)~t)), x=expression(hat(alpha)~+~hat(beta)~t))
autoplot(pacf(m$residuals, lag=600, main="", plot=FALSE)) + theme_minimal() + labs(y="Partial ACF")


library(lmtest); library(strucchange)
f = log(BCT1)~uts
Run969 <- slice(Run969, .sub[1]:.sub[2])

harvtest(f, data=Run969)
raintest(f, data=Run969)
sctest(f, data=Run969, type="aveF")
sctest(f, data=Run969, type="ME")
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
  labs(x="Off-cycle", y=expression(hat(sigma)[0]~"(a.u.)"))

cs0mb%>%group_by(Soundness, Closeness) %>% 
  dplyr::summarise(
    NUM = n(),
    MEAN = mean(Estimate),
    W.MEAN = weighted.mean(Estimate),
    SD = sd(Estimate),
    SE = SD/sqrt(NUM)
  ) 

ggplot(cs0mb, aes(Estimate, col=Closeness)) +
  facet_grid(Soundness~., scale="free_y", space="free_y") +
  scale_color_manual(breaks=c("Close","Far"), values=c("black","red")) + 
  geom_density(trim=TRUE, kernel="rect",bw=10) + 
  geom_segment(aes(x=Estimate, xend=Estimate, y=0, yend=.005, col=Closeness)) +
  theme_minimal() + labs(x=expression(hat(sigma)[0]~"(a.u.)"), y="Kernel density") + theme(legend.position="top")

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
  labs(y=expression(hat(beta)))

slopes%>%group_by(Targ) %>% 
  dplyr::summarise(
    NUM = n(),
    MEAN = mean(Estimate),
    W.MEAN = weighted.mean(Estimate),
    SD = sd(Estimate),
    SE = SD/sqrt(NUM)
  )
