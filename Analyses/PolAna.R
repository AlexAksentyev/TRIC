library(dplyr); library(plyr)
library(ggplot2); library(reshape2)
library(lattice)
library(mosaic)

# right now I don't have the data to do the qu15/qu26 analysis
# only polarization life-time analysis

source("EstiStats.R")

getPolData <- function(){
  read.table("./Stats/Pol_data.txt", sep = "\t", quote="")[,1:4] -> poldata
  names(poldata) <- c("Run","Ring", "P0", "SEP0")
  poldata%>%mutate(rSEP0 = SEP0/P0, Run = factor(Run), Ring=factor(Ring))
}
getCRData <- function(){
  read.table("./Stats/CR_data.txt", sep = "\t", quote="")[,1:6] -> crdata
  names(crdata) <- c("Run","Ring", "CR0", "SECR0", "CRB", "SECRB")
  crdata%>%mutate(CRLT = -1/CRB, SECRLT = CRLT^2*SECRB, rSECRLT = SECRLT/CRLT, rSECR0 = SECR0/CR0, Run = factor(Run), Ring=factor(Ring))
}

#### metadata prep ####
run = 7080:7102 %>% c(7105:7121)
mqu15 = c(20, 30, 40, 50, 50, 60, 70, 80, 30, 40, 50, 60, 70, 80, 90, 20, 30, 40, 50, 60, 70, 90, 20) %>%
  c(rep(50, 17))
mqu26 = c(10, 6, 3, 0, 10, 20, 30, 40, 10, 6, 3, 10, 20, 30, 40, 6, 3, 0, 20, 30, 40, 40, 10) %>%
  c(c(-5:5, -4.5:.5))
data.frame(Run = run, MQU15 = mqu15, MQU26 = mqu26) -> MD

#### Polarizarion data ####
poldata <- getPolData()

if(FALSE){
  poldata %>% .markOutliers("SEP0") %>% filter(FOut == "F") %>%
    ggplot(aes(Run, P0, col=Ring)) + 
    geom_pointrange(aes(ymin=P0-SEP0, ymax=P0+SEP0)) +
    ggtitle("Polarization") + 
    theme_bw() + theme(axis.text.x=element_text(angle=90))
  
  join(poldata, MD) -> poldata
  
  poldata%>%filter(MQU26<1, Ring%in%c(9:12))%>%ggplot(aes(MQU26/10, P0, col=Ring)) + geom_pointrange(aes(ymin=P0-SEP0, ymax=P0+SEP0)) + theme_bw() +
    facet_grid(Ring~.)
  
  poldata%>%filter(Ring%in%9:12, !is.na(MQU26))->dat
  levelplot(P0 ~ MQU15*MQU26|Ring, data=dat, pretty=TRUE, drape=TRUE)
}


#### Cross ratio data ####
crdata <- getCRData()

if(FALSE){
  crdata %>% filter(Run %in% 7105:7121, Ring%in%9:12) %>%
    .markOutliers("CRLT") %>% filter(FOut == "F") %>%
    ggplot(aes(Run,CRLT, col=Ring)) + geom_pointrange(aes(ymin=CRLT-SECRLT, ymax=CRLT+SECRLT)) + ggtitle("Cross-Ratio") + scale_y_log10() +
    theme_bw() + theme(axis.text.x=element_text(angle=90)) +
    facet_grid(Ring~.)
}


#### lifetime analysis ####
run=7129:7136
btime = seq(3, 4.75, .25)
utime = rep(2, 8)

MD = data.frame(Run=run, BTime=btime, UTime=utime)
join(crdata,MD)%>%filter(!is.na(BTime)) -> x

mutate(x, G = derivedFactor(
  "D" = BTime <= 4,
  .default = "U"
)) -> x

library(lme4)
lmer(CR0~BTime + (BTime|Ring) + (BTime|G), data=x) -> m3
summary(m3)
coef(summary(m3))[2,1:2] -> b
lt = -1/b[1]; selt = lt^2*b[2]

mutate(x, Run=as.numeric(Run), Ring=as.numeric(Ring)) ->x
kmfit <- kmeans(x[,1:3], 8)
mutate(x, KClus = as.factor(kmfit$cluster))->x
library(cluster)
clusplot(dplyr::select(x,Run,Ring, CR0), kmfit$cluster, color=TRUE, shade=TRUE, lines=0, stan=TRUE)

prcomp(dplyr::select(x,Run,Ring, CR0)) -> xpca

xyplot(CR0~BTime|Ring, data=filter(x, Ring%in%9:15))
xyplot(CR0~Ring|Run, data=mutate(x, BTime = as.factor(BTime)))#, groups=KClus)
