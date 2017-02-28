library(dplyr); library(plyr)
library(ggplot2); library(reshape2)
library(ggfortify)
library(lattice)
library(mosaic)
library(cluster)

source("EstiStats.R")

getPolData <- function(){
  read.table("./Stats/Pol_data.txt", sep = "\t", quote="")[,1:6] -> poldata
  names(poldata) <- c("Run","Ring", "P0", "SEP0", "Chi2", "NDF")
  poldata%>%mutate(Run = factor(Run))
}
getCRData <- function(){
  read.table("./Stats/CR_data.txt", sep = "\t", quote="")[,1:8] -> crdata
  names(crdata) <- c("Run","Ring", "CR0", "SECR0", "CRB", "SECRB", "Chi2", "NDF")
  crdata%>%mutate(CRLT = -1/CRB, SECRLT = CRLT^2*SECRB, 
                  rSECRLT = SECRLT/CRLT, rSECR0 = SECR0/CR0,
                  Chi2red = Chi2/NDF, Run = factor(Run))
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

## QU15, QU26 analysis
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

## QU15, QU26 analysis
if(FALSE){
  crdata %>% filter(Run %in% 7105:7121, Ring%in%9:12) %>%
    .markOutliers("CRLT") %>% filter(FOut == "F") %>%
    ggplot(aes(Run,CRLT, col=Ring)) + geom_pointrange(aes(ymin=CRLT-SECRLT, ymax=CRLT+SECRLT)) + ggtitle("Cross-Ratio") + scale_y_log10() +
    theme_bw() + theme(axis.text.x=element_text(angle=90)) +
    facet_grid(Ring~.)
}


#### lifetime analysis ####
## runs 7129--36 are not useful for that!
if(FALSE){
  run=7129:7136
  btime = seq(3, 4.75, .25)
  utime = rep(2, 8)
  
  MD = data.frame(Run=run, BTime=btime, UTime=utime)
  join(crdata,MD)%>%filter(!is.na(BTime)) -> x
  
  mutate(x, G = derivedFactor(
    "D" = BTime <= 4,
    .default = "U"
  )) %>% mutate(G=as.numeric(G)) -> x
  
  cmts = c("CR0","CRB","SECR0","SECRB")
  ddply(x, "Run", function(s) 
    mutate(s, 
           KClus = as.factor(kmeans(s[,cmts],3)$cluster)
    )
  ) -> x
  autoplot(kmfit, data=x, frame=TRUE) + theme_bw()
  
  prcomp(x[,cmts]) -> xpca; summary(xpca)
  
  xyplot(CRLT~BTime|factor(Ring, labels=9:15), data=filter(x), groups = KClus)
  xyplot(CRLT~Ring|Run, data=mutate(x, BTime = as.factor(BTime)), groups=KClus)
  
}

data = join(crdata,poldata, c("Run", "Ring"))
melt(data, id.vars = c("Run","Ring")) -> mdata
ggplot(filter(mdata,variable%in%c("CRLT","PLT"))) + 
  geom_histogram(aes(value), fill="white",col="black") + 
  facet_grid(.~variable, scales = "free_x") + theme_bw()

wch = "P0"
filter(data,!Ring%in%c(10,15)) %>% 
  ggplot(aes_string("Ring", wch)) + 
  geom_pointrange(aes_string(ymin=paste0(wch,"-",paste0("SE",wch)), 
                             ymax=paste0(wch,"+",paste0("SE",wch))
                  )
  ) + 
  theme_bw()

mutate(data, Ring=as.factor(Ring)) -> data
filter(data, !Ring%in%c(10,15)) %>%
  ggplot(aes(CR0,CRLT, col=Ring)) + geom_point() +
  theme_bw() + theme(legend.position="") +
  labs(x=expression(CR[0]), y=expression(hat(tau)[CR]))+
  geom_errorbar(aes(ymin=CRLT-SECRLT,ymax=CRLT+SECRLT)) +
  geom_errorbarh(aes(xmin=CR0-SECR0,xmax=CR0+SECR0)) -> corplot

filter(data, !Ring%in%c(10,15)) %>% ggplot(aes(CRLT)) + labs(x=expression(hat(tau)[CR])) +
  geom_density(kernel="gaus") +geom_rug(aes(col=Ring)) + theme_bw() +theme(legend.position="top") -> dplot

library(cowplot)
plot_grid(corplot,dplot,nrow = 2)
