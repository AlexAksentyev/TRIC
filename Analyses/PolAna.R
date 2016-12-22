library(dplyr); library(plyr)
library(ggplot2); library(plotly); library(reshape2)
library(lattice)
library(mosaic)

source("EstiStats.R")

if(TRUE){
  read.table("./Stats/FitRes.txt", sep = "\t", quote="")[,1:5] -> crcomdata
  
  names(crcomdata) <- c("Run", "P0", "SEP0", "B", "SEB")
  
  mutate(crcomdata, LT = -1/B, SELT = LT^2*SEB, rSEP0 = SEP0/P0, rSEB = SEB/B, rSELT = SELT/LT) %>% 
    filter(!Run %in% c(7103, 7104, 7122)) -> crcomdata; crcomdata
  
  run = 7080:7102 %>% c(7105:7121)
  mqu15 = c(20, 30, 40, 50, 50, 60, 70, 80, 30, 40, 50, 60, 70, 80, 90, 20, 30, 40, 50, 60, 70, 90, 20) %>%
    c(rep(50, 17))
  mqu26 = c(10, 6, 3, 0, 10, 20, 30, 40, 10, 6, 3, 10, 20, 30, 40, 6, 3, 0, 20, 30, 40, 40, 10) %>%
    c(c(-5:5, -4.5:.5))
  data.frame(Run = run, MQU15 = mqu15, MQU26 = mqu26) -> MD
  
  dat = filter(crcomdata, Run %in% MD$Run)%>%cbind(MD%>%dplyr::select(-Run))%>%filter(Run!=7101)
  dat %>% mutate() %>% dplyr::arrange(abs(rSEP0)) %>% filter(abs(rSEP0)<1)
  
  acast(dat, MQU15~MQU26, value.var = "LT", fun.aggregate = mean) -> dat_xyz
  plot_ly(z=dat_xyz, type="surface")
  wireframe(LT ~ MQU15*MQU26, data=dat)
}



read.table("./Stats/Pol_data.txt", sep = "\t", quote="")[,1:4] -> poldata
names(poldata) <- c("Run","Ring", "P0", "SEP0")
poldata%>%mutate(rSEP0 = SEP0/P0, Run = factor(Run), Ring=factor(Ring)) -> poldata

poldata %>% .markOutliers("SEP0") %>% filter(FOut == "F") %>%
ggplot(aes(Run, P0, col=Ring)) + 
  geom_pointrange(aes(ymin=P0-SEP0, ymax=P0+SEP0)) +
  ggtitle("Polarization") + 
  theme_bw() + theme(axis.text.x=element_text(angle=90))
  


read.table("./Stats/CR_data.txt", sep = "\t", quote="")[,1:6] -> crdata
names(crdata) <- c("Run","Ring", "CR0", "SECR0", "CRB", "SECRB")
crdata%>%mutate(CRLT = -1/CRB, SECRLT = CRLT^2*SECRB, rSECRLT = SECRLT/CRLT, rSECR0 = SECR0/CR0, Run = factor(Run), Ring=factor(Ring)) -> crdata

crdata %>% 
  filter(Run %in% 7105:7121, Ring%in%9:12) %>%
  .markOutliers("CRLT") %>% filter(FOut == "F") %>%
  ggplot(aes(Run,CRLT, col=Ring)) + geom_pointrange(aes(ymin=CRLT-SECRLT, ymax=CRLT+SECRLT)) + ggtitle("Cross-Ratio") + scale_y_log10() +
  theme_bw() + theme(axis.text.x=element_text(angle=90)) +
  facet_grid(Ring~.)

join(poldata, MD) -> Data

Data%>%filter(MQU26<1, Ring%in%c(9:12))%>%ggplot(aes(MQU26/10, P0, col=Ring)) + geom_pointrange(aes(ymin=P0-SEP0, ymax=P0+SEP0)) + theme_bw() +
  facet_grid(Ring~.)

Data%>%filter(Ring%in%9:12, !is.na(MQU26))->dat
acast(dat, MQU15~MQU26, value.var = "P0", fun.aggregate = mean) -> dat_xyz
# plot_ly(z=dat_xyz, type="surface")
levelplot(P0 ~ MQU15*MQU26|Ring, data=dat, pretty=TRUE, drape=TRUE)
# persp(as.numeric(rownames(dat_xyz)), as.numeric(colnames(dat_xyz)), dat_xyz, phi=30, theta=-120)
