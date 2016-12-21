library(dplyr); library(plyr)
library(ggplot2); library(plotly); library(reshape2)
library(lattice)
library(mosaic)


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

ggplot(poldata, aes(Run, P0, col=Ring)) + 
  geom_pointrange(aes(ymin=P0-SEP0, ymax=P0+SEP0)) +
  ggtitle("Polarization") + 
  theme_bw() + theme(axis.text.x=element_text(angle=60))
  


read.table("./Stats/CR_data.txt", sep = "\t", quote="")[,1:6] -> crdata
names(crdata) <- c("Run","Ring", "CR0", "SECR0", "CRB", "SECRB")
crdata%>%mutate(CRLT = -1/CRB, SECRLT = CRLT^2*SECRB, rSECRLT = SECRLT/CRLT, rSECR0 = SECR0/CR0, Run = factor(Run), Ring=factor(Ring)) -> crdata

ggplot(filter(crdata, Run==7114, Ring%in%c(4:15)), aes(Ring, CR0)) + geom_point() + ggtitle("Cross-Ratio")

join(crdata, poldata) -> Data

Data %>%
  filter(!Ring %in% c(15:11)) %>% #Run%in%c(7107, 7114)) %>%
  ggplot(aes(CR0, P0, col=Ring)) + geom_point() + geom_text(aes(label=Ring, vjust=-.5))
