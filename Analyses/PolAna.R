library(dplyr); library(plyr)
library(ggplot2); library(plotly); library(reshape2)
library(lattice)

read.table("./Stats/FitRes.txt", sep = "\t", quote="")[,1:5] -> poldata

names(poldata) <- c("Run", "P0", "SEP0", "B", "SEB")

mutate(poldata, LT = -1/B, SELT = LT^2*SEB) %>% 
  filter(!Run %in% c(7103, 7104, 7122)) -> poldata

run = 7080:7102 %>% c(7105:7121)
mqu15 = c(20, 30, 40, 50, 50, 60, 70, 80, 30, 40, 50, 60, 70, 80, 90, 20, 30, 40, 50, 60, 70, 90, 20) %>%
  c(rep(50, 17))
mqu26 = c(10, 6, 3, 0, 10, 20, 30, 40, 10, 6, 3, 10, 20, 30, 40, 6, 3, 0, 20, 30, 40, 40, 10) %>%
  c(c(-5:5, -4.5:.5))
data.frame(Run = run, MQU15 = mqu15, MQU26 = mqu26) -> MD


dat = filter(poldata, Run %in% MD$Run)%>%cbind(MD%>%dplyr::select(-Run))%>%filter(Run!=7101, LT>0)
dat%>%mutate(RSELT = SELT/LT * 100) %>% dplyr::arrange(RSELT)

acast(dat, MQU15~MQU26, value.var = "LT", fun.aggregate = mean) -> dat_xyz
plot_ly(z=dat_xyz, type="surface")
wireframe(LT ~ MQU15*MQU26, data=dat)
