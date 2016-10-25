source("DataAna.R")


library(dplyr); library(plyr); library(splitstackshape)
library(parallel); library(doParallel)


source("Parameters.R")
thick = 6.92e13

#### 2016 data ####
get2016Data() -> Data16
slopes16 <- Data16$Slopes; Data16 <- Data16$Data

#### 2012 data ####
get2012Data() -> Data12
slopes12 <- Data12$Slopes; Data12 <- Data12$Data

#### merging slopes ####
slopes12 %>% mutate(Year = 2012) -> slopes12
slopes16 %>% filter(B.Spin == "N") %>%  
  mutate(Year = 2016, Run = NA, Quality = "ok", Bunch = "RF", eCool.I = 45, Size = "Small", Targ = derivedFactor(
    "On" = FABS == "T", "Chopper" = FABS == "F"
  )) %>% dplyr::select(-T.Spin, -B.Spin, -FABS) %>% 
  rbind(slopes12) -> slopes

slopes %>% ddply(.(Bunch, Size, eCool.I), function(.sub) mutate(.sub, Weight = SE^(-2)/max(SE^(-2)))) %>%
  mutate(
    I0 = cut(I0, breaks=4, labels=letters[1:4]), 
    Year = factor(Year)
  ) -> slopes

rm(slopes12, slopes16)

if(FALSE){
  #### Slope-CS0 Analysis ##
  ## big cycles ##
  library(cowplot)
  slopes12.big = filter(slopes, Year==2012, Size=="Big", FOut=="F") 
  slopes12.big %>% dlply("Targ") -> l.big; 
  
  l.big[c("Chopper","On")] %>% dbeta.(.fields=c("Run","Clock")) %>% 
    dplyr::select(-Unit.On,-Unit.Chopper) %>% 
    mutate(Estimate = Estimate/nu/thick*1e27, SE = SE/nu/thick*1e27, Weight = SE^(-2)/max(SE^(-2))) -> cs0mb.big
  
  l.big %>% ldply(.id="Targ") %>% mutate(Group = derivedFactor("5--6th" = Targ=="Valve", .default="7th")) %>%
    ggplot(aes(y=Estimate, x=Clock, col=Targ, shape=I0)) + 
    geom_point() + geom_errorbar(aes(ymin=Estimate-SE, ymax=Estimate+SE), width=0) +
    facet_grid(~ Group, scales = "free_x") +
    geom_smooth(method=lm) + 
    theme_minimal() + ggtitle("slope; 2012-10-07") -> slp.big.plot
  
  
  cs0mb.big %>% ggplot(aes(y=Estimate, x=Clock.On, col=Run.Chopper)) + 
    geom_point() + #geom_errorbar(aes(ymin=Estimate-SE, ymax=Estimate+SE), width=.1) +
    # geom_smooth(method=lm) +
    theme_minimal() + ggtitle("cs0; 2012-10-07") -> cs0mb.big.plot
  
  plot_grid(slp.big.plot, cs0mb.big.plot, nrow=2)
  
  ## small cycles ##
  slopes12.small = filter(slopes, Year==2012, Size=="Small", FOut=="F", Quality=="ok")
  slopes12.small %>% filter(Bunch=="Off", FOut == "F", Run == 935) %>%
    mutate(uts = (as.numeric(Clock)-as.numeric(Clock)[1])/3600) %>%  dlply("Targ") -> l.small
  l.small[c("Valve","On")] %>% dbeta.(.fields=c("Run","Clock","uts")) %>% 
    mutate(
      Estimate = Estimate/nu/thick*1e27, 
      SE = SE/nu/thick*1e27, 
      Weight=SE^(-2)/max(SE^(-2)),
      HDiff = uts.Valve-uts.On
    ) -> cs0mb.small
  
  l.small %>% ldply(.id="Targ") %>% 
    ggplot(aes(y=Estimate, x=Clock, shape=I0)) + 
    geom_point(aes(col=eCool.I)) + geom_errorbar(aes(ymin=Estimate-SE, ymax=Estimate+SE, col=eCool.I), width=0) +
    geom_smooth(aes(linetype=Targ), method=lm) +
    theme_minimal() + ggtitle("slope; 2012-10-02") -> slp.small.plot
  
  cs0mb.small  %>% ggplot(aes(y=Estimate, x=HDiff, col=as.factor(Clock.On))) + 
    geom_point() + geom_errorbar(aes(ymin=Estimate-SE, ymax=Estimate+SE), width=0) + 
    geom_smooth(method=lm) + 
    theme_minimal() + 
    labs(x="Valve-On time difference") + ggtitle("cs0; 2012-10-02") -> cs0mb.small.plot
  
  plot_grid(slp.small.plot, cs0mb.small.plot, nrow=2)
}

#### CLUSTERING ####
library(mclust)
slopes %>%filter(Year == 2012, FOut=="F", Quality=="ok") %>% 
  ggplot(aes(x=Estimate,y=SE, col=Targ, shape=Size)) + geom_point(size=2) + scale_y_log10()

slopes %>% filter(Year == 2012, Size=="Big", Targ != "Valve") %>% mutate(Weight = cut(Weight,breaks=4)) %>%
  ggplot(aes(x=Run, y=Estimate, shape=Targ, col=Weight)) + 
  geom_point(size=2) + geom_errorbar(aes(ymin=Estimate-SE,ymax=Estimate+SE)) +
  geom_smooth() + 
  ggtitle("slope; 2012")

slopes %>% filter(FOut=="F",Quality=="ok", Size=="Small", Year==2016) %>% 
  ddply(.(Bunch), function(.sub) .sub %>% mutate(Weight = SE^(-2)/max(SE^(-2)), RSE = abs(SE/Estimate))) -> slp

slp %>% dplyr::select(Estimate, Weight) %>% scale(center=FALSE) %>% Mclust() %>% plot(what="class")

slp %>% mutate(Estimate = scale(Estimate, center=FALSE), Weight = scale(Weight, center=FALSE)) %>%
  ggplot(aes(x=Estimate, y=Weight, col=Targ, shape=I0)) + 
  geom_point() + 
  # scale_x_log10() +
  theme_minimal() + theme(legend.position="top")

#### ANOMALIES ####
.anomalies <- function(.data, neigh = 5, categ = 3){
  vars = names(.data)
  .data %>% lof(neigh) %>% cut(breaks=categ, labels=1:categ) -> lof.w
  
  .data %>% ggplot(aes_string(x=vars[1], y=vars[2], col="lof.w")) + geom_point()
}

slopes16 %>% mutate(Unit = rep(1:(length(Unit)/2), each=2)) -> slopes16

slopes16 %>% ggplot(aes(Clock, Estimate, col=FABS)) + 
  geom_pointrange(aes(ymin=Estimate-SE, ymax=Estimate+SE)) + 
  geom_smooth(aes(weight = slopes16$SE^(-2)), span=.5) -> slopes16.plot

slopes16 %>% ggplot(aes(Clock, I0, col=FABS)) + 
  geom_point() + 
  geom_smooth(aes(weight = slopes16$SE^(-2)), span=.5) -> I016.plot

plot_grid(slopes16.plot, I016.plot, nrow=2)


slopes16 %>% filter(FABS=="T") %>% mutate(Weight = SE^(-2)%>%scale(center=FALSE)) %>% 
  dplyr::select(I0, Estimate) %>% .anomalies(5,5)
