source("EstiStats.R")
source("DataPrep.R")
source("DataAna.R")

library(Rcpp)
library(dplyr); library(plyr)
library(Rlof)

sourceCpp("./src/gas_calk.cc")
source("./Data/Pressures.R")

.calc.conc <- function(Pressures, Temp = 27) density_mb2Natccm(mean(Pressures), 273.15 + Temp)

.calc.conc2 <- function(Pressure, Temp = 27, Vmol = 24){mean(Pressure)/1013.25 * 273.15/(273.15 + Temp) * 6.02e23/Vmol * 1e-3}

### constants ####
MHat = 1.008
MDat = 2.014

FTL = 10; FTD = .96;
CellLl = 20; CellLr = 20; CellL = CellLl+CellLr
CellD = .96;
XtrtL = 38; XtrtD = .96

### variables ####
TN = 273.15 + 27
Fabs = 3.3e16
MGAt = MDat; #Molar mass [g/mol] of used Gas

#### computation ####

Ca = Conduct_SCT2(FTL,FTD,TN,MGAt);
Cb = Conduct_SCT2(CellLl,CellD,TN,MGAt);
Cc = Conduct_SCT2(CellLr,CellD,TN,MGAt);
Cd = Conduct_SCT2(XtrtL,XtrtD,TN,MGAt)
Ctot = Ca+Cb+Cc+Cd

th.cell = CellThickness(Fabs,Ctot,CellL);


if(FALSE){
  
  # list("12" = Pressures12, "16" = Pressures16) %>% llply(function(e) e[which(names(e) == "PAX")]) %>%
  #   laply(function(e) .calc.conc(e) * 1e2*.5e2*.5e2/(40*pi*10^2/4)) %>%
  #   ( function(x) print( c(formatC(x, format="e"), round(x[2]/x[1], 1)) ) )
  

  Vvc = 1e2*.5e2*.5e2; Vcell = CellL*pi*CellD^2/4; Lcosy=18400
  Pressures16[which(names(Pressures16)=="PAX")] -> Ppax
    .calc.conc(Ppax)*Vvc/Vcell*CellL -> th.cell.pres
  Pressures16[which(names(Pressures16)!="PAX")] -> Pcosy
    .calc.conc(Pressures16)*Lcosy -> th.cosy
    
  cat(paste0("Cell thickness ", formatC(th.cell,2,format="e"), ", density ", formatC(th.cell/CellL,2,format="e"), "\n"))
  cat(paste0("COSY thickness ", formatC(th.cosy,2,format="e"), ", density ", formatC(th.cosy/Lcosy,2,format="e"), "\n"))
  cat(paste("Cell thickness/Cosy thickness", round(th.cell/th.cosy, 2), "\n"))
  cat(paste("Cell thickness (calk) / Cell thickness (mine)", round(th.cell/th.cell.pres, 2), "\n"))
}


if(TRUE){
  #### ABS data ####
  list.files(path = "./Data", pattern = "ABS2") %>% ldply(getABS) %>%
    mutate(State = cut(Flux, breaks = 2, labels = c("F", "T")), uts = UTS - UTS[1]) -> ABS
  
  ABS %>% filter(State=="T", Flux > 8e4) %>% mutate(Clock = as.POSIXct(UTS, origin="1970-1-1")) -> ABS1
  lm(Flux~uts, data = ABS1) -> m
  .collect.stats(m) -> .stats; .stats <- .stats[1:3]; .stats[1:2] <- .stats[1:2]#*3600
  
  ABS1 %>% ggplot(aes(x=Clock, y=Flux)) + geom_point() + geom_smooth(method=lm) +
    ggtitle("ABS Flux 2016, July 1--2") + 
    annotate(
      "text", 
      label = paste("Slope =", paste0(round(.stats[1:2],2), collapse ="+-")), 
      x=mean(ABS1$Clock), y=min(ABS1$Flux), col="red"
    ) + theme_minimal() -> flux.plot
  
  ABS1 %>% slice(30:130) %>% ggplot(aes(Clock, Flux)) + geom_line() + theme_minimal() + ggtitle("ABS 100 secs")
  
  ## comparing with slopes
  get2016Data() -> Data16
  Data16$Slopes %>% mutate(
    I0 = cut(I0, breaks=4, labels=letters[1:4]), 
    Targ = derivedFactor(
      "On" = FABS == "T", "Off" = FABS == "F"
    )
  ) %>% dplyr::select(-FABS) -> slopes16
  Data16$Data %>% mutate(
    Unit = (1:length(rle(Unit)$values)) %>% rep(times=rle(Unit)$lengths),
    I0 = cut(BCT2, breaks=4, labels=letters[1:4]),
    Datetime = as.POSIXct(paste(Date, Time), format = "%Y-%m-%d %H:%M:%S")
  ) -> Data16
  
  slopes16 %>% filter(FOut == "F") %>% 
    mutate(Estimate = -1/Estimate, SE = SE*Estimate^2) -> tau
  
  tau %>% ggplot(aes(x=Clock, y=Estimate, col=B.Spin)) + 
    geom_point() + geom_errorbar(aes(ymin=Estimate-SE,ymax=Estimate+SE)) + 
    geom_smooth(method="loess", se=FALSE, mapping=aes(weight = tau$SE^(-2))) + 
    facet_grid(Targ~.) + 
    ggtitle("tau; 2016, July 1") + theme_minimal() -> tau.plot
  
  library(cowplot)
  
  plot_grid(flux.plot, tau.plot, nrow=2)
  
  #### looking for the extra, non-linear scattering ####
  
  Data16 %>% ddply(.(Unit, FABS), function(u){
    u %>% mutate(DlnI = c(NA, diff(log(BCT2)))) -> u
    (lm(DlnI~1)%>% summary %>% coef)[1:2] -> b
    u %>% mutate(Xtra = DlnI-b[1], XSE = b[2]) -> u
  }) -> Data16
  
  Data16%>% filter(!is.na(DlnI), FABS=="T", T.Spin=="1") %>% 
    ggplot(aes(x=UTS, y=Xtra)) + 
    geom_point() #+ 
  # facet_grid(I0~T.Spin, scales="free_x", space="free_x")
  
  Data16 %>% filter(Unit==5, FABS=="T") -> u5
  u5 %>% ggplot(aes(x=1:350, y=Xtra)) + geom_line()
  
  u5 %>% mutate(Frame = cut(UTS, breaks=24, labels=1:24)) -> u5
  u5 %>% ggplot(aes(x=Frame, y=Xtra)) + geom_violin(scale="count", draw_quantiles = c(0.5), na.rm=TRUE) 
  Xmed = median(Xtra~Frame, data=u5, na.rm = TRUE)
  
  
  
  Data16 %>% ddply(.(Unit, FABS), function(u){
    mutate(u, Frame = cut(UTS, breaks=24, labels=1:24))[-1,] %>%
      group_by(Frame) %>% 
      dplyr::summarise(Xmed = median(Xtra, na.rm=TRUE), Num = n())
  }) -> Xtra  
}


