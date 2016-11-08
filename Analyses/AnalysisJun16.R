# library(lmtest) #tests linear regression models
library(nleqslv) #non-linear system of equations solver
library(dplyr); library(plyr)
library(parallel); library(doParallel)

source("Utilities.R")
source("EstiStats.R")
source("./Tests/DataPrep.R")
source("Parameters.R")
source("DataAna.R")


registerDoParallel(detectCores())

get2016Data() -> Data16
slopes.uc <- Data16$Slopes; Data <- Data16$Data

###########################################################################################

a = ExpParameters["Rev.f"]*ExpParameters["Pol.targ"]*ExpParameters["Targ.len"]*ExpParameters["Targ.dens.on"]; names(a)<-NULL
dP = -diff(Pb)
beta <- filter(slopes, FABS == "T"); boff <- filter(slopes, FABS == "F", B.Spin == "N")
bU = filter(beta, B.Spin == "U"); bD = filter(beta, B.Spin == "D"); bN = filter(beta, B.Spin == "N");
bUD = filter(beta, B.Spin != "N")

d1 = dbeta.(list("D" = bD, "U" = bU)) %>% mutate(Estimate = Estimate/dP/a, SE = SE/dP/a)
d2 = dbeta.(list("N" = bN, "UD" = bUD), .fields = "B.Spin") %>% 
  mutate(Estimate = Estimate/Pb[B.Spin.UD]/a, SE = abs(SE/Pb[B.Spin.UD]/a))
d3 = dbeta.(list("Off" = boff, "N" = bN)) %>% mutate(Estimate = Estimate/nu, SE = SE/nu)
d4 = dbeta.(list("Off" = boff, "UD" = bUD), .fields = "B.Spin") %>% 
  mutate(Estimate = Estimate/Pb[B.Spin.UD]/a, SE = abs(SE/Pb[B.Spin.UD]/a))

###################### Way 2: CS0 FROM DATA ###########################################################
## just estimating cs0 first assuming knowledge of Dthickness (1.1e14) ##
cs0 = d3 %>% mutate(
  Estimate = Estimate/ExpParameters["Targ.len"]/(ExpParameters["Targ.dens.on"]- ExpParameters["Targ.dens.off"]),
  SE = SE/ExpParameters["Targ.len"]/(ExpParameters["Targ.dens.on"]- ExpParameters["Targ.dens.off"]),
  N.Cat = .smartCut(d3, which="Estimate.N"),
  Off.Cat = derivedFactor("Up" = Estimate.Off > -1.6e-4, .default = "Norm")
) %>% mutate(N.Cat = derivedFactor("Low" = N.Cat == "F", "Mid" = N.Cat == "T" & Unit.N != 28, .default = "Up"))

cs0 %>% mutate(Estimate = Estimate*1e27, SE = SE*1e27) %>% group_by(N.Cat, Off.Cat) %>% 
  dplyr::summarise(NUM = n(), 
                   MEAN = mean(Estimate), 
                   W.MEAN = weighted.mean(Estimate, SE^(-2)), 
                   SD = sd(Estimate),
                   SE = SD/sqrt(NUM)) -> cs0.stats; cs0.stats

write.table(x=cs0.stats, file='summary_table_2016.ods',sep='\t',row.names=F)


cs0.est = filter(cs0.stats, N.Cat == "Mid", Off.Cat == "Norm")
cs0mb = filter(cs0, N.Cat == "Mid", Off.Cat == "Norm") %>% mutate(Estimate = Estimate*1e27, SE = SE*1e27)

Ayy <- d1 %>% mutate(
  Estimate = Estimate/mean(cs0.est$WMN), SE = SE/mean(cs0.est$WMN),
  U.Cat = factor(.smartCut(d1, "Estimate.U"), labels = c("Low", "High"))
)

Ayy %>% group_by(U.Cat) %>%
  dplyr::summarise(NUM = n(), 
                   RMN = mean(Estimate), WMN = weighted.mean(Estimate, SE^(-2)), 
                   SE = sd(Estimate)/sqrt(n())) -> Ayy.stats; Ayy.stats


## pics
par(mfrow=c(3,1))
filter(Ayy, U.Cat == "Low") %>% .QQNorm("Ayy")
filter(Ayy, U.Cat == "High")%>% .QQNorm("Ayy")
Ayy %>% .QQNorm("Ayy")

cs0mb %>% .QQNorm(plot.it = TRUE)

#####################################################################
library(beanplot)

tau = transmute(slopes, Estimate = -1/Estimate, SE = SE*Estimate^2, B.Spin, T.Spin, FABS, Unit)
par(mfrow=c(1,1))
beanplot(
  Estimate~FABS*B.Spin, data = tau, kernel = "rect", 
  side = "both", method = "stack",
  col = list("lightblue","red"),
  what = c(0,1,1,1),
  overallline = "median",
  main = "Bean plot of lifetimes",
  xlab = "Beam Spin", ylab = bquote(hat(tau)~"sec"),
  names = c("Up","Down","Null")
)
legend(
  "topleft", fill=c("lightblue","red"),
  legend = paste("Target", c("off", "on"))
)

beans <- function(x, whch) bwplot(
  Estimate~B.Spin|Grp*FABS, data = x,
  panel=whch, kernel="rect", adjust=2
)

if(FALSE){ ## what if Spin Up--Down classification changed?
  g = 2; tau %>% arrange(Unit) %>% mutate(Grp = rep(1:g, each=72/g) %>% factor) %>% beans(panel.bwplot)
  slopes %>% arrange(Unit) %>% mutate(Grp = rep(1:g, each=72/g) %>% factor) %>% 
    mutate(
      B.Spin = derivedFactor(
        "U" = Grp==2 & B.Spin=="D",
        "D" = Grp==2 & B.Spin=="U",
        "U" = Grp==1 & B.Spin=="U",
        "D" = Grp==1 & B.Spin=="D",
        "N" = B.Spin=="N",
        .method = "last"
      )
    ) -> slopes.1
  slopes.1 %>% transmute(Estimate = -1/Estimate, SE = SE*Estimate^2, B.Spin, T.Spin, FABS, Unit) -> tau.inv
  
  beanplot(
    Estimate~FABS*B.Spin, data = tau.inv, kernel = "rect", 
    side = "both", method = "stack",
    col = list("lightblue","red"),
    what = c(0,1,1,1),
    overallline = "median",
    main = "Bean plot of lifetimes",
    xlab = "Beam Spin", ylab = bquote(hat(tau)~"sec"),
    names = c("Up","Down","Null")
  )
  legend(
    "topleft", fill=c("lightblue","red"),
    legend = paste("Target", c("off", "on"))
  )
  

}

# cat(WMN(d3)/117e-27) #@Berkley, 98MeV (440Mev/c)
# cat(WMN(d3)/50e-27) #@Particle Data Group, pd total collision cs, 135MeV (521Mev/c)
# cat(WMN(d3)/80e-27)#@PDG, pd, 98MeV