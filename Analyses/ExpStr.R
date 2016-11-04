# I'll try to reconstruct the experimental model from data

source("DataAna.R")
library(ggfortify)

get2012Data() -> Data
slopes <- Data$Slopes; Data <- Data$Data

tRun = filter(Data, Run==969) %>% mutate(uts=UTS-UTS[1])
tRun%>%ggplot(aes(Clock, log(BCT2))) + geom_line()

n=nrow(tRun); lm(log(BCT2)~uts, data=tRun) -> tm
x = tm$model$uts; y = tm$coef[1] - tm$coef[2]*x

enorm = rnorm(n,sd=summary(tm)$sigma) ##random normal error model
earima = arima.sim(list(1,0,0),n, sd=summary(tm)$sigma) ##random arima error model

df = data.frame(Y = log(tRun$BCT2), X = y+enorm); lm(Y~X, data=df) -> icgm
df %>% ggplot(aes(X,Y)) + geom_line() + geom_abline(slope=1,intercept=0, col="red") + 
  geom_smooth(method="loess", se=FALSE, col="red", linetype=2, size=.75) +
  labs(x="model",y="real") + theme_minimal()


confint(icgm,parm="X", level=.99)


loess(Y~X, data=df) -> icgml
plot(icgml$fitted~icgml$x, type="l")
abline(icgm, col="red")
plot(icgml$fitted~icgm$fitted.values, type="l"); abline(0,1)

optimize(
  f=function(lam) (lm(BoxCox(BCT2, lam)~uts, data=tRun)$residuals %>% moments::jarque.test())$statistic,
  interval= c(-1,1)
)$minimum -> lamb
## apparently, the Box-Cox transformation with lambda=-.19 is better than with lambda=0 (ln transform)
## in terms of the model residuals

lm(BoxCox(BCT2, lambda=lamb)~uts, data=tRun) -> tm_19
autoplot(tm,1); autoplot(tm_19,1)


get2016Data() -> Data16
slopes16 <- Data16$Slopes; Data16 <- Data16$Data

#### PANEL ANALYSIS ####
# Suppose the offset does indeed depend on beam current. We have multiple cycles with varying beam intensities;
# presumably, the offsets for each are different. If we panel-analyze them, we should see a difference in 

##### MACHINE LEARNING #####  
run = filter(Data16, Unit == 1, FABS=="T", FSgl=="T", eCool=="Acc") %>% mutate(uts=UTS-UTS[1])
f = log(BCT2) ~ uts
run%>%with(plot(log(BCT2)~uts))
nnet(f, data=run, size=3) -> x
lm(f, data=run) -> y
with(y, points(fitted.values~run$uts, col="red"))
with(x, points(fitted.values~run$uts, col="blue"))
