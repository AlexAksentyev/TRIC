library(dplyr)
source("Utilities.R")

Clean.Spin <- function(.Data){
  sb.i = which(substr(names(.Data),1,2) == "SB")
  
  sapply(sb.i, function(i, X) cut(X[,i], breaks = 2, labels = c("0", "1")), .Data) -> .Data[,sb.i]
  
  .Data$B.Spin = strtoi(do.call(paste0, .Data[,sb.i[order(sb.i, decreasing = TRUE)]]), base = 2)
  
  return(subset(.Data, select = -sb.i))
}



loadData(from = "2016-6-30 22:50:00", to = "2016-7-1 04:00:00") %>% mutate(UTS = UTS - UTS[1]) -> Data



with(Data, which(diff(UTS) == 2))->i.skip; ind = seq(1, dim(Data)[1], length.out = 6)

Data %>% mutate(Skip = "F")-> Data; Data[i.skip, ]$Skip = "T"

Clean.Spin(Data) -> Data

# rep(i.skip,3)+rep(c(-1,0,1),each = length(i.skip))->i; View( Data[i[order(i)], ])

Data[i.skip,] -> df; d = diff(df$UTS)
par(mfrow=c(2,1))
t0 = 0:(dim(Data)[1]-1); with(Data, plot((UTS-t0)~UTS, type="l", main = "Overtime", xlab="Run Time", xaxt="n")); #abline(v=Data$UTS[i.skip], col="red")
with(Data[ind,], axis(1, at = UTS, labels = deciTots(UTS/3600, "h:m")))
plot(density(d, kernel = "rect", bw=2), main = "EPDF of Smooth Run Period")
legend("topright", ncol=2, legend = c(c("Mean","SD"), round(c(mean(d), sd(d)),1)))

