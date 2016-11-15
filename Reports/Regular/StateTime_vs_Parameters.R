# rm(list= ls())
# cat("\014")

source("Parameters.R")
source("EstiStats.R")
source("Utilities.R")

StateTimeWrap <- function(.rslsd, QIPrec){
  IntPeriod = Tint
  EstCoef = C
  BeamTime = H*Duty
  
  return(StateTime(IntPeriod, EstCoef, .rslsd, QIPrec, BeamTime))
}

BeamTimeWrap <- function(.rslsd, QIPrec){
  IntPeriod = Tint
  EstCoef = C
  ResolStD = .rslsd
  h = 15*60
  
  return(BeamTime(.rslsd, IntPeriod, EstCoef, h, QIPrec))
}

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}

rsd = c(1%o% 10^(-4:-6)); names(rsd) <- as.character(rsd)

var.beta = 1e-15

Precision <- c(1:10 %o% 10^(-6:-4))

QIprec.trg = 1e-6; trg.i = which(QIprec.trg %in% Precision)

acK = 1 #HAC correction

library(plyr)
ldply(Precision, function(p) BeamTimeWrap(rsd, p)) %>% melt(variable.name="ResSD", value.name="H") %>% 
  cbind(data.frame("Prec" = rep(Precision, times=length(rsd)))) -> Hplot; head(hplot)

ggplot(Hplot, aes(Prec, H/86400, col=ResSD)) + geom_line() +
  scale_x_log10(name=expression(sigma[EA[yy]]~"(a.u.)"), minor_breaks=Precision, labels=fancy_scientific) + 
  scale_y_log10(name="H (days)", minor_breaks=(c(1:10) %o% 10^(-4:4))) + 
  guides(col=guide_legend(title="Resolution")) + 
  theme_minimal() +
  theme(panel.grid.minor = element_line(colour = "gray", linetype = "dashed"), legend.position="top")

x = unlist(hplot); min(x)->minst; max(x) -> maxst;
