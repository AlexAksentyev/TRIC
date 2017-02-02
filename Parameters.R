source("TestParameters.R")

Pb = c(U = .4, D = -.4); Pt = .88;
d = 40; L = 18400 #target and accelerator lengths in cm
schottky = 1.1e14
Rho.on = schottky/d; rate = 1.157e-6 # target density dissipation rate
  Rho.off = 0#Rho.on / 2
nu = 793345;

days = 30; H = days*24*3600
Duty = .7

## WATCH OUT FOR TESTING FACTOR
Sigma0 = 70e-27#70e-27 *Factor$CS;
Tint = 1 *Factor$Tint;
MeanSlope = Sigma0*Rho.on*d*nu *Factor$Slope 
p1Dt = MeanSlope*Tint *Factor$LoseProb;
ErrorStD = 1e-7 *Factor$ErrorStD; I0 = 1e-3 *Factor$Current; 
rStD = ErrorStD/I0
N0 = I0/(1.6e-19*nu);

EstCoef <- function(CSec, BeamPol, TargPol, TargDens, TargLen, RevFreq){
  C = -1/(CSec*diff(BeamPol)*TargPol*TargDens*TargLen*RevFreq)
  return(C)
}

C = EstCoef(Sigma0, Pb, Pt, Rho.on, d, nu);

#### SUMMARY #### 
ExpParameters = c(
  Pol = c(beam = Pb, targ = Pt), 
  Rev.f = nu, 
  Targ = c(
    dens = c(on = Rho.on, off = Rho.off), 
    len = d
  )
)