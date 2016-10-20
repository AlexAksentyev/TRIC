/*
 Ideal gas caclulations oriented on PIT cell

   started 17.05.2014 by SMkrtch
*/

#include <stdio.h>
#include <math.h> // add in Makefile "-lm "
#include <Rcpp.h> //for use in R

using namespace Rcpp;

// Common Constants in SI
static const float NA = 6.0221367E23;//Avogadro´s number [molecules/mol]
static const float RG = 8.3144621;   //universal gas constant [J/(mol K)]
static const float kB = 1.3806E-23 ;//Boltzmann´s const kB = RG/NA [J/K ]

static const float VM = 22.414; //molar volume L/mol at 0 °C

// [[Rcpp::export]]
static float density_Pa2gl(float PrsPa, float TK, float Mgmol)
{ //Prms: PrsPa  Pressure [Pa], TK temperature in Kelvin,
  //      Mgmol Molecular weight im g/mol, in SI kg/mol = 10^3 g/mol
  //return Density [kg/m^3] = [g/liter]

  return (PrsPa * Mgmol * 1.0E-3) / (RG * TK);
}

// [[Rcpp::export]]
static float density_Pa2gccm(float PrsPa, float TK, float Mgmol)
{
  return (PrsPa * Mgmol) / (RG * TK);
}

// [[Rcpp::export]]
static float density_mb2gccm(float Prsmbar, float TK, float Mgmol)
{
  return density_Pa2gl(Prsmbar*100, TK, Mgmol) * 1.0E-3;
}

// [[Rcpp::export]]
static float density_mb2Natccm(float Prsmbar, float TK)
{ // mbar=1.0e2Pa, cm=1.0e-6ccm -> 1.0e-4
  return  (Prsmbar /(kB * TK)) * 1.0e-4;// * 10e-0;
}

// [[Rcpp::export]]
static float massFlow_mbls2gs(float Qmbls, float TK, float Mgmol)
{ // Q_mbarliter/s = 0.1 Q_Pam^3/s 
  return (Qmbls * 0.1) * Mgmol /(RG * TK);
}

// [[Rcpp::export]]
static float massFlow_mbls2Nats(float Qmbls, float TK)
{ // mass Flow Q[mBar*l/s] ->  N[at/s]
  return (Qmbls * 0.1) / (kB * TK);
}

// [[Rcpp::export]]
static float massFlow_Nats2mbls(float Nats, float TK)
{ // mass Flow N[at/s] -> Q[mBar*l/s] 
  return Nats * (kB * TK) / 0.1;
}

// Conductance in [Liter/sec] of the gas molecular flow 
// for circular tube with D: diameter [cm], L: length [cm]
// :
// [[Rcpp::export]]
static float Conduct_LCT(float L, float D, float TK, float Mgmol) 
{//Long tube 
  return 3.81 * sqrt(TK/Mgmol) *(D*D*D / L);
}

// [[Rcpp::export]]
static float Conduct_SCT1(float L, float D, float TK, float Mgmol) 
{//Short tube formula 1
  return 3.81 * sqrt(TK/Mgmol) *D*D*D/(L + 1.33*D);
}

// [[Rcpp::export]]
static float Conduct_SCT2(float L, float D, float TK, float Mgmol) 
{//Short tube formula 2
  return 2.85 * sqrt(TK/Mgmol) *D*D/(1 + 0.75*L/D);
}

// [[Rcpp::export]]
static float CellThickness(float Iabs, float Ctot, float Lb)
{// CellThickness in [Nat/cm^2], Ctot [l/s], Lbeam [cm]
  float dens = (Iabs /(2 * Ctot)) * 1.0e-3; //[Nat/l] * 10^-3 = [Nat/ccm]
  //printf("Dens=%.2e[Nt/ccm]\n",dens);
  return dens*Lb;
} 

static float PZeroCalc()
{//for test cell with 6 probes at X1..X6
  float 
    X0 =  0,   P0 = 2.0e-7, //P[mbar] in target chamber: SVP622
    X1 =  1.0, P1 = 0.0181e-3,//V1 op, Baratron
    X2 =  6.0, P2 = 0.0540e-3,
    X3 = 11.0, P3 = 0.3970e-3,
    X4 = 16.0, P4 = 0.0831e-3,
    X5 = 25.0, P5 = 0.0508e-3,
    X6 = 38.0, P6 = 0.0180e-3,
    X7 = 39.0, P7 = P0,       

    Xz = X3,   Pbz = -0.012e-3,//Batron Zero Pressure: V1-V6 cl,V7 op
    dP, PzL12,PzL02,PzL, PzR40,PzR50,PzR46,PzR,Pz;
/*
  dP = P2 - P1; 
  PzL12 = P1 + (P2 - P1)*((Xz-X1)/(X2-X1));// - P0;
  PzL02 = P0 + (P2 - P0)*((Xz-X0)/(X2-X0));// - P0;
  PzL = (PzL12 + PzL02)/2;
  PzR40 = P0 + (P4 - P0)*((Xz-X7)/(X4-X7));// - P0;
  PzR50 = P0 + (P5 - P0)*((Xz-X7)/(X5-X7));// - P0;
  PzR46 = P6 + (P4 - P6)*((Xz-X6)/(X4-X6));// - P0;
  PzR = (PzR40 + PzR50 + PzR46)/3;
  Pz = (PzL + PzR)/2;
  printf("\nP in the center of cell.\n");
  printf("PzL 02=%.2e, 12=%.2e -> %.2e |",PzL02,PzL12,PzL);
  printf("PzR 40=%.2e, 50=%.2e, 46=%.2e -> %.2e\n",PzR40,PzR50,PzR46,PzR);
  printf(" P0 = %.2e\n",Pz); 
*/

  printf("\nP in the center of cell 29.04.2014.\n");
  PzL12 = P1 + (P2 - P1)*((Xz-X1)/(X2-X1)) - P1;
  PzL02 = P0 + (P2 - P0)*((Xz-X0)/(X2-X0)) - P1;
  PzL = (PzL12 + PzL02)/2;
  PzR40 = P0 + (P4 - P0)*((Xz-X7)/(X4-X7)) - P6;
  PzR50 = P0 + (P5 - P0)*((Xz-X7)/(X5-X7)) - P6;
  PzR46 = P6 + (P4 - P6)*((Xz-X6)/(X4-X6)) - P6;
  PzR = (PzR40 + PzR50 + PzR46)/3;

  Pz = (PzL12 + PzR46)/2;
  printf("PoL1_2=%.2e, Po4_6=%.2e -> Po = %.2e\n",PzL12,PzR46,Pz);
  
  printf("\nPo in the center of cell 23.06.2014.\n");
  //  dFT=0.8;    dBb=dBc=0.98; //cm
  //   L=13.71;  11.0 28.0
  P0 = 1.7e-7; //
  P1 = 0.0384e-3; // (P1+P6)/2 = 0.0443e-3;
  P2 = 0.0782e-3;
  P3 = 0.5180e-3;
  P4 = 0.1032e-3;
  P5 = 0.0796e-3;
  P6 = 0.0503e-3;
  P7 = P0;       

  PzL12 = P1 + (P2 - P1)*((Xz-X1)/(X2-X1)) - P1;
  PzR46 = P6 + (P4 - P6)*((Xz-X6)/(X4-X6)) - P6;
  PzR = (PzR40 + PzR50 + PzR46)/3;

  Pz = (PzL12 + PzR46)/2;
  printf("PoL1_2=%.2e, PoR4_6=%.2e -> Po = %.2e\n",PzL12,PzR46,Pz);

  float C,Po,Q; //Ca,Cb,Cc,

  C = 9.39;//l/s FT 0.8
  Po = 7.03e-5;//HJul01//8.53e-5;//HMay//mbar
  Q = massFlow_mbls2Nats(C*Po,300);
  printf("\"H_1/2\" FT d80: Ctot = %.2e, Po = %.2e, ABS Flux=%.2e\n",C,Po,Q);
/*
  C = 11.23;//l/s FT 0.8
  Po = 7.58e-5;//mbar
  Q = massFlow_mbls2Nats(C*Po,300);
  printf("\"H_1/2\" FT d98: Ctot = %.2e, Po = %.2e, ABS Flux=%.2e\n",C,Po,Q);
*/
  C = 11.23;// FTd0.98; //9,39;//l/s FTd0.80
  Po = 10.17e-5;//8.53e-5;//mbar
  Q = massFlow_mbls2Nats(C*Po,300);
  printf("\"H_1/2\" FT d98: Ctot = %.2e, Po = %.2e, ABS Flux=%.2e\n",C,Po,Q);

  //C = 11.23 / 1.41421356; // Deiterium//l/s FT 0.98
  C =  9.39 / 1.41421356; // Deiterium//l/s FT 0.8
  Po =10.12e-5;//mbar
  Q = massFlow_mbls2Nats(C*Po,300);
  printf("\"D_2/6\" FT d80: Ctot = %.2e, Po = %.2e, ABS Flux=%.2e\n",C,Po,Q);



/*

Po in the center of cell 23.06.2014.

  //  dFT=0.8;    dBb=dBc=0.98; //cm
  //   L=13.71;  11.0 28.0
  P0 = 1.7e-7; //
  P1 = 0.0384e-3; // (P1+P6)/2 = 0.0443e-3;
  P2 = 0.0782e-3;
  P3 = 0.5180e-3;
  P4 = 0.1032e-3;
  P5 = 0.0796e-3;
  P6 = 0.0503e-3;
  P7 = P0;       

PoL1_2=6.78e-05, PoR4_6=7.23e-05 -> Po = 7.00e-05

P in the center of cell 29.04.2014.
PoL1_2=7.18e-05, Po4_6=7.99e-05 -> Po = 7.58e-05


FT: L=13.71     D=0.98  C1=4.11
Bb: L=11.00     D=0.98  C1=5.01
Bc: L=28.00     D=0.98  C1=2.11
Conductance [l/s]: Tot= 11.23; BC= 7.12
ABS flux=4.00e+16 [at/s], M=1.0080; Cell at T=300.0 °K, t_BC=6.95e+13 [at/cm^2]

FT: L=13.71     D=0.80  C1=2.27
Bb: L=11.00     D=0.98  C1=5.01
Bc: L=28.00     D=0.98  C1=2.11
Conductance [l/s]: Tot= 9.39; BC= 7.12
ABS flux=4.00e+16 [at/s], M=1.0080; Cell at T=300.0 °K, t_BC=8.31e+13 [at/cm^2]

FT: L=13.71     D=0.60  C1=0.98
Bb: L=11.00     D=0.98  C1=5.01
Bc: L=28.00     D=0.98  C1=2.11
Conductance [l/s]: Tot= 8.09; BC= 7.12
ABS flux=4.00e+16 [at/s], M=1.0080; Cell at T=300.0 °K, t_BC=9.64e+13 [at/cm^2]


*/
  return 0;
}

struct tcell {// in cm for Cell geometry
  float FTd;
  float FTl;
  float Bbd;
  float Bbl;
  float Bcd;
  float Bcl;
  float xP[6]; //Position of the Mes. points
};

tcell CellGeom = {0.98, 13.71, 0.98, 11.00, 0.98, 28.00,
		  //   x1   x2    x3    x4    x5    x6:  Position of the P sensors
		     {1.0, 6.0, 11.0, 16.0, 25.0, 38.0}};

int  MesOrd[6] = {1,6,5,2,4,3};

struct tBarMes {
  float Pb[14]; //Baratron data: 0,2,4..= P7, 1,3,5..= Px
  float Pbkg;   //background Averaged over 6 sensors with no gas load
  //tcell cell;   
  float FTd;
  char Comment[128];
};

tBarMes Jun04 = { // x e-7
  //p7    p1     p7    p6     p7    p5    p7    p2     p7    p4      p7     p3    p7
  {-3.0, 342.0, -5.0, 265.0, -2.0, 709.0, 4.0, 787.0, -4.0, 1000.0, -5.0, 8040.0, 6.0},
  220.0, // Pbkg
  0.8,   //  FTd;
  "Jun04. ABS Pos: X=+2; Z=0" //Comment
};

static void PoCalc(tBarMes *Dat) 
{ 
  //float Px[6];
  int i,j;
  printf("\nComment: %s\n",Dat->Comment);
  for (i = 0; i<7; i++) {
    j = i*2;
    //Px[MesOrd[i]] = (Dat->Pb[j+1] -((Dat->Pb[j] + Dat->Pb[j+2])/2) - Dat->Pbkg)*1.0e-7;
    printf("[%d/%d]: Px = %.3e, Pofs(%.3e/%.3e),Pbk = %.3e \n", 
	   i,j,Dat->Pb[j+1],Dat->Pb[j], Dat->Pb[j+2],Dat->Pbkg);   

    //printf("P[%d] = %.3e  ",i, Px[MesOrd[i]]); 
  }

  printf("\n");
}

int main(int argc,char **argv)
{
  float 
    PNmbX = 1013.25,//mbar
    TN    = 273.15, //K -> 0°C
  //Atomic MHat but Molecular MHat*2
    MHat  = 1.008,  //Molar mass [g/mol] Hydrogen Atoms
    MDat  = 2.014,  //Molar mass [g/mol] Deuterium Atoms
    MCH4at= 16.04,  //Molar mass [g/mol] Methane 
    VM    = 22.414;    //molar volume [L/mol]
  //(0 °C, 101.325 kPa) 0.08988 g/L

#ifdef GAS_CALC_COMMON
  printf("Molecular Hydrogen density (M=%.3f) at T=%.2f°K, P=%.2fmbar\n %e g/l, ",
	 MHat, TN, PNmbX,
	 density_Pa2gl(PNmbX*100,TN, 2*MHat));

  printf(" %e g/ccm, %e Nat/ccm.\n",
	 density_mb2gccm(PNmbX,TN, 2*MHat), 
	 density_mb2Natccm(PNmbX,TN));

  printf("Molar volume = %f L/mol\n\n", NA/density_mb2Natccm(PNmbX,TN)/1000);

  printf("Atomic Hydrogen density (M=%.3f) at T=%.2f°K, P=%.2fmbar\n %e g/l, ",
	 MHat, TN, PNmbX,
	 density_Pa2gl(PNmbX*100,TN, MHat));
  printf(" %e g/ccm, %e Nat/ccm.\n",
	 density_mb2gccm(PNmbX,TN, MHat), 
	 density_mb2Natccm(PNmbX,TN));

  float Pc,
    P0 = 7.6e-5, //P [mbar] in central cell point
    lbc = 39,   // Cell length 39 cm
    dbc = 1.0;  // Cell diam   1.0 cm
  Pc = P0/2;
  printf("\nCell: L=%f[cm], D=%.2f[cm]\n", lbc, dbc);
  printf("Atomic Hydrogen density (M=%.3f) at T=%.2f°K, P=%.2embar\n %e g/l, ",
	 MHat, TN, Pc,
	 density_Pa2gl(Pc*100,TN, MHat));
  printf(" %e g/ccm, %e Nat/ccm. \n thickness %e at/cm^2\n",
	 density_mb2gccm(Pc,TN, MHat), 
	 density_mb2Natccm(Pc,TN),
	 density_mb2Natccm(Pc,TN)*lbc);
#endif //#ifdef GAS_CALC

  float 
    MGAt,       //Molar mass [g/mol] of used Gas
    Iabs,       //Feeding flux from ABS [at/s]
    //Tubes  A: Feeding, B:, C: beam tubes
    La, Da, Ca, //Tube A: Len., Dia. Conduct.
    Lb, Db, Cb, //Tube B  
    Lc, Dc, Cc, //Tube C   
    C3,T3 ;

  Iabs = 3.0e16;//4.0e16;
  La = 13.71; 
  Lb = 11.00; 
  Lc = 28.00; 
  TN = 300; MGAt = MHat; //Molar mass [g/mol] of used Gas

  Da = 1.5;//1.0;
  Db = 1.5;//0.96;
  Dc = Db;

  Ca = Conduct_SCT2(La,Da,TN,MGAt);
  Cb = Conduct_SCT2(Lb,Db,TN,MGAt);
  Cc = Conduct_SCT2(Lc,Dc,TN,MGAt);
  C3 = Ca+Cb+Cc;
  T3 = CellThickness(Iabs,C3,Lb+Lc);
  printf("\nFT: L=%.2f\tD=%.2f\tC1=%.2f\n", La,Da,Ca);
  printf("Bb: L=%.2f\tD=%.2f\tC1=%.2f\n", Lb,Db,Cb);
  printf("Bc: L=%.2f\tD=%.2f\tC1=%.2f\n", Lc,Dc,Cc);
  printf("Conductance [l/s]: Tot= %.2f; BC= %.2f\n", C3, Cb+Cc);
  printf("ABS flux=%.2e [at/s], M=%.4f; Cell at T=%.1f °K, t_BC=%.2e [at/cm^2]\n", Iabs,MGAt,TN,T3);


/*
  Da = 0.8;
  Db = 0.98;
  Dc = Db;
float  Po =10.12e-5;//mbar
  Ca = Conduct_SCT2(La,Da,TN,MDat);
  Cb = Conduct_SCT2(Lb,Db,TN,MDat);
  Cc = Conduct_SCT2(Lc,Dc,TN,MDat);
float  C = (Ca+Cb+Cc);// / 1.41421356;
float  Q = massFlow_mbls2Nats(C*Po,300);
  printf("\"D_2/6\" FT d98: Ctot = %.2e, Po = %.2e, ABS Flux=%.2e\n",C,Po,Q);
*/
  //  PZeroCalc();

  //PoCalc(&Jun04);
  return 0;
}

/*
FT: L=13.71     D=1.18  C1=7.05
Bb: L=11.00     D=1.18  C1=8.57
Bc: L=28.00     D=1.18  C1=3.64
Conductance [l/s]: Tot= 19.26; BC= 12.21
ABS flux=4.00e+16 [at/s], M=1.0080; Cell at T=300.0 °K, t_BC=4.05e+13 [at/cm^2]
FT: L=13.71     D=0.80  C1=2.27
Bb: L=11.00     D=1.18  C1=8.57
Bc: L=28.00     D=1.18  C1=3.64
Conductance [l/s]: Tot= 14.48; BC= 12.21
ABS flux=4.00e+16 [at/s], M=1.0080; Cell at T=300.0 °K, t_BC=5.39e+13 [at/cm^2]

Target thicknes t_bc as f(Da): 
   FT Da = 1.18              4.05  e+13  gain in thickness
           1.1       4.05 -> 4.34 e+13   +9%  
           1.0       4.05 -> 4.71 e+13   +18%
           0.9       4.05 -> 5.06 e+13   +25% 
           0.8       4.05 -> 5.39 e+13   +30% 
           0.7       4.05 -> 5.68 E+13   +40%
           0.6       4.05 -> 5.92 e+13   +45%

ANKE ABS paper in Nucl. Instr. Meth. A 721, p. 83-89, (2013)
Hydrogen H1:   7.5⋅10^16 atoms/s (two hyperfine states) 
deuterium D1:  3.9⋅10^16 atoms/s (three hyperfine states)

Exp.#219 pp-> K Lamda request:
p-pol_beam:  5⋅10^9   1/s      stacking: 2⋅10^10 1/s
PIT:         3⋅10^13  at/cm^2

Proposal Exp.#219: Cell Beam Tube D ~ 10mm (stochastic cooled p-beam)
Nov.2014 test:
 Openable Cell Stainless Steel wall thickness 0.1mm: 
  Da=11.8mm, La=130mm, Db=Dc=11.8mm, Lb=110, Lc=280mm, Lx=-38 (Baratron) 
 Obtained PIT thickness ~ 3⋅10^13  at/cm^2 for H1 and D1
 Baratron Pm = (0.08 - 0.12)⋅10^-3 mbar, Pzero = 0.05⋅10^-3 mbar 
Lets calculate target thickness.
 so, in cell Px = <Pm> - Pzero = (0.1 - 0.05)⋅10^-3 = 5.0⋅10^-5 mbar 
 Pressure in target chamber Pch = 2⋅10^-7 mbar. (can be neglected!)
Lb - Lx = 11.0 - 3.8 = 7.2 cm dP = Px - Pch ~= Px.
Pressure in centrum of the cell:
     Po = Px⋅Lb/(Lb-Lx) = 5.0⋅10^-5 1.528 = 7.64⋅10^-5 mbar
The same value one can get by calculations using conductances Cbx and Cx.
Averaged pressure in the cell is <P>=0.5⋅Po⋅(Lb+Lc)/(Lb+Lc)=0.5⋅Po.
So, P = 3.8⋅10^-5 mbar;
Density at T=300°K: dn[at/cm^2] = P[mbar]⋅10^-5/kT[°K] = P⋅7.2432⋅10^18/T = 
              = 9.2⋅10^11 at/cm^2.
Target thickness t = dn⋅(Lb+Lc) = 9.2⋅10^11⋅39 = 
              = 3.59⋅10^13 at/cm^2
ABS flux Iabs estimation:
  Averaged gas debsity in cell is:
       dn[at/cm^3] = Iabs[at/s]/(10^3⋅Ctot[l/s]), 
 where Ctot is total conductance of the cell tubes Ctot = Ca+Cb+Cc = 
       Ctot = 7.05 + 8.57 + 3.64 = 19.26 l/s.
 So, ABS flux is
  Iabs[at/s] = dn[at/cm^3]⋅10^3⋅Ctot[l/s] = 9.2 ⋅10^11⋅ 10^3 ⋅19.26 = 177⋅10^14=
             = 1.77⋅10^16 at/s

 */

/*
FT: L=13.71     D=1.50  C1=14.08
Bb: L=11.00     D=1.50  C1=17.02
Bc: L=28.00     D=1.50  C1=7.38
Conductance [l/s]: Tot= 38.48; BC= 24.39
ABS flux=4.00e+16 [at/s], M=1.0080; Cell at T=300.0 °K, t_BC=2.03e+13 [at/cm^2]

FT: L=13.71     D=1.00  C1=4.36
Bb: L=11.00     D=1.50  C1=17.02
Bc: L=28.00     D=1.50  C1=7.38
Conductance [l/s]: Tot= 28.75; BC= 24.39
ABS flux=4.00e+16 [at/s], M=1.0080; Cell at T=300.0 °K, t_BC=2.71e+13 [at/cm^2]

FT: L=13.71     D=0.80  C1=2.27
Bb: L=11.00     D=1.50  C1=17.02
Bc: L=28.00     D=1.50  C1=7.38
Conductance [l/s]: Tot= 26.67; BC= 24.39
ABS flux=4.00e+16 [at/s], M=1.0080; Cell at T=300.0 °K, t_BC=2.93e+13 [at/cm^2]

FT: L=13.71     D=0.60  C1=0.98
Bb: L=11.00     D=1.50  C1=17.02
Bc: L=28.00     D=1.50  C1=7.38
Conductance [l/s]: Tot= 25.37; BC= 24.39
ABS flux=4.00e+16 [at/s], M=1.0080; Cell at T=300.0 °K, t_BC=3.07e+13 [at/cm^2]

Cell Dbc = 1.5 cm 
     DFT      Dens.*e+13 [at/cm^2] 
     1.5      2.03                   1.00
     1.0      2.71                   1.33
     0.8      2.93                   1.44
     0.6      3.07                   1.55
*/

/*
FT: L=13.71     D=1.18  C1=7.05
Bb: L=11.00     D=1.18  C1=8.57
Bc: L=28.00     D=1.18  C1=3.64
Conductance [l/s]: Tot= 19.26; BC= 12.21
ABS flux=4.00e+16 [at/s], M=1.0080; Cell at T=300.0 °K, t_BC=4.05e+13 [at/cm^2]

FT: L=13.71     D=1.00  C1=4.36
Bb: L=11.00     D=1.18  C1=8.57
Bc: L=28.00     D=1.18  C1=3.64
Conductance [l/s]: Tot= 16.57; BC= 12.21
ABS flux=4.00e+16 [at/s], M=1.0080; Cell at T=300.0 °K, t_BC=4.71e+13 [at/cm^2]

FT: L=13.71     D=0.80  C1=2.27
Bb: L=11.00     D=1.18  C1=8.57
Bc: L=28.00     D=1.18  C1=3.64
Conductance [l/s]: Tot= 14.48; BC= 12.21
ABS flux=4.00e+16 [at/s], M=1.0080; Cell at T=300.0 °K, t_BC=5.39e+13 [at/cm^2]

FT: L=13.71     D=0.60  C1=0.98
Bb: L=11.00     D=1.18  C1=8.57
Bc: L=28.00     D=1.18  C1=3.64
Conductance [l/s]: Tot= 13.18; BC= 12.21
ABS flux=4.00e+16 [at/s], M=1.0080; Cell at T=300.0 °K, t_BC=5.92e+13 [at/cm^2]

Cell Dbc = 1.18 cm 
     DFT      Dens.*e+13 [at/cm^2] 
     1.18     4.05                   1.00
     1.0      4.71                   1.16
     0.8      5.39                   1.33
     0.6      5.92                   1.46
*/

/*
FT: L=13.71     D=1.00  C1=4.36
Bb: L=11.00     D=0.96  C1=4.72
Bc: L=28.00     D=0.96  C1=1.98
Conductance [l/s]: Tot= 11.06; BC= 6.70
ABS flux=4.00e+16 [at/s], M=1.0080; Cell at T=300.0 °K, t_BC=7.05e+13 [at/cm^2]

FT: L=13.71     D=0.96  C1=3.87
Bb: L=11.00     D=0.96  C1=4.72
Bc: L=28.00     D=0.96  C1=1.98
Conductance [l/s]: Tot= 10.57; BC= 6.70
ABS flux=4.00e+16 [at/s], M=1.0080; Cell at T=300.0 °K, t_BC=7.38e+13 [at/cm^2]

FT: L=13.71     D=0.80  C1=2.27
Bb: L=11.00     D=0.96  C1=4.72
Bc: L=28.00     D=0.96  C1=1.98
Conductance [l/s]: Tot= 8.98; BC= 6.70
ABS flux=4.00e+16 [at/s], M=1.0080; Cell at T=300.0 °K, t_BC=8.69e+13 [at/cm^2]

FT: L=13.71     D=0.60  C1=0.98
Bb: L=11.00     D=0.96  C1=4.72
Bc: L=28.00     D=0.96  C1=1.98
Conductance [l/s]: Tot= 7.68; BC= 6.70
ABS flux=4.00e+16 [at/s], M=1.0080; Cell at T=300.0 °K, t_BC=1.02e+14 [at/cm^2]

Cell Dbc = 0.96 cm 
     DFT      Dens.*e+13 [at/cm^2] 
     1.0      7.05                   1.00
     0.96     7.38                   1.05
     0.8      8.69                   1.23
     0.6     10.20                   1.45

*/
#ifdef NOCALK
Atomic Hydrogen M=1.0080; Cell at T=300.0 °K
ABS flux=4.00e+16 [at/s], 

Cell Dbc = 1.5 cm                  
     DFT      Dens.*e+13 [at/cm^2]  rDens  rIabs  rDens*rIabs
     1.5      2.03                   1.00  1.00*    1.00      rDens
     1.0      2.71                   1.33  1.00     1.32       1.00  1.00
     0.8      2.93                   1.44  0.92     1.32       1.08  0.99
     0.6      3.07                   1.55  0.84     1.30       1.13  0.95

Cell Dbc = 1.18 cm 
     DFT      Dens.*e+13 [at/cm^2] 
     1.18     4.05                   1.00  1.00*    1.00
     1.0      4.71                   1.16  1.00     1.16       1.00  1.00
     0.8      5.39                   1.33  0.92     1.22       1.14  1.05
     0.6      5.92                   1.46  0.84     1.22       1.26  1.06

Cell Dbc = 0.96 cm 
     DFT      Dens.*e+13 [at/cm^2] 
     1.0      7.05                   1.00  1.00     1.00
     0.96     7.38                   1.05  0.99*    1.04
     0.8      8.69                   1.23  0.92     1.13
     0.6     10.20                   1.45  0.84     1.21

ABS beam profile has FWHM = 6.28mm = 2.355 sigma; 
Sigma S = 2.667 mm. Int(6S) = 1.0;
FT 
D[cm]   S    +/-S   F(U)/2   F(U)  rIabs
1.0   3.73   1.87   0.469   0.938  1.00
0.9   3.37   1.68   0.462   0.924  0.98
0.8   3.00   1.50   0.433   0.866  0.92
0.7   2.62   1.31   0.405   0.810  0.86
0.6   2.25   1.12   0.394   0.788  0.84
0.5   1.87   0.94   0.325   0.650  0.69
0.4   1.50   0.75   0.223   0.446  0.47
   
#endif

#ifdef BT_D_15

     Cell BT D=1.5cm; t_BC = f(FT_D)

FT: L=13.71     D=1.50  C1=14.08
Bb: L=11.00     D=1.50  C1=17.02
Bc: L=28.00     D=1.50  C1=7.38
Conductance [l/s]: Tot= 38.48; BC= 24.39
ABS flux=3.00e+16 [at/s], M=1.0080; Cell at T=300.0 °K, t_BC=1.52e+13 [at/cm^2]

FT: L=13.71     D=1.20  C1=7.40
Bb: L=11.00     D=1.50  C1=17.02
Bc: L=28.00     D=1.50  C1=7.38
Conductance [l/s]: Tot= 31.79; BC= 24.39
ABS flux=3.00e+16 [at/s], M=1.0080; Cell at T=300.0 °K, t_BC=1.84e+13 [at/cm^2]

FT: L=13.71     D=1.00  C1=4.36
Bb: L=11.00     D=1.50  C1=17.02
Bc: L=28.00     D=1.50  C1=7.38
Conductance [l/s]: Tot= 28.75; BC= 24.39
ABS flux=3.00e+16 [at/s], M=1.0080; Cell at T=300.0 °K, t_BC=2.03e+13 [at/cm^2]

FT: L=13.71     D=0.98  C1=4.11
Bb: L=11.00     D=1.50  C1=17.02
Bc: L=28.00     D=1.50  C1=7.38
Conductance [l/s]: Tot= 28.50; BC= 24.39
ABS flux=3.00e+16 [at/s], M=1.0080; Cell at T=300.0 °K, t_BC=2.05e+13 [at/cm^2]

FT: L=13.71     D=0.80  C1=2.27
Bb: L=11.00     D=1.50  C1=17.02
Bc: L=28.00     D=1.50  C1=7.38
Conductance [l/s]: Tot= 26.67; BC= 24.39
ABS flux=3.00e+16 [at/s], M=1.0080; Cell at T=300.0 °K, t_BC=2.19e+13 [at/cm^2]


FT: L=13.71     D=1.00  C1=4.36
Bb: L=11.00     D=0.80  C1=2.78
Bc: L=28.00     D=0.80  C1=1.15
Conductance [l/s]: Tot= 8.29; BC= 3.94
ABS flux=3.00e+16 [at/s], M=1.0080; Cell at T=300.0 °K, t_BC=7.05e+13 [at/cm^2]

FT: L=13.71     D=1.00  C1=4.36
Bb: L=11.00     D=1.00  C1=5.32
Bc: L=28.00     D=1.00  C1=2.23
Conductance [l/s]: Tot= 11.91; BC= 7.55
ABS flux=3.00e+16 [at/s], M=1.0080; Cell at T=300.0 °K, t_BC=4.91e+13 [at/cm^2]


FT: L=13.71     D=1.20  C1=7.40
Bb: L=11.00     D=1.50  C1=17.02
Bc: L=28.00     D=1.50  C1=7.38
Conductance [l/s]: Tot= 31.79; BC= 24.39
ABS flux=3.00e+16 [at/s], M=1.0080; Cell at T=300.0 °K, t_BC=1.84e+13 [at/cm^2]

FT: L=13.71     D=1.00  C1=4.36
Bb: L=11.00     D=1.50  C1=17.02
Bc: L=28.00     D=1.50  C1=7.38
Conductance [l/s]: Tot= 28.75; BC= 24.39
ABS flux=3.00e+16 [at/s], M=1.0080; Cell at T=300.0 °K, t_BC=2.03e+13 [at/cm^2]


#endif //#ifdef BT_D_15

#ifdef METHANE
  Dcell = 10mm; 
  Gas: METHANE CH4
FT: L=50.00     D=0.10  C1=0.00
Bb: L=11.00     D=1.00  C1=1.33
Bc: L=28.00     D=1.00  C1=0.56
Conductance [l/s]: Tot= 1.89; BC= 1.89
ABS flux=3.00e+16 [at/s]-> 1.24e-03 [mBar*l/s], 
M=16.0400; Cell at T=300.0 °K, t_BC=3.09e+14 [at/cm^2]


FT: L=50.00     D=0.10  C1=0.00
Bb: L=11.00     D=1.00  C1=3.76
Bc: L=28.00     D=1.00  C1=1.58
Conductance [l/s]: Tot= 5.34; BC= 5.34
ABS flux=3.00e+16 [at/s]-> 1.24e-03 [mBar*l/s], 
M=1.0080; Cell at T=300.0 °K, 
	    H2:    t_BC=1.10e+14 [at/cm^2]
	  * atoms: t_BC=7.75e+13 [at/cm^2]


  Dcell = 20mm;
  Gas: METHANE CH4
FT: L=50.00     D=0.10  C1=0.00
Bb: L=11.00     D=2.00  C1=9.62
Bc: L=28.00     D=2.00  C1=4.29
Conductance [l/s]: Tot= 13.91; BC= 13.91
ABS flux=3.00e+16 [at/s]-> 1.24e-03 [mBar*l/s], 
M=16.0400; Cell at T=300.0 °K, t_BC=4.21e+13 [at/cm^2]


FT: L=50.00     D=0.10  C1=0.00
Bb: L=11.00     D=2.00  C1=27.13
Bc: L=28.00     D=2.00  C1=12.09
Conductance [l/s]: Tot= 39.23; BC= 39.23
ABS flux=3.00e+16 [at/s] -> 1.24e-03 [mBar*l/s], 
M=2.0160; Cell at T=300.0 °K, 
	    H2:    t_BC=1.49e+13 [at/cm^2]
	  * atoms: t_BC=1.05e+13 [at/cm^2]
				 
     L = 1.0e+10 * 1.0e+6 * t_BC [at/cm^2] ~ 1.0e+29 

#endif//METHANE
