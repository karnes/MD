int tdpoint,npoints;
int nDDC,nw,nEr;
#define DDCsites 12
ljcon WDlj[3][DDCsites];// 3 H2O X 1 Er = 3 combinations
ljcon EDlj[DDCsites], EWlj[3];
ljcon DDlj[DDCsites][DDCsites];

double DDCkstr,DDCreq,DDCkbend,DDCtheq;
double pTors[3], factor14; 
double bondE,torsE,bendE,nbE; // DDC intramolecular potentials

#define rmin 3.204 //Er3+ - SPC/E parameters from Li,Song,& Merz;JCPB;2014
#define EWeps 0.08034231 // kcal/mol

double binSize;
#define binRDF 0.05
#define RDFbins 300
double grSOL[2][RDFbins]; // Er-solvent radial distribution function
// 0 = water O
// 1 = water H

#define WMass 18.0153
#define H2Odensity 0.0334
#define DDCmass 170.34
#define DDCdensity  0.0026497 // based on 0.7495 g/cm3 and MW=170.34 
#define ErMass 167.259

#define rShellO_1 3.50
#define rShellO_2 5.35
int shellO_1,shellO_2;

tripd *H2Oden, *DDCden;

double syscomz; //system CoM on z-axis
double Erz; // Er3+ z position

double DDCBondE, DDCbendE, DDCtorsE;

double V_EW, V_ED, V_WD;
double EW_C, ED_C, WD_C;
double DDCNB, INTRADDC,DD_C;
double H2OC;

double Erz; //z-position of Er ion

double BIAS_W, BIAS_A, BIAS_B, BIAS_C, BIAS_D, BIAS_G;
int windowOn, biasingOn, KillFinger;
double center_w, width_w, pot_w, power_w;
double Gibbs_surf;
double V_w,V_bias;
double osysCmz;

/*
 * Solvent solute switching parameters
 */
double swSoMax;
double swSoMin;
double swSocoef[4];
double pswSocoef[3];

/*
 * Switching charges variables and types
 */

typedef struct	{
 	double	qa;		/* Charge */
 	double	qr;		/* Rate */
}	qswitch;



/*
 * Water finger coordinate
 */
double wcor;
