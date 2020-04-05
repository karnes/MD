int	nNIT;		/* number of nitrobenzene (NIT) molecules*/
int	tdpoint;
int	neqFlag;
int	Estate;		/* flag for solute electronic state	*/
double	NITV;   	/* intramolecular potential energy of NIT*/
double	NITNB;   	/* intermolecular potential energy of NIT*/
double	pTors[8][3];	/* torsion potential parameters		*/
double	kBend[5], eqBend[5];	/* NIT bending potential parameters	*/
double	kStr[5], eqBond[5];	/* NIT stretching potential parameters	*/
double  bondE;		/* bond stretch potential energy	*/
double  bendE;		/* angle bending potential energy	*/
double  torsE;		/* angle torsion potential energy	*/
double  impE;		/* 1-4 non-bond energy			*/
double  Q14, A14, C14;	/* parameters for non-bond interactions */
ljcon   NITlj[19][19];	/* lj - coul for NIT and NIT/solute     */

double WNITC;		/* water-nit coulombic */
double NITC;		/* NIT-NIT coulombic */
double H2OC;		/* water-water coulombic */
double VINTH2O;         /* ion-water interaction */
double VCH2O[2];	/* 1st and 2nd H2O shell coulombic interactions */
double VINT_NIT;  	/* ion-nitrobenzene interaction */
double VC_NIT;
double grSOL[8][300];/* on the fly solute-solvent rdfs*/
double EXFIELD;
double H2ONITV;		/* water - nitrobenzene potential energy */
double V_teth;		/* energy used to tether fixShel water to ion */
//double tteth1, tteth2, length;  /* tethering tests */
//int bump;
//double totMass;

#define NITMass 123.1094
#define WMass 18.0153
#define ionMass 35.453
#define H2Odensity 0.0334
#define NITdensity 0.005865
double HQ, OQ;  //spc water charges.


double wcor; // w (Morita) coordinate
/* 
 * densitiies 
 */
tripd *H2Oden;
tripd *NITden;

int	npoints;
double	binSize;

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


double Es, Eb, Et, Eimt, Eintranb;
/*
 * For the electric potential
  */
int EFieldOn;             /* If 1 add constant electric potential */
double EForce[19];
double EVNIT;             /* interactions of NIT with electric field*/
double EVH2O;             /* interactions of H2O with electric field*/
double EVION[2];             /* interaction of Ion with electric field*/
int EFieldQ[3];           /* turn on E field on the water, NIT and ion*/

/* For windowing*/
double center_w;                /* center of window (A)         */
double width_w;                 /* width of window (A)          */
double power_w;                 /* power of window potential    */
double pot_w;                   /* window potential strength    */
double osysCmz;                 /* to check change in sys c.o.m */
double V_w;                 /* windowing energy */
/*On the fly Rfd  orientations and hydrogen-bonds calculations*/
 
#define tcData 5001 // numDP +1 
#define maxW 1001
#define binRDF 0.05
double H2ONITRdf[9][500];
int biasingOn;
int KillFinger;
int fixShel; // number of water molecules to tether to the ion
int fixw[10]; //indices of waters tethered to ion
int fixed; //number of water molecules tethered to ion
double V_bias;
double VINTH2O_s; /*total interaction with 1st shell water*/
int     Nshel; /*number of water molecules in first shell*/
double  Cshel[maxW]; /*orient of water dipole wr water-ion*/
double Rshel; /*radius of first hydraton shell*/
tripd poshist[tcData];
double  sqDiv[3][tcData];
double  sqDivNorm[tcData];
double  pwdip[101];
int     wdipNorm;
#define rODrad 10.0
#define rODbinw 0.2
#define rODbins ((int)(rODrad/rODbinw))+1
#define zODdist 20.0
#define zODbinw 0.5
#define zODbins ((int)(zODdist/zODbinw))+1
#define NITzODdist 20.0
#define NITzODbinw 0.5
#define NITzODbins ((int)(NITzODdist/NITzODbinw))+1
double  zODdip[zODbins][101];
double     zdipNorm[zODbins];
double  NITzOD[NITzODbins][101];
double     NITzNorm[NITzODbins];
double  pwdip2[rODbins][101];
int     wdipNorm2[rODbins];
int     Inshel[maxW][tcData];
double  dipmat[maxW][tcData];
double  CorSh[tcData];
double  CorNorm[tcData];
double  avCorSh[tcData];
double  Cordip[tcData];
double  avCordip[tcData];
double  hhCor[3][tcData]; /* water HH correlations*/
tripd   hhV[maxW][tcData];/*tine history - HH vector for water i at time t*/
int outSh[maxW];
// for snapshots
double snap;
#define snap_w 0.7
