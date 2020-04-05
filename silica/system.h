int	nSi;		/* number of Si atoms bonded to OH*/
int nCH3OH;		/* number of methanol molecules*/
int nCH3CN;		/* number of acetonitrile molecules*/
int	tdpoint;
int 	npoints;
double kStr[3], eqBond[3], kBend[2], eqBend[2];/*SiOH and SiOO vibrational parameters*/
double binSize; /*for density profie*/
double VINTSA, VINTSM, VINTSS;
double INTRAV_S, VNB_S;	/*silica intramolecular and intermolecular energies*/
double INTRAV_M, VNB_M;	/*methanol intramolecular and intermolecular energies*/
double INTRAV_A, VNB_A;	/*acetonitrile intramolecular and intermolecular energies*/
double VMA;		/*acetonitrile-methanol interaction*/
double VMS;		/*SiOH-methanol interaction*/
double VAS;		/*SiOH-acetonitrile interaction*/

//#define SiSkip 6   // zero every nth silica's partial charges (turn off by making this > nSi )

ljcon Silj[6][6], liqlj[6][6];
ljcon sollj[5][9]; /*assume 5 different solute atoms*/

/*On the fly Rfd and diff calculations*/
#define binRDF 0.05
double SiLiqRdf[6][3][400];/*liquid-SiOH Rdfs*/
double ccRDF[400];
tripd *Mdensity[3], *Adensity[3];  /* M: 0=O,1=H,2=(CH3). A: 0=C,1=N,0=(CH3) */
double *CMden[2]; /* center of mass density. 0 = CH3OH, 1=CH3CN  */
double *Qden[2]; /* charge density. 0 = CH3OH, 1=CH3CN  */
int HBdt;   /* time step skips for H-bond lifetime  calculations */
//#define HBdt 10
int **SiHB; /* hydrogen bonds */
int **SiHB2;
int *HBcount[6];
int doubleSiA;
int doubleSiD;
/* variables to track non-equilibrium parameters
 * for one methanol molecule runs */
/*double SiHtoMO;  // silica H to MeOH O
double cosMa;  // cosine of angle (if calculated) acceptor MeOH
double SiOtoMOa;  // silica O to MeOH O (M acceptor)
double MHtoSiO;  // silica O to MeOH H
double cosMd;  // cosine of angle (if calculated) donor MeOH
double SiOtoMOd;  // silica O to MeOH O (M donor)
double SiDtoA;   // silica (donating to MeOH) H to ACN N
double SiAtoA;   // silica (accept from MeOH) H to ACN N
int SiOMO; // index of Si with O closest to MeOH O
int SiHMO; // index of Si with H closest to MeOH O
int SiOMH; // index of Si with O closest to MeOH H
int MeOHd; // flag-- is MeOH donating to HB with SiOH
double durMeOHd; // duration of HB with MeOH donor
int MeOHa; // flag-- is MeOH accepting HB from SiOH
double durMeOHa; // duration of HB with MeOH acceptor
double MeOHtau; // min time (ps) req'd, each MeOH HB, to end trajectory
int MeQ; // flag to determine end of trajectory: 2 MeOH HB for MeOHtau time
*/
//int uHB; /* number of unique ACNs h-bonded to SiOH */
/*************************
 *       wxtr.c      
*************************
int *SHBTCF[6];
int *CHBTCF[6];
int normSHB[6];
int normCHB[6];
double *SiDonor;
double *SiAnyM;
double *MtoSi;
double *SitoM;
double *SitoA;
double *SiTot;
double *SiDonAcc;
int *snormOH[NumSBins], *snormOC[NumSBins], *snormCN[NumSBins];
double *sTCFOH[NumSBins], *sTCFOC[NumSBins], *sTCFCN[NumSBins];
int *surM[NumSBins], normSurM[NumSbins];
int *surA[NumSBins], normSurA[NumSbins];
************************/
#define CosBinSize 0.05
#define ZBinSize 0.5  /* orientational distribution binning */
#define HBZBinSize 0.1 /* HB count Z bin size */
//#define TCFZbinsize 1.0
int HBSitoM;
int HBSitoA;
int HBMtoSi;
int **odOC;
int **odOH;
int **odCN;
#define NumSBins 4
int *sodOC[NumSBins], *sodOH[NumSBins], *sodCN[NumSBins];
int TCFdt;   /* time step skips for TCFs and ODs */ 
//#define TCFdt 10
tripfz **MOHvec;
tripfz **MOCvec;
tripfz **ACNvec;
double sbins[2 * NumSBins];
#define taud 500.0  /* max time allowed H-bond to temp separate (in fs) */ 
/* window parameters */
double center_w;
double width_w;
double power_w;
double pot_w;
double osysCmz;
tripd sysCm;
double V_w;
#define CH3OHMass 32.04237
#define CH3CNMass 41.05285
#define OOdist2 (3.4*3.4)  /* square of O-O h-bond max r */
#define NOdist2 (3.5*3.5)  /* square of N-O h-bond max r */
#define cosHOO 0.866025    /* cos of H-O...O h-bond max angle */
#define cosNOH 0.866025    /* cos of N...O-H h-bond max angle */
/*Globals for the electric potential profile calculation*/
int EPdt;
//#define EPdt 10
#define NX 8
#define NY 8
#define NZ 100
double pGrid[4][NZ];
double pGrid2[4][NZ];
double gridQ2[9]; /*coulomb parameters*/
tripd potGrid[NX*NY*NZ];
/*window parameters to keep CH3OH in a slab*/
int setWindow;/*flag set to 1 (2) if the methanol (acetonitrile)  window is activated*/
double W1, W2;/*window boundaries*/
