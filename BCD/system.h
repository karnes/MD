#define BrOs 9
#define BCDs 147
int nBCD, nBrO;
int tdpoint; 
int npoints;

#define factor14 0.5
#define BCDfactorlj 0.5
#define BCDfactorq 0.8333
ljcon   slj[BrOs+3][BrOs+3]; /* lj - coul for water - Br octane */   
ljcon	BCDlj[21][21];	/* lj - coul for BCD     */
ljcon   wBCDlj[21][3];  /* lj - coul for water-BCD */
ljcon   bBCDlj[21][BrOs]; /* lj - coul for water-BRB */

double BeqBond[4];	/* Equilibrium bond length for BrO*/
double BkStr[4];	/* Stretching force constants for BrO*/
double BeqBend[4];	/* Equilibrium bond angles for BrO*/
double BkBend[4];	/* Bending force constants for BrO*/
double BpTors[4][4];	/* Bending force constants for BrO*/

double eqBond[6];	/* Equilibrium bond length for BCD*/
double kStr[6];		/* Stretching force constants for BCD*/
double eqBend[11];	/* Equilibrium bond angles for BCD*/
double kBend[11];	/* Bending force constants for BCD*/
double pTors[9][3];	/* Bending force constants for BCD*/
int pair[63][2];	/* 1-4 L-J pairs (intra-BCD) 
			   (63 per sugar unit) */
int pairNB[9849][2];	/* 1-5 or further L-Jq 'non-bonded' pairs (intra-BCD)*/

double sn2BCDrad;       /* BCD-sn2 com-com distance */
tripd BCDcom;		/* BCD center of mass */
tripd BCDz;		/* BCD molecular vector (secondary to primary) */
tripd BCDx; 		/* vector from CoM to center of glucose unit 1 */
tripd BCDy;		/* cross product of BCDx and BCDz */
tripd pBrO;		/* C8 to Br vector of pore BrO (if c8 in pore; else 0,0,0) */ 
double cosBCDpBr;		/* cosine of pBrO ~ BCDz angle */
double capR,Ms,Ml;
int guestBiasOn;
double gbA, gbB, gbO,gBias;
double V_gAng;

/* On the fly RDF and diff calculations */
#define BrOdensity 0.0034861
#define H2Odensity 0.0334
#define WMass 18.0153
#define BMass 193.128
#define binRDF 0.05
#define RDFbins 400
#define binD 0.1
double binSize;
tripd *BCDden, *BrOden[3], *H2Oden;

double U_teth;
double prevBCDcos;
double guestZ,gaC,gC8,g_ang; //guest CoM,aC,C8 pos on BCDz (p) axis, ac-C8 vs p angle
/*
 * hydrogen bonding calculations
 */
#define rOOmax 3.4
#define cosHOO 0.8660254
// max time HB broken to be 'really' broken (ps)
#define HBtau 1.0
#define iHBtau 0.0
// minimum HB bond length (ps)
#define mnHBtau 0.0
#define mniHBtau 0.0
int priBCDHB, iuBCDHB, BCDHB, BCDWHB1, BCDWHB2,  WWHB; /* number of hydrogen bonds detected */
int priD,priA,secD,secA,watHBtot;
#define nBCDHB (3*7*2+3*7)
#define HBhistlen 2000
int HBhist[5][HBhistlen],HBhistn[5];
int *HB[nBCDHB]; // BCD hydrogen bonding 

/*
 * interaction potentials
 */
double INTRABCD, BCDNB;
double VCH2O, VINTH2O, VINTH2Ooh; /* BCD-H2O potential */
double H2OBrOV, BrONB, BrOV;
double VC_BrO, VINT_BrO; 
double VINT_BCD; // BCD-solvent interaction
double bondE, bendE, torsE, nb14E,nb15E;
double BrObondE,BrOtorsE,BrObendE,BrOnb14E,BrOnb15E;

/* on the fly solute-solvent rdfs*/
double grSOL[8][RDFbins]; /* BCD CoM - solvent RDFs: 
		0 = water CoM
		1 = Br-octane CoM
		2 = Br-octane Br atom  */
double grCl[3][4][RDFbins]; /* SN2 Cl rdfs
   0 = CH3	0 = aC
   1 = Cl(1)	1 = Br
   2 = Cl(2)	2 = bCD OH
		3 = water O  */
#define rWatCl 3.85
#define rOHCl 3.85
#define rBrCl 6.20
#define raCCl 5.20
int Clsol[3][4]; // tally of solvent in solvation shell, same indexing
#define contRad 12.0
#define contBin 0.2
#define xyMapx 12.0
#define xyMapy 12.0
#define xyMapz 8.0
int BWmap[2*(int)(contRad/contBin)+2][2*(int)(contRad/contBin)+2];
int BWmapxy[2*(int)(contRad/contBin)+2][2*(int)(contRad/contBin)+2];
int BaCmap[2*(int)(contRad/contBin)+2][2*(int)(contRad/contBin)+2];
int BaCmapxy[2*(int)(contRad/contBin)+2][2*(int)(contRad/contBin)+2];
int BC8map[2*(int)(contRad/contBin)+2][2*(int)(contRad/contBin)+2];
int BC8mapxy[2*(int)(contRad/contBin)+2][2*(int)(contRad/contBin)+2];
int BBrmap[2*(int)(contRad/contBin)+2][2*(int)(contRad/contBin)+2];
int BBrmapxy[2*(int)(contRad/contBin)+2][2*(int)(contRad/contBin)+2];
int BBcommap[2*(int)(contRad/contBin)+2][2*(int)(contRad/contBin)+2];
int BBcommapxy[2*(int)(contRad/contBin)+2][2*(int)(contRad/contBin)+2];

double *sqDiv[4]; /* 0 = x,y,x
		     1 = x
		     2 = x,y
		     3 = norm const. */

/*
 * orientational distributions
 */
#define nODbin 0.05
#define nODrange 12.0
#define nODbinz 1.0
#define cradMax 15.0
int nOD[5][(int)(nODrange/nODbinz)*2][(int)(2.0/nODbin)];
int nODnorm[5][(int)(nODrange/nODbinz)*2];
#define BCDODbin 0.05  // bin width cos(BCD vector vs z axis)
#define BCDODrange 4.0
#define BCDODbinz 0.5
#define sODbin 0.05
#define sODrange 12.0
#define sODbinz 1.0
#define sOOPrange 20.0
#define sOOPbinz 0.2
int BCDOD[(int)(BCDODrange/BCDODbinz)*2][(int)(2.0/BCDODbin)]; //[z bin][cosine]
int BCDODnorm[(int)(BCDODrange/BCDODbinz)*2];
int BCDODt[(int)(2.0/BCDODbin)];
int BCDODtnorm;
int sOD[5][(int)(sODrange/sODbinz)*2][(int)(2.0/sODbin)];
int sODnorm[5][(int)(sODrange/sODbinz)*2];
double sOOP[5][(int)(sOOPrange/sOOPbinz)*2]; // orientational order parameter
int sOOPnorm[5][(int)(sOOPrange/sOOPbinz)*2];
int sODt[5][(int)(2.0/sODbin)];
int sODtnorm[5]; // 0:water dip,1:Br-aC,2:C7-C8,3:C1-C8,4:C4-C5
#define maxBrComr 12.0
double cosgz;   // cosine of guest C8-->Br and BCDz vectors
int BrOgOD[(int)(2.0/sODbin)],BrOgODnorm;
int aCBCDr[(int)(2.0*maxBrComr/0.5)];

double radGyr[(int)(sOOPrange/sOOPbinz)*2];
int radGyrn[(int)(sOOPrange/sOOPbinz)*2];
/*
 * b-CD pore guests
 */
#define maxPoreWat 11
#define maxPoreB 4
#define pwtau 2.0
#define mnpwtau 0.0
#define pbtau 2.0
#define mnpbtau 0.0
int poreWat,poreW[maxPoreWat+1];
int probPoreWat[maxPoreWat+2];
int poreBrO[4],poreB[4][maxPoreB+1];
int probPoreBrO[4][maxPoreB+2];
tripfz *poreWdat[maxPoreWat+1];
double topH,botH,poreRad;
int poreCoM,poreBr,poreC8,poreaC;
int gst[3],ph2o;
double raC[3],rCoM[3],rC8[3];
/*
 * orientational dynamics
 */
tripd *comhist, *zhist; 
tripfzp **watDipVec,**BrOVec[3];
tripfz **watHHVec;
double watTCFbins[6],BrOTCFbins[6];
double *watHHTCF[3],*watHHP2TCF[3],*watDipTCF[3],*BrOTCF[3][3];
int *watHHTCFn[3],*watDipTCFn[3],*BrOTCFn[3][3];
double *BCDOTCF[2];
int *wTCF,*wTCFn;
/*
 * Solvent solute switching parameters
 */
double swSoMax;
double swSoMin;
double swSocoef[4];
double pswSocoef[3];
int KillFinger;
double EXFIELD;
int windowOn;
double center_w,width_w,power_w,pot_w;
double V_window;
double zGH,wGH; // host-guest distance, window
/*
 * Switching charges variables and types
 */
typedef struct	{
 	double	qa;		/* Charge */
 	double	qr;		/* Rate */
}	qswitch;

/*
 * EVB 
 */
int evbFlag;
tripd *fevb; /*To store temporary forces for EVB*/
double Qa, Qb, Qc, Qd;
double Dmors, amors, reqmors;
double Zion, rstar, EA;
double bid, epsid, sigid, nid;
double mua, mub, Qbeta;
double C1sq, RC1, RC2, cosTheta;
double VINT_WI1, VINT_WI2, VINT_BI1, VINT_BI2;
double VINT_WD1, VINT_WD2, VINT_BD1, VINT_BD2;
double VCL_WI1, VCL_WI2, VCL_BI1, VCL_BI2;
double VCL_WD1, VCL_WD2, VCL_BD1, VCL_BD2;
double VINT_CDI1, VINT_CDI2, VCL_CDI1, VCL_CDI2;
double VINT_CDD1, VINT_CDD2, VCL_CDD1, VCL_CDD2;
double VINT_oCDI1, VINT_oCDI2, VCL_oCDI1, VCL_oCDI2; // SN2 - bCD hydroxyl interactions
double VINT_oCDD1, VINT_oCDD2, VCL_oCDD1, VCL_oCDD2; // "
double VSN2_BCDoh; // SN2 - bCD hydroxyl total interaction
double VSN2_BCD; // SN2 - bCD total interaction
double VINT_EVB; // reactive system - solvent interaction
double s_coord;
double V_teth,EKVIB,VVIB,rVIB; // vibrational energy and radius of the sn2 vibration
ljcon	wslj[3][3];  /* lj - coul for water-solute (ICN)     */
ljcon	bslj[BrOs][3];  /* lj - coul for Br-octane-solute (ICN)     */
ljcon	cdslj[21][3];  /* lj - coul for b-CD - solute (ICN)     */
double cosccz,rBCDsol,rBCDz;
double sn2z,V_sn2w;
int Nshel[2]; // waters in solvation shell of Chlorine molecules
/*
 * Switching charges variables
 */
double ksw,kswcon,oe_k2; //jjk not sure what this is
double soluteD[3][400];		/*solvent and solute density profiles*/

/* EVB window parameters */
double parabA,parabB;
double Ecenter_w;
double Ewidth_w;
double Epower_w;
double Epot_w;
double solBCD_c,solBCD_w; // solute CoM-BCD CoM center of window, width of window
double osysCmz;
double V_wBSN2; // bCD-SN2 CoM tethering potential
/*reaction coordinate constraints*/
double kRC, x1RC, x2RC, VRC, FRC, NW;
double Vbias, Abias, Bbias, alpha_bias, beta_bias;
double V_solB; // solute-BCD biasing
/*dynamic reactive  flux calculations*/
int fluxQ;/*flag to signal to sysForce to calculate flux correlation function*/
int corflux[1001];/*array to save the correlation function*/
double totc1sq[1001]; /* sum of c1^2 */
double totc[1001];
double cprod[1001];
double endC;
int flux_t;/*time point for flux calculations*/ 
double init_dir;/*initial direction of motion from the TS*/
double init_RC;/*exact location of TS, should be near zero*/

//int tagwat;
tripd btop,bbot;
