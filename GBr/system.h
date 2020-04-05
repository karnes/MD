int tdpoint;
int nGLY,nTHA,nBr,nCl2,nTS;
#define GLYsites 14
ljcon GLYlj[GLYsites+1][GLYsites+1];//+1 is for Br-
//indices of pairs
//GLY force field
#define nGLYstr 13
int GLYstr[nGLYstr][3];
#define nGLYbend 21
int GLYbend[nGLYbend][4];
#define nGLYtors 27
int GLYtors[nGLYtors][5];
#define nGLY15 30
int GLY15[nGLY15][2];
#define g14qscale (1.0/1.2)
#define g14ljscale (1.0/2.0)
double gkstr[4],greq[4];
double gkbend[6],geqbend[6];
double gtors[7][3];

//THA+ force field
#define THAsites 77
ljcon THAlj[THAsites+1][THAsites+1];//+1 is for Br-
#define nTHAstr 19
int THAstr[nTHAstr][3];
#define nTHAbend1 36 //intra
#define nTHAbend2 6 //inter
int THAbend1[nTHAbend1][4],THAbend2[nTHAbend2][4];
#define nTHAtors1 37 //intra
#define nTHAtors2 36 //inter
int THAtors1[nTHAtors1][5],THAtors2[nTHAtors2][5];
#define t14qscale (1.0/1.2)
#define t14ljscale (1.0/2.0)
#define nTHA15 2516
int THA15[nTHA15][2];
double tkstr[3],treq[3];
double tkbend[6],teqbend[6];
double ttors[7][3];

//GLY-solute and THA-solute
//6 = M,M0,Cl(Cl2),Cl(TScenter),Cl(TS),Br(TS)
ljcon solGLYlj[GLYsites][6];
ljcon solTHAlj[THAsites][6];
double Clkstr,Clreq;
double TSkstr[2],TSreq[2];
double TSkbend,TSeqbend;
double rClM;

double syscomz; //system CoM on z-axis
double BrZForce;

#define grdt 40
int gGOX[3],gGCX[3],gGHX[3];
double gGOXmax[3],gGCXmax[3],gGHXmax[3];
//ensure no double counting of g_GLY-X(r)
int gflag[GLYsites];

double cosTHAX; //cosine of angle between +z and r_THA-halogen
double cosXz; //cosine of angle between halogen and z (Cl2 & TS)
double Xz; // z position of halogen
// (X to THA com vector)

tripd THAcent, NTHAvec;
double kappa;
#define kappaNorm 2.0

tripd THAcom;
double RgTHA,THANz,rXN;//rad gyration,N's z,Br-N distance
double rXCom;
double XGLYV,XGLYC,XTHAV,XTHAC;
double XGLYs[4]; // first solvation shell interaction only: C,O,H(oh),H
double fXR,fTR;

// BIASING (THA - Br)
int biasingOn;
double center_w,width_w,pot_w,power_w;
double BIAS_W,BIAS_A,BIAS_B,BIAS_C,BIAS_D,BIAS_G;
double V_bias,V_w,osysCmz;

ljcon GTlj[GLYsites][THAsites];
ljcon WTlj[3][THAsites];

double gbondE,gbendE,gtorsE,g14E,g15E;
double tbondE,tbendE,ttorsE,t14E,t15E;
double Cl2bondE,TSbondE,TSbendE;

double INTRA_GLY,INTER_GLY,GLYC;
double INTRA_THA,INTER_THA,THAC;
double GLYTHAV,GLYTHAC;
double WATERV,VWATTS;

int windowOn,XWindowOn;
double Xmin;
double center_r,window_w;
double V_window;
double V_Xwindow;

int kappaWindowOn,kappaBiasingOn;
double kap_c,kap_w; //cavity coordinate center, width
double VkapWin;

/* 
 * densitiies 
 */
tripd *GLYden[3];
int	npoints;
double	binSize;
#define grGLYtypes 16//10
/* 0: C-C
 * 1: O-H(hydroxyl)
 * 2: C-O
 * 3: C-H
 * 4: O-H
 * 5: H-H
 * 6: O-O
 * 7: C-X1 //Br
 * 8: O-X1 //Br
 * 9: H(oh)-X1 //Br
 *10: C-X2 //
 *11: O-X2 //
 *12: H(oh)-X2 //
 *10: C-X3 //
 *11: O-X3 //
 *12: H(oh)-X3 //
 */
double grGLY[grGLYtypes][300];
#define GLYdensity 0.008246
#define binRDF 0.05

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
