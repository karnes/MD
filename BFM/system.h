int	nCFM;			/* number of CFM molecules */
int	tdpoint;
double	INTRAV;			/* intramolecular potential energy of CFM*/
double	CFMNB;			/* intermolecular potential energy of CFM*/
double	VCOUL;			/* Coulomb potential energy of solute-CFM*/
ljcon	CFMlj[5][5];		/* lj - coul for CFM     */

double Req[5];		/* Equilibrium bond length for chloroform*/
double ks[5];		/* Streatching force constants for chloroform*/
double Beq[5][5];	/* Equilibrium bond angles for chloroform*/
double kb[5][5];	/* Bending force constants for chloroform*/
double QCFM[5];		/* charges on chloroform atoms*/

//#define CFMMass 119.38
//#define CFMMass 252.731
double CFMMass;
double CFMden;
/* On the fly Rfd and diff calculations */
#define binRDF 0.05
double CFMRdf[7][600];	/*chloroform atomic radial distributions*/
//#define binD 0.1
#define FS1 3.7
#define FS2 6.0 

#define TCFdt 10
tripd **dipVec; // dipole vector for orientational TCF
tripd **c3Vec; // normal to C3 axis vector, for orienational TCF
tripd **chVec; // C-H axis vector, for orienational TCF
int **pstack; // polar stacking vector
int **pstack2; // polar stacking vector
//int printStx[2000]; //indices to print (for testing)
// Kirkwood correlation factor
#define gKrad 22.0
#define gKbin 0.2
double gK[(int)(gKrad/gKbin)];
double gKn[(int)(gKrad/gKbin)];
/*
 *    angluar dependency plots
 */
double CHmax2,CHmin2;
double rCCmax,rCCmin;
double avgDipDip, avgDipDip2, avgDipDipn;
double avgDipDip2n, totInHCC;
#define contRad 10.0 // max radius to look for gCOM-COM
//#define contRbin 0.05 // size (Å) of bin
#define dAng 2.5
#define dAngsm 1.0  // sm = small
#define planAng 0.25
double planCos;
//#define contAbin 0.5
#define contBin 0.1
//int dipMap[5][2*(int)(contRad/contBin)+1][2*(int)(contRad/contBin)+1];
double dipMap[5][2*(int)(contRad/contBin)+1][2*(int)(contRad/contBin)+1];
double dipMapsm[5][2*(int)(contRad/contBin)+1][2*(int)(contRad/contBin)+1];
int CHmap[2*(int)(contRad/contBin)+1][2*(int)(contRad/contBin)+1];
int CClmap[2*(int)(contRad/contBin)+1][2*(int)(contRad/contBin)+1];
int CCmap[2*(int)(contRad/contBin)+1][2*(int)(contRad/contBin)+1];
#define MAXPOLAR 10000
//double CHPolar[MAXPOLAR][2],CClPolar[MAXPOLAR][2],dipPolar[5][MAXPOLAR][2]; // [0] is radius, [1] is angle in degrees
//int dipPolarCount[5],CHPolarCount,CClPolarCount;
//
int **stkLife[4];
// Added 24-June-2017: tabulate actual dipole-dipole angle of 'stacked' CFM
#define dStkCos 0.01
int stackCos[(int)(2.0/dStkCos)];
#define cosApollo 0.8660254
int apollo,nApollo;
