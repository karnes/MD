int tdpoint;
int nBr,nCl2;
ljcon clj[3];// 3 Cl X 1 Br = 3 combinations

double Clkstr,Clreq;

double syscomz; //system CoM on z-axis

double Cl2BondE;

double INTER_X,X_C;
double XXdist;
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
