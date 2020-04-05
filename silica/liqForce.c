#include	<md.h>
#include        <system.h>
liqForce(pos, force)
	tripd	*pos;
	tripd	*force;
{
//printf("enter liqForce.c\n");
double Reflect(), ER = 0.;
double window(), EWin = 0.;
INTRAV_S = VNB_S = 0.;
SiForce();  
VLIQ = INTRAV_S + VNB_S;
tc++;/* global integer time for correlation function*/

// reset non-equilib and HB data if required
if(tc%HBdt==0){
	HBSitoM = HBSitoA = HBMtoSi = 0;
/*	SiHtoMO = MHtoSiO = SiDtoA = SiAtoA = SiOtoMOd = SiOtoMOa = 99.9;
	SiOMO = SiHMO = SiOMH = 0;
	MeOHd = MeOHa = 0;
	cosMa = cosMd = 0.0;*/
}

static int liqAlloc = 1;

//printf("after static declaration, before init, liqAlloc = %d\n",liqAlloc);
if(tc==0) initDataMatrices(liqAlloc);
liqAlloc = 0;
//printf("init'd matrices\n");
//printf("tc = %d, numDP_x_dataRate = %d, nCH3CN = %d, nCH3OH = %d\n",tc,numDP_x_dataRate,nCH3CN,nCH3OH);
//printf("HBdt = %d, TCFdt = %d, EPdt = %d\n",HBdt, TCFdt, EPdt);
if (nCH3OH > 0){
	INTRAV_M = VNB_M = 0.;
	CH3OHForce();
	VLIQ += INTRAV_M + VNB_M;
	VMS = 0.;
	MSiForce();
	VLIQ += VMS;
}
//printf("liqForce.c -- after CH3OHForce\n");
if (nCH3CN > 0){
	INTRAV_A = VNB_A = 0.;
	CH3CNForce();
	VLIQ += INTRAV_A + VNB_A;
	VAS = 0.;
	ASiForce();
	VLIQ += VAS;
}
//printf("liqForce.c -- after CH3CNForce\n");
if (nCH3OH > 0 && nCH3CN > 0){
	VMA = 0.;
	MAForce();
	VLIQ += VMA;
}
//printf("liqForce.c -- after forces\n");
/* add soft reflecting wall at Zwall-5 to prevent evaporating
 * liquid from entering the other side
 */
if (nSi > 0 && nCH3OH+nCH3CN > 0){
	ER = Reflect();
	VLIQ += ER;
	if (ER > 0) fprintf(stderr," ER = %f\n",ER*KCAL);
//	printf("setWindow = %d, W1 = %f, W2 = %f\n",setWindow, W1, W2);
	if (setWindow){/*modified by ilan 9/9/14*/
		EWin = window(setWindow,W1,W2);
		VLIQ += EWin;
		if (EWin > 0) fprintf(stderr," EWin = %f\n",EWin*KCAL);
	}
}
//printf("liqForce.c -- after reflect & window\n");
if (tc % EPdt == 0)
{
	EpotProf();/* Calculate the electrostatic potential at different grid points and average*/
}
//printf("liqForce.c -- after EpotProf\n");
}
#define ExCon 0.5
double
Reflect()
{
int i,j,nm;
double cmposZ,zext;
double RFIELD;
RFIELD = 0.;
nm = nCH3OH*3;
zext = zwall -5;/*reflecting wall 5A from the box edge*/

for     (i = 0; i < nm; i=i+3){
        cmposZ = 0.;
        for (j=0;j<3;j++)
                cmposZ += mass[i+j]*pos[i+j].fz;
        cmposZ /= CH3OHMass;
        if (cmposZ > zext) {
                RFIELD += ExCon*(cmposZ-zext)*(cmposZ-zext);
                for (j=0;j<3;j++)
                    force[i+j].fz -= 2*ExCon*mass[i+j]*(cmposZ-zext)/CH3OHMass;
	}
}
for	(i = 0; i < nCH3CN*3; i=i+3){
	cmposZ = 0.;
	for (j=0;j<3;j++)
		cmposZ += mass[nm+i+j]*pos[nm+i+j].fz;
        cmposZ /= CH3CNMass;
        if (cmposZ > zext) {
                RFIELD += ExCon*(cmposZ-zext)*(cmposZ-zext);
                for (j=0;j<3;j++)
                   force[nm+i+j].fz -= 2*ExCon*mass[nm+i+j]*(cmposZ-zext)/CH3CNMass;
	}
}
return(RFIELD);
}
/* restrict all the methanol molecules to a slab between w1 and w2*/
double
window(sol,w1,w2)
int sol;/*sol = 1 for methanol , 2 for acetonitrile*/
double w1, w2;
{
int i,j,nm;
double cmposZ;
double Wpot;
Wpot = 0.;
nm = nCH3OH*3;
if (sol == 1){
   for     (i = 0; i < nm; i=i+3){
        cmposZ = 0.;
        for (j=0;j<3;j++)
                cmposZ += mass[i+j]*pos[i+j].fz;
        cmposZ /= CH3OHMass;
        if (cmposZ > w2) {
                Wpot += ExCon*(cmposZ-w2)*(cmposZ-w2);
                for (j=0;j<3;j++)
                    force[i+j].fz -= 2*ExCon*mass[i+j]*(cmposZ-w2)/CH3OHMass;
	}
        if (cmposZ < w1) {
                Wpot += ExCon*(cmposZ-w1)*(cmposZ-w1);
                for (j=0;j<3;j++)
                    force[i+j].fz -= 2*ExCon*mass[i+j]*(cmposZ-w1)/CH3OHMass;
	}
   }
}
if (sol == 2){
   for     (i = 0; i < nCH3CN*3; i=i+3){
        cmposZ = 0.;
        for (j=0;j<3;j++)
                cmposZ += mass[nm+i+j]*pos[nm+i+j].fz;
        cmposZ /= CH3CNMass;
        if (cmposZ > w2) {
                Wpot += ExCon*(cmposZ-w2)*(cmposZ-w2);
                for (j=0;j<3;j++)
                    force[nm+i+j].fz -= 2*ExCon*mass[nm+i+j]*(cmposZ-w2)/CH3CNMass;
	}
        if (cmposZ < w1) {
                Wpot += ExCon*(cmposZ-w1)*(cmposZ-w1);
                for (j=0;j<3;j++)
                    force[nm+i+j].fz -= 2*ExCon*mass[nm+i+j]*(cmposZ-w1)/CH3CNMass;
	}
   }
} 
return(Wpot);
}

void initDataMatrices(int liqAlloc)
{
	int i,j,k;
	int zbins = (int)(zwall / ZBinSize); /* z bin for orientational distributions */
	int HBzbins = (int)(zwall / HBZBinSize);

/*	durMeOHa = durMeOHd = 0.0;	
*/	
// temp bins for survival probability run
/*	sbins[0] = 0.0;
	sbins[1] = 6.2;
	sbins[2] = 0.0;
	sbins[3] = 6.0;
	sbins[4] = 0.0;
	sbins[5] = 7.2;
	sbins[6] = 25.0;
	sbins[7] = 34.0;  */
	sbins[0] = 0.0;
	sbins[1] = 4.5;
	sbins[2] = 4.5;
	sbins[3] = 6.0;
	sbins[4] = 0.0;
	sbins[5] = 7.2;
	sbins[6] = 25.0;
	sbins[7] = 34.0; 

	int cosbins = (int)(2.0 / CosBinSize);
	/*TCFs, H-Bonds, orientational distributions */
	if(liqAlloc){
	if((MOHvec = (tripfz **)calloc(nCH3OH, sizeof(tripfz *))) == NULL
	|| (MOCvec = (tripfz **)calloc(nCH3OH, sizeof(tripfz *))) == NULL
	|| (ACNvec = (tripfz **)calloc(nCH3CN, sizeof(tripfz *))) == NULL
	|| (SiHB  = (int **)calloc(2 * nSi, sizeof(int *))) == NULL
	|| (SiHB2 = (int **)calloc(2 * nSi, sizeof(int *))) == NULL
	|| (odOH = (int **)calloc(zbins, sizeof(int *))) == NULL
	|| (odOC = (int **)calloc(zbins, sizeof(int *))) == NULL
	|| (odCN = (int **)calloc(zbins, sizeof(int *))) == NULL)
		ERROR((stderr,"liqForce.c: out of core\n"),exit);
	}
//	printf("liqForce.c arrays --- nSi+nCH3OH = %d, int(numDP_x_dataRate/HBdt) + 1=%d\n",nSi+nCH3OH,(int)(numDP_x_dataRate / HBdt) + 1);
//printf("before HBcount/n");	
	if(liqAlloc){
	for(i=0;i<6;i++)
	{
		if((HBcount[i] = (int *)calloc(HBzbins, sizeof(int))) == NULL) 
			ERROR((stderr,"liqForce.c: out of core\n"),exit);
	}
	}
	for(i=0;i<6;i++)
	{
		for(j=0;j<HBzbins;j++)
		{
			HBcount[i][j] = 0;
		}
	}
//printf("after HBcount\n");
	if(liqAlloc){
	for(i=0; i < zbins; i++)
	{
		if((odOH[i] = (int *)calloc(cosbins, sizeof(int))) == NULL
		|| (odOC[i] = (int *)calloc(cosbins, sizeof(int))) == NULL
		|| (odCN[i] = (int *)calloc(cosbins, sizeof(int))) == NULL)
			ERROR((stderr,"liqForce.c: out of core\n"),exit);
	}
	}
	for(i=0; i < zbins; i++)
	{
		for(j=0; j < cosbins; j++)
		{
			odOH[i][j] = odOC[i][j] = 0;
		}
	}
//printf("odOH done\n");
	for(i=0; i < zbins; i++)
	{
		for(j=0; j < cosbins; j++)
		{
			odCN[i][j] = 0;
		}
	}
	if(liqAlloc){
	for(i=0; i < NumSBins; i++)
	{
		if((sodOH[i] = (int *)calloc(cosbins, sizeof(int))) == NULL
		|| (sodOC[i] = (int *)calloc(cosbins, sizeof(int))) == NULL
		|| (sodCN[i] = (int *)calloc(cosbins, sizeof(int))) == NULL)
			ERROR((stderr,"liqForce.c: out of core\n"),exit);

	}
	}
	for(i=0; i < NumSBins; i++)
	{
		for(j=0; j < cosbins; j++)
		{
			sodOH[i][j] = sodOC[i][j] = sodCN[i][j] = 0;
		}
	}
//printf("odOC done\n");
	if(liqAlloc){
	for(i=0; i < nCH3OH; i++)
	{
		if((MOHvec[i] = (tripfz *)calloc((int)(numDP_x_dataRate / TCFdt) + 1, sizeof(tripfz))) == NULL
		|| (MOCvec[i] = (tripfz *)calloc((int)(numDP_x_dataRate / TCFdt) + 1 , sizeof(tripfz))) == NULL)
			ERROR((stderr,"liqForce.c: out of core\n"),exit);
	}
	for(i=0; i < nCH3CN; i++)
	{
		if((ACNvec[i] = (tripfz *)calloc((int)(numDP_x_dataRate / TCFdt)  + 1, sizeof(tripfz))) == NULL)
			ERROR((stderr,"liqForce.c: out of core\n"),exit);
	}
	}
	/* manually initialize arrays to avoid "optimistic" segfaults */ 
	for(i=0; i < nCH3OH; i++)
	{
		for(j=0; j < (int)(numDP_x_dataRate / TCFdt) + 1; j++)
		{
			MOHvec[i][j].fx = MOHvec[i][j].fy = MOHvec[i][j].fz = MOHvec[i][j].zpos = 0.0;
			MOCvec[i][j].fx = MOCvec[i][j].fy = MOCvec[i][j].fz = MOCvec[i][j].zpos = 0.0;
		}
	}
	for(i=0; i < nCH3CN; i++)
	{
		for(j=0; j < (int)(numDP_x_dataRate / TCFdt) + 1; j++)
		{
			ACNvec[i][j].fx = ACNvec[i][j].fy = ACNvec[i][j].fz = ACNvec[i][j].zpos = 0.0;
		}
	}
	if(liqAlloc){
	for(i=0; i < (2 * nSi); i++)
	{
		if((SiHB[i]  = (int *)calloc((int)(numDP_x_dataRate / HBdt) + 1, sizeof(int))) == NULL
		|| (SiHB2[i] = (int *)calloc((int)(numDP_x_dataRate / HBdt) + 1, sizeof(int))) == NULL)
			ERROR((stderr,"liqForce.c: out of core\n"),exit);
	}
	}
	/* initialize H-Bond matrix to -1  **/
//	printf("initialize H-Bond matrix to -1's\n");
	for(i=0; i < (2 * nSi); i++)
	{
		for(j = 0; j < (int)(numDP_x_dataRate / HBdt) + 1; j++)
		{
			SiHB[i][j]  = -1;
			SiHB2[i][j] = -1;
		}
	}
//	printf("liqForce.c -- 999'd HBinitDataMatrices\n");
	doubleSiA = doubleSiD = 0;

	for(i=0; i < 4; i++)
	{
		for(j=0; j < NZ; j++)
		{
			pGrid[i][j] = pGrid2[i][j] = 0.0;
		}
	}
//printf("liqForce - init'd matrices.\n");

}
