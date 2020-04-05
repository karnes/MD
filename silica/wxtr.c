#include	<md.h>
#include	<system.h>
#include	<math.h>
#define ZS 4.5 /* top location of SiOH H atoms*/
#define MeDen 0.014882
#define AcDen 0.011531

wxtr(fp, initQ,pFreq)
	FILE	*fp;
	int	initQ;
	int pFreq;
{
int	i,j,k,l,index;	/* dummy indexes			*/
double z, r, gfactor;
double cmZ;
tripd	n_factor;	/* normalization factor to convert to **
			**  reduced densities */
int maxz;    /*  z of highest particle  */


int *SHBTCF[6]; /* SiOH-related H-bond lifetimes TCF*/
int *CHBTCF[6]; /* SiOH is: 0=donor(to M) 1=accept(from M) 2=donor(to A) */
int *normSHB[6];
int *normCHB[6];
int *vec;
double *SiDonor;
double *SiAnyM;
double *MtoSi;
double *SitoM;
double *SitoA;
double *SiTot;
double *SiDonAcc;
int *snormOH[NumSBins], *snormOC[NumSBins], *snormCN[NumSBins];
double *sTCFOH[NumSBins], *sTCFOC[NumSBins], *sTCFCN[NumSBins];
int *surM[NumSBins], *normSurM[NumSBins];
int *surA[NumSBins], *normSurA[NumSBins];

//printf("wxtr.c\n");
if (initQ == 1) {
	binSize = 0.2;
	npoints = (int) (zwall/binSize);
	if ((CMden[0]=(double *)calloc(npoints,sizeof(double))) == NULL
	|| (CMden[1]=(double *)calloc(npoints,sizeof(double))) == NULL
	|| (Qden[0]=(double *)calloc(npoints,sizeof(double))) == NULL
	|| (Qden[1]=(double *)calloc(npoints,sizeof(double))) == NULL
	|| (Mdensity[0]=(tripd *)calloc(npoints,sizeof(tripd))) == NULL
	|| (Mdensity[1]=(tripd *)calloc(npoints,sizeof(tripd))) == NULL
	|| (Mdensity[2]=(tripd *)calloc(npoints,sizeof(tripd))) == NULL
	|| (Adensity[0]=(tripd *)calloc(npoints,sizeof(tripd))) == NULL
	|| (Adensity[1]=(tripd *)calloc(npoints,sizeof(tripd))) == NULL
	|| (Adensity[2]=(tripd *)calloc(npoints,sizeof(tripd))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	/* find new bin size so we have whole bins*/
	binSize = zwall/npoints; /* exact bin size*/
	for	(i = 0; i < npoints; i++){
	     CMden[0][i] = CMden[1][i] = 0.;
	     Qden[0][i] = Qden[1][i] = 0.;/*added by ilan 8/22*/
	     Mdensity[0][i].fx=Mdensity[0][i].fy=Mdensity[0][i].fz = 0.;
	     Mdensity[1][i].fx=Mdensity[1][i].fy=Mdensity[1][i].fz = 0.;
	     Mdensity[2][i].fx=Mdensity[2][i].fy=Mdensity[2][i].fz = 0.;
	     Adensity[0][i].fx=Adensity[0][i].fy=Adensity[0][i].fz = 0.;
	     Adensity[1][i].fx=Adensity[1][i].fy=Adensity[1][i].fz = 0.;
	     Adensity[2][i].fx=Adensity[2][i].fy=Adensity[2][i].fz = 0.;
        }
	for (i = 0; i<6; i++)
		for (j = 0; j<3; j++)
			for (k = 0; k < 400; k++){
				SiLiqRdf[i][j][k] = 0.;
				ccRDF[k] = 0.;
			}
				
	tdpoint = 0;
}
	tdpoint++;	/* total number of data points	*/
	maxz = 0;
	for(i=0; i < natoms; i++)
	{
		if((int)(pos[i].fz + 1.0) > maxz) maxz = (int)(pos[i].fz + 1.0);
	}
	maxz++;
	/***	Calculating the CH3OH density profile	***/
	for	(i = 0; i < 3*nCH3OH; i = i + 3){
		cmZ = 0;
		for (j = 0; j<3; j++){
			index = (int) ((pos[i+j].fx + xwall+5)/(binSize));
			Mdensity[j][index].fx += 1.;
			index = (int) ((pos[i+j].fy + ywall+5)/(binSize));
			Mdensity[j][index].fy += 1.;
			index = (int) ((pos[i+j].fz)/(binSize));
			Mdensity[j][index].fz += 1.;
			cmZ += mass[i+j]*pos[i+j].fz;
			Qden[0][index] += E2*gridQ2[j];/*added by ilan 8/22*/
		}
		cmZ /= CH3OHMass;/*methanol mass is defined in system.h*/
		index = (int) (cmZ/binSize);
		CMden[0][index] += 1.;
	}
	/***	Calculating the CH3CN density profile	***/
	for	(i = 0; i < 3*nCH3CN; i = i + 3){
		cmZ = 0;
		for (j = 0; j<3; j++){
			index = (int) ((pos[3*nCH3OH+i+j].fx + xwall+5)/(binSize));
			Adensity[j][index].fx += 1.;
			index = (int) ((pos[3*nCH3OH+i+j].fy + ywall+5)/(binSize));
			Adensity[j][index].fy += 1.;
			index = (int) ((pos[3*nCH3OH+i+j].fz)/(binSize));
			Adensity[j][index].fz += 1.;
			cmZ += mass[3*nCH3OH+i+j]*pos[3*nCH3OH+i+j].fz;
			Qden[1][index] += E2*gridQ2[j+3];/*added by ilan 8/22*/
		}
		cmZ /= CH3CNMass;/*defined in system.h*/
		index = (int) (cmZ/binSize);
		CMden[1][index] += 1.;
	}
/*********************************************************
*** 	    	print distributions                    ***
**********************************************************/ 
	if((tdpoint-1) % pFreq == 0 && tdpoint > 1 ){
/**********************************************************
  	SOLVENT DETAILS
***********************************************************/
	/*** print solvent details for external calculations ***/
		fprintf(fp,"nSi\tnCH3OH\tnCH3CN\n");
		fprintf(fp,"%d\t%d\t%d\n\n", nSi, nCH3OH, nCH3CN);
/**********************************************************
  	ELECTROSTATIC POTENTIAL PROFILE
***********************************************************/
	/*** electrostatic potential calculations ***/
		fprintf(fp,"\n Electrostatic Potential Profile\n");
		fprintf(fp,"Z\t  MeOH\t SigM\t ACN\t SigA\t Si\t SigSi \tTot \tSigTot\n");
		double sigma;
		for(i=0; i<NZ; i++)
		{
			fprintf(fp, "%4.2f\t", 2.6+i*0.25);
			for(j=0; j<4; j++)
			{
				pGrid[j][i] /= (numDP_x_dataRate / EPdt) + 1;
				pGrid2[j][i] /= (numDP_x_dataRate / EPdt) + 1;
				sigma = sqrt(pGrid2[j][i] - pGrid[j][i] * pGrid[j][i]);
				fprintf(fp,"% f\t % f\t",EVOLT*pGrid[j][i],EVOLT*sigma);
			}
			fprintf(fp,"\n");
		}	
		
/**********************************************************
  	CENTER OF MASS DENSITIES, RADIAL DISTRIBUTIONS
***********************************************************/
		fprintf(fp,"\nDensity Profiles\n");
		fprintf(fp,"M(O) M(H) M(M) M(c) A(C) A(N) A(M) A(c) QM QA\n");
		n_factor.fx = 1./(2*ywall*(zwall-ZS)*binSize*tdpoint); 
		n_factor.fy = 1./(2*xwall*(zwall-ZS)*binSize*tdpoint); 
		n_factor.fz = 1./(4*xwall*ywall*binSize*tdpoint); 
		for	(i = 0; i < npoints; i++){
			z =  (i+0.5)*binSize;
			  fprintf(fp,"%f %f %f %f %f %f %f %f %f %f %f\n",z,
				Mdensity[0][i].fz*n_factor.fz,
				Mdensity[1][i].fz*n_factor.fz,
				Mdensity[2][i].fz*n_factor.fz,
				CMden[0][i],
				Adensity[0][i].fz*n_factor.fz,
				Adensity[1][i].fz*n_factor.fz,
				Adensity[2][i].fz*n_factor.fz,
				CMden[1][i],
				Qden[0][i]*n_factor.fz,
				Qden[1][i]*n_factor.fz);
		}
		if (nCH3OH > 0){
			fprintf(fp,"CH3OH/SiOH radial distributions\n");
			fprintf(fp,"r      OO    OS    OH    HO    HS    HH    MO    MS    MH  CH3CH3\n");
			gfactor = (tdpoint-1)*dataRatex*4*PI*binRDF*binRDF*binRDF*MeDen*nSi;
			for (k = 0; k< 400; k++){
				fprintf(fp,"%-6.2f",k*binRDF);
                		for (i=0;i<3;i++)
                		    for (j=0;j<3;j++)
                  			fprintf(fp,"%-6.3f",SiLiqRdf[i][j][k]/(gfactor*(1.0/3+k*(k+1))));
                  			fprintf(fp,"%-6.3f",ccRDF[k]/(gfactor*(1.0/3+k*(k+1))));
				fprintf(fp,"\n");
			}
		}
		if (nCH3CN > 0){
			fprintf(fp,"CH3CN/SiOH radial distributions\n");
			fprintf(fp,"r      CO    CS    CH    NO    NS    NH    MO    MS    MH\n");
			gfactor = (tdpoint-1)*dataRatex*4*PI*binRDF*binRDF*binRDF*AcDen*nSi;
			for (k = 0; k< 400; k++){
				fprintf(fp,"%-6.2f",k*binRDF);
                		for (i=0;i<3;i++)
                		    for (j=0;j<3;j++)
                  			fprintf(fp,"%-6.2f",SiLiqRdf[i+3][j][k]/(gfactor*(1.0/3+k*(k+1))));
				fprintf(fp,"\n");
			}
		}
//	  /*original end of "if((tdpoint-1) % pFreq... loop*/
/*****************************************************
  	ORIENTATIONAL DISTRIBUTIONS
******************************************************/
	/***  print orientational distributions  ***/
	fprintf(fp,"\nORIENTATIONAL DISTRIBUTIONS     (timestep = %d)\n",tc);
	fprintf(fp,"cosine of angle from +z axis\n");
	fprintf(fp,"z \tcos \t#(O->H)\t#(O->C)\t#(C->N)\n");
	for(j=0; j < (int)((float)maxz/ZBinSize); j++)
	{
		for(i=0; i < (int)(2.0/CosBinSize); i++){
			fprintf(fp,"%4.1f \t%-3.2f \t%d \t %d \t %d\n",(double)(j)*ZBinSize,
			(double)(i)*(CosBinSize) - 1.0, odOH[j][i],odOC[j][i],odCN[j][i]);
		}
	}
	fprintf(fp,"\nORIENTATIONAL DISTRIBUTIONS   'smart bins' \n");
	fprintf(fp,"cosine of angle from +z axis\n");
	for(j=0; j < NumSBins; j++)
	{
		fprintf(fp, "\nbin %d: %3.1f to %3.1f A\n",j,sbins[2*j],sbins[2*j+1]);
		fprintf(fp,"bin \t cos\t#(O->H)\t#(O->C)\t#(C->N)\n");
		for(i=0; i < (int)(2.0/CosBinSize); i++){
			fprintf(fp,"%d \t %-3.2f \t%d \t %d \t %d\n",j,
			(double)(i)*(CosBinSize) - 1.0, sodOH[j][i],sodOC[j][i],sodCN[j][i]);
		}
	}
/********************************************************
  	ORIENTATIONAL TIME CORRELATION FUNCTIONS
*********************************************************/
	/*** setup matrices for calculating TCFs **/
//	initTCFmatrices();
	for(i=0; i < NumSBins; i++)
	{
		if((sTCFOH[i] = (double *)calloc((int)(numDP_x_dataRate / TCFdt) + 1, sizeof(double))) == NULL
		|| (sTCFOC[i] = (double *)calloc((int)(numDP_x_dataRate / TCFdt) + 1, sizeof(double))) == NULL
		|| (sTCFCN[i] = (double *)calloc((int)(numDP_x_dataRate / TCFdt) + 1, sizeof(double))) == NULL
		|| (snormOH[i] = (int *)calloc((int)(numDP_x_dataRate / TCFdt) + 1, sizeof(int))) == NULL
		|| (snormOC[i] = (int *)calloc((int)(numDP_x_dataRate / TCFdt) + 1, sizeof(int))) == NULL
		|| (snormCN[i] = (int *)calloc((int)(numDP_x_dataRate / TCFdt) + 1, sizeof(int))) == NULL)
			ERROR((stderr,"wxtr: out of core\n"), exit);
	}
	for(i=0; i< NumSBins; i++)
	{
		for(j=0; j < (int)(numDP_x_dataRate / TCFdt) + 1; j++)
		{
			snormOH[i][j] = snormOC[i][j] = 0;
			snormCN[i][j] = 0;
			sTCFOH[i][j] = sTCFOC[i][j] = 0.0;
			sTCFCN[i][j] = 0.0;
		}
	}
	/*** calculate orientational time correlation functions ***/
	tcfs(sTCFOH,sTCFOC,sTCFCN,snormOH,snormOC,snormCN);
//	 printf("wxtr: back from the tcf()\n");
	/*** print TCFs  ****/
	fprintf(fp,"\nORIENTATIONAL TIME CORRELATION FUNCTIONS - 'smart bins'\n");
	for(i=0; i < NumSBins; i++)
	{
		fprintf(fp, "\nbin %d: %3.1f to %3.1f A\n",i,sbins[2*i],sbins[2*i+1]);
		fprintf(fp,"bin \ttau(fs)\tcos(O->H) \tcos(O->C) \tcos(C->N)\n");
		for(j=0; j < (int)(tc / TCFdt) -1; j++)
		{
			fprintf(fp,"%d\t %d\t %-4.3f\t %-4.3f\t %-4.3f\n",i,(int)((double)(j * TCFdt) * 0.5),sTCFOH[i][j],sTCFOC[i][j],sTCFCN[i][j]);
		}
	}
//	printf("wxtr: Finished printing orientational TCFs\n");
/***************************************************
  	HYDROGEN BONDING
****************************************************/
	/* calculate probability of Si-related H-bonds existing */
	if((MtoSi = (double *)calloc((int)(numDP_x_dataRate / HBdt) + 1, sizeof(double))) == NULL
	|| (SitoM = (double *)calloc((int)(numDP_x_dataRate / HBdt) + 1, sizeof(double))) == NULL
	|| (SitoA = (double *)calloc((int)(numDP_x_dataRate / HBdt) + 1, sizeof(double))) == NULL
	|| (SiAnyM = (double *)calloc((int)(numDP_x_dataRate / HBdt) + 1, sizeof(double))) == NULL
	|| (SiDonor = (double *)calloc((int)(numDP_x_dataRate / HBdt) + 1, sizeof(double))) == NULL
	|| (SiTot = (double *)calloc((int)(numDP_x_dataRate / HBdt) + 1, sizeof(double))) == NULL
	|| (SiDonAcc = (double *)calloc((int)(numDP_x_dataRate / HBdt) + 1, sizeof(double))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"),exit);
	for(i=0; i < (int)(numDP_x_dataRate / TCFdt) + 1; i++)
	{
		SiDonor[i] = SiAnyM[i] = SiDonAcc[i] = 0.0;
		SitoM[i] = SitoA[i] = MtoSi[i] = 0.0;
	}
	hbondProbs(MtoSi,SitoM,SitoA,SiAnyM,SiDonor,SiTot,SiDonAcc);
//	printf("wxtr: done with H bond probablities\n");
	for(j=0; j<6; j++)
	{
		if((SHBTCF[j]  = (int *)calloc((int)(tc / HBdt) + 1, sizeof(int))) == NULL
		|| (CHBTCF[j]  = (int *)calloc((int)(tc / HBdt) + 1, sizeof(int))) == NULL
		|| (normCHB[j]  = (int *)calloc((int)(tc / HBdt) + 1, sizeof(int))) == NULL
		|| (normSHB[j]  = (int *)calloc((int)(tc / HBdt) + 1, sizeof(int))) == NULL
		|| (    vec     = (int *)calloc((int)(tc / HBdt) + 1, sizeof(int))) == NULL)
			ERROR((stderr, "wxtr: out of core\n"), exit);
	}
	for(j=0;j<6;j++)
	{
		for(i=0; i < (int)(tc / HBdt) + 1; i++)
		{
			normSHB[j][i] = SHBTCF[j][i] = 0;
			normCHB[j][i] = CHBTCF[j][i] = 0;
		}
	}
//	printf("wxtr: initialized H-bond matrix \n");
	makeSHBTCF(SHBTCF,normSHB);
//	printf("wxtr: finished making S_HB matrices\n");
	makeCHBTCF(CHBTCF,normCHB,vec);
//	printf("wxtr: finished making C_HB matrices\n");
	/* print hydrogen bonding TCFs */
	fprintf(fp,"\nHYDROGEN BOND LIFETIME CORRELATIONS\n");
	if(5 !=0)//do always. fix this jjk
	{
		fprintf(fp,"t(fs) \tC_SM \tC_MS \tC_SMMS \tC_SA \tC_SMA \t C_AllSi \tS_SM \tS_MS \tS_SMMS \tS_SA \tS_SMA \t S_AllSi \n");
		for(i=0; i < (int)(tc / HBdt); i++)
		{
			fprintf(fp,"%d \t", (int)((i)*HBdt*0.5));
			for(j=0;j<6;j++)
			{
				fprintf(fp,"%4.3f \t",(double)(CHBTCF[j][i])/(double)(normCHB[j][i]));
			}
			for(j=0;j<6;j++)
			{
				fprintf(fp,"%4.3f \t",(double)(SHBTCF[j][i])/(double)(normSHB[j][i]));
			}
			fprintf(fp,"\n");
		}
//		fprintf(fp,"\n%d total bonds detected\n",normSHB);
		fprintf(fp,"Text line for axgxtr.c detected CHB_SM[0]:%d,CHB_SM[500]:%d,CHB_SA:%dSHB_SM:%d,SHB_MS:%d,SHB_SA:%d\n",normCHB[0][0],normCHB[0][50],normCHB[0][75],normCHB[0][100],normSHB[1][0],normSHB[3][0]);
	}
	else
	{
		fprintf(fp, "No hydrogen bonds detected\n\n");
	}
//	printf("wxtr: printed h-bonds\n");
	/* print silanol H-bond probabilities */
	fprintf(fp,"\nHYDROGEN BONDING: probability of SiOH in H-bond\n");
	fprintf(fp,"t(fs)\tSitoM\tSitoA\tMtoSi\tSiDon\tSianyM\tSiD&A\n");
	for(i=0; i < (int)(numDP_x_dataRate / HBdt) + 1; i++)
	{
		fprintf(fp,"%d \t%4.3f \t%4.3f \t%4.3f \t%4.3f \t%4.3f \t%4.3f \n", (int)(i*HBdt*0.5),SitoM[i], SitoA[i], MtoSi[i], SiDonor[i], SiAnyM[i], SiDonAcc[i]);
 	}
	/* print hbond count vs z bin */
	fprintf(fp,"\nHYDROGEN BONDS: count per z-bin\n");
	for(i=0; i < (int)((double)maxz/HBZBinSize); i++)
	{
		fprintf(fp,"%4.1f \t", (double)i*HBZBinSize);
		for(j=0;j<6;j++)
		{
			fprintf(fp,"%5.2f \t", ((double) HBcount[j][i]) / (((double)tc / (double)HBdt) + 1.0));
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"\nInstances of >1 H-bond to an SiOH as donor: %d as acceptor: %d\n",doubleSiD, doubleSiA);
/*************************************************
****           SURVIVABILITY PROFILE	      ****
**************************************************/
//initSurvival();
	for(i=0; i < NumSBins; i++)
	{
		if((surM[i] = (double *)calloc((int)(numDP_x_dataRate / TCFdt) + 1, sizeof(double))) == NULL
		|| (surA[i] = (double *)calloc((int)(numDP_x_dataRate / TCFdt) + 1, sizeof(double))) == NULL
		|| (normSurM[i] = (double *)calloc((int)(numDP_x_dataRate / TCFdt) + 1, sizeof(double))) == NULL
		|| (normSurA[i] = (double *)calloc((int)(numDP_x_dataRate / TCFdt) + 1, sizeof(double))) == NULL)
			ERROR((stderr,"wxtr: out of core\n"), exit);
	}
	
	for(i=0; i < NumSBins; i++)
	{
		for(j=0;j < (int)(numDP_x_dataRate / TCFdt) + 1; j++)
		{
			surM[i][j] = surA[i][j] = 0;
			normSurM[i][j] = normSurA[i][j] = 0;
		}
	}	
popSurvival(surM,surA,normSurM,normSurA,vec);

fprintf(fp,"\nSURVIVAL PROBABILITY\n");
fprintf(fp,"t(fs)\tM(0)\tM(1)\tM(2)\tM(3)\tA(0)\tA(1)\tA(2)\tA(3)\tT(0)\tT(1)\tT(2)\tT(3)\n");
for(j=0;j < (int)(tc/TCFdt); j++)
{
	fprintf(fp,"%f\t",j*TCFdt*0.5);
	for(i=0; i<NumSBins; i++)
	{
		fprintf(fp,"%4.3f\t",((double)surM[i][j])/((double)normSurM[i][j]));
	}
	for(i=0; i<NumSBins; i++)
	{
		fprintf(fp,"%4.3f\t",((double)surA[i][j])/((double)normSurA[i][j]));
	}
	for(i=0; i<NumSBins; i++)
	{
		fprintf(fp,"%4.3f\t",((double)(surM[i][j]+surA[i][j]))/((double)(normSurM[i][j]+normSurA[i][j])));
	}
	fprintf(fp,"\n");
}
fprintf(fp,"\nSURVIVAL PROBABILITY: total particle appearances in bin... this is broken.\n");
fprintf(fp,"M(0)\tM(1)\tM(2)\tM(3)\tA(0)\tA(1)\tA(2)\tA(3)\tT(0)\tT(1)\tT(2)\tT(3)\n");
	for(i=0; i<NumSBins; i++)
	{
		fprintf(fp,"%d\t",333/*normSurM[i][0]*/);
	}
	for(i=0; i<NumSBins; i++)
	{
		fprintf(fp,"%d\t",333/*normSurA[i][0]*/);
	}
	for(i=0; i<NumSBins; i++)
	{
		fprintf(fp,"%d\t",333/*normSurM[i]+normSurA[i][0]*/);
	}
	fprintf(fp,"\n");
	


//printf("before free\n");
for(i=0; i<NumSBins; i++)
{
	free(sTCFOH[i]);
	free(sTCFOC[i]);
	free(sTCFCN[i]);
	free(snormOH[i]);
	free(snormOC[i]);
	free(snormCN[i]);
	free(surM[i]);
	free(surA[i]);
	free(normSurM[i]);
	free(normSurA[i]);
}
free(vec);
free(MtoSi);
free(SitoM);
free(SitoA);
free(SiAnyM);
free(SiDonor);
free(SiTot);
free(SiDonAcc);
for(i=0; i<6; i++)
{
	free(SHBTCF[i]);
	free(CHBTCF[i]);
	free(normCHB[i]);
	free(normSHB[i]);
	
}

} /*end of "(tdpoint-1) % pFreq == 0 && tdpoint > 1 ) loop */
fflush(fp);
//printf("exit wxtr\n");
}




void tcfs(double **sTCFOH,double **sTCFOC,double **sTCFCN, int **snormOH,int **snormOC,int **snormCN){
int i, j, k, tau;
double dot;
//printf("wxtr: in tcfs() function\n");
//printf("tc=%d\n",tc);
if(nCH3OH > 0){
	for(i = 0; i<nCH3OH; i++){
		for(j = 0; j < (int)(tc / TCFdt); j++){ /*step through all times */
			for(tau = 0; tau < j; tau ++){
			//	zbin = (int)(MOHvec[i][j].zpos / TCFZbinsize);
				dot =   MOHvec[i][j].fx * MOHvec[i][j-tau].fx +
					MOHvec[i][j].fy * MOHvec[i][j-tau].fy +
					MOHvec[i][j].fz * MOHvec[i][j-tau].fz; 
		//		printf("dot=%f\n",dot);
				for(k=0; k < NumSBins; k++)
				{
					if(MOHvec[i][j].zpos > sbins[2*k] && MOHvec[i][j].zpos < sbins[2*k+1])
					{
						sTCFOH[k][tau] +=dot;
						snormOH[k][tau]++;
					}
				}
				dot =   MOCvec[i][j].fx * MOCvec[i][j-tau].fx +
					MOCvec[i][j].fy * MOCvec[i][j-tau].fy +
					MOCvec[i][j].fz * MOCvec[i][j-tau].fz; 
				for(k=0; k < NumSBins; k++)
				{
					if(MOCvec[i][j].zpos > sbins[2*k] && MOCvec[i][j].zpos < sbins[2*k+1])
					{
						sTCFOC[k][tau] +=dot;
						snormOC[k][tau]++;
					}
				}
			}
		}
	}
//	printf("after tcf CH3OH\n");
	/* normalize and prevent division by zero  */
	for(i=0; i < NumSBins; i++)
	{
		for(j=0; j < (int)(tc / TCFdt); j++)
		{
			sTCFOH[i][j] = (snormOH[i][j] == 0)? 0.0: sTCFOH[i][j] / (float)snormOH[i][j];
			sTCFOC[i][j] = (snormOC[i][j] == 0)? 0.0: sTCFOC[i][j] / (float)snormOC[i][j];
		}
	}
}
//printf("wxtr: tcf() after CH3OH TCFs\n"); 
if(nCH3CN > 0){
//	printf("wxtr.c --- in the CH3CN part of tcfs()\n");
	for(i = 0; i < nCH3CN; i++){
		for(j = 0; j < (int)(tc / TCFdt); j++){ /*step through all times */
			for(tau = 0; tau < j; tau ++){
			//	zbin = (int)(ACNvec[i][j].zpos / TCFZbinsize);
				dot =   ACNvec[i][j].fx * ACNvec[i][j-tau].fx +
					ACNvec[i][j].fy * ACNvec[i][j-tau].fy +
					ACNvec[i][j].fz * ACNvec[i][j-tau].fz; 
				for(k=0; k < NumSBins; k++)
				{
					if(ACNvec[i][j].zpos > sbins[k*2] && ACNvec[i][j].zpos < sbins[(2*k)+1])
					{
						sTCFCN[k][tau] +=dot;
						snormCN[k][tau]++;
					}
				}
			}
		}
	}
	/* normalize and prevent division by zero */
	for(i=0; i < NumSBins; i++)
	{
		for(j=0; j < (int)(tc / TCFdt); j++)
		{
			sTCFCN[i][j] = (snormCN[i][j] == 0)? 0.0: sTCFCN[i][j] / (float)snormCN[i][j];
		}
	}
}
//printf("wxtr: tcf() after CH3CN TCFs\n"); 
}

void makeSHBTCF(int **SHBTCF,int **normSHB)
{
	int a,i,j,k,loop,start;
	int partner, duration;
	for(a=0;a<2;a++)
	{
		for(i=a*nSi; i < (a+1)*nSi; i++)  /*this loop looks for SiOH HB to MeOH*/
		{
			partner = -1;
			duration = start = 0;
			for(j=0; j < (int)(tc / HBdt); j++)
			{
				if(SiHB[i][j] == -1 && partner == -1)
				{
				}
				else if(SiHB[i][j] != -1 && SiHB[i][j] < nCH3OH && partner == -1) 
				{
					partner = SiHB[i][j];
					start = j;
					duration++;
				}
				else if((SiHB[i][j] != -1 && SiHB[i][j] < nCH3OH && SiHB[i][j] == partner) 
					|| (SiHB2[i][j] != -1 && SiHB[i][j] < nCH3OH && SiHB2[i][j] ==partner))
				{
					duration++;
				}
				else if(SiHB[i][j] != -1 && SiHB[i][j] != partner && SiHB2[i][j] != partner)
				{
					for(k=0; k < duration; k++)
					{
						SHBTCF[a][k]++;
					}
					for(k=0; k < (int)(tc / HBdt) - start; k++)
					{
						normSHB[a][k]++;
					}
					if(SiHB[i][j] < nCH3OH)
					{
						partner = SiHB[i][j];
						start = j;
						duration = 1;
					}
					else
					{
						partner = -1;
						duration = start = 0;
					}
				}
				else if(SiHB[i][j] == -1 && partner != -1)
				{	
					for(k=0; k < duration; k++)
					{	
						SHBTCF[a][k]++;
					}
					for(k=0; k < (int)(tc / HBdt) - start; k++)
					{
						normSHB[a][k]++;
					}
					partner = -1;
					duration = start = 0;
				}
			}
			if(duration > 0)
			{
				for(k=0; k < duration; k++)
				{
					SHBTCF[a][k]++;
				}
				for(k=0; k < (int)(tc / HBdt) - start; k++)
				{
					normSHB[a][k]++;
				}
				start = 0;
			}
		}
	}
	/*make total SiOH matrices*/
	for(i=0; i < (int)(tc / HBdt);i++)
	{
		SHBTCF[2][i] = SHBTCF[0][i] + SHBTCF[1][i];
		normSHB[2][i] = normSHB[0][i] + normSHB[1][i];
	}
	for(i=0; i < nSi; i++)  /*this loop looks for SiOH HB to ACN*/
	{
		a=3;
		partner = -1;
		duration = start = 0;
		for(j=0; j < (int)(tc / HBdt); j++)
		{
			if(SiHB[i][j] == -1 && partner == -1)
			{
			}
			else if(SiHB[i][j] != -1 && SiHB[i][j] >= nCH3OH && partner == -1) 
			{
				partner = SiHB[i][j];
				duration++;
				start = j;
			}
			else if((SiHB[i][j] != -1 && SiHB[i][j] >= nCH3OH && SiHB[i][j] == partner) 
				|| (SiHB2[i][j] != -1 && SiHB[i][j] >= nCH3OH && SiHB2[i][j] ==partner))
			{
				duration++;
			}
			else if(SiHB[i][j] != -1 && SiHB[i][j] != partner && SiHB2[i][j] != partner)
			{
				for(k=0; k < duration; k++)
				{
					SHBTCF[a][k]++;
				}
				for(k=0; k < (int)(tc / HBdt) - start; k++)
				{
					normSHB[a][k]++;
				}
				if(SiHB[i][j] >= nCH3OH)
				{
					partner = SiHB[i][j];
					duration = 1;
					start = j;
				}
				else
				{
					partner = -1;
					duration = start = 0;
				}
			}
			else if(SiHB[i][j] == -1 && partner != -1)
			{	
				for(k=0; k < duration; k++)
				{	
					SHBTCF[a][k]++;
				}
				partner = -1;
				for(k=0; k < (int)(tc / HBdt) - start; k++)
				{
					normSHB[a][k]++;
				}
				duration = start = 0;
			}
		}
		if(duration > 0)
		{
			for(k=0; k < duration; k++)
			{
				SHBTCF[a][k]++;
			}
			for(k=0; k < (int)(tc / HBdt) - start; k++)
			{
				normSHB[a][k]++;
			}
			start = duration = 0;
		}
	}
	for(a=4;a<6;a++)
	{
		if(nCH3OH > 0 && a==5) /*this loop gets the total S_HB for all types */
		{
			loop = 2 * nSi;
		}
		else /* SiOH as donor only */
		{
			loop = nSi;
		} 
		for(i=0; i < loop; i++)
		{
			partner = -1;
			duration = start = 0;
			for(j=0; j < (int)(tc / HBdt); j++)
			{
				if(SiHB[i][j] == -1 && partner == -1)
				{
				}
				else if(SiHB[i][j] != -1 && partner == -1) 
				{
					partner = SiHB[i][j];
					duration++;
				}
				else if((SiHB[i][j] != -1 && SiHB[i][j] == partner) 
					|| (SiHB2[i][j] != -1 && SiHB2[i][j] ==partner))
				{
					duration++;
				}
				else if(SiHB[i][j] != -1 && SiHB[i][j] != partner && SiHB2[i][j] != partner)
				{
					for(k=0; k < duration; k++)
					{
						SHBTCF[a][k]++;
					}
					for(k=0; k < (int)(tc / HBdt) - start; k++)
					{
						normSHB[a][k]++;
					}
					partner = SiHB[i][j];
					duration = 1;
					start = j;
				}
				else if(SiHB[i][j] == -1 && partner != -1)
				{	
					for(k=0; k < duration; k++)
					{	
						SHBTCF[a][k]++;
					}
					for(k=0; k < (int)(tc / HBdt) - start; k++)
					{
						normSHB[a][k]++;
					}
					partner = -1;
					duration = start = 0;
				}
				
			}
			if(duration > 0)
			{
				for(k=0; k < duration; k++)
				{
					SHBTCF[a][k]++;
				}
				for(k=0; k < (int)(tc / HBdt) - start; k++)
				{
					normSHB[a][k]++;
				}
			}
		}
	}
}

void makeCHBTCF(int **CHBTCF,int **normCHB, int *vec)
{
	int i,j,k,tau,t,l,m;
	int *partner;
	if((partner = (int *)calloc(nCH3OH + nCH3CN, sizeof(int))) == NULL)
			ERROR((stderr, "wxtr: out of core\n"), exit);
/*
//  count total StoM and StoA bonds per time step (possibly for normalization)
	for(i=0;i<nSi;i++){
		for(j=0;j<(int)(numDP_x_dataRate/HBdt)+1;j++){
		   k=SiHB[i][j];
		   if(k!=-1){
		       if(k<nCH3OH){
		          countHB[0][j]++;
		       }
		       else if(k>=nCH3OH && k < nCH3OH+nCH3CN){
			       countHB[2][j]++;
		       }
		   }
		   k=SiHB2[i][j];
		   if(k!=-1){
		       if(k<nCH3OH){
		          countHB[0][j]++;
		       }
		       else if(k>=nCH3OH && k < nCH3OH+nCH3CN){
			       countHB[2][j]++;
		       }
		   }
		}
	}

//  count total MtoS bonds per time step (possibly for normalization)
	for(i=Si;i<2*nSi;i++){
		for(j=0;j<(int)(numDP_x_dataRate/HBdt)+1;j++){
		   k=SiHB[i][j];
		   if(k!=-1){
		       if(k<nCH3OH){
		          countHB[1][j]++;
		       }
		   }
		   k=SiHB2[i][j];
		   if(k!=-1){
		       if(k<nCH3OH){
		          countHB[1][j]++;
		       }
		   }
		}
	}
*/
// make one HB vector per step, tabulate TCF.
	for(i=0;i<nSi;i++){
	   for(j=0;j<nCH3OH;j++){
	      for(k=0;k<(int)(numDP_x_dataRate/HBdt)+1;k++){
	         vec[k]=0; // wipe vector clean while populating 
		 if(SiHB[i][k]==j || SiHB2[i][k]==j){
		     vec[k]=1; //if i is HB'd to j at time k, set to 1
		 }
	      }
	      // HERE CALCULATE TCF~~~~~~~~~~~
	      for(tau=0;tau<(int)(numDP_x_dataRate/HBdt);tau++){
		 if(vec[tau]==1){
		    for(t=0;t<(int)(numDP_x_dataRate/HBdt)-tau;t++){
		       CHBTCF[0][t]+= vec[tau]*vec[tau+t];
		       normCHB[0][t]++;
		    }
		 }
	      }
	   }
	}
// make one HB vector per step, tabulate TCF.
	for(i=nSi;i<2*nSi;i++){
	   for(j=0;j<nCH3OH;j++){
	      for(k=0;k<(int)(numDP_x_dataRate/HBdt)+1;k++){
	         vec[k]=0; // wipe vector clean while populating 
		 if(SiHB[i][k]==j || SiHB2[i][k]==j){
		     vec[k]=1; //if i is HB'd to j at time k, set to 1
		 }
	      }
	      // HERE CALCULATE TCF~~~~~~~~~~~
	      for(tau=0;tau<(int)(numDP_x_dataRate/HBdt);tau++){
		 if(vec[tau]==1){
		    for(t=0;t<(int)(numDP_x_dataRate/HBdt)-tau;t++){
		       CHBTCF[1][t]+= vec[tau]*vec[tau+t];
		       normCHB[1][t]++;
		    }
		 }
	      }
	   }
	}
// make one HB vector per step, tabulate TCF.
	for(i=0;i<nSi;i++){
	   for(j=nCH3OH;j<nCH3OH+nCH3CN;j++){
	      for(k=0;k<(int)(numDP_x_dataRate/HBdt)+1;k++){
	         vec[k]=0; // wipe vector clean while populating 
		 if(SiHB[i][k]==j || SiHB2[i][k]==j){
		     vec[k]=1; //if i is HB'd to j at time k, set to 1
		 }
	      }
	      // HERE CALCULATE TCF~~~~~~~~~~~
	      for(tau=0;tau<(int)(numDP_x_dataRate/HBdt);tau++){
		 if(vec[tau]==1){
		    for(t=0;t<(int)(numDP_x_dataRate/HBdt)-tau;t++){
		       CHBTCF[2][t]+= vec[tau]*vec[tau+t];
		       normCHB[2][t]++;
		    }
		 }
	      }
	   }
	}
}

void hbondProbs(double *MtoSi,double *SitoM,double *SitoA,double *SiAnyM,double *SiDonor,double *SiTot,double *SiDonAcc)
{
	int i,j;

	for(i=0; i < nSi; i++)
	{
		for(j=0; j < (int)(tc / HBdt) + 1; j++)
		{
			if(SiHB[i][j] != -1 && SiHB[i][j] < nCH3OH)
			{
				SitoM[j] += 1.0 / (double)nSi;
				SiAnyM[j] += 1.0 / (double)nSi;
				SiDonor[j] += 1.0 / (double)nSi;
				SiTot[j] += 1.0 / (double)nSi;
			}
			else if(SiHB[i][j] != -1 && SiHB[i][j] >= nCH3OH)
			{
				SitoA[j] += 1.0 / (double)nSi;
				SiDonor[j] += 1.0 / (double)nSi;
				SiTot[j] += 1.0 / (double)nSi;
			}
			if(SiHB[i + nSi][j] != -1)
			{
				MtoSi[j] += 1.0 / (double)nSi;
				SiAnyM[j] += 1.0 / (double)nSi;
				SiTot[j] = 1.0 / (double)nSi;

			}
			if(SiHB[i][j] != -1 && SiHB[i + nSi][j] != -1)
			{
				SiDonAcc[j] += 1.0 / (double)nSi;
			}
		}
	}
}

void popSurvival(int **surM,int **surA,int **normSurM,int **normSurA,int *vec)
{
	int i,j,k,bin;
	int tau,t;

	for(i=0; i<nCH3OH; i++)
	{
           for(bin=0; bin < NumSBins; bin++)
	   {
	      for(j=0; j < (int)(tc/TCFdt)+1; j++)
	      {
	         vec[j]=0;
		 if(MOCvec[i][j].zpos > sbins[2*bin] && MOCvec[i][j].zpos < sbins[(2*bin)+1]){
			vec[j]=1;
		 }
              }
	      // HERE CALCULATE TCF~~~~~~~~~~~
	      for(tau=0;tau<(int)(numDP_x_dataRate/HBdt);tau++){
		  if(vec[tau]==1){
		     for(t=0;t<(int)(numDP_x_dataRate/HBdt)-tau;t++){
		        surM[bin][t]+= vec[tau]*vec[tau+t];
		        normSurM[bin][t]++;
		     }
		  }
	       }
	   }
	}
	for(i=0; i<nCH3CN; i++)
	{
           for(bin=0; bin < NumSBins; bin++)
	   {
	      for(j=0; j < (int)(tc/TCFdt)+1; j++)
	      {
	         vec[j]=0;
		 if(ACNvec[i][j].zpos > sbins[2*bin] && ACNvec[i][j].zpos < sbins[(2*bin)+1]){
			vec[j]=1;
		 }
	      }
	      // HERE CALCULATE TCF~~~~~~~~~~~
	      for(tau=0;tau<(int)(numDP_x_dataRate/HBdt);tau++){
		 if(vec[tau]==1){
		    for(t=0;t<(int)(numDP_x_dataRate/HBdt)-tau;t++){
		       surA[bin][t]+= vec[tau]*vec[tau+t];
		       normSurA[bin][t]++;
		    }
		 }
	      }
	    }
	 }
	

}
