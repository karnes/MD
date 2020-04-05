#include	<md.h>
#include	<system.h>
#include	<math.h>

wxtr(fp, initQ, pFreq)
	FILE	*fp;
	int	initQ;	/*initialization flag	*/
	int	pFreq;	/* frequency to print output	*/
{
void offsetBCD(),getBrOAng(),getWatAng();
int i,j,k,l,m,t,bin,wip,newWat,hbtype;
int lines, length, inp, maxGap,maxiGap,timeout, oe;
int hbp,newHB, hbpList[2000];
int nw,watList[2000];
int index,**poreRes;
int **wcrTCF,minbond,minres,minibond;
int *prdTCF,*praTCF,*prdTCFn,*praTCFn;
int *sdTCF,*saTCF,*sdTCFn,*saTCFn,*BBsdTCF,*BBsdTCFn;
int **hbRes,**hbrTCF,*HB_TCF[nBCDHB],*HB_TCFn[nBCDHB];
double *pwTCF[3],*pbTCF[3][3];
double gyr,r2,r,gfactor[2],idens[12],zavg;
double g2dwfactor,g2dbfactor;
tripd n_factor,dip,havg;
nw = natoms - nBrO*BrOs - nBCD*BCDs - nsolute;
//fprintf(stderr,"in wxtr.c  tdpoint=%d pFreq=%d\n",tdpoint,pFreq);
if (initQ == 1){
	iuBCDHB = priBCDHB = 0; // intra-unit and primary intra-BCD HB's
	pBrO.fx = pBrO.fy = pBrO.fz = 0.0;
	tdpoint = 0;
	binSize = 0.5;
	npoints = (int)(2*zwall/(binSize));
	/**********************************/
	/* allocate space and zero arrays */
	/**********************************/
	/* densities */
	for(i=0;i<2*(int)(contRad/contBin)+2;i++){
	   for(j=0;j<2*(int)(contRad/contBin)+2;j++){
	      BWmap[i][j] = 0;
	      BWmapxy[i][j] = 0;
	      BaCmap[i][j] = 0;
	      BaCmapxy[i][j] = 0;
	      BC8map[i][j] = 0;
	      BC8mapxy[i][j] = 0;
	      BBcommap[i][j] = 0;
	      BBcommapxy[i][j] = 0;
	      BBrmap[i][j] = 0;
	      BBrmapxy[i][j] = 0;
	   }
	}
	if((H2Oden = (tripd *) calloc(npoints,sizeof(tripd))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	for(i=0;i<3;i++){
	   if((BrOden[i] = (tripd *) calloc(npoints,sizeof(tripd))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	}
	if((BCDden = (tripd *) calloc(npoints,sizeof(tripd))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	for(i=0;i<npoints;i++){
		H2Oden[i].fz = 0.0;
		for(j=0;j<3;j++){
		   BrOden[j][i].fz = 0.0;
		}
		BCDden[i].fz = 0.0;
	}
	/* BCD center of mass and BCD molecular vector */
	if((comhist = (tripd *) calloc(numDPx+1,sizeof(tripd))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	if((zhist = (tripd *) calloc(numDPx+1,sizeof(tripd))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	for(i=0;i<numDPx;i++){
	   comhist[i].fx = comhist[i].fy = comhist[i].fz = -999.9;
	   zhist[i].fx = zhist[i].fy = zhist[i].fz = -999.9;
	}
	/* solvent molecular vectors */
	if((watDipVec = (tripfzp **) calloc(nw/3,sizeof(tripfzp *))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	if((watHHVec = (tripfz **) calloc(nw/3,sizeof(tripfz *))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	for(i=0;i<nw/3;i++){
	   if((watDipVec[i] = (tripfzp *) calloc(numDPx+1,sizeof(tripfzp))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	   if((watHHVec[i] = (tripfz *) calloc(numDPx+1,sizeof(tripfz))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	}
	for(j=0;j<3;j++){
	   if((BrOVec[j] = (tripfzp **) calloc(nBrO,sizeof(tripfzp *))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	   for(i=0;i<nBrO;i++){
	      if((BrOVec[j][i] = (tripfzp *) calloc(numDPx+1,sizeof(tripfzp))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	   }
	}
	for(i=0;i<nw/3;i++){
	   for(j=0;j<numDPx+1;j++){
	      watDipVec[i][j].fx = watDipVec[i][j].fy = watDipVec[i][j].fz = 0.0;
	      watDipVec[i][j].zpos = -999.9;
	      watDipVec[i][j].inp = -1;
	      watHHVec[i][j].fx  = watHHVec[i][j].fy  = watHHVec[i][j].fz  = 0.0;
	      watDipVec[i][j].zpos = -999.9;
	   }
	}
	for(i=0;i<3;i++){
	   for(j=0;j<nBrO;j++){
	      for(k=0;k<numDPx+1;k++){
		 BrOVec[i][j][k].fx = BrOVec[i][j][k].fy = BrOVec[i][j][k].fz = 0.0;
		 BrOVec[i][j][k].zpos = -999.9;
		 BrOVec[i][j][k].inp = -1;
	      }
	   }
	}
	/* pore water molecules */
	for(i=0;i<maxPoreWat+1;i++){
	   if((poreWdat[i] = (tripfz *) calloc(numDPx+1,sizeof(tripfz))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	}
	for(i=0;i<maxPoreWat+1;i++){
	   for(j=0;j<numDPx+1;j++){
	      poreWdat[i][j].fx = poreWdat[i][j].fy = poreWdat[i][j].fz = 0.0;
	      poreWdat[i][j].zpos = -999.9;
	   }
	}	
	/* square distance, BCD */
	for(i=0;i<4;i++){
	   if((sqDiv[i] = (double *) calloc(numDPx+1,sizeof(double))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	}
	for(i=0;i<4;i++){
	   for(j=0;j<numDPx+1;j++){
	      sqDiv[i][j] = 0.0;
	   }
	}
	/* BCD orientational TCFs */
	for(i=0;i<2;i++){
	   if((BCDOTCF[i] = (double *) calloc(numDPx+1,sizeof(double))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	}
	for(i=0;i<2;i++){
	   for(j=0;j<numDPx+1;j++){
	      BCDOTCF[i][j] = 0.0;
	   }
	}
	if((wTCF = (int *) calloc(numDPx+1,sizeof(int))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	if((wTCFn = (int *) calloc(numDPx+1,sizeof(int))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	for(i=0;i<numDPx+1;i++){
	   wTCF[i] = wTCFn[i] = 0;
	}
	for(i=0;i<nBCDHB;i++){
	   if((HB[i] = (int *) calloc(numDPx+1,sizeof(int))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	}
	for(i=0;i<nBCDHB;i++){
	   for(j=0;j<numDPx+1;j++){
	      HB[i][j] = -1;
	   }
	}
	/* find new bin size so we have whole bins*/
	binSize = 2*zwall/npoints; /* exact bin size*/
	/* initialize matrices to zero */
	for(i=0;i<(int)(BCDODrange/BCDODbinz)*2;i++){
	   for(j=0;j<(int)(2.0/BCDODbin);j++){
	       BCDOD[i][j]=0;
	       BCDODt[j]=0;
	   }
	   BCDODnorm[i]=0;
	}
	BCDODtnorm = 0;
	for(k=0;k<5;k++){
	   for(i=0;i<(int)(nODrange/nODbinz)*2;i++){
	      for(j=0;j<(int)(2.0/nODbin);j++){
	       nOD[k][i][j]=0;
	      }
	      nODnorm[k][i]=0;
	   }
	}
	for(k=0;k<5;k++){
	   for(i=0;i<(int)(sODrange/sODbinz)*2;i++){
	      for(j=0;j<(int)(2.0/sODbin);j++){
	       BrOgOD[j] = 0;
	       sOD[k][i][j]=0;
	       sODt[k][j]=0;
	      }
	      sODnorm[k][i]=0;
	   }
	   for(i=0;i<(int)(sOOPrange/sOOPbinz)*2;i++){
	      sOOP[k][i] = 0.0;
	      sOOPnorm[k][i] = 0;
	      radGyr[i] = 0.0;
	      radGyrn[i] = 0;
	   }
	   sODtnorm[k] = 0;
	}
	BrOgODnorm = 0;
	for(i=0;i<(int)(2.0*maxBrComr/0.5);i++){
	   aCBCDr[i] = 0;
	}
	for(i=0;i<8;i++){
   	   for(j=0;j<RDFbins;j++){
		   grSOL[i][j] = 0.0;
	   }
	}
	for(i=0;i<3;i++){
	   for(j=0;j<4;j++){
	      for(k=0;k<RDFbins;k++){
		 grCl[i][j][k] = 0.0;
	      }
	   }
	}
	for(i=0;i<maxPoreWat+2;i++){
	   probPoreWat[i]=0;
	}
	for(j=0;j<4;j++){
	   for(i=0;i<maxPoreB+2;i++){
	      probPoreBrO[j][i]=0;
	   }
	}
	/* hydrogen bond duration histograms */
	for(i=0;i<5;i++){
	   HBhistn[i] = 0;
	   for(j=0;j<HBhistlen;j++){
	      HBhist[i][j] = 0;
	   }
	}

	/* solvent bins in z-space */
        watTCFbins[0] = -6.0; //min
        watTCFbins[1] = 0.5; //max
        watTCFbins[2] = -25.0;
        watTCFbins[3] = -10.0;
        watTCFbins[4] = -zwall; 
        watTCFbins[5] = zwall;
	BrOTCFbins[0] = -0.5; //min
	BrOTCFbins[1] = 7.0; //max
	BrOTCFbins[2] = 20.0;
	BrOTCFbins[3] = 50.0;
	BrOTCFbins[4] = -zwall;
	BrOTCFbins[5] = zwall;

	comhist[0].fx = BCDcom.fx; 
	comhist[0].fy = BCDcom.fy; 
	comhist[0].fz = BCDcom.fz;
	zhist[0].fx = BCDz.fx; 
	zhist[0].fy = BCDz.fy; 
	zhist[0].fz = BCDz.fz; 

//fprintf(stderr,"wxtr.c: allocated and initialized arrays. tc = %d, tdpoint = %d\n",tc,tdpoint);
}  //   **********   END INITIALIZATION LOOP ********
/* get solvent orientational data */
if(nBrO>0){
   for(i=0;i<nBrO;i++){
      k = nw + i*BrOs;
      getBrOAng(k,1,0,1,tdpoint);//head, tail, index of OD array
      getBrOAng(k,7,8,2,tdpoint);
      getBrOAng(k,4,5,3,tdpoint);
      getBrOAng(k,1,8,4,tdpoint);
      dip.fx = dip.fy = dip.fz = 0.0;
      for(j=0;j<BrOs;j++){
	 dip.fx += pos[k+j].fx;//*mass[k+j];
	 dip.fy += pos[k+j].fy;//*mass[k+j];
	 dip.fz += pos[k+j].fz;//*mass[k+j];
      }  // use mean position, not CoM
      dip.fx/=BrOs;
      dip.fy/=BrOs;
      dip.fz/=BrOs;
      if(fabs(dip.fz)<sOOPrange){
         gyr = 0.0;
         for(j=0;j<BrOs;j++){
	 gyr += sq(dip.fx - pos[k+j].fx) + sq(dip.fy - pos[k+j].fy) 
		 + sq(dip.fz - pos[k+j].fz);
         }
	 radGyr[(int)((dip.fz+sOOPrange)/sOOPbinz)]+=(gyr/(float)(BrOs));
	 radGyrn[(int)((dip.fz+sOOPrange)/sOOPbinz)]++;
      }
   }
}
//fprintf(stderr,"wxtr.c: after getBrO.\n");
if(nw>0){
   for(i=0;i<nw;i+=3){
      getWatAng(i,tdpoint);
   }
}
//fprintf(stderr,"wxtr.c: after getWat.\n");
/* detect hydrogen bonds */
if(nBCD>0 || nw>0)
   getHB();
//fprintf(stderr,"wxtr.c: after getHB.\n");
/* calculate BCD orientational distributions */
if(fabs(BCDcom.fz) < BCDODrange){
   BCDOD[(int)((BCDcom.fz+BCDODrange)/BCDODbinz)][(int)((BCDz.fz+1.0)/BCDODbin)]++;
   BCDODnorm[(int)((BCDcom.fz+BCDODrange)/BCDODbinz)]++;
}
BCDODt[(int)((BCDz.fz+1.0)/BCDODbin)]++;
BCDODtnorm++;
/* calculate pore water diploes */
for(i=0;i<maxPoreWat;i++){
   if(poreW[i]!=-1){
      poreWdat[i][tdpoint].zpos = (float)(poreW[i]);
      havg.fx = (pos[i+1].fx + pos[i+2].fx)/2.0;
      havg.fy = (pos[i+1].fy + pos[i+2].fy)/2.0;
      havg.fz = (pos[i+1].fz + pos[i+2].fz)/2.0;
      dip.fx = havg.fx - pos[i].fx;
      dip.fy = havg.fy - pos[i].fy;
      dip.fz = havg.fz - pos[i].fz;
      r = sqrt(dip.fx*dip.fx + dip.fy*dip.fy + dip.fz*dip.fz);
      poreWdat[i][tdpoint].fx = (float)(dip.fx/r);
      poreWdat[i][tdpoint].fy = (float)(dip.fy/r);
      poreWdat[i][tdpoint].fz = (float)(dip.fz/r);
   }
}
// calculate BCD displacement
// and molecular vector orientational TCF 
if(nBCD == 1 && tdpoint>0){
   comhist[tdpoint].fx = BCDcom.fx+comhist[tdpoint-1].fx;	
   comhist[tdpoint].fy = BCDcom.fy+comhist[tdpoint-1].fy;
   if(nw>0 && nBrO >0){
	comhist[tdpoint].fz = BCDcom.fz;
   }
   else{
	comhist[tdpoint].fz = BCDcom.fz+comhist[tdpoint-1].fz;
   }
   zhist[tdpoint].fx = BCDz.fx;
   zhist[tdpoint].fy = BCDz.fy;
   zhist[tdpoint].fz = BCDz.fz;
   for(k=0;k<tdpoint+1;k++){
	dip.fx = comhist[tdpoint].fx - comhist[k].fx;
	dip.fy = comhist[tdpoint].fy - comhist[k].fy;
	dip.fz = comhist[tdpoint].fz - comhist[k].fz;
	sqDiv[0][tdpoint-k] += sq(dip.fx)+sq(dip.fy)+sq(dip.fz);
	sqDiv[1][tdpoint-k] += sq(dip.fz);
	sqDiv[2][tdpoint-k] += sq(dip.fx)+sq(dip.fy);
	sqDiv[3][tdpoint-k] += 1.0; // normalization constant
	BCDOTCF[0][tdpoint-k] += BCDz.fx*zhist[k].fx + BCDz.fy*zhist[k].fy +  BCDz.fz*zhist[k].fz;
	BCDOTCF[1][tdpoint-k] += 1.0;
   }
}
tdpoint++;
if(nBCD>0)
   offsetBCD();
/* Calculate BCD z-probability and the density profiles */
index = (int)((BCDcom.fz+zwall)/(binSize));
if(index >= npoints)
	ERROR((stderr,"wxtr: out of range\n"), exit);
BCDden[index].fz += 1.;
for(i=0;i<nw;i=i+3)
{
	index = (int) ((pos[i].fz + zwall)/(binSize));
	if(index >= npoints)
		ERROR((stderr,"wxtr: out of range\n"), exit);
	H2Oden[index].fz += 1.;
}
if(nw>0){
for(i=0;i<nBrO;i++)
{
//fprintf(stderr,"nw = %d, i = %d\n",nw,i);
	zavg = 0.0;
	for(j=0;j<BrOs;j++){
	   zavg += pos[nw+i*BrOs+j].fz*mass[nw+i*BrOs+j];
	}
	zavg /= BMass;   
	index = (int) ((zavg + zwall)/binSize);
	if (index >= npoints)
		ERROR((stderr,"wxtr: BrO com out of range, index = %d i = %d\n",index,i), exit);
	BrOden[0][index].fz += 1.;
	index = (int) ((pos[nw+BrOs*i].fz + zwall)/(binSize));
	if (index >= npoints)
		ERROR((stderr,"wxtr: BrO aC out of range, index = %d i = %d\n",index,i), exit);
	BrOden[1][index].fz += 1.;
	index = (int)((pos[nw+BrOs*i+1].fz + zwall)/(binSize));
	if (index >= npoints){
		ERROR((stderr,"wxtr: BrO Br out of range, index = %d i = %d,atom = %d\n",index,i,nw+i*BrOs+1), exit);
	}
	BrOden[2][index].fz += 1.;
}
}
//fprintf(stderr,"wxtr.c: completed dataPoint updates.\n");
if ((tdpoint-1) % pFreq == 0 && tdpoint > 1) {
	/* print run information */
	fprintf(fp,"nBCD\tnBrO\tnsolute\tnatoms\n");
	fprintf(fp,"%d\t%d\t%d\t%d\n",nBCD,nBrO,nsolute,natoms);
        fprintf(fp,"\n");
	/* print density profiles */
	n_factor.fz = 1./(4*xwall*ywall*binSize*tdpoint);
	fprintf(fp,"density profs / BCD CoM prob\n\n");
	fprintf(fp,"z\tH2O\tBrO(CoM)\tBrO(aC)\tBrO(Br)\tBCD\n");
	for(i=0;i<npoints-1;i++)
	{
	   	fprintf(fp,"%f\t%f\t%f\t%f\t%f\t%f\n",
                   -zwall+(i+0.5)*binSize,H2Oden[i].fz*n_factor.fz/H2Odensity/*29.91565*/,
	           BrOden[0][i].fz*n_factor.fz/BrOdensity/*286.861733*/,
	           BrOden[1][i].fz*n_factor.fz/BrOdensity/*286.861733*/,
	           BrOden[2][i].fz*n_factor.fz/BrOdensity/*286.861733*/,
		   BCDden[i].fz/tdpoint);
	}
	/* print radial distribution functions */
	gfactor[0]=(tdpoint-1)*dataRatex*4*PI*binRDF*binRDF*binRDF*H2Odensity;
	gfactor[1]=(tdpoint-1)*dataRatex*4*PI*binRDF*binRDF*binRDF*BrOdensity;
	fprintf(fp,"\nradial density functions - vs BCD CoM\n\n");
	fprintf(fp,"r\tH2O(O)\tiO\tH2O(CoM)\tiCoM\tBrO(alphC)\tiaC\tBrO(CoM)\tiCoM\tBrO(Br)\tiBr\tBrO(C8)\tiC8\n");
	for(i=0;i<9;i++) 
		idens[i] = 0.0;
	for(i=0;i<RDFbins;i++){
		fprintf(fp,"%f\t",i*binRDF);
		for(k=0;k<2;k++){
			idens[k]+=(grSOL[k][i]/(dataRatex*(tdpoint-1))); //((gfactor[0]/H2Odensity)*(1.0/3+i*(i+1)));
		//	idens[k]+=grSOL[k][i]/((gfactor[0]/H2Odensity)*(1.0/3+i*(i+1)));
			fprintf(fp,"%f\t%f\t",grSOL[k][i]/(gfactor[0]*(1.0/3.0+i*(i+1))),idens[k]);
		}
		for(k=3;k<7;k++){
			idens[k]+=(grSOL[k][i]/(dataRatex*(tdpoint-1)));//((gfactor[1]/BrOdensity)*(1.0/3+i*(i+1)));
		//	idens[k]+=grSOL[k][i]/((gfactor[1]/BrOdensity)*(1.0/3+i*(i+1)));
			fprintf(fp,"%f\t%f\t",grSOL[k][i]/(gfactor[1]*(1.0/3.0+i*(i+1))),idens[k]);
		}
		fprintf(fp,"\n");
	}
	/* print pore guest molecule probabilities */
	fprintf(fp,"\nPore guest molecule probabilities: water\n");
	fprintf(fp,"\nnWat\tP(nWat)\tcounts=\t%d\n",probPoreWat[maxPoreWat+1]);
	for(i=0;i<maxPoreWat+1;i++){
	   fprintf(fp,"%d\t%5.4f\n",i,(double)(probPoreWat[i])/(double)(probPoreWat[maxPoreWat+1]));
	}
	/* print pore guest molecule probabilities */
	fprintf(fp,"\nPore guest molecule probabilities: Br-oct CoM,a-C,Br,C8\n");
	fprintf(fp,"\nnBrO\tP(CoM)\tP(a-C)\tP(Br)\tP(C8)\tcounts=\t%d\n",probPoreBrO[maxPoreB+1]);
	for(i=0;i<maxPoreB+1;i++){
	   fprintf(fp,"%d\t",i);
	   for(j=0;j<4;j++){
	      fprintf(fp,"%5.4f\t",(double)(probPoreBrO[j][i])/(double)(probPoreBrO[j][maxPoreB+1]));
	   }
	   fprintf(fp,"\n");
	}
	/* print BCD OTCD, MSD */
	fprintf(fp,"\nBCD orientational TCF and mean square displacement\n\n");
	fprintf(fp,"t\tOTCF\t3D\tz\tx-y\n");
	for(i=0;i<numDPx;i++)
	{
	   	fprintf(fp,"%f\t%f\t%f\t%f\t%f\n"
                  ,i*dataRatex*h/1000.0,BCDOTCF[0][i]/BCDOTCF[1][i],sqDiv[0][i]/sqDiv[3][i],sqDiv[1][i]/sqDiv[3][i],sqDiv[2][i]/sqDiv[3][i]);
	}
	/* print BCD ODs */
	fprintf(fp,"\nBCD orientational distributions\n");
	fprintf(fp,"cos(th)\t");
	for(i=0;i<2*(BCDODrange/BCDODbinz);i++){
	   fprintf(fp,"%4.2f\t",(double)(i)*BCDODbinz-BCDODrange);
	}
	fprintf(fp,"all_z\tBrOguest\n");
	for(i=0;i<(int)(2.0/BCDODbin);i++){
	   fprintf(fp,"%4.3f\t",(double)(i)*BCDODbin-1.0);
	   for(j=0;j<2*(BCDODrange/BCDODbinz);j++){
	      fprintf(fp,"%4.3f\t",(double)(BCDOD[j][i])/(double)(BCDODnorm[j]));
	   }
	   fprintf(fp,"%4.3f\t%4.3f\n",(double)(BCDODt[i])/(double)(BCDODtnorm),(double)(BrOgOD[i])/(double)(BrOgODnorm));
	}
	/* print solvent ODs */
	for(k=0;k<5;k++){
	   if(k==0)
	      fprintf(fp,"\n water dipole orientational distributions\n");
	   else if(k==1)
	      fprintf(fp,"\n aC->Br vec orientational distributions\n");
	   else if(k==2)
	      fprintf(fp,"\n C8->C7 vec orientational distributions\n");
	   else if(k==3)
	      fprintf(fp,"\n C5->C4 vec orientational distributions\n");
	   else if(k==4)
	      fprintf(fp,"\n C8->C1 vec orientational distributions\n");
	   fprintf(fp,"cos(th)\t");
	   for(i=0;i<2*(sODrange/sODbinz);i++){
	      fprintf(fp,"%4.2f\t",(double)(i)*sODbinz-sODrange);
	   }
	   fprintf(fp,"all_z\n");
	   for(i=0;i<(int)(2.0/sODbin);i++){
	      fprintf(fp,"%4.3f\t",(int)(i)*sODbin-1.0);
	      for(j=0;j<2*(sODrange/sODbinz);j++){
	         fprintf(fp,"%4.3f\t",(double)(sOD[k][j][i])/(double)(sODnorm[k][j]));
	      }
	      fprintf(fp,"%4.3f\n",(double)(sODt[k][i])/(double)(sODtnorm[k]));
	   }
	}
	/* print orientational order parameter & radius of gyration*/
	fprintf(fp,"\nOrientational order parameters.\n");
	fprintf(fp,"z\twatDip\taC->Br\tC8->C7\tC5->C4\tC8->C1\tradGyr\n");
	for(i=0;i<2*(sOOPrange/sOOPbinz);i++){
	   fprintf(fp,"%4.2f\t",(double)(i)*sOOPbinz-sOOPrange);
	   for(j=0;j<5;j++){
	      fprintf(fp,"%5.3f\t",sOOP[j][i]/(double)(sOOPnorm[j][i]));
	   }
	   fprintf(fp,"%5.3f\n",radGyr[i]/(float)(radGyrn[i]));
	}

//fprintf(stderr,"before water pore res time stuff\n");
	/* calculate water pore residence times and dipole correlation functions */
	/* count water molecules that appear in pore */
	wip = 0;  // "waters in pore"
	for(i=0;i<2000;i++){
	   watList[i]=-1;
	}
	for(j=0;j<numDPx;j++){
	   for(i=0;i<maxPoreWat+1;i++){
	      newWat = 0;
	      if((int)(poreWdat[i][j].zpos) != -1){
		 newWat = 1; // assume it's new to the list
		 for(k=0;k<wip;k++){
		    if((int)(poreWdat[i][j].zpos) == watList[k]){
		       newWat=0;
		    }
		 }
	      }
	      if(newWat == 1 && wip < 1999){
	         watList[wip]=(int)(poreWdat[i][j].zpos);
	         wip++;
//	 fprintf(stderr,"timestep = %d, row = %d, water = %d.\n",j,i,(int)(poreWdat[i][j].zpos));
	      }
	   }
	}
	if(wip > 1998){
	     fprintf(stderr,"wxtr.c: WARNING: max number of waters that moved through pore exceeded. Ignoring n > 2000.\n");
	}
	/* construct array of in/out of pore */
	if((poreRes = (int **) calloc(wip,sizeof(int *))) == NULL)
	   ERROR((stderr,"wxtr: out of core\n"), exit);
	for(i=0;i<wip;i++){
		if((poreRes[i] = (int *) calloc(numDPx,sizeof(int))) == NULL)
		   ERROR((stderr,"wxtr: out of core\n"), exit);
	}
	for(i=0;i<wip;i++){
	   for(j=0;j<numDPx;j++){
	      poreRes[i][j] = 0;
	      for(k=0;k<maxPoreWat+1;k++){
		 if((int)(poreWdat[k][j].zpos) == watList[i])
		    poreRes[i][j] = 1;
	      }
	   }
	}
	/* sweep array to remove 'short exits' */
	maxGap = (int)(pwtau*1000.0/(h*(float)(dataRatex)));
	for(i=0;i<wip;i++){
	   inp = 0;
	   timeout = 0;
	   for(j=0;j<numDPx;j++){
	      if(poreRes[i][j]==1){
		inp = 1;
		if(timeout>0 && timeout<maxGap){
		   for(k=j-timeout;k<j;k++){
		      poreRes[i][k] = 1;
		   }
		}
		timeout = 0;
	      }
	      else if(poreRes[i][j]==0){
		 if(inp==1){
		   timeout = 1;
		 }
		 else if(timeout>0){
	           timeout++;
		 }
		 inp=0;
	      }
	   }
	}
	/* sweep array to count number of occupancy events */
	oe = 0;
	for(i=0;i<wip;i++){
	   inp = 0;
	   if(poreRes[i][0]==1){
	     oe++;
	   }
	   for(j=1;j<numDPx;j++){
	      if(poreRes[i][j]==1 && poreRes[i][j-1]==0){
		      oe++;
	      }
	   }
	}
/*
fprintf(stderr,"wxtr.c: waters appeared in pore: %d. occupancy events: %d\n",wip,oe);
for(i=0;i<wip;i++){
   fprintf(stderr,"%d  ",watList[i]);
}
fprintf(stderr,"\n");
*/
	/* create new array of individual occupancy events */
	/* water cavity residence time correlation function */
	minres = (int)(mnpwtau*1000.0/(h*(float)(dataRatex)));
	if((wcrTCF = (int **) calloc(oe,sizeof(int *))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	for(i=0;i<oe;i++){
		if((wcrTCF[i] = (int *) calloc(numDPx,sizeof(int))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	}
	/* set zeroes */
	for(i=0;i<oe;i++){
	   for(j=0;j<numDPx;j++){
	      wcrTCF[i][j] = 0;
	   }
	}
	/* read and write occupancies */
	lines = 0;
	for(i=0;i<wip;i++){
	   length = 0;
	   for(j=0;j<numDPx;j++){
	      if(poreRes[i][j]==1){
		 length++;
		 if(j==numDPx-1){ // end of row
		    for(k=0;k<length;k++){ //here insert minimun residence check
		       if(length>=minres)
			  wcrTCF[lines][k]=1;
		    }
		    lines++;
		 }
	      }
	      else if(poreRes[i][j]==0 && j!=0){
		 if(poreRes[i][j-1]== 1){ //end of residence
		    for(k=0;k<length;k++){
		       if(length>=minres)
			  wcrTCF[lines][k]=1;
		    }
		    lines++;
		 }
		 length=0;
	      }
	   }
	}
	/* calculate wcrTCF */
	for(i=0;i<lines;i++){
	   for(j=0;j<numDPx;j++){
	      for(k=0;k<numDPx-j;k++){
		 wTCF[k]+=wcrTCF[i][j]*wcrTCF[i][j+k];
		 wTCFn[k]+=wcrTCF[i][j];
	      }
	   }
	}
	/* solvent molecular vector OTCFs */
	for(j=0;j<3;j++){
	   if((watDipTCF[j] = (double *) calloc(numDPx+1,sizeof(double))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	   if((watDipTCFn[j] = (int *) calloc(numDPx+1,sizeof(int))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	   if((watHHTCF[j] = (double *) calloc(numDPx+1,sizeof(double))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	   if((watHHTCFn[j] = (int *) calloc(numDPx+1,sizeof(int))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	   if((watHHP2TCF[j] = (double *) calloc(numDPx+1,sizeof(double))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	   for(k=0;k<3;k++){
	      if((BrOTCF[k][j] = (double *) calloc(numDPx+1,sizeof(double))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	      if((BrOTCFn[k][j] = (int *) calloc(numDPx+1,sizeof(int))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	   }
	}
	for(j=0;j<numDPx+1;j++){
	   for(i=0;i<3;i++){
	      for(k=0;k<3;k++){
		 watDipTCF[i][j] = watHHTCF[i][j] = watHHP2TCF[i][j] = BrOTCF[k][i][j] = 0.0;
		 watDipTCFn[i][j] = watHHTCFn[i][j] = BrOTCFn[k][i][j] = 0;
	      }
	   }
	}
	/* calculate solvent OTCFs */
	for(i=0;i<nw/3;i++){
	    for(j=0;j<3;j++){  // bins in z-space
		for(k=0;k<numDPx;k++){
		    for(t=0;t<numDPx-k;t++){
	                if(watDipVec[i][k].zpos > watTCFbins[j*2] && watDipVec[i][k].zpos < watTCFbins[j*2+1] &&
	                   watDipVec[i][k+t].zpos > watTCFbins[j*2] && watDipVec[i][k+t].zpos < watTCFbins[j*2+1]){
			    watDipTCF[j][t] += watDipVec[i][k].fx*watDipVec[i][k+t].fx +
			       watDipVec[i][k].fy*watDipVec[i][k+t].fy + watDipVec[i][k].fz*watDipVec[i][k+t].fz;
			    watDipTCFn[j][t]++;
		        }
	                if(watHHVec[i][k].zpos > watTCFbins[j*2] && watHHVec[i][k].zpos < watTCFbins[j*2+1] &&
	                   watHHVec[i][k+t].zpos > watTCFbins[j*2] && watHHVec[i][k+t].zpos < watTCFbins[j*2+1]){
			    watHHTCF[j][t] += watHHVec[i][k].fx*watHHVec[i][k+t].fx +
			       watHHVec[i][k].fy*watHHVec[i][k+t].fy + watHHVec[i][k].fz*watHHVec[i][k+t].fz;
			    watHHTCFn[j][t]++;
			    watHHP2TCF[j][t] += (3.0*(sq(watHHVec[i][k].fx*watHHVec[i][k+t].fx +
			       watHHVec[i][k].fy*watHHVec[i][k+t].fy + watHHVec[i][k].fz*watHHVec[i][k+t].fz))-1.0)/2.0;
			}
		    }
	        }
	    }
	}
	/*
	j = 2;
	for(i=0;i<nw/3;i++){
	 for(k=0;k<numDPx;k++){
	    for(t=0;t<numDPx-k;t++){
		watDipTCF[j][t] += watDipVec[i][k].fx*watDipVec[i][k+t].fx +
		     watDipVec[i][k].fy*watDipVec[i][k+t].fy + watDipVec[i][k].fz*watDipVec[i][k+t].fz;
		watDipTCFn[j][t]++;
		watHHTCF[j][t] += watHHVec[i][k].fx*watHHVec[i][k+t].fx +
		     watHHVec[i][k].fy*watHHVec[i][k+t].fy + watHHVec[i][k].fz*watHHVec[i][k+t].fz;
		watHHTCFn[j][t]++;
		watHHP2TCF[j][t] += (3.0*(sq(watHHVec[i][k].fx*watHHVec[i][k+t].fx +
		     watHHVec[i][k].fy*watHHVec[i][k+t].fy + watHHVec[i][k].fz*watHHVec[i][k+t].fz))-1.0)/2.0;
	    }
	 }
	}
	*/
	for(l=0;l<3;l++){
	   for(i=0;i<nBrO;i++){
	     for(j=0;j<3;j++){  // bins in z-space
		 for(k=0;k<numDPx;k++){
		    for(t=0;t<numDPx-k;t++){
	                if(BrOVec[l][i][k].zpos > BrOTCFbins[j*2] && BrOVec[l][i][k].zpos < BrOTCFbins[j*2+1] &&
	                   BrOVec[l][i][k+t].zpos > BrOTCFbins[j*2] && BrOVec[l][i][k+t].zpos < BrOTCFbins[j*2+1]){
			     BrOTCF[l][j][t] += BrOVec[l][i][k].fx*BrOVec[l][i][k+t].fx +
			       BrOVec[l][i][k].fy*BrOVec[l][i][k+t].fy + BrOVec[l][i][k].fz*BrOVec[l][i][k+t].fz;
			     BrOTCFn[l][j][t]++;
		        }
		    }
	         }
	     }
	   }
	}
	/* print water correlation TCF */
	fprintf(fp,"\npore water res, solvent orientation TCFs\n");
	fprintf(fp,"t=0 norms: wat: i: %d b: %d BrO(C-Br): i: %d b: %d\n",watDipTCFn[0][0],watDipTCFn[1][0],BrOTCFn[0][0],BrOTCFn[1][0]);
	fprintf(fp,"t\tpwrTCF\twDip,i\twDip,b\twDip,all\twHH,i\t,wHH,b\twHH,all\twHHP2,i\twHHP2,b\twHHP2,all\tC-Br,i\tC-Br,b\tC8-C7,i\tC8-C7,b\tC5-C4,i\tC4-C5,b\n");
	for(i=0;i<numDPx;i++){
	   fprintf(fp,"%6.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\n",
		   (double)(i*dataRatex)*h/1000.0,(double)(wTCF[i])/(double)(wTCFn[i]),
		   watDipTCF[0][i]/(float)(watDipTCFn[0][i]),
		   watDipTCF[1][i]/(float)(watDipTCFn[1][i]),
		   watDipTCF[2][i]/(float)(watDipTCFn[2][i]),
		   watHHTCF[0][i]/(float)(watHHTCFn[0][i]),
		   watHHTCF[1][i]/(float)(watHHTCFn[1][i]),
		   watHHTCF[2][i]/(float)(watHHTCFn[2][i]),
		   watHHP2TCF[0][i]/(float)(watHHTCFn[0][i]),
		   watHHP2TCF[1][i]/(float)(watHHTCFn[1][i]),
		   watHHP2TCF[2][i]/(float)(watHHTCFn[2][i]),
		   BrOTCF[0][0][i]/(float)(BrOTCFn[0][0][i]),
		   BrOTCF[0][1][i]/(float)(BrOTCFn[0][1][i]),
		   BrOTCF[1][0][i]/(float)(BrOTCFn[1][0][i]),
		   BrOTCF[1][1][i]/(float)(BrOTCFn[1][1][i]),
		   BrOTCF[2][0][i]/(float)(BrOTCFn[2][0][i]),
		   BrOTCF[2][1][i]/(float)(BrOTCFn[2][1][i]));
	}
	for(i=0;i<3;i++){
	   free(watDipTCF[i]);
	   free(watDipTCFn[i]);
	   free(watHHTCF[i]);
	   free(watHHP2TCF[i]);
	   free(watHHTCFn[i]);
	   for(j=0;j<3;j++){
	      free(BrOTCF[i][j]);
	      free(BrOTCFn[i][j]);
	   }
	}
/* hydrogen bond lifetime correlation funcions */
	/* create array forHBs */
	for(i=0;i<nBCDHB-21;i++){
		if((HB_TCF[i] = (int *) calloc(numDPx,sizeof(int))) == NULL)
		   ERROR((stderr,"wxtr: out of core\n"), exit);
		if((HB_TCFn[i] = (int *) calloc(numDPx,sizeof(int))) == NULL)
		   ERROR((stderr,"wxtr: out of core\n"), exit);
	}
	for(i=0;i<nBCDHB-21;i++){
	   for(j=0;j<numDPx;j++){
	      HB_TCF[i][j] = HB_TCFn[i][j] = 0;
	   }
	}
/*** sweep vector of each site to find HB partners ***/
	for(m=0;m<nBCDHB-21;m++){
	   hbp = 0; // h-bond partners
	   for(i=0;i<2000;i++){
	      hbpList[i] = -1;
	   }
/******* loop around each bonding site ******/
	   for(j=0;j<numDPx;j++){
	      newHB = 0;
	      if(HB[m][j]!=-1 && HB[m][j]<nw){ //only consider water HBs
		 newHB = 1; //asume new to list
		 for(k=0;k<hbp;k++){ // sweep current list, see if already there
		     if(HB[m][j]==hbpList[k]){
			newHB = 0;
		     }
		 }
		 if(newHB==1 && hbp < 1998){
			 hbpList[hbp] = HB[m][j];
			 hbp++;
		 }
		 if(m>20){  // BCD acceptor, check for secondary.
	            if(HB[m+21][j]!=-1 && HB[m][j]<nw){ //only consider water HBs
		 	newHB = 1; //asume new to list
		 	for(k=0;k<hbp;k++){
		     	   if(HB[m][j]==hbpList[k]){
			      newHB = 0;
		     	   }
			}
		    }
		 
		    if(newHB==1 && hbp < 1998){
			 hbpList[hbp] = HB[m][j];
			 hbp++;
		    }
		 }   
	   }
	}
	if(hbp > 1998){
	     fprintf(stderr,"wxtr.c: WARNING: max number of HB partners exceeded. Ignoring n > 2000.\n");
	}
/** construct array of bonded / not bonded **/
	if((hbRes = (int **) calloc(hbp,sizeof(int *))) == NULL)
	   ERROR((stderr,"wxtr: out of core\n"), exit);
	for(i=0;i<hbp;i++){
		if((hbRes[i] = (int *) calloc(numDPx,sizeof(int))) == NULL)
		   ERROR((stderr,"wxtr: out of core\n"), exit);
	}
	for(i=0;i<hbp;i++){
	   for(j=0;j<numDPx;j++){
	      hbRes[i][j] = 0;
	      if(HB[m][j] == hbpList[i])
		  hbRes[i][j] = 1;
	      else if(m>20){
		  if(HB[m+21][j] == hbpList[i])
		     hbRes[i][j] = 1;
	      }
	   }
	}
	/* sweep array to remove 'short exits' */
	maxGap = (int)(HBtau*1000.0/(h*(float)(dataRatex)));
	for(i=0;i<hbp;i++){
	   inp = 0;
	   timeout = 0;
	   for(j=0;j<numDPx;j++){
	      if(hbRes[i][j]==1){
		inp = 1;
		if(timeout>0 && timeout<maxGap){
		   for(k=j-timeout;k<j;k++){
		      hbRes[i][k] = 1;
		   }
		}
		timeout = 0;
	      }
	      else if(hbRes[i][j]==0){
		 if(inp==1){
		   timeout = 1;
		 }
		 else if(timeout>0){
	           timeout++;
		 }
		 inp=0;
	      }
	   }
	}
	/* sweep array to count number of h-bonding events */
	oe = 0;
	for(i=0;i<hbp;i++){
	   inp = 0;
	   if(hbRes[i][0]==1){
	     oe++;
	   }
	   for(j=1;j<numDPx;j++){
	      if(hbRes[i][j]==1 && hbRes[i][j-1]==0){
		      oe++;
	      }
	   }
	}
	/* create new array of individual occupancy events */
	if((hbrTCF = (int **) calloc(oe,sizeof(int *))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	for(i=0;i<oe;i++){
		if((hbrTCF[i] = (int *) calloc(numDPx,sizeof(int))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	}
	/* set zeroes */
	for(i=0;i<oe;i++){
	   for(j=0;j<numDPx;j++){
	      hbrTCF[i][j] = 0;
	   }
	}
	/* read and write occupancies */
	minbond = (int)(mnHBtau*1000.0/(h*(float)(dataRatex)));
	lines = 0;
	for(i=0;i<hbp;i++){
	   length = 0;
	   for(j=0;j<numDPx;j++){
	      if(hbRes[i][j]==1){
		 length++;
		 if(j==numDPx-1){ // end of row
		    for(k=0;k<length;k++){
		       if(length>=minbond) // <-- lines to length
			  hbrTCF[lines][k]=1;
		    }
		    if(length<HBhistlen){
		       if(m<7)
			 hbtype=0;
		       else if(m<21)
			 hbtype=1;
		       else if(m<28)
			 hbtype=2;
		       else if(m<42)
			 hbtype=3;
		       HBhist[hbtype][length]++;
		       HBhistn[hbtype]++;
		    }
		    lines++;
		 }
	      }
	      else if(hbRes[i][j]==0 && j!=0){
		 if(hbRes[i][j-1]== 1 ){ //end of residence
		    for(k=0;k<length;k++){ //here insert minimum bond duration check
		       if(length>=minbond) 
			  hbrTCF[lines][k]=1;
		    }
		    if(length<HBhistlen){
		       if(m<7)
			 hbtype=0;
		       else if(m<21)
			 hbtype=1;
		       else if(m<28)
			 hbtype=2;
		       else if(m<42)
			 hbtype=3;
		       HBhist[hbtype][length]++;
		       HBhistn[hbtype]++;
		    }
		    lines++;
		 }
		 length=0;
	      }
	   }
	}
	/* create new array of individual occupancy events */
/*	for(i=0;i<nBCDHB;i++){
		if((HB_TCF[i] = (int *) calloc(numDPx,sizeof(int))) == NULL)
		   ERROR((stderr,"wxtr: out of core\n"), exit);
		if((HB_TCFn[i] = (int *) calloc(numDPx,sizeof(int))) == NULL)
		   ERROR((stderr,"wxtr: out of core\n"), exit);
	}
	for(i=0;i<nBCDHB;i++){
	   for(j=0;j<numDPx;j++){
	      HB_TCF[i][j] = HB_TCFn[i][j] = 0;
	   }
	}*/
	/* calculate hbTCF for this site */
	for(i=0;i<lines;i++){
	   for(j=0;j<numDPx;j++){
	      for(k=0;k<numDPx-j;k++){
		 HB_TCF[m][k]+=hbrTCF[i][j]*hbrTCF[i][j+k];
		 HB_TCFn[m][k]+=hbrTCF[i][j];
	      }
	   }
	}
	/* free memory from this internal loop */
	for(i=0;i<hbp;i++){
	   free(hbRes[i]);
	}
	free(hbRes);
	for(i=0;i<oe;i++){
	   free(hbrTCF[i]);
	}
	free(hbrTCF);

	}
	/* set up TCF summation matrices */
	if((prdTCF = (int *) calloc(numDPx,sizeof(int))) == NULL)
	   ERROR((stderr,"wxtr: out of core\n"), exit);
	if((prdTCFn = (int *) calloc(numDPx,sizeof(int))) == NULL)
	   ERROR((stderr,"wxtr: out of core\n"), exit);
	if((praTCF = (int *) calloc(numDPx,sizeof(int))) == NULL)
	   ERROR((stderr,"wxtr: out of core\n"), exit);
	if((praTCFn = (int *) calloc(numDPx,sizeof(int))) == NULL)
	   ERROR((stderr,"wxtr: out of core\n"), exit);
	if((sdTCF = (int *) calloc(numDPx,sizeof(int))) == NULL)
	   ERROR((stderr,"wxtr: out of core\n"), exit);
	if((sdTCFn = (int *) calloc(numDPx,sizeof(int))) == NULL)
	   ERROR((stderr,"wxtr: out of core\n"), exit);
	if((saTCF = (int *) calloc(numDPx,sizeof(int))) == NULL)
	   ERROR((stderr,"wxtr: out of core\n"), exit);
	if((saTCFn = (int *) calloc(numDPx,sizeof(int))) == NULL)
	   ERROR((stderr,"wxtr: out of core\n"), exit);
	for(i=0;i<numDPx;i++){
	   prdTCF[i] = praTCF[i] = sdTCF[i] = saTCF[i] = 0;	
	   prdTCFn[i] = praTCFn[i] = sdTCFn[i] = saTCFn[i] = 0;	
	}
	/* sum for type of h-bond */
	for(i=0;i<7;i++){ // primary hydroxyl is donor
	   for(j=0;j<numDPx;j++){
	      prdTCF[j] += HB_TCF[i][j];
	      prdTCFn[j]+= HB_TCFn[i][j];
	   }
	}
	for(i=7;i<21;i++){ // sec.hydroxyl is donor
	   for(j=0;j<numDPx;j++){
	      sdTCF[j] += HB_TCF[i][j];
	      sdTCFn[j]+= HB_TCFn[i][j];
	   }
	}
	for(i=21;i<28;i++){ // primary hydroxyl is acceptor
	   for(j=0;j<numDPx;j++){
	      praTCF[j] += HB_TCF[i][j];
	      praTCFn[j]+= HB_TCFn[i][j];
//	      praTCF[j+21] += HB_TCF[i][j];
//	      praTCFn[j+21]+= HB_TCFn[i][j];
	   }
	}
	for(i=28;i<42;i++){ // sec.hydroxyl is acceptor
	   for(j=0;j<numDPx;j++){
	      saTCF[j] += HB_TCF[i][j];
	      saTCFn[j]+= HB_TCFn[i][j];
//	      saTCF[j+21] += HB_TCF[i][j];
//	      saTCFn[j+21]+= HB_TCFn[i][j];
	   }
	}
/* now consider intra-BCD hydrogen bonds 
 * note: only need to consider donor vectors */
	/* wipe out old arrays */
	for(i=0;i<nBCDHB-21;i++){
	   for(j=0;j<numDPx;j++){
	      HB_TCF[i][j] = HB_TCFn[i][j] = 0;
	   }
	}
/*** sweep vector of each site to find HB partners ***/
	for(m=7;m<21;m++){
	   hbp = 0; // h-bond partners
	   for(i=0;i<2000;i++){
	      hbpList[i] = -1;
	   }
/******* loop around each bonding site ******/
	   for(j=0;j<numDPx;j++){
	      newHB = 0;
	      if(HB[m][j]!=-1 && HB[m][j]>nw){ //only consider BCD-BCD HBs
		 newHB = 1; //assume new to list
		 for(k=0;k<hbp;k++){ // sweep current list, see if already there
		     if(HB[m][j]==hbpList[k]){
			newHB = 0;
		     }
		 }
		 if(newHB==1 && hbp < 1998){
			 hbpList[hbp] = HB[m][j];
			 hbp++;
		 }
	      }
	   }
	   if(hbp > 1998){
	      fprintf(stderr,"wxtr.c: WARNING: max number of HB partners exceeded. Ignoring n > 2000.\n");
	   }
/** construct array of bonded / not bonded **/
	   if((hbRes = (int **) calloc(hbp,sizeof(int *))) == NULL)
	      ERROR((stderr,"wxtr: out of core\n"), exit);
	   for(i=0;i<hbp;i++){
		if((hbRes[i] = (int *) calloc(numDPx,sizeof(int))) == NULL)
		   ERROR((stderr,"wxtr: out of core\n"), exit);
	   }   
	   for(i=0;i<hbp;i++){
	      for(j=0;j<numDPx;j++){
	         hbRes[i][j] = 0;
	         if(HB[m][j] == hbpList[i])
		  hbRes[i][j] = 1;
	      }
	   }
	   /* sweep array to remove 'short exits' */
	   maxiGap = (int)(iHBtau*1000.0/(h*(float)(dataRatex)));
	   for(i=0;i<hbp;i++){
	      inp = 0;
	      timeout = 0;
	      for(j=0;j<numDPx;j++){
	         if(hbRes[i][j]==1){
		   inp = 1;
		   if(timeout>0 && timeout<maxiGap){
		      for(k=j-timeout;k<j;k++){
		         hbRes[i][k] = 1;
		      }
		   }
		   timeout = 0;
	         }
	         else if(hbRes[i][j]==0){
		    if(inp==1){
		      timeout = 1;
		    }
		    else if(timeout>0){
	              timeout++;
		    }
		    inp=0;
	         } 
	      }
	   }
	   /* sweep array to count number of h-bonding events */
	   oe = 0;
	   for(i=0;i<hbp;i++){
	      inp = 0;
	      if(hbRes[i][0]==1){
	         oe++;
	      }
	      for(j=1;j<numDPx;j++){
	         if(hbRes[i][j]==1 && hbRes[i][j-1]==0){
		      oe++;
	         }
	      }
	   }
	   /* create new array of individual occupancy events */
	   if((hbrTCF = (int **) calloc(oe,sizeof(int *))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	   for(i=0;i<oe;i++){
		if((hbrTCF[i] = (int *) calloc(numDPx,sizeof(int))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	   }
	   /* set zeroes */
	   for(i=0;i<oe;i++){
	      for(j=0;j<numDPx;j++){
	         hbrTCF[i][j] = 0;
	      }
	   }
	   /* read and write occupancies */
	   minibond = (int)(mniHBtau*1000.0/(h*(float)(dataRatex)));
	   lines = 0;
	   for(i=0;i<hbp;i++){
	      length = 0;
	      for(j=0;j<numDPx;j++){
	         if(hbRes[i][j]==1){
		    length++;
		    if(j==numDPx-1){ // end of row
		       for(k=0;k<length;k++){
		          if(length>=minibond)
			     hbrTCF[lines][k]=1;
		       }
		       if(length<HBhistlen){
			  hbtype=4;
		          HBhist[hbtype][length]++;
		          HBhistn[hbtype]++;
		       }
		       lines++;
		    }
	         }
	         else if(hbRes[i][j]==0 && j!=0){
		    if(hbRes[i][j-1]== 1){ //end of residence
		       for(k=0;k<length;k++){
		          if(length>=minibond)
			     hbrTCF[lines][k]=1;
		       }
		       if(length<HBhistlen){
			  hbtype=4;
		          HBhist[hbtype][length]++;
		          HBhistn[hbtype]++;
		       }
		       lines++;
		    }
		    length=0;
	         }
	      }
 	   }
	   /* create new array of individual occupancy events */
/*	for(i=0;i<nBCDHB;i++){
		if((HB_TCF[i] = (int *) calloc(numDPx,sizeof(int))) == NULL)
		   ERROR((stderr,"wxtr: out of core\n"), exit);
		if((HB_TCFn[i] = (int *) calloc(numDPx,sizeof(int))) == NULL)
		   ERROR((stderr,"wxtr: out of core\n"), exit);
	}
	for(i=0;i<nBCDHB;i++){
	   for(j=0;j<numDPx;j++){
	      HB_TCF[i][j] = HB_TCFn[i][j] = 0;
	   }
	}*/
	/* calculate hbTCF for this site */
	   for(i=0;i<lines;i++){
	      for(j=0;j<numDPx;j++){
	         for(k=0;k<numDPx-j;k++){
		    HB_TCF[m][k]+=hbrTCF[i][j]*hbrTCF[i][j+k];
		    HB_TCFn[m][k]+=hbrTCF[i][j];
	         }
	      }
	   }
	   /* free memory from this internal loop */
	   for(i=0;i<hbp;i++){
	      free(hbRes[i]);
	   }
	   free(hbRes);
	   for(i=0;i<oe;i++){
	      free(hbrTCF[i]);
	   }
	   free(hbrTCF);

	}
	/* set up TCF summation matrices */
	if((BBsdTCF = (int *) calloc(numDPx,sizeof(int))) == NULL)
	   ERROR((stderr,"wxtr: out of core\n"), exit);
	if((BBsdTCFn = (int *) calloc(numDPx,sizeof(int))) == NULL)
	   ERROR((stderr,"wxtr: out of core\n"), exit);
	for(i=0;i<numDPx;i++){
	   BBsdTCF[i] = BBsdTCFn[i] = 0;	
	}
	for(i=7;i<21;i++){ // sec.hydroxyl is donor
	   for(j=0;j<numDPx;j++){
	      BBsdTCF[j] += HB_TCF[i][j];
	      BBsdTCFn[j]+= HB_TCFn[i][j];
	   }
	}

	/* print HBlifetime correlation TCF */
	fprintf(fp,"\nHB lifetime TCFs. maxGap: %d, minBond: %d,maxiGap: %d,minibond: %d\n",maxGap,minbond,maxiGap,minibond);
	fprintf(fp,"\n<intra-unit BCDHB> = %f, <primary intraBCD HB> = %f\n",(float)(iuBCDHB)/(float)(tdpoint),(float)(priBCDHB)/(float)(tdpoint));
	fprintf(fp,"t\tpriDon\tpriAcc\tsecDon\tsecAcc\tintraBCD\n");
	for(i=0;i<numDPx;i++){
	   fprintf(fp,"%6.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\n",
		   (double)(i*dataRatex)*h/1000.0,
		   (float)(prdTCF[i])/(float)(prdTCFn[i]),
		   (float)(praTCF[i])/(float)(praTCFn[i]),
		   (float)(sdTCF[i])/(float)(sdTCFn[i]),
		   (float)(saTCF[i])/(float)(saTCFn[i]),
		   (float)(BBsdTCF[i])/(float)(BBsdTCFn[i]));
	}
	/* print HB duration histograms */
	fprintf(fp,"\nHB duration histograms\n");
	fprintf(fp,"t\tpriDon\tpriAcc\tsecDon\tsecAcc\tintraBCD\n");
	for(i=0;i<HBhistlen;i++){
	   fprintf(fp,"%6.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\n",
		   (double)(i*dataRatex)*h/1000.0,
		   (float)(HBhist[0][i])/(float)(HBhistn[0]),
		   (float)(HBhist[1][i])/(float)(HBhistn[1]),
		   (float)(HBhist[2][i])/(float)(HBhistn[2]),
		   (float)(HBhist[3][i])/(float)(HBhistn[3]),
		   (float)(HBhist[4][i])/(float)(HBhistn[4]));
	}
	/*****************************/
	/* pore water residence time */
	/*****************************/
	for(i=0;i<3;i++){
	   if((pwTCF[i] = (double *) calloc(numDPx+1,sizeof(double))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	}
	for(i=0;i<3;i++){
	   for(j=0;j<numDPx+1;j++){
	      pwTCF[i][j] = 0.0;
	   }
	}
if(nBCD==1 && nw > 1){
	/* sweep out short exits */
	maxGap = (int)((pwtau*1000.0)/(h*dataRatex));
//	fprintf(stderr,"pore water res time maxGap = %d\n",maxGap); 
	for(i=0;i<nw/3;i++){
	   timeout = 0;
	   t = 0;
	   for(j=0;j<numDPx+1;j++){
	      if(watDipVec[i][j].inp==1){
		if(t==1){ // flag to detect entry into pore
		  if(timeout>0 && timeout<maxGap){
		     for(k=1;k<=timeout;k++){
			watDipVec[i][j-k].inp = 1;
		     }
		  }
		}
		t = 1;
		timeout = 0;
	      }
	      else if(watDipVec[i][j].inp==0){
		 if(t==1){
		     timeout++;
		 }
	      }
	      else if(watDipVec[i][j].inp==-1){
		 fprintf(stderr,"wxtr.c: problem with pore water residences... detected -1. i=%d,j=%d\n,",i,j);
		 exit(1);
	      }
	      else{
		 fprintf(stderr,"wxtr.c: problem with pore water residences... detected unknown value: %d. i=%d,j=%d\n,",watDipVec[i][j].inp,i,j);
	         exit(1);
	      }
	   }
	}	
	/* calculate pore water residence TCF */
	for(i=0;i<nw/3;i++){
	   oe = 0;
	   for(j=0;j<numDPx+1;j++){
	      if(watDipVec[i][j].inp==1){
		 oe = 1;  // water shows up in pore
	      }
	      else if(watDipVec[i][j].inp==0 && oe == 1){ //water has left pore. delete history.
		 for(k=0;k<=j;k++){
		     watDipVec[i][k].inp = 0;
		 }
		 oe = 0;
	      }
	      for(t=0;t<j;t++){ //calculate residence TCF, going backwards in time
		 pwTCF[0][t] += (watDipVec[i][j].fx*watDipVec[i][j-t].fx +
				 watDipVec[i][j].fy*watDipVec[i][j-t].fy +
				 watDipVec[i][j].fz*watDipVec[i][j-t].fz)*(double)(watDipVec[i][j-t].inp*watDipVec[i][j].inp); 
		 pwTCF[1][t] += watDipVec[i][j].inp*watDipVec[i][j-t].inp;
		 pwTCF[2][t] += watDipVec[i][j].inp*watDipVec[i][j].inp;
	      }
	   }
	}
}
	for(i=0;i<3;i++){
	   for(j=0;j<3;j++){
	      if((pbTCF[i][j] = (double *) calloc(numDPx+1,sizeof(double))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	   }
	}
	for(i=0;i<3;i++){
	   for(k=0;k<3;k++){
	      for(j=0;j<numDPx+1;j++){
	         pbTCF[i][k][j] = 0.0;
	      }
	   }
	}
if(nBCD==1 && nBrO > 1){
	/* sweep out short exits */
   maxGap = (int)((pbtau*1000.0)/(h*dataRatex));
   //	fprintf(stderr,"pore water res time maxGap = %d\n",maxGap); 
   for(m=0;m<3;m++){
	for(i=0;i<nBrO;i++){
	   timeout = 0;
	   t = 0;
	   for(j=0;j<numDPx+1;j++){
	      if(BrOVec[m][i][j].inp==1){
		if(t==1){ // flag to detect entry into pore
		  if(timeout>0 && timeout<maxGap){
		     for(k=1;k<=timeout;k++){
			BrOVec[m][i][j-k].inp = 1;
		     }
		  }
		}
		t = 1;
		timeout = 0;
	      }
	      else if(BrOVec[m][i][j].inp==0){
		 if(t==1){
		     timeout++;
		 }
	      }
	      else if(BrOVec[m][i][j].inp==-1){
		 fprintf(stderr,"wxtr.c: problem with pore BrO residences... detected -1. i=%d,j=%d\n,",i,j);
		 exit(1);
	      }
	      else{
		 fprintf(stderr,"wxtr.c: problem with pore BrO residences... detected unknown value: %d. i=%d,j=%d\n,",BrOVec[m][i][j].inp,i,j);
	         exit(1);
	      }
	   }
	}	
	/* calculate pore BrO residence TCF */
	for(i=0;i<nBrO;i++){
	   oe = 0;
	   for(j=0;j<numDPx+1;j++){
	      if(BrOVec[m][i][j].inp==1){
		 oe = 1;  // water shows up in pore
	      }
	      else if(BrOVec[m][i][j].inp==0 && oe == 1){ //water has left pore. delete history.
		 for(k=0;k<=j;k++){
		     BrOVec[m][i][k].inp = 0;
		 }
		 oe = 0;
	      }
	      for(t=0;t<j;t++){ //calculate residence TCF, going backwards in time
		 pbTCF[m][0][t] += (BrOVec[m][i][j].fx*BrOVec[m][i][j-t].fx +
				 BrOVec[m][i][j].fy*BrOVec[m][i][j-t].fy +
				 BrOVec[m][i][j].fz*BrOVec[m][i][j-t].fz)*(double)(BrOVec[m][i][j-t].inp*BrOVec[m][i][j].inp); 
		 pbTCF[m][1][t] += BrOVec[m][i][j].inp*BrOVec[m][i][j-t].inp;
		 pbTCF[m][2][t] += BrOVec[m][i][j].inp*BrOVec[m][i][j].inp;
	      }
	   }
	}
   }
}
	/* print */
	fprintf(fp,"pore guest residence TCF, OTCF\n");
	fprintf(fp,"t(ps)\tOTCF\tresTCF\tnorm");
	fprintf(fp,"\taC_O\taCres\tnorm");
	fprintf(fp,"\tC8_O\tC8res\tnorm");
	fprintf(fp,"\tbCoM_O\tCoMres\tnorm\n");
	for(i=0;i<numDPx;i++){
	   fprintf(fp,"%f\t%f\t%f\t%d",i*h*dataRatex/1000.0,pwTCF[0][i]/pwTCF[2][i],
			   pwTCF[1][i]/pwTCF[2][i],(int)(pwTCF[2][i]));
	   for(m=0;m<3;m++){
	      fprintf(fp,"\t%f\t%f\t%d",pbTCF[m][0][i]/pbTCF[m][2][i],
			   pbTCF[m][1][i]/pbTCF[m][2][i],(int)(pbTCF[m][2][i]));
	   }
	   fprintf(fp,"\n");
	}
	
	fprintf(fp,"reaction center - BCDcom distance\n");
	fprintf(fp,"r(Å)\traC-BCDcom\n");
	for(i=0;i<(int)(2.0*maxBrComr/0.5);i++){
	   fprintf(fp,"%f\t%f\n",(i+0.5)*0.5-maxBrComr,(double)(aCBCDr[i])/(double)(BrOgODnorm));
	}
	
	fprintf(fp,"avg number of HBs\n");
	fprintf(fp,"priD\tpriA\tsecD\tsecA\twater/nw\n");
	fprintf(fp,"%f\t%f\t%f\t%f\t%f\n",
			(float)(priD)/(float)(tdpoint),
			(float)(priA)/(float)(tdpoint),
			(float)(secD)/(float)(tdpoint),
			(float)(secA)/(float)(tdpoint),
			(float)(watHBtot)/(float)(tdpoint)/(float)((natoms-nBCD*BCDs-nsolute-nBrO*BrOs)/3));
	g2dwfactor = (numDP_x_dataRate*H2Odensity); 
	g2dbfactor = (numDP_x_dataRate*BrOdensity); 
	fprintf(fp,"BCD-water map BCDz axis: g(r,theta)\n");
	for(i=0;i<2*(int)(contRad/contBin)+1;i++){
	   for(j=0;j<(int)(contRad/contBin)+1;j++){
   	      fprintf(fp,"%f\t",BWmap[i][j]/( g2dwfactor* ( sq((j+1)*contBin) - sq(j*contBin) ) * contBin * PI));
	   }
	   fprintf(fp,"\n");
	}	   
	fprintf(fp,"BCD-water map XY-plane: g(r,theta)\n");
	for(i=0;i<2*(int)(xyMapy/contBin)+1;i++){
	   for(j=0;j<2*(int)(xyMapx/contBin)+1;j++){
   	      fprintf(fp,"%f\t",BWmapxy[i][j]/(contBin*contBin*2.0*xyMapz*g2dwfactor) );
	   }
	   fprintf(fp,"\n");
	}	   
	fprintf(fp,"BCD-BrO:CoM map BCDz axis: g(r,theta)\n");
	for(i=0;i<2*(int)(contRad/contBin)+1;i++){
	   for(j=0;j<(int)(contRad/contBin)+1;j++){
   	      fprintf(fp,"%f\t",BBcommap[i][j]/( g2dbfactor* ( sq((j+1)*contBin) - sq(j*contBin) ) * contBin * PI));
	   }
	   fprintf(fp,"\n");
	}	   
	fprintf(fp,"BCD-BrO:CoM map XY-plane: g(r,theta)\n");
	for(i=0;i<2*(int)(xyMapy/contBin)+1;i++){
	   for(j=0;j<2*(int)(xyMapx/contBin)+1;j++){
   	      fprintf(fp,"%f\t",BBcommapxy[i][j]/(contBin*contBin*2.0*xyMapz*g2dbfactor) );
	   }
	   fprintf(fp,"\n");
	}	   
	fprintf(fp,"BCD-BrO:aC map BCDz axis: g(r,theta)\n");
	for(i=0;i<2*(int)(contRad/contBin)+1;i++){
	   for(j=0;j<(int)(contRad/contBin)+1;j++){
   	      fprintf(fp,"%f\t",BaCmap[i][j]/( g2dbfactor* ( sq((j+1)*contBin) - sq(j*contBin) ) * contBin * PI));
	   }
	   fprintf(fp,"\n");
	}	   
	fprintf(fp,"BCD-BrO:aC map XY-plane: g(r,theta)\n");
	for(i=0;i<2*(int)(xyMapy/contBin)+1;i++){
	   for(j=0;j<2*(int)(xyMapx/contBin)+1;j++){
   	      fprintf(fp,"%f\t",BaCmapxy[i][j]/(contBin*contBin*2.0*xyMapz*g2dbfactor) );
	   }
	   fprintf(fp,"\n");
	}	   
	fprintf(fp,"BCD-BrO:C8 map BCDz axis: g(r,theta)\n");
	for(i=0;i<2*(int)(contRad/contBin)+1;i++){
	   for(j=0;j<(int)(contRad/contBin)+1;j++){
   	      fprintf(fp,"%f\t",BC8map[i][j]/( g2dbfactor*( sq((j+1)*contBin) - sq(j*contBin) ) * contBin * PI));
	   }
	   fprintf(fp,"\n");
	}	   
	fprintf(fp,"BCD-BrO:C8 map XY-plane: g(r,theta)\n");
	for(i=0;i<2*(int)(xyMapy/contBin)+1;i++){
	   for(j=0;j<2*(int)(xyMapx/contBin)+1;j++){
   	      fprintf(fp,"%f\t",BC8mapxy[i][j]/(contBin*contBin*2.0*xyMapz*g2dbfactor) );
	   }
	   fprintf(fp,"\n");
	}	   
	fprintf(fp,"BCD-BrO:Br map BCDz axis: g(r,theta)\n");
	for(i=0;i<2*(int)(contRad/contBin)+1;i++){
	   for(j=0;j<(int)(contRad/contBin)+1;j++){
   	      fprintf(fp,"%f\t",BBrmap[i][j]/( g2dbfactor*( sq((j+1)*contBin) - sq(j*contBin) ) * contBin * PI));
	   }
	   fprintf(fp,"\n");
	}	   
	fprintf(fp,"BCD-BrO:Br map XY-plane: g(r,theta)\n");
	for(i=0;i<2*(int)(xyMapy/contBin)+1;i++){
	   for(j=0;j<2*(int)(xyMapx/contBin)+1;j++){
   	      fprintf(fp,"%f\t",BBrmapxy[i][j]/(contBin*contBin*2.0*xyMapz*g2dbfactor) );
	   }
	   fprintf(fp,"\n");
	}	   
	/* print BCD neighboring solvent ODs */
	fprintf(fp,"BCD 'neighbor solvent' ODs\n");
	for(k=1;k<5;k++){
	   if(k==0)
	      fprintf(fp,"\n water dipole orientational distributions\n");
	   else if(k==1)
	      fprintf(fp,"\n aC->Br vec orientational distributions\n");
	   else if(k==2)
	      fprintf(fp,"\n C8->C7 vec orientational distributions\n");
	   else if(k==3)
	      fprintf(fp,"\n C5->C4 vec orientational distributions\n");
	   else if(k==4)
	      fprintf(fp,"\n C8->C1 vec orientational distributions\n");
	   fprintf(fp,"cos(th)\t");
	   for(i=0;i<2*(nODrange/nODbinz);i++){
	      fprintf(fp,"%4.2f\t",(double)(i)*nODbinz-nODrange);
	   }
	   fprintf(fp,"\n");
	   for(i=0;i<(int)(2.0/nODbin);i++){
	      fprintf(fp,"%4.3f\t",(int)(i)*nODbin-1.0);
	      for(j=0;j<2*(nODrange/nODbinz);j++){
	         fprintf(fp,"%4.3f\t",(double)(nOD[k][j][i])/(double)(nODnorm[k][j]));
	      }
	      fprintf(fp,"\n");
	   }
	}

	/* print EVB radial distribution functions */
	gfactor[1]=(tdpoint-1)*dataRatex*4*PI*binRDF*binRDF*binRDF*H2Odensity;
	gfactor[0]=(tdpoint-1)*dataRatex*4*PI*binRDF*binRDF*binRDF*BrOdensity;
	fprintf(fp,"\nradial density functions - vs SN2 Cl\n\n");
	fprintf(fp,"r\tC-aC\ti\tCl1-aC\ti\tCl2-aC\ti\tC-Br\ti\tCl1-Br\ti\tCl2-Br\ti\tC-bOH\ti\tCl1-bOH\ti\tCl2-bOH\ti\tC-H2O\ti\tCl1-H2O\ti\tCl2-H2O\ti\n");
	for(i=0;i<12;i++) 
	   idens[i] = 0.0;
	for(i=0;i<RDFbins;i++){
	   fprintf(fp,"%f\t",i*binRDF);
	   for(j=0;j<4;j++){
	      for(k=0;k<3;k++){
	         idens[j*3+k]+=(grCl[k][j][i]/(dataRatex*(tdpoint-1))); //((gfactor[0]/H2Odensity)*(1.0/3+i*(i+1)));
		 fprintf(fp,"%f\t%f\t",grCl[k][j][i]/(gfactor[(int)(j/2)]*(1.0/3.0+i*(i+1))),idens[j*3+k]);
	      }
	   }
	   fprintf(fp,"\n");
	}
//	fprintf(stderr,"\n...printing complete\n");
	      
//fprintf(stderr,"before free BCDOTCF loop\n");
for(i=0;i<2;i++){
	free(BCDOTCF[i]);
}
for(i=0;i<4;i++){
	free(sqDiv[i]);
}
for(i=0;i<maxPoreWat+1;i++){
	free(poreWdat[i]);
}
for(i=0;i<wip;i++){
	free(poreRes[i]);
}
for(i=0;i<oe;i++){
	free(wcrTCF[i]);
}
free(wcrTCF);
free(poreRes);
free(H2Oden);
for(i=0;i<3;i++){
   free(BrOden[i]);
}
free(BCDden);
free(comhist);
free(zhist);
free(wTCF);
free(wTCFn);
for(i=0;i<nw/3;i++){
	free(watDipVec[i]);
	free(watHHVec[i]);
}
free(watDipVec);
free(watHHVec);
for(j=0;j<3;j++){
   free(pwTCF[j]);
   for(i=0;i<3;i++){
      free(BrOVec[j][i]);
      free(pbTCF[i][j]);
   }
}
for(i=0;i<3;i++){
   free(BrOVec[i]);
}
for(i=0;i<nBCDHB;i++){
   free(HB[i]);
}
}
fflush(fp);
}

void getWatAng(int k, int td){ // water dipole for OD
   tripd vec,hvec;
   double r,costh,zed;
   hvec.fx = (pos[k+1].fx + pos[k+2].fx)/2.0;
   hvec.fy = (pos[k+1].fy + pos[k+2].fy)/2.0;
   hvec.fz = (pos[k+1].fz + pos[k+2].fz)/2.0;
   vec.fx = hvec.fx - pos[k].fx;
   vec.fy = hvec.fy - pos[k].fy;
   vec.fz = hvec.fz - pos[k].fz;
   r = sqrt(vec.fx*vec.fx + vec.fy*vec.fy + vec.fz*vec.fz);
   costh = vec.fz/r;
   zed = (pos[k].fz + hvec.fz)/2.0;
   watDipVec[k/3][td].fx = (float)(vec.fx/r);
   watDipVec[k/3][td].fy = (float)(vec.fy/r);
   watDipVec[k/3][td].fz = (float)(vec.fz/r);
   watDipVec[k/3][td].zpos = (float)(zed);
   watHHVec[k/3][td].zpos = (float)(hvec.fz);
   hvec.fx = pos[k+1].fx - pos[k+2].fx;
   hvec.fy = pos[k+1].fy - pos[k+2].fy;
   hvec.fz = pos[k+1].fz - pos[k+2].fz;
   r = sqrt(hvec.fx*hvec.fx + hvec.fy*hvec.fy + hvec.fz*hvec.fz);
   watHHVec[k/3][td].fx = (float)(hvec.fx/r);
   watHHVec[k/3][td].fy = (float)(hvec.fy/r);
   watHHVec[k/3][td].fz = (float)(hvec.fz/r);
   
   if(fabs(zed)<sODrange){
      sOD[0][(int)((zed+sODrange)/sODbinz)][(int)((1.0+costh)/sODbin)]++;
      sODnorm[0][(int)((zed+sODrange)/sODbinz)]++;
   }
   if(fabs(zed)<sOOPrange){
      sOOPnorm[0][(int)((zed+sOOPrange)/sOOPbinz)]++;
      sOOP[0][(int)((zed+sOOPrange)/sOOPbinz)]+=(3.0*costh*costh - 1.0)/2.0;
   }
   sODt[0][(int)((1.0+costh)/sODbin)]++;
   sODtnorm[0]++;
}

void getBrOAng(int k, int l, int m, int j, int td){
//if(tc==0) fprintf(stderr,"get BrO angle: k = %d.\n",k);
   tripd vec;
   int nw,nb;
   double r,costh,zed,crad;
   nw = natoms - nBCD*BCDs - nBrO*BrOs - nsolute;
   nb = (k - nw)/BrOs; //index of BrO (0 - nBrO)
//if(tc==0) fprintf(stderr,"get BrO angle: nb = %d.\n",nb);
   vec.fx = pos[k+l].fx - pos[k+m].fx;
   vec.fy = pos[k+l].fy - pos[k+m].fy;
   vec.fz = pos[k+l].fz - pos[k+m].fz;
   r = sqrt(vec.fx*vec.fx + vec.fy*vec.fy + vec.fz*vec.fz);
   costh = vec.fz/r;
   zed = (pos[k+l].fz + pos[k+m].fz)/2.0;
   if(j>=1 && j <=3){  // 0-2 are BrO vector indices
      BrOVec[j-1][nb][td].fx = (float)(vec.fx/r);
      BrOVec[j-1][nb][td].fy = (float)(vec.fy/r);
      BrOVec[j-1][nb][td].fz = (float)(vec.fz/r);
      BrOVec[j-1][nb][td].zpos = (float)(zed);
   }

   if(fabs(zed)<sODrange){ // 1-4 are solvent OD indices for BrO (water is 0)
      sOD[j][(int)((zed+sODrange)/sODbinz)][(int)((1.0+costh)/sODbin)]++;
      sODnorm[j][(int)((zed+sODrange)/sODbinz)]++;
      if(nBCD==1){
         crad = sqrt(sq(vec.fx-BCDcom.fx) + sq(vec.fy-BCDcom.fy));
         if(crad < cradMax){
	    nOD[j][(int)((zed+nODrange)/nODbinz)][(int)((1+costh)/nODbin)]++;
	    nODnorm[j][(int)((zed+nODrange)/nODbinz)]++;
	 }
      }
   }
   if(fabs(zed)<sOOPrange){
      sOOPnorm[j][(int)((zed+sOOPrange)/sOOPbinz)]++;
      sOOP[j][(int)((zed+sOOPrange)/sOOPbinz)]+=(3.0*costh*costh - 1.0)/2.0;
   }
   sODt[j][(int)((1.0+costh)/sODbin)]++;
   sODtnorm[j]++;
}

/* offsetBCD:
 *	This routine offsets all atoms by an appropriate amount so as to bring
 *	BCD CoM to the center of the box (0,0,0) 
 */

void offsetBCD()
{
//	fprintf(stderr,"in offsetBCD\n");
	double sn2mass;
	int	i,j,nw;
	tripd	image, off, sn2com;
	nw = natoms - nBCD*BCDs - nBrO*BrOs - nsolute;
	if (nBCD != 1) return;
 	if(nsolute==3){
	   sn2mass = mass[natoms-3] + mass[natoms-2] + mass[natoms-1];
   	   sn2com.fx = pos[natoms-3].fx*mass[natoms-3] + pos[natoms-2].fx*mass[natoms-2] + pos[natoms-1].fx*mass[natoms-1]; 
   	   sn2com.fy = pos[natoms-3].fy*mass[natoms-3] + pos[natoms-2].fy*mass[natoms-2] + pos[natoms-1].fy*mass[natoms-1]; 
  	   sn2com.fz = pos[natoms-3].fz*mass[natoms-3] + pos[natoms-2].fz*mass[natoms-2] + pos[natoms-1].fz*mass[natoms-1];
   	   sn2com.fx/=sn2mass;
   	   sn2com.fy/=sn2mass;
   	   sn2com.fz/=sn2mass;
   	   sn2BCDrad = sqrt(sq(sn2com.fx-BCDcom.fx) + sq(sn2com.fy-BCDcom.fy) + sq(sn2com.fz-BCDcom.fz));
	}
	off.fx = BCDcom.fx;
	off.fy = BCDcom.fy;
	BCDcom.fx = 0.0;
	BCDcom.fy = 0.0;
        if (nw > 0 && nBrO > 1) // l/l interface 
	   off.fz = 0.0;
	else{
	   off.fz = BCDcom.fz;
	   BCDcom.fz = 0.0;
	}
	for (i = 0; i < natoms; i++) {
		pos[i].fx -= off.fx;
		pos[i].fy -= off.fy;
		pos[i].fz -= off.fz;

		if ((atom[i].flags & A_MAJOR) && atom[i].param1) {
			image.fx = -pos[i].fx;
			image.fy = -pos[i].fy;
			image.fz = -pos[i].fz;
			mvimage(&pos[i]);
			image.fx += pos[i].fx;
			image.fy += pos[i].fy;
			image.fz += pos[i].fz;
			for (j=i+1; j <= i + atom[i].param1; j++) {
				pos[j].fx += image.fx;
				pos[j].fy += image.fy;
				pos[j].fz += image.fz;
			}
		}
		else if (atom[i].flags & A_MAJOR)
			mvimage(&pos[i]);
	}
}
void offsetSN2()
{
//	fprintf(stderr,"in offsetBCD\n");
	if(nsolute!=3) return;
	int	i,j,nw;
	tripd	image, off;
	nw = natoms - nBCD*BCDs - nBrO*BrOs - nsolute;
	off.fx = pos[natoms-3].fx;
	off.fy = pos[natoms-3].fy;
        if (nw > 0 && nBrO > 1) // l/l interface 
	   off.fz = 0.0;
	else{
	   off.fz = pos[natoms-3].fz;
	}
	for (i = 0; i < natoms; i++) {
		pos[i].fx -= off.fx;
		pos[i].fy -= off.fy;
		pos[i].fz -= off.fz;

		if ((atom[i].flags & A_MAJOR) && atom[i].param1) {
			image.fx = -pos[i].fx;
			image.fy = -pos[i].fy;
			image.fz = -pos[i].fz;
			mvimage(&pos[i]);
			image.fx += pos[i].fx;
			image.fy += pos[i].fy;
			image.fz += pos[i].fz;
			for (j=i+1; j <= i + atom[i].param1; j++) {
				pos[j].fx += image.fx;
				pos[j].fy += image.fy;
				pos[j].fz += image.fz;
			}
		}
		else if (atom[i].flags & A_MAJOR)
			mvimage(&pos[i]);
	}
}
