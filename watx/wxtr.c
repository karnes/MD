#include	<md.h>
#include	<system.h>
#include	<math.h>

wxtr(fp, initQ, pFreq)
	FILE	*fp;
	int	initQ;	/*initialization flag	*/
	int	pFreq;	/* frequency to print output	*/
{
int i,j,k,nw,/*ion,*/bin;
int	index;
double gfactor[2];
tripd	n_factor,vec,dip;

nw = natoms-nsolute-14*nNIT;
//ion = natoms - nsolute;

if (initQ == 1 ) {
	tdpoint = 0;
	binSize = 0.5;
	npoints = (int) (2*zwall/(binSize));
	if((H2Oden = (tripd *) calloc(npoints,sizeof(tripd))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	if((NITden = (tripd *) calloc(npoints,sizeof(tripd))) == NULL)
		ERROR((stderr,"wxtr: out of core\n"), exit);
	/* find new bin size so we have whole bins*/
	binSize = 2*zwall/npoints; /* exact bin size*/
	for(i=0;i<npoints;i++){
		H2Oden[i].fz = 0.;
		NITden[i].fz = 0.;
	}
	for(i=0;i<8;i++){
		for(bin=0;bin<300;bin++){
			grSOL[i][bin]=0.0;
		}
	}
	poshist[0].fx=watx.fx;
	poshist[0].fy=watx.fy;
	poshist[0].fz=watx.fz;
	for(k=0;k<tcData;k++){
		for(i=0;i<3;i++){
			sqDiv[i][k]=sqDivNorm[k]=0.0;
		}
		CorSh[k]=CorNorm[k]=0.0;
		avCorSh[k]=0.0;
		avCordip[k]=0.0;
		hhCor[0][k]=hhCor[1][k]=hhCor[2][k]=0.0;
		for(i=0;i<maxW;i++){
			hhV[i][k].fx=hhV[i][k].fy=hhV[i][k].fz=0.0;
			outSh[i]=0;
		}
	}
	for(i=0;i<rODbins;i++){
		for(k=0;k<101;k++){
			pwdip[k]=0.0;
			pwdip2[i][k]=0.0;
		}
		wdipNorm2[i]=0.0;
	}
	wdipNorm=0.0;
	for(i=0;i<zODbins;i++){
		for(k=0;k<101;k++){
			zODdip[i][k]=0.0;
		}
		zdipNorm[i]=0.0;
	}
	for(i=0;i<NITzODbins;i++){
		for(k=0;k<101;k++){
			NITzOD[i][k]=0.0;
		}
		NITzNorm[i]=0.0;
	}
}
/***	Calculate the density profiles	***/
for(i = 0; i < nw; i = i + 3)
{
	index = (int) ((pos[i].fz + zwall)/(binSize));
	if(index >= npoints)
		ERROR((stderr,"wxtr: out of range\n"), exit);
	H2Oden[index].fz += 1.;
}
for(i = 0; i < nNIT; i++)
{
	vec.fz = pos[nw+14*i].fz;   /* use image carbon */
	index = (int) ((vec.fz + zwall)/(binSize));
	if (index >= npoints)
		ERROR((stderr,"wxtr: out of range, index = %d i = %d\n",index,i), exit);
	NITden[index].fz += 1.;
}

analyze_cshel(tdpoint);
gethhCor(tdpoint);

for(k=0;k<tdpoint;k++){
	dip.fx = watx.fx - poshist[k].fx;
	dip.fy = watx.fy - poshist[k].fy;
	dip.fz = watx.fz - poshist[k].fz;
	pimage(&dip);
	sqDiv[0][tdpoint-k] += sq(dip.fx)+sq(dip.fy)+sq(dip.fz);
	sqDiv[1][tdpoint-k] += sq(dip.fz);
	sqDiv[2][tdpoint-k] += sq(dip.fx)+sq(dip.fy);
	sqDivNorm[tdpoint-k] += 1.0;
}
poshist[tdpoint].fx = watx.fx;
poshist[tdpoint].fy = watx.fy;
poshist[tdpoint].fz = watx.fz;
tdpoint++;
//fprintf(stderr,"tdpoint=%d pFreq=%d\n",tdpoint,pFreq);
if((tdpoint-1) % pFreq == 0 && tdpoint > 1) {
	/* print density profiles */
	n_factor.fz = 1./(4*xwall*ywall*binSize*tdpoint);
	fprintf(fp,"density profiles\n\n");
	fprintf(fp,"z\tNIT\tH2O\n");
	for(i=0;i<npoints-1;i++)
	{
	   	fprintf(fp,"%f\t%f\t%f\n"
                  ,-zwall+(i+0.5)*binSize,NITden[i].fz*n_factor.fz*216.2554,
	           H2Oden[i].fz*n_factor.fz*29.91565);
	}
	/* print radial distribution functions */
	gfactor[0]=(tdpoint-1)*dataRatex*4*PI*binRDF*binRDF*binRDF*H2Odensity;
	gfactor[1]=(tdpoint-1)*dataRatex*4*PI*binRDF*binRDF*binRDF*NITdensity;
	fprintf(fp,"\nradial density functions\n\n");
	fprintf(fp,"z\tO(w)\tH(w)\tC(a)\tC(o)\tC(m)\tC(p)\tN\tO(n)\n");
	for(i=20;i<300;i++){
		fprintf(fp,"%f\t",i*binRDF);
		for(k=0;k<2;k++){
			fprintf(fp,"%f\t",grSOL[k][i]/(gfactor[0]*(1.0/3.0+i*(i+1))));
		}
		for(k=2;k<8;k++){
			fprintf(fp,"%f\t",grSOL[k][i]/(gfactor[1]*(1.0/3.0+i*(i+1))));
		}
		fprintf(fp,"\n");
	}
//fprintf(stderr,"\nradial density functions complete\n\n");
	/* print time-dependent mean square dist. */
	fprintf(fp,"\ntime-dependent MSD, corr\n\n");
	fprintf(fp,"t(ps)\tMSD\tMSD(z)\tMSD(xy)\tCorSh\tavCorSh\tCordip\tavCordip\tCorResTCF\n");
        fprintf(fp,"0.000000\t0.000000\t0.000000\t0.000000\t%f\t%f\t%f\t%f\t%f\n",CorSh[0]/(double)(tcData),
		avCorSh[0],Cordip[0]/(double)(tcData),avCordip[0],CorSh[0]/CorNorm[0]);
	for(i=1;i<tdpoint;i++){
        	fprintf(fp,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",dataRatex*i/2000.0,
			sqDiv[0][i]/sqDivNorm[i],sqDiv[1][i]/sqDivNorm[i],sqDiv[2][i]/sqDivNorm[i],
			CorSh[i]/((double)(tcData-i)),
			avCorSh[i],Cordip[i]/((double)(tcData-i)),avCordip[i],CorSh[i]/CorNorm[i]);
	}
//fprintf(stderr,"\ntime-dependent MSD, corr complete\n\n");
	/* print hyrdation shell water dipole orientation wrt h2o->ion vector  */
	fprintf(fp,"\nshell water dipole orientation wrt water->ion vector\n\n");
	fprintf(fp,"cos(th)\tP(cos(th))\n");
	for(i=0;i<100;i++){
		fprintf(fp,"%f\t%f\n",i*0.02-0.99,pwdip[i]/(wdipNorm*1.0));
	}
//fprintf(stderr,"\nshell water dipole orientation wrt water->ion vector complete\n\n");
	/* print radius-binned water dipole orientation wrt h2o->ion vector  */
	fprintf(fp,"\nr-binned water dipole orientation wrt water->ion vector\n\n");
	fprintf(fp,"cos(th)\t");
	for(i=2;i<rODbins;i++){	
		fprintf(fp,"%3.2f\t",i*rODbinw);
	}
	fprintf(fp,"\n");
	for(i=0;i<100;i++){
		fprintf(fp,"%f\t",i*0.02-0.99);
		for(k=2;k<rODbins;k++){
			fprintf(fp,"%f\t",pwdip2[k][i]/(wdipNorm2[k]*1.0));
		}
		fprintf(fp,"\n");
	}
//fprintf(stderr,"\nr-binned water dipole orientation wrt water->ion vector complete\n\n");
	/* print water HH vector orientational TCF */
	fprintf(fp,"\nwater orientational TCF\n\n");
	fprintf(fp,"t(ps)\tHH_OTCF\thhCor[1]\thhCor[2]\n");
	for(i=0;i<tdpoint;i++){
		fprintf(fp,"%f\t%f\t%f\t%f\n",dataRatex*i/2000.0,
			hhCor[0][i]/hhCor[2][i]/*((double)(tcData-i))*/,hhCor[1][i]/hhCor[2][i]/*((double)(tcData-i))*/,hhCor[2][i]);
	}
//fprintf(stderr,"\nwater orientational TCF complete\n\n");
	/* print water dipole vector ODs */
	fprintf(fp,"\nwater dipole ODs by z\n\n");
	fprintf(fp,"cos(th)\t");
	for(i=0;i<zODbins;i++){	
		fprintf(fp,"%3.2f\t",-1.0*i*zODbinw);
	}
	fprintf(fp,"\n");
	for(i=0;i<100;i++){
		fprintf(fp,"%f\t",i*0.02-0.99);
		for(k=0;k<zODbins;k++){
			fprintf(fp,"%f\t",zODdip[k][i]/zdipNorm[k]);
		}
		fprintf(fp,"\n");
	}
	/* print nitrobenzene dipole vector ODs */
	fprintf(fp,"\nnitrobenzene dipole ODs by z\n\n");
	fprintf(fp,"cos(th)\t");
	for(i=0;i<NITzODbins;i++){	
		fprintf(fp,"%3.2f\t",i*NITzODbinw);
	}
	fprintf(fp,"\n");
	for(i=0;i<100;i++){
		fprintf(fp,"%f\t",i*0.02-0.99);
		for(k=0;k<NITzODbins;k++){
			fprintf(fp,"%f\t",NITzOD[k][i]/NITzNorm[k]);
		}
		fprintf(fp,"\n");
	}
}	

fflush(fp);

}

void analyze_cshel(t)
int t;
{
/*
 * Cshel[i] contains -2 for any water molecule that at the current time t
 * is not in the first hydration shell and the actual cos of the angle between
 * the water dipole and the water ion vector if this water molecule is in
 * the first hydration shell.
*/
int i,j, k, nw;
nw = (natoms - nsolute - 14*nNIT)/3;
for (i=0;i<nw;i++){
	if (Cshel[i] >= -1){
		Inshel[i][t] = 1;
		dipmat[i][t] = Cshel[i];
	}
	else {/*reset to zero all histories of a molecule that exited*/
		for (k = 0; k <=t;k++){
			Inshel[i][k] = 0;
			dipmat[i][k] = 0;
		}
	}
	avCorSh[t] += Inshel[i][t];
	avCordip[t] += dipmat[i][t];
	for (j=0;j<=t;j++){
		CorSh[t-j] += Inshel[i][t]*Inshel[i][j];
		Cordip[t-j] += dipmat[i][t]*dipmat[i][j];
		if(Inshel[i][t] == 1)
			CorNorm[j]+=1.0;
	}
}
}
/*calculate the water H-H correlation function for the water in the hydration shell*/
#define tau 0 //time we allow water to be outside shell (fs)
void gethhCor(t)
int t;
{
int i, k, nw;
double uv, c1, c2, sqrt();
tripd r;
nw = (natoms - nsolute - 14*nNIT)/3;
/* initialze all molecules that are in
 * the shell at t=0 to have a tolarence of tau fs  */
if(t == 0){
        for (i=0;i<nw;i++){
                if (Cshel[i] >= -1) outSh[i] = tau/dataRatex/2;
        }
}
for(i=0;i<nw;i++){
 if(Cshel[i] >= -1 || outSh[i] > 0){/*this water is in shell or out briefly*/
        if (Cshel[i] >= -1) outSh[i] = 1 + tau/dataRatex/2;/*reset tolerance if entered*/
	r.fx = pos[3*i+2].fx-pos[3*i+1].fx;
	r.fy = pos[3*i+2].fy-pos[3*i+1].fy;
	r.fz = pos[3*i+2].fz-pos[3*i+1].fz;
	uv = sqrt(r.fx*r.fx+r.fy*r.fy+r.fz*r.fz);
	hhV[i][t].fx = r.fx/uv;
	hhV[i][t].fy = r.fy/uv;
	hhV[i][t].fz = r.fz/uv;
	for (k=0; k<=t; k++){
		c1 = hhV[i][k].fx*hhV[i][t].fx 
		   + hhV[i][k].fy*hhV[i][t].fy 
		   + hhV[i][k].fz*hhV[i][t].fz; 
		c2 = (3*c1*c1-1)/2.0;
		hhCor[0][t-k] += c1;
		hhCor[1][t-k] += c2;
		if(c1!=0.0)
			hhCor[2][t-k] += 1.0;
	}
   }
   if(Cshel[i] < -1) outSh[i]--;
   if(outSh[i] <= 0){ /*this water remains outside more than 2 ps*/
        for (k=0; k<=t; k++){ /* delete history */
                hhV[i][k].fx = hhV[i][k].fy = hhV[i][k].fz = 0.0;
        }
     }
}
}
