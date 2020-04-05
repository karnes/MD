#include	<md.h>
#include	<system.h>
#include	<math.h>

wxtr(fp, initQ, pFreq)
    FILE	*fp;
    int	initQ;	/*initialization flag	*/
    int	pFreq;	/* frequency to print output	*/
{
    int i,j,index;
    int nw;
    double gfactor[2];
    double integr[2] = {0.0};
    tripd n_factor,vec;
    nw = natoms-nEr-nDDC*DDCsites;

    if(initQ == 1){
	tdpoint = 0;
	binSize = 0.5;
	npoints = (int) (2*zwall/(binSize));
	if((H2Oden = (tripd *) calloc(npoints,sizeof(tripd))) == NULL)
	    ERROR((stderr,"wxtr: out of core\n"), exit);
	if((DDCden = (tripd *) calloc(npoints,sizeof(tripd))) == NULL)
	    ERROR((stderr,"wxtr: out of core\n"), exit);
	/* find new bin size so we have whole bins*/
	binSize = 2*zwall/npoints; /* exact bin size*/
	for(i=0;i<npoints;i++){
	    H2Oden[i].fz = 0.;
	    DDCden[i].fz = 0.;
	}
	for(i=0;i<2;i++){
	    for(j=0;j<RDFbins;j++){
		grSOL[i][j]=0.0;
	    }
	}
    }

    /***    Calculate the density profiles  ***/
    for(i=0;i<nw;i=i+3)
    {
	index = (int) ((pos[i].fz + zwall)/(binSize));
	if(index >= npoints)
	    ERROR((stderr,"wxtr: out of range\n"), exit);
	H2Oden[index].fz += 1.;
    }
    for(i=0;i<nDDC;i++)
    {
	vec.fz = pos[nw+DDCsites*i].fz;   /* use image carbon */
	index = (int) ((vec.fz + zwall)/(binSize));
	if (index >= npoints)
	    ERROR((stderr,"wxtr: out of range, index = %d i = %d\n",index,i), exit);
	DDCden[index].fz += 1.;
    }   
    tdpoint++;
    if((tdpoint-1) % pFreq == 0 && tdpoint > 1) {
	/* print density profiles */
	n_factor.fz = 1./(4*xwall*ywall*binSize*tdpoint);
	fprintf(fp,"density profiles\n\n");
	fprintf(fp,"z\tDDC\tH2O\n");
	for(i=0;i<npoints-1;i++)
	{
	    fprintf(fp,"%f\t%f\t%f\n"
		    ,-zwall+(i+0.5)*binSize,DDCden[i].fz*n_factor.fz/DDCdensity,
		    H2Oden[i].fz*n_factor.fz/H2Odensity);
	}
	/* print radial distribution functions */
	gfactor[0]=(tdpoint-1)*dataRatex*4*PI*binRDF*binRDF*binRDF*H2Odensity;
	gfactor[1]=(tdpoint-1)*dataRatex*4*PI*binRDF*binRDF*binRDF*DDCdensity;
	fprintf(fp,"\nradial density functions\n\n");
	fprintf(fp,"z\tO(w)\tH(w)\tiO(w)\tiH(w)\n");
	for(i=20;i<RDFbins;i++){
	    fprintf(fp,"%f\t",i*binRDF);
	    for(j=0;j<2;j++){ // Er - water
		fprintf(fp,"%f\t", grSOL[j][i]/(gfactor[0]*(1.0/3.0+i*(i+1))));
	    }
	    for(j=0;j<2;j++){
		integr[j] += grSOL[j][i]/(dataRatex*(tdpoint-1));
		fprintf(fp,"%f\t",(j+1.0)*integr[j]);
	    }	
	    fprintf(fp,"\n");
	}
    }
    fflush(fp);
}
