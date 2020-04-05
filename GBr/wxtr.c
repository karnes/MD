#include	<md.h>
#include	<system.h>
#include	<math.h>

wxtr(fp, initQ, pFreq)
	FILE	*fp;
	int	initQ;	/*initialization flag	*/
	int	pFreq;	/* frequency to print output	*/
{

int i,j,k,l,m,n,bin;
double delfx,delfy,delfz;
int index;
double gfactor[2],r2;
tripd n_factor,sdl,image;

if(initQ == 1){
   tdpoint = 0;
   binSize = 0.2;
   npoints = (int) (2*zwall/(binSize));
   for(j=0;j<3;j++){
      if((GLYden[j] = (tripd *) calloc(npoints,sizeof(tripd))) == NULL)
         ERROR((stderr,"wxtr: out of core\n"), exit);
   }
   /* find new bin size so we have whole bins*/
   binSize = 2*zwall/npoints; /* exact bin size*/
   for(j=0;j<3;j++){
      for(i=0;i<npoints;i++){
         GLYden[j][i].fz = 0.;
      }
   }
   for(i=0;i<grGLYtypes;i++){
      // 0 = C-C
      // 1 = O - hydroxyl H
      for(bin=0;bin<300;bin++){
      grGLY[i][bin]=0.0;
      }
   }
}
/***	Calculate the density profiles	***/
for(i=0;i<nGLY*GLYsites;i=i+GLYsites)
{
   index = (int)((pos[i].fz + zwall)/(binSize));
   if(index >= npoints)
      ERROR((stderr,"wxtr: out of range\n"), exit);
   GLYden[0][index].fz += 1.;
   //carbon
   GLYden[1][index].fx += 1.0/3.0;
   index = (int)((pos[i+4].fz + zwall)/(binSize));
   GLYden[1][index].fx += 1.0/3.0;
   index = (int)((pos[i+9].fz + zwall)/(binSize));
   GLYden[1][index].fx += 1.0/3.0;
   //oxygen
   index = (int)((pos[i+1].fz + zwall)/(binSize));
   GLYden[2][index].fx += 1.0/3.0;
   index = (int)((pos[i+5].fz + zwall)/(binSize));
   GLYden[2][index].fx += 1.0/3.0;
   index = (int)((pos[i+10].fz + zwall)/(binSize));
   GLYden[2][index].fx += 1.0/3.0;
}

tdpoint++;

if((tdpoint-1) % pFreq == 0 && tdpoint > 1) {
	/* print density profiles */
	n_factor.fz = 1./(4.*xwall*ywall*binSize*tdpoint);
	fprintf(fp,"density profiles\n\n");
	fprintf(fp,"z\tGLY\n");
	for(i=625;i<npoints-1-625;i++)
	{
	   	fprintf(fp,"%f\t",-zwall+(i+0.5)*binSize);
		for(j=0;j<3;j++){
		   fprintf(fp,"%f\t",GLYden[j][i].fz*n_factor.fz*(1./GLYdensity));
		}
	   	fprintf(fp,"\n");
	}
	/* print radial distribution functions */
	gfactor[0] = /*tdpoint*/(int)(numDP_x_dataRate/grdt)*nGLY*4*PI*binRDF*binRDF*binRDF*GLYdensity;
//	gfactor[0] = numDP_x_dataRate*4*PI*binRDF*binRDF*binRDF*GLYdensity;
	fprintf(fp,"\nradial density functions\n\n");
	fprintf(fp,"z\tC-C\tO-H(OH)\tC-O\tC-H\tO-H(any)\tH-H\tO-O\t");
	for(i=0;i<3;i++){
	   fprintf(fp,"C-X%d\tO-X%d\tH(oh)-X%d\t",i,i,i);
	}
	fprintf(fp,"\n");
	for(i=20;i<300;i++){
		fprintf(fp,"%f\t",i*binRDF);
		for(k=0;k<grGLYtypes;k++){
			fprintf(fp,"%f\t",grGLY[k][i]/(gfactor[0]*(1.0/3.0+i*(i+1))));
		}
		fprintf(fp,"\n");
	}
}
return;

}
