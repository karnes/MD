/*
 *	wdat:
 *
 *	This routine writes the .dat file, 
 *	currently using to write all positions of all CFMs
 */

#include	<md.h>
#include	<system.h>
//tripd		shpos[500];

wdat(fp,initQ)
	FILE	*fp;
	int initQ;
{
	int j,i;
	int CFMprint = 0;

/* print all stacks n>=3 */
	if(tc<numDPx*dataRatex) return (0);
//fprintf(stderr,"in wdat.c\n");
/*	for(i=0;i<2000;i++){
	   if(printStx[i] != -1){
	      CFMprint++;
	   }
	}
	fprintf(fp,"%d\n\n",CFMprint*5);
	for(j=0;j<2000;j++){
	   if(printStx[j] != -1){
	      i = printStx[j];
	      fprintf(fp,"C %f %f %f\n",pos[i*5+0].fx,pos[i*5+0].fy,pos[i*5+0].fz);
	      fprintf(fp,"Cl %f %f %f\n",pos[i*5+1].fx,pos[i*5+1].fy,pos[i*5+1].fz);
	      fprintf(fp,"Cl %f %f %f\n",pos[i*5+2].fx,pos[i*5+2].fy,pos[i*5+2].fz);
	      fprintf(fp,"Cl %f %f %f\n",pos[i*5+3].fx,pos[i*5+3].fy,pos[i*5+3].fz);
	      fprintf(fp,"H %f %f %f\n",pos[i*5+4].fx,pos[i*5+4].fy,pos[i*5+4].fz);
	   }
	}

*/

/* print everything */
	/*
	fprintf(fp,"%d\n\n",nCFM*5);

	for(i=0;i<nCFM;i++){
	   fprintf(fp,"C %f %f %f\n",pos[i*5+0].fx,pos[i*5+0].fy,pos[i*5+0].fz);
	   fprintf(fp,"Cl %f %f %f\n",pos[i*5+1].fx,pos[i*5+1].fy,pos[i*5+1].fz);
	   fprintf(fp,"Cl %f %f %f\n",pos[i*5+2].fx,pos[i*5+2].fy,pos[i*5+2].fz);
	   fprintf(fp,"Cl %f %f %f\n",pos[i*5+3].fx,pos[i*5+3].fy,pos[i*5+3].fz);
	   fprintf(fp,"H %f %f %f\n",pos[i*5+4].fx,pos[i*5+4].fy,pos[i*5+4].fz);
	}
*/


}
