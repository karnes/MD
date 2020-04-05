/*
 *	wdat:
 *	This routine writes the .dat file 
 */
#include	<md.h>
#include	<system.h>
wdat(fp,initQ,pFreq)
	FILE	*fp;
	int initQ,pFreq;
{
   int i,j,k;
   char gly[14] = {'C','O','H','H',
		   'C','O','H','H','H',
		   'C','O','H','H','H'};
   char tha[19] = {'C','H','H','H',
		   'C','H','H',
		   'C','H','H',
		   'C','H','H',
		   'C','H','H',
		   'C','H','H'};
   if(nTHA==1)
	centerTHA();
   fprintf(fp,"%d\n",natoms);
   fprintf(fp,"timeStep = %d\n",tc);
   for(i=0;i<nGLY;i++){
      for(j=0;j<GLYsites;j++){ 
         fprintf(fp,"%c %f %f %f\n",gly[j],pos[i*GLYsites+j].fx,pos[i*GLYsites+j].fy,pos[i*GLYsites+j].fz);
      }
   }
   if(nTHA==1){
      k = nGLY*GLYsites;
      fprintf(fp,"%c %f %f %f\n",'N',pos[k].fx,pos[k].fy,pos[k].fz);
      for(i=1;i<THAsites;i++){
	 k = nGLY*GLYsites + i;
         fprintf(fp,"%c %f %f %f\n",tha[(i-1)%19],pos[k].fx,pos[k].fy,pos[k].fz);
      }
   }
   if(nBr==1){
      k = nGLY*GLYsites+nTHA*THAsites;
      fprintf(fp,"Br %f %f %f\n",pos[k].fx,pos[k].fy,pos[k].fz);
   }
   return(0);
}
