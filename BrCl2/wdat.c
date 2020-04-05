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
   fprintf(fp,"Br %f %f %f\n",pos[0].fx,pos[0].fy,pos[0].fz);
   fprintf(fp,"Cl %f %f %f\n",pos[1].fx,pos[1].fy,pos[1].fz);
   fprintf(fp,"Cl %f %f %f\n",pos[2].fx,pos[2].fy,pos[2].fz);
   return(0);
}
