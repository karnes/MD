/*
 *	This routine will write a status line to the file pointed to by fp,
 */

#include	<md.h>
#include	<system.h>
#include	<water.h>

wstat(fp)
FILE	*fp;
{
// adding solvation shell populations
   fprintf(fp, "%-19s%-12s%-12s%-12s%-12s%-12s","TIME", "TEMP", "ETOT", "INTER_X","X_C","Br-b-dist");
   fprintf(fp, "\n");
   fprintf(fp,"%-18.1f%-12.2f%-12.2f%-12.2f%-12.2f%-12.3f",etime,temp,E*KCAL,INTER_X*KCAL,X_C,XXdist);
   fprintf(fp, "\n");
}



