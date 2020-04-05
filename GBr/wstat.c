/*
 *	This routine will write a status line to the file pointed to by fp,
 */

#include	<md.h>
#include	<system.h>
#include	<water.h>

wstat(fp)
FILE	*fp;
{
   int i;
   centerBr();
   //debug GLYForce.c
/*
   fprintf(fp, "%-10s%-17s%-17s%-17s%-17s%-17s%-17s%-17s%-17s%-17s","TIME", "TEMP", "ETOT", "GLYV", "GLYNB","gBond","gBend","gTors","g14","g15");
   fprintf(fp, "\n");
   fprintf(fp,"%-9.1f%-17.2f%-17.2f%-17.2f%-17.2f%-17.2f%-17.2f%-17.2f%-17.2f%-17.2f",etime,temp,E*KCAL,INTRA_GLY*KCAL,INTER_GLY*KCAL,gbondE*KCAL,gbendE*KCAL,gtorsE*KCAL,g14E*KCAL,g15E*KCAL);
   fprintf(fp, "\n");
*/
// adding solvation shell populations
   fprintf(fp, "%-10s%-12s%-12s%-12s%-12s%-12s%-12s%-12s%-12s%-10s%-10s%-10s","TIME", "TEMP", "ETOT", "GLYV", "GLYNB","XGLYV","V_vindow","V_bias","cosXz","Xz","sysCoM","BrZFrc");
   fprintf(fp, "\n");
   fprintf(fp,"%-9.1f%-12.2f%-12.2f%-12.2f%-12.2f%-12.2f%-12.2f%-12.2f%-12.2f%-10.3f%-10.3f%-10.3f",etime,temp,E*KCAL,INTRA_GLY*KCAL,INTER_GLY*KCAL,XGLYV*KCAL,V_w*KCAL,V_bias*KCAL,cosXz,Xz,syscomz,BrZForce*KCAL);
   fprintf(fp, "\n");
}

centerBr() //shift system so that THA N is at x=0,y=0
{
   int i,k;
   tripd cent;
   if(nBr!=1)
      return 0;
   k=nGLY*GLYsites;
   cent.fx = pos[k].fx;
   cent.fy = pos[k].fy;
   if(xwall==ywall && xwall==zwall)
      cent.fz = pos[k].fz;
   for(i=0;i<natoms;i++){
      pos[i].fx-=cent.fx;
      pos[i].fy-=cent.fy;
      if(xwall==ywall && xwall==zwall)
         pos[i].fz-=cent.fz;
   }
   return 1;
}

