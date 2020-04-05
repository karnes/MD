/*
 *	This routine will write a status line to the file pointed to by fp,
 */

#include	<md.h>
#include	<system.h>
#include	<water.h>

wstat(fp)
FILE	*fp;
{
int nw;
nw = natoms - 14*nNIT - 1;
offsetWatx(natoms-1);

if(nw>12){
   fprintf(fp, "%-10s%-7s%-9s%-9s%-11s","TIME", "TEMP", "ETOT", "WATERV", "H2ONB");
   fprintf(fp, "%-10s%-7s%-9s","NITV", "NITNB", "H2ONIT");
   fprintf(fp, "%-8s%-8s%-10s%-10s%-10s%-8s","VINT","VINTH2O","VINTH2O_s","VCH2O[0]","VCH2O[1]","Nshel");
   fprintf(fp, "%-9s%-7s%-10s%-8s%-8s%-8s%-8s","VINT_NIT", "VC_NIT", "V_Window","Ion_z","w_cor","V_bias","zCoM");
   fprintf(fp, "\n");
   fprintf(fp,"%-9.1f%-7.2f%-9.2f%-9.2f%-11.2f",etime,temp,E*KCAL,WATERV*KCAL,VWATTS*KCAL);
   fprintf(fp, "%-8.2f%-9.2f%-10.2f",NITV*KCAL,NITNB*KCAL,H2ONITV*KCAL);
   fprintf(fp, "%-9.2f%-9.2f%-10.2f%-10.2f%-10.2f%-7d",
	VINT*KCAL,VINTH2O*KCAL,VINTH2O_s*KCAL,VCH2O[0]*KCAL,VCH2O[1]*KCAL,Nshel);
   fprintf(fp, "%-8.2f%-9.2f%-8.2f%-8.3f%-8.3f%-8.3f%-8.3f",VINT_NIT*KCAL,VC_NIT*KCAL, V_w*KCAL,pos[natoms-1].fz,wcor,V_bias*KCAL,osysCmz);
   fprintf(fp, "\n");
//fprintf(fp, "%-8s%-8s%-10s%-10s%-10s%-8s","VINT","VINTH2O","VINTH2O_s","VCH2O[0]","VCH2O[1]","Nshel");
//fprintf(fp, "%-10s%-7s","VINT_NIT", "VC_NIT");
//fprintf(fp, "\n");
//fprintf(fp, "%-8.2f%-8.2f%-10.2f%-10.2f%-10.2f%-9d",
//	VINT*KCAL,VINTH2O*KCAL,VINTH2O_s*KCAL,VCH2O[0]*KCAL,VCH2O[1]*KCAL,Nshel);
//	fprintf(fp, "%-8.2f%-8.2f",VINT_NIT*KCAL,VC_NIT*KCAL);
//fprintf(fp, "\n");
}
else{
   fprintf(fp, "%-10s%-7s%-9s%-9s%-11s","TIME", "TEMP", "ETOT", "WATERV", "H2ONB");
   fprintf(fp, "%-10s%-7s%-9s","NITV", "NITNB", "H2ONIT");
   fprintf(fp, "%-8s%-8s%-10s%-10s%-10s%-8s","VINT","VINTH2O","VINTH2O_s","VCH2O[0]","VCH2O[1]","Nshel");
   fprintf(fp, "%-9s%-7s%-10s%-8s%-8s%-8s%-8s","VINT_NIT", "VC_NIT", "  WW_C","Ion_z","NBNB_C","V_bias"," WNITC");
   fprintf(fp, "\n");
   fprintf(fp,"%-9.1f%-7.2f%-9.2f%-9.2f%-11.2f",etime,temp,E*KCAL,WATERV*KCAL,VWATTS*KCAL);
   fprintf(fp, "%-8.2f%-9.2f%-10.2f",NITV*KCAL,NITNB*KCAL,H2ONITV*KCAL);
   fprintf(fp, "%-9.2f%-9.2f%-10.2f%-10.2f%-10.2f%-7d",
	VINT*KCAL,VINTH2O*KCAL,VINTH2O_s*KCAL,VCH2O[0]*KCAL,VCH2O[1]*KCAL,Nshel);
   fprintf(fp, "%-8.2f%-9.2f%-8.2f%-7.3f%-9.2f%-8.3f%-8.3f",VINT_NIT*KCAL,VC_NIT*KCAL, H2OC*KCAL,pos[natoms-1].fz,NITC*KCAL,V_bias*KCAL,WNITC*KCAL);
   fprintf(fp, "\n");



}
}
