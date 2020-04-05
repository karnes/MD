/*
 *	This routine will write a status line to the file pointed to by fp,
 */

#include	<md.h>
#include	<system.h>
#include	<water.h>

wstat(fp)
FILE *fp;
{

if(nBCD==1){ 
   offsetBCD();
} 

fprintf(fp,"%-8s%-10s%-10s%-11s%-10s%-12s%-9s%-9s%-11s%-10s%-9s%-9s","TIME","TEMP","ETOT","INTRABCD","WATERV","H2ONB","VINT_W","VINT_B","VLIQ","H2OBrOV","BrOV","BrONB");
fprintf(fp,"%-7s%-7s%-7s%-6s%-6s%-9s  %-5s%-6s%-5s%-5s%-5s%-5s\n","BCD_z","BCDHB","BCDWHB","WWHB","POREW","V_window","cos","rBCDz","paC","pC8","pCoM","pBr");

fprintf(fp,"%-8.1f %-8.2f%-11.2f%-11.2f%-9.2f%-12.2f %-9.2f%-9.2f%-10.2f%-10.2f%-9.2f%-7.2f",etime,temp,E*KCAL,INTRABCD*KCAL,WATERV*KCAL,VWATTS*KCAL,(VINTH2O+VINTH2Ooh)*KCAL,VINT_BrO*KCAL,VLIQ*KCAL,H2OBrOV*KCAL,BrOV*KCAL,BrONB*KCAL);
fprintf(fp,"   %-6.2f  %-4d  %-3d   %-5d   %-3d   %-6.3f  %-4.3f  %-4.2f %-5d%-5d%-5d%-5d\n",BCDcom.fz,BCDHB,BCDWHB1+BCDWHB2,WWHB,poreWat,V_window*KCAL,cosccz,rBCDz,poreaC,poreC8,poreCoM,poreBr);

if( (natoms-nsolute-nBCD*BCDs) > 0 && nBCD == 1){
checkGuests();
fprintf(fp,"%-9s %-7s%-7s%-7s%-6s%-6s%-7s%-7s%-7s%-6s%-6s %-6s %-6s%-6s%-7s%-7s%-7s%-7s  %-7s%-7s%-7s%-7s%-7s%-7s\n","TIME","G1","rC8","rCoM","raC","G2","rC8","rCoM","raC","G3","rC8","rCoM","raC","pWat","BCDth","U_PMF","gCoMp","delcos","ETOT","gaCp","gC8p","g_ang","gBias","U_ang");
fprintf(fp,"%-9.1f %-6d%-5.3f %-5.3f %-5.3f   %-4d%-5.3f %-5.3f %-5.3f   %-5d%-5.3f %-5.3f %-5.3f   %-4d  %-6.3f %-7.2f %-5.2f %-6.3f  %-8.2f %-7.2f%-7.2f%-7.2f%-7.2f%-7.2f\n",etime,gst[0],rC8[0],rCoM[0],raC[0],gst[1],rC8[1],rCoM[1],raC[1],gst[2],rC8[2],rCoM[2],raC[2],ph2o,BCDz.fz,U_teth,guestZ,prevBCDcos-BCDz.fz,E*KCAL,gaC,gC8,g_ang,gBias*KCAL,V_gAng);
prevBCDcos = BCDz.fz;

U_teth = 0;
}

if(nsolute==3){
   if(nBCD==0) offsetSN2();
   fprintf(fp,"%-8s%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s","TIME","VINT_WI1","VINT_WI2","VINT_WD1","VINT_WD2","VINT_BI1","VINT_BI2","VINT_BD1","VINT_BD2");
   fprintf(fp,"%-10s%-10s%-10s%-10s%-9s%-9s%-9s%-7s%-6s%-6s%-9s%-9s\n","VINT_CDI1","VINT_CDI2","VINT_CDD1","VINT_CDD2","VINT_EVB","VINT_BCD","V_teth","rBCDsol","Nshel1","Nshel2","VSN2_bCD","VSN2_BCDoh");
   fprintf(fp,"%-8.1f%-9.3f%-9.3f%-9.3f%-9.3f%-9.3f%-9.3f%-9.3f%-9.3f",etime,VINT_WI1*KCAL,VINT_WI2*KCAL,VINT_WD1*KCAL,VINT_WD2*KCAL,VINT_BI1*KCAL,VINT_BI2*KCAL,VINT_BD1*KCAL,VINT_BD2*KCAL);
   fprintf(fp,"%-10.3f%-10.3f%-10.3f%-10.3f%-9.3f%-9.3f%-9.3f%-8.3f%-6d%-6d%-9.2f%-9.2f\n",VINT_CDI1*KCAL,VINT_CDI2*KCAL,VINT_CDD1*KCAL,VINT_CDD2*KCAL,VINT_EVB*KCAL,VINT_BCD*KCAL,V_teth*KCAL,rBCDsol/* *(float)(solBCDside)*/,Nshel[0],Nshel[1],VSN2_BCD*KCAL,VSN2_BCDoh*KCAL);
   fprintf(fp,"%-9s%-9s%-9s%-9s%-9s%-11s%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-6s%-7s%-7s%-7s%-7s\n","TIME","VSYS","r1","r2","C1sq","cosTheta","z0","z1","z2","VRC","FRC","Vbias","VVIB","EKVIB","rVIB","rSN2B","sn2Z","Vsn2w","V_solB");
   if(NW > 0) FRC /= (NW*1.0);
   fprintf(fp,"%-8.1f%-9.2f%-9.3f%-9.3f%-9.3f%-9.2f%-9.2f%-9.2f%-9.2f%-9.2f%-9.2f%-9.3f%-9.3f%-9.3f%-7.3f%-8.3f%-8.3f%-8.3f%-8.3f\n",etime,VSYS*KCAL,RC1,RC2,C1sq,cosTheta,pos[natoms-3].fz, pos[natoms-2].fz, pos[natoms-1].fz,VRC*KCAL/*rxn coord windowing potential*/, FRC*KCAL, Vbias*KCAL,VVIB*KCAL,EKVIB*KCAL,rVIB,sn2BCDrad,sn2z,KCAL*V_sn2w,KCAL*V_solB);
   FRC = 0.0;
   NW = 0.0;
   fprintf(fp,"%-8s%-8s%-8s%-8s%-8s%-8s%-8s%-8s%-8s%-8s\n","TIME","Cl1-aC","Cl2-aC","Cl1-Br","Cl2-Br","Cl1-bOH","Cl2-bOH","Cl1-H2O","Cl2-H2O","s_coord");
   fprintf(fp,"%-8.1f%-8d%-8d%-8d%-8d%-8d%-8d%-8d%-8d%-8.3f\n",etime,Clsol[1][0],Clsol[2][0],Clsol[1][1],Clsol[2][1],Clsol[1][2],Clsol[2][2],Clsol[1][3],Clsol[2][3],KCAL*s_coord);
}
}

checkGuests(){

gst[0] = gst[1] = gst[2] = -1;
rC8[0] = rC8[1] = rC8[2] = 0.0;
raC[0] = raC[1] = raC[2] = 0.0;
rCoM[0] = rCoM[1] = rCoM[2] = 0.0;

float   Bpos[800][2];
int i,j,k,nw;
tripd sdl,com,wat;
double r,totm;

nw = natoms - nBCD*BCDs - nBrO*BrOs - nsolute;

for (i=0;i<nBrO;i++){
   Bpos[i][0] = (float)i;
   sdl.fx = pos[nw+i*BrOs+8].fx - BCDcom.fx;
   sdl.fy = pos[nw+i*BrOs+8].fy - BCDcom.fy;
   sdl.fz = pos[nw+i*BrOs+8].fz - BCDcom.fz;
   mvimage(&sdl);
   r = sqrt(sdl.fx*sdl.fx + sdl.fy*sdl.fy + sdl.fz*sdl.fz);
   Bpos[i][1] = r;
}
for (i= nBrO; i < 800; i++){
   Bpos[i][0] = (float)i;
   Bpos[i][1] = 888.0;
}
quicksort(Bpos,0,800-1);

for(j=0;j<3;j++){
   if(j>=nBrO) continue;
   gst[j] = (int)(Bpos[j][0]);
   rC8[j] = Bpos[j][1];
//sign reveals side of bCD
   sdl.fx = pos[nw+gst[j]*BrOs+8].fx - BCDcom.fx;
   sdl.fy = pos[nw+gst[j]*BrOs+8].fy - BCDcom.fy;
   sdl.fz = pos[nw+gst[j]*BrOs+8].fz - BCDcom.fz;
   mvimage(&sdl);
   rC8[j] *= sgn(sdl.fx*BCDz.fx + sdl.fy*BCDz.fy + sdl.fz*BCDz.fz);
   
   k = nw+gst[j]*BrOs;
   sdl.fx = pos[k].fx - BCDcom.fx;
   sdl.fy = pos[k].fy - BCDcom.fy;
   sdl.fz = pos[k].fz - BCDcom.fz;
   mvimage(&sdl);
   raC[j] = sqrt(sdl.fx*sdl.fx + sdl.fy*sdl.fy + sdl.fz*sdl.fz);
//sign reveals side of bCD
   raC[j] *= sgn(sdl.fx*BCDz.fx + sdl.fy*BCDz.fy + sdl.fz*BCDz.fz);
   com.fx = com.fy = com.fz = totm = 0.0;
   for(i=0;i<BrOs;i++){
      k = nw+gst[j]*BrOs+i;
      com.fx += pos[k].fx*mass[k];
      com.fy += pos[k].fy*mass[k];
      com.fz += pos[k].fz*mass[k];
      totm += mass[k];
   }
   com.fx/=totm;
   com.fy/=totm;
   com.fz/=totm;
   sdl.fx = com.fx - BCDcom.fx;
   sdl.fy = com.fy - BCDcom.fy;
   sdl.fz = com.fz - BCDcom.fz;
   mvimage(&sdl);
   rCoM[j] = sqrt(sdl.fx*sdl.fx + sdl.fy*sdl.fy + sdl.fz*sdl.fz);
//sign reveals side of bCD
   rCoM[j] *= sgn(sdl.fx*BCDz.fx + sdl.fy*BCDz.fy + sdl.fz*BCDz.fz);
}
/* count water in pore */
ph2o = 0;
for(i=0;i<nw;i+=3){
   wat.fx = wat.fy = wat.fz = 0.0;
   /* base on water CoM */
   for(j=0;j<3;j++){
      wat.fx+=pos[i+j].fx*mass[i+j];
      wat.fy+=pos[i+j].fy*mass[i+j];
      wat.fz+=pos[i+j].fz*mass[i+j];
   }
   wat.fx/=WMass;
   wat.fy/=WMass;
   wat.fz/=WMass;
   sdl.fx = wat.fx - BCDcom.fx;
   sdl.fy = wat.fy - BCDcom.fy;
   sdl.fz = wat.fz - BCDcom.fz;
   mvimage(&sdl);
   r = sqrt(sdl.fx*sdl.fx + sdl.fy*sdl.fy + sdl.fz*sdl.fz);
   if(r<3.5) ph2o++;
}
   
}
