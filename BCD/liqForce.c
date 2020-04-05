#include	<md.h>
#include	<system.h>
#include        <water.h>
#include	<string.h>

liqForce(pos, force)
	tripd	*pos;
	tripd	*force;
{
	int i,j,nw;
	double ER, Reflect();
	tc++;
// FOR MOVIE MAKING ONLY
//	if(tc==500){
//	  kRC=0.0;
//	} 
// END MOVIE
	if(tc==0){
	   tdpoint = 0;
//	 tagwat = -1;
	}
	for(i=0;i<3;i++){
	   for(j=0;j<4;j++){
	      Clsol[i][j] = 0;
	   }
	}	   
	
//fprintf(stderr,"liqForce.c: wBCDljq.a = %f, b = %f, q = %f\n",wBCDlj[0][0].a*KCAL, wBCDlj[0][0].b*KCAL,wBCDlj[0][0].q*E2);
//fprintf(stderr,"tdpoint = %d, dt = %d\n",tdpoint,dt);
//	solBCDside = 0;
	poreRad = 3.0;
	topH = botH = 3.9;
	nw = natoms - nBrO*BrOs - nBCD*BCDs - nsolute;
//fprintf(stderr,"liqForce. npoints = %d, tdpoint = %d, natoms = %d, nBCD = %d, nBrO = %d, nw = %d\n",npoints,tdpoint, natoms, nBCD, nBrO, nw);
	Nshel[0] = Nshel[1] = 0;
	INTRABCD = BCDNB = 0.;
	V_window = WATERV = VWATTS = H2OBrOV = BrONB = BrOV = 0.;
	if(nBCD>0){
	   BCDForce();
	}
//fprintf(stderr,"liqForce. after BCDForce call.xwall = %f ywall = %f \n", xwall, ywall);
	if(nw>0){
//	   fprintf(stderr,"liqForce. calling waterForce. nw = %d\n", nw);
	   waterForce(pos, force, nw);
	}
//fprintf(stderr,"liqForce. after waterForce call. H2ONB = %f, H2OV = %f\n", VWATTS*KCAL, WATERV*KCAL);
	if(nBrO>0){
	   BrOForce();
	}
//fprintf(stderr,"liqForce. after BrOForce call. BrONB = %f, BrOV = %f\n", BrONB*KCAL, BrOV*KCAL);
	if(nBrO>0 && nw>0){
	  waterBrO();
	}
//fprintf(stderr,"liqForce. after waterBrO call. H2OBrOV = %f\n", H2OBrOV*KCAL);

	VLIQ = INTRABCD + VWATTS + WATERV + H2OBrOV + BrONB + BrOV;
/* add soft reflecting wall at Zwall-5 to prevent evaporating
 * liquid from entering the other side
 */
	if (nw + nBrO > 0 && pbcType[0] != 'O'){
	   ER = Reflect();
	   VLIQ += ER;
	   if (ER > 0) fprintf(stderr," ER = %f\n",ER*KCAL);
	}
        if(KillFinger==1){
	   extForce(nw,nBrO*BrOs);
	   VLIQ += EXFIELD;
	}
//fprintf(stderr,"liqForce. after VLIQ sum.VLIQ = %f\n", VLIQ*KCAL);
	if(nBCD>0 && action[1]!='E')
	   getRDF(); //separate since relative to BCD CoM

//fprintf(stderr,"liqForce.c: wBCDljq.a = %f, b = %f, q = %f\n",wBCDlj[0][0].a*KCAL, wBCDlj[0][0].b*KCAL,wBCDlj[0][0].q*E2);
}

void getRDF(void){
   int i,j,k,nw,index,index2,xbin,ybin;
   tripd image, sdist, com, rcs;
   double r,r2,dx,dy,dz,totm,costh;
   nw = natoms - nBrO*BrOs - nBCD*BCDs - nsolute;
   if(nw>0){
   for(i=0;i<nw;i+=3){
     // water O - BCD CoM dist 
     image.fx = -(sdist.fx = BCDcom.fx - pos[i].fx);
     image.fy = -(sdist.fy = BCDcom.fy - pos[i].fy);
     image.fz = -(sdist.fz = BCDcom.fz - pos[i].fz);
     mvimage(&sdist);
     image.fx += sdist.fx;
     image.fy += sdist.fy;
     image.fz += sdist.fz;
     r2 = sdist.fx*sdist.fx + sdist.fy*sdist.fy + sdist.fz*sdist.fz;
     r = sqrt(r2);
     index = (int)(r/binRDF);
     if(index < RDFbins){
        grSOL[0][index] += 1.0; //BCD com - water O
        com.fx = com.fy = com.fz = totm = 0.0;
        for(j=0;j<3;j++){
 	   com.fx += pos[i+j].fx * mass[i+j];
	   com.fy += pos[i+j].fy * mass[i+j];
	   com.fz += pos[i+j].fz * mass[i+j];
	   totm += mass[i+j];
        }
        com.fx /= totm;
        com.fy /= totm;
        com.fz /= totm;
        rcs.fx = BCDcom.fx - com.fx + image.fx;
        rcs.fy = BCDcom.fy - com.fy + image.fy;
        rcs.fz = BCDcom.fz - com.fz + image.fz;
        r2 = rcs.fx*rcs.fx + rcs.fy*rcs.fy + rcs.fz*rcs.fz;
        r = sqrt(r2);
        index2 = (int)(r/binRDF);
        if(index2 < RDFbins)
	     grSOL[1][index2] += 1.0; //BCD com - water com
	// g(r,theta)
	if(r < contRad){
	   costh = BCDz.fx*rcs.fx+BCDz.fy*rcs.fy+BCDz.fz*rcs.fz;
	   costh /= r;
	   dy = r*costh;
	   dx = r*sqrt(1-(costh*costh));
	   ybin = (int)(round((dy+contRad)/contBin));
	   xbin = (int)(dx/contBin);
	   BWmap[ybin][xbin]++;
/*	   costh = BCDxy.fx*rcs.fx+BCDxy.fy*rcs.fy+BCDxy.fz*rcs.fz;
	   costh /= r;
	   dy = r*costh;
	   dx = r*sqrt(1-(costh*costh));
	   ybin = (int)(round((dy+contRad)/contBin));
	   xbin = (int)(dx/contBin);
	   BWmapxy[ybin][xbin]++;*/
	}
	dz = rcs.fx*BCDz.fx + rcs.fy*BCDz.fy + rcs.fz*BCDz.fz;
	if(fabs(dz) < xyMapz){
	   dx = rcs.fx*BCDx.fx + rcs.fy*BCDx.fy + rcs.fz*BCDx.fz;
	   if(fabs(dx) < xyMapx){
	      dy = rcs.fx*BCDy.fx + rcs.fy*BCDy.fy + rcs.fz*BCDy.fz;
	      if(fabs(dy) < xyMapy){
		  BWmapxy[(int)(round((dy+xyMapy)/contBin))][(int)(round((dx+xyMapx)/contBin))]++;
	      }
	   }
	}

     }
   }
   }
   if(nBrO > 0){
   for(k=0;k<nBrO;k++){
     i = nw + k*BrOs;
     // BrO C - BCD CoM dist 
     image.fx = -(sdist.fx = BCDcom.fx - pos[i].fx);
     image.fy = -(sdist.fy = BCDcom.fy - pos[i].fy);
     image.fz = -(sdist.fz = BCDcom.fz - pos[i].fz);
     mvimage(&sdist);
     image.fx += sdist.fx;
     image.fy += sdist.fy;
     image.fz += sdist.fz;
     r2 = sdist.fx*sdist.fx + sdist.fy*sdist.fy + sdist.fz*sdist.fz;
     r = sqrt(r2);
     // g(r,theta)
     if(r < contRad){
	costh = BCDz.fx*sdist.fx+BCDz.fy*sdist.fy+BCDz.fz*sdist.fz;
	costh /= r;
	dy = r*costh;
	dx = r*sqrt(1-(costh*costh));
	ybin = (int)(round((dy+contRad)/contBin));
	xbin = (int)(dx/contBin);
	BaCmap[ybin][xbin]++;
     }
     rcs.fx = sdist.fx; 
     rcs.fy = sdist.fy;
     rcs.fz = sdist.fz;
     dz = rcs.fx*BCDz.fx + rcs.fy*BCDz.fy + rcs.fz*BCDz.fz;
     if(fabs(dz) < xyMapz){
	   dx = rcs.fx*BCDx.fx + rcs.fy*BCDx.fy + rcs.fz*BCDx.fz;
	   if(fabs(dx) < xyMapx){
	      dy = rcs.fx*BCDy.fx + rcs.fy*BCDy.fy + rcs.fz*BCDy.fz;
	      if(fabs(dy) < xyMapy){
		  BaCmapxy[(int)(round((dy+xyMapy)/contBin))][(int)(round((dx+xyMapx)/contBin))]++;
	      }
	   }
     }
     index = (int)(r/binRDF);
     if(index < RDFbins)
        grSOL[3][index] += 1.0; //BCD com - BrO alpha C
     com.fx = com.fy = com.fz = totm = 0.0;
     for(j=0;j<BrOs;j++){
        com.fx += pos[i+j].fx * mass[i+j];
	com.fy += pos[i+j].fy * mass[i+j];
	com.fz += pos[i+j].fz * mass[i+j];
	totm += mass[i+j];
     }
     com.fx /= totm;
     com.fy /= totm;
     com.fz /= totm;
     rcs.fx = BCDcom.fx - com.fx + image.fx;
     rcs.fy = BCDcom.fy - com.fy + image.fy;
     rcs.fz = BCDcom.fz - com.fz + image.fz;
     r2 = rcs.fx*rcs.fx + rcs.fy*rcs.fy + rcs.fz*rcs.fz;
     r = sqrt(r2);
     // g(r,theta)
     if(r < contRad){
	costh = BCDz.fx*rcs.fx+BCDz.fy*rcs.fy+BCDz.fz*rcs.fz;
	costh /= r;
	dy = r*costh;
	dx = r*sqrt(1-(costh*costh));
	ybin = (int)(round((dy+contRad)/contBin));
	xbin = (int)(dx/contBin);
	BBcommap[ybin][xbin]++;
     }
     dz = rcs.fx*BCDz.fx + rcs.fy*BCDz.fy + rcs.fz*BCDz.fz;
     if(fabs(dz) < xyMapz){
	   dx = rcs.fx*BCDx.fx + rcs.fy*BCDx.fy + rcs.fz*BCDx.fz;
	   if(fabs(dx) < xyMapx){
	      dy = rcs.fx*BCDy.fx + rcs.fy*BCDy.fy + rcs.fz*BCDy.fz;
	      if(fabs(dy) < xyMapy){
		  BBcommapxy[(int)(round((dy+xyMapy)/contBin))][(int)(round((dx+xyMapx)/contBin))]++;
	      }
	   }
     }
     index2 = (int)(r/binRDF);
     if(index2 < RDFbins)
	  grSOL[4][index2] += 1.0; //BCD com - BrO com
     rcs.fx = BCDcom.fx - pos[i+1].fx + image.fx;
     rcs.fy = BCDcom.fy - pos[i+1].fy + image.fy;
     rcs.fz = BCDcom.fz - pos[i+1].fz + image.fz;
     r2 = rcs.fx*rcs.fx + rcs.fy*rcs.fy + rcs.fz*rcs.fz;
     r = sqrt(r2);
     // g(r,theta)
     if(r < contRad){
	costh = BCDz.fx*rcs.fx+BCDz.fy*rcs.fy+BCDz.fz*rcs.fz;
	costh /= r;
	dy = r*costh;
	dx = r*sqrt(1-(costh*costh));
	ybin = (int)(round((dy+contRad)/contBin));
	xbin = (int)(dx/contBin);
	BBrmap[ybin][xbin]++;
     }
     dz = rcs.fx*BCDz.fx + rcs.fy*BCDz.fy + rcs.fz*BCDz.fz;
     if(fabs(dz) < xyMapz){
	   dx = rcs.fx*BCDx.fx + rcs.fy*BCDx.fy + rcs.fz*BCDx.fz;
	   if(fabs(dx) < xyMapx){
	      dy = rcs.fx*BCDy.fx + rcs.fy*BCDy.fy + rcs.fz*BCDy.fz;
	      if(fabs(dy) < xyMapy){
		  BBrmapxy[(int)(round((dy+xyMapy)/contBin))][(int)(round((dx+xyMapx)/contBin))]++;
	      }
	   }
     }
     index2 = (int)(r/binRDF);
     if(index2 < RDFbins)
	  grSOL[5][index2] += 1.0; //BCD com  -BrO Br
     rcs.fx = BCDcom.fx - pos[i+BrOs-1].fx + image.fx;
     rcs.fy = BCDcom.fy - pos[i+BrOs-1].fy + image.fy;
     rcs.fz = BCDcom.fz - pos[i+BrOs-1].fz + image.fz;
     r2 = rcs.fx*rcs.fx + rcs.fy*rcs.fy + rcs.fz*rcs.fz;
     r = sqrt(r2);
     // g(r,theta)
     if(r < contRad){
	costh = BCDz.fx*rcs.fx+BCDz.fy*rcs.fy+BCDz.fz*rcs.fz;
	costh /= r;
	dy = r*costh;
	dx = r*sqrt(1-(costh*costh));
	ybin = (int)(round((dy+contRad)/contBin));
	xbin = (int)(dx/contBin);
	BC8map[ybin][xbin]++;
     }
     dz = rcs.fx*BCDz.fx + rcs.fy*BCDz.fy + rcs.fz*BCDz.fz;
     if(fabs(dz) < xyMapz){
	   dx = rcs.fx*BCDx.fx + rcs.fy*BCDx.fy + rcs.fz*BCDx.fz;
	   if(fabs(dx) < xyMapx){
	      dy = rcs.fx*BCDy.fx + rcs.fy*BCDy.fy + rcs.fz*BCDy.fz;
	      if(fabs(dy) < xyMapy){
		  BC8mapxy[(int)(round((dy+xyMapy)/contBin))][(int)(round((dx+xyMapx)/contBin))]++;
	      }
	   }
     }
     index2 = (int)(r/binRDF);
     if(index2 < RDFbins)
          grSOL[6][index2] += 1.0; //BCD com - BrO C8(tail)
     
   }
   }
     
 
}
/* reflecting wall to keep solvents from crossing over to other sides of box */
#define ExCon 0.5
double
Reflect()
{
int i,j,nw;
double cmposZ,ztop,zbot;
double RFIELD;
RFIELD = 0.;
nw = natoms - nBCD*BCDs - nBrO*BrOs - nsolute;
zbot = -zwall + 5.0;/*reflecting wall 5A from the box edge*/
ztop =  zwall - 5.0;/*reflecting wall 5A from the box edge*/

for(i = 0; i < nw; i=i+3){
   cmposZ = 0.;
   for(j=0;j<3;j++)
      cmposZ += mass[i+j]*pos[i+j].fz;
   cmposZ /= WMass;
   if(cmposZ < zbot) {
      fprintf(stderr,"water %d is above. zcom = %f.\n",i,cmposZ);
      RFIELD += ExCon*(cmposZ-zbot)*(cmposZ-zbot);
      for(j=0;j<3;j++)
         force[i+j].fz -= 2*ExCon*mass[i+j]*(cmposZ-zbot)/WMass;
   }
   else if(cmposZ > ztop) {
      fprintf(stderr,"water %d is below. zcom = %f.\n",i,cmposZ);
      RFIELD += ExCon*(cmposZ-ztop)*(cmposZ-ztop);
      for(j=0;j<3;j++)
         force[i+j].fz -= 2*ExCon*mass[i+j]*(cmposZ-ztop)/WMass;
   }
}
for(i = 0; i < nBrO*BrOs; i=i+BrOs){
   cmposZ = 0.;
   for(j=0;j<BrOs;j++)
      cmposZ += mass[nw+i+j]*pos[nw+i+j].fz;
   cmposZ /= BMass;
   if(cmposZ > ztop){
      fprintf(stderr,"Br-oct %d is above. zcom = %f.\n",i,cmposZ);
      RFIELD += ExCon*(cmposZ-ztop)*(cmposZ-ztop);
      for(j=0;j<BrOs;j++)
         force[nw+i+j].fz -= 2*ExCon*mass[nw+i+j]*(cmposZ-ztop)/BMass;
   }
   else if(cmposZ < zbot){
      fprintf(stderr,"Br-oct %d is below. zcom = %f.\n",i,cmposZ);
      RFIELD += ExCon*(cmposZ-zbot)*(cmposZ-zbot);
      for(j=0;j<BrOs;j++)
         force[nw+i+j].fz -= 2*ExCon*mass[nw+i+j]*(cmposZ-zbot)/BMass;
   }
}
return(RFIELD);
}
/* calculate the external forces and energy that remove surface roghness *
 * We assume that the water is in the z<0 region and the Br-octane	 *
 * in the z > 0 region.							 */
extForce(nw,nb)
	int nw;/* number of water atoms*/
	int nb;/* number of Br-octane atoms*/
{
int i,j;
double cmposZ;
EXFIELD = 0.;

for(i=0;i<nw;i=i+3){
    cmposZ = 0.;
    for (j=0;j<3;j++)
 	cmposZ += mass[i+j]*pos[i+j].fz;
    if (cmposZ > 0) {/* water in the z>0 side, force it to the z<0 side */
		cmposZ /= WMass;
		EXFIELD += ExCon*cmposZ*cmposZ;
		for (j=0;j<3;j++)
			force[i+j].fz -= 2*ExCon*mass[i+j]*cmposZ/WMass;
    }
    
}
for(i=0;i<nb;i+=BrOs){
	cmposZ = 0.;
	for (j=0;j<BrOs;j++)
		cmposZ += mass[nw+i+j]*pos[nw+i+j].fz;
	if (cmposZ < 0) {/* NIT in the z<0 side, force it to the z>0 side */
		cmposZ /= BMass;
		EXFIELD += ExCon*cmposZ*cmposZ;
		for (j=0;j<BrOs;j++)
			force[nw+i+j].fz -= 2*ExCon*mass[nw+i+j]*cmposZ/BMass;
	}
}
}
