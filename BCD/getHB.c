#include <md.h>
#include <system.h>
#include <math.h>
#include <water.h>

//int inPore(tripd);

getHB(){
int i,j,k,l,m,nw,ns;
int o1, o2,over;
double totm,r2,r,h,rcc,raC_8p;
double doti1, doti2, dotj1, dotj2;
tripd sdist, image, wat, wat2, C8Br;
nw = natoms - nBCD*BCDs - nBrO*BrOs - nsolute;
ns = natoms - nBCD*BCDs - nsolute;
//fprintf(stderr,"in getHB.c ... ns = %d\n",ns);
BCDHB = WWHB = BCDWHB1 = BCDWHB2 = 0;
if(tdpoint==0){
   priD = priA = secD = secA = watHBtot = 0;
}
cosgz = 2.0;
rcc = 99.0;
/* water-water HBs */
if(nw > 0){
for(i=0;i<nw;i+=3){
   for(j=i+3;j<nw;j+=3){
      /* get O - O distance */
      image.fx = -(sdist.fx = pos[i].fx - pos[j].fx);
      image.fy = -(sdist.fy = pos[i].fy - pos[j].fy);
      image.fz = -(sdist.fz = pos[i].fz - pos[j].fz);
      mvimage(&sdist);
      image.fx += sdist.fx;
      image.fy += sdist.fy;
      image.fz += sdist.fz;
      r2 = sdist.fx*sdist.fx + sdist.fy*sdist.fy + sdist.fz*sdist.fz;
      if(r2 < rOOmax*rOOmax){
         checkWWHB(i,j,image,sdist,sqrt(r2));
      }
   }
}
}
//fprintf(stderr,"done with WW...\n");
if(nBCD>0){
/* BCD-BCD HBs */
for(l=0;l<7;l++){
   if(l==6)
	over = 0;
   else
	over = (l+1)*21;
   /* intra-glucose unit */
   /* get O - O distance */
/* o1 = ns + l*21 + 19;
   o2 = ns + l*21 + 17;
   image.fx = -(sdist.fx = pos[o1].fx - pos[o2].fx);
   image.fy = -(sdist.fy = pos[o1].fy - pos[o2].fy);
   image.fz = -(sdist.fz = pos[o1].fz - pos[o2].fz);
   mvimage(&sdist);
   image.fx += sdist.fx;
   image.fy += sdist.fy;
   image.fz += sdist.fz;
   r2 = sdist.fx*sdist.fx + sdist.fy*sdist.fy + sdist.fz*sdist.fz;
   if(r2 < rOOmax*rOOmax){
      checkBCDHB(o1,o2,image,sdist,sqrt(r2),3);
   }
*/
   /* inter-glucose unit */
   /* get O-O distance (secondary)*/
   o1 = ns + over + 19;
   o2 = ns + l*21 + 17;
   image.fx = -(sdist.fx = pos[o1].fx - pos[o2].fx);
   image.fy = -(sdist.fy = pos[o1].fy - pos[o2].fy);
   image.fz = -(sdist.fz = pos[o1].fz - pos[o2].fz);
   mvimage(&sdist);
   image.fx += sdist.fx;
   image.fy += sdist.fy;
   image.fz += sdist.fz;
   r2 = sdist.fx*sdist.fx + sdist.fy*sdist.fy + sdist.fz*sdist.fz;
   if(r2 < rOOmax*rOOmax){
      checkBCDHB(o1,o2,image,sdist,sqrt(r2),2);
   }
   /* get O-O distance (primary)*/
/* o1 = ns + over + 14;
   o2 = ns + l*21 + 14;
   image.fx = -(sdist.fx = pos[o1].fx - pos[o2].fx);
   image.fy = -(sdist.fy = pos[o1].fy - pos[o2].fy);
   image.fz = -(sdist.fz = pos[o1].fz - pos[o2].fz);
   mvimage(&sdist);
   image.fx += sdist.fx;
   image.fy += sdist.fy;
   image.fz += sdist.fz;
   r2 = sdist.fx*sdist.fx + sdist.fy*sdist.fy + sdist.fz*sdist.fz;
   if(r2 < rOOmax*rOOmax){
      checkBCDHB(o1,o2,image,sdist,sqrt(r2),1);
   }
   */
}
//fprintf(stderr,"done with BCD-BCD...\n");
if(nw>0){
/* BCD - water HBs */
for(l=0;l<7;l++){
   for(m=14;m<20;m+=2){
      if(m==16) m=17;
      o1 = ns + (l*21) + m;
      for(i=0;i<nw;i+=3){
	 image.fx = -(sdist.fx = pos[o1].fx - pos[i].fx);
	 image.fy = -(sdist.fy = pos[o1].fy - pos[i].fy);
	 image.fz = -(sdist.fz = pos[o1].fz - pos[i].fz);
	 mvimage(&sdist);
         image.fx += sdist.fx;
         image.fy += sdist.fy;
         image.fz += sdist.fz;
         r2 = sdist.fx*sdist.fx + sdist.fy*sdist.fy + sdist.fz*sdist.fz;
	 if(r2 < rOOmax*rOOmax){
	    checkBCDWHB(o1,i,image,sdist,sqrt(r2),m);
	 }
      }
   }
}
//fprintf(stderr,"done with BCD-W...\n");

/* count water in pore */
poreWat = 0;
for(i=0;i<maxPoreWat+1;i++){
   poreW[i] = -1;
}
for(i=0;i<nw;i+=3){
   wat.fx = wat.fy = wat.fz = totm = 0.0;
   /* base on water CoM */
   for(j=0;j<3;j++){
      wat.fx+=pos[i+j].fx*mass[i+j];
      wat.fy+=pos[i+j].fy*mass[i+j];
      wat.fz+=pos[i+j].fz*mass[i+j];
   }
   wat.fx/=WMass;
   wat.fy/=WMass;
   wat.fz/=WMass;
// if(i==0) fprintf(stderr,"inPore(water %d) = %d, tc = %d, dataRate = %d, tdpoint = %d\n",i/3,inPore(wat),tc,dataRatex,tdpoint);
   watDipVec[i/3][tdpoint].inp = inPore(wat);   
   if(watDipVec[i/3][tdpoint].inp==1){
// NEW FOR TESTING
//     if(tagwat==-1){ tagwat = i; fprintf(stderr,"tagged water %d.\n",i); }
//END NEW FOR TESTING

//   fprintf(stderr,"water in pore. (water %d)\n",i/3);
     if(poreWat<=maxPoreWat){
        poreW[poreWat]=i/3;
        poreWat++;
     }
     else
     {
	  fprintf(stderr,"Warning: detected more than %d water molecules in pore.\n",maxPoreWat);
     }
   }
}
probPoreWat[poreWat]++;
probPoreWat[maxPoreWat+1]++; //total for normalization
}//end 'if nw>0' loop
if(nBrO>0){
poreaC = poreC8 = poreCoM = poreBr = 0; 
/* count BrO in pore */
for(j=0;j<4;j++){
   poreBrO[j] = 0;
   for(i=0;i<maxPoreB+1;i++){
     poreB[j][i] = -1;
   }
}
for(k=0;k<nBrO;k++){
   i=nw+k*BrOs;
   wat.fx = wat.fy = wat.fz = totm = 0.0;
   /* check BrO CoM */
   for(j=0;j<BrOs;j++){
      wat.fx+=pos[i+j].fx*mass[i+j];
      wat.fy+=pos[i+j].fy*mass[i+j];
      wat.fz+=pos[i+j].fz*mass[i+j];
   }
   wat.fx/=BMass;
   wat.fy/=BMass;
   wat.fz/=BMass;
   if(BrOVec[2][k][tdpoint].inp = inPore(&wat)){
        if(poreBrO[0]<=maxPoreB){
	  poreB[0][poreBrO[0]]=i;
	  poreBrO[0]++;
	  poreCoM++;
	}
	else
	{
	  fprintf(stderr,"Warning: detected more than %d Br-octane CoM in pore.\n",maxPoreB);
	}
   }
   /* check alpha C atom */
   if(BrOVec[0][k][tdpoint].inp = inPore(pos[i])){
        if(poreBrO[1]<=maxPoreB){
	  poreB[1][poreBrO[1]]=i;
	  poreBrO[1]++;
	  poreaC++;
	}
	else
	{
	  fprintf(stderr,"Warning: detected more than %d Br-octane alpha-C in pore.\n",maxPoreB);
	}
   }
   /* check Br atom */
   if(inPore(pos[i+1])){
        if(poreBrO[2]<=maxPoreB){
	  poreB[2][poreBrO[2]]=i;
	  poreBrO[2]++;
	  poreBr++;
	}
	else
	{
	  fprintf(stderr,"Warning: detected more than %d Br-octane Br in pore.\n",maxPoreB);
	}
   }
   /* check C8 atom */
   if(BrOVec[1][k][tdpoint].inp = inPore(pos[i+8])){
        if(poreBrO[3]<=maxPoreB){
	  poreB[3][poreBrO[3]]=i;
	  poreBrO[3]++;
	  poreC8++;
	  // is this closest C8 to CoM?
	  r = 0.0;
	  r += sq(pos[i+8].fx - BCDcom.fx);
	  r += sq(pos[i+8].fy - BCDcom.fy);
	  r += sq(pos[i+8].fz - BCDcom.fz);
	  r = sqrt(r);
	  if(fabs(rcc) > r){
	     rcc = r;
	     r = 0.0;
	     r += sq(pos[i+0].fx - BCDcom.fx);
	     r += sq(pos[i+0].fy - BCDcom.fy);
	     r += sq(pos[i+0].fz - BCDcom.fz);
	     raC_8p = sqrt(r);
	     C8Br.fx = pos[i+0].fx - pos[i+8].fx;
	     C8Br.fy = pos[i+0].fy - pos[i+8].fy;
	     C8Br.fz = pos[i+0].fz - pos[i+8].fz;
	     mvimage(&C8Br);
	     r = sqrt(sq(C8Br.fx)+sq(C8Br.fy)+sq(C8Br.fz));
	     C8Br.fx /= r;
	     C8Br.fy /= r;
	     C8Br.fz /= r;
	     cosgz = BCDz.fx*C8Br.fx + BCDz.fy*C8Br.fy + BCDz.fz*C8Br.fz;
	     raC_8p *= sgn(cosgz);
	  }
	}
	else
	{
	  fprintf(stderr,"Warning: detected more than %d Br-octane C8 in pore.\n",maxPoreB);
	}
   }
}
// if a C8 in pore, tally statistics
if(cosgz != 2.0){
   BrOgOD[(int)((1.0+cosgz)/sODbin)]++;
   BrOgODnorm++;
   aCBCDr[(int)((raC_8p+maxBrComr)/0.5)]++;
}

for(i=0;i<4;i++){
   probPoreBrO[i][poreBrO[i]]++;
   probPoreBrO[i][maxPoreB+1]++; //total for normalization
}
}//end 'if BrO>0' loop

}//end 'if BCD>0' loop
}


int inPore(tripd wat){

double r2,r,h,guest;
tripd wat2,sdist;
guest = 0;

wat.fx-=BCDcom.fx;
wat.fy-=BCDcom.fy;
wat.fz-=BCDcom.fz;
mvimage(&wat);
/* "height" from CoM along BCD vector*/
h = fabs(wat.fx*BCDz.fx + wat.fy*BCDz.fy + wat.fz*BCDz.fz);
if(h<topH || h<botH){
   /* "radius" from BCDz vector */
   wat2.fx = -(BCDz.fx - wat.fx);
   wat2.fy = -(BCDz.fy - wat.fy);
   wat2.fz = -(BCDz.fz - wat.fz);
   sdist.fx = wat.fy*wat2.fz - wat.fz*wat2.fy;
   sdist.fy = wat.fz*wat2.fx - wat.fx*wat2.fz;
   sdist.fz = wat.fx*wat2.fy - wat.fy*wat2.fx;
   r = sqrt(sdist.fx*sdist.fx + sdist.fy*sdist.fy + sdist.fz*sdist.fz);    
   if(r<poreRad){
      guest = 1;
   }
}
return guest;
}



checkBCDWHB(int o1, int i, tripd image, tripd oodist, double rOO,int m){
   double dotB, dot1, dot2, rohB, roh1, roh2;
   tripd ohB, oh1, oh2, oBo1, o1oB;
   int ns,unit,index;

   ns = natoms - nBCD*BCDs - nsolute;
   unit = (o1-ns-m)/21;
   if(m==14){ // primary
      index = unit;
   }
   else if(m==17){
      index = unit + 7;
   }
   else if(m==19){ // secondary, either side
      index = unit + 14;
   }
   else{
      fprintf(stderr, "getHB: bad index passed to HB detection routine.\n");
   }
   /* BCD O -> BCD H */
   ohB.fx = pos[o1+1].fx - pos[o1].fx;
   ohB.fy = pos[o1+1].fy - pos[o1].fy;
   ohB.fz = pos[o1+1].fz - pos[o1].fz;
      
   /* BCD O -> w O */
   oBo1.fx = - oodist.fx;
   oBo1.fy = - oodist.fy;
   oBo1.fz = - oodist.fz;

   /* w O -> 1 H */
   oh1.fx = pos[i+1].fx - pos[i].fx;
   oh1.fy = pos[i+1].fy - pos[i].fy;
   oh1.fz = pos[i+1].fz - pos[i].fz;
      
   /* w O -> 2 H */
   oh2.fx = pos[i+2].fx - pos[i].fx;
   oh2.fy = pos[i+2].fy - pos[i].fy;
   oh2.fz = pos[i+2].fz - pos[i].fz;
      
   /* w O -> BCD O */
   o1oB.fx = oodist.fx;
   o1oB.fy = oodist.fy;
   o1oB.fz = oodist.fz;

   dotB = oBo1.fx*ohB.fx + oBo1.fy*ohB.fy + oBo1.fz*ohB.fz;
   rohB = sqrt(ohB.fx*ohB.fx + ohB.fy*ohB.fy + ohB.fz*ohB.fz);
   dot1 = o1oB.fx*oh1.fx + o1oB.fy*oh1.fy + o1oB.fz*oh1.fz;
   roh1 = sqrt(oh1.fx*oh1.fx + oh1.fy*oh1.fy + oh1.fz*oh1.fz);
   dot2 = o1oB.fx*oh2.fx + o1oB.fy*oh2.fy + o1oB.fz*oh2.fz;
   roh2 = sqrt(oh2.fx*oh2.fx + oh2.fy*oh2.fy + oh2.fz*oh2.fz);

   if(dotB/(rOO*rohB) > cosHOO){
	/* BCD H is donor */
	if(m==14){
		BCDWHB1++; // (1 or 2 denote primary or secondary hydroxyl of BCD)
		priD++;
		HB[index][tdpoint] = i/3;
	}
	else{
		BCDWHB2++;
		secD++;
		HB[index][tdpoint] = i/3;
	}
   }
   if(dot1/(rOO*roh1) > cosHOO){
	/* water H1 is donor */
	if(m==14){
		BCDWHB1++;
		priA++;
		if(HB[index+21][tdpoint] == -1) HB[index+21][tdpoint] = i/3;
		else HB[index+21+21][tdpoint] = i/3;
	}
	else{
		BCDWHB2++;
		secA++;
		if(HB[index+21][tdpoint] == -1) HB[index+21][tdpoint] = i/3;
		else HB[index+21+21][tdpoint] = i/3;
	}
   }
   if(dot2/(rOO*roh2) > cosHOO){
        /* water H2 is donor */
	if(m==14){
		BCDWHB1++;
		priA++;
		if(HB[index+21][tdpoint] == -1) HB[index+21][tdpoint] = i/3;
		else HB[index+21+21][tdpoint] = i/3;
	}
	else{
		BCDWHB2++;
		secA++;
		if(HB[index+21][tdpoint] == -1) HB[index+21][tdpoint] = i/3;
		else HB[index+21+21][tdpoint] = i/3;
	}
   }
   

}

checkBCDHB(int o1, int o2, tripd image, tripd oodist, double rOO,int deg){
   double dot1, dot2, roh1, roh2;
   tripd oh1, oh2, o1o2, o2o1;
   int ns,unit1,unit2,index1,index2,m,n;
   ns = natoms - nBCD*BCDs - nsolute;
   unit1 = (o1 - ns)/21;
   unit2 = (o2 - ns)/21;
   m = o1 - ns - 21*unit1;
   n = o2 - ns - 21*unit2;
   if(m==14){ // primary
      index1 = unit1;
   }
   else if(m==17){
      index1 = unit1 + 7;
   }
   else if(m==19){ // secondary, ether side
      index1 = unit1 + 14;
   }
   else{
      fprintf(stderr, "getHB: bad index passed to HB detection routine.\n");
   }
   if(n==14){ // primary
      index2 = unit2;
   }
   else if(n==17){
      index2 = unit2 + 7;
   }
   else if(n==19){ // secondary, ether side
      index2 = unit2 + 14;
   }
   else{
      fprintf(stderr, "getHB: bad index passed to HB detection routine.\n");
   }

   /* 1 O -> H */
   oh1.fx = pos[o1+1].fx - pos[o1].fx;
   oh1.fy = pos[o1+1].fy - pos[o1].fy;
   oh1.fz = pos[o1+1].fz - pos[o1].fz;

   /* 1 O -> 2 O */
   o1o2.fx = -oodist.fx;
   o1o2.fy = -oodist.fy;
   o1o2.fz = -oodist.fz;

   /* 2 O -> H */
   oh2.fx = pos[o2+1].fx - pos[o2].fx;
   oh2.fy = pos[o2+1].fy - pos[o2].fy;
   oh2.fz = pos[o2+1].fz - pos[o2].fz;
      
   /* 2 O -> 1 O */
   o2o1.fx = oodist.fx;
   o2o1.fy = oodist.fy;
   o2o1.fz = oodist.fz;

   dot1 = oh1.fx*o1o2.fx + oh1.fy*o1o2.fy + oh1.fz*o1o2.fz;
   roh1 = sqrt(oh1.fx*oh1.fx + oh1.fy*oh1.fy + oh1.fz*oh1.fz);
   dot2 = oh2.fx*o2o1.fx + oh2.fy*o2o1.fy + oh2.fz*o2o1.fz;
   roh2 = sqrt(oh2.fx*oh2.fx + oh2.fy*oh2.fy + oh2.fz*oh2.fz);

   if(dot1/(rOO*roh1) > cosHOO){
      /* hydrogen o1 hydroxyl is donor */
      if(deg==2){
	 BCDHB++;
         HB[index1][tdpoint] = o2;
         if(HB[index2+21][tdpoint] == -1) HB[index2+21][tdpoint] = o1;
         else HB[index2+21+21][tdpoint] = o1;
      }
      else if(deg==1){
	 priBCDHB++;
      }
      else if(deg==3){
	 iuBCDHB++;
      }

   }
   if(dot2/(rOO*roh2) > cosHOO){
      /* hydrogen o2 hydroxyl is donor */
      if(deg==2){
	 BCDHB++;
         if(HB[index1+21][tdpoint] == -1) HB[index1+21][tdpoint] = o2;
         else HB[index1+21+21][tdpoint] = o2;
         HB[index2][tdpoint] = o1;
      }
      else if(deg==1){
	 priBCDHB++;
      }
      else if(deg==3){
	 iuBCDHB++;
      }
   }
}

checkWWHB(int i, int j, tripd image, tripd oodist, double rOO){
   double doti1, doti2, dotj1, dotj2;
   double ri1, ri2, rj1, rj2;
   tripd iOjO,jOiO,i1,i2,j1,j2;
   
   /* i O -> H */
   i1.fx = pos[i+1].fx - pos[i].fx;
   i1.fy = pos[i+1].fy - pos[i].fy;
   i1.fz = pos[i+1].fz - pos[i].fz;
   i2.fx = pos[i+2].fx - pos[i].fx;
   i2.fy = pos[i+2].fy - pos[i].fy;
   i2.fz = pos[i+2].fz - pos[i].fz;
   
   /* i O -> j O */
   iOjO.fx = - oodist.fx;	
   iOjO.fy = - oodist.fy;	
   iOjO.fz = - oodist.fz;	

   /* j O -> H */
   j1.fx = pos[j+1].fx - pos[j].fx;
   j1.fy = pos[j+1].fy - pos[j].fy;
   j1.fz = pos[j+1].fz - pos[j].fz;
   j2.fx = pos[j+2].fx - pos[j].fx;
   j2.fy = pos[j+2].fy - pos[j].fy;
   j2.fz = pos[j+2].fz - pos[j].fz;

   /* j O -> i O */
   jOiO.fx = oodist.fx;	
   jOiO.fy = oodist.fy;	
   jOiO.fz = oodist.fz;

   doti1 = i1.fx*iOjO.fx + i1.fy*iOjO.fy + i1.fz*iOjO.fz;
   ri1 = sqrt(i1.fx*i1.fx + i1.fy*i1.fy + i1.fz*i1.fz);
   doti2 = i2.fx*iOjO.fx + i2.fy*iOjO.fy + i2.fz*iOjO.fz;
   ri2 = sqrt(i2.fx*i2.fx + i2.fy*i2.fy + i2.fz*i2.fz);

   dotj1 = j1.fx*jOiO.fx + j1.fy*jOiO.fy + j1.fz*jOiO.fz;
   rj1 = sqrt(j1.fx*j1.fx + j1.fy*j1.fy + j1.fz*j1.fz);
   dotj2 = j2.fx*jOiO.fx + j2.fy*jOiO.fy + j2.fz*jOiO.fz;
   rj2 = sqrt(j2.fx*j2.fx + j2.fy*j2.fy + j2.fz*j2.fz);
   if(doti1/(rOO*ri1) > cosHOO){
      /* hydrogen i1 in W-W HB */
      WWHB++;
      watHBtot++;
   }
   if(doti2/(rOO*ri2) > cosHOO){
      /* hydrogen i2 in W-W HB */
      WWHB++;
      watHBtot++;
   }
   if(dotj1/(rOO*rj1) > cosHOO){
      /* hydrogen j1 in W-W HB */
      WWHB++;
      watHBtot++;
   }
   if(dotj2/(rOO*rj2) > cosHOO){
      /* hydrogen j2 in W-W HB */
      WWHB++;
      watHBtot++;
   }

}
