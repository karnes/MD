#include	<md.h>
#include	<system.h>
#include	<water.h>
#include	<math.h>
/* non-bonded forces between the b-cyclodextrin and the H2O/bromooctane (BrO)
 * liquids */

intForce(pos, force)
	tripd	*pos;
	tripd	*force;
{
	int i,j,l,m,b,nw,nb,nsolvent;
	double /*c,*/ pe; //, ElecPot();
//return (0);
	VINT = VINTH2O = VINTH2Ooh = VCH2O = 0.;/*total interaction with all water*/
	VC_BrO = VINT_BrO = 0.;
	VINT_EVB = VINT_BCD = VSN2_BCD = VSN2_BCDoh = 0.0;
	rBCDsol = EKVIB = VVIB = 0.0;
	if (nsolute==3) solForce(pos,force);
//	if (0 && tagwat!=-1) watEC();
//fprintf(stderr,"intForce.c after solForce call. natoms = %d, nsolute = %d, nw = %d swr2max = %f \n",natoms, nsolute, nw, swr2max);
	VINT += VINT_EVB;	
	if (nBCD != 1) return;
	nb = nBrO*BrOs;
	nsolvent = natoms - nBCD*BCDs - nsolute; // BCDs = # atom sites in BCD
	nw = nsolvent - nb; //BrOs is # BrO atom sites
        
//fprintf(stderr,"in intForce.c natoms = %d, nsolute = %d, nw = %d swr2max = %f \n",natoms, nsolute, nw, swr2max);
        for(l=0;l<7;l++){ /* loop each glucose unit */
	   for(m=0;m<15;m++){ /* loop non hydroxyl atoms */
	      if(m==14) m=16; /* include ether oxygen (index = 16) */
	      b = l*21 + m + nsolvent; /* index of BCD atom to pass to force routine */   
	      for(i=0;i<nw;i+=3) {/* loop over water molecules*/
	         pe = 0.;
//	         c = -2;
//if(l==6&&i==0) fprintf(stderr,"i = %d, b = %d, m = %d\n",i,b,m);
	         h2oljc(&pos[i],&pos[b],&force[i],&force[b],&pe,/*&c,*/i/3,m);
//		 if(fabs(KCAL*pe)>10.0) fprintf(stderr,"liqForce.c: large water-BCD (l=%d, m=%d,water=%d) pe = %f\n",l,m,i,pe*KCAL);
	         VINTH2O += pe;
	      }
	
	      for(j=0;j<nBrO;j++){/* loop over BRB molecules*/
		 i = nw+ j*BrOs;
	         pe = 0.;
	         brolj(&pos[i],&pos[b],&force[i],&force[b],&pe,m);
	         VINT_BrO += pe;
	         pe = 0.;
		 broq(&pos[i],&pos[b],&force[i],&force[b],&pe,m);
	         VINT_BrO += pe;
	      }
	   }
	   for(m=14;m<20;m+=2){ /* loop over BCD hydroxyl groups */
	      if(m==16) m=17;
	      b = l*21 + m + nsolvent; /* index of BCD atom to pass */
	      for(i=0;i<nw;i+=3){ /* loop over waters */
		 pe = 0.;
//		 c = -2;
//if(l==6&&i==0) fprintf(stderr,"i = %d, b = %d, m = %d\n",i,b,m);
		 h2oohljc(&pos[i],&pos[b],&force[i],&force[b],&pe,/*&c,*/i/3,m);
		 VINTH2Ooh += pe;
	      }
	      for(j=0;j<nBrO;j++) {/* loop over BrO molecules*/
		 i = nw+ j*BrOs;
	         pe = 0.;
	         broohlj(&pos[i],&pos[b],&force[i],&force[b],&pe,m);
	         VINT_BrO += pe;
	         pe = 0.;
	         broohq(&pos[i],&pos[b],&force[i],&force[b],&pe,m);
	         VINT_BrO += pe;
	      }
	   }
	}
//fprintf(stderr, "exiting liqForce.c: VINTH2O = %f, VINTH2Ooh = %f, VINT_BrO = %f \n", VINTH2O*KCAL,VINTH2Ooh*KCAL,VINT_BrO*KCAL);
	if(wGH!=999.9/* && tc%50>3 && tc%50<48*/){
	   hgTeth();
	}
	else if(nBCD==1 && nBrO > 50 && nw > 50)
	   BCDang();
	VINT_BCD += VINT_BrO + VINTH2O + VINTH2Ooh;
	VINT += VINT_BCD;

}

h2oljc(ri,rj,fi,fj,pe,/*c,*/ind,m)
tripd *ri;/* ri[0], ri[1] and ri[2] are position of the O H1 and H2 of water*/
tripd *rj;/* rj[0] is the position of the solute atom*/
tripd *fi;/* fi[0], fi[1] and fi[2] are forces on the O H1 and H2 of water*/
tripd *fj;/* force on solute atom*/
double *pe;/* interaction energy between one water and the solute*/
//double *c;/* water dipole orienataion */
int ind; /* index of water O (i/3) */
int m; /* 0-20 index of BCD atom */
{
double sqrt(), s, sp, f, a, aa, b, bb, q2, r2, r, ir, ir6, EC, rOI,shelMax;
tripd image, rij, sdist;
int i, j;//, rindex, w_index;

/***	Determine Oxygen-BCD site image vector			***/

image.fx = -(sdist.fx = ri[0].fx - rj[0].fx);
image.fy = -(sdist.fy = ri[0].fy - rj[0].fy);
image.fz = -(sdist.fz = ri[0].fz - rj[0].fz);
mvimage(&sdist);
image.fx += sdist.fx;
image.fy += sdist.fy;
image.fz += sdist.fz;
r2 = sdist.fx*sdist.fx +  sdist.fy*sdist.fy +  sdist.fz*sdist.fz;
rOI = sqrt(r2);
//put here to avoid editing water routines in ../water
if ( r2 >= swr2max )
	return;
if ( r2 <= swr2min ){
	s = 1.0;
	sp = 0.0;
 }
else
	swtch(r2 - swr2min, &s, &sp);

/***	Loop over atoms in water molecule			***/

for (i=0; i< 3; i++){
//	if (i == 0) w_index = 2;
//      if (i > 0) w_index = 3;
	rij.fx = ri[i].fx - rj[0].fx + image.fx;
	rij.fy = ri[i].fy - rj[0].fy + image.fy;
	rij.fz = ri[i].fz - rj[0].fz + image.fz;
	r2 = rij.fx*rij.fx +  rij.fy*rij.fy +  rij.fz*rij.fz ;
       	a = wBCDlj[m][i].a;
	aa = 12.*a;
       	b = wBCDlj[m][i].b;
	bb = 6.*b;
       	q2 = wBCDlj[m][i].q;
//if(i==0&&m==0) fprintf(stderr,"b = %f, a = %f, q2 = %f, r2 = %f\n",b,a,q2,r2);
	ir = 1. / r2;
	ir6 = ir * ir * ir;
	r = sqrt(r2);
/*	rindex = r/binRDF;
	if(rindex<300)
		grSOL[w_index-2][rindex] += 1.0/(w_index-1);   
*/
//	if (neqFlag == 2) q2 = -q2;
	EC = q2/r;
	
	*pe +=  ( a * ir6 - b ) * ir6 + EC;
	f = s*((aa * ir6 - bb ) * ir6 * ir + EC/r2); /* -dV/dr/r	*/
	fi[i].fx += (rij.fx *= f);
	fi[i].fy += (rij.fy *= f);
	fi[i].fz += (rij.fz *= f);

	fj[0].fx -= rij.fx;
	fj[0].fy -= rij.fy;
	fj[0].fz -= rij.fz;
}
if (sp != 0.0){
	f =sp*(*pe);
	fi[0].fx -= (sdist.fx *= f);
	fi[0].fy -= (sdist.fy *= f);
	fi[0].fz -= (sdist.fz *= f);

	fj[0].fx += sdist.fx;
	fj[0].fy += sdist.fy;
	fj[0].fz += sdist.fz;
	*pe *= s;
}
//if(fabs(*pe*KCAL)>10.0){fprintf(stderr,"r = %f\n",r);
//fprintf(stderr,"intForce.c: wBCDljq.a = %f, b = %f, q = %f\n",wBCDlj[0][0].a*KCAL, wBCDlj[0][0].b*KCAL,wBCDlj[0][0].q*E2);
//}
}
//////
h2oohljc(ri,rj,fi,fj,pe,/*c,*/ind,m)
tripd *ri;/* ri[0], ri[1] and ri[2] are position of the O H1 and H2 of water*/
tripd *rj;/* rj[0] is the position of the BCD O*/
tripd *fi;/* fi[0], fi[1] and fi[2] are forces on the O H1 and H2 of water*/
tripd *fj;/* force on solute atom*/
double *pe;/* interaction energy between one water and the solute*/
//double *c;/* water dipole orienataion */
int ind; /* index of water O (i/3) */
int m; /* 0-20 index of BCD atom */
{
double sqrt(), s, sp, f, a, aa, b, bb, q2, r2, r, ir, ir6, EC, rOI,shelMax;
tripd image, rij, sdist;
int i, j, rindex, w_index;
double pe2[2];
int ps;

if(m==14) ps=0;
else if(m==17 || m==19) ps=1;

pe2[0] = pe2[1] = 0.0;

/***	Determine water-oxygen - BCD-OH-oxygen-site image vector			***/

image.fx = -(sdist.fx = ri[0].fx - rj[0].fx);
image.fy = -(sdist.fy = ri[0].fy - rj[0].fy);
image.fz = -(sdist.fz = ri[0].fz - rj[0].fz);
mvimage(&sdist);
image.fx += sdist.fx;
image.fy += sdist.fy;
image.fz += sdist.fz;
r2 = sdist.fx*sdist.fx +  sdist.fy*sdist.fy +  sdist.fz*sdist.fz;
rOI = sqrt(r2);
if ( r2 >= swr2max )
	return;
if ( r2 <= swr2min ){
	s = 1.0;
	sp = 0.0;
 }
else
	swtch(r2 - swr2min, &s, &sp);

/***	Loop over atoms in water molecule			***/

for(j=0;j<2;j++){ //loop over BCD O and H
   for (i=0; i< 3; i++){ //loop over water atoms
//	if (i == 0) w_index = 2;
//      if (i > 0) w_index = 3;
	rij.fx = ri[i].fx - rj[j].fx + image.fx;
	rij.fy = ri[i].fy - rj[j].fy + image.fy;
	rij.fz = ri[i].fz - rj[j].fz + image.fz;
	r2 = rij.fx*rij.fx +  rij.fy*rij.fy +  rij.fz*rij.fz ;
       	a = wBCDlj[m+j][i].a;
	aa = 12.*a;
       	b = wBCDlj[m+j][i].b;
	bb = 6.*b;
       	q2 = wBCDlj[m+j][i].q;
	ir = 1. / r2;
	ir6 = ir * ir * ir;
	r = sqrt(r2);
/*	rindex = r/binRDF;
	if (rindex < 300){
	   if(j==0){
		if (i == 0) 
		   grBCDOH[ps][0][rindex] += 1.0;   
		else
		   grBCDOH[ps][1][rindex] += 0.5;   
	   }
	}
*/
//	rindex = r/binRDF;
//	if(rindex<300)
//		grSOL[w_index-2][rindex] += 1.0/(w_index-1);   
//	if (neqFlag == 2) q2 = -q2;
	EC = q2/r;
//	*pe +=  ( a * ir6 - b ) * ir6 + EC;
	pe2[j] +=  ( a * ir6 - b ) * ir6 + EC;
	f = s*((aa * ir6 - bb ) * ir6 * ir + EC/r2); /* -dV/dr/r	*/
	fi[i].fx += (rij.fx *= f);
	fi[i].fy += (rij.fy *= f);
	fi[i].fz += (rij.fz *= f);

	fj[j].fx -= rij.fx;
	fj[j].fy -= rij.fy;
	fj[j].fz -= rij.fz;
   }
   if(sp != 0.0){
	f =sp*(pe2[j]);
	fi[0].fx -= sdist.fx * f;
	fi[0].fy -= sdist.fy * f;
	fi[0].fz -= sdist.fz * f;

	fj[0].fx += sdist.fx * f;
	fj[0].fy += sdist.fy * f;
	fj[0].fz += sdist.fz * f;
/*	fi[0].fx -= (sdist.fx *= f);
	fi[0].fy -= (sdist.fy *= f);
	fi[0].fz -= (sdist.fz *= f);

	fj[j].fx += sdist.fx;
	fj[j].fy += sdist.fy;
	fj[j].fz += sdist.fz;
*/
	pe2[j] *= s;
   }
}

*pe += pe2[0] + pe2[1];

}
//////
brolj(ri,rj,fi,fj,pe,m)
tripd *ri;/* ri[0], ri[1], etc are positions of BrO atoms */
tripd *rj;/* rj[0] is the position of the BCD atom*/
tripd *fi;/* same indices as above for the forces   */
tripd *fj;/* force on solute atom*/
double *pe;/* interaction energy between one NIT and the solute*/
int m;
{
double sqrt(), s, sp, f, a, aa, b, bb, r2, ir, ir6,eg;
tripd rij;
int k;
/*****	Loop over atoms in BrO molecules ******/
for(k = 0; k < BrOs; k++){
/*****	Determine image vector ******/
   rij.fx = ri[k].fx - rj[0].fx;
   rij.fy = ri[k].fy - rj[0].fy;
   rij.fz = ri[k].fz - rj[0].fz;
   mvimage(&rij);
   r2 = rij.fx*rij.fx + rij.fy*rij.fy + rij.fz*rij.fz;
   if(r2 >= swr2max)
	continue;
   if(r2 <= swr2min)
   {
	s = 1.0;
	sp = 0.0;
   }
   else
	swtch(r2-swr2min,&s,&sp);

/*****	Get (1/r dV/dr) ******/
//   dedr = ttraljq(r2,m,n,&eg);
   a = bBCDlj[m][k].a;
   aa = 12.*a;
   b = bBCDlj[m][k].b;
   bb = 6.*b;
//   q2 = bBCDlj[m][3+i].q;
   ir = 1. / r2;
   ir6 = ir * ir * ir;
//   r = sqrt(r2);
//   EC = q2/r;
//   VC_BrO += EC*s;
   eg =  (a * ir6 - b ) * ir6/* + *EC*/;
   f = s*(aa * ir6 - bb) * ir6 * ir/* + s*EC/r2*/; /* -dV/dr/r	*/
   fi[k].fx += (rij.fx *= (f + eg*sp));
   fi[k].fy += (rij.fy *= (f + eg*sp));
   fi[k].fz += (rij.fz *= (f + eg*sp));
   fj[0].fx -= rij.fx;
   fj[0].fy -= rij.fy;
   fj[0].fz -= rij.fz;
   *pe += eg*s;
}
}
//////
broq(ri,rj,fi,fj,pe,m)
tripd *ri;/* ri[0], ri[1] are positions of Br-oct Br,C head */
tripd *rj;/* rj[0] is the position of the BCD atom*/
tripd *fi;/* same indices as above for the forces   */
tripd *fj;/* force on solute atom*/
double *pe;/* q interaction energy between BCD atom and the headgroup*/
int m;
{

int i, r2, q2, r, EC, f, s, sp;
tripd image, rij, sdist;

/***	Determine alpha Carbon-BCD atom image vector	***/

image.fx = -(sdist.fx = ri[0].fx - rj[0].fx);
image.fy = -(sdist.fy = ri[0].fy - rj[0].fy);
image.fz = -(sdist.fz = ri[0].fz - rj[0].fz);
mvimage(&sdist);
image.fx += sdist.fx;
image.fy += sdist.fy;
image.fz += sdist.fz;
r2 = sdist.fx*sdist.fx +  sdist.fy*sdist.fy +  sdist.fz*sdist.fz;
if ( r2 >= swr2max )
	return;
if ( r2 <= swr2min ){
	s = 1.0;
	sp = 0.0;
 }
else
	swtch(r2 - swr2min, &s, &sp);

/***	Loop over atoms in Br-octane headgroup 		***/

for (i=0; i< 2; i++){
	rij.fx = ri[i].fx - rj[0].fx + image.fx;
	rij.fy = ri[i].fy - rj[0].fy + image.fy;
	rij.fz = ri[i].fz - rj[0].fz + image.fz;
	r2 = rij.fx*rij.fx +  rij.fy*rij.fy +  rij.fz*rij.fz ;
//	a = bBCDlj[m][3+i].a;
//	aa = 12.*a;
//	b = bBCDlj[m][3+i].b;
//	bb = 6.*b;
	q2 = bBCDlj[m][i].q;
//	ir = 1. / r2;
//	ir6 = ir * ir * ir;
	r = sqrt(r2);
/*	rindex = r/binRDF;
	if (rindex < 300){
		if (i == 0) grSOL[2][rindex] += 1.0;   
		if (i == 1 || i == 5) grSOL[3][rindex] += 0.5;   
		if (i == 2 || i == 4) grSOL[4][rindex] += 0.5;   
		if (i == 3 ) grSOL[5][rindex] += 1.0;   
		if (i == 6 ) grSOL[6][rindex] += 1.0;   
		if (i == 7 || i == 8 ) grSOL[7][rindex] += 0.5;   
	}
*/
//	if (neqFlag == 2) q2 = -q2;
	EC = q2/r;
//	VC_BrO += EC*s;
//	if (neqFlag == 1)
		/* run with no charge on the ion, but note that VCNIT include *
		 * the coulomb interaction that would be if the charges were  *
		 * on but with the solvent dynamics without the charges       */
//		 EC = 0.;
	
	*pe += /* ( a * ir6 - b ) * ir6 +*/ EC;
	f = s*(/*(aa * ir6 - bb ) * ir6 * ir +*/ EC/r2); /* -dV/dr/r	*/
	fi[i].fx += (rij.fx *= f);
	fi[i].fy += (rij.fy *= f);
	fi[i].fz += (rij.fz *= f);

	fj[0].fx -= rij.fx;
	fj[0].fy -= rij.fy;
	fj[0].fz -= rij.fz;
}
if (sp != 0.0){
	f =sp*(*pe);
	fi[0].fx -= (sdist.fx *= f);
	fi[0].fy -= (sdist.fy *= f);
	fi[0].fz -= (sdist.fz *= f);

	fj[0].fx += sdist.fx;
	fj[0].fy += sdist.fy;
	fj[0].fz += sdist.fz;
	*pe *= s;
}

}

broohlj(ri,rj,fi,fj,pe,m)
tripd *ri;/* ri[0], ri[1] are positions of BrO atoms */
tripd *rj;/* rj[0] is the position of the BCD OH atom*/
tripd *fi;/* same indeces as above for the forces		   */
tripd *fj;/* force on solute atom*/
double *pe;/* interaction energy between one NIT and the solute*/
int m;
{
   int i,j;
   tripd image, sdist, rij;
   double a,aa,b,bb,f,ir,ir6,r2,s,sp,eg;

   for(i=0;i<BrOs;i++){
   
   /***	 Determine BrO atom-BCD OH oxygen atom image vector	***/

      image.fx = -(sdist.fx = ri[i].fx - rj[0].fx);
      image.fy = -(sdist.fy = ri[i].fy - rj[0].fy);
      image.fz = -(sdist.fz = ri[i].fz - rj[0].fz);
      mvimage(&sdist);
      image.fx += sdist.fx;
      image.fy += sdist.fy;
      image.fz += sdist.fz;
      r2 = sdist.fx*sdist.fx +  sdist.fy*sdist.fy +  sdist.fz*sdist.fz;
      if ( r2 >= swr2max )
	continue;
      if ( r2 <= swr2min ){
	s = 1.0;
	sp = 0.0;
       }
      else
	swtch(r2 - swr2min, &s, &sp);
      eg = 0.0;
      for(j=0;j<2;j++){
	rij.fx = ri[i].fx - rj[j].fx + image.fx;
	rij.fy = ri[i].fy - rj[j].fy + image.fy;
	rij.fz = ri[i].fz - rj[j].fz + image.fz;
	r2 = rij.fx*rij.fx +  rij.fy*rij.fy +  rij.fz*rij.fz ;
       	a = bBCDlj[m+j][i].a;
	aa = 12.*a;
       	b = bBCDlj[m+j][i].b;
	bb = 6.*b;
	ir = 1. / r2;
	ir6 = ir * ir * ir;
        eg +=  ( a * ir6 - b ) * ir6/* + *EC*/;
        f = s*(aa * ir6 - bb) * ir6 * ir/* + s*EC/r2*/; /* -dV/dr/r	*/
	fi[i].fx += (rij.fx *= f);
	fi[i].fy += (rij.fy *= f);
	fi[i].fz += (rij.fz *= f);

	fj[0].fx -= rij.fx;
	fj[0].fy -= rij.fy;
	fj[0].fz -= rij.fz;
      }
      if (sp != 0.0){
	f =sp*(eg);
	fi[i].fx -= (sdist.fx *= f);
	fi[i].fy -= (sdist.fy *= f);
	fi[i].fz -= (sdist.fz *= f);

	fj[0].fx += sdist.fx;
	fj[0].fy += sdist.fy;
	fj[0].fz += sdist.fz;
	eg *= s;
      }
      *pe += eg;
   }
}
broohq(ri,rj,fi,fj,pe,m)
tripd *ri;/* ri[0], ri[1] are positions of the alpha C, Br   */
tripd *rj;/* rj[0], rj[1] are positions of the BCD OH O,H*/
tripd *fi;/* same indeces as above for the forces		   */
tripd *fj;/* force on BCD atom*/
double *pe;/* interaction energy between Br-C headgroup and BCD OH group*/
int m;
{
   int i,j,d;
   tripd image,sdist,rij,f[4],frc;
   double dedr,r2,r,q2,ec,s,sp,eg;

/***	Determine alpha Carbon-BCD OH O atom image vector	***/

image.fx = -(sdist.fx = ri[0].fx - rj[0].fx);
image.fy = -(sdist.fy = ri[0].fy - rj[0].fy);
image.fz = -(sdist.fz = ri[0].fz - rj[0].fz);
mvimage(&sdist);
image.fx += sdist.fx;
image.fy += sdist.fy;
image.fz += sdist.fz;
r2 = sdist.fx*sdist.fx +  sdist.fy*sdist.fy +  sdist.fz*sdist.fz;
if ( r2 >= swr2max )
   return;
if ( r2 <= swr2min ){
	s = 1.0;
	sp = 0.0;
}
else
	swtch(r2 - swr2min, &s, &sp);

/***	Loop over atoms in Br-octane headgroup 		***/
for(d=0;d<4;d++){
   f[d].fx = f[d].fy = f[d].fz = 0.0;
}
eg = 0.0;
for(i=0; i< 2; i++){
   for (j=0; j< 2; j++){
	rij.fx = ri[i].fx - rj[j].fx + image.fx;
	rij.fy = ri[i].fy - rj[j].fy + image.fy;
	rij.fz = ri[i].fz - rj[j].fz + image.fz;
	r2 = rij.fx*rij.fx +  rij.fy*rij.fy +  rij.fz*rij.fz ;
//	a = bBCDlj[m][3+i].a;
//	aa = 12.*a;
//	b = bBCDlj[m][3+i].b;
//	bb = 6.*b;
	q2 = bBCDlj[m+j][i].q;
//	ir = 1. / r2;
//	ir6 = ir * ir * ir;
	r = sqrt(r2);
	ec = q2/r;
	eg += ec;
	dedr = -ec/r2;
	f[j+2].fx += (rij.fx *= dedr);
	f[j+2].fy += (rij.fy *= dedr);
	f[j+2].fz += (rij.fz *= dedr);
	f[i].fx -= rij.fx;
	f[i].fy -= rij.fy;
	f[i].fz -= rij.fz;
   }
}

if(sp != 0.0)
{
   for(d=0;d<4;d++)
   {
      f[d].fx *= s;
      f[d].fy *= s;
      f[d].fz *= s;
   }
   f[0].fx -= (frc.fx = sp*eg*sdist.fx);
   f[0].fy -= (frc.fy = sp*eg*sdist.fy);
   f[0].fz -= (frc.fz = sp*eg*sdist.fz);
   f[2].fx+= frc.fx;
   f[2].fy+= frc.fy;
   f[2].fz+= frc.fz;
   eg *= s;
}

fi[0].fx += f[0].fx;
fi[0].fy += f[0].fy;
fi[0].fz += f[0].fz;
fi[1].fx += f[1].fx;
fi[1].fy += f[1].fy;
fi[1].fz += f[1].fz;

fj[0].fx += f[2].fx;
fj[0].fy += f[2].fy;
fj[0].fz += f[2].fz;
fj[1].fx += f[3].fx;
fj[1].fy += f[3].fy;
fj[1].fz += f[3].fz;

*pe += eg;
}
	  
#define HGTETH 1.0
#define HGANG 50.0

hgTeth(){

int i,j,nsolvent,nw;
double MBO,MCD,deriv;
double k,zeta,gamma;
double U,cosGH,ccr;
double Wc,W_width,rBrO;
tripd sdl,image,c,b,p,cb,cb_hat,a8;

c.fx = BCDcom.fx;
c.fy = BCDcom.fy;
c.fz = BCDcom.fz;
p.fx = BCDz.fx;
p.fy = BCDz.fy;
p.fz = BCDz.fz;

nsolvent = natoms - nsolute - BCDs*nBCD;
// nw = number of water atoms
nw = nsolvent - nBrO*BrOs;
U = 0.0;
MBO = MCD = 0.0;
Wc = zGH;
W_width = wGH;

for(i=0;i<BCDs;i++){
   j = nsolvent + i;
   MCD += mass[j];
}
// b is guest molecule center of mass
b.fx = b.fy = b.fz = 0.0;
for(i=0;i<BrOs;i++){
   b.fx += pos[nw+i].fx * mass[nw+i];
   b.fy += pos[nw+i].fy * mass[nw+i];
   b.fz += pos[nw+i].fz * mass[nw+i];
   MBO += mass[nw+i];
}
b.fx /= MBO;
b.fy /= MBO;
b.fz /= MBO;

a8.fx = pos[nw].fx - pos[nw+8].fx;
a8.fy = pos[nw].fy - pos[nw+8].fy;
a8.fz = pos[nw].fz - pos[nw+8].fz;
rBrO = sqrt(sq(a8.fx)+sq(a8.fy)+sq(a8.fz));
a8.fx/=rBrO;
a8.fy/=rBrO;
a8.fz/=rBrO;
g_ang = p.fx*a8.fx+p.fy*a8.fy+p.fz*a8.fz;

image.fx = -(sdl.fx = pos[nw].fx - c.fx);
image.fy = -(sdl.fy = pos[nw].fy - c.fy);
image.fz = -(sdl.fz = pos[nw].fz - c.fz);
mvimage(&sdl);
image.fx += sdl.fx;
image.fy += sdl.fy;
image.fz += sdl.fz;

cb.fx = b.fx - c.fx + image.fx;
cb.fy = b.fy - c.fy + image.fy;
cb.fz = b.fz - c.fz + image.fz;

gC8 = p.fx*(pos[nw+8].fx - c.fx + image.fx) + p.fy*(pos[nw+8].fy - c.fy + image.fy) + p.fz*(pos[nw+8].fz - c.fz + image.fz);
gaC = p.fx*sdl.fx + p.fy*sdl.fy + p.fz*sdl.fz;
// dz is position of gcom on BCDz
gamma = p.fx*cb.fx + p.fy*cb.fy + p.fz*cb.fz;
// apply biasing potential
if(guestBiasOn){
  gBias = -(gbA*(gamma-gbO)*(gamma-gbO)-gbB)/KCAL;
  deriv = (2.0*gbA*(gamma-gbO))/KCAL;
  for(i=0;i<BrOs;i++){
     force[nw+i].fx += deriv*p.fx*mass[nw+i]/MBO;
     force[nw+i].fy += deriv*p.fy*mass[nw+i]/MBO;
     force[nw+i].fz += deriv*p.fz*mass[nw+i]/MBO;
  }
  for(i=0;i<BCDs;i++){
     j = nsolvent + i;
     if(i%21==3 || i%21==5){
         force[j].fx -= deriv*(p.fx*mass[j]/MCD + cb.fx*mass[j]/(capR*Ml) - p.fx*gamma*mass[j]/(capR*Ml));
	 force[j].fy -= deriv*(p.fy*mass[j]/MCD + cb.fy*mass[j]/(capR*Ml) - p.fy*gamma*mass[j]/(capR*Ml));
	 force[j].fz -= deriv*(p.fz*mass[j]/MCD + cb.fz*mass[j]/(capR*Ml) - p.fz*gamma*mass[j]/(capR*Ml));
     }
     else if(i%21==9 || i%21==16){
	 force[j].fx -= deriv*(p.fx*mass[j]/MCD - cb.fx*mass[j]/(capR*Ms) + p.fx*gamma*mass[j]/(capR*Ms));
	 force[j].fy -= deriv*(p.fy*mass[j]/MCD - cb.fy*mass[j]/(capR*Ms) + p.fy*gamma*mass[j]/(capR*Ms));
	 force[j].fz -= deriv*(p.fz*mass[j]/MCD - cb.fz*mass[j]/(capR*Ms) + p.fz*gamma*mass[j]/(capR*Ms));
     }
     else{
	 force[j].fx -= deriv*(p.fx*mass[j]/MCD);
	 force[j].fy -= deriv*(p.fy*mass[j]/MCD);
	 force[j].fz -= deriv*(p.fz*mass[j]/MCD);
     }
   }
   VINT += gBias;
}

//ccr is comcom radius
ccr = sqrt(sq(cb.fx)+sq(cb.fy)+sq(cb.fz));
//gamma = ccr;
// comcom is BCDcom-BrOcom direction vector
cb_hat.fx = cb.fx/ccr;
cb_hat.fy = cb.fy/ccr;
cb_hat.fz = cb.fz/ccr;
//cosine BCDcom - guest CoM
cosGH = cb_hat.fx*p.fx + cb_hat.fy*p.fy + cb_hat.fz*p.fz;

zeta = fabs(gamma - Wc) - (W_width/2.0); 
//zeta = gamma - Wc - (W_width/2.0); 

if(zeta > 0.0){
   if(zeta > 2.0*W_width){
      fprintf(stderr,"ERROR: BCDz window-- guest too far outside window. zeta = %f\n",zeta);
      fprintf(stderr,"ERROR: ... G-H radius = %f, gamma = %f\n",ccr,gamma);
      fprintf(stderr,"ERROR: ... c.fx = %f, c.fy = %f, c.fy = %f\n",c.fx,c.fy,c.fz);
      fprintf(stderr,"ERROR: ... aC.fx = %f, aC.fy = %f, aC.fy = %f\n",pos[nw].fx,pos[nw].fy,pos[nw].fz);
      exit(1);
   }
   else{
      U = HGTETH*pow(zeta,3.);
      deriv = -3.0*pow(zeta,2.)*HGTETH*sgn(gamma-Wc);

      for(i=0;i<BrOs;i++){
	 force[nw+i].fx += deriv*p.fx*mass[nw+i]/MBO;
	 force[nw+i].fy += deriv*p.fy*mass[nw+i]/MBO;
	 force[nw+i].fz += deriv*p.fz*mass[nw+i]/MBO;
      }

      for(i=0;i<BCDs;i++){
	j = nsolvent + i;
	if(i%21==3 || i%21==5){
	 force[j].fx -= deriv*(p.fx*mass[j]/MCD + cb.fx*mass[j]/(capR*Ml) - p.fx*gamma*mass[j]/(capR*Ml));
	 force[j].fy -= deriv*(p.fy*mass[j]/MCD + cb.fy*mass[j]/(capR*Ml) - p.fy*gamma*mass[j]/(capR*Ml));
	 force[j].fz -= deriv*(p.fz*mass[j]/MCD + cb.fz*mass[j]/(capR*Ml) - p.fz*gamma*mass[j]/(capR*Ml));
	}
	else if(i%21==9 || i%21==16){
	 force[j].fx -= deriv*(p.fx*mass[j]/MCD - cb.fx*mass[j]/(capR*Ms) + p.fx*gamma*mass[j]/(capR*Ms));
	 force[j].fy -= deriv*(p.fy*mass[j]/MCD - cb.fy*mass[j]/(capR*Ms) + p.fy*gamma*mass[j]/(capR*Ms));
	 force[j].fz -= deriv*(p.fz*mass[j]/MCD - cb.fz*mass[j]/(capR*Ms) + p.fz*gamma*mass[j]/(capR*Ms));
	}
	else{
	 force[j].fx -= deriv*(p.fx*mass[j]/MCD);
	 force[j].fy -= deriv*(p.fy*mass[j]/MCD);
	 force[j].fz -= deriv*(p.fz*mass[j]/MCD);
	}
      }
   }
}
VINT += U;
U_teth = 0.0;
U_teth += U*KCAL;
guestZ = gamma;

// Now impose angular window if guest-windowed l/l system

V_gAng = 0.0;
U = 0.0;
//if(1){
if(nw > 0 && nBrO > 1){
   gamma = p.fz;
   zeta = 0.95 - fabs(gamma);
   if(zeta > 0.0){
      U = HGANG*pow(zeta,3.);
      deriv = 3.*HGANG*pow(zeta,2.)*sgn(gamma);
      
      for(i=0;i<BCDs;i++){
	j = nsolvent + i;
	if(i%21==3 || i%21==5){ // large ring
	 force[j].fx += deriv*gamma*(mass[j]/(capR*Ml))*p.fx;
	 force[j].fy += deriv*gamma*(mass[j]/(capR*Ml))*p.fy;
	 force[j].fz += deriv*(mass[j]/(capR*Ml))*(gamma*gamma-1.);
	}
	else if(i%21==9 || i%21==16){ // small ring
	 force[j].fx -= deriv*gamma*(mass[j]/(capR*Ms))*p.fx;
	 force[j].fy -= deriv*gamma*(mass[j]/(capR*Ms))*p.fy;
	 force[j].fz -= deriv*(mass[j]/(capR*Ms))*(gamma*gamma-1.);
	}
      }
      V_gAng = U*KCAL;
      VINT += U;
   }
}


}

BCDang(){

int i,j,nsolvent,nw;
double MBO,MCD,deriv;
double k,zeta,gamma;
double U,cosGH,ccr;
double Wc,W_width,rBrO;
tripd sdl,image,c,b,p,cb,cb_hat,a8;

c.fx = BCDcom.fx;
c.fy = BCDcom.fy;
c.fz = BCDcom.fz;
p.fx = BCDz.fx;
p.fy = BCDz.fy;
p.fz = BCDz.fz;

nsolvent = natoms - nsolute - BCDs*nBCD;
// nw = number of water atoms
nw = nsolvent - nBrO*BrOs;
U = 0.0;
MBO = MCD = 0.0;
Wc = zGH;
W_width = wGH;

for(i=0;i<BCDs;i++){
   j = nsolvent + i;
   MCD += mass[j];
}
// b is guest molecule center of mass
/*
b.fx = b.fy = b.fz = 0.0;
for(i=0;i<BrOs;i++){
   b.fx += pos[nw+i].fx * mass[nw+i];
   b.fy += pos[nw+i].fy * mass[nw+i];
   b.fz += pos[nw+i].fz * mass[nw+i];
   MBO += mass[nw+i];
}
b.fx /= MBO;
b.fy /= MBO;
b.fz /= MBO;

a8.fx = pos[nw].fx - pos[nw+8].fx;
a8.fy = pos[nw].fy - pos[nw+8].fy;
a8.fz = pos[nw].fz - pos[nw+8].fz;
rBrO = sqrt(sq(a8.fx)+sq(a8.fy)+sq(a8.fz));
a8.fx/=rBrO;
a8.fy/=rBrO;
a8.fz/=rBrO;
g_ang = p.fx*a8.fx+p.fy*a8.fy+p.fz*a8.fz;

image.fx = -(sdl.fx = pos[nw].fx - c.fx);
image.fy = -(sdl.fy = pos[nw].fy - c.fy);
image.fz = -(sdl.fz = pos[nw].fz - c.fz);
mvimage(&sdl);
image.fx += sdl.fx;
image.fy += sdl.fy;
image.fz += sdl.fz;

cb.fx = b.fx - c.fx + image.fx;
cb.fy = b.fy - c.fy + image.fy;
cb.fz = b.fz - c.fz + image.fz;

gC8 = p.fx*(pos[nw+8].fx - c.fx + image.fx) + p.fy*(pos[nw+8].fy - c.fy + image.fy) + p.fz*(pos[nw+8].fz - c.fz + image.fz);
gaC = p.fx*sdl.fx + p.fy*sdl.fy + p.fz*sdl.fz;
// dz is position of gcom on BCDz
gamma = p.fx*cb.fx + p.fy*cb.fy + p.fz*cb.fz;
// apply biasing potential
if(guestBiasOn){
  gBias = -(gbA*(gamma-gbO)*(gamma-gbO)-gbB)/KCAL;
  deriv = (2.0*gbA*(gamma-gbO))/KCAL;
  for(i=0;i<BrOs;i++){
     force[nw+i].fx += deriv*p.fx*mass[nw+i]/MBO;
     force[nw+i].fy += deriv*p.fy*mass[nw+i]/MBO;
     force[nw+i].fz += deriv*p.fz*mass[nw+i]/MBO;
  }
  for(i=0;i<BCDs;i++){
     j = nsolvent + i;
     if(i%21==3 || i%21==5){
         force[j].fx -= deriv*(p.fx*mass[j]/MCD + cb.fx*mass[j]/(capR*Ml) - p.fx*gamma*mass[j]/(capR*Ml));
	 force[j].fy -= deriv*(p.fy*mass[j]/MCD + cb.fy*mass[j]/(capR*Ml) - p.fy*gamma*mass[j]/(capR*Ml));
	 force[j].fz -= deriv*(p.fz*mass[j]/MCD + cb.fz*mass[j]/(capR*Ml) - p.fz*gamma*mass[j]/(capR*Ml));
     }
     else if(i%21==9 || i%21==16){
	 force[j].fx -= deriv*(p.fx*mass[j]/MCD - cb.fx*mass[j]/(capR*Ms) + p.fx*gamma*mass[j]/(capR*Ms));
	 force[j].fy -= deriv*(p.fy*mass[j]/MCD - cb.fy*mass[j]/(capR*Ms) + p.fy*gamma*mass[j]/(capR*Ms));
	 force[j].fz -= deriv*(p.fz*mass[j]/MCD - cb.fz*mass[j]/(capR*Ms) + p.fz*gamma*mass[j]/(capR*Ms));
     }
     else{
	 force[j].fx -= deriv*(p.fx*mass[j]/MCD);
	 force[j].fy -= deriv*(p.fy*mass[j]/MCD);
	 force[j].fz -= deriv*(p.fz*mass[j]/MCD);
     }
   }
   VINT += gBias;
}

//ccr is comcom radius
ccr = sqrt(sq(cb.fx)+sq(cb.fy)+sq(cb.fz));
//gamma = ccr;
// comcom is BCDcom-BrOcom direction vector
cb_hat.fx = cb.fx/ccr;
cb_hat.fy = cb.fy/ccr;
cb_hat.fz = cb.fz/ccr;
//cosine BCDcom - guest CoM
cosGH = cb_hat.fx*p.fx + cb_hat.fy*p.fy + cb_hat.fz*p.fz;

zeta = fabs(gamma - Wc) - (W_width/2.0); 
//zeta = gamma - Wc - (W_width/2.0); 

if(zeta > 0.0){
   if(zeta > 2.0*W_width){
      fprintf(stderr,"ERROR: BCDz window-- guest too far outside window. zeta = %f\n",zeta);
      fprintf(stderr,"ERROR: ... G-H radius = %f, gamma = %f\n",ccr,gamma);
      fprintf(stderr,"ERROR: ... c.fx = %f, c.fy = %f, c.fy = %f\n",c.fx,c.fy,c.fz);
      fprintf(stderr,"ERROR: ... aC.fx = %f, aC.fy = %f, aC.fy = %f\n",pos[nw].fx,pos[nw].fy,pos[nw].fz);
      exit(1);
   }
   else{
      U = HGTETH*pow(zeta,3.);
      deriv = -3.0*pow(zeta,2.)*HGTETH*sgn(gamma-Wc);

      for(i=0;i<BrOs;i++){
	 force[nw+i].fx += deriv*p.fx*mass[nw+i]/MBO;
	 force[nw+i].fy += deriv*p.fy*mass[nw+i]/MBO;
	 force[nw+i].fz += deriv*p.fz*mass[nw+i]/MBO;
      }

      for(i=0;i<BCDs;i++){
	j = nsolvent + i;
	if(i%21==3 || i%21==5){
	 force[j].fx -= deriv*(p.fx*mass[j]/MCD + cb.fx*mass[j]/(capR*Ml) - p.fx*gamma*mass[j]/(capR*Ml));
	 force[j].fy -= deriv*(p.fy*mass[j]/MCD + cb.fy*mass[j]/(capR*Ml) - p.fy*gamma*mass[j]/(capR*Ml));
	 force[j].fz -= deriv*(p.fz*mass[j]/MCD + cb.fz*mass[j]/(capR*Ml) - p.fz*gamma*mass[j]/(capR*Ml));
	}
	else if(i%21==9 || i%21==16){
	 force[j].fx -= deriv*(p.fx*mass[j]/MCD - cb.fx*mass[j]/(capR*Ms) + p.fx*gamma*mass[j]/(capR*Ms));
	 force[j].fy -= deriv*(p.fy*mass[j]/MCD - cb.fy*mass[j]/(capR*Ms) + p.fy*gamma*mass[j]/(capR*Ms));
	 force[j].fz -= deriv*(p.fz*mass[j]/MCD - cb.fz*mass[j]/(capR*Ms) + p.fz*gamma*mass[j]/(capR*Ms));
	}
	else{
	 force[j].fx -= deriv*(p.fx*mass[j]/MCD);
	 force[j].fy -= deriv*(p.fy*mass[j]/MCD);
	 force[j].fz -= deriv*(p.fz*mass[j]/MCD);
	}
      }
   }
}
VINT += U;
U_teth = 0.0;
U_teth += U*KCAL;
guestZ = gamma;
*/

// Now impose angular window if guest-windowed l/l system

V_gAng = 0.0;
U = 0.0;
//if(1){
if(nw > 0 && nBrO > 1){
   gamma = p.fz;
   zeta = 0.95 - fabs(gamma);
   if(zeta > 0.0){
      U = HGANG*pow(zeta,3.);
      deriv = 3.*HGANG*pow(zeta,2.)*sgn(gamma);
      
      for(i=0;i<BCDs;i++){
	j = nsolvent + i;
	if(i%21==3 || i%21==5){ // large ring
	 force[j].fx += deriv*gamma*(mass[j]/(capR*Ml))*p.fx;
	 force[j].fy += deriv*gamma*(mass[j]/(capR*Ml))*p.fy;
	 force[j].fz += deriv*(mass[j]/(capR*Ml))*(gamma*gamma-1.);
	}
	else if(i%21==9 || i%21==16){ // small ring
	 force[j].fx -= deriv*gamma*(mass[j]/(capR*Ms))*p.fx;
	 force[j].fy -= deriv*gamma*(mass[j]/(capR*Ms))*p.fy;
	 force[j].fz -= deriv*(mass[j]/(capR*Ms))*(gamma*gamma-1.);
	}
      }
      V_gAng = U*KCAL;
      VINT += U;
   }
}


}

