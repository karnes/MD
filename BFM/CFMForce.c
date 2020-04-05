#include	<md.h>
#include	<system.h>
#include	<math.h>
#define cosHCHmax (-0.866025)
#define cosHCCmax (-0.866025)
CFMForce()
{
int	i, j, k, l, m, n, ir, io;
double headTail(),comcom;
double	r2, r, dedr, delfx, delfy, delfz, dotp;
double		eg, s, sp;
tripd		frc, image, sdl, f[10],CFdip[800];
double		ttraljq();
double contX,contY,dx,dy,rCC;
tripd iCH,ijCH,ijCCl,ijCC;
double ilen, ijlen,cosHCH,cosHCCl,cosalp,alpha,cosHCC;
int zz,angBin,angBinsm;
double ang,angMod;
double circR;
int binR,xbin,ybin;
double dotp1, dotp2,radd,dipdip;
tripd jl2,ijdip,icom,jcom,ref;
tripd A1,A2,B1,B2;
tripd A2A1,B2B1;

if(tc%TCFdt==0){
 if(dataRatex > 0){
   for(i = 0; i < nCFM; i++){
      k = i*5;
      getc3vec(k);
   }
   for(i = 0; i < nCFM; i++){
      k = i*5;
      n = 4;
/*****	Get unit vector in the direction of electric dipole  for i'th CFM.  ******/
//      getc3vec(k);
     /*
       *  POLAR STACKING
       *  "To determine the percentages of molecules 
       *  that take part in stacks of CHCl3 molecules with 
       *  approximately collinear dipoles we define a C–H distance 
       *  range from 2 to 4.2 Å and an H–C...H angle range 
       *  from theta = 150 to 210 as the condition for polar stacking"
       *  from Chem. Comm., 2015,51, 4770-4773
       */
      /** loop to look for stacking **/
      for(j=0;j<nCFM;j++){
         if(i!=j){
            n = 4; // hydrogen
            l = j*5;
	    k = i*5;
            /** Determine image vector for Cj -> Ci **/
	    image.fx = -(sdl.fx = pos[k].fx - pos[l].fx);
	    image.fy = -(sdl.fy = pos[k].fy - pos[l].fy);
	    image.fz = -(sdl.fz = pos[k].fz - pos[l].fz);
	    mvimage(&sdl);
	    image.fx += (delfx = sdl.fx);
	    image.fy += (delfy = sdl.fy);
	    image.fz += (delfz = sdl.fz);
	    rCC = sqrt(sdl.fx*sdl.fx + sdl.fy*sdl.fy + sdl.fz*sdl.fz);
	    /* ith C -> jth C unit vector */
	    ijCC.fx = -sdl.fx/rCC;
	    ijCC.fy = -sdl.fy/rCC;
	    ijCC.fz = -sdl.fz/rCC;
	    /* ith C -> jth H distance  */
	    delfx = -(pos[k].fx - pos[l+n].fx + image.fx);
	    delfy = -(pos[k].fy - pos[l+n].fy + image.fy);
	    delfz = -(pos[k].fz - pos[l+n].fz + image.fz);
	    r2 = delfx*delfx + delfy*delfy + delfz*delfz;
	    /* get ith C->H vector */
	    iCH.fx = pos[k+n].fx - pos[k].fx;
	    iCH.fy = pos[k+n].fy - pos[k].fy;
	    iCH.fz = pos[k+n].fz - pos[k].fz;
	    ilen = sqrt(iCH.fx*iCH.fx + iCH.fy*iCH.fy + iCH.fz*iCH.fz);
	    /* get HC...C cosine */
	    cosHCC = (ijCC.fx*iCH.fx + ijCC.fy*iCH.fy + ijCC.fz*iCH.fz)/ilen;
	    /* get ith C-> jth H vector */
	    ijCH.fx = pos[l+n].fx - pos[k].fx + image.fx;
	    ijCH.fy = pos[l+n].fy - pos[k].fy + image.fy;
	    ijCH.fz = pos[l+n].fz - pos[k].fz + image.fz;
	    ijlen = sqrt(ijCH.fx*ijCH.fx + ijCH.fy*ijCH.fy + ijCH.fz*ijCH.fz);
	    /* get cos */
	    cosHCH = (iCH.fx*ijCH.fx + iCH.fy*ijCH.fy + iCH.fz*ijCH.fz)/(ilen*ijlen);
	    if(r2 < CHmax2 && r2 > CHmin2){
	       if(cosHCH<cosHCHmax){
	          // then i is above j in a 'polar stack'
	          // assign value
		  dipdip = dipVec[i][(int)(tc/TCFdt)].fx*dipVec[j][(int)(tc/TCFdt)].fx + 
		          dipVec[i][(int)(tc/TCFdt)].fy*dipVec[j][(int)(tc/TCFdt)].fy + 
		          dipVec[i][(int)(tc/TCFdt)].fz*dipVec[j][(int)(tc/TCFdt)].fz;
//		  if(-dipdip < cosHCCmax){
		  avgDipDip += dipdip;
		  avgDipDipn+= 1.0;
	          pstack[i][(int)(tc/TCFdt)] = j;
		  
		  stackCos[(int)((dipdip+1.0)/dStkCos)]++;
		  nApollo++;
		  if(dipdip>cosApollo)
		     apollo++;
		  
	       }
	    }
	    if(rCC < rCCmax && rCC > rCCmin){
	       if(cosHCC < cosHCCmax){
		  totInHCC += 1.0;
		  dipdip = dipVec[i][(int)(tc/TCFdt)].fx*dipVec[j][(int)(tc/TCFdt)].fx + 
		          dipVec[i][(int)(tc/TCFdt)].fy*dipVec[j][(int)(tc/TCFdt)].fy + 
		          dipVec[i][(int)(tc/TCFdt)].fz*dipVec[j][(int)(tc/TCFdt)].fz;
		  if(-dipdip < cosHCCmax){
		     pstack2[i][(int)(tc/TCFdt)] = j;
		     avgDipDip2 += dipdip;
		     avgDipDip2n+= 1.0;
		  }
	       }
	    }		     
	    if(r2 < (contRad*contRad)){
		r = sqrt(r2);
		dy = r*cosHCH;
		dx = r*sqrt(1-(cosHCH*cosHCH));
		ybin = (int)(floor((dy+contRad)/contBin));
		xbin = (int)(dx/contBin);
		CHmap[ybin][xbin]++;
	    }
	    for(n=1;n<4;n++){ // loop around chlorine atoms of jth CFM
		 k=i*5;
		 l=j*5;
		 ijCCl.fx = pos[l+n].fx - pos[k].fx + image.fx;	
		 ijCCl.fy = pos[l+n].fy - pos[k].fy + image.fy;	
		 ijCCl.fz = pos[l+n].fz - pos[k].fz + image.fz;
		 r = sqrt(sq(ijCCl.fx)+sq(ijCCl.fy)+sq(ijCCl.fz));
		 if(r<contRad){
			cosHCCl = (iCH.fx*ijCCl.fx + iCH.fy*ijCCl.fy + iCH.fz*ijCCl.fz)/(ilen*r);
			dy = r*cosHCCl;
			dx = r*sqrt(1-(cosHCCl*cosHCCl));
			ybin = (int)(floor((dy+contRad)/contBin));
			xbin = (int)(dx/contBin);
			CClmap[ybin][xbin]++;
		 }
	    }
	    if(rCC < contRad){
		dy = rCC*cosHCC;
		dx = rCC*sqrt(1-(cosHCC*cosHCC));
		ybin = (int)(floor((dy+contRad)/contBin));
		xbin = (int)(dx/contBin);
		CCmap[ybin][xbin]++;
	    }
	 }
      }
   }

   for(i=0;i<nCFM;i++){
      k=i*5;
      /* loop to look for co-planar CFMs for contour plots */
      /* find CoM */
      icom.fx = icom.fy = icom.fz = 0.0;
      for(n=0;n<5;n++){
         icom.fx += pos[k+n].fx*mass[k+n];
         icom.fy += pos[k+n].fy*mass[k+n];
         icom.fz += pos[k+n].fz*mass[k+n];
      }
      icom.fx /= CFMMass;
      icom.fy /= CFMMass;
      icom.fz /= CFMMass;
      /* get i dipole endpoints */
      A1.fx = (pos[k].fx*QCFM[0] + pos[k+4].fx*QCFM[4])/(QCFM[0]+QCFM[4]);
      A1.fy = (pos[k].fy*QCFM[0] + pos[k+4].fy*QCFM[4])/(QCFM[0]+QCFM[4]);
      A1.fz = (pos[k].fz*QCFM[0] + pos[k+4].fz*QCFM[4])/(QCFM[0]+QCFM[4]);
      A2.fx = (pos[k+1].fx + pos[k+2].fx + pos[k+3].fx)/3.0;
      A2.fy = (pos[k+1].fy + pos[k+2].fy + pos[k+3].fy)/3.0;
      A2.fz = (pos[k+1].fz + pos[k+2].fz + pos[k+3].fz)/3.0;
      /* get ith dipole vector */
      A2A1.fx = A1.fx - A2.fx;
      A2A1.fy = A1.fy - A2.fy;
      A2A1.fz = A1.fz - A2.fz;
      r = sqrt(sq(A2A1.fx)+sq(A2A1.fy)+sq(A2A1.fz));
      A2A1.fx/=r;
      A2A1.fy/=r;
      A2A1.fz/=r;
      for(j=0;j<nCFM;j++){
       if(j!=i){
  	 l=j*5;
         /** Determine image vector for Ci -> Cj **/
	 image.fx = -(sdl.fx = pos[k].fx - pos[l].fx);
	 image.fy = -(sdl.fy = pos[k].fy - pos[l].fy);
	 image.fz = -(sdl.fz = pos[k].fz - pos[l].fz);
	 mvimage(&sdl);
	 image.fx += (delfx = sdl.fx);
	 image.fy += (delfy = sdl.fy);
	 image.fz += (delfz = sdl.fz);
	 jcom.fx = jcom.fy = jcom.fz = 0.0;
	 // get CoM of jth CFM
	 for(n=0;n<5;n++){
	    jcom.fx += pos[l+n].fx*mass[l+n];
	    jcom.fy += pos[l+n].fy*mass[l+n];
	    jcom.fz += pos[l+n].fz*mass[l+n];
	 }
	 jcom.fx /= CFMMass;
	 jcom.fy /= CFMMass;
	 jcom.fz /= CFMMass;
	 sdl.fx = icom.fx - jcom.fx + image.fx;
	 sdl.fy = icom.fy - jcom.fy + image.fy;
	 sdl.fz = icom.fz - jcom.fz + image.fz;
	 // r2 is com-com distance
	 r2 = sdl.fx*sdl.fx + sdl.fy*sdl.fy + sdl.fz*sdl.fz;
	 comcom = sqrt(r2);
	 if(comcom < 29.0){
	    CFMRdf[6][(int)(comcom/binRDF)] += 0.5; 
	 }
	 if(comcom < gKrad){
	    gK[(int)(comcom/gKbin)] += dipVec[i][(int)(tc/TCFdt)].fx*dipVec[j][(int)(tc/TCFdt)].fx +
	    			dipVec[i][(int)(tc/TCFdt)].fy*dipVec[j][(int)(tc/TCFdt)].fy +
	    			dipVec[i][(int)(tc/TCFdt)].fz*dipVec[j][(int)(tc/TCFdt)].fz;
	    gKn[(int)(comcom/gKbin)]+=0.5;
	 }
	 if(r2<(contRad*contRad)){ // within radius. check theta.
//	    fprintf(stderr,"In range!\n");
	    /* if ith C->H -- CoM angle is OK */
	    angBin = getAngBin(A2A1,icom,jcom,image); 
	    angBinsm = getAngBinsm(A2A1,icom,jcom,image); 
	    // determines if the i - j angle
	    // corresponds to a plot we wish to construct
	    // {0º,45º,...,180º} +/- dAng 
	    if(angBin != 5){
//	    fprintf(stderr,"In an angle bin!\n");
	       /* check if planar */
	       /* get dipole endpoints */
	       B1.fx = (pos[l].fx*QCFM[0] + pos[l+4].fx*QCFM[4])/(QCFM[0]+QCFM[4]);
	       B1.fy = (pos[l].fy*QCFM[0] + pos[l+4].fy*QCFM[4])/(QCFM[0]+QCFM[4]);
	       B1.fz = (pos[l].fz*QCFM[0] + pos[l+4].fz*QCFM[4])/(QCFM[0]+QCFM[4]);
	       B2.fx = (pos[l+1].fx + pos[l+2].fx + pos[l+3].fx)/3.0;
      	       B2.fy = (pos[l+1].fy + pos[l+2].fy + pos[l+3].fy)/3.0;
      	       B2.fz = (pos[l+1].fz + pos[l+2].fz + pos[l+3].fz)/3.0;
	       B1.fx+=image.fx;
	       B1.fy+=image.fy;
	       B1.fz+=image.fz;
	       B2.fx+=image.fx;
	       B2.fy+=image.fy;
	       B2.fz+=image.fz;
	       /* check if planar */
	       if(angBin==0 || angBin==4 || inPlane(A1,A2,B1,B2)){
	          angMod = headTail(A1,A2,B1,B2); // determ if angle is 0-180(0) or 180-360(360)
		  B2B1.fx = B1.fx - B2.fx;
		  B2B1.fy = B1.fy - B2.fy;
		  B2B1.fz = B1.fz - B2.fz;
		  r = sqrt(sq(B2B1.fx)+sq(B2B1.fy)+sq(B2B1.fz));
		  B2B1.fx/=r;
		  B2B1.fy/=r;
		  B2B1.fz/=r;

		  cosalp = A2A1.fx*B2B1.fx + A2A1.fy*B2B1.fy + A2A1.fz*B2B1.fz;
		  if(comcom < contRad){
		     dy = comcom*cosalp;
		     dx = comcom*sqrt(1.0-sq(cosalp));
		     ybin = (int)(floor((dy+contRad)/contBin) );
		     if(angBin != 0 && angBin != 4){
		        if(angMod==360.0){
			   xbin = (int)(floor((-dx+contRad)/contBin));
		        }
		        else{
			   xbin = (int)(floor((dx+contRad)/contBin));
		        }
		     }
		     else
			   xbin = (int)(dx/contBin);
// new print
		     if(angBin != 0 && angBin != 4){
		        dipMap[angBin][ybin][xbin] += 1.0;
			if(angBinsm > 0 && angBinsm < 4)
		           dipMapsm[angBin][ybin][xbin] += 1.0;
		     }
		     else if(angBin == 0 || angBin == 4){
			if(fabs(cosalp) > 0.99985) 
			   cosalp = sgn(cosalp)*0.99985;
			dipMap[angBin][ybin][xbin] += 1.0/sin(acos(cosalp));
		     	if(angBinsm == 0 || angBinsm == 4)
			   dipMapsm[angBin][ybin][xbin] += 1.0/sin(acos(cosalp));
		     }
// end new print
//		     fprintf(stderr,"%d, ybin= %d, xbin= %d\n",angBin*45,ybin,xbin);
		  }
	       }
	    }
	 }
       }
      }
   }
 }
}
//fprintf(stderr,"CFMForce.c: after polar stacking\n");
for	(i = 0; i < nCFM; i++){
	k = i*5;
/*****	Get intramolecular forces for i'th CFM.  ******/
	intraCFM(&pos[k],&force[k]);
/*****	Get intermolecular forces for i - j CFM-CFM interactions.  ******/
	for	(j = i+1; j < nCFM; j++,k = i*5){
		l = j*5;
		eg = 0.0;

/*****	Determine image vector for Ci -> Cj.  ******/
		image.fx = -(sdl.fx = pos[k].fx - pos[l].fx);
		image.fy = -(sdl.fy = pos[k].fy - pos[l].fy);
		image.fz = -(sdl.fz = pos[k].fz - pos[l].fz);
		mvimage(&sdl);
		image.fx += (delfx =sdl.fx);
		image.fy += (delfy =sdl.fy);
		image.fz += (delfz =sdl.fz);
		r2 = delfx*delfx + delfy*delfy + delfz*delfz;
		if	(r2 >= swr2max)
			continue;
		if	(r2 <= swr2min)
			{
			s = 1.0;
			sp = 0.0;
			}
		else
			swtch(r2-swr2min,&s,&sp);

		for	(m = 0; m < 10; m++)
			f[m].fx = f[m].fy = f[m].fz = 0.;

/*****	Loop over atoms in CFM molecules ******/
		for	(m = 0; m < 5; m++)
			for	(n = 0; n < 5; n++)
				{
/*****	Determine image vector ******/
				delfx = pos[k+m].fx - pos[l+n].fx + image.fx;
				delfy = pos[k+m].fy - pos[l+n].fy + image.fy;
				delfz = pos[k+m].fz - pos[l+n].fz + image.fz;
				r2 = delfx*delfx + delfy*delfy + delfz*delfz;

/*****	Get (1/r dV/dr) ******/
				dedr = ttraljq(r2,m,n,&eg);

/*****	Resolve forces on atoms.  ******/
				f[n+5].fx += (delfx *= dedr);
				f[n+5].fy += (delfy *= dedr);
				f[n+5].fz += (delfz *= dedr);
				f[m].fx -= delfx;
				f[m].fy -= delfy;
				f[m].fz -= delfz;
				}
		if	(sp != 0.0)
			{
			for	(m = 0; m < 10; m++)
				{
				f[m].fx *= s;
				f[m].fy *= s;
				f[m].fz *= s;
				}
			f[0].fx -= (frc.fx = sp*eg*sdl.fx);
			f[0].fy -= (frc.fy = sp*eg*sdl.fy);
			f[0].fz -= (frc.fz = sp*eg*sdl.fz);
			f[5].fx += frc.fx;
			f[5].fy += frc.fy;
			f[5].fz += frc.fz;
			eg *= s;
			}
		k = i*5;
		l = j*5;
		force[k].fx += f[0].fx;
		force[k].fy += f[0].fy;
		force[k].fz += f[0].fz;
		force[++k].fx += f[1].fx;
		force[k].fy += f[1].fy;
		force[k].fz += f[1].fz;
		force[++k].fx += f[2].fx;
		force[k].fy += f[2].fy;
		force[k].fz += f[2].fz;
		force[++k].fx += f[3].fx;
		force[k].fy += f[3].fy;
		force[k].fz += f[3].fz;
		force[++k].fx += f[4].fx;
		force[k].fy += f[4].fy;
		force[k].fz += f[4].fz;

		force[l].fx += f[5].fx;
		force[l].fy += f[5].fy;
		force[l].fz += f[5].fz;
		force[++l].fx += f[6].fx;
		force[l].fy += f[6].fy;
		force[l].fz += f[6].fz;
		force[++l].fx += f[7].fx;
		force[l].fy += f[7].fy;
		force[l].fz += f[7].fz;
		force[++l].fx += f[8].fx;
		force[l].fy += f[8].fy;
		force[l].fz += f[8].fz;
		force[++l].fx += f[9].fx;
		force[l].fy += f[9].fy;
		force[l].fz += f[9].fz;
		CFMNB += eg;
		}
	}
}

double
ttraljq(r2,m,n,eg)
double r2, *eg;
int m,n;
{
int bin,k;
double der,r,ir,ir6, ec,a,b,q,aa,bb;

	ir = 1. / r2;
	ir6 = ir * ir * ir;
	a = CFMlj[m][n].a;
	b = CFMlj[m][n].b;
	q = CFMlj[m][n].q;
	aa =12*a;
	bb = 6*b;
	r = sqrt(r2);
	bin = r/binRDF;
	if (m*n == 0){
		k = 0;
		if (m+n  > 0) k = 1 + (m+n)/4;
	}
	else {
		k = 3;
		if (n > 3 || m > 3) k =4;
		if (n ==4 && m == 4) k =5;
	}
	CFMRdf[k][bin] += 1.0;
	ec = q/r;
	*eg +=  ( a * ir6 - b ) * ir6 +ec;
	der = (bb - aa * ir6) * ir6 * ir -ec/r2;
	return(der);
}

/***	Calculate the intramolecular vibrational contribution using the
 ***	most general quadratic force field.
 ***/	


intraCFM(p,f)
tripd *p,*f;
{
tripd r[5],g[3];
double dr[5],df,r1r2,n1n2,da;
int i,j;
/* stretchings*/
for (i=1;i<5;i++){
	r[i].fx = p[0].fx - p[i].fx;
	r[i].fy = p[0].fy - p[i].fy;
	r[i].fz = p[0].fz - p[i].fz;
	dr[i] = sqrt(sq(r[i].fx) + sq(r[i].fy) + sq(r[i].fz));
	INTRAV += 0.5*ks[i]*sq(dr[i]-Req[i]);
	df = -ks[i]*(dr[i]-Req[i])/dr[i];
	f[i].fx -= df*r[i].fx;
	f[i].fy -= df*r[i].fy;
	f[i].fz -= df*r[i].fz;
	f[0].fx += df*r[i].fx;
	f[0].fy += df*r[i].fy;
	f[0].fz += df*r[i].fz;
}
/* bending */
for (i=1;i<4;i++){
   for (j=i+1;j<5;j++){
	r1r2 = r[i].fx*r[j].fx + r[i].fy*r[j].fy + r[i].fz*r[j].fz;
	n1n2 = r1r2 / (dr[i]*dr[j]);
	da = acos(n1n2) - Beq[i][j];
	INTRAV += 0.5*kb[i][j]*da*da;
	df = kb[i][j]*da/(sqrt(1.-n1n2*n1n2)*dr[i]*dr[j]);
	f[i].fx -= (g[1].fx = df*(r[j].fx - r1r2*r[i].fx/(dr[i]*dr[i])));
	f[i].fy -= (g[1].fy = df*(r[j].fy - r1r2*r[i].fy/(dr[i]*dr[i])));
	f[i].fz -= (g[1].fz = df*(r[j].fz - r1r2*r[i].fz/(dr[i]*dr[i])));
	f[j].fx -= (g[2].fx = df*(r[i].fx - r1r2*r[j].fx/(dr[j]*dr[j])));
	f[j].fy -= (g[2].fy = df*(r[i].fy - r1r2*r[j].fy/(dr[j]*dr[j])));
	f[j].fz -= (g[2].fz = df*(r[i].fz - r1r2*r[j].fz/(dr[j]*dr[j])));
	f[0].fx += g[1].fx + g[2].fx;
	f[0].fy += g[1].fy + g[2].fy;
	f[0].fz += g[1].fz + g[2].fz;
   }
}
}
getc3vec(int k){ // C-H & vector perpendicular to C3 axis

tripd Clcent,CHcent,r;
double d;
int i;
	tripd a,b;
	double ab, bb;
	double dx,dy,dz;
	double len;
	b.fx = pos[k+4].fx-pos[k].fx; // C->H
	b.fy = pos[k+4].fy-pos[k].fy;
	b.fz = pos[k+4].fz-pos[k].fz;
	a.fx = pos[k+1].fx-pos[k].fx; // C->Cl_1
	a.fy = pos[k+1].fy-pos[k].fy;
	a.fz = pos[k+1].fz-pos[k].fz;

	ab = a.fx*b.fx + a.fy*b.fy + a.fz*b.fz;
	bb = b.fx*b.fx + b.fy*b.fy + b.fz*b.fz;

	dx = a.fx - (ab/bb)*b.fx;
	dy = a.fy - (ab/bb)*b.fy;
	dz = a.fz - (ab/bb)*b.fz;
	len = sqrt(dx*dx+dy*dy+dz*dz);
	c3Vec[k/5][(int)(tc/TCFdt)].fx = dx/len;
	c3Vec[k/5][(int)(tc/TCFdt)].fy = dy/len;
	c3Vec[k/5][(int)(tc/TCFdt)].fz = dz/len;

	len = sqrt(bb);
//	len = bb;

	chVec[k/5][(int)(tc/TCFdt)].fx = b.fx/len;
	chVec[k/5][(int)(tc/TCFdt)].fy = b.fy/len;
	chVec[k/5][(int)(tc/TCFdt)].fz = b.fz/len;
//CHECK
/*
if(k==0)
   fprintf(stderr,"CFMForce.c: %f should be zero.\n",
	c3Vec[k/5][(int)(tc/TCFdt)].fx*chVec[k/5][(int)(tc/TCFdt)].fx +
	c3Vec[k/5][(int)(tc/TCFdt)].fy*chVec[k/5][(int)(tc/TCFdt)].fy +
	c3Vec[k/5][(int)(tc/TCFdt)].fz*chVec[k/5][(int)(tc/TCFdt)].fz);
*/
	r.fx=r.fy=r.fz=0.;
/*
Clcent.fx = (pos[k+1].fx + pos[k+2].fx + pos[k+3].fx)/3.;
Clcent.fy = (pos[k+1].fy + pos[k+2].fy + pos[k+3].fy)/3.;
Clcent.fz = (pos[k+1].fz + pos[k+2].fz + pos[k+3].fz)/3.;
// charges for 'pseudo-dipole' of CCl4 set in sysInit
CHcent.fx = pos[k].fx*QCFM[0]/(QCFM[0]+QCFM[4]) + pos[k+4].fx*QCFM[4]/(QCFM[0]+QCFM[4]);
CHcent.fy = pos[k].fy*QCFM[0]/(QCFM[0]+QCFM[4]) + pos[k+4].fy*QCFM[4]/(QCFM[0]+QCFM[4]);
CHcent.fz = pos[k].fz*QCFM[0]/(QCFM[0]+QCFM[4]) + pos[k+4].fz*QCFM[4]/(QCFM[0]+QCFM[4]);

r.fx = CHcent.fx - Clcent.fx;
r.fy = CHcent.fy - Clcent.fy;
r.fz = CHcent.fz - Clcent.fz;
*/
for(i=0;i<5;i++){
   r.fx += QCFM[i]*pos[k+i].fx;
   r.fy += QCFM[i]*pos[k+i].fy;
   r.fz += QCFM[i]*pos[k+i].fz;
}
d = sqrt(r.fx*r.fx+r.fy*r.fy+r.fz*r.fz);

dipVec[k/5][(int)(tc/TCFdt)].fx = r.fx/d;
dipVec[k/5][(int)(tc/TCFdt)].fy = r.fy/d;
dipVec[k/5][(int)(tc/TCFdt)].fz = r.fz/d;

}

/* check if dipoles are planar */
int inPlane(A1,A2,B1,B2)
tripd A1,A2,B1,B2;
{
   /* B1, B2 are imaged before calling inPlane */
   double r,dot,costh;
   double dotP;
   tripd B2A1,B2A2,B2B1,cross;
   tripd A1A2;
   B2A1.fx = A1.fx - B2.fx;
   B2A1.fy = A1.fy - B2.fy;
   B2A1.fz = A1.fz - B2.fz;
   B2A2.fx = A2.fx - B2.fx;
   B2A2.fy = A2.fy - B2.fy;
   B2A2.fz = A2.fz - B2.fz;
   A1A2.fx = A2.fx - A1.fx;
   A1A2.fy = A2.fy - A1.fy;
   A1A2.fz = A2.fz - A1.fz;
   r = sqrt(A1A2.fx*A1A2.fx + A1A2.fy*A1A2.fy + A1A2.fz*A1A2.fz);
   A1A2.fx /= r;
   A1A2.fy /= r;
   A1A2.fz /= r;

   B2B1.fx = B1.fx - B2.fx;
   B2B1.fy = B1.fy - B2.fy;
   B2B1.fz = B1.fz - B2.fz;
   r = sqrt(B2B1.fx*B2B1.fx + B2B1.fy*B2B1.fy + B2B1.fz*B2B1.fz);
   B2B1.fx /= r;
   B2B1.fy /= r;
   B2B1.fz /= r;
   /* get cross product A1B1 X A1B2 */
   cross.fx = B2A1.fy*B2A2.fz - B2A1.fz*B2A2.fy;
   cross.fy = B2A1.fz*B2A2.fx - B2A1.fx*B2A2.fz;
   cross.fz = B2A1.fx*B2A2.fy - B2A1.fy*B2A2.fx;
   r = sqrt(sq(cross.fx)+sq(cross.fy)+sq(cross.fz));
// if(r < 0.087156)
//    return(1);
   cross.fx/=r;
   cross.fy/=r;
   cross.fz/=r;
   
   dot = cross.fx*B2B1.fx + cross.fy*B2B1.fy + cross.fz*B2B1.fz;
   if(fabs(dot) < planCos/*0.087156*//*17365*/){ // if in plane, would be orthogonal to cross product
	   		    // (cosine = 0)
			    // 0.17365 is cos 80º, -0.17 is cos 100º
	return(1);
   }
   else{
	return(0);
   }
} 
   

int getAngBinsm(tripd A2A1, tripd icom, tripd jcom, tripd image)
{
   int i,j;
   tripd ij;
   double r,theta,costh;
   /* icom --> jcom vector */
   ij.fx = jcom.fx - icom.fx + image.fx;
   ij.fy = jcom.fy - icom.fy + image.fy;
   ij.fz = jcom.fz - icom.fz + image.fz;
   r = sqrt(sq(ij.fx)+sq(ij.fy)+sq(ij.fz));
   ij.fx/=r;
   ij.fy/=r;
   ij.fz/=r;

   costh = A2A1.fx*ij.fx + A2A1.fy*ij.fy + A2A1.fz*ij.fz;
   theta = acos(costh)*180/PI;
   if(costh > cos((dAngsm+0.0)*PI/180.0)){
      return 0;
   }
   else if(costh < cos((45.0-dAngsm)*PI/180.0) && costh  > cos((45.0+dAngsm)*PI/180.0)){
      return 1;
   }
   else if(costh < cos((90.0-dAngsm)*PI/180.0) && costh > cos((90.0+dAngsm)*PI/180.0)){
      return 2;
   }
   else if(costh < cos((135.0- dAngsm)*PI/180.0) && costh > cos((135.0+dAngsm)*PI/180.0)){
      return 3;
   }
   else if(costh < cos((180.0-dAngsm)*PI/180.0)){
      return 4;
   }
   else{
      return 5;
   }
}
int getAngBin(tripd A2A1, tripd icom, tripd jcom, tripd image)
{
   int i,j;
   tripd ij;
   double r,theta,costh;
   /* icom --> jcom vector */
   ij.fx = jcom.fx - icom.fx + image.fx;
   ij.fy = jcom.fy - icom.fy + image.fy;
   ij.fz = jcom.fz - icom.fz + image.fz;
   r = sqrt(sq(ij.fx)+sq(ij.fy)+sq(ij.fz));
   ij.fx/=r;
   ij.fy/=r;
   ij.fz/=r;

   costh = A2A1.fx*ij.fx + A2A1.fy*ij.fy + A2A1.fz*ij.fz;
   theta = acos(costh)*180/PI;
   if(costh > cos((dAng+0.0)*PI/180.0)){
      return 0;
   }
   else if(costh < cos((45.0-dAng)*PI/180.0) && costh  > cos((45.0+dAng)*PI/180.0)){
      return 1;
   }
   else if(costh < cos((90.0-dAng)*PI/180.0) && costh > cos((90.0+dAng)*PI/180.0)){
      return 2;
   }
   else if(costh < cos((135.0- dAng)*PI/180.0) && costh > cos((135.0+dAng)*PI/180.0)){
      return 3;
   }
   else if(costh < cos((180.0-dAng)*PI/180.0)){
      return 4;
   }
   else{
      return 5;
   }
}


double headTail(tripd A1,tripd A2, tripd B1, tripd B2)
{
   /* 
    * uses formula for distance of point (x0) to 
    * line that passes through 2 points (x1,x2):
    * d = |(x1x2 X x0x1 )| / |x1x2|
    * where x1x2 is line segment x1 to x2
    */
   double rB1, rB2, rd, rc;
   tripd d1,d2,d3,crossp;
   
   d1.fx = B1.fx - A1.fx; 
   d1.fy = B1.fy - A1.fy;
   d1.fz = B1.fz - A1.fz;
   d2.fx = B1.fx - A2.fx;
   d2.fy = B1.fy - A2.fy;
   d2.fz = B1.fz - A2.fz;
   d3.fx = A2.fx - A1.fx;
   d3.fy = A2.fy - A1.fy;
   d3.fz = A2.fz - A1.fz;
   rd = sqrt(sq(d3.fx)+sq(d3.fy)+sq(d3.fz));

   crossp.fx = d1.fy*d2.fz - d1.fz*d2.fy;
   crossp.fy = d1.fz*d2.fx - d1.fx*d2.fz;
   crossp.fz = d1.fx*d2.fy - d1.fy*d2.fx;
   rc = sqrt(sq(crossp.fx)+sq(crossp.fy)+sq(crossp.fz));

   rB1 = rc/rd;
   
   d1.fx = B2.fx - A1.fx;
   d1.fy = B2.fy - A1.fy;
   d1.fz = B2.fz - A1.fz;
   d2.fx = B2.fx - A2.fx;
   d2.fy = B2.fy - A2.fy;
   d2.fz = B2.fz - A2.fz;

   crossp.fx = d1.fy*d2.fz - d1.fz*d2.fy;
   crossp.fy = d1.fz*d2.fx - d1.fx*d2.fz;
   crossp.fz = d1.fx*d2.fy - d1.fy*d2.fx;
   rc = sqrt(sq(crossp.fx)+sq(crossp.fy)+sq(crossp.fz));

   rB2 = rc/rd;

   if(rB2>rB1){
      return 360.0;
   }
   else{
      return 0.0;
   }
}

