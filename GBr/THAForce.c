#include	<md.h>
#include	<system.h>
#include	 <math.h>
#define EXCON 0.1
#define KAP_ALPHA 0.01

THAForce()
{
int i, k;
double r;
double ttraljq();
int THAkappa(),centerTHA();

int nw;
nw = natoms - nGLY*GLYsites - nBr - 2*nCl2 - 3*nTS - nTHA*THAsites;

tbondE = tbendE = ttorsE = t14E = t15E = 0.;
for(i=0;i<nTHA;i++){
   k = nGLY*GLYsites + nw;

   intraTHA(k);
   THAkappa();
}
//calculate the cavity constant kappa

/*if(nTHA==1){
   k = nGLY*GLYsites+nw;
   THAcent.fx = THAcent.fy = THAcent.fz = 0.0;
   kappa = 0.0;
   NTHAvec.fx = NTHAvec.fy = NTHAvec.fz = 0.0;
   for(i=0;i<THAsites;i++){
      THAcent.fx += pos[k+i].fx;
      THAcent.fy += pos[k+i].fy;
      THAcent.fz += pos[k+i].fz;
   }

   THAcent.fx /= (double)(THAsites);
   THAcent.fy /= (double)(THAsites);
   THAcent.fz /= (double)(THAsites);
   NTHAvec.fx = THAcent.fx - pos[k].fx;
   NTHAvec.fy = THAcent.fy - pos[k].fy;
   NTHAvec.fz = THAcent.fz - pos[k].fz;

   r = sqrt(NTHAvec.fx*NTHAvec.fx + NTHAvec.fy*NTHAvec.fy + NTHAvec.fz*NTHAvec.fz);

   NTHAvec.fx /= r;
   NTHAvec.fy /= r;
   NTHAvec.fz /= r;

   kappa = r/kappaNorm;
}
*/
}

intraTHA(k)
int k;
{
int i,j,n, n1, n2, n3, m;
double dx, dy, dz, r2;
double tcalcstrch(), tcalcbend(), tcalctorq(),ttraljq14(), ttraljq15();  
double tbond_e, ttors_e, tbend_e, t15_e, t14_e;
double e14,c14,e15,c15;
double dedr;

/* loop over all 1-5+ pairs */
t15_e = 0.;
for(i=0;i<nTHA15;i++){
   e15 = c15 = 0.;
   dx=pos[k+THA15[i][0]].fx-pos[k+THA15[i][1]].fx;
   dy=pos[k+THA15[i][0]].fy-pos[k+THA15[i][1]].fy;
   dz=pos[k+THA15[i][0]].fz-pos[k+THA15[i][1]].fz;
   r2=dx*dx + dy*dy + dz*dz;
   dedr = ttraljq15(r2,THA15[i][0],THA15[i][1],&e15,&c15);
// ****	Resolve forces on atoms.  ******
   force[k+THA15[i][1]].fx += (dx *= dedr);
   force[k+THA15[i][1]].fy += (dy *= dedr);
   force[k+THA15[i][1]].fz += (dz *= dedr);
   force[k+THA15[i][0]].fx -= dx;
   force[k+THA15[i][0]].fy -= dy;
   force[k+THA15[i][0]].fz -= dz;
//   if(tc==0)
//      fprintf(stderr,"1-5 i=%02d (%d,%d) U = %f\n",i,THA15[i][0],THA15[i][1],e15*KCAL);
   t15_e += e15;
}

/* loop over all 1-4 pairs */
t14_e = 0.;
// intra-unit 1-4's
for(j=0;j<4;j++){
   for(i=0;i<nTHAtors1;i++){
      n1 = (THAtors1[i][0]==0)?0:j*19+THAtors1[i][0];
      n2 = (THAtors1[i][3]==0)?0:j*19+THAtors1[i][3];
      dx=pos[k+n1].fx-pos[k+n2].fx;
      dy=pos[k+n1].fy-pos[k+n2].fy;
      dz=pos[k+n1].fz-pos[k+n2].fz;
      r2=dx*dx + dy*dy + dz*dz;
      e14 = c14 = 0.;
      dedr = ttraljq14(r2,THAtors1[i][0],THAtors1[i][3],&e14,&c14);
// ****	Resolve forces on atoms.  ******
      force[k+n2].fx += (dx *= dedr);
      force[k+n2].fy += (dy *= dedr);
      force[k+n2].fz += (dz *= dedr);
      force[k+n1].fx -= dx;
      force[k+n1].fy -= dy;
      force[k+n1].fz -= dz;
//   if(tc==0)
//      fprintf(stderr,"1-4 i=%02d (%d,%d) U = %f\n",i,THAtors[i][0],THAtors[i][3],e14*KCAL);
      t14_e += e14;
   }
}

// inter-unit 1-4's
for(i=0;i<nTHAtors2;i++){
   n1 = THAtors2[i][0];
   n2 = THAtors2[i][3];
   dx=pos[k+n1].fx-pos[k+n2].fx;
   dy=pos[k+n1].fy-pos[k+n2].fy;
   dz=pos[k+n1].fz-pos[k+n2].fz;
   r2=dx*dx + dy*dy + dz*dz;
   e14 = c14 = 0.;
   dedr = ttraljq14(r2,THAtors2[i][0],THAtors2[i][3],&e14,&c14);
// ****	Resolve forces on atoms.  ******
   force[k+n2].fx += (dx *= dedr);
   force[k+n2].fy += (dy *= dedr);
   force[k+n2].fz += (dz *= dedr);
   force[k+n1].fx -= dx;
   force[k+n1].fy -= dy;
   force[k+n1].fz -= dz;
//   if(tc==0)
//      fprintf(stderr,"1-4 i=%02d (%d,%d) U = %f\n",i,THAtors[i][0],THAtors[i][3],e14*KCAL);
      t14_e += e14;
}

/* loop over all bond stretches*/
tbond_e = 0.;
for(j=0;j<4;j++){
   for(i=0;i<nTHAstr;i++){
      n = THAstr[i][0]==0?k:j*19+k+THAstr[i][0];
      m = THAstr[i][1]==0?k:j*19+k+THAstr[i][1];
      tbond_e += tcalcstrch(n,m,THAstr[i][2]);
   }
}

/* loop over all torsions*/
ttors_e=0.;
/*i-j-k-l*/
//intra-unit torsions
for(j=0;j<4;j++){
   for(i=0;i<nTHAtors1;i++){
      n = (THAtors1[i][0]==0)?k:j*19+k+THAtors1[i][0];//i
      n1= (THAtors1[i][1]==0)?k:j*19+k+THAtors1[i][1];//j
      n2= (THAtors1[i][2]==0)?k:j*19+k+THAtors1[i][2];//k
      n3= (THAtors1[i][3]==0)?k:j*19+k+THAtors1[i][3];//l
      ttors_e += tcalctorq(n,n1,n2,n3,THAtors1[i][4]);
   }
}
//inter-unit torsions
for(i=0;i<nTHAtors2;i++){
   n = k+THAtors2[i][0];//i
   n1= k+THAtors2[i][1];//j
   n2= k+THAtors2[i][2];//k
   n3= k+THAtors2[i][3];//l
   ttors_e += tcalctorq(n,n1,n2,n3,THAtors2[i][4]);
}

/* loop over all angle bends*/
tbend_e = 0.;
/*i-j-k*/
//intra-unit bends
for(j=0;j<4;j++){
   for(i=0;i<nTHAbend1;i++){
      n = (THAbend1[i][0]==0)?k:j*19+k+THAbend1[i][0];//i
      n1= (THAbend1[i][1]==0)?k:j*19+k+THAbend1[i][1];//j
      n2= (THAbend1[i][2]==0)?k:j*19+k+THAbend1[i][2];//k
      tbend_e += tcalcbend(n,n1,n2,THAbend1[i][3]);
   }
}
//intra-unit bends
for(i=0;i<nTHAbend2;i++){
   n = k+THAbend2[i][0];//i
   n1= k+THAbend2[i][1];//j
   n2= k+THAbend2[i][2];//k
   tbend_e += tcalcbend(n,n1,n2,THAbend2[i][3]);
}

/* total sum for all molecules*/
tbondE += tbond_e;
tbendE += tbend_e;
ttorsE += ttors_e;
t14E += t14_e;
t15E += t15_e;
/* total intramolecular potential energy*/
INTRA_THA += tbond_e + ttors_e + tbend_e + t14_e + t15_e;
}

double tcalcstrch(i,j,k)
int i,j,k;
{
tripd 	d; 
double	bond, dedr, pe;

d.fx = pos[i].fx - pos[j].fx;
d.fy = pos[i].fy - pos[j].fy;
d.fz = pos[i].fz - pos[j].fz;
bond = sqrt(d.fx*d.fx+d.fy*d.fy+d.fz*d.fz);
dedr = tkstr[k]*(bond-treq[k])/bond;
pe = 0.5*tkstr[k]*sq(bond-treq[k]);
force[i].fx -= (d.fx = dedr*d.fx);
force[i].fy -= (d.fy = dedr*d.fy);
force[i].fz -= (d.fz = dedr*d.fz);
force[j].fx += d.fx;
force[j].fy += d.fy;
force[j].fz += d.fz;
return(pe);

}
double tcalctorq (i1,i2,i3,i4,m)
int i1,i2,i3,i4,m;
{
int i,j,k,index,ats[4];
double dmat[3][3], 
      cmat[3][3],
      grad[4][3], 
      d[3][3], 
      da,cosa, 
      K1,K10,K11,K12,K13,K21,K22,K23, 
      D1201, C1100, 
      dudcos, 
      pe,V1,V2;
ats[0] = i1;
ats[1] = i2;
ats[2] = i3;
ats[3] = i4;
    for (i=0;  i<3;  i++){/* get vectors for atoms */
        d[i][0]=pos[ats[i+1]].fx-pos[ats[i]].fx;
        d[i][1]=pos[ats[i+1]].fy-pos[ats[i]].fy;
        d[i][2]=pos[ats[i+1]].fz-pos[ats[i]].fz;
    }
    for (i=0; i<3;  i++){
        for (j=0;  j<3;  j++){
	    cmat[i][j]=0.0;
            for (k=0;  k<3;  k++){
                cmat[i][j] += d[i][k]*d[j][k];
            }
        }
    }
    for (i=0; i<3;  i++){
        for (j=0;  j<3;  j++){
            dmat[i][j] = cmat[i][i]*cmat[j][j]-cmat[i][j]*cmat[i][j];
        }
    }

    
    K1 = cmat[0][1]*cmat[1][2] - cmat[0][2]*cmat[1][1];
    D1201 = sqrt (dmat[1][2]*dmat[0][1]);

    if (D1201 == 0){
	fprintf(stderr,"No torsion defined\n");
	exit(1);
    }


    /*
    **	    minus gradient at atom1 with respect to the cosine of the angle
    */

    for (i=0;  i<3;  i++){
        grad[0][i]=( ( (K1* (cmat[0][1] * d[1][i]-cmat[1][1] * d[0][i]) )
	/dmat[0][1]) + cmat[1][2]*d[1][i] - cmat[1][1]*d[2][i])/D1201;
    }

    /*
    **	    minus gradient at atom2 with respect to the cosine of the angle
    */

    K11 = cmat[0][1] + cmat[1][1];
    K12 = cmat[1][2] + 2*cmat[0][2];
    K13 = cmat[0][0] + cmat[0][1];
   
    for (i=0;  i<3;  i++){
	V1 = cmat[1][2]*d[2][i] - cmat[2][2]*d[1][i];
	V2 = K11*d[0][i] - K13*d[1][i];
        grad[1][i] = (K11*d[2][i] - K12*d[1][i] + cmat[1][2]*d[0][i] + 
			K1*(V1/dmat[1][2] + V2/dmat[0][1]))/D1201;
    }

    /*
    **	   minus gradient at atom3 with respect to the cosine of the angle
    */
    K21= cmat[0][1] + 2*cmat[0][2];
    K22= cmat[1][2] + cmat[1][1];
    K23= cmat[2][2] + cmat[1][2];
    for (i=0;  i<3;  i++){
        V1 = cmat[0][0]*d[1][i] - cmat[0][1]*d[0][i];
	V2 = K23*d[1][i] - K22*d[2][i];
	grad[2][i] = (K21*d[1][i] - K22*d[0][i] - cmat[0][1]*d[2][i] +
			K1*(V1/dmat[0][1] + V2/dmat[1][2]))/D1201;
    }

    /*
    **	    minus gradient at atom4 with respect to the cosine of the angle
    */
    for (i=0;  i<3;  i++){
        grad[3][i] = (-cmat[0][1]*d[1][i] + cmat[1][1]*d[0][i] + 
	    K1*(cmat[1][1]*d[2][i] - cmat[1][2]*d[1][i])/dmat[1][2])/D1201;
    }
    /*
    **	    calculate cos angle
    */
    cosa = K1/D1201;
    /*IMPA[m] = cosa;*/
/*
    index = acos(cosa)*36/PI;
    gTors[l][index]++;
*/
    /*
    **	    evaluate dU/dcos(a)
    **  OPLS-style
    **  U = pTors[l][0]*(1+cos(a))+pTors[l][1]*(1-cos(2*a))+pTors[l][2]*(1+cos(3*a))
    */
//    pe = (1+cosa)*(gtors[m][0]+2*gtors[m][1]*(1-cosa)+gtors[m][2]*sq(2*cosa-1));
//    dudcos = gtors[m][0] - 4*gtors[m][1]*cosa + 3*gtors[m][2]*(4*cosa*cosa-1);
    /*
    **	evaluate dU/dcos(a)
    **  AMBER-style
    **  U = pTors[l][0]*(1+cos(a-phi))+pTors[l][1]*(1+cos(2*a-phi))+pTors[l][2]*(1+cos(3*a-phi))
    **  (let phi = 0)
    */
    //pe = ttors[m][1]*2.0*cosa*cosa + (1+cosa)*(ttors[m][0]+ttors[m][2]*sq(2*cosa-1));
    //dudcos = ttors[m][0] + 4*ttors[m][1]*cosa + 3*ttors[m][2]*(4*cosa*cosa-1);
    if(m!=2){ // phi1 = phi2 = phi3 = 0º
       pe = ttors[m][0]*(1.+cosa) + 2.*ttors[m][1]*cosa*cosa + ttors[m][2]*(cosa+1.)*(2.*cosa-1.)*(2.*cosa-1.);
       dudcos = ttors[m][0] + 4.*ttors[m][1]*cosa + ttors[m][2]*(12.*cosa*cosa-3.);
    }
    else if(m==2){ // phi1 = phi2 = 180º, phi3 = 0º
       pe = ttors[m][0]*(1.-cosa) + 2.*ttors[m][1]*(1.-cosa*cosa) + ttors[m][2]*(cosa+1.)*(2.*cosa-1.)*(2.*cosa-1.);
       dudcos = -ttors[m][0] - 4.*ttors[m][1]*cosa + ttors[m][2]*(12.*cosa*cosa-3.);
    }
    /*
    **	    add resulting force to atoms for angle
    */
    for(i=0;i<4;i++){
        force[ats[i]].fx += dudcos*grad[i][0];
	force[ats[i]].fy += dudcos*grad[i][1];
	force[ats[i]].fz += dudcos*grad[i][2];
    }
    return(pe);
}

double tcalcbend(i,j,k,l)
int i,j,k,l;
{
int i0,i1,i2;
tripd r[2], grad[3];
double pe,r1,r2,r1r2,n1n2,da,dedd;

i0 = j;/*center atom*/
i1 = i;
i2 = k;

r[0].fx = pos[i0].fx - pos[i1].fx;
r[0].fy = pos[i0].fy - pos[i1].fy;
r[0].fz = pos[i0].fz - pos[i1].fz;
r[1].fx = pos[i0].fx - pos[i2].fx;
r[1].fy = pos[i0].fy - pos[i2].fy;
r[1].fz = pos[i0].fz - pos[i2].fz;

r1 = sqrt(r[0].fx*r[0].fx + r[0].fy*r[0].fy + r[0].fz*r[0].fz);
r2 = sqrt(r[1].fx*r[1].fx + r[1].fy*r[1].fy + r[1].fz*r[1].fz);
r1r2 =   (r[0].fx*r[1].fx + r[0].fy*r[1].fy + r[0].fz*r[1].fz);
n1n2 = r1r2 / (r1*r2);
da = acos(n1n2) - teqbend[l];
pe = 0.5*tkbend[l]*da*da;
dedd = tkbend[l]*da/ (sqrt(1.-n1n2*n1n2)*r1*r2);
force[i1].fx -= (grad[1].fx = dedd*(r[1].fx - r1r2*r[0].fx/(r1*r1)));
force[i1].fy -= (grad[1].fy = dedd*(r[1].fy - r1r2*r[0].fy/(r1*r1)));
force[i1].fz -= (grad[1].fz = dedd*(r[1].fz - r1r2*r[0].fz/(r1*r1)));
force[i2].fx -= (grad[2].fx = dedd*(r[0].fx - r1r2*r[1].fx/(r2*r2)));
force[i2].fy -= (grad[2].fy = dedd*(r[0].fy - r1r2*r[1].fy/(r2*r2)));
force[i2].fz -= (grad[2].fz = dedd*(r[0].fz - r1r2*r[1].fz/(r2*r2)));
force[i0].fx += grad[1].fx + grad[2].fx;
force[i0].fy += grad[1].fy + grad[2].fy;
force[i0].fz += grad[1].fz + grad[2].fz;
    return(pe);
}
double
ttraljq15(r2,m,n,eg,nc)
double r2, *eg, *nc;
int m,n;
{
double der,r,ir,ir6, ec,a,b,q,aa,bb,sqrt();

	ir = 1. / r2;
	ir6 = ir * ir * ir;
	a = THAlj[m][n].a;
	b = THAlj[m][n].b;

	q = THAlj[m][n].q;
	aa =12*a;
	bb = 6*b;

	r = sqrt(r2);
	ec = q/r;

	*eg +=  ( a * ir6 - b ) * ir6 + ec;
	*nc += ec;
	der = (bb - aa * ir6) * ir6 * ir - ec/r2;
	return(der);
}


double
ttraljq14(r2,m,n,eg,nc)
double r2, *eg, *nc;
int m,n;
{
double der,r,ir,ir6, ec,a,b,q,aa,bb,sqrt();

	ir = 1. / r2;
	ir6 = ir * ir * ir;
	a = THAlj[m][n].a;
	b = THAlj[m][n].b;

	q = THAlj[m][n].q;
	aa =12*a;
	bb = 6*b;

	r = sqrt(r2);
	ec = q/r;

	*eg += (t14ljscale)*(( a * ir6 - b ) * ir6) + (t14qscale)*ec;
	*nc += (t14qscale)*ec;
	der  = (t14ljscale)*((bb - aa * ir6) * ir6 * ir) - (t14qscale)*ec/r2;
	return(der);
}

THAkappa()  //calculate CoM, radius of gyration, etc for THA
{
   int i,j,k,l,nw;
   double totMass,r2;
   double kx,ky,kz;
   double gamma,deriv;
   RgTHA = 0.0;
   THAcom.fx = THAcom.fy = THAcom.fz = totMass = 0.0;
   THANz = 0.0;
   nw = natoms - nGLY*GLYsites - nBr - 2*nCl2 - 3*nTS - nTHA*THAsites;
   kappa = 0.0;
   if(nTHA !=1)
      return 0;
// radius of gyration
   for(i=0;i<4;i++){
      for(j=i+1;j<4;j++){
	 k = nw+nGLY*GLYsites+1+19*i;
	 l = nw+nGLY*GLYsites+1+19*j;
	 r2 = sq(pos[k].fx-pos[l].fx) + sq(pos[k].fy-pos[l].fy) + sq(pos[k].fz-pos[l].fz);
         RgTHA += (1./32.)*r2;
      }
   }
// Center of Mass, kappa
   for(i=0;i<THAsites;i++){
      k=nGLY*GLYsites+i;
      totMass+=mass[k];
      THAcom.fx += pos[k].fx*mass[k];
      THAcom.fy += pos[k].fy*mass[k];
      THAcom.fz += pos[k].fz*mass[k];
   }
   THAcom.fx/=totMass;
   THAcom.fy/=totMass;
   THAcom.fz/=totMass;

   THANz = pos[nGLY*GLYsites+nw].fz;

   kx = THAcom.fx - pos[k].fx;
   ky = THAcom.fy - pos[k].fy;
   kz = THAcom.fz - pos[k].fz;
   kappa = sqrt(kx*kx + ky*ky + kz*kz);

// kappa window potential
   if(kappaWindowOn==1){
   gamma = fabs(kappa - kap_c) - kap_w;
   if(gamma > 0.0){
      VkapWin += KAP_ALPHA*gamma*gamma*gamma;
      deriv = -3.0*KAP_ALPHA*gamma*gamma*sgn(kappa - kap_c);
      k = nw + nGLY*GLYsites;
      //force on N
      force[k].fx+=deriv*(mass[k]/totMass - 1.0)*(THAcom.fx-pos[k].fx);
      force[k].fy+=deriv*(mass[k]/totMass - 1.0)*(THAcom.fy-pos[k].fy);
      force[k].fz+=deriv*(mass[k]/totMass - 1.0)*(THAcom.fz-pos[k].fz);
      //force on non-N atoms
      for(i=1;i<THAsites;i++){
	      force[k+i].fx+=deriv*(mass[k+i]/totMass)*(THAcom.fx-pos[k+i].fx);
	      force[k+i].fy+=deriv*(mass[k+i]/totMass)*(THAcom.fy-pos[k+i].fy);
	      force[k+i].fz+=deriv*(mass[k+i]/totMass)*(THAcom.fz-pos[k+i].fz);
      }
   }
   }
// kappa biasing potential
   if(kappaBiasingOn==1){
   } 	   
}

centerTHA() //shift system so that THA N is at x=0,y=0
{
   int i,k;
   tripd cent;
   if(nTHA!=1 && nBr+nCl2+nTS !=1)
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
