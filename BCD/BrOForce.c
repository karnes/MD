#include <md.h>
#include <system.h>
#include <math.h>

BrOForce()
{
int	i, j, k, l, m, n,nw,d;
double	r2, dedr, delfx, delfy, delfz;
double	eg, s, sp;
tripd	del, frc, image, sdl, f[2*BrOs];
double	BrOlj();
double	BrOq();
//fprintf(stderr,"swr2max = %f, swr2min = %f, swsomin = %f, swsomax = %f\n",swr2min,swr2max,swSoMin,swSoMax);
nw = natoms - nBCD*BCDs - nBrO*BrOs - nsolute;
BrONB = BrObondE = BrOtorsE = BrObendE = BrOnb14E = BrOnb15E = 0.0;

for(i = 0; i < nBrO; i++){
	k = nw+i*BrOs;
/*****	Get intramolecular forces for i'th bromo-octane.  ******/
	intraBrO(k);
/*****	Get intermolecular forces for i - j BrO-BrO interactions.  ******/
	for(j =i+1; j < nBrO; j++){
		k = nw+i*BrOs;
		l = nw+j*BrOs;
/*****	Resolve L-J forces ******/
/*****	Loop over atoms in BrOctane molecules ******/
		for(m = 0; m < BrOs; m++)
		{
		   for(n = 0; n < BrOs; n++)
		   {
		      eg = 0.0;
/*****	Determine image vector ******/
		      del.fx = pos[k+m].fx - pos[l+n].fx;
		      del.fy = pos[k+m].fy - pos[l+n].fy;
		      del.fz = pos[k+m].fz - pos[l+n].fz;
		      mvimage(&del);
		      r2 = del.fx*del.fx + del.fy*del.fy + del.fz*del.fz;
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
		      dedr = BrOlj(r2,m,n,&eg);
/*****	Resolve forces on atoms.  ******/
		      force[l+n].fx += (del.fx *= (dedr*s+sp*eg));
		      force[l+n].fy += (del.fy *= (dedr*s+sp*eg));
		      force[l+n].fz += (del.fz *= (dedr*s+sp*eg));
		      force[k+m].fx -= del.fx;
		      force[k+m].fy -= del.fy;
		      force[k+m].fz -= del.fz;
		      BrONB += eg*s;
		   }
		}
		k = nw+i*BrOs;
		l = nw+j*BrOs;
		eg = 0.0;
/*****	Determine coulombic forces between Br-C headgroups.  ******/
/*****	Determine image vector for i - j headgroups.  ******/
		image.fx = -(sdl.fx = pos[k].fx - pos[l].fx);
		image.fy = -(sdl.fy = pos[k].fy - pos[l].fy);
		image.fz = -(sdl.fz = pos[k].fz - pos[l].fz);
		mvimage(&sdl);
		image.fx += (delfx =sdl.fx);
		image.fy += (delfy =sdl.fy);
		image.fz += (delfz =sdl.fz);
		r2 = delfx*delfx + delfy*delfy + delfz*delfz;
		if(r2 >= swr2max)
		   continue;
		if(r2 <= swr2min)
		{
		   s = 1.0;
		   sp = 0.0;
		}
		else
		   swtch(r2-swr2min,&s,&sp);
		for(d = 0; d < 4; d++)
		   f[d].fx = f[d].fy = f[d].fz = 0.;
/*****	Loop over atoms in BrO head-groups (Br-C) ******/
		for(m = 0; m < 2; m++)
		{
			for(n = 0; n < 2; n++)
			{
/*****	Determine image vector ******/
			   delfx = pos[k+m].fx - pos[l+n].fx + image.fx;
			   delfy = pos[k+m].fy - pos[l+n].fy + image.fy;
			   delfz = pos[k+m].fz - pos[l+n].fz + image.fz;
			   r2 = delfx*delfx + delfy*delfy + delfz*delfz;

/*****	Get (1/r dV/dr) ******/
	 		   dedr = BrOq(r2,m,n,&eg);
/*****	Resolve forces on atoms.  ******/
			   f[2+n].fx += (delfx *= dedr);
			   f[2+n].fy += (delfy *= dedr);
			   f[2+n].fz += (delfz *= dedr);
			   f[m].fx -= delfx;
			   f[m].fy -= delfy;
			   f[m].fz -= delfz;
			}
		}
		if(sp != 0.0){
		   for(m = 0; m < 4; m++)
		   {
			f[m].fx *= s;
			f[m].fy *= s;
			f[m].fz *= s;
		   }
		   f[0].fx -= (frc.fx = sp*eg*sdl.fx);
		   f[0].fy -= (frc.fy = sp*eg*sdl.fy);
		   f[0].fz -= (frc.fz = sp*eg*sdl.fz);
		   f[2].fx += frc.fx;
		   f[2].fy += frc.fy;
		   f[2].fz += frc.fz;
		   eg *= s;
		}
		k = nw+i*BrOs;
		l = nw+j*BrOs;
		force[k].fx += f[0].fx;
		force[k].fy += f[0].fy;
		force[k].fz += f[0].fz;
		force[++k].fx += f[1].fx;
		force[k].fy += f[1].fy;
		force[k].fz += f[1].fz;

		force[l].fx += f[2].fx;
		force[l].fy += f[2].fy;
		force[l].fz += f[2].fz;
		force[++l].fx += f[3].fx;
		force[l].fy += f[3].fy;
		force[l].fz += f[3].fz;
		BrONB += eg;
	}
		
}

}

double
BrOlj(r2,m,n,eg)
double r2, *eg;
int m,n;
{
double q2,der,ir,ir6,a,b,aa,bb;

	ir = 1. / r2;
	ir6 = ir * ir * ir;
	a = slj[3+m][3+n].a;
	b = slj[3+m][3+n].b;
	aa =12*a;
	bb = 6*b;
	*eg +=  ( a * ir6 - b ) * ir6;
	der = (bb - aa * ir6) * ir6 * ir;
	return(der);
}
double
BrOq(r2,m,n,eg)
double r2, *eg;
int m,n;
{
double q2,r,der,/*ir,ir6, a,b,aa,bb,*/EC;

//	ir = 1. / r2;
//	ir6 = ir * ir * ir;
//	a = slj[3+m][3+n].a;
//	b = slj[3+m][3+n].b;
//	aa =12*a;
//	bb = 6*b;
	q2 = slj[3+m][3+n].q;
	r = sqrt(r2);
	EC = q2/r;
	*eg += EC;
	der = -EC/r2;
	return(der);
}


/***	Calculate the contributions of stretch, bend, torsion and 1-4
 ***	van der Waals to forces and energy for a 6 atoms chain molecule.
 ***/	


intraBrO(fa)
int fa;
{
double BrOcalcstrch(), BrOcalcbend(), BrOcalctorq(), calcnb(),BrOnb15();  
double bond_e, tors_e, bend_e, nb14_e, nb15_e;
int i;
/* connectivity table:
  C1------C2------C3------C4------C5------C6
 fa+4    fa+2     fa      fa+1   fa+3    fa+5
  Br---C1-----C2-----C3-----C4-----C5----C6-----C7-----C8
 fa+1  fa    fa+2   fa+3   fa+4   fa+5  fa+6   fa+7   fa+8
*/
tors_e = nb14_e = nb15_e = bond_e = bend_e = 0.0;
/* loop over all bond stretches*/
    	bond_e = BrOcalcstrch(fa+1,fa,0);
    	bond_e += BrOcalcstrch(fa,fa+2,1);
    	bond_e += BrOcalcstrch(fa+2,fa+3,2);
    	bond_e += BrOcalcstrch(fa+3,fa+4,2);
    	bond_e += BrOcalcstrch(fa+4,fa+5,2);
    	bond_e += BrOcalcstrch(fa+5,fa+6,2);
    	bond_e += BrOcalcstrch(fa+6,fa+7,2);
    	bond_e += BrOcalcstrch(fa+7,fa+8,3);
/* loop over all torsions*/

    	tors_e = BrOcalctorq(fa+1,fa,fa+2,fa+3,0);
    	tors_e += BrOcalctorq(fa,fa+2,fa+3,fa+4,1);
    	tors_e += BrOcalctorq(fa+2,fa+3,fa+4,fa+5,2);
   	tors_e += BrOcalctorq(fa+3,fa+4,fa+5,fa+6,2);
    	tors_e += BrOcalctorq(fa+4,fa+5,fa+6,fa+7,2);
    	tors_e += BrOcalctorq(fa+5,fa+6,fa+7,fa+8,3);

/* loop over all angle bends*/
    	bend_e = BrOcalcbend(fa+1,fa,fa+2,0);
//fprintf(stderr,"bend_e = %f\n",bend_e*KCAL);
    	bend_e += BrOcalcbend(fa,fa+2,fa+3,1);
//fprintf(stderr,"bend_e = %f\n",bend_e*KCAL);
    	bend_e += BrOcalcbend(fa+2,fa+3,fa+4,2);
//fprintf(stderr,"bend_e = %f\n",bend_e*KCAL);
    	bend_e += BrOcalcbend(fa+3,fa+4,fa+5,2);
//fprintf(stderr,"bend_e = %f\n",bend_e*KCAL);
    	bend_e += BrOcalcbend(fa+4,fa+5,fa+6,2);
//fprintf(stderr,"bend_e = %f\n",bend_e*KCAL);
    	bend_e += BrOcalcbend(fa+5,fa+6,fa+7,2);
//fprintf(stderr,"bend_e = %f\n",bend_e*KCAL);
    	bend_e += BrOcalcbend(fa+6,fa+7,fa+8,3);
//fprintf(stderr,"bend_e = %f\n",bend_e*KCAL);
/* loop over 1-4 non-bonded energies*/
    	nb14_e = calcnb(fa+1,fa+3,1,3);
//    	nb14_e += calcnb(fa,fa+4,0,4);
//    	nb14_e += calcnb(fa+2,fa+5,2,5);
//    	nb14_e += calcnb(fa+3,fa+6,3,6);
//    	nb14_e += calcnb(fa+4,fa+7,4,7);
//    	nb14_e += calcnb(fa+5,fa+8,5,8);
/* loop over 1-5+ non-bonded energies*/
	for(i=4;i<9;i++)
    	    nb15_e += BrOnb15(fa+1,fa+i,1,i);
	for(i=5;i<9;i++)
    	    nb15_e += BrOnb15(fa,fa+i,0,i);
	for(i=6;i<9;i++)
    	    nb15_e += BrOnb15(fa+2,fa+i,2,i);
	for(i=7;i<9;i++)
    	    nb15_e += BrOnb15(fa+3,fa+i,3,i);
    	nb15_e += BrOnb15(fa+4,fa+8,4,8);

/* total sum for all molecules*/
	BrObondE += bond_e;
	BrObendE += bend_e;
	BrOtorsE += tors_e;
	BrOnb14E += nb14_e;
	BrOnb15E += nb15_e;
/* total intramolecular potential energy*/
    	BrOV += bond_e + tors_e + bend_e + nb14_e + nb15_e;
//fprintf(stderr,"fa = %d, BrOV = %f, %f, %f, %f, %f, %f\n",fa,BrOV*KCAL,bond_e*KCAL,bend_e*KCAL,tors_e*KCAL,nb14_e*KCAL,nb15_e*KCAL);
}

double BrOcalcstrch(i,j,k)
int i,j,k;
{
tripd 	d; 
double	bond, dedr, pe/*,dr*/;

d.fx = pos[i].fx - pos[j].fx;
d.fy = pos[i].fy - pos[j].fy;
d.fz = pos[i].fz - pos[j].fz;
bond = sqrt(d.fx*d.fx+d.fy*d.fy+d.fz*d.fz);
//dr = bond - BeqBond[k];
//if(fabs(dr)>0.2) fprintf(stderr,"bond stretched!! %d--%d,(%d-%d) bond. dr = %f\n",(i-(natoms-BrOs*nBrO-nBCD*BCDs))%BrOs,(j-(natoms-BrOs*nBrO-nBCD*BCDs))%BrOs,i,j,dr); 

dedr = BkStr[k]*(bond-BeqBond[k])/bond;
pe = 0.5*BkStr[k]*sq(bond-BeqBond[k]);
force[i].fx -= (d.fx = dedr*d.fx);
force[i].fy -= (d.fy = dedr*d.fy);
force[i].fz -= (d.fz = dedr*d.fz);
force[j].fx += d.fx;
force[j].fy += d.fy;
force[j].fz += d.fz;
return(pe);

}
double BrOcalctorq (i1,i2,i3,i4,m)
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
      dudcos,n1n2, 
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
/*	for angles very close to 180 degrees, n1n2 = -1 + eps,
 *	where eps is very small and can be both positive and
 *	negative. If its negative acos will blow up. If its
 *	positive but very small, da will be very small for a linear
 *	equilibrium state and dedr will be the ratio between two
 *	very small quantities. We therefore use a switch	*/
//    n1n2 = K1/D1201;
//    if (n1n2 + 1 < 1.0e-8)  n1n2 = -1;
//    cosa = n1n2;
    cosa = K1/D1201;
//    index = acos(cosa)*36/PI;
//    gTors[m][index]++;
    /*
    **	    evaluate dU/dcos(a)
    **  U = pTors[m][0]*(1+cos(a))+pTors[m][1]*(1-cos(2*a))+pTors[m][2]*(1+cos(3*a))
    */
    pe = (1+cosa)*(BpTors[m][0]+2*BpTors[m][1]*(1-cosa)+BpTors[m][2]*sq(2*cosa-1));
    dudcos = BpTors[m][0] - 4*BpTors[m][1]*cosa + 3*BpTors[m][2]*(4*cosa*cosa-1);
    /*
    **	    add resulting force to atoms for angle
    */
    for (i=0;  i<4;  i++){
        force[ats[i]].fx += dudcos*grad[i][0];
	force[ats[i]].fy += dudcos*grad[i][1];
	force[ats[i]].fz += dudcos*grad[i][2];
    }
    return(pe);
}

double BrOcalcbend(i,j,k,m)
int i,j,k,m;
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
/*	for angles very close to 180 degrees, n1n2 = -1 + eps,
 *	where eps is very small and can be both positive and
 *	negative. If its negative acos will blow up. If its
 *	positive but very small, da will be very small for a linear
 *	equilibrium state and dedr will be the ratio between two
 *	very small quantities. We therefore use a switch	*/
if (n1n2+1 < 1.0e-8) n1n2 = -1;
da = acos(n1n2) - BeqBend[m];
pe = 0.5*BkBend[m]*da*da;
dedd = BkBend[m]*da/ (sqrt(1.-n1n2*n1n2)*r1*r2);
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

double calcnb(i,j,l,m)
int i,j,l,m;
{
double pe, ftmp, r2, ir2, ir6, a,b,aa,bb;
tripd rij;

rij.fx = pos[i].fx - pos[j].fx;
rij.fy = pos[i].fy - pos[j].fy;
rij.fz = pos[i].fz - pos[j].fz;

r2 = rij.fx*rij.fx + rij.fy*rij.fy + rij.fz*rij.fz;
ir2 = 1. / r2;
ir6 = ir2*ir2*ir2;
a = slj[3+l][3+m].a;
b = slj[3+l][3+m].b;
aa = 12*a;
bb = 6*b;
pe = factor14*((a*ir6-b) * ir6);
ftmp = factor14*(ir6*(aa*ir6 - bb)) * ir2; // -dV/dr
force[i].fx += (rij.fx *= ftmp);
force[i].fy += (rij.fy *= ftmp);
force[i].fz += (rij.fz *= ftmp);
force[j].fx -= rij.fx;
force[j].fy -= rij.fy;
force[j].fz -= rij.fz;
return(pe);
}

double BrOnb15(i,j,m,n)
int i,j,m,n;
{
double pe,ftmp,/*q2,*/r2,/*r,*/der,ir2,ir6,a,b,aa,bb/*,EC*/;
tripd rij;

rij.fx = pos[i].fx - pos[j].fx;
rij.fy = pos[i].fy - pos[j].fy;
rij.fz = pos[i].fz - pos[j].fz;

r2 = rij.fx*rij.fx + rij.fy*rij.fy + rij.fz*rij.fz;
ir2 = 1. / r2;
ir6 = ir2 * ir2 * ir2;
a = slj[3+m][3+n].a;
b = slj[3+m][3+n].b;
aa =12*a;
bb = 6*b;
//q2 = slj[3+m][3+n].q;
//r = sqrt(r2);
//EC = q2/r;
pe = (a*ir6-b)*ir6/* + EC*/;
ftmp = (ir6*(aa*ir6 - bb))*ir2/* + EC/r2*/; //-dV/dr
force[i].fx += (rij.fx *= ftmp);
force[i].fy += (rij.fy *= ftmp);
force[i].fz += (rij.fz *= ftmp);
force[j].fx -= rij.fx;
force[j].fy -= rij.fy;
force[j].fz -= rij.fz;
return(pe);
}
/*
void getBrOAng(int k, int l, int m, int j){
   tripd vec;
   double r,costh,zed;
   vec.fx = pos[k+l].fx - pos[k+m].fx;
   vec.fy = pos[k+l].fy - pos[k+m].fy;
   vec.fz = pos[k+l].fz - pos[k+m].fz;
   r = sqrt(vec.fx*vec.fx + vec.fy*vec.fy + vec.fz*vec.fz);
   costh = vec.fz/r;
   zed = (pos[k+l].fz + pos[k+m].fz)/2.0;
   if(fabs(zed)<sODrange){
      sOD[j][(int)((zed+sODrange)/sODbinz)][(int)((1.0+costh)/sODbin)]++;
      sODnorm[j][(int)((zed+sODrange)/sODbinz)]++;
      sOOP[j][(int)((zed+sODrange)/sODbinz)]+=(3.0*costh*costh - 1.0)/2.0;
   }
   sODt[j][(int)((1.0+costh)/sODbin)]++;
   sODtnorm[j]++;
}*/
