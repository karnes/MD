#include	<md.h>
#include	<system.h>
#include	 <math.h>

HEXForce()
{
int	i, j, k, l, m, n;
double	r2, dedr, delfx, delfy, delfz;
double		eg, s, sp;
tripd		frc, image, sdl, f[12];
double		ttraljq();
bondE = torsE = bendE = nb14E = 0.0;
for	(i = 0; i < nHEX; i++){
	k = natoms-nsolute-6*nHEX+i*6;

/*****	Get intramolecular forces for i'th HEX.  ******/
	intraHEX(k);
/*****	Get intermolecular forces for i - j HEX-HEX interactions.  ******/
	for	(j = i+1; j < nHEX; j++,k = natoms-nsolute-6*nHEX+i*6){
		l = natoms-nsolute-6*nHEX+j*6;
		eg = 0.0;

/*****	Determine image vector for primeC i - primeC j.  ******/
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

		for	(m = 0; m < 12; m++)
			f[m].fx = f[m].fy = f[m].fz = 0.;

/*****	Loop over atoms in HEX molecules ******/
		for	(m = 0; m < 6; m++)
			for	(n = 0; n < 6; n++)
				{
/*****	Determine image vector ******/
				delfx = pos[k+m].fx - pos[l+n].fx + image.fx;
				delfy = pos[k+m].fy - pos[l+n].fy + image.fy;
				delfz = pos[k+m].fz - pos[l+n].fz + image.fz;
				r2 = delfx*delfx + delfy*delfy + delfz*delfz;

/*****	Get (1/r dV/dr) ******/
				dedr = ttraljq(r2,m,n,&eg);

/*****	Resolve forces on atoms.  ******/
				f[n+6].fx += (delfx *= dedr);
				f[n+6].fy += (delfy *= dedr);
				f[n+6].fz += (delfz *= dedr);
				f[m].fx -= delfx;
				f[m].fy -= delfy;
				f[m].fz -= delfz;
				}
		if	(sp != 0.0)
			{
			for	(m = 0; m < 12; m++)
				{
				f[m].fx *= s;
				f[m].fy *= s;
				f[m].fz *= s;
				}
			f[0].fx -= (frc.fx = sp*eg*sdl.fx);
			f[0].fy -= (frc.fy = sp*eg*sdl.fy);
			f[0].fz -= (frc.fz = sp*eg*sdl.fz);
			f[6].fx += frc.fx;
			f[6].fy += frc.fy;
			f[6].fz += frc.fz;
			eg *= s;
			}
		k = natoms-nsolute-6*nHEX+i*6;
		l = natoms-nsolute-6*nHEX+j*6;
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
		force[++k].fx += f[5].fx;
		force[k].fy += f[5].fy;
		force[k].fz += f[5].fz;

		force[l].fx += f[6].fx;
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
		force[++l].fx += f[10].fx;
		force[l].fy += f[10].fy;
		force[l].fz += f[10].fz;
		force[++l].fx += f[11].fx;
		force[l].fy += f[11].fy;
		force[l].fz += f[11].fz;
		HEXNB += eg;
		}
	}
}

double
ttraljq(r2,m,n,eg)
double r2, *eg;
int m,n;
{
double der,ir,ir6, a,b,aa,bb;

	ir = 1. / r2;
	ir6 = ir * ir * ir;
	a = HEXlj[m][n].a;
	b = HEXlj[m][n].b;
	aa =12*a;
	bb = 6*b;
	*eg +=  ( a * ir6 - b ) * ir6;
	der = (bb - aa * ir6) * ir6 * ir;
	return(der);
}

/***	Calculate the contributions of stretch, bend, torsion and 1-4
 ***	van der Waals to forces and energy for a 6 atoms chain molecule.
 ***/	


intraHEX(fa)
int fa;
{
double calcstrch(), calcbend(), calctorq(), calcnb();  
double bond_e, tors_e, bend_e, nb14_e;
/* connectivity table:
  C1------C2------C3------C4------C5------C6
 fa+4    fa+2     fa      fa+1   fa+3    fa+5
*/
/* loop over all bond stretches*/
    	bond_e = calcstrch(fa+4,fa+2,0);
    	bond_e += calcstrch(fa+2,fa,1);
    	bond_e += calcstrch(fa,fa+1,2);
    	bond_e += calcstrch(fa+1,fa+3,3);
    	bond_e += calcstrch(fa+3,fa+5,4);
/* loop over all torsions*/
    	tors_e = calctorq(fa+4,fa+2,fa,fa+1,0);
    	tors_e += calctorq(fa+2,fa,fa+1,fa+3,1);
    	tors_e += calctorq(fa,fa+1,fa+3,fa+5,2);
/* loop over all angle bends*/
    	bend_e = calcbend(fa+4,fa+2,fa,0);
    	bend_e += calcbend(fa+2,fa,fa+1,1);
    	bend_e += calcbend(fa,fa+1,fa+3,2);
    	bend_e += calcbend(fa+1,fa+3,fa+5,3);
/* loop over 1-4 non-bonded energies*/
    	nb14_e = calcnb(fa+4,fa+1,0);
    	nb14_e += calcnb(fa+2,fa+3,1);
    	nb14_e += calcnb(fa,fa+5,2);
/* total sum for all molecules*/
	bondE += bond_e;
	bendE += bend_e;
	torsE += tors_e;
	nb14E += nb14_e;
/* total intramolecular potential energy*/
    	HEXV += bond_e + tors_e + bend_e + nb14_e;
}

double calcstrch(i,j,k)
int i,j,k;
{
tripd 	d; 
double	bond, dedr, pe;

d.fx = pos[i].fx - pos[j].fx;
d.fy = pos[i].fy - pos[j].fy;
d.fz = pos[i].fz - pos[j].fz;
bond = sqrt(d.fx*d.fx+d.fy*d.fy+d.fz*d.fz);
dedr = kStr[k]*(bond-eqBond[k])/bond;
pe = 0.5*kStr[k]*sq(bond-eqBond[k]);
force[i].fx -= (d.fx = dedr*d.fx);
force[i].fy -= (d.fy = dedr*d.fy);
force[i].fz -= (d.fz = dedr*d.fz);
force[j].fx += d.fx;
force[j].fy += d.fy;
force[j].fz += d.fz;
return(pe);

}
double calctorq (i1,i2,i3,i4,m)
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
    index = acos(cosa)*36/PI;
    gTors[m][index]++;
    /*
    **	    evaluate dU/dcos(a)
    **  U = pTors[m][0]*(1+cos(a))+pTors[m][1]*(1-cos(2*a))+pTors[m][2]*(1+cos(3*a))
    */
    pe = (1+cosa)*(pTors[m][0]+2*pTors[m][1]*(1-cosa)+pTors[m][2]*sq(2*cosa-1));
    dudcos = pTors[m][0] - 4*pTors[m][1]*cosa + 3*pTors[m][2]*(4*cosa*cosa-1);
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

double calcbend(i,j,k,m)
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
da = acos(n1n2) - eqBend[m];
pe = 0.5*kBend[m]*da*da;
dedd = kBend[m]*da/ (sqrt(1.-n1n2*n1n2)*r1*r2);
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

double calcnb(i,j,m)
int i,j,m;
{
double pe, ftmp, r2, ir2, ir6, a,b,aa,bb;
tripd rij;

rij.fx = pos[i].fx - pos[j].fx;
rij.fy = pos[i].fy - pos[j].fy;
rij.fz = pos[i].fz - pos[j].fz;

r2 = rij.fx*rij.fx + rij.fy*rij.fy + rij.fz*rij.fz;
ir2 = 1. / r2;
ir6 = ir2*ir2*ir2;
a = HEXlj[2*m+1][4-2*m].a;
b = HEXlj[2*m+1][4-2*m].b;
aa = 12*a;
bb = 6*b;
pe = factor14*((a*ir6 - b) * ir6);
ftmp = factor14*(ir6*(aa*ir6 - bb)) * ir2;
force[i].fx += (rij.fx *= ftmp);
force[i].fy += (rij.fy *= ftmp);
force[i].fz += (rij.fz *= ftmp);
force[j].fx -= rij.fx;
force[j].fy -= rij.fy;
force[j].fz -= rij.fz;
return(pe);
}
