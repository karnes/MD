#include	<md.h>
#include	<system.h>
#include	 <math.h>
#define EXCON 0.1


NITROForce()
{
int	i, j, k, l, m, n;
double	r2, dedr, delfx, delfy, delfz;
double		eg, nc, s, sp;
tripd		frc, image, sdl, f[28];
double		ttraljq();

bondE = bendE = torsE = impE = 0.;
for	(i = 0; i < nNIT; i++){
	k = natoms - nsolute-14*nNIT + i*14;
	intraNIT(k);
	if(pos[k].fz > 0.0 && pos[k].fz < NITzODdist){
		NIT_OD(k);
	}	
/*****	Get intermolecular forces for i - j NIT-NIT interactions.  ******/
	for	(j = i+1; j < nNIT; j++){
		k = natoms - nsolute-14*nNIT + i*14;

		l = natoms-nsolute-14*nNIT+j*14;
		eg = nc = 0.0;

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

		for	(m = 0; m < 28; m++)
			f[m].fx = f[m].fy = f[m].fz = 0.;

/*****	Loop over atoms in NIT molecules ******/
		for	(m = 0; m < 14; m++)
			for	(n = 0; n < 14; n++)
				{
/*****	Determine image vector ******/
				delfx = pos[k+m].fx - pos[l+n].fx + image.fx;
				delfy = pos[k+m].fy - pos[l+n].fy + image.fy;
				delfz = pos[k+m].fz - pos[l+n].fz + image.fz;
				r2 = delfx*delfx + delfy*delfy + delfz*delfz;
//fprintf(stderr," %6.2f \n ",sqrt(r2));
/*****	Get (1/r dV/dr) ******/
				dedr = ttraljq(r2,m,n,&eg,&nc);
/*****	Resolve forces on atoms.  ******/
				f[n+14].fx += (delfx *= dedr);
				f[n+14].fy += (delfy *= dedr);
				f[n+14].fz += (delfz *= dedr);
				f[m].fx -= delfx;
				f[m].fy -= delfy;
				f[m].fz -= delfz;
				}
		if	(sp != 0.0)
			{
			for	(m = 0; m < 28; m++)
				{
				f[m].fx *= s;
				f[m].fy *= s;
				f[m].fz *= s;
				}
/* make sure that f[0] and f[14] are the center of mass atoms */
			f[0].fx -= (frc.fx = sp*eg*sdl.fx);
			f[0].fy -= (frc.fy = sp*eg*sdl.fy);
			f[0].fz -= (frc.fz = sp*eg*sdl.fz);
			f[14].fx += frc.fx;
			f[14].fy += frc.fy;
			f[14].fz += frc.fz;
			eg *= s;
			nc *= s;
			}
		k = natoms - nsolute-14*nNIT + i*14;
		l = natoms-nsolute-14*nNIT+j*14;
/* This code has been modified for octanol  */
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
		force[++k].fx += f[6].fx;
		force[k].fy += f[6].fy;
		force[k].fz += f[6].fz;
		force[++k].fx += f[7].fx;
		force[k].fy += f[7].fy;
		force[k].fz += f[7].fz;
		force[++k].fx += f[8].fx;
		force[k].fy += f[8].fy;
		force[k].fz += f[8].fz;
		force[++k].fx += f[9].fx;
		force[k].fy += f[9].fy;
		force[k].fz += f[9].fz;
		force[++k].fx += f[10].fx;
		force[k].fy += f[10].fy;
		force[k].fz += f[10].fz;
		force[++k].fx += f[11].fx;
		force[k].fy += f[11].fy;
		force[k].fz += f[11].fz;
		force[++k].fx += f[12].fx;
		force[k].fy += f[12].fy;
		force[k].fz += f[12].fz;
		force[++k].fx += f[13].fx;
		force[k].fy += f[13].fy;
		force[k].fz += f[13].fz;
		
		force[l].fx += f[14].fx;
		force[l].fy += f[14].fy;
		force[l].fz += f[14].fz;
		force[++l].fx += f[15].fx;
		force[l].fy += f[15].fy;
		force[l].fz += f[15].fz;
		force[++l].fx += f[16].fx;
		force[l].fy += f[16].fy;
		force[l].fz += f[16].fz;
		force[++l].fx += f[17].fx;
		force[l].fy += f[17].fy;
		force[l].fz += f[17].fz;
		force[++l].fx += f[18].fx;
		force[l].fy += f[18].fy;
		force[l].fz += f[18].fz;
		force[++l].fx += f[19].fx;
		force[l].fy += f[19].fy;
		force[l].fz += f[19].fz;
		force[++l].fx += f[20].fx;
		force[l].fy += f[20].fy;
		force[l].fz += f[20].fz;
		force[++l].fx += f[21].fx;
		force[l].fy += f[21].fy;
		force[l].fz += f[21].fz;
		force[++l].fx += f[22].fx;
		force[l].fy += f[22].fy;
		force[l].fz += f[22].fz;
		force[++l].fx += f[23].fx;
		force[l].fy += f[23].fy;
		force[l].fz += f[23].fz;
		force[++l].fx += f[24].fx;
		force[l].fy += f[24].fy;
		force[l].fz += f[24].fz;
		force[++l].fx += f[25].fx;
		force[l].fy += f[25].fy;
		force[l].fz += f[25].fz;
		force[++l].fx += f[26].fx;
		force[l].fy += f[26].fy;
		force[l].fz += f[26].fz;
		force[++l].fx += f[27].fx;
		force[l].fy += f[27].fy;
		force[l].fz += f[27].fz;
		
		NITNB += eg;
		NITC += nc;
		}
	}

}

double
ttraljq(r2,m,n,eg,nc)
double r2, *eg, *nc;
int m,n;
{
double der,r,ir,ir6, ec,a,b,q,aa,bb,sqrt();

	ir = 1. / r2;
	ir6 = ir * ir * ir;
	a = NITlj[m][n].a;
	b = NITlj[m][n].b;

	q = NITlj[m][n].q;
	aa =12*a;
	bb = 6*b;

	r = sqrt(r2);
	ec = q/r;

	*eg +=  ( a * ir6 - b ) * ir6 + ec;
	*nc += ec;
	der = (bb - aa * ir6) * ir6 * ir - ec/r2;
	return(der);
}

intraNIT(k)
int k;
{
int i,n, n1, n2, n3, m;
double calcstrch(), calcbend(), calctorq();  
double bond_e, tors_e, bend_e, imp_e;
/* loop over all bond stretches*/
	bond_e = 0.;
/*fprintf(stderr,"%d ",k);*/
	/*carbon-carbon*/
	for (i=0;i<6;i++){
		n = k+i;
		m = k + (i+1)%6;
    		bond_e += calcstrch(n,m,0);
/*fprintf(stderr,"%f ",bond_e*KCAL);*/
	}
	/*carbon-hydrogen*/
	for (i=1;i<6;i++){
		n = k + i;
		m = k + i + 8;
    		bond_e += calcstrch(n,m,1);
/*fprintf(stderr,"%f ",bond_e*KCAL);*/
	}
	/*carbon-nitrogen*/
		n = k;
		m = k + 6;
    		bond_e += calcstrch(n,m,2);
/*fprintf(stderr,"%f ",bond_e*KCAL);*/
	/*oxygen-nitrogen*/
		n = k + 6;
		m = k + 7;
    		bond_e += calcstrch(n,m,3);
/*fprintf(stderr,"%f ",bond_e*KCAL);*/
		m = k + 8;
    		bond_e += calcstrch(n,m,3);
/*fprintf(stderr,"%f ",bond_e*KCAL);*/
/*fprintf(stderr,"\n");*/
/* loop over all torsions*/
	tors_e=0.;
	/*H-C-C-H*/
	for (i=1;i<5;i++){
		n = k + i + 8;/*H*/
		n1 = k + i;/*C*/
		n2 = n1 + 1;/*C*/
		n3 = n + 1;/*H*/
    		tors_e += calctorq(n,n1,n2,n3,0);
	}
	/*N-C-C-H*/
	n = k + 6;/*N*/
	n1 = k;/*C*/
	n2 = k + 1;/*C*/
	n3 = k + 9;/*H*/
    	tors_e += calctorq(n,n1,n2,n3,1);
	n2 = k + 5;/*C*/
	n3 = k + 13;/*H*/
    	tors_e += calctorq(n,n1,n2,n3,1);
	/*N-C-C-C*/
	n = k + 6;/*N*/
	n1 = k;/*C*/
	n2 = k + 1;/*C*/
	n3 = k + 2;/*C*/
    	tors_e += calctorq(n,n1,n2,n3,2);
	n2 = k + 5;/*C*/
	n3 = k + 4;/*C*/
    	tors_e += calctorq(n,n1,n2,n3,2);
	/*O-N-C-C*/
	n = k + 7;/*O*/
	n1 = k + 6;/*N*/
	n2 = k;/*C*/
	n3 = k + 1;/*C*/
    	tors_e += calctorq(n,n1,n2,n3,3);
	n3 = k + 5;/*C*/
    	tors_e += calctorq(n,n1,n2,n3,3);
	n = k + 8;/*O*/
	n3 = k + 1;/*C*/
    	tors_e += calctorq(n,n1,n2,n3,3);
	n3 = k + 5;/*C*/
    	tors_e += calctorq(n,n1,n2,n3,3);
/* loop over all angle bends*/
	bend_e = 0.;
	/*H-C-C*/
	for (i=1;i<6;i++){
		n = k+8+i;/*H*/
		n1 = k + i;/*C*/
		n2 = k + (i+1)%6;/*C*/
    		bend_e += calcbend(n,n1,n2,1);
		n2 = k + (i+5)%6;/*C*/
    		bend_e += calcbend(n,n1,n2,1);
	}
	/*C-C-C*/
	for (i=0;i<6;i++){
		n = k+i;
		n1 = k + (i+1)%6;
		n2 = k + (i+2)%6;
    		bend_e += calcbend(n,n1,n2,0);
	}
	/*N-C-C*/
	n = k+6;/*N*/
	n1 = k;/*C*/
	n2 = k + 1;/*C*/
    	bend_e += calcbend(n,n1,n2,2);
	n2 = k + 5;/*C*/
    	bend_e += calcbend(n,n1,n2,2);
	/*O-N-C*/
	n = k + 7;/*O*/
	n1 = k + 6;/*N*/
	n2 = k;/*C*/
    	bend_e += calcbend(n,n1,n2,3);
	n = k + 8;/*O*/
    	bend_e += calcbend(n,n1,n2,3);
	/*O-N-O*/
	n = k + 7;/*O*/
	n1 = k + 6;/*N*/
	n2 = k + 8;/*O*/
    	bend_e += calcbend(n,n1,n2,4);
/*improper torsions to keep the ring in a plan*/
	imp_e = 0.;
	imp_e += calctorq(k+1,k,k+3,k+4,4);
	imp_e += calctorq(k,k+5,k+2,k+3,4);
	imp_e += calctorq(k+2,k+1,k+4,k+5,4);
/*improper torsions to keep the H in the ring plan*/
	for (i=1;i<6;i++){
		n = k + i;/*i*/
		n1 = k + (i+1)%6;/*j*/
		n2 = k + (i+5)%6;/*k*/
		n3 = k + i + 8;/*l*/
    		imp_e += calctorq(n,n1,n2,n3,5);
	}
/*improper torsions to keep the N in the ring plan*/
	n = k;
	n1 = k + 1;
	n2 = k + 5;
	n3 = k + 6;
    	imp_e += calctorq(n,n1,n2,n3,6);
/*improper torsions to keep the NO2 in the ring plan*/
	n = k + 6;
	n1 = k + 7;
	n2 = k + 8;
	n3 = k;
    	imp_e += calctorq(n,n1,n2,n3,7);
/* total sum for all molecules*/
	bondE += bond_e;
	bendE += bend_e;
	torsE += tors_e;
	impE += imp_e;
/* total intramolecular potential energy*/
    	NITV += bond_e + tors_e + bend_e + imp_e;
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
    /*IMPA[m] = cosa;*/
/*
    index = acos(cosa)*36/PI;
    gTors[l][index]++;
*/
    /*
    **	    evaluate dU/dcos(a)
    **  U = pTors[l][0]*(1+cos(a))+pTors[l][1]*(1-cos(2*a))+pTors[l][2]*(1+cos(3*a))
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

double calcbend(i,j,k,l)
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
da = acos(n1n2) - eqBend[l];
pe = 0.5*kBend[l]*da*da;
dedd = kBend[l]*da/ (sqrt(1.-n1n2*n1n2)*r1*r2);
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

void NIT_OD(int i)
{
	tripd d;
	double di,c;
	int inx;

//	s.fx = sd.fx/rOI;
//	s.fy = sd.fy/rOI;
//	s.fz = sd.fz/rOI;

	d.fx = 2*pos[i+3].fx-pos[i+7].fx-pos[i+8].fx;
	d.fy = 2*pos[i+3].fy-pos[i+7].fy-pos[i+8].fy;
	d.fz = 2*pos[i+3].fz-pos[i+7].fz-pos[i+8].fz;
	di=sqrt(d.fx*d.fx+d.fy*d.fy+d.fz*d.fz);
	
	c = d.fz /= di;

//	d.fx /= di;
//	d.fy /= di;
//	d.fz /= di;

//	c=d.fx*0.0+d.fy*0.0+d.fz*1.0;
	inx = (1+c)/0.02;
	NITzOD[(int)(pos[i].fz/NITzODbinw)][inx] += 1.0;
	NITzNorm[(int)(pos[i].fz/NITzODbinw)] += 1.0;
	return;
}
