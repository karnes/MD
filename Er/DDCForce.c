#include	<md.h>
#include	<system.h>
#include	 <math.h>

DDCForce()
{
    int	i, j, k, l, m, n;
    double	r2, dedr;

    tripd del;
    double		eg, s, sp;
    double		DDljq();
    bondE = torsE = bendE = nbE = 0.0;
    for	(i = 0; i < nDDC; i++){
	k = natoms-nEr-nDDC*DDCsites+i*DDCsites;

	/*****	Get intramolecular forces for i'th CHAIN.  ******/
	intraDDC(k); 
	/*****	Get intermolecular forces for i - j CHAIN-CHAIN interactions.  ******/
	for	(j = i+1; j < nDDC; j++){
	    l = natoms-nEr-nDDC*DDCsites+j*DDCsites;
	    /*****	Loop over atoms in CHA molecules ******/
	    for	(m = 0; m < DDCsites; m++)
		for	(n = 0; n < DDCsites; n++){
		    eg = 0.0;
		    /*****	Determine image vector ******/
		    del.fx = pos[k+m].fx - pos[l+n].fx;
		    del.fy = pos[k+m].fy - pos[l+n].fy;
		    del.fz = pos[k+m].fz - pos[l+n].fz;
		    mvimage(&del);
		    r2 = del.fx*del.fx + del.fy*del.fy + del.fz*del.fz;
		    if	(r2 >= swr2max)
			continue;
		    if	(r2 <= swr2min)
		    {
			s = 1.0;
			sp = 0.0;
		    }
		    else
			swtch(r2-swr2min,&s,&sp);

		    /*****	Get (1/r dV/dr) ******/
		    dedr = DDljq(r2,m,n,&eg);

		    /*****	Resolve forces on atoms.  ******/
		    force[l+n].fx += (del.fx *= (dedr*s+sp*eg));
		    force[l+n].fy += (del.fy *= (dedr*s+sp*eg));
		    force[l+n].fz += (del.fz *= (dedr*s+sp*eg));
		    force[k+m].fx -= del.fx;
		    force[k+m].fy -= del.fy;
		    force[k+m].fz -= del.fz;
		    DDCNB += eg*s;
		}
	}
    }
}

    double
DDljq(r2,m,n,eg)
    double r2, *eg;
    int m,n;
{
    double der,r,ir,ir6, a,b,aa,bb;
    //int mm,nn;
    ir = 1. / r2;
    ir6 = ir * ir * ir;
    //	if (m == 0 || m == nCarbon-1) mm = 0; else mm = 1;
    //	if (n == 0 || n == nCarbon-1) nn = 0; else nn = 1;
    a = DDlj[m][n].a;
    b = DDlj[m][n].b;
    aa =12*a;
    bb = 6*b;
    *eg +=  ( a * ir6 - b ) * ir6;
    der = (bb - aa * ir6) * ir6 * ir;
    return(der);
}

/***	Calculate the contributions of stretch, bend, torsion and 1-4
 ***	van der Waals to forces and energy for a 6 atoms chain molecule.
 ***/	

intraDDC(fa)
    int fa;
{
    int i,j,n;
    double calcstrch(), calcbend(), calctorq(), calcnb();  
    double bond_e, tors_e, bend_e, nb_e;
    /* loop over all bond stretches*/
    bond_e = calcstrch(fa,fa+1,0);/*CH3-CH2 bond*/
    for (i=1;i<DDCsites-2;i++){
	n = fa+i;
	bond_e += calcstrch(n,n+1,1);/*CH2-CH2 bonds*/
    }
    bond_e += calcstrch(fa+DDCsites-2,fa+DDCsites-1,0);/*CH3-CH2 bond*/
    /* loop over all torsions*/
    tors_e=0.;
    for (i=0;i<DDCsites-3;i++){
	n = fa+i;
	tors_e += calctorq(n,n+1,n+2,n+3,i);
    }
    /* loop over all angle bends*/
    bend_e = 0.;
    for (i=0;i<DDCsites-2;i++){
	n = fa+i;
	bend_e += calcbend(n,n+1,n+2);
    }
    /* loop over 1-5+ non-bonded energies*/
    nb_e = 0.;
    for (i=0;i<DDCsites-/*3*/4;i++){
	for (j=i+/*3*/4;j<DDCsites;j++){
	    nb_e += calcnb(fa+i,fa+j);
	}
    }
    /* total sum for all molecules*/
    bondE += bond_e;
    bendE += bend_e;
    torsE += tors_e;
    nbE += nb_e;
    /* total intramolecular potential energy*/
    INTRADDC += bond_e + tors_e + bend_e + nb_e;
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
    dedr = DDCkstr*(bond-DDCreq)/bond;
    pe = 0.5*DDCkstr*sq(bond-DDCreq);
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
    int i,j,k,ats[4];
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
    //    index = acos(cosa)*36/PI;
    //    gTors[m][index]++;
    /*
     **	    evaluate dU/dcos(a)
     **  U = pTors[0]*(1+cos(a))+pTors[1]*(1-cos(2*a))+pTors[2]*(1+cos(3*a))
     */
    pe = (1+cosa)*(pTors[0]+2*pTors[1]*(1-cosa)+pTors[2]*sq(2*cosa-1));
    dudcos = pTors[0] - 4*pTors[1]*cosa + 3*pTors[2]*(4*cosa*cosa-1);
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

double calcbend(i,j,k)
    int i,j,k;
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
    da = acos(n1n2) - DDCtheq;
    pe = 0.5*DDCkbend*da*da;
    dedd = DDCkbend*da/ (sqrt(1.-n1n2*n1n2)*r1*r2);
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

double calcnb(i,j)
    int i,j;
{
    double pe, ftmp, r2, r, ir2, ir6, EC,a,b,q,aa,bb;
    tripd rij;

    rij.fx = pos[i].fx - pos[j].fx;
    rij.fy = pos[i].fy - pos[j].fy;
    rij.fz = pos[i].fz - pos[j].fz;

    r2 = rij.fx*rij.fx + rij.fy*rij.fy + rij.fz*rij.fz;
    ir2 = 1. / r2;
    ir6 = ir2*ir2*ir2;
    a = DDlj[2][2].a;
    b = DDlj[2][2].b;
    aa = 12*a;
    bb = 6*b;
    pe = /*factor14**/((a*ir6 - b) * ir6);
    ftmp = /*factor14**/(ir6*(aa*ir6 - bb)) * ir2;
    force[i].fx += (rij.fx *= ftmp);
    force[i].fy += (rij.fy *= ftmp);
    force[i].fz += (rij.fz *= ftmp);
    force[j].fx -= rij.fx;
    force[j].fy -= rij.fy;
    force[j].fz -= rij.fz;
    return(pe);
}
