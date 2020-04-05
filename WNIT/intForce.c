#include	<md.h>
#include	<system.h>
#include	<water.h>
#include	<math.h>
/* non-bonded forces between the system atoms and the H2O/nitrobenzene (NIT)
 * liquids */

intForce(pos, force)
    tripd	*pos;
    tripd	*force;
{
    int	i, k, nw;
    double c, pe, ElecPot();

    if (nsolute == 0) return;
    k = natoms - nsolute;
    nw = k - 14*nNIT;
    VCH2O[0] = 0.;/* first shell water eletrostatic*/
    VCH2O[1] = 0.;/* all other water electrostatic*/
    VINTH2O = 0.;/*total interaction with all water*/
    VINTH2O_s = 0.;/*total interaction with 1st shell water*/
    V_teth = 0.0; /*total energy tethering fixShel waters */
    Nshel = 0; /*number of water molecules in first shell*/
    for	(i = 0; i < nw; i = i+3) {/* loop over water molecules*/
	pe = 0.;
	c = -2;
	h2oljc(&pos[i],&pos[k],&force[i],&force[k],&pe,&c,i/3);
	VINTH2O += pe;
	Cshel[i/3]=c;
	if(fixShel==3){
	    if(i/3==fixw[0] || i/3==fixw[1] || i/3==fixw[2]){
		teth(&pos[i],&pos[k],&force[i],&force[k]);
	    }
	}
    }

    VC_NIT = VINT_NIT = 0.;
    for	(i = 0; i < nNIT*14; i = i+14) {/* loop over NIT molecules*/
	pe = 0.;
	nitljc(&pos[i+nw],&pos[k],&force[i+nw],&force[k],&pe);
	VINT_NIT += pe;
    }
    VINT = VINT_NIT + VINTH2O + V_teth;
    if (EFieldOn) 
	VINT += ElecPot(pos,force);
    /***	Restrict the center of mass of the solute to a region.***/
    if(width_w > zwall) /* no window */
	return;
    if(width_w < 0.) /* if the size of the region is negative*/
	freez_cm(); /* freez center of mass	*/
    else /* use windowing potential centered at the closest integer
	  * and width width_w outside which there is a power law
	  * potential of power power_w and coeficients pot_w
	  */
	window();
    if(tc==0){
	if(fixShel==3){
	    if(fixed!=3){
		fprintf(stderr, "intForce.c: ion lost its fixed hydration shell. terminating run.\n");
		exit(1);
	    }
	}
    }

}
#define	EXBIAS	0.5
#define BIAS_W 11.0
#define BIAS_A 0.18
//#define BIAS_M 3.0
//#define BIAS_B 5.0
//#define BIAS_a 2.9
//#define BIAS_b -40.8
//#define BIAS_c 189.4
//#define BIAS_d -267
#define Gibbs_surf 7
#define BIAS_MOD 1.0 // multiply by this factor to scale the biasing force
window()
{
    double totMass,zeta,deriv,pow(),fabs(),zpos;
    tripd sysCm;
    int i;
    double bias_com = 0.0;
    static int printCOM = 1;
    V_w = 0;
    sysCm.fx = sysCm.fy = sysCm.fz = totMass = 0.;
    for(i = 0; i < natoms; i++) {
	sysCm.fx += mass[i] * pos[i].fx;
	sysCm.fy += mass[i] * pos[i].fy;
	sysCm.fz += mass[i] * pos[i].fz;
	totMass += mass[i];
    }
    sysCm.fx /= totMass;
    sysCm.fy /= totMass;
    sysCm.fz /= totMass;
    if (fabs(sysCm.fz - osysCmz) > 0.001 || printCOM){
	fprintf(stderr,"C.O.M = %f \n",sysCm.fz);
	printCOM = 0;
    }
    osysCmz = sysCm.fz;
    if(nsolute == 1){
	/*** Add biasing potential for ion. Form of potential is *** 
	 *** V_bias = BIAS_W (1 + tanh(BIAS_A*x) )/2             ***
	 ***   where x = pos[ion].fz - Gibbs_surf - bias_com  ***/
	if(biasingOn==1){
	    bias_com = osysCmz;
	    zpos = pos[natoms-1].fz - Gibbs_surf + bias_com;
	    //fprintf(stderr,"zpos = %f, posIon = %f, G_s = %f, bias_com = %f\n",zpos,pos[natoms-1].fz,Gibbs_surf,bias_com);
	    //  hyperbolic tangent form here
	    V_bias = - (BIAS_MOD * ( BIAS_W*(1+tanh(BIAS_A*zpos)) / 2.0 ) / KCAL );
	    deriv = - (BIAS_MOD * ( (BIAS_W*BIAS_A/2.0) / 
			( cosh(BIAS_A*zpos) * cosh(BIAS_A*zpos) ) ) / KCAL);

	    //  linear form here
	    //	    V_bias = - ((BIAS_M * zpos) + BIAS_B ) / KCAL;
	    //	    deriv = - BIAS_M / KCAL;

	    //  quadratic form here
	    /*	    V_bias = - (BIAS_a *(zpos - BIAS_b) * (zpos - BIAS_b) + BIAS_c ) / KCAL;
		    deriv = - 2 * BIAS_a * (zpos - BIAS_b) / KCAL;
		    */
	    /*  cubic form here
		V_bias = - (BIAS_a * zpos * zpos * zpos +  BIAS_b * zpos * zpos 
		+ BIAS_c * zpos + BIAS_d ) / KCAL;
		deriv = - (3.0 * BIAS_a * zpos * zpos + 2.0 *BIAS_b * zpos
		+ BIAS_c ) / KCAL;
		*/
	    force[natoms-1].fz -= deriv;
	    for(i=0;i<natoms;i++){
		force[i].fz+=deriv*mass[i]/totMass;
	    }
	    VINT+=V_bias;
	}
	/*** Add the window potential if solute outside the window ***/
	if(KillFinger==0){
	    zpos = pos[natoms-1].fz - sysCm.fz - center_w;
	}
	else{
	    zpos = pos[natoms-1].fz - 9.8 - center_w;
	}
	zeta = fabs(zpos)-width_w/2.;
	if(zeta > 0.){/* outside window, calculate energy */
	    if (zeta > 1.5*width_w){
		fprintf(stderr,"window: too much outside the window\n");
		exit(1);
	    }
	    else { /* 0 < zeta < width_w*/
		V_w = pot_w*pow(zeta,power_w-1);
		deriv = V_w*power_w*sgn(zpos);
		V_w *= zeta;
		force[natoms-1].fz -= deriv;
		if(KillFinger==0){
		    for (i=0;i<natoms;i++){
			force[i].fz += (deriv * mass[i]/totMass);
		    }
		    VINT += V_w;
		}
		else{
		    VINT += V_w;
		}
	    }
	}
    }
    else{
	if(nsolute == 2){
	    int solA, solB;
	    double zpos, muA, muB;
	    solA = natoms-2;
	    solB = natoms-1;
	    muA = mass[solA]/(mass[solA] + mass[solB]);
	    muB = 1. - muA;
	    zpos = muA * pos[solA].fz + muB * pos[solB].fz - center_w - sysCm.fz;
	    zeta = fabs(zpos)-width_w/2.;
	    if (zeta < 0.) /*inside window, no interactions*/
		return;
	    else {/*zeta > 0.; outside window, calculate energy */
		if (zeta > width_w){
		    fprintf(stderr,"window:too much outside the window\n");
		    exit(1);
		}
		else { /* 0 < zeta < width_w*/
		    V_w = pot_w*pow(zeta,power_w-1);
		    deriv = V_w*power_w*sgn(zpos);
		    V_w *= zeta;
		    force[solA].fz -= muA*deriv;
		    force[solB].fz -= muB*deriv;
		    for (i=0;i<natoms;i++)
			force[i].fz += (deriv * mass[i]/totMass);
		    VINT += V_w;
		}
	    }
	}
	else{
	    fprintf(stderr,"No window setup for polyatomic solute\n");
	    exit(1);
	}
    }
}

freez_cm()
{
    if (nsolute == 1){
	force[natoms-1].fz = 0.;
    }
    else	{
	if (nsolute == 2){
	    int solA, solB;
	    double zCMfrc, muA, muB, zCMvel;
	    solA = natoms-2;
	    solB = natoms -1;
	    muA = mass[solA]/(mass[solA] + mass[solB]);
	    muB = 1. - muA;
	    zCMfrc = force[solA].fz + force[solB].fz;
	    force[solA].fz -= muA*zCMfrc;
	    force[solB].fz -= muB*zCMfrc;
	}
    }
}

h2oljc(ri,rj,fi,fj,pe,c,ind)
    tripd *ri;/* ri[0], ri[1] and ri[2] are position of the O H1 and H2 of water*/
    tripd *rj;/* rj[0] is the position of the solute atom*/
    tripd *fi;/* fi[0], fi[1] and fi[2] are forces on the O H1 and H2 of water*/
    tripd *fj;/* force on solute atom*/
    double *pe;/* interaction energy between one water and the solute*/
    double *c;/* water dipole orienataion */
    int ind; /* index of water O (i/3) */
{
    double sqrt(), s, sp, f, a, aa, b, bb, q2, r2, r, ir, ir6, EC, rOI,shelMax;
    tripd image, rij, sdist;
    int i, j, rindex, w_index;
    double wsdip();
    void wsdip2();
    void dipOD();

    /***	Determine Oxygen-solute image vector			***/

    image.fx = -(sdist.fx = ri[0].fx - rj[0].fx);
    image.fy = -(sdist.fy = ri[0].fy - rj[0].fy);
    image.fz = -(sdist.fz = ri[0].fz - rj[0].fz);
    mvimage(&sdist);
    image.fx += sdist.fx;
    image.fy += sdist.fy;
    image.fz += sdist.fz;
    //if(ind==fixw[0] || ind==fixw[1] || ind==fixw[2]){
    //fprintf(stderr, "befr: fx = %f, fy = %f, fz = %f\n",sdist.fx,sdist.fy,sdist.fz);}
    r2 = sdist.fx*sdist.fx +  sdist.fy*sdist.fy +  sdist.fz*sdist.fz;
    rOI = sqrt(r2);
    //put here to avoid editing water routines in ../water
    if (ri[0].fz < 0.0 && ri[0].fz > -zODdist){
	dipOD(ri);
    }
    if ( r2 >= swSoMax )
	return;
    if ( r2 <= swSoMin ){
	s = 1.0;
	sp = 0.0;
    }
    else
	swtchSolu(r2 - swSoMin, &s, &sp);
    //rOI = sqrt(r2);
    if (rOI < Rshel){
	*c = wsdip(rOI,sdist,ri);
	Nshel++;/* number of first shell water*/
    }
    //fixShel: assign indices at tc = 0
    /*if(tc==0){
      if(action[1]=='E')
      shelMax = Rshel + 0.10;
      else
      shelMax = Rshel - 0.2;
      if(rOI < shelMax){
      if(fixed < fixShel){
      fixw[fixed] = ind;
      fixed++;
      fprintf(stderr,"liqForce.c - tether water %d to ion. (%d/%d)\n",ind,fixed,fixShel);
      }
      }
      }*/
    if (rOI < rODrad){
	wsdip2(rOI,sdist,ri);
    }

    /***	Loop over atoms in water molecule			***/

    for (i=0; i< 3; i++){
	if (i == 0) w_index = 2;
	if (i > 0) w_index = 3;
	rij.fx = ri[i].fx - rj[0].fx + image.fx;
	rij.fy = ri[i].fy - rj[0].fy + image.fy;
	rij.fz = ri[i].fz - rj[0].fz + image.fz;
	r2 = rij.fx*rij.fx +  rij.fy*rij.fy +  rij.fz*rij.fz ;
	a = NITlj[i+14][17].a;
	aa = 12.*a;
	b = NITlj[i+14][17].b;
	bb = 6.*b;
	q2 = NITlj[i+14][17].q;
	ir = 1. / r2;
	ir6 = ir * ir * ir;
	r = sqrt(r2);
	rindex = r/binRDF;
	if(rindex<300)
	    grSOL[w_index-2][rindex] += 1.0/(w_index-1);   
	if (neqFlag == 2) q2 = -q2;
	EC = q2/r;
	if (rOI <Rshel){/*this water is in the first hydration shell*/
	    VCH2O[0] += EC*s;
	    VINTH2O_s +=  (( a * ir6 - b ) * ir6 + EC)*s;
	} else VCH2O[1] += EC*s;/*electrostatic energy of all other water*/

	if (neqFlag == 1)
	    /* run with no charge on the ion, but note that VCH2O include *
	     * the coulomb interaction that would be if the charges were  *
	     * on but with the solvent dynamics without the charges       */
	    EC = 0.;

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
}

nitljc(ri,rj,fi,fj,pe)
    tripd *ri;/* ri[0], ri[2] are position of the two carbons of NIT   *
	       * ri[1], ri[3] are position of the two chlorines of NIT */
    tripd *rj;/* rj[0] is the position of the solute atom*/
    tripd *fi;/* same indeces as above for the forces		   */
    tripd *fj;/* force on solute atom*/
    double *pe;/* interaction energy between one NIT and the solute*/
{
    double sqrt(), s, sp, f, a, aa, b, bb, q2, r2, r, ir, ir6, EC;
    tripd image, rij, sdist;
    int i, rindex;

    /***	Determine first Carbon-solute image vector		***/

    image.fx = -(sdist.fx = ri[0].fx - rj[0].fx);
    image.fy = -(sdist.fy = ri[0].fy - rj[0].fy);
    image.fz = -(sdist.fz = ri[0].fz - rj[0].fz);
    mvimage(&sdist);
    image.fx += sdist.fx;
    image.fy += sdist.fy;
    image.fz += sdist.fz;
    r2 = sdist.fx*sdist.fx +  sdist.fy*sdist.fy +  sdist.fz*sdist.fz;
    if ( r2 >= swSoMax )
	return;
    if ( r2 <= swSoMin ){
	s = 1.0;
	sp = 0.0;
    }
    else
	swtchSolu(r2 - swSoMin, &s, &sp);

    /***	Loop over atoms in NIT molecule			***/

    for (i=0; i< 14; i++){
	rij.fx = ri[i].fx - rj[0].fx + image.fx;
	rij.fy = ri[i].fy - rj[0].fy + image.fy;
	rij.fz = ri[i].fz - rj[0].fz + image.fz;
	r2 = rij.fx*rij.fx +  rij.fy*rij.fy +  rij.fz*rij.fz ;
	a = NITlj[i][17].a;
	aa = 12.*a;
	b = NITlj[i][17].b;
	bb = 6.*b;
	q2 = NITlj[i][17].q;
	ir = 1. / r2;
	ir6 = ir * ir * ir;
	r = sqrt(r2);
	rindex = r/binRDF;
	if (rindex < 300){
	    if (i == 0) grSOL[2][rindex] += 1.0;   
	    if (i == 1 || i == 5) grSOL[3][rindex] += 0.5;   
	    if (i == 2 || i == 4) grSOL[4][rindex] += 0.5;   
	    if (i == 3 ) grSOL[5][rindex] += 1.0;   
	    if (i == 6 ) grSOL[6][rindex] += 1.0;   
	    if (i == 7 || i == 8 ) grSOL[7][rindex] += 0.5;   
	}
	if (neqFlag == 2) q2 = -q2;
	EC = q2/r;
	VC_NIT += EC*s;
	if (neqFlag == 1)
	    /* run with no charge on the ion, but note that VCNIT include *
	     * the coulomb interaction that would be if the charges were  *
	     * on but with the solvent dynamics without the charges       */
	    EC = 0.;

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
}
/*
 *	swtchSolu:
 *
 *	This subroutine generates the switching function and its derivative
 *	given the value of zz.  it is normalized to the size of the box.
 */
swtchSolu(zz,s,sp)
    double	zz, *s, *sp;
{
    double	zz2;
    zz2 = zz * zz;
    *s = swSocoef[0]+zz*zz2*(swSocoef[1]+zz*swSocoef[2]+zz2*swSocoef[3]);
    *sp = zz2 * (pswSocoef[0] + zz*pswSocoef[1] + zz2*pswSocoef[2]);
}

double wsdip(double rOI,tripd sd,tripd *ri)
{
    tripd s,d;
    double di,c;
    int inx;

    s.fx = sd.fx/rOI;
    s.fy = sd.fy/rOI;
    s.fz = sd.fz/rOI;

    d.fx = ri[1].fx+ri[2].fx-2*ri[0].fx;
    d.fy = ri[1].fy+ri[2].fy-2*ri[0].fy;
    d.fz = ri[1].fz+ri[2].fz-2*ri[0].fz;
    di=sqrt(d.fx*d.fx+d.fy*d.fy+d.fz*d.fz);

    d.fx /= di;
    d.fy /= di;
    d.fz /= di;

    c=d.fx*s.fx+d.fy*s.fy+d.fz*s.fz;
    inx = (1+c)/0.02;
    pwdip[inx] += 1.0;
    wdipNorm += 1.0;
    return(c);
}

void wsdip2(double rOI,tripd sd,tripd *ri)
{
    tripd s,d;
    double di,c;
    int inx;

    s.fx = sd.fx/rOI;
    s.fy = sd.fy/rOI;
    s.fz = sd.fz/rOI;

    d.fx = ri[1].fx+ri[2].fx-2*ri[0].fx;
    d.fy = ri[1].fy+ri[2].fy-2*ri[0].fy;
    d.fz = ri[1].fz+ri[2].fz-2*ri[0].fz;
    di=sqrt(d.fx*d.fx+d.fy*d.fy+d.fz*d.fz);

    d.fx /= di;
    d.fy /= di;
    d.fz /= di;

    c=d.fx*s.fx+d.fy*s.fy+d.fz*s.fz;
    inx = (1+c)/0.02;
    pwdip2[(int)(rOI/rODbinw)][inx] += 1.0;
    wdipNorm2[(int)(rOI/rODbinw)] += 1.0;
    return;
}

void dipOD(tripd *ri)
{
    tripd d;
    double di,c;
    int inx;

    //	s.fx = sd.fx/rOI;
    //	s.fy = sd.fy/rOI;
    //	s.fz = sd.fz/rOI;

    d.fx = ri[1].fx+ri[2].fx-2*ri[0].fx;
    d.fy = ri[1].fy+ri[2].fy-2*ri[0].fy;
    d.fz = ri[1].fz+ri[2].fz-2*ri[0].fz;
    di=sqrt(d.fx*d.fx+d.fy*d.fy+d.fz*d.fz);

    c = d.fz /= di;

    //	d.fx /= di;
    //	d.fy /= di;
    //	c = d.fz /= di;

    //	c=d.fx*0.0+d.fy*0.0+d.fz*1.0;
    inx = (1+c)/0.02;
    zODdip[(int)(-ri[0].fz/zODbinw)][inx] += 1.0;
    zdipNorm[(int)(-ri[0].fz/zODbinw)] += 1.0;
    return;
}

/* tethering forces */
teth(ri,rj,fi,fj)
    tripd *ri;/* ri[0], ri[1] and ri[2] are position of the O H1 and H2 of water*/
    tripd *rj;/* rj[0] is the position of the solute atom*/
    tripd *fi;/* fi[0], fi[1] and fi[2] are forces on the O H1 and H2 of water*/
    tripd *fj;/* force on solute atom*/
{
    double sqrt(),r2,rOI,deriv;
    double dx,dy,dz;
    tripd image, sdist;

    /***	Determine Oxygen-solute image vector			***/

    image.fx = -(sdist.fx = ri[0].fx - rj[0].fx);
    image.fy = -(sdist.fy = ri[0].fy - rj[0].fy);
    image.fz = -(sdist.fz = ri[0].fz - rj[0].fz);
    mvimage(&sdist);
    image.fx += sdist.fx;
    image.fy += sdist.fy;
    image.fz += sdist.fz;
    r2 = sdist.fx*sdist.fx +  sdist.fy*sdist.fy +  sdist.fz*sdist.fz;
    rOI = sqrt(r2);

    //fixShel tethering forces
#define TethCon 1.0 
    if(rOI > Rshel){
	V_teth += TethCon*(rOI-Rshel)*(rOI-Rshel);
	deriv = -2.0*TethCon*(rOI-Rshel);
	dx=(ri[0].fx-rj[0].fx+image.fx)/rOI/*(sdist.fx/rOI)*/;
	dy=(ri[0].fy-rj[0].fy+image.fy)/rOI/*(sdist.fy/rOI)*/;
	dz=(ri[0].fz-rj[0].fz+image.fz)/rOI/*(sdist.fz/rOI)*/;
	fi[0].fx += deriv*dx;
	fi[0].fy += deriv*dy;
	fi[0].fz += deriv*dz;
	fj[0].fx -= deriv*dx;
	fj[0].fy -= deriv*dy;
	fj[0].fz -= deriv*dz;
    }

}
/*
   void getTeth(int nw){

   tripd image, sdist;
   int i,j;
   double rOI, r2;
   double watr[nw/3][2]={0.0};

   for(i=0;i<nw;i+=3){
   image.fx = -(sdist.fx = pos[i].fx - pos[natoms-1].fx);
   image.fy = -(sdist.fy = pos[i].fy - pos[natoms-1].fy);
   image.fz = -(sdist.fz = pos[i].fz - pos[natoms-1].fz);
   mvimage(&sdist);
   image.fx += sdist.fx;
   image.fy += sdist.fy;
   image.fz += sdist.fz;
   r2 = sdist.fx*sdist.fx +  sdist.fy*sdist.fy +  sdist.fz*sdist.fz;
   rOI = sqrt(r2);
   watr[i][0]=(float)i/3;
   watr[i][1]=rOI;
   }

   quicksort(watr,0,nw/3);

   for(i=0;i<fixShel;i++){
   fixw[i]=(int)(watr[i][0]);
   fixed++;
   fprintf(stderr,"liqForce.c - tether water %d to ion. (%d/%d)\n",fixw[i],fixed,fixShel);
   }

   }*/
