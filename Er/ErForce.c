#include	<md.h>
#include	<system.h>
#include	 <math.h>

ErForce()
{
    int	i, j, k, l, m;
    int mm;
    double totMass, r, r2, dedr, delfx, delfy, delfz;
    double eg, wnc, s, sp;
    int nw;
    tripd com, frc, image, sdl, f[3+1];
    double ewljq();
    double edljq();

    if(nEr!=1){
	fprintf(stderr,"ErForce.c -- ErForce was called but # of Er ions != 1\n.");
	exit(0);
    }

    nw = (natoms - nEr - nDDC*DDCsites)/3;

    for(i=0;i<nDDC;i++){
	k = i*DDCsites + nw*3;
	l = nDDC*DDCsites + nw*3; //index of Er3+ ion
	eg = wnc = 0.0;

	for(m=0;m<DDCsites;m++){ /*Loop over atoms in DDC molecules*/
	    for(mm=0;mm<4; mm++) //clear force array
		f[mm].fx = f[mm].fy = f[mm].fz = 0.;
	    /*****	Determine image vector, DDC site - Er3+.  ******/
	    sdl.fx = pos[k+m].fx - pos[l].fx;
	    sdl.fy = pos[k+m].fy - pos[l].fy;
	    sdl.fz = pos[k+m].fz - pos[l].fz;
	    mvimage(&sdl);
	    delfx = sdl.fx;
	    delfy = sdl.fy;
	    delfz = sdl.fz;
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
	    /*****	Get (1/r dV/dr) ******/
	    dedr = edljq(r2,m,&eg,&wnc);
	    /*****	Resolve forces on atoms.  ******/
	    f[1].fx += (delfx *= dedr);
	    f[1].fy += (delfy *= dedr);
	    f[1].fz += (delfz *= dedr);
	    f[0].fx -= delfx;
	    f[0].fy -= delfy;
	    f[0].fz -= delfz;
	}
	if(sp != 0.0){
	    for(m=0;m<4;m++){
		f[m].fx *= s;
		f[m].fy *= s;
		f[m].fz *= s;
	    }
	    f[0].fx -= (frc.fx = sp*eg*sdl.fx);
	    f[0].fy -= (frc.fy = sp*eg*sdl.fy);
	    f[0].fz -= (frc.fz = sp*eg*sdl.fz);
	    f[1].fx += frc.fx;
	    f[1].fy += frc.fy;
	    f[1].fz += frc.fz;
	    eg *= s;
	    wnc *= s;
	}
	k=nw*3+i*DDCsites+m;
	l = nDDC*DDCsites+nw*3;
	/* DDC */
	force[k].fx += f[0].fx;
	force[k].fy += f[0].fy;
	force[k].fz += f[0].fz;

	/* Er */
	force[l].fx += f[1].fx;
	force[l].fy += f[1].fy;
	force[l].fz += f[1].fz;
	
	V_ED += eg;
	ED_C += wnc;
    }

    for(i=0;i<nw;i++){
	k = i*3;
	l = nw*3 + nDDC*DDCsites; //index of Er3+ ion
	eg = wnc = 0.0;

	/*****	Determine image vector, water O parent - Er3+.  ******/
	image.fx = -(sdl.fx = pos[k].fx - pos[l].fx);
	image.fy = -(sdl.fy = pos[k].fy - pos[l].fy);
	image.fz = -(sdl.fz = pos[k].fz - pos[l].fz);
	mvimage(&sdl);
	image.fx += (delfx =sdl.fx);
	image.fy += (delfy =sdl.fy);
	image.fz += (delfz =sdl.fz);
	r2 = delfx*delfx + delfy*delfy + delfz*delfz;
	if(r2 >= swSoMax)
	    continue;
	if(r2 <= swSoMin)
	{
	    s = 1.0;
	    sp = 0.0;
	}
	else
	    swtch(r2-swSoMin,&s,&sp);
	for(m=0;m<4; m++) //clear full array (same for Er/water and Er/DDC) 
	    f[m].fx = f[m].fy = f[m].fz = 0.;
	for(m=0;m<3;m++){ /* Loop over atoms in water */
	    /*****	Determine image vector ******/
	    delfx = pos[k+m].fx - pos[l].fx + image.fx;
	    delfy = pos[k+m].fy - pos[l].fy + image.fy;
	    delfz = pos[k+m].fz - pos[l].fz + image.fz;
	    r2 = delfx*delfx + delfy*delfy + delfz*delfz;
	    /*****	Get (1/r dV/dr) ******/
	    dedr = ewljq(r2,m,&eg,&wnc);
	    /*****	Resolve forces on atoms.  ******/
	    f[3].fx += (delfx *= dedr);
	    f[3].fy += (delfy *= dedr);
	    f[3].fz += (delfz *= dedr);
	    f[m].fx -= delfx;
	    f[m].fy -= delfy;
	    f[m].fz -= delfz;
	}
	if(sp != 0.0){
	    for(m=0;m<4;m++){
		f[m].fx *= s;
		f[m].fy *= s;
		f[m].fz *= s;
	    }
	    f[0].fx -= (frc.fx = sp*eg*sdl.fx);
	    f[0].fy -= (frc.fy = sp*eg*sdl.fy);
	    f[0].fz -= (frc.fz = sp*eg*sdl.fz);
	    f[3].fx += frc.fx;
	    f[3].fy += frc.fy;
	    f[3].fz += frc.fz;
	    eg *= s;
	    wnc *= s;
	}
	k=i*3;
	l = nw*3+nDDC*DDCsites;
	/* H2O */
	force[k].fx += f[0].fx;
	force[k].fy += f[0].fy;
	force[k].fz += f[0].fz;
	force[++k].fx += f[1].fx;
	force[k].fy += f[1].fy;
	force[k].fz += f[1].fz;
	force[++k].fx += f[2].fx;
	force[k].fy += f[2].fy;
	force[k].fz += f[2].fz;

	/* Er3+ */
	force[l].fx += f[3].fx;
	force[l].fy += f[3].fy;
	force[l].fz += f[3].fz;
	
	V_EW += eg;
	EW_C += wnc;
    }


}

    double
edljq(r2,m,eg,wnc)
    double r2, *eg,*wnc;
    int m; /* m is the DDC atom index */
{
    double der,r,ir,ir6, ec,a,b,q,aa,bb;
    int bin;
    double en1;
    int index;
    ir = 1. / r2;
    ir6 = ir * ir * ir;
    a = EDlj[m].a;
    b = EDlj[m].b;
    q = EDlj[m].q;
    aa =12*a;
    bb = 6*b;
    r = sqrt(r2);
    ec = q/r;
    *wnc += ec;
    en1 =  ( a * ir6 - b ) * ir6 + ec;
    *eg +=  ( a * ir6 - b ) * ir6 + ec;
    der = (bb - aa * ir6) * ir6 * ir - ec/r2 ;

    return(der);
}
    double
ewljq(r2,m,eg,wnc) 
    double r2, *eg,*wnc;
    int m; /* m is the water atom index */
{
    double der,r,ir,ir6, ec,a,b,q,aa,bb;
    int w_index;
    int rindex;
    if(m == 0)
	w_index = 2;
    else if (m > 0) 
	w_index = 3;
    ir = 1. / r2;
    ir6 = ir * ir * ir;
    a = EWlj[m].a;
    b = EWlj[m].b;
    q = EWlj[m].q;
    aa = 12*a;
    bb = 6*b;
    r = sqrt(r2);
    rindex = r/binRDF;
    if(rindex<RDFbins){
	grSOL[w_index-2][rindex] += 1.0/(w_index-1);
	if(m==0)
	    if(r<rShellO_2)
		if(r<rShellO_1)
		    shellO_1++;
		else
		    shellO_2++;
    }
    ec = q/r;
    *wnc += ec;
    *eg +=  ( a * ir6 - b ) * ir6 + ec;
    der = (bb - aa * ir6) * ir6 * ir - ec/r2;
    return(der);
}
