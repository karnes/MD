#include	<md.h>
#include	<system.h>
#include	 <math.h>
#include	<water.h>

waterDDC(pos,force)
    tripd	*pos;
    tripd	*force;
{
    int	i, j, k, l, m, n,nw;
    int dd,mm;
    double	r2, dedr, delfx, delfy, delfz;
    double		eg, s, sp;
    tripd		frc, image, sdl, f[4];
    double		wdljq();
    
    nw = natoms-nEr-DDCsites*nDDC;
    for(i=0;i<nDDC;i++){
	k=nw+i*DDCsites;

	for(j=0;j<nw/3;j++){
	    l = j*3;
	    k = nw+i*DDCsites;

	    /*****	Determine image vector for C i - oxygen j.  ******/
	    for(m = 0; m < DDCsites; m++){ /*Loop over atoms in hex molecules*/
		dd = m + k;
	        l = j*3;
		image.fx = -(sdl.fx = pos[dd].fx - pos[l].fx);
		image.fy = -(sdl.fy = pos[dd].fy - pos[l].fy);
		image.fz = -(sdl.fz = pos[dd].fz - pos[l].fz);
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

	        eg = 0.0;
		for(mm = 0; mm < 4; mm++)
		    f[mm].fx = f[mm].fy = f[mm].fz = 0.;

		for(n = 0; n < 3; n++)/*Loop over atoms in H2O*/
		{
		    /*****	Determine image vector ******/
		    delfx = pos[dd].fx - pos[l+n].fx + image.fx;
		    delfy = pos[dd].fy - pos[l+n].fy + image.fy;
		    delfz = pos[dd].fz - pos[l+n].fz + image.fz;
		    r2 = delfx*delfx + delfy*delfy + delfz*delfz;

		    /*****	Get (1/r dV/dr) ******/
		    dedr = wdljq(r2,n,m,&eg);

		    /*****	Resolve forces on atoms.  ******/
		    f[n+1].fx += (delfx *= dedr);
		    f[n+1].fy += (delfy *= dedr);
		    f[n+1].fz += (delfz *= dedr);
		    f[0].fx -= delfx;
		    f[0].fy -= delfy;
		    f[0].fz -= delfz;
		}
		if	(sp != 0.0)
		{
		    for	(mm = 0; mm < 4; mm++)
		    {
			f[mm].fx *= s;
			f[mm].fy *= s;
			f[mm].fz *= s;
		    }
		    f[0].fx -= (frc.fx = sp*eg*sdl.fx);
		    f[0].fy -= (frc.fy = sp*eg*sdl.fy);
		    f[0].fz -= (frc.fz = sp*eg*sdl.fz);
		    f[1].fx += frc.fx;
		    f[1].fy += frc.fy;
		    f[1].fz += frc.fz;
		    eg *= s;
		}
		l = j*3;
		// DODECANE
		force[dd].fx += f[0].fx;
		force[dd].fy += f[0].fy;
		force[dd].fz += f[0].fz;

		//WATER
		force[l].fx += f[1].fx;
		force[l].fy += f[1].fy;
		force[l].fz += f[1].fz;
		force[++l].fx += f[2].fx;
		force[l].fy += f[2].fy;
		force[l].fz += f[2].fz;
		force[++l].fx += f[3].fx;
		force[l].fy += f[3].fy;
		force[l].fz += f[3].fz;
		V_WD += eg;
	    }
	}
    }
}

    double
wdljq(r2,n,m,eg)
    double r2;
    int m; /* 0<=m<= DDC atom*/
    int n; /* 0<=n<= 2 - water Oxygen and Hydrogen.*/
    double *eg;
{
    double der,r,ir,ir6, a,b,aa,bb;
    ir = 1. / r2;
    ir6 = ir * ir * ir;
    //water-DDC cross terms: [water][DDC]
    a = WDlj[n][m].a;
    b = WDlj[n][m].b;
    aa =12*a;
    bb = 6*b;
    *eg +=  ( a * ir6 - b ) * ir6;
    der = (bb - aa * ir6) * ir6 * ir;
    return(der);
}
