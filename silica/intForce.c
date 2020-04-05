#include	<md.h>
#include	<system.h>

/* calculate the non-bonded LJ forces between the system atoms the solvents and the silica*/

intForce(pos, force)
	tripd	*pos;
	tripd	*force;
{
int	i, j,k,l;
double pe,vc;
       
if	(nsolute == 0)
	return;
/*solute-CH3OH interactions*/
VINTSM = 0.;/*solute-mehanol interactions*/
l = 0;
for	(i =0;i<nCH3OH*3;i=i+3) {/* loop over methanol molecules*/
	for	(j = 0; j < nsolute; j++) {
		pe = vc = 0.;
		liqljc(&pos[i],&pos[natoms-nsolute+j],&force[i],&force[natoms-nsolute+j],&pe,&vc,j,l);
		VINTSM += pe;
		VINT += pe;
	}
}

/*solute-CH3CN interactions*/
VINTSA = 0.;/*solute-acetonitrile interactions*/
l = 1;
for	(i =0;i<nCH3CN*3;i=i+3) {/* loop over acetonitrile molecules*/
	k = nCH3OH*3+i;
	for	(j = 0; j < nsolute; j++) {
		pe = vc = 0.;
		liqljc(&pos[k],&pos[natoms-nsolute+j],&force[k],&force[natoms-nsolute+j],&pe,&vc,j,l);
		VINTSA += pe;
		VINT += pe;
	}
}
/*solute-silica interactions*/
VINTSS = 0.;/*solute-silica interactions*/
l = 2;
for	(i =0;i<nSi*3;i=i+3) {/* loop over SiOH surface units*/
	k = nCH3OH*3 + nCH3CN*3+i;
	for	(j = 0; j < nsolute; j++) {
		pe = vc = 0.;
		liqljc(&pos[k],&pos[natoms-nsolute+j],&force[k],&force[natoms-nsolute+j],&pe,&vc,j,l);
		VINTSS += pe;
		VINT += pe;
	}
}
	
}

liqljc(ri,rj,fi,fj,pe,vc,jsol,l)
tripd *ri;/* ri[0], ri[1] and ri[2] are position of the liquids or SiOH atoms*/
tripd *rj;/* rj[0] is the position of the solute atom*/
tripd *fi;/* fi[0], fi[1] and fi[2] are forces on the liquid or SiOH atoms*/
tripd *fj;/* force on solute atom*/
double *pe;/* interaction energy between one water and the solute*/
double *vc;/* coulomb interaction energy between one water and the solute*/
int jsol;/*solute atom in the pos array*/
int l;/*l = 0 for CH3OH, l = 1 for CH3CN and l = 2 for SiOH*/
{
double sqrt(), s, sp, f, a, aa, b, bb, q2, r2, r, ir, ir6, EC;
tripd image, rij, sdist;
int i,j;

/***	Determine image vector			***/

image.fx = -(sdist.fx = ri[0].fx - rj[0].fx);
image.fy = -(sdist.fy = ri[0].fy - rj[0].fy);
image.fz = -(sdist.fz = ri[0].fz - rj[0].fz);
pimage(&sdist);
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

/***	Loop over the three atoms in liquid or SiOH			***/

for (i=0; i< 3; i++){
	rij.fx = ri[i].fx - rj[0].fx + image.fx;
	rij.fy = ri[i].fy - rj[0].fy + image.fy;
	rij.fz = ri[i].fz - rj[0].fz + image.fz;
	r2 = rij.fx*rij.fx +  rij.fy*rij.fy +  rij.fz*rij.fz ;
	j = jsol;/*get atom type for lj parameters from the location in the array*/
	a = sollj[j][i+l*3].a;
	aa = 12.*a;
	b = sollj[j][i+l*3].b;
	bb = 6.*b;
	q2 = sollj[j][i+l*3].q;
	ir = 1. / r2;
	ir6 = ir * ir * ir;
	r = sqrt(r2);
	EC = q2/r;
	*vc += EC*s;
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
