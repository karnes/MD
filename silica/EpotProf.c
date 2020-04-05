#include	<md.h>
#include	<system.h>

/* calculate the average and fluctuations in the electrostatic potential in slabs parallel to the surface*/

EpotProf()
{
int	i, j, k, l, jx, jy, jz;
double vc, VC, stepTot[NZ];
       
/*CH3OH contribution*/
l = 0;
for	(jz = 0; jz < NZ; jz++) {
	VC = 0.;
	for (jx = 0; jx < NX; jx++) {
	    for (jy = 0; jy < NY; jy++) {
		j = NX*NY*jz+NY*jx+jy;
		for (i =0;i<nCH3OH*3;i=i+3) {
			vc = 0.;
			Epot(&pos[i],&potGrid[j],&vc,l);
//if(jx==0 && jy==0 && jz==3) fprintf(stderr,"MeOH\t%f\n",vc);
			VC += vc;
		}
	    }
	}
	VC = VC/(NX*NY);
	pGrid[0][jz] += VC;
	pGrid2[0][jz] += VC*VC;
	stepTot[jz] = VC;
}
/*CH3CN contribution*/
l = 1;
for	(jz = 0; jz < NZ; jz++) {
	VC = 0.;
	for (jx = 0; jx < NX; jx++) {
	    for (jy = 0; jy < NY; jy++) {
		j = NX*NY*jz+NY*jx+jy;
		for (i =0;i<nCH3CN*3;i=i+3) {
			k = nCH3OH*3+i;
			vc = 0.;
			Epot(&pos[k],&potGrid[j],&vc,l);
			VC += vc;
		}
	    }
	}
	VC = VC/(NX*NY);
	pGrid[1][jz] += VC;
	pGrid2[1][jz] += VC*VC;
	stepTot[jz] += VC;
}
/*Silica contribution*/
l = 2;
for	(jz = 0; jz < NZ; jz++) {
	VC = 0.;
	for (jx = 0; jx < NX; jx++) {
	    for (jy = 0; jy < NY; jy++) {
		j = NX*NY*jz+NY*jx+jy;
		for (i =0;i<nSi*3;i=i+3) {
			k = nCH3OH*3 + nCH3CN*3 +i;
			vc = 0.;
//			if( (i/3) % SiSkip != 0 || SiSkip > nSi){
				Epot(&pos[k],&potGrid[j],&vc,l);
//			}
//if(jx==0 && jy==0 && jz==3) fprintf(stderr,"SiOH\t%f\n",vc);
			VC += vc;
		}
	    }
	}
	VC = VC/(NX*NY);
	pGrid[2][jz] += VC;
	pGrid2[2][jz] += VC*VC;
	stepTot[jz] += VC;
}

for(jz=0; jz < NZ; jz++)
{
	pGrid[3][jz] += stepTot[jz];
	pGrid2[3][jz] += stepTot[jz]  * stepTot[jz];
}

}

Epot(ri,rj,vc,l)
tripd *ri;/* ri[0], ri[1] and ri[2] are position of the liquids or SiOH atoms*/
tripd *rj;/* rj[0] is the position of the grid point*/
double *vc;/* coulomb interaction energy with one liquid or SiOH molecule*/
int l;/*l = 0 for CH3OH, l = 1 for CH3CN and l = 2 for SiOH*/
{
double sqrt(), s, sp, r2, EC;
tripd image, rij, sdist;
int i;

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
	EC = gridQ2[i+l*3]/sqrt(r2);
	*vc += EC*s;
}
}
