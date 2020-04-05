#include	<md.h>
#include	<system.h>
#include	 <math.h>
#include	<water.h>

waterHEX(pos,force)
	tripd	*pos;
	tripd	*force;
{
int	i, j, k, l, m, n,nw;
double	r2, dedr, delfx, delfy, delfz;
double		eg, s, sp;
tripd		frc, image, sdl, f[9];
double		wdljq();

nw = (natoms-nsolute-6*nHEX)/3;

for	(i = 0; i < nHEX; i++){
	k = nw*3+i*6;

	for	(j = 0; j < nw; j++,k = nw*3+i*6){
		l = j*3;
		eg = 0.0;

/*****	Determine image vector for C i - oxygen j.  ******/
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

		for	(m = 0; m < 9; m++)
			f[m].fx = f[m].fy = f[m].fz = 0.;


		for	(m = 0; m < 6; m++)/*Loop over atoms in hex molecules*/
			for	(n = 0; n < 3; n++)/*Loop over atoms in H2O*/
				{
/*****	Determine image vector ******/
				delfx = pos[k+m].fx - pos[l+n].fx + image.fx;
				delfy = pos[k+m].fy - pos[l+n].fy + image.fy;
				delfz = pos[k+m].fz - pos[l+n].fz + image.fz;
				r2 = delfx*delfx + delfy*delfy + delfz*delfz;

/*****	Get (1/r dV/dr) ******/
				dedr = wdljq(r2,m,n,&eg);

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
			for	(m = 0; m < 9; m++)
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
		k = nw*3+i*6;
		l = j*3;
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
		H2OHEXV += eg;
		}
	}
}

double
wdljq(r2,m,n,eg)
double r2;
int m; /* 0<=m<= 5 hexane atoms*/
int n; /* 0<=n<= 2 - Oxygen and Hydrogen.*/
double *eg;
{
double der,r,ir,ir6, a,b,aa,bb;
int k,bin;
	ir = 1. / r2;
	ir6 = ir * ir * ir;
	a = HEXlj[m][6+n].a;
	b = HEXlj[m][6+n].b;
	aa =12*a;
	bb = 6*b;
	*eg +=  ( a * ir6 - b ) * ir6;
	der = (bb - aa * ir6) * ir6 * ir;
	return(der);
}
