/*methanol-acetonitrile interactions*/
#include	<md.h>
#include	<system.h>
#include	 <math.h>

MAForce()
{
int	i, j, k, l, m, n;
double	r2, dedr, delfx, delfy, delfz;
double		eg, s, sp;
tripd		frc, image, sdl, f[6];
double		ljq12();

for	(i = 0; i < nCH3OH; i++){
	k = i*3;

	for	(j = 0; j < nCH3CN; j++,k = i*3){
		l = nCH3OH*3+j*3;
		eg = 0.0;

/*****	Determine image vector for O i - C j.  ******/
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

		for	(m = 0; m < 6; m++)
			f[m].fx = f[m].fy = f[m].fz = 0.;

/*****	Loop over atoms ******/
		for	(m = 0; m < 3; m++)
			for	(n = 0; n < 3; n++)
				{
/*****	Determine image vector ******/
				delfx = pos[k+m].fx - pos[l+n].fx + image.fx;
				delfy = pos[k+m].fy - pos[l+n].fy + image.fy;
				delfz = pos[k+m].fz - pos[l+n].fz + image.fz;
				r2 = delfx*delfx + delfy*delfy + delfz*delfz;
/***** Routine for determining H-bonds *******/
				if(r2 < OOdist2 && tc%HBdt == 0 && m == 0 && n == 1){
					getMAHB(k,l,delfx,delfy,delfz);
				}

/*****	Get (1/r dV/dr) ******/
				dedr = ljq12(r2,m,n,&eg);

/*****	Resolve forces on atoms.  ******/
				f[n+3].fx += (delfx *= dedr);
				f[n+3].fy += (delfy *= dedr);
				f[n+3].fz += (delfz *= dedr);
				f[m].fx -= delfx;
				f[m].fy -= delfy;
				f[m].fz -= delfz;
				}
		if	(sp != 0.0)
			{
			for	(m = 0; m < 6; m++)
				{
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
			}
		k = i*3;
		l = j*3+nCH3OH*3;
		force[k].fx += f[0].fx;
		force[k].fy += f[0].fy;
		force[k].fz += f[0].fz;
		force[++k].fx += f[1].fx;
		force[k].fy += f[1].fy;
		force[k].fz += f[1].fz;
		force[++k].fx += f[2].fx;
		force[k].fy += f[2].fy;
		force[k].fz += f[2].fz;

		force[l].fx += f[3].fx;
		force[l].fy += f[3].fy;
		force[l].fz += f[3].fz;
		force[++l].fx += f[4].fx;
		force[l].fy += f[4].fy;
		force[l].fz += f[4].fz;
		force[++l].fx += f[5].fx;
		force[l].fy += f[5].fy;
		force[l].fz += f[5].fz;
		VMA += eg;
		}
	}
}

double
ljq12(r2,m,n,eg)
double r2, *eg;
int m,n;
{
double der,r,ir,ir6, ec,a,b,q,aa,bb;

	ir = 1. / r2;
	ir6 = ir * ir * ir;
	a = liqlj[m][n+3].a;
	b = liqlj[m][n+3].b;
	q = liqlj[m][n+3].q;
	aa =12*a;
	bb = 6*b;
	r = sqrt(r2);
	ec = q/r;
	*eg +=  ( a * ir6 - b ) * ir6 +ec;
	der = (bb - aa * ir6) * ir6 * ir -ec/r2;
	return(der);
}

void getMAHB(int k, int l, double Ox, double Oy, double Oz)
{
	/* Oxyz is the N -> M(O) vector */
	tripd MOHvec;
	double rON, rOH, dot;
	/* get O->H vector */
	MOHvec.fx = pos[k+1].fx - pos[k].fx;
	MOHvec.fy = pos[k+1].fy - pos[k].fy;
	MOHvec.fz = pos[k+1].fz - pos[k].fz;

	rON = sqrt(Ox*Ox + Oy*Oy + Oz*Oz);
	rOH = sqrt(MOHvec.fx*MOHvec.fx + MOHvec.fy*MOHvec.fy + MOHvec.fz*MOHvec.fz);
	dot = Ox*MOHvec.fx + Oy*MOHvec.fy + Oz*MOHvec.fz;
	if(dot/(rON*rOH) > cosNOH)
	{
		/* don't do anything */
	}
}
		



