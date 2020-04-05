#include	<md.h>
#include	<system.h>
#include	 <math.h>

CH3CNForce()
{
int	i, j, k, l, m, n;
double	r2,dedr,delfx,delfy,delfz;
double		eg, s, sp;
tripd		frc, image, sdl, f[6];
double		ljq2();
for	(i = 0; i < nCH3CN; i++)
	{
	k = nCH3OH*3+ i*3;

/*****	Get intramolecular forces for i'th CH3CN.  ******/
	intraCH3CN(&pos[k],&force[k]);

/****	Determine the cos of C-N vector (versus +z vector normal to silica)
	and store normalized vector for TCF  ***/
	if(tc % TCFdt == 0) CNangles(k,i);
/*****	Get intermolecular forces for i - j CH3CN-CH3CN interactions.  ******/
	for	(j = i+1; j < nCH3CN; j++, k = nCH3OH*3+i*3)
		{
		l = nCH3OH*3+j*3;
		eg = 0.0;

/*****	Determine image vector for Carbon i - Carbon j.  ******/
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

/*****	Loop over atoms in CH3CN molecules ******/
		for	(m = 0; m < 3; m++)
			for	(n = 0; n < 3; n++)
				{
/*****	Determine image vector ******/
				delfx = pos[k+m].fx - pos[l+n].fx + image.fx;
				delfy = pos[k+m].fy - pos[l+n].fy + image.fy;
				delfz = pos[k+m].fz - pos[l+n].fz + image.fz;
				r2 = delfx*delfx + delfy*delfy + delfz*delfz;

/*****	Get (1/r dV/dr) ******/
				dedr = ljq2(r2,m,n,&eg);

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
		k = nCH3OH*3+i*3;
		l = nCH3OH*3+j*3;
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
		VNB_A += eg;
		}
	}
}
/* intramolecular potential for linear triatomic*/
/* pos[0] = C, pos[1] = N, pos[2] = CH3*/
intraCH3CN(pos,force)
tripd	*pos,*force;
{
double		dedr[3],r[3];
tripdouble	r0,r1,r2,grad[2];
/*CN bond*/
r0.fx = pos[1].fx - pos[0].fx;
r0.fy = pos[1].fy - pos[0].fy;
r0.fz = pos[1].fz - pos[0].fz;
/*CC bond*/
r1.fx = pos[2].fx - pos[0].fx;
r1.fy = pos[2].fy - pos[0].fy;
r1.fz = pos[2].fz - pos[0].fz;

r[0] = sqrt(r0.fx*r0.fx + r0.fy*r0.fy + r0.fz*r0.fz);
r[1] = sqrt(r1.fx*r1.fx + r1.fy*r1.fy + r1.fz*r1.fz);
r[2] = r0.fx*r1.fx+r0.fy*r1.fy+r0.fz*r1.fz;

/* calculate molecular axis distribution*/
//AcOr[(int)((r0.fz/r[0]+1)/0.04)] += 1.0;
/* calculate molecular reorientation tcf*/
//ox[tc] = r0.fy/sqrt(r0.fx*r0.fx + r0.fy*r0.fy);/* tc is a global time integer*/
//oy[tc] = r0.fx/sqrt(r0.fx*r0.fx + r0.fy*r0.fy);
//getAcOr();

getharmfd(r,dedr);
force[1].fx -= dedr[0]*r0.fx;
force[1].fy -= dedr[0]*r0.fy;
force[1].fz -= dedr[0]*r0.fz;
force[2].fx -= dedr[1]*r1.fx;
force[2].fy -= dedr[1]*r1.fy;
force[2].fz -= dedr[1]*r1.fz;
force[0].fx += dedr[0]*r0.fx + dedr[1]*r1.fx;
force[0].fy += dedr[0]*r0.fy + dedr[1]*r1.fy;
force[0].fz += dedr[0]*r0.fz + dedr[1]*r1.fz;
force[1].fx += (grad[0].fx = dedr[2]*(r1.fx-r[2]*r0.fx/(r[0]*r[0])));
force[1].fy += (grad[0].fy = dedr[2]*(r1.fy-r[2]*r0.fy/(r[0]*r[0])));
force[1].fz += (grad[0].fz = dedr[2]*(r1.fz-r[2]*r0.fz/(r[0]*r[0])));
force[2].fx += (grad[1].fx = dedr[2]*(r0.fx-r[2]*r1.fx/(r[1]*r[1])));
force[2].fy += (grad[1].fy = dedr[2]*(r0.fy-r[2]*r1.fy/(r[1]*r[1])));
force[2].fz += (grad[1].fz = dedr[2]*(r0.fz-r[2]*r1.fz/(r[1]*r[1])));
force[0].fx -= grad[1].fx + grad[0].fx;
force[0].fy -= grad[1].fy + grad[0].fy;
force[0].fz -= grad[1].fz + grad[0].fz;
}

/*
 *	**** function to calculate a triatom harmonic potential		****
 *	**** with respect to the internal coordinates r.      		****
 *      **** r[0] and r[1] are the bond length and r[2]/(r[0]*r[1]) is 	****
 *	**** cos the bond angle						****
 */
#define eqCCN PI
#define eqCN 1.17
#define eqCC 1.46
#define kCCN (43.7/KCAL)
#define kCN  (2508.0/KCAL)
#define kCC (755.0/KCAL)
getharmfd(r,dedr)
double	*r,*dedr;
{
double n1n2, da, q1,q2;
	n1n2 = r[2] /(r[0]*r[1]);
/*	for angles very close to 180 degrees, n1n2 = -1 + eps,
 *	where eps is very small and can be both positive and
 *	negative. If its negative acos will blow up. If its
 *	positive but very small, da will be very small for a linear
 *	equilibrium state and dedr will be the ratio between two
 *	very small quantities. We therefore use a switch	*/
	if (n1n2+1 < 1.0e-8) n1n2 = -1;
	da = acos(n1n2) - eqCCN;
	q1 = r[0]-eqCN;
	q2 = r[1]-eqCC;
	INTRAV_A += 0.5*(kCN*q1*q1+kCC*q2*q2+kCCN*da*da);
	dedr[0] = kCN*q1/r[0];
	dedr[1] = kCC*q2/r[1];
	if (n1n2+1 < 1.0e-8)
		dedr[2] = -kCCN/(r[0]*r[1]);
	else
		dedr[2] = kCCN*da/ (sqrt(1.-n1n2*n1n2)*r[0]*r[1]);
}
double
ljq2(r2,m,n,eg)
double r2, *eg;
int m,n;
{
double der,r,ir,ir6, ec,a,b,q,aa,bb;

	ir = 1. / r2;
	ir6 = ir * ir * ir;
	a = liqlj[m+3][n+3].a;
	b = liqlj[m+3][n+3].b;
	q = liqlj[m+3][n+3].q;
	aa =12*a;
	bb = 6*b;
	r = sqrt(r2);
	ec = q/r;
	*eg +=  ( a * ir6 - b ) * ir6 +ec;
	der = (bb - aa * ir6) * ir6 * ir -ec/r2;
	return(der);
}
//getAcOr()
//{
//int t;
//for (t=0;t<=tc;t++)
//	corrAC[t] = ox[tc]*ox[tc-t]+oy[tc]*oy[tc-t];
//}

/****	Determine the cos of C-N vector (versus +z vector normal to silica)
	and store normalized vector for TCF  ***/
void CNangles(int k, int i)
{
	double	r, r2, delfx, delfy, delfz;
	double cosCN; 
	delfx = pos[k+1].fx - pos[k].fx;
	delfy = pos[k+1].fy - pos[k].fy;
	delfz = pos[k+1].fz - pos[k].fz;
	r2 = delfx*delfx + delfy*delfy + delfz*delfz;
	r = sqrt(r2);
	/* normalized C-N vector and z-position*/
	ACNvec[i][(int)(tc/TCFdt)].fx = (float)(delfx/r);
	ACNvec[i][(int)(tc/TCFdt)].fy = (float)(delfy/r);
	ACNvec[i][(int)(tc/TCFdt)].fz = (float)(delfz/r);
	ACNvec[i][(int)(tc/TCFdt)].zpos = (float)(pos[k].fz);//jjk - plan to store and calc TCFs in wxtr.c
	cosCN = delfz/r;
	odCN[(int)(pos[k].fz/ZBinSize)][(int)((cosCN + 1.0) / CosBinSize)]++;
	int z;
	for(z=0; z < NumSBins; z++)
	{
		if(pos[k].fz > sbins[2 * z] && pos[k].fz < sbins[2 * z + 1])
		{
			sodCN[z][(int)((cosCN + 1.0) / CosBinSize)]++;
		}
	}
}
