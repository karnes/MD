#include	<md.h>
#include	<water.h>
#include	<system.h>

waterForceWNIT(pos, force, na)
tripd	*pos;
tripd	*force;
int	na;
{
register int	i, j, k, l, m, n;
register double	r2, dedr, delfx, delfy, delfz;
double		eg, wc, s, sp;
tripd		frc, image, sdl, f[6];
double		spc(), sqrt();
double		watts();
/*static double	*addr[3][3] = {	{ootbl, ohtbl, ohtbl},
				{ohtbl, hhtbl, hhtbl},
				{ohtbl, hhtbl, hhtbl}};
*/
na /= 3;
for	(i = 0; i < na; i++)
	{
	k = i*3;

/*****	Get intramolecular forces for i'th water.  ******/
	intrawater(&pos[k],&force[k]);

/*****	Get intermolecular forces for i - j water-water interactions.  ******/
	for	(j = i+1; j < na; j++, k = i*3)
		{
		l = j*3;
		eg = wc = 0.0;

/*****	Determine image vector for OXYGEN i - OXYGEN j.  ******/
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

/*****	Loop over atoms in water molecules ******/
		for	(m = 0; m < 3; m++)
			for	(n = 0; n < 3; n++)
				{
/*****	Determine image vector ******/
				delfx = pos[k+m].fx - pos[l+n].fx + image.fx;
				delfy = pos[k+m].fy - pos[l+n].fy + image.fy;
				delfz = pos[k+m].fz - pos[l+n].fz + image.fz;
				r2 = delfx*delfx + delfy*delfy + delfz*delfz;

/*****	Get (1/r dV/dr) ******/
				if (waterP[0] == 'W'){ /* watts water */
				//	dedr = watts(r2,addr[m][n],&eg);
				}
				else if (waterP[0] == 'S'){ /* spc water */
					dedr = spc(r2,m,n,&eg,&wc);
				}
				else {
					fprintf(stderr, "waterForce: water model undefined\n");
					exit(1);
				}

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
			wc *= s;
			}
		k = i*3;
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

		force[l].fx += f[3].fx;
		force[l].fy += f[3].fy;
		force[l].fz += f[3].fz;
		force[++l].fx += f[4].fx;
		force[l].fy += f[4].fy;
		force[l].fz += f[4].fz;
		force[++l].fx += f[5].fx;
		force[l].fy += f[5].fy;
		force[l].fz += f[5].fz;
		VWATTS += eg;
		H2OC += wc;
		}
	}
}
#ifdef NEWSPC
#define	RE	1.000
#define	AE	1.9106332
#endif

intrawater(pos,force)
register tripd	*pos,*force;
{
double		acos(), sqrt();
double		dedd[3];
register double	d1, d2, da;
register double	r1, r2;
register double	r1r2, n1n2;
tripd		r[2], grad[3];

/******	first we need the vector differences ******/
r[0].fx = pos[0].fx - pos[1].fx;
r[0].fy = pos[0].fy - pos[1].fy;
r[0].fz = pos[0].fz - pos[1].fz;
r[1].fx = pos[0].fx - pos[2].fx;
r[1].fy = pos[0].fy - pos[2].fy;
r[1].fz = pos[0].fz - pos[2].fz;

/******	next we need r1, r2, and alpha (which requires a dot product) ******/
r1 = sqrt(r[0].fx*r[0].fx + r[0].fy*r[0].fy + r[0].fz*r[0].fz);
r2 = sqrt(r[1].fx*r[1].fx + r[1].fy*r[1].fy + r[1].fz*r[1].fz);
r1r2 =   (r[0].fx*r[1].fx + r[0].fy*r[1].fy + r[0].fz*r[1].fz);
n1n2 = r1r2 / (r1*r2);

/******	next come the useful variables (as they're defined in the pot.) ******/
d1 = r1 - RE;
d2 = r2 - RE;
da = RE * (acos(n1n2) - AE);

/******	we need the derivatives w.r.t. d1, d2, da of WATERV ******/
intrapot(d1, d2, da, dedd);

/******	now we need to calculate the derivs of d1, d2, da ******/

/******	first d1 ******/
dedd[0] /= r1;
force[1].fx += (grad[0].fx = dedd[0] * r[0].fx);
force[1].fy += (grad[0].fy = dedd[0] * r[0].fy);
force[1].fz += (grad[0].fz = dedd[0] * r[0].fz);
force[0].fx -= grad[0].fx;
force[0].fy -= grad[0].fy;
force[0].fz -= grad[0].fz;

/******	next d2 ******/
dedd[1] /= r2;
force[2].fx += (grad[0].fx = dedd[1] * r[1].fx);
force[2].fy += (grad[0].fy = dedd[1] * r[1].fy);
force[2].fz += (grad[0].fz = dedd[1] * r[1].fz);
force[0].fx -= grad[0].fx;
force[0].fy -= grad[0].fy;
force[0].fz -= grad[0].fz;

/******	finally da ******/
dedd[2] *= RE / (sqrt(1.-n1n2*n1n2)*r1*r2);
force[1].fx -= (grad[1].fx = dedd[2]*(r[1].fx - r1r2*r[0].fx/(r1*r1)));
force[1].fy -= (grad[1].fy = dedd[2]*(r[1].fy - r1r2*r[0].fy/(r1*r1)));
force[1].fz -= (grad[1].fz = dedd[2]*(r[1].fz - r1r2*r[0].fz/(r1*r1)));
force[2].fx -= (grad[2].fx = dedd[2]*(r[0].fx - r1r2*r[1].fx/(r2*r2)));
force[2].fy -= (grad[2].fy = dedd[2]*(r[0].fy - r1r2*r[1].fy/(r2*r2)));
force[2].fz -= (grad[2].fz = dedd[2]*(r[0].fz - r1r2*r[1].fz/(r2*r2)));
force[0].fx += grad[1].fx + grad[2].fx;
force[0].fy += grad[1].fy + grad[2].fy;
force[0].fz += grad[1].fz + grad[2].fz;
}

#define	K1	( 1213.61 / KCAL)
#define	K2	(-14.4990 / KCAL)
#define	K3	( 109.245 / KCAL)
#define	K4	( 32.7304 / KCAL)
#define	K5	(-1370.95 / RE / KCAL)
#define	K6	(-45.937 / RE / KCAL)
#define	K7	( 22.969 / RE / KCAL)
#define	K8	(-94.746 / RE / KCAL)
#define	K9	( 21.533 / RE / KCAL)
#define	K10	(-20.0976 / RE / KCAL)
#define	K11	( 2210.74 / RE / RE / KCAL)
#define	K12	( 114.84 / RE / RE / KCAL)
#define	K13	( 186.62 / RE / RE / KCAL)
#define	K14	(-244.04 / RE / RE / KCAL)
#define	K15	(-71.777 / RE / RE / KCAL)
#define	U0	(13.25 / KCAL)

intrapot(d1,d2,da,dedd)
register double	d1,d2,da;
double		*dedd;
{
register double	d12,d22,da2,d13,d23,d1d2,d12d22,dd12;

d12 = d1*d1;
d22 = d2*d2;
dd12 = d1*d2;
da2 = da*da;
d13 = d1*d12;
d23 = d2*d22;
d1d2 = d1 + d2;
d12d22 = d12 + d22;

WATERV += (0.5*((K1*d12d22) + (2.0*K2*dd12) + (K3*da2) + (2.0*K4*d1d2*da))
		+ (K5*(d13 + d23)) + (K6*d1d2*dd12) + (K7*d12d22*da)
		+ (K8*dd12*da) + (K9*d1d2*da2)
		+ (K10*da2*da) + (K11*(d12*d12 + d22*d22))
		+ (K12*d12d22*dd12) + (K13*d12*d22)
		+ (K14*d12d22*da2) + (K15*dd12*da2));

dedd[0] = ((K1*d1) + (K2*d2) + (K4*da) + (3.0*K5*d12) + (K6*(2.0*dd12 + d22))
		+ (2.0*K7*d1*da) + (K8*d2*da) + (K9*da2) + (4.0*K11*d13)
		+ (K12*(3.0*d12*d2 + d23)) + (2.0*K13*d1*d22)
		+ (2.0*K14*d1*da2) + (K15*d2*da2));

dedd[1] = ((K1*d2) + (K2*d1) + (K4*da) + (3.0*K5*d22) + (K6*(2.0*dd12 + d12))
		+ (2.0*K7*d2*da) + (K8*d1*da) + (K9*da2) + (4.0*K11*d23)
		+ (K12*(3.0*d22*d1 + d13)) + (2.0*K13*d2*d12)
		+ (2.0*K14*d2*da2) + (K15*d1*da2));

dedd[2] = ((K3*da) + (K4*d1d2) + (K7*d12d22) + (K8*dd12)
		+ (2.0*K9*d1d2*da) + (3.0*K10*da2)
		+ (2.0*K14*d12d22*da) + (2.0*K15*dd12*da));
}

int watts_inc, watts_cum;

double
watts(r2, buf, eg)
register double	*buf, r2, *eg;
{
register int	id;
register double	dz, der;

watts_inc++;
id =  1;
dz =  r2 - buf[id];
while	(dz < 0.)
	{
	watts_cum++;
	id += 12;
	dz =  r2 - buf[id];
	}
der = buf[id+1] + dz *(buf[id+2] + dz *(buf[id+3]
		 + dz *(buf[id+4] + dz * buf[id+5])));
*eg += buf[id+6] + dz *(buf[id+7] + dz *(buf[id+8] + dz *(buf[id+9]
		 + dz *(buf[id+10] + dz * buf[id+11]))));
return(der);
}

double
watts_calc(id, r2, buf, eg)
register int	id;
register double	*buf, r2, *eg;
{
register double	dz, der;

dz =  r2 - buf[id];
der = buf[id+1] + dz *(buf[id+2] + dz *(buf[id+3]
		 + dz *(buf[id+4] + dz * buf[id+5])));
*eg += buf[id+6] + dz *(buf[id+7] + dz *(buf[id+8] + dz *(buf[id+9]
		 + dz *(buf[id+10] + dz * buf[id+11]))));
return(der);
}

watts_seek(r2, buf, eg)
register double	*buf, r2, *eg;
{
register int	id;
register double	dz;

watts_inc++;
id =  1;
dz =  r2 - buf[id];
while	(dz < 0.)
	{
	watts_cum++;
	id += 12;
	dz =  r2 - buf[id];
	}
return(id);
}

#define	OO_SIGMA2	(3.16554*3.16554)
#define	OO_EPSIL4	(4*0.1554/KCAL)
#define HQ	0.41
#define	OQ	(-0.82)
double
spc(r2,m,n,eg,wc)
double r2, *eg, *wc;
int m,n;
{
int i;
double der,ir,ir2,ir6, EC,sqrt();

ir = 1./sqrt(r2);
if (m == 0){ /*Oxygen*/
	if (n == 0){ /*O-O*/
		EC = OQ*OQ*ir/E2;
		ir2 = OO_SIGMA2/r2;
		ir6 = ir2*ir2*ir2;
		*wc += EC;
		*eg += OO_EPSIL4*ir6*(ir6-1.) + EC;
		der = OO_EPSIL4*ir6*(6.-12*ir6)/r2 - EC/r2;
		i = 20.0/ir;
		waterRdf[i] += 1.0;
	}
	else { /*O-H*/
		EC = OQ*HQ*ir/E2;
		*eg += EC;
		*wc += EC;
		der = -EC/r2;
	}
}
else { /*Hydrogen*/
	if (n == 0){ /*H-O*/
		EC = OQ*HQ*ir/E2;
		*eg +=  EC;
		*wc += EC;
		der = -EC/r2;
	}
	else { /*H-H*/
		EC = HQ*HQ*ir/E2;
		*eg += EC;
		*wc += EC;
		der = -EC/r2;
	}
}
return(der);
}
