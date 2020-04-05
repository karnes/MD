#include	<stdio.h>
#include	<math.h>
#define PI              3.14159265358979323846
#define E2 7.197587876
#define KCAL (1.e4/4.184)
typedef struct  {
	double  fx;
	double  fy;
	double  fz;
	}       tripd;  
double	swcoef[4];
double	pswcoef[3];
double swr2max,swr2min;

main(argc,argv)
	int argc;
	char	*argv[];
{
int i;
double	dz, dz3, dz4, dz5, bc_factor,xwall,d,t,p,Fi,Fo,pe,R;
tripd pos[4];
t = atof(argv[1])*PI/180.0;
p = atof(argv[2])*PI/180.0;
	xwall = 12.0;
	bc_factor = 1.; /* for cubic */
/*	bc_factor = 0.75; for truncated octahedron */
	swr2max = bc_factor * xwall * xwall;
	/* This is the square of 1/2 the distance between the hexagonal faces
	 * of the pto box, or the sq of 1/2 the small length of the cube */
	swr2min = bc_factor * (xwall - 1.0) * (xwall - 1.0);
	
	dz = swr2max - swr2min;   /* r**2 (min max) for switching region */
	dz3 = (dz*dz*dz);
	dz4 = dz3*dz;
	dz5 = dz4*dz;
	swcoef[0] =    1.0;       /* switching parameters for group cut-offs */
	swcoef[1] =  -10.0 / dz3;
	swcoef[2] =   15.0 / dz4;
	swcoef[3] =   -6.0 / dz5;
	pswcoef[0] = -60.0 / dz3;	/* 2 times swcoef's derivs */
	pswcoef[1] = 120.0 / dz4;
	pswcoef[2] = -60.0 / dz5;

	pos[0].fx = pos[0].fy = pos[0].fz = 0.0;
	pos[1].fx = -0.76;pos[1].fy = -0.47; pos[1].fz = 0.09;
	pos[2].fx = 0.69;pos[2].fy = -0.47; pos[2].fz = 0.35;
	for (i = 0; i<101; i++){
		d = 2+i*(xwall-2)/100.0;
		pos[3].fx = d*sin(t)*cos(p);
		pos[3].fy = d*sin(t)*sin(p);
		pos[3].fz = d*cos(t);
		pe = 0.0;
		h2oljc(&pos[0],&pos[3],&Fi,&Fo,&pe,&R);
		printf("%f %f %f %f\n",R,Fi*KCAL,Fo*KCAL,pe*KCAL);
	}
}

h2oljc(ri,rj,Fi,Fo,pe,R)
tripd *ri;/* ri[0], ri[1] and ri[2] are position of the O H1 and H2 of water*/
tripd *rj;/* rj[0] is the position of the solute atom*/
double *Fi;/* fi[0], fi[1] and fi[2] are forces on the O H1 and H2 of water*/
double *Fo;/* force on solute atom*/
double *pe;
double *R;
{
double s, sp, f, R2, q2, r2, r, ir, ir6, EC,qW;
tripd  rij, sdist,fi[3],fj[3];
int i;

/***	Determine Oxegen-solute image vector			***/

sdist.fx = ri[0].fx - rj[0].fx;
sdist.fy = ri[0].fy - rj[0].fy;
sdist.fz = ri[0].fz - rj[0].fz;
R2 = sdist.fx*sdist.fx +  sdist.fy*sdist.fy +  sdist.fz*sdist.fz;
*R = sqrt(R2);
if ( R2 >= swr2max )
	return;
if ( R2 <= swr2min ){
	s = 1.0;
	sp = 0.0;
 }
else
	swtch(R2 - swr2min, &s, &sp);

/***	Loop over atoms in water molecule			***/
for (i=0; i< 3; i++)
	fi[i].fx = fi[i].fy = fi[i].fz = 0.0;
fj[0].fx = fj[0].fy = fj[0].fz = 0.0;
for (i=0; i< 3; i++){
	if (i == 0) qW = -0.82;
	if (i > 0) qW = 0.41;
	rij.fx = ri[i].fx - rj[0].fx;
	rij.fy = ri[i].fy - rj[0].fy;
	rij.fz = ri[i].fz - rj[0].fz;
	r2 = rij.fx*rij.fx +  rij.fy*rij.fy +  rij.fz*rij.fz ;
       	q2 = qW/E2;
	ir = 1. / r2;
	r = sqrt(r2);
	EC = q2/r;
	*pe +=  EC;
	f = s*(EC/r2); /* -dV/dr/r	*/
	fi[i].fx += (rij.fx *= f);
	fi[i].fy += (rij.fy *= f);
	fi[i].fz += (rij.fz *= f);

	fj[0].fx -= rij.fx;
	fj[0].fy -= rij.fy;
	fj[0].fz -= rij.fz;
}
if (sp != 0.0){
	f =sp*(*pe);
	fi[0].fx -= (sdist.fx * f);
	fi[0].fy -= (sdist.fy * f);
	fi[0].fz -= (sdist.fz * f);

	fj[0].fx += (sdist.fx * f);
	fj[0].fy += (sdist.fy * f);
	fj[0].fz += (sdist.fz * f);
	*pe *= s;
	}
	*Fi = (fj[0].fx*sdist.fx+fj[0].fy*sdist.fy+fj[0].fz*sdist.fz)/sqrt(R2);
	*Fo = (fi[0].fx*sdist.fx+fi[0].fy*sdist.fy+fi[0].fz*sdist.fz)/sqrt(R2);
}
swtch(zz,s,sp)
double	zz, *s, *sp;
	{
	register double	zz2;
	zz2 = zz * zz;
	*s = swcoef[0] + zz * zz2  * (swcoef[1] + zz*swcoef[2] + zz2*swcoef[3]);
	*sp = zz2 * (pswcoef[0] + zz*pswcoef[1] + zz2*pswcoef[2]);
	}
