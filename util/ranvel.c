/*
 *	ranvel:
 *
 *	This routine chooses the initial velocity configuration for the atoms
 *	from a MB distribution at temperature temp.  If rmcom is nonzero, 
 *	then the center of mass velocity is subtracted and the number of
 *	degrees of freedom is reduced by 3.
 */

#include	<md.h>
#include	<math.h>
#include	<sys/time.h>

#define	RMAX	2147483648.
#define	EQTIME	50.

ranvel(temp, rmcom, vScaleQ, newVQ)
double	temp;
int	rmcom;		/* if == 1 remove center of mass velocity */
int 	vScaleQ;	/* if == 1 scale velocity so temp = Teq	  */
int 	newVQ;		/* if == 0 do not select random velocities*/
{
int	i, j;
double	x1, x2;		/* uniform random variates */
double	sig;		/* normalization factor */
double	im;		/* for sqrt(imass[i]) */
union 	{
	tripd	v[2];	/* corresp. to vel[i] & vel[i+1] */
	double	y[6];	/* for generating gauss. RV's in pairs */
} r;
static int	init = 1;

/* set up arrays for orientational distributions */
void setupOD();  //jjk

if	(init) {
	struct timeval t;
	struct timezone tz;
	gettimeofday(&t, &tz);
	srandom((int) (t.tv_sec + t.tv_usec)); 
	init = random() * 100. / RMAX;
	for	(i = 0; i < init; i++)
		random();	/* warm it up */
	init = 0;
}

sig = sqrt(2. * KBOLTZ * temp);
h2 = h1 = 0.;		/* initialize last-integration variable */
g = 0;			/* number of degrees of freedom */
	
if	(newVQ) {
	for	(i = 0; i < natoms;) {
		for	(j = 0; j < 6; j += 2) {
			x1 = (double)random() / RMAX;
			x2 = (double)random() * 2. * PI / RMAX;
			x1 = sqrt(-log(1.-x1)) * sig;
			r.y[j] = x1 * cos(x2);
			r.y[j+1] = x1 * sin(x2);
		}
		for	(j = 0; j < 2; j++) {
			while	(i < natoms && (atom[i].flags & A_FIXED)) {
				vel[i].fx = vel[i].fy = vel[i].fz = 0.;
				i++;
			}
			if	(i == natoms)
				break;
			im = sqrt(imass[i]);
			vel[i].fx = r.v[j].fx * im;
			vel[i].fy = r.v[j].fy * im;
			vel[i].fz = r.v[j].fz * im;
			i++;
			g++;
		}
	}
}
else 	{ /* newVQ = 0. use the current velocities */	
	for 	(i = 0; i< natoms;i++) {
		g++;
		if	(atom[i].flags & A_FIXED) {
			if(vel[i].fx != 0. || vel[i].fy != 0. || vel[i].fz !=0.)
				fprintf(stderr,"ranvel:Warning- atom %d is killed, but its velocity is not zero.\n Velocity = (%f %f %f)\n",i,vel[i].fx,vel[i].fy,vel[i].fz);
			g--;
		}
	}
}
		
		

/****** if rmcom is set and nothing is killed, remove COM velocity ******/
if	(rmcom && (g == natoms)) {
	rmcm();
	g--;
}
	
/****** if something is killed, then COM velocity is zero ******/
else if	(g != natoms){
	if (g != nsolute)/* do not reduce g if all solvent is frozen*/
		g--;
}
g *= 3;			/* total degrees of freedom */
/****** remove velocity component along constraint bond	*****/
constrain('v');
g -= cons_num;	/* reduce number of degrees of freedom	*
		 * by the number of constraints		*/
etime = 0.;
tc = -1;/*time counter is increased in liqForce.c*/
s = 1.;
vs = 0.;
Teq = temp;

/****** Q = mass for EQTIME fs period of oscillation ******/
Q = KBOLTZ * EQTIME * EQTIME * .5 / (PI * PI) * Teq * g;
if (vScaleQ)
	{ /* scale velocities so that the temperature is Teq */
		double kinet, scale_factor;
		kinet = 0.; /* used to accumulate kinetic energy */
		/* calculate (twice the) kinetic energy: */
		for (i=0; i< natoms; i++)
			{
			kinet += mass[i] * (vel[i].fx * vel[i].fx +
				vel[i].fy * vel[i].fy +
				vel[i].fz * vel[i].fz);
			}
		temp = kinet/ (g * KBOLTZ);   /* current temperature */
		scale_factor = sqrt(Teq /temp);
		/* calculate new velocities */
		for (i=0; i< natoms; i++)
			{
				vel[i].fx *= scale_factor;
				vel[i].fy *= scale_factor;
				vel[i].fz *= scale_factor;
			}
	}
#ifdef	LEPS
j = natoms - 3;
getleps_r(&pos[j]);
if	(nbias == 2)
	{
	double	chqas();
	chqas(leps_r, &vel[j], &mass[j], 0., 'v', 2);
	g--;	/* constraint removes one degree of freedom */
	}
else
	{
	double	chqas();
	x1 = chqas(leps_r, &vel[j], &mass[j], 0., 'v', 1);
	x1 = abs(x1);
	chqas(leps_r, &vel[j], &mass[j], x1, 'v', 2);
	}
#endif
}

/*
 *	Remove Center-of-Mass velocity
 */
rmcm()
{
	double	comfx, comfy, comfz, mtot;
	int	i;

	comfx = comfy = comfz = mtot = 0.;

	for	(i = 0; i < natoms; i++) {
		comfx += mass[i] * vel[i].fx;
		comfy += mass[i] * vel[i].fy;
		comfz += mass[i] * vel[i].fz;
		mtot += mass[i];
	}
	comfx /= mtot;
	comfy /= mtot;
	comfz /= mtot;
	for	(i = 0; i < natoms; i++) {
		vel[i].fx -= comfx;
		vel[i].fy -= comfy;
		vel[i].fz -= comfz;
	}
}
