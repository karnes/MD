/*
 *	Verlet (second order) integrator
 *
 *	Standard Variables:
 *		h	time increment
 *		h1	previous integration step
 *		h2	previous elapsed real time step
 *		etime	elapsed real time
 *		pos	positions
 *		vel	(pos(t+h1) - pos(t))/h1
 *		tvel	true velocity
 *		force	forces
 *		K	instantaneous kinetic energy
 *		V	potential energy
 *		E	total energy
 *		H	Hamiltonian
 *
 *	Constant Temperature Algorithm Variables:
 *		g	number of degrees of freedom
 *		Teq	equilibrium temperature
 *		s	scaling variable
 *		vs	scaling variable's equivalent to "vel"
 *		Q	scaling variable's "mass"
 *
 *	Reference for Verlet algorithm:
 *		Swope, Anderson, ...
 *
 *	Reference for constant temperature algorithm:
 *		Suichi Nose, "A molecular dynamics method for simulations
 *		in the canonical ensemble," Molecular Physics, 52(2),255(1984)
 *
 */
#include	<md.h>
#include	<math.h>

integ(ngo, fix_T, fix_h)
int	ngo, fix_T, fix_h;
/*
 *	ngo	is number of time steps to integrate
 *	fix_T 	is nonzero if the constant temperature algorithm is used
 *	fix_h	is nonzero if the elapsed time increment is to be held fixed
 */
{
//printf("enter verlet\n");
double	alpha, beta, h_2, gkT, is2, sdot, fs;
int	i, j, k;
tripdouble	image;
double		log();
gkT = (double) (g + 1) * KBOLTZ * Teq;

while 	(ngo--)
	{
/******	Increment elapsed time, initialize other variables ******/
	etime += h2;
	H = K = KLIQ = 0.0;
	h_2 = h1 * .5;

/****** increment positions ******/
	for 	(i = 0; i < natoms; i++)
		{
		pos[i].fx += h1 * vel[i].fx;
		pos[i].fy += h1 * vel[i].fy;
		pos[i].fz += h1 * vel[i].fz;
		}
	if 	(natoms > nsolute) 
	        {
/******	move through periodic boundary conditions ******/

		for	(i = 0; i < natoms; i++)
			{
			if	((atom[i].flags & A_MAJOR) && atom[i].param1)
				{
				image.fx = -pos[i].fx;
				image.fy = -pos[i].fy;
				image.fz = -pos[i].fz;
				mvimage(&pos[i]);
				image.fx += pos[i].fx;
				image.fy += pos[i].fy;
				image.fz += pos[i].fz;
/*
				if (abs (image.fx) > 0.1)
				fprintf(stderr,"im.fx = %f\n",image.fx);
				if (abs (image.fy) > 0.1)
				fprintf(stderr,"im.fy = %f\n",image.fy);
				if (abs (image.fz) > 0.1)
				fprintf(stderr,"im.fz = %f\n",image.fz);
*/
				for 	(j=i+1; j <= i + atom[i].param1; j++)
					{
					pos[j].fx += image.fx;
					pos[j].fy += image.fy;
					pos[j].fz += image.fz;
					}
				}
			else if	(atom[i].flags & A_MAJOR)
				mvimage(&pos[i]);
			}
		}
		constrain('r'); /* correct positions	*/
/****** get force (-dV/dr) ******/
	getforce(pos, force);
/*	if (!ngo) wfsolute(stderr,1);*/
	for 	(i = 0; i < natoms; i++)
		if 	(atom[i].flags & A_FIXED)
			force[i].fx = force[i].fy = force[i].fz = 0.;


	if	(!fix_T)	/* normal verlet algorithm */
		{
		if	(ngo)	/* not returning yet */
			{
			for	(i = 0; i < natoms; i++)
				{
/****** just calculate accelerations ******/
				force[i].fx *= imass[i];
				force[i].fy *= imass[i];
				force[i].fz *= imass[i];
				}
			}
		else		/* will be returning */
			{
			for	(i = 0; i < natoms; i++)
				{
/****** calculate true velocities and accelerations ******/
				tvel[i].fx = vel[i].fx +
					h_2 * (force[i].fx *= imass[i]);
				tvel[i].fy = vel[i].fy +
					h_2 * (force[i].fy *= imass[i]);
				tvel[i].fz = vel[i].fz +
					h_2 * (force[i].fz *= imass[i]);

/****** Kinetic energy calculated from true velocities ******/
				K += mass[i] * (tvel[i].fx * tvel[i].fx +
						tvel[i].fy * tvel[i].fy +
						tvel[i].fz * tvel[i].fz);
				if (i == natoms - nsolute -1)
				   /* we just finish the accumulation of
				    * the kinetic energy of the solvent atoms*/
				    KLIQ = K/2.;
				}
		constrain('t');/* adjust kinetic energy. tvel is not adjusted */
		
			}
		h1 = h;		/* next time step */
		h_2 += h1 * .5;
		h2 = h;
		}
	else			/* constant temperature */
		{
		is2 = 1./ (s*s);
		alpha = s * s / (s + h1*vs);

		for 	(i = 0; i < natoms; i++)
			{
			force[i].fx *= imass[i] * is2;
			force[i].fy *= imass[i] * is2;
			force[i].fz *= imass[i] * is2;
			tvel[i].fx = (vel[i].fx + h_2 * force[i].fx) * alpha;
			tvel[i].fy = (vel[i].fy + h_2 * force[i].fy) * alpha;
			tvel[i].fz = (vel[i].fz + h_2 * force[i].fz) * alpha;
			K += mass[i] * (tvel[i].fx * tvel[i].fx +
					tvel[i].fy * tvel[i].fy +
					tvel[i].fz * tvel[i].fz);
			}

		sdot = vs + h_2 * (K - gkT) / (Q * s);
		alpha =  (s + h1 * vs) / (s + h1 * sdot);
		K *= alpha * alpha;
		fs = (K - gkT) / (Q * s);
		beta = 2. * sdot * is2;
		 
		 if 	(!ngo)
			H = .5 * Q * sdot * sdot + gkT * log(s);

		for	(i = 0; i < natoms; i++)
			{
			force[i].fx -= beta * (tvel[i].fx *= alpha);
			force[i].fy -= beta * (tvel[i].fy *= alpha);
			force[i].fz -= beta * (tvel[i].fz *= alpha);
			}

		if	(fix_h)	/* elapsed time step is fixed */
			{
			alpha = 2. * h * s / (1+sqrt(1 - 2. *h *sdot));
			beta = s * s - alpha * (sdot * s * .5 +
				alpha * (2. * sdot * sdot - fs * s) / 3.);
			h1 = h * s * s * s / beta;
			}
		else	/* integration step is fixed */
			h1 = h;

/****** calculate elapsed time ******/
		beta = s * s - h1 * (sdot * s * .5 +
			h1 * (2. * sdot * sdot - fs * s) / 3.);
		h2 = h1 * beta / (s * s * s);

		h_2 += h1 * .5;

/****** increment scale parameter ******/
		vs += h_2 * fs;
		s += h1 * vs;
		}

/****** both algorithms join here ******/

/****** increment velocities ******/
	for 	(i = 0; i < natoms; i++)
		{
		vel[i].fx += h_2 * force[i].fx;
		vel[i].fy += h_2 * force[i].fy;
		vel[i].fz += h_2 * force[i].fz;
		}
	constrain('v');/* correct velocites */
	}      /***** end_while ******/
temp = K / (g * KBOLTZ);	/* instantaneous temperature */
K *= .5;
H += (E = K + V);
//printf("exit verlet\n");
}
