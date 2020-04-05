#include	<md.h>

/*
 *	This routine does the periodic imaging of any given
 *	vector for the parallelopiped or the truncated octahedral
 *      boundary conditions.
 */

mvimage(del)
	tripdouble	*del;
{
	double	fx, fy, fz;
	if (natoms == nsolute) /* no solvent */
	     return;
	fx = del->fx * iperiod.fx;  /* Normalizing the positions */
	fy = del->fy * iperiod.fy;  /* |fx|, |fy|, |fz| <= 0.5 for the  */
	fz = del->fz * iperiod.fz;  /* atom to be in the box */
	fx -= (int) fx;
	fy -= (int) fy;   /* bringing the particle to the nearest box */
	fz -= (int) fz;
	fx -= (int) (2 * fx);
	fy -= (int) (2 * fy);  /* bringing the particle to the box */
	fz -= (int) (2 * fz);
	if (pbcType[0] == 'O'){
	/* truncated octahedral boundaries */
		if	((abs(fx) + abs(fy) + abs(fz)) > 0.75) {
			fx -= hsgn(fx);
			fy -= hsgn(fy);
			fz -= hsgn(fz);
		}
	}
	else if (pbcType[0] != 'C')
		ERROR((stderr,"mvimage is set for either OCT or CUBE\n"), exit);

	del->fx = fx * period.fx;
	del->fy = fy * period.fy;
	del->fz = fz * period.fz;
}
