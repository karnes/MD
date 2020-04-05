/*
 *
 *	This routine calculates the contribution to total potential energy
 *	and the forces from the intramolecular vibrations of a diatomic liq.
 */
#include	<md.h>
#include	<system.h>

diaVib(n, pos, frc)
	int	n; /*	number of liquid atoms 	*/
	tripd	*pos;
	tripd	*frc;
{
	int	i;
	tripd	r;
	double	r1,r2,f, sqrt();
	VIBLIQ = 0.;
	for	(i = 0; i < n; i += 2){
	/*	atom i is bonded to atom i+1 */
		r.fx = (pos[i+1].fx - pos[i].fx);
		r.fy = (pos[i+1].fy - pos[i].fy);
		r.fz = (pos[i+1].fz - pos[i].fz);
		mvimage(&r);
		r2 = r.fx * r.fx + r.fy * r.fy + r.fz * r.fz;
		r1 = sqrt (r2);
		VIBLIQ += 0.5*vibF*(r1-equB)*(r1-equB);
		f =  vibF * (equB - r1)/r1;
		frc[i+1].fx += (r.fx *= f);
		frc[i+1].fy += (r.fy *= f);
		frc[i+1].fz += (r.fz *= f);
		frc[i].fx -= r.fx;
		frc[i].fy -= r.fy;
		frc[i].fz -= r.fz;
	}
}
