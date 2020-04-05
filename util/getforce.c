/*
 *	getforce:
 *
 *	calls the system dependent force routines.
 */

#include <md.h>

getforce(pos, force)
	tripd	*pos;
	tripd	*force;
{
//fprintf("getforce.c: width_w = %f\n",width_w);
	int	i;

	V = VSYS = VLIQ = VINT = 0.;
	for	(i = 0; i < natoms; i++) {
		force[i].fx = 0.0;
		force[i].fy = 0.0;
		force[i].fz = 0.0;
	}
//fprintf(stderr,"getforce.c\n");
/***	Do not change the order of these calls		***/
	liqForce(pos, force);	/* get liq-liq forces */
	intForce(pos, force);	/* get interaction liq-sys forces */
	sysForce(pos, force);	/* get system forces */
	V = VSYS + VLIQ + VINT;
//fprintf(stderr,"at end of getforce.c\n");
}
