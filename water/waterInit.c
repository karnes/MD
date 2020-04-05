#include	<md.h>
#include	<math.h>
#include	<water.h>
#include	<atomtypes.h>
#define 	en_prec 1.0E-5

/*
 *	this subroutine initializes the polynomials for the water routines.
 */
waterInit()
{
if (waterP[0] == 'W'){ /* watts water */
	fprintf(stderr, "waterInit.c: Watts water disabled. exiting...\n");
	exit(1);
/*	oo_poly();
	oh_poly();
	hh_poly();
*/
}
else if (waterP[0] == 'S'){ /* spc water */
	return;
}
else {
		fprintf(stderr, "waterForce: water model undefined\n");
		exit(1);
}
}
