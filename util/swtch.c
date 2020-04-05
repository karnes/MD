/*
 *	swtch:
 *
 *	This subroutine generates the switching function and its derivative
 *	given the value of zz.  it is normalized to the size of the box.
 */

#include	<md.h>

swtch(zz,s,sp)
register double	zz, *s, *sp;
	{
	register double	zz2;
	zz2 = zz * zz;
	*s = swcoef[0] + zz * zz2  * (swcoef[1] + zz*swcoef[2] + zz2*swcoef[3]);
	*sp = zz2 * (pswcoef[0] + zz*pswcoef[1] + zz2*pswcoef[2]);
	}
