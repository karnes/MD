/*
 *	aafe612:
 *
 *	This routine calculates the system-liquid interaction, assuming
 *	Lennard-Jones potentials.
 */

#include	<md.h>

aafe612(ri, rj, fi, fj, ti, tj)
	tripd	*ri, *rj;
	tripd	*fi, *fj;
	int	ti, tj;
{
	register double	r2, ir2, ir6;
	register double	ftmp;
	register double	E;
	double	sw, swp;
	tripd	rij;
	ljcon	ljp;

	rij.fx = ri->fx - rj->fx;
	rij.fy = ri->fy - rj->fy;
	rij.fz = ri->fz - rj->fz;
	mvimage(&rij);

	r2 = rij.fx*rij.fx + rij.fy*rij.fy + rij.fz*rij.fz;

	if	(r2 >= swr2max)
		return;
	if	(r2 <= swr2min)
		{
		sw = 1.;
		swp = 0.;
		}
	else
		swtch(r2 - swr2min, &sw, &swp);
	
	ljp = lj[ti][tj];
	ir2 = 1. / r2;
	ir6 = ir2*ir2*ir2;
	E = (ljp.a * ir6 - ljp.b) * ir6;

	ftmp = (12. * ljp.a * ir6 - 6. * ljp.b) * ir6 * ir2 * sw - E * swp;
	fi->fx += (rij.fx *= ftmp);
	fi->fy += (rij.fy *= ftmp);
	fi->fz += (rij.fz *= ftmp);
	fj->fx -= rij.fx;
	fj->fy -= rij.fy;
	fj->fz -= rij.fz;

	VINT += E * sw;	/* this is the total Vnb over all solvent atoms */
}
