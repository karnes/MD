/*
 *	mmlj612:
 *
 *	This routine calculates the liquid-liquid interactions, assuming
 *	Lennard-Jones potentials.
 */
#include	<md.h>

mmlj612(n, pos, frc)
	register int	n;
	register tripd	*pos;
	register tripd	*frc;
{
	register int	i, j, ki, kj;
	register double	fx, fy, fz, fxy, fz2;
	register double	ir, ir6, a, aa, b, bb;
	register double	r2, pe, f;
	double	s, sp, sqrt();

	a = lj[0][0].a;
	aa = 12.*a;
	b = lj[0][0].b;
	bb = 6.*b;
	for	(i = 0; i < n-1; i++) {
		ki = (int) ((pos[i].fz + zwall)/pslabSize);
		pTensor[ki].fx += 1.;
		for	(j = i+1; j < n; j++) {
			/*
			 *	this next part is the imaging.
			 */
			kj = (int) (pos[j].fz + zwall)/pslabSize;
			fx = (pos[i].fx - pos[j].fx) * iperiod.fx;
			fy = (pos[i].fy - pos[j].fy) * iperiod.fy;
			fz = (pos[i].fz - pos[j].fz) * iperiod.fz;

			fx -= (int) (2 * fx);
			fy -= (int) (2 * fy);
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
				ERROR((stderr,"mmlj612:imaging error\n"), exit);
			fx *= period.fx;
			fy *= period.fy;
			fz *= period.fz;

			fz2 = fz*fz;
			fxy = fx*fx + fy*fy;
			r2 =fxy + fz2;

			if	(r2 >= swr2max)
				continue;
/***	Radial distribution function calculated here
	This is an expensive addition, so the use of the ifdef	***/
#ifdef RDF
			rdf[(int) (sqrt(r2)*20)] += 2.;
#endif RDF
			ir = 1. / r2;
			ir6 = ir * ir * ir;

			if	(r2 <= swr2min) {
				VLIQ += ( a * ir6 - b ) * ir6;
				f = (aa * ir6 - bb ) * ir6 * ir;
			}
			else {
				pe = r2-swr2min;
				f = pe * pe;
				s = swcoef[0] + pe * f *
					(swcoef[1]+pe*swcoef[2]+f*swcoef[3]);
				sp = f *
					(pswcoef[0]+pe*pswcoef[1]+f*pswcoef[2]);
				pe = ( a * ir6 - b ) * ir6;
				f = s * (aa * ir6 - bb ) * ir6 * ir - sp * pe;
				VLIQ += s * pe;
			}

			pTensor[ki].fy += f*fxy/4.;
			pTensor[kj].fy += f*fxy/4.;
			pTensor[ki].fz += f*fz2/2.;
			pTensor[kj].fz += f*fz2/2.;
			frc[i].fx += (fx *= f);
			frc[i].fy += (fy *= f);
			frc[i].fz += (fz *= f);

			frc[j].fx -= fx;
			frc[j].fy -= fy;
			frc[j].fz -= fz;
		}
	}
}
