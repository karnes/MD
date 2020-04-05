/*
 *
 *	This routine calculates the liquid-liquid interactions for a polar
 *	or nonpolar diatomic Lennard-Jonesium.
 */
#include	<md.h>

dialj(n, pos, frc)
	register int	n; /*	number of liquid atoms 	*/
	register tripd	*pos;
	register tripd	*frc;
{
	register int	i, j, ii, jj;
	register double	fx, fy, fz;
	register double	ir, ir6, a, aa, b, bb, q2;
	register double	r2, pe, f;
	double	s, sp;
	double add_coul();

	a = lj[0][0].a;
	aa = 12.*a;
	b = lj[0][0].b;
	bb = 6.*b;
	q2 = lj[0][0].q;
	for	(i = 0; i < n-2; i += 2)
		for	(j = i+2; j < n; j += 2) {
		/* i and i+1 are the two diatoms and so are j and j+1 	*
		 * if the atom type is 101 we have + charge on it	*
		 * if the atom type is 102 we have - charge on it	*
		 * if the atom type is 103 we have no charge on it	*
		 * we need to consider the 4 interactions :		*
		 * i with j, i with j+1, i+1 with j and i+1 with j+1	*/
			for (ii =0; ii<2; ii++)
			    for (jj =0; jj<2; jj++) {
				/*
				 *	this next part is the imaging.
			 	*/
				fx = (pos[i+ii].fx - pos[j+jj].fx) * iperiod.fx;
				fy = (pos[i+ii].fy - pos[j+jj].fy) * iperiod.fy;
				fz = (pos[i+ii].fz - pos[j+jj].fz) * iperiod.fz;

				fx -= (int) (2 * fx);
				fy -= (int) (2 * fy);
				fz -= (int) (2 * fz);
#ifdef PTO /* truncated octahedral boundaries */
				if	((abs(fx) + abs(fy) + abs(fz)) > 0.75) {
					fx -= hsgn(fx);
					fy -= hsgn(fy);
					fz -= hsgn(fz);
				}
#endif PTO
				fx *= period.fx;
				fy *= period.fy;
				fz *= period.fz;

				r2 =fx*fx + fy*fy + fz*fz;

				if	(r2 >= swr2max)
					continue;
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
					f = s*(aa * ir6 - bb )*ir6*ir - sp * pe;
					VLIQ += s * pe;
				}
				f += add_coul(r2,q2,i+ii,j+jj);
				frc[i+ii].fx += (fx *= f);
				frc[i+ii].fy += (fy *= f);
				frc[i+ii].fz += (fz *= f);

				frc[j+jj].fx -= fx;
				frc[j+jj].fy -= fy;
			        frc[j+jj].fz -= fz;
                            }
			}
}
/* calculate coulomb interactions */
double add_coul(r2,q2,i,j)
double r2,q2;
int i,j;
{
double f, r, pe, s, sp, sqrt();
	if	(atom[i].type == 103 ||  atom[j].type == 103 ) 	
		return(0.);
	/* The following is specific to a system where all the
	 * site are identical except for the sign of the
	 * charge. It is assumed that atom type 101 carry
	 * positive charge and atom 102 negative. The 
	 * lj parameters in potentials.h contains the 
	 * correct charge for the general case
	 */
	q2 *= (203 - 2*atom[i].type) * (203 - 2*atom[j].type);
	r = sqrt(r2);
	if	(r2 <= swr2min) {
		f = q2/r;
		VLIQ += f;
		f /= r2;
		}
		else {
			pe = r2-swr2min;
			f = pe * pe;
			s = swcoef[0] + pe * f *
		        (swcoef[1]+pe*swcoef[2]+f*swcoef[3]);
			sp = f *
			(pswcoef[0]+pe*pswcoef[1]+f*pswcoef[2]);
			pe = q2/r; 
			f = s*pe/r2 - sp * pe;
			VLIQ += s * pe;
		}		       
	return(f);
}
