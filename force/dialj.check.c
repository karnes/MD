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
	register double	r, ir, ir6, a, aa, b, bb, q2;
	register double	r2, pe, f;
	double	s, sp;
	double sqrt();
	double colint[256][256], sumi;

	a = lj[0][0].a;
	aa = 12.*a;
	b = lj[0][0].b;
	bb = 6.*b;
	for	(i = 0; i < n-2; i += 2)
		{
		for	(j = i+2; j < n; j += 2) {
		/* i and i+1 are the two diatoms and so are j and j+1 	*
		 * if the atom type is 101 we have + charge on it	*
		 * if the atom type is 102 we have - charge on it	*
		 * if the atom type is 103 we have no charge on it	*
		 * we need to consider the 4 interactions :		*
		 * i with j, i with j+1, i+1 with j and i+1 with j+1	*/
			colint[i/2][j/2] = colint[j/2][i/2] = 0.;
			for (ii =0; ii<2; ii++)
			    {
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
					if(atom[i+ii].type != 103 && atom[j+jj].type != 103) 	
					{
					q2 = lj[0][0].q *
					(203 - 2*atom[i+ii].type) * (203 - 2*atom[j+jj].type);
					r = sqrt(r2);
					pe = q2/r;
					VLIQ += pe;
					colint[i/2][j/2] += pe;
					f += pe/r2;
					}
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
					if(atom[i+ii].type != 103 && atom[j+jj].type != 103) 	
					{
					q2 = lj[0][0].q *
					(203 - 2*atom[i+ii].type) * (203 - 2*atom[j+jj].type);
					r = sqrt(r2);
					pe = q2/r; 
					f += s*pe/r2 - sp * pe;
					VLIQ += s * pe;
					colint[i/2][j/2] += pe*s;
					}
				}
				frc[i+ii].fx += (fx *= f);
				frc[i+ii].fy += (fy *= f);
				frc[i+ii].fz += (fz *= f);

				frc[j+jj].fx -= fx;
				frc[j+jj].fy -= fy;
			        frc[j+jj].fz -= fz;
                            }
			}
		    }
		}
	for (i=0; i <n/2; i++)
	    {
		sumi = 0.;
		for (j=0; j <n/2; j++)
			sumi += (colint[i][j] + colint[j][i]);
	        fprintf(stderr," %d %f\n",i,sumi*KCAL);
	     }
}
