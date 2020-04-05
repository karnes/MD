/*
 *
 *	This routine calculates the liquid-liquid interactions for a polar
 *	or nonpolar diatomic Lennard-Jonesium.
 */

#include	<md.h>
#include	<system.h>

dialj(n, pos, frc)
	register int	n; /*	number of liquid atoms 	*/
	register tripd	*pos;
	register tripd	*frc;
{
	register int	i, j, ii, jj;
	register double	fx, fy, fz;
	register double	r, ir, ir6, a, aa, b, bb, q2;
	register double	r2, pe, f, ec;
	tripd vecAB, S;
	double	s, sp;
	double sqrt();

	a = ljDiaA;
	aa = 12.*a;
	b =  ljDiaB;
	bb = 6.*b;
	/* double loop over molecules	*/
	for	(i = 0; i < n-2; i += 2)
		for	(j = i+2; j < n; j += 2) {
		 /* The truncation of the interaction between any two    *
		 * molecules is determined according to the distance of *
		 * center of mass between the two molecules. All the 4  *
		 * pair interactions are smoothed using the same factor. *
		 * The diaCM is calculated in liqForce.c		*/

		 vecAB.fx = diaCM[i/2].fx - diaCM[j/2].fx;
		 vecAB.fy = diaCM[i/2].fy - diaCM[j/2].fy;
		 vecAB.fz = diaCM[i/2].fz - diaCM[j/2].fz;
		 mvimage(&vecAB);
		 r2 = vecAB.fx*vecAB.fx+vecAB.fy*vecAB.fy+vecAB.fz*vecAB.fz; 
		 if ( r2 >= swr2max )
		 	continue; /* to the next diatom pair*/
		 if ( r2 <= swr2min ){
			s = 1.;
			sp = 0.;
		 }
		else 
			swtch(r2 - swr2min, &s, &sp);
		pe = 0; /* pe is the total interaction energy between two
			 * diatomic molecules. We need it to calculate the
			 * forces due to the switching function
			 */
		/* Double loop over atoms in molecules:			*
		 * i and i+1 are the two diatoms and so are j and j+1 	*
		 * if the atom type is 101 we have + charge on it	*
		 * if the atom type is 102 we have - charge on it	*
		 * if the atom type is 103 we have no charge on it	*
		 * we need to consider the 4 interactions :		*
		 * i with j, i with j+1, i+1 with j and i+1 with j+1	*/
			for (ii =0; ii<2; ii++)
			    for (jj =0; jj<2; jj++) {
/*
 *	This next part is the imaging on the atom-atom distances.
 */
				fx = (pos[i+ii].fx - pos[j+jj].fx) * iperiod.fx;
				fy = (pos[i+ii].fy - pos[j+jj].fy) * iperiod.fy;
				fz = (pos[i+ii].fz - pos[j+jj].fz) * iperiod.fz;

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

				r2 =fx*fx + fy*fy + fz*fz;
				ir = 1. / r2;
				ir6 = ir * ir * ir;
				pe +=  ( a * ir6 - b ) * ir6;
				f = (aa * ir6 - bb ) * ir6 * ir *s;
		/***	Add coulombic interactions	*/
				if(atom[i+ii].type!=103&&atom[j+jj].type!=103){	
					q2=sq(sCharge)*(203-2*atom[i+ii].type)*
						(203 - 2*atom[j+jj].type);
					r = sqrt(r2);
					ec = q2/(E2*r);
					pe += ec;
					f += s*ec/r2;
				}
				/* add forces due to inter atom interactios*/
				frc[i+ii].fx += (fx *= f);
				frc[i+ii].fy += (fy *= f);
				frc[i+ii].fz += (fz *= f);

				frc[j+jj].fx -= fx;
				frc[j+jj].fy -= fy;
			        frc[j+jj].fz -= fz;
                            }/* end loop over this diatom pair*/
			/* add forces due to switching function */
			S.fx = -0.5*sp*pe*vecAB.fx; /* the factor of 0.5 here*/
			S.fy = -0.5*sp*pe*vecAB.fy; /* is because of the equal*/
			S.fz = -0.5*sp*pe*vecAB.fz; /* masses of all atoms    */
			frc[i].fx += S.fx;
			frc[i+1].fx += S.fx;
			frc[j].fx -= S.fx;
			frc[j+1].fx -= S.fx;
			frc[i].fy += S.fy;
			frc[i+1].fy += S.fy;
			frc[j].fy -= S.fy;
			frc[j+1].fy -= S.fy;
			frc[i].fz += S.fz;
			frc[i+1].fz += S.fz;
			frc[j].fz -= S.fz;
			frc[j+1].fz -= S.fz;
			VLIQ += pe*s;
			}
}
