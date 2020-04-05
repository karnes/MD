/*
 *	This procedure rotates the molecule a random amount for placement
 *	by addmol.
 */


#include	<math.h>
#include	"typedefs.h"

#define RMAX	2147483648.0	/* random() returns in range [0..RMAX-1] */
#define TWO_PI	6.283185/RMAX	/* 2 PI / RMAX */
#define TWO	2.0/RMAX	/* 2 / RMAX */

rotcoord(coord,n)
	tripfloat	*coord;
	int		n;
{
	int		i;
	tripfloat	temp;
	float		tmat[3][3];
	float		phi, theta, psi;
	float		cphi,ctheta,cpsi,sphi,stheta,spsi;

	phi = random() * TWO_PI;
	theta = acos(random() * TWO - 1.0);
	psi = random() * TWO_PI;
	cphi = cos(phi);
	ctheta = cos(theta);
	cpsi = cos(psi);
	sphi = sin(phi);
	stheta = sin(theta);
	spsi = sin(psi);
	tmat[0][0] = ctheta * cphi;
	tmat[0][1] = ctheta * sphi;
	tmat[0][2] = -stheta;
	tmat[1][0] = spsi * stheta * cphi - cpsi * sphi;
	tmat[1][1] = spsi * stheta * sphi + cpsi * cphi;
	tmat[1][2] = ctheta * spsi;
	tmat[2][0] = cpsi * stheta * cphi + spsi * sphi;
	tmat[2][1] = cpsi * stheta * sphi - spsi * cphi;
	tmat[2][2] = ctheta * cpsi;
/*
 * perform rotation
 */
	while	(n--) {
		temp.fx = coord->fx * tmat[0][0] +
			  coord->fy * tmat[0][1] +
			  coord->fz * tmat[0][2];
		temp.fy = coord->fx * tmat[1][0] +
			  coord->fy * tmat[1][1] +
			  coord->fz * tmat[1][2];
		temp.fz = coord->fx * tmat[2][0] +
			  coord->fy * tmat[2][1] +
			  coord->fz * tmat[2][2];
		coord->fx = temp.fx;
		coord->fy = temp.fy;
		coord->fz = temp.fz;
		coord++;
	}
}
