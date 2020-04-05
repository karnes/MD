/*
 *	Given two atom numbers, i and j, say, dist(i,j) gives the magnitude
 *	of the distance between the vectors pos[i] and pos[j], using the
 *	minimum image convention.
 */
#include	<md.h>

double dist(atom1, atom2)
	int	atom1;
	int	atom2;
{
	double	r, sqrt();
	tripd	rij;
	rij.fx = pos[atom1].fx - pos[atom2].fx;
	rij.fy = pos[atom1].fy - pos[atom2].fy;
	rij.fz = pos[atom1].fz - pos[atom2].fz;
	mvimage(&rij);

	r = sqrt(rij.fx*rij.fx + rij.fy*rij.fy + rij.fz*rij.fz);
	return (r);
}
