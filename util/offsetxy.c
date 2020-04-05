/*
 *	offset:
 *	
 *	This routine offsets all atoms by an appropriate amount so as to bring
 *	the center of the last two atoms in the list to x=0, y=0
 */

#include	<md.h>

offset()
{
	int	i, j;
	tripd	image, off;
	off.fx = 0.5*(pos[natoms-1].fx+pos[natoms-2].fx);
	off.fy = 0.5*(pos[natoms-1].fy+pos[natoms-2].fy);

	for (i = 0; i < natoms; i++) {
		pos[i].fx -= off.fx;
		pos[i].fy -= off.fy;

		if ((atom[i].flags & A_MAJOR) && atom[i].param1) {
			image.fx = -pos[i].fx;
			image.fy = -pos[i].fy;
			image.fz = -pos[i].fz;
			mvimage(&pos[i]);
			image.fx += pos[i].fx;
			image.fy += pos[i].fy;
			image.fz += pos[i].fz;
			for (j=i+1; j <= i + atom[i].param1; j++) {
				pos[j].fx += image.fx;
				pos[j].fy += image.fy;
				pos[j].fz += image.fz;
			}
		}
		else if (atom[i].flags & A_MAJOR)
			mvimage(&pos[i]);
	}
}
