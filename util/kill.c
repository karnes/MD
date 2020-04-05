/*
 *	kill:
 *
 *	This routine freezes the atoms numbered from index to (index + n - 1).
 */

#include	<md.h>

kill(index,n)
	int	index,n;
	{
	register int	i;

	if	(index + n > natoms)
		{
		fprintf(stderr,"KILL: Not that many atoms\n");
		exit(1);
		}
	for	(i = index; i < index + n; i++)
		{
		atom[i].flags |= A_FIXED;
		vel[i].fx = 0.0;
		vel[i].fy = 0.0;
		vel[i].fz = 0.0;
		}
	}
