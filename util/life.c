/*
 *	life:
 *
 *	This routine unfreezes the atoms numbered from index to (index+n-1).
 */

#include	<md.h>

life(index,n)
	int	index,n;
	{
	register int	i;

	if	(index + n > natoms)
		{
		fprintf(stderr,"LIFE: Not that many atoms\n");
		exit(1);
		}
	for	(i = index; i < index + n; i++)
		{
		atom[i].flags &= ~A_FIXED;
		}
	}
