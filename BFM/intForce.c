#include	<md.h>
#include	<system.h>

intForce(pos, force)
tripd	*pos;
tripd	*force;
{
int k;
k = natoms - nsolute;
if	(k == 0 || nsolute == 0)
	return;
else{
	fprintf(stderr,"no solute routine\n");
	exit(1);
}
}
