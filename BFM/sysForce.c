#include	<md.h>
#include	<system.h>
 
sysForce(pos, force)
	tripd	*pos;
	tripd	*force;
{
if (nsolute == 0) return;
else{
	fprintf(stderr,"nsolute = %d . Wrong number of solute atoms\n",nsolute);
	exit(1);
}
}
