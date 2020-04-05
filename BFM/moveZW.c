#include	<md.h>
#include	<system.h>
#define 	FinalZwall	29.846323 /* Final value of zwall */
/*	This routine move the Z wall to a new position, calculate *
 *	the density of CFM, adjust the value of zwall, period.fz   *
 *	and iperiod.fz
 */
	
moveZW(done)
int *done;
{
/*	find the closest distance from a CFM molecule and zwall */
int index,i;
double minDis;
double delta, diagonal;

minDis = zwall;
for(i=0;i<nCFM*5;i = i+5){
	diagonal = abs(pos[i].fx)+ abs(pos[i].fy) + abs(pos[i].fz);
      	if ((1.5*zwall - diagonal) < minDis)
     		 minDis =1.5*zwall - diagonal;
}
delta = 0.9*minDis; /* we move the walls by this distance*/

if (zwall-delta < FinalZwall){
	delta = zwall - FinalZwall;
	xwall = ywall = zwall = FinalZwall;
	*done = 1;
}
else	{
	zwall -= delta;
	xwall = ywall = zwall;
	*done = 0;
}

	period.fz = 2.0 * zwall;	/* box period (boundary conditions) */
	iperiod.fz = 1.0 / period.fz;		/* inverse period */
	period.fy = 2.0 * ywall;	/* box period (boundary conditions) */
	iperiod.fy = 1.0 / period.fy;		/* inverse period */
	period.fx = 2.0 * xwall;	/* box period (boundary conditions) */
	iperiod.fx = 1.0 / period.fx;		/* inverse period */

fprintf(stderr,"%f    %f    \n",zwall,delta );
}
