#include	<stdio.h>

/*
 *	This routine loops over randomizations to equilibrate the system.
 */

equil(numRand, numSteps, kelvin, fixTQ, prStatQ, vScaleQ)
	int	numRand;	/* # of equilibration steps/randomizations */
	int	numSteps;	/* # of timesteps per randomization */
	double	kelvin;		/* temperature in Kelvin */
	int	fixTQ;		/* !0 == use constant temperature algorithm */
	int	prStatQ;	/* !0 == print status line */
	int	vScaleQ;	/* !0 == scale initial velocity */
{
//printf("enter equil \n");
	int	i;

	for	(i = 0; i < numRand; i++) {
		ranvel(kelvin, 1, vScaleQ,1);
		if	(prStatQ) {
			integ(1, fixTQ, 0);
			wstat(stderr);
			integ(numSteps-1, fixTQ, 0);
			wstat(stderr);
		}
		else
			integ(numSteps, fixTQ, 0);
	}
//printf("exit equil\n");
}

