/*
 *	init:
 *
 *	This routine does the system independent initializations (feathering
 *	parameters, elapsed time, period, and inverse period) and then
 *	calls sysInit for the system dependent initializations.
 */

#include	<md.h>
#include	<math.h>

init(conFile)
	char	conFile[];
{
     double	w, dz, dz3, dz4, dz5, bc_factor;
     if (natoms > nsolute) 
	{
	if (pbcType[0] == 'C'){
			bc_factor = 1.; /* for cubic */
	}
	else if (pbcType[0] == 'O'){
			bc_factor = 0.75; /* for truncated octahedron */
	}
	else
		ERROR((stderr,"init: periodic boundaries not defined\n"), exit);

	/*Use the smallest cubic side for the cuttoff disatnce*/
	w = xwall;
	if (ywall < w) w = ywall;
	if (zwall < w) w = zwall;

	swr2max = bc_factor * w * w;
	swr2min = bc_factor * (w - 1.0) * (w - 1.0);

	/* This is the square of 1/2 the distance between the hexagonal faces
	 * of the pto box, or the sq of 1/2 the small length of the cube */
	
	dz = swr2max - swr2min;   /* r**2 (min max) for switching region */
	dz3 = (dz*dz*dz);
	dz4 = dz3*dz;
	dz5 = dz4*dz;
	swcoef[0] =    1.0;       /* switching parameters for group cut-offs */
	swcoef[1] =  -10.0 / dz3;
	swcoef[2] =   15.0 / dz4;
	swcoef[3] =   -6.0 / dz5;
	pswcoef[0] = -60.0 / dz3;	/* 2 times swcoef's derivs */
	pswcoef[1] = 120.0 / dz4;
	pswcoef[2] = -60.0 / dz5;

/*
 *	Set up other constants
 */
	period.fx = 2.0 * xwall;	/* box period (boundary conditions) */
	period.fy = 2.0 * ywall;	/* box period (boundary conditions) */
	period.fz = 2.0 * zwall;	/* box period (boundary conditions) */
	iperiod.fx = 1.0 / period.fx;		/* inverse period */
	iperiod.fy = 1.0 / period.fy;		/* inverse period */
	iperiod.fz = 1.0 / period.fz;		/* inverse period */
	
	pslabSize = 2.;
	nslabs = (int) (2*zwall/(pslabSize));
	if (
	   ( pTensor = (tripd *) calloc(nslabs,sizeof(tripd))) == NULL
	   )
		ERROR((stderr,"init: out of core\n"), exit);
	pslabSize = 2*zwall/nslabs;
        }
/*
 *	Call the system dependent initialization routine
 */
	sysInit(conFile);
}
