/*
 *	rvel:
 *
 *	This routine reads one time step from a velocity file (.vel or .rsv),
 *      assuming that space allocation has been made for the velocities arays.
 *	The calling program must know how many time steps this file contain.
 */

#include	<md.h>
#include	<sys/types.h>

rvel(file,initQ)
char	*file;
int	initQ;
{
FILE		*fp;
if	(initQ) {
/*	Open input file ...						*/

	if	((fp = fopen(file,"r")) == NULL) {
		fprintf(stderr,"rvel:  can't open %s\n",file);
		exit(1);
	}

/*	Read filetype string and warn if not ".vel" ...		*/

	fread(filetype, sizeof(char), 5, fp);
	if (strcmp(filetype, ".vel") != 0)
		fprintf(stderr, "rvel:  file is not a .vel file");

/*	Read header information ...					*/

	fread(&datestamp, sizeof(datestamp), 1, fp);
}
/*	Read one time step				*/
fread(vel, sizeof vel[0], natoms, fp);
}
