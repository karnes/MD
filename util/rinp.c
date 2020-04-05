/*
 *	rinp:
 *
 *	This routine reads the input files (formerly fill-files) and stores
 *	the data in appropriate globals.  Space is also allocated for
 *	the arrays, if it's the first time through.  The mass arrays (mass,
 *	imass) are calculated.
 */

#include	<md.h>
#include	<sys/types.h>

rinp(file)
	char	*file;
{
	static int	allocate = 1;
	FILE		*fp;
	int		i;
	float		getmass();

/*	Open input file ...						*/

	if	((fp = fopen(file,"r")) == NULL) {
		fprintf(stderr,"rinp:  can't open %s\n",file);
		exit(1);
	}

/*	Read filetype string and warn if not ".inp" ...		*/

	fread(filetype, sizeof(char), 5, fp);
	if (strcmp(filetype, ".inp") != 0)
		fprintf(stderr, "rinp:  file is not a .inp file");

/*	Read header information ...					*/

	fread(&datestamp, sizeof(datestamp), 1, fp);
	fread(status, sizeof(char), 4, fp);	/* status = NEW, MIN, EQU */
	fread(&natoms, sizeof(natoms), 1, fp);
	fread(&nsolute, sizeof(nsolute), 1, fp);
	fread(&xwall, sizeof(xwall), 1, fp);
	fread(&ywall, sizeof(ywall), 1, fp);
	fread(&zwall, sizeof(zwall), 1, fp);
	fread(&EqTemp, sizeof(EqTemp), 1, fp);	/* if status = EQU */
	fread(&DEqTemp, sizeof(DEqTemp), 1, fp); 	/* '' */
	fread(&EqPress, sizeof(EqPress), 1, fp);
	fread(&DEqPress, sizeof(DEqPress), 1, fp);
		
/*	Allocate space for arrays, if this is the first time through ...*/

	if (allocate &&
	(  (atom	= (parts  *) calloc(natoms,sizeof(parts))) == NULL
	|| (mass	= (double *) calloc(natoms,sizeof(double))) == NULL
	|| (imass	= (double *) calloc(natoms,sizeof(double))) == NULL
	|| (force	= (tripd  *) calloc(natoms,sizeof(tripd))) == NULL
	|| (pos		= (tripd  *) calloc(natoms,sizeof(tripd))) == NULL
	|| (vel		= (tripd  *) calloc(natoms,sizeof(tripd))) == NULL
	|| (tvel	= (tripd  *) calloc(natoms,sizeof(tripd))) == NULL
	|| (itable	= (int  *) calloc(natoms,sizeof(int))) == NULL
	)) {
		ERROR((stderr,"rinp: out of core\n"), exit);
	}
/*	if (allocate){
	   for(i=0;i<natoms;i++){
		force[i].fx = force[i].fy = force[i].fz = 0.0;
		pos[i].fx = pos[i].fy = pos[i].fz = 0.0;
		tvel[i].fx = tvel[i].fy = tvel[i].fz = 0.0;
		itable[i] = 0;
		vel[i].fx = vel[i].fy = vel[i].fz = 0.0;
		mass[i] = imass[i] = 0.0;
	   }
	   fprintf(stderr,"rinp.c -- wrote zeroes to arrays.\n");
	}*/

	allocate = 0;	/* so we won't allocate space the next time */

/*	Read the atom specifications and the positions ...		*/

	fread(atom,sizeof atom[0],natoms,fp);
	fread(&xtrInQ, sizeof(xtrInQ), 1, fp);
	fread(pos, sizeof pos[0], natoms, fp);

/*	Initialize the mass arrays appropriately ...			*/

	for	(i = 0; i < natoms; i++) {
		mass[i] = (double) getmass(atom[i].type);
		imass[i] = 1.0 / mass[i];
	}

/*	All done ...							*/

	fclose(fp);
}
