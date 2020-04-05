/*
 *	wrsv:
 *
 *	This routine writes out the velocities for a restart trajectories.
 *	It has the same format as the .vel file (we could not used wvel
 *	since wvel does not close the file).
 */

#include	<md.h>
#include	<sys/types.h>

wrsv(file)
	char	*file;
{
	char	filetype[5];
	FILE	*fp;

/*	Get the time of day to stamp file ...				*/

	time(&datestamp);

/*	Open file and give error message if there is a problem ...	*/

	if	((fp = fopen(file,"w")) == NULL) {
		fprintf(stderr,"wrsv:  can't open %s\n",file);
		exit(1);
	}

/*	Write the pertinent stuff ...					*/

	sprintf(filetype, ".vel");
	fwrite(filetype, sizeof(char), 5, fp);
	fwrite(&datestamp, sizeof(datestamp), 1, fp);
	fwrite(tvel, sizeof tvel[0], natoms, fp);
	fclose(fp);
}
