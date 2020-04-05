/*
 *	winp:
 *
 *	This routine writes out the appropriate globals to a file with the
 *	same format as an input file.  Typically for use after minimization
 *	or equilibration when you want to write a new input file with the
 *	minimized/equilibrated atom positions.
 */

#include	<md.h>
#include	<sys/types.h>

winp(file)
	char	*file;
{
	char	filetype[5];
	FILE	*fp;

/*	Get the time of day to stamp file ...				*/

	time(&datestamp);

/*	Open file and give error message if there is a problem ...	*/

	if	((fp = fopen(file,"w")) == NULL) {
		fprintf(stderr,"winp:  can't open %s\n",file);
		exit(1);
	}

/*	Write the pertinent stuff ...					*/

	sprintf(filetype, ".inp");
	fwrite(filetype, sizeof(char), 5, fp);
	fwrite(&datestamp, sizeof(datestamp), 1, fp);
	fwrite(status, sizeof(char), 4, fp);	/* must define status */
	fwrite(&natoms, sizeof(natoms), 1, fp);
	fwrite(&nsolute, sizeof(nsolute), 1, fp);
	fwrite(&xwall, sizeof(xwall), 1, fp);
	fwrite(&ywall, sizeof(ywall), 1, fp);
	fwrite(&zwall, sizeof(zwall), 1, fp);
	fwrite(&EqTemp, sizeof(EqTemp), 1, fp);
	fwrite(&DEqTemp, sizeof(DEqTemp), 1, fp);
	fwrite(&EqPress, sizeof(EqPress), 1, fp);
	fwrite(&DEqPress, sizeof(DEqPress), 1, fp);
	fwrite(atom,sizeof atom[0],natoms,fp);
	fwrite(&xtrInQ, sizeof(xtrInQ), 1, fp);
	fwrite(pos, sizeof pos[0], natoms, fp);

/*	All done ...							*/

	fclose(fp);
}
