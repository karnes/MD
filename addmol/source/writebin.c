/*
 *	writebin.c:	This program converts ascii text input files
 *			created by readbin back into binary form so
 *			that they may be used for the dynamics code.
 *
 *	To call:  writebin <ascii_file> <binary_file>
 *
 */

#include	<stdio.h>
#include	<stdlib.h>
#include	<typedefs.h>
#include	<globals.h>
#include	<time.h>

main(argc, argv)
	int	argc;
	char	*argv[];
{
	int	anum, type, flags, parent, param1, param2;
	char	sbuf[256];
	int	i;
	FILE	*fp;

/*
 *  argv[1] contains the ascii inp-file, and argv[2] the binary name
 */
	if (argc != 3) {
		fprintf(stderr, "usage: %s <ascii_file> <binary_file>\n",
			argv[0]);
		exit(0);
	}

/*
 *  open the ascii input file and read the header input
 */
	if ((fp = fopen(argv[1], "r")) == NULL) {
		fprintf(stderr, "writebin: can't open ascii input file\n");
	}

	fgets(sbuf, 256, fp);		/* Input file name -- ignore */

	fgets(sbuf, 256, fp);		/* File type */
	sscanf(sbuf, "%*s %*s %s", filetype);
	if (strcmp(filetype, ".inp") != 0) {
		fprintf(stderr, "writebin: file not of type .inp\n");
		exit(1);
	}

	fgets(sbuf, 256, fp);		/* datestamp -- ignore, will recalc. */

	fgets(sbuf, 256, fp);		/* status -- NEW, MIN, or EQU */
	sscanf(sbuf, "%*s %*s %s", status);

	fgets(sbuf, 256, fp);		/* natoms */
	sscanf(sbuf, "%*s %*s %*s %*s %*s %d", &natoms);

	fgets(sbuf, 256, fp);		/* nsolute */
	sscanf(sbuf, "%*s %*s %*s %*s %*s %*s %d", &nsolute);

	fgets(sbuf, 256, fp);		/* xwall */
	sscanf(sbuf, "%*s %*s %*s %*s %*s %lf", &xwall);

	fgets(sbuf, 256, fp);		/* ywall */
	sscanf(sbuf, "%*s %*s %*s %*s %*s %lf", &ywall);

	fgets(sbuf, 256, fp);		/* zwall */
	sscanf(sbuf, "%*s %*s %*s %*s %*s %lf", &zwall);

	fgets(sbuf, 256, fp);		/* Equil. Temp. */
	sscanf(sbuf, "%*s %*s %*s %*s %lf", &EqTemp);

	fgets(sbuf, 256, fp);		/* St. Dev. of equil. temp. */
	sscanf(sbuf, "%*s %*s %*s %*s %*s %*s %*s %lf", &DEqTemp);

	fgets(sbuf, 256, fp);		/* Equil. Press */
	sscanf(sbuf, "%*s %*s %*s %*s %lf", &EqPress);

	fgets(sbuf, 256, fp);		/* St. Dev. of equil. press. */
	sscanf(sbuf, "%*s %*s %*s %*s %*s %*s %*s %lf", &DEqPress);

	fgets(sbuf, 256, fp);		/* Extra file flag */
	sscanf(sbuf, "%*s %*s %*s %*s %*s %d", &xtrInQ);

/*
 *  allocate space for the arrays and then read them
 */
	if ((atom = (parts *) calloc(natoms, sizeof(parts))) == NULL ||
	    (pos  = (tripd *) calloc(natoms, sizeof(tripd))) == NULL) {
		fprintf(stderr, "writebin: out of memory\n");
		exit(1);
	}

	fgets(sbuf, 256, fp);		/* header line -- ignore */
	for (i=0; i<natoms; i++) {
		fgets(sbuf, 256, fp);	/* pos and atom info */
		sscanf(sbuf, "%d %d %lf %lf %lf %o %d %d %d", &anum, &type,
			&pos[i].fx, &pos[i].fy, &pos[i].fz,
			&flags, &parent, &param1, &param2);
		atom[i].type   = type;
		atom[i].flags  = flags;
		atom[i].parent = parent;
		atom[i].param1 = param1;
		atom[i].param2 = param2;
	}

/*
 *  close the ascii file
 */
	fclose(fp);

/*
 *  Write output to binary input file
 */
	if ((fp = fopen(argv[2], "w")) == NULL) {
		fprintf(stderr, "writebin: can't open binary inp-file\n");
		exit(1);
	}

	fwrite(filetype, sizeof(char), 5, fp);
	time(&datestamp);
	fwrite(&datestamp, sizeof(datestamp), 1, fp);
	fwrite(status, sizeof(char), 4, fp);
	fwrite(&natoms, sizeof(natoms), 1, fp);
	fwrite(&nsolute, sizeof(nsolute), 1, fp);
	fwrite(&xwall, sizeof(xwall), 1, fp);
	fwrite(&ywall, sizeof(ywall), 1, fp);
	fwrite(&zwall, sizeof(zwall), 1, fp);
	fwrite(&EqTemp, sizeof(EqTemp), 1, fp);
	fwrite(&DEqTemp, sizeof(DEqTemp), 1, fp);
	fwrite(&EqPress, sizeof(EqPress), 1, fp);
	fwrite(&DEqPress, sizeof(DEqPress), 1, fp);
	fwrite(atom, sizeof(parts), natoms, fp);
	fwrite(&xtrInQ, sizeof(xtrInQ), 1, fp);
	fwrite(pos, sizeof(tripd), natoms, fp);

/*
 *  close the binary input file
 */
	fclose(fp);
}
