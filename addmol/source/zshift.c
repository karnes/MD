/*
 * move all atom along Z by Z and change the running number by N
 * if RFlag = -1, first multiply all zpos by -1 (reflection)
 */

#include	<stdio.h>
#include	<typedefs.h>
#include	<globals.h>
#include	<time.h>

main (argc, argv)
	int	argc;
	char	*argv[];
{
	int	i, N, RFlag;
	double atof(),Z;
	FILE	*fp;

/*
 *  argv[1] and argv[2] contain the input and output file names
 */
	if (argc != 6) {
		fprintf(stderr, "usage: %s <binary_file> <ascii_file> <N> <Z> <RFlag>\n",
			argv[0]);
		exit(0);
	}

N = atoi(argv[3]);
Z = atof(argv[4]);
RFlag = atoi(argv[5]);
/*
 *  open the binary input file
 */
	if ((fp = fopen(argv[1], "r")) == NULL) {
		fprintf(stderr, "readbin: can't open binary input-file\n");
		exit(1);
	}

/*
 *  Read the header data, allocate space for the arrays, and read them too.
 */
	fread(filetype, sizeof(char), 5, fp);
	fread(&datestamp, sizeof(datestamp), 1, fp);
	fread(status, sizeof(char), 4, fp);
	fread(&natoms, sizeof(natoms), 1, fp);
	fread(&nsolute, sizeof(nsolute), 1, fp);
	fread(&xwall, sizeof(xwall), 1, fp);
	fread(&ywall, sizeof(ywall), 1, fp);
	fread(&zwall, sizeof(zwall), 1, fp);
	fread(&EqTemp, sizeof(EqTemp), 1, fp);
	fread(&DEqTemp, sizeof(DEqTemp), 1, fp);
	fread(&EqPress, sizeof(EqPress), 1, fp);
	fread(&DEqPress, sizeof(DEqPress), 1, fp);

	if ((atom = (parts *) calloc(natoms, sizeof(atom[0]))) == NULL ||
	    (pos  = (tripd *) calloc(natoms, sizeof(pos[0]))) == NULL) {
		fprintf(stderr, "readbin: out of memory\n");
		exit(1);
	}

	fread(atom, sizeof(parts), natoms, fp);
	fread(&xtrInQ, sizeof(xtrInQ), 1, fp);
	fread(pos, sizeof(tripd), natoms, fp);

/*
 *  Close the binary file
 */
	fclose(fp);

/*
 *  Write the output to the ascii input-file in argv[2]
 */
	if ((fp = fopen(argv[2], "w")) == NULL) {
		fprintf(stderr, "readbin: can't open ascii input-file\n");
		exit(1);
	}

	fprintf(fp,"Input File -- \t%s\n",argv[1]);
	fprintf(fp,"Filetype  -- \t.inp\n");

/*  ctime() inserts a newline character so we don't need to do it */
	fprintf(fp,"Datestamp -- \t%s", ctime(&datestamp));

	fprintf(fp, "Status -- \t%s\n", status);
	fprintf(fp,"Number of atoms (natoms) -- \t%d\n", natoms);
	fprintf(fp,"Number of solute atoms (nsolute) -- \t%d\n", nsolute);
	fprintf(fp,"X box size (xwall) -- \t%f\n", xwall);
	fprintf(fp,"Y box size (ywall) -- \t%f\n", ywall);
	fprintf(fp,"Z box size (zwall) -- \t%f\n", zwall);
	fprintf(fp, "Equilibration Temperature (EqTemp) -- \t%f\n",
		EqTemp);
	fprintf(fp, "St. Dev. of Equil. Temp. (DEqTemp) -- \t%f\n",
		DEqTemp);
	fprintf(fp, "Equilibration Pressure (EqTemp) -- \t%f\n",
		EqPress);
	fprintf(fp, "St. Dev. of Equil. Press. (DEqPress) --\t%f\n",
		DEqPress);
	fprintf(fp, "Extra file flag (xtrInQ) -- \t%d\n", xtrInQ);

	fprintf(fp,
"   #  type  x position    y position    z position  flags parent param1 param2\n");
for (i=0; i<natoms; i++) {
	fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
i+N,atom[i].type,pos[i].fx,pos[i].fy,RFlag*pos[i].fz+Z,atom[i].flags,atom[i].parent+N,atom[i].param1,atom[i].param2);
	}

/*
 *  Close the ascii file, and deallocate the memory
 */
	fclose(fp);
	cfree(atom);
	cfree(pos);
}
