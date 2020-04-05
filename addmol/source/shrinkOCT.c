/*
 *	shrinkOCT.c:reads in binary OCT file
 *			shrinks system geometrically by input value
 *
 *	To call:	shrink <OCT binary> <output .as> <ratio (0-1)>
 *
 */
#include "typedefs.h"
#include "atomtypes.h"
#include "globals.h"
#include <stdio.h>
#include <time.h>
#include <math.h>
#define abs(x)          ((x)>0.?(x):-(x))
#define sq(x)	((x)*(x))
#define hsgn(x)         ((x)<0.?-.5:.5)
main (argc, argv)
	int	argc;
	char	*argv[];
{
	int	i,j,k;
	double sqrt(), atof(), r; // z1, dz, rdisk;
        char    sbuf[256];
	FILE	*fp;
	FILE	*fpo;

/*
 *  argv[1] and argv[2] contain the two input files names
 */
	if (argc != 4) {
		fprintf(stderr, "usage: %s <OCT binary> <output .as> <ratio to shrink (0-1)>\n",
			argv[0]);
		exit(0);
	}
	
	r = atof(argv[3]);

/*
 *  open the l/l binary input file
 */
	if ((fp = fopen(argv[1], "r")) == NULL) {
		fprintf(stderr, "can't open OCT binary input-file\n");
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
	    (pos  = (tripd *) calloc(natoms, sizeof(pos[0]))) == NULL)
	 {
		fprintf(stderr, "out of memory\n");
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
 *  Write the output to the ascii output-file in argv[3]
 */
	if ((fpo = fopen(argv[2], "w")) == NULL) {
		fprintf(stderr, "can't open ascii output-file %s\n",argv[2]);
		exit(1);
	}

	fprintf(fpo,"Input File -- \t%s\n",argv[1]);
	fprintf(fpo,"Filetype  -- \t.inp\n");

/*  ctime() inserts a newline character so we don't need to do it */
	fprintf(fpo,"Datestamp -- \t%s", ctime(&datestamp));

	fprintf(fpo, "Status -- \t%s\n", status);
	fprintf(fpo,"Number of atoms (natoms) -- \t%d\n", natoms);
	fprintf(fpo,"Number of solute atoms (nsolute) -- \t%d\n", 0);
	fprintf(fpo,"X box size (xwall) -- \t%f\n", xwall*r);
	fprintf(fpo,"Y box size (ywall) -- \t%f\n", ywall*r);
	fprintf(fpo,"Z box size (zwall) -- \t%f\n", zwall*r);
	fprintf(fpo, "Equilibration Temperature (EqTemp) -- \t%f\n",
		EqTemp);
	fprintf(fpo, "St. Dev. of Equil. Temp. (DEqTemp) -- \t%f\n",
		DEqTemp);
	fprintf(fpo, "Equilibration Pressure (EqTemp) -- \t%f\n",
		EqPress);
	fprintf(fpo, "St. Dev. of Equil. Press. (DEqPress) --\t%f\n",
		DEqPress);
	fprintf(fpo, "Extra file flag (xtrInQ) -- \t%d\n", xtrInQ);
	fprintf(fpo,
"   #  type  x position    y position    z position  flags parent param1 param2\n");
	/*print the atoms */
	for (j=0; j<natoms; j++) {
		fprintf(fpo,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
		j,atom[j].type,pos[j].fx*(r-0.01),pos[j].fy*(r-0.01),pos[j].fz*(r-0.01),atom[j].flags,atom[j].parent,atom[j].param1,atom[j].param2);
	}
/*
 *  Close the file, and deallocate the memory
 */
	fclose(fpo);
/*
 *  deallocate the memory
 */
	cfree(atom);
	cfree(pos);

}
