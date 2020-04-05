/*
 *	This inserts a water into another configuration file
 *	at ~0,0,0 
 *
 *	To call:	insertH2O <binary_file> <ascii_file>
 *
 */

#include "typedefs.h"
#include "atomtypes.h"
#include "globals.h"
#include <stdio.h>
#include <time.h>
#define abs(x)          ((x)>0.?(x):-(x))
#define sq(x)   ((x)*(x))
#define hsgn(x)         ((x)<0.?-.5:.5)

main (argc, argv)
	int	argc;
	char	*argv[];
{
	int	i,k,j,m,nCF,clst[260], atype, aflag, pa, p1, p2;
	double r2[260], rC, sqrt();
	tripd posO, posH[2];
	FILE	*fp;
	posO.fx = posO.fy = posO.fz = 0.0;
	posH[0].fx = 0.38205;
	posH[0].fy = -0.30234;
	posH[0].fz = -0.889344;
	posH[1].fx = -0.101586;
	posH[1].fy = 1.019924;
	posH[1].fz = -0.153315;

/*
 *  argv[1] and argv[2] contain the input and output file names
 */
	if (argc != 3) {
		fprintf(stderr, "usage: %s <binary_file> <ascii_file>\n",
			argv[0]);
		exit(0);
	}

/*
 *  open the binary input file which contains one sn2 system and chloroform molecules
 */
	if ((fp = fopen(argv[1], "r")) == NULL) {
		fprintf(stderr, "can't open binary input-file\n");
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

fprintf(stderr,"read binary...\n");
fprintf(stderr,"natoms = %d\n",natoms);
/*
 *  Write the output to the ascii input-file in argv[2].
 *  Note that the number of atoms is 3 more than the original file
 */
	if ((fp = fopen(argv[2], "w")) == NULL) {
		fprintf(stderr, "readbin: can't open ascii input-file\n");
		exit(1);
	}

fprintf(stderr,"opened file...\n");
	fprintf(fp,"Input File -- \t%s\n",argv[1]);
	fprintf(fp,"Filetype  -- \t.inp\n");

fprintf(stderr,"printed 2 lines...\n");
/*  ctime() inserts a newline character so we don't need to do it */
	fprintf(fp,"Datestamp -- \t%s", ctime(&datestamp));
fprintf(stderr,"printed time/date...\n");

	fprintf(fp, "Status -- \t%s\n", status);
	fprintf(fp,"Number of atoms (natoms) -- \t%d\n", natoms+3);
	fprintf(fp,"Number of solute atoms (nsolute) -- \t%d\n", nsolute);
	fprintf(fp,"X box size (xwall) -- \t%f\n", xwall);
	fprintf(fp,"Y box size (ywall) -- \t%f\n", ywall);
	fprintf(fp,"Z box size (zwall) -- \t%f\n", zwall);
fprintf(stderr,"printed to walls...\n");
	fprintf(fp, "Equilibration Temperature (EqTemp) -- \t%f\n",
		EqTemp);
	fprintf(fp, "St. Dev. of Equil. Temp. (DEqTemp) -- \t%f\n",
		DEqTemp);
	fprintf(fp, "Equilibration Pressure (EqTemp) -- \t%f\n",
		EqPress);
	fprintf(fp, "St. Dev. of Equil. Press. (DEqPress) --\t%f\n",
		DEqPress);
	fprintf(fp, "Extra file flag (xtrInQ) -- \t%d\n", xtrInQ);

	fprintf(fp,"   #  type  x position    y position    z position  flags parent param1 param2\n");
/*print the water oxygen*/
fprintf(stderr,"printed header...\n");
	i =0; atype = 0 ; aflag = 0026; pa = 0; p1 = 2; p2 = 3;
	fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",i,atype,posO.fx,posO.fy,posO.fz,aflag,pa,p1,p2);
/*print the water hydrogen*/
	i++; atype = 1 ; aflag = 0004; pa = 0; p1 = 0; p2 = 0;
	fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",i,atype,posH[0].fx,posH[0].fy,posH[0].fz,aflag,pa,p1,p2);
	i++;
	fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",i,atype,posH[1].fx,posH[1].fy,posH[1].fz,aflag,pa,p1,p2);
fprintf(stderr,"printed water...\n");
/*print the rest of the atoms*/
	for (j=0; j<natoms; j++) {
	    fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",j+3,atom[j].type,pos[j].fx,pos[j].fy,pos[j].fz,atom[j].flags,3+atom[j].parent,atom[j].param1,atom[j].param2);
	}
/*
 *  Close the ascii file, and deallocate the memory
 */
	fclose(fp);
	free(atom);
	free(pos);
}
