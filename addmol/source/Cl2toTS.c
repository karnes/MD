/*
 *	Cl2toTS.c:	reads in binary file where
 *			Cl2 molecule is at end.. inserts (Cl-Cl-Br)- in its place
 *
 *	To call:	cl2ts <water binary> <output .as> 
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
	double rClCl,cx,cy,cz;
	tripd Clvec;
	

/*
 *  argv[1] and argv[2] contain the two input files names
 */
	if (argc != 3) {
		fprintf(stderr, "usage: %s <GLY/Cl2 binary> <output .as>\n",
			argv[0]);
		exit(0);
	}

/*
 *  open the water binary input file
 */
	if ((fp = fopen(argv[1], "r")) == NULL) {
		fprintf(stderr, "can't open water binary input-file\n");
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

//find Cl-Cl unit vector
rClCl = sqrt(sq(pos[natoms-1].fx-pos[natoms-2].fx)+sq(pos[natoms-1].fy-pos[natoms-2].fy)+sq(pos[natoms-1].fz-pos[natoms-2].fz));

Clvec.fx = 0.7*(pos[natoms-1].fx-pos[natoms-2].fx)/rClCl;
Clvec.fy = 0.7*(pos[natoms-1].fy-pos[natoms-2].fy)/rClCl;
Clvec.fz = 0.7*(pos[natoms-1].fz-pos[natoms-2].fz)/rClCl;

cx = (pos[natoms-1].fx+pos[natoms-2].fx)/2.0;
cy = (pos[natoms-1].fy+pos[natoms-2].fy)/2.0;
cz = (pos[natoms-1].fz+pos[natoms-2].fz)/2.0;

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
	fprintf(fpo,"Number of atoms (natoms) -- \t%d\n", natoms+1);
	fprintf(fpo,"Number of solute atoms (nsolute) -- \t%d\n", 0);
	fprintf(fpo,"X box size (xwall) -- \t%f\n", xwall);
	fprintf(fpo,"Y box size (ywall) -- \t%f\n", ywall);
	fprintf(fpo,"Z box size (zwall) -- \t%f\n", zwall);
	fprintf(fpo, "Equilibration Temperature (EqTemp) -- \t%f\n",
		EqTemp);
	fprintf(fpo, "St. Dev. of Equil. Temp. (DEqTemp) -- \t%f\n",
		DEqTemp);
	fprintf(fpo, "Equilibration Pressure (EqTemp) -- \t%f\n",
		EqPress);
	fprintf(fpo, "St. Dev. of Equil. Press. (DEqPress) --\t%f\n",
		DEqPress);
	fprintf(fpo, "Extra file flag (xtrInQ) -- \t%d\n", xtrInQ);
/*print the water -- no change*/
	fprintf(fpo,
"   #  type  x position    y position    z position  flags parent param1 param2\n");
/*print all atoms except the Br */
for (i=0;i<natoms-2;i++) {/*loop over all the water molecules*/
		fprintf(fpo,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
i,atom[i].type,pos[i].fx,pos[i].fy,pos[i].fz,atom[i].flags,atom[i].parent,atom[i].param1,atom[i].param2);
}
// print the center Cl2
		fprintf(fpo,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
i,67,cx,cy,cz,atom[i].flags,atom[i].parent,2,3);
		fprintf(fpo,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
i+1,67,cx+Clvec.fx,cy+Clvec.fy,cz+Clvec.fz,0004,atom[i].parent,0,0);
		fprintf(fpo,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
i+2,94,cx-Clvec.fx,cy-Clvec.fy,cz-Clvec.fz,0004,atom[i].parent,0,0);
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
