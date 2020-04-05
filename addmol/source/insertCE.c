/*
 *	insertCE.c:	This program read the xyz positions of the coumarine
 *                      and insert this in the center of an ethanol OCT box
 *                      after all the ethanol molecules inside a disk centered
 *                      at 0,0,0 with radius rdisk and half thickness dz
*                       normal to the z axis are removed.
 *
 *	To call:	insCE <ETH binary> <Coumarine xyz> <ascii_file> <rdisk> <dz>
 *
 */
#include "typedefs.h"
#include "atomtypes.h"
#include "globals.h"
#include <stdio.h>
#include <time.h>
#define abs(x)          ((x)>0.?(x):-(x))
#define sq(x)	((x)*(x))
#define hsgn(x)         ((x)<0.?-.5:.5)
#define NSOL 42
main (argc, argv)
	int	argc;
	char	*argv[];
{
	int	i,j,k,l,nE,nER, inC[1100];
	int aflag,  pa1, pa2, ctype[NSOL];
	double sqrt(), atof(), r, z1, dz, rdisk;
	tripd pc[NSOL];
        char    sbuf[256];
	FILE	*fp;
	FILE	*fpo;

/*
 *  argv[1] and argv[2] contain the two input files names
 */
	if (argc != 6) {
		fprintf(stderr, "usage: %s <ETH binary> <coumarine xyz_file> <ascii_file> <rdisk> <dz>\n",
			argv[0]);
		exit(0);
	}

rdisk = atof(argv[4]);/* cavity radius*/
dz = atof(argv[5]);/* half thickness of the cylindrecial cavity)*/
/*
 *  open the ethanol binary input file
 */
	if ((fp = fopen(argv[1], "r")) == NULL) {
		fprintf(stderr, "can't open ethanol binary input-file\n");
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
if (nsolute != 0){
	fprintf(stderr,"nsolute = %d must be equal to 0\n",nsolute);
	exit(1);
}
nE = natoms/4;/* number of ethanol molecules in the original file*/
nER = 0;
for (i=0;i<nE;i = i+1){
	inC[i] = 0;/* assume this ethanol is outside the cylinderical cavity*/
	r = sqrt(sq(pos[i*4].fx) + sq(pos[i*4].fy));
	z1 = sqrt(sq(pos[i*4].fz));
	if (r < rdisk && z1 < dz){
		 inC[i] = 1;/* this ethanol will be removed*/
		nER++;
	}
}	

fprintf(stderr," %d ethanol molecules will be removed\n",nER);

natoms = natoms - nER*4  + NSOL;
/*
 *  Write the output to the ascii output-file in argv[3]
 */
	if ((fpo = fopen(argv[3], "w")) == NULL) {
		fprintf(stderr, "can't open ascii output-file %s\n",argv[3]);
		exit(1);
	}

	fprintf(fpo,"Input File -- \t%s\n",argv[1]);
	fprintf(fpo,"Filetype  -- \t.inp\n");

/*  ctime() inserts a newline character so we don't need to do it */
	fprintf(fpo,"Datestamp -- \t%s", ctime(&datestamp));

	fprintf(fpo, "Status -- \t%s\n", status);
	fprintf(fpo,"Number of atoms (natoms) -- \t%d\n", natoms);
	fprintf(fpo,"Number of solute atoms (nsolute) -- \t%d\n", NSOL);
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
/*print the ETH except those that were removed*/
i = 0;/*index run over atoms of ETH that are kept*/
for (j=0; j<nE; j++) {/*loop over all the ETH molecules*/
	if (inC[j] == 0){/*if this ETH is outside, print it*/
		l = 4*j;
		fprintf(fpo,
"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
i,atom[l].type,pos[l].fx,pos[l].fy,pos[l].fz,atom[l].flags,i,atom[l].param1,atom[l].param2);
i++;
		fprintf(fpo,
"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
i,atom[l+1].type,pos[l+1].fx,pos[l+1].fy,pos[l+1].fz,atom[l+1].flags,i-1,atom[l+1].param1,atom[l+1].param2);
i++;
		fprintf(fpo,
"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
i,atom[l+2].type,pos[l+2].fx,pos[l+2].fy,pos[l+2].fz,atom[l+2].flags,i-2,atom[l+2].param1,atom[l+2].param2);
i++;
		fprintf(fpo,
"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
i,atom[l+3].type,pos[l+3].fx,pos[l+3].fy,pos[l+3].fz,atom[l+3].flags,i-1,atom[l+3].param1,atom[l+3].param2);
i++;
	}
}
/*
 *  deallocate the memory
 */
	cfree(atom);
	cfree(pos);

/*
 *  open the Coumarine xyz input file
 */
	if ((fp = fopen(argv[2], "r")) == NULL) {
		fprintf(stderr, "can't open coumarine xyz input-file\n");
		exit(1);
	}
for (i=0;i<NSOL;i++){
        fgets(sbuf,256,fp);
        if (sscanf(sbuf,"%lf%lf%lf%d",&pc[i].fx, &pc[i].fy, &pc[i].fz, &ctype[i]) != 4) {
                fprintf(stderr, "error reading line %d coumarine file\n",i);
                exit(1);
        }
}



/*
 *  Close the coumarine file
 */
	fclose(fp);
j=(nE-nER)*4;
for (i=0; i<NSOL; i++) {
	if (i == 0){
	aflag = 0026;
	pa1 = NSOL-1;
	pa2 = NSOL;
	fprintf(fpo,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
	j+i,ctype[i],pc[i].fx,pc[i].fy,pc[i].fz,aflag,j,pa1,pa2);
	}
	if (i > 0){
	aflag = 0004;
	pa1 = 0;
	pa2 = 0;
	fprintf(fpo,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
	j+i,ctype[i],pc[i].fx,pc[i].fy,pc[i].fz,aflag,j,pa1,pa2);
	}
}
/*
 *  Close the file, and deallocate the memory
 */
	fclose(fpo);
}
