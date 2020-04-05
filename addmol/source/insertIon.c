/*
 *	insertIon.c:	This program read the ion in water file, keep the
 *                      hydrated ion and insert this complex in the center
 *			of the DCE box after a given number of DCE were
 *			removed
 *
 *	To call:	readbin <h2o binary_file> <DCE binary> <ascii_file> <hydration number> <DCE to remove>
 *
 */
#include "typedefs.h"
#include "atomtypes.h"
#include "globals.h"
#include <stdio.h>
#include <sys/time.h>
#define abs(x)          ((x)>0.?(x):-(x))
#define sq(x)	((x)*(x))
#define hsgn(x)         ((x)<0.?-.5:.5)
#define SIZE 1500
main (argc, argv)
	int	argc;
	char	*argv[];
{
	int	i,j,k,nw,nh;
	double sqrt(), atof(), rshel, rdce, r;
	parts s,Natom[SIZE];
	tripd spos;
	tripd npos[SIZE];
	FILE	*fp;

/*
 *  argv[1] and argv[2] contain the two input files names
 */
	if (argc != 7) {
		fprintf(stderr, "usage: %s <H2O binary_file> <DCE binary> <ascii_file> <nh> <rshel> <rdce>\n",
			argv[0]);
		exit(0);
	}
nh = atoi(argv[4]);
rshel = atof(argv[5]);
rdce = atof(argv[6]);
/*
 *  open the first binary input file
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
	    (pos  = (tripd *) calloc(natoms, sizeof(pos[0]))) == NULL) {
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
if (nsolute != 1){
	fprintf(stderr,"nsolute = %d must be equal to 1\n",nsolute);
	exit(1);
}
offset(natoms-1);/*center the box on the ion*/
nw = natoms-nsolute;
while (j < nh*3){/* we want to keep exactly nh water molecules*/
j = 0;
for (i=0;i<nw;i = i+3){
    r = sqrt(sq(pos[i].fx-pos[natoms-1].fx) + sq(pos[i].fy-pos[natoms-1].fy)+
       sq(pos[i].fz-pos[natoms-1].fz));
    if (r < rshel && j<nh*3){
	npos[j].fx = pos[i].fx;npos[j].fy = pos[i].fy;npos[j].fz = pos[i].fz;
	Natom[j].type = atom[i].type;
	Natom[j].flags = atom[i].flags | A_FIXED;
	Natom[j].parent = j;
	Natom[j].param1 = atom[i].param1;
	Natom[j].param2 = atom[i].param2;
	j++;
	npos[j].fx = pos[i+1].fx;npos[j].fy = pos[i+1].fy;npos[j].fz = pos[i+1].fz;
	Natom[j].type = atom[i+1].type;
	Natom[j].flags = atom[i+1].flags | A_FIXED;
	Natom[j].parent = j-1;
	Natom[j].param1 = atom[i+1].param1;
	Natom[j].param2 = atom[i+1].param2;
	j++;
	npos[j].fx = pos[i+2].fx;npos[j].fy = pos[i+2].fy;npos[j].fz = pos[i+2].fz;
	Natom[j].type = atom[i+2].type;
	Natom[j].flags = atom[i+2].flags | A_FIXED;
	Natom[j].parent = j-2;
	Natom[j].param1 = atom[i+2].param1;
	Natom[j].param2 = atom[i+2].param2;
	j++;
   }
}	
rshel *= 1.02;
}
spos.fx = spos.fy = spos.fz = 0;
s.type = atom[natoms-1].type;
s.flags = atom[natoms-1].flags | A_FIXED;
s.param1 = atom[natoms-1].param1;
s.param2 = atom[natoms-1].param2;
nh = j/3;
printf("hydration number = %d  ",nh);
	cfree(atom);
	cfree(pos);
/*
 *  open the second binary input file
 */
	if ((fp = fopen(argv[2], "r")) == NULL) {
		fprintf(stderr, "can't open DCE binary input-file\n");
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
if(natoms+nh*3+1 > SIZE){
	fprintf(stderr,"SIZE too small. Increase to %d\n",natoms+nh*3+1);
	exit(1);
}
for (i=0;i<natoms;i = i+4){
    r = sqrt(sq(pos[i].fx) + sq(pos[i].fy)+ sq(pos[i].fz));
    if (r > rdce){
	for (k=0;k<4;k++){
		npos[j+k].fx = pos[i+k].fx;
		npos[j+k].fy = pos[i+k].fy;
		npos[j+k].fz = pos[i+k].fz;
		Natom[j+k].type = atom[i+k].type;
		Natom[j+k].flags = atom[i+k].flags;
		Natom[j+k].parent = j;
		if (k == 3 ) Natom[j+k].parent = j+2;
		Natom[j+k].param1 = atom[i+k].param1;
		Natom[j+k].param2 = atom[i+k].param2;
	}
	j = j +4;
    }
}
printf("DCE removed = %d\n",(natoms-(j-3*nh))/4);
nsolute = 1;
/*
 *  Write the output to the ascii input-file in argv[3]
 */
	if ((fp = fopen(argv[3], "w")) == NULL) {
		fprintf(stderr, "readbin: can't open ascii input-file\n");
		exit(1);
	}

	fprintf(fp,"Input File -- \t%s\n",argv[1]);
	fprintf(fp,"Filetype  -- \t.inp\n");

/*  ctime() inserts a newline character so we don't need to do it */
	fprintf(fp,"Datestamp -- \t%s", ctime(&datestamp));

	fprintf(fp, "Status -- \t%s\n", status);
	fprintf(fp,"Number of atoms (natoms) -- \t%d\n", j+1);
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
	for (i=0; i<j; i++) {
		fprintf(fp,
"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
i,Natom[i].type,npos[i].fx,npos[i].fy,npos[i].fz,Natom[i].flags,Natom[i].parent,Natom[i].param1,Natom[i].param2);
	}
		fprintf(fp,
"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
j,s.type,spos.fx,spos.fy,spos.fz,s.flags,j,s.param1,s.param2);

/*
 *  Close the ascii file, and deallocate the memory
 */
	fclose(fp);
	cfree(atom);
	cfree(pos);
}
/*
 *	offset:
 *	
 *	This routine offsets all atoms by an appropriate amount so as to bring
 *	atom n to the center of the box (0,0,0).
 */

offset(n)
	int	n;	/* atom number from which to offset positions */
{
	int	i, j;
	tripd	image, off;
	if (n > natoms-1 || n < 0) return;
	off.fx = pos[n].fx;
	off.fy = pos[n].fy;
	off.fz = pos[n].fz;

	for (i = 0; i < natoms; i++) {
		pos[i].fx -= off.fx;
		pos[i].fy -= off.fy;
		pos[i].fz -= off.fz;

		if ((atom[i].flags & A_MAJOR) && atom[i].param1) {
			image.fx = -pos[i].fx;
			image.fy = -pos[i].fy;
			image.fz = -pos[i].fz;
			mvimage(&pos[i]);
			image.fx += pos[i].fx;
			image.fy += pos[i].fy;
			image.fz += pos[i].fz;
			for (j=i+1; j <= i + atom[i].param1; j++) {
				pos[j].fx += image.fx;
				pos[j].fy += image.fy;
				pos[j].fz += image.fz;
			}
		}
		else if (atom[i].flags & A_MAJOR)
			mvimage(&pos[i]);
	}
}

/*
 *	This routine does the periodic imaging of any given
 *	vector for the parallelopiped or the truncated octahedral
 *      boundary conditions.
 */

mvimage(del)
	tripdouble	*del;
{
	double	fx, fy, fz;
	tripd period, iperiod;
	if (natoms == nsolute) /* no solvent */
	     return;
	period.fx = 2.0 * xwall;	/* box period (boundary conditions) */
	period.fy = 2.0 * ywall;	/* box period (boundary conditions) */
	period.fz = 2.0 * zwall;	/* box period (boundary conditions) */
	iperiod.fx = 1.0 / period.fx;		/* inverse period */
	iperiod.fy = 1.0 / period.fy;		/* inverse period */
	iperiod.fz = 1.0 / period.fz;		/* inverse period */
	fx = del->fx * iperiod.fx;  /* Normalizing the positions */
	fy = del->fy * iperiod.fy;  /* |fx|, |fy|, |fz| <= 0.5 for the  */
	fz = del->fz * iperiod.fz;  /* atom to be in the box */
	fx -= (int) fx;
	fy -= (int) fy;   /* bringing the particle to the nearest box */
	fz -= (int) fz;
	fx -= (int) (2 * fx);
	fy -= (int) (2 * fy);  /* bringing the particle to the box */
	fz -= (int) (2 * fz);
	/* truncated octahedral boundaries */
		if	((abs(fx) + abs(fy) + abs(fz)) > 0.75) {
			fx -= hsgn(fx);
			fy -= hsgn(fy);
			fz -= hsgn(fz);
		}

	del->fx = fx * period.fx;
	del->fy = fy * period.fy;
	del->fz = fz * period.fz;
}

/*
 *	Set up other constants
 */
