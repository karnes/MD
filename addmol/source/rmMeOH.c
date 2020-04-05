/*
 *	rmMeOH.c:	remove all the CH3OH molecules between z1 and z2 (z1 < z2) leaving only nM
 *	To call:	rmMeOH <binary_file> <ascii_file> <initial nCH3OH> <z1> <z2> <nM>
 *
 */

#include	<stdio.h>
#include	<typedefs.h>
#include	<globals.h>
#include	<time.h>

main (argc, argv)
	int	argc;
	char	*argv[];
{
	int	i,j,n,nCH3OH,nM;
	FILE	*fp;
	double t,z1,z2,atof();
/*
 *  argv[1] and argv[2] contain the input and output file names
 */
	if (argc != 7) {
		fprintf(stderr, "usage: %s <binary_file> <ascii_file> <nCH3OH> <z1> <z2> <nM>\n",
			argv[0]);
		exit(0);
	}
nCH3OH = atoi(argv[3]);
z1 = atof(argv[4]);
z2 = atof(argv[5]);
if (z2 < z1){/*make sure that z1 is less than z2*/
t = z1;
z1 = z2;
z2 = t;
}
nM = atoi(argv[6]);
/*
 *  open the binary input file
 */
	if ((fp = fopen(argv[1], "r")) == NULL) {
		fprintf(stderr, "mix: can't open binary input-file\n");
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
	if (nsolute >0){
		fprintf(stderr, "mix: program work if nsolute = 0\n");
		exit(1);
	}
	if ((atom = (parts *) calloc(natoms, sizeof(atom[0]))) == NULL ||
	    (pos  = (tripd *) calloc(natoms, sizeof(pos[0]))) == NULL) {
		fprintf(stderr, "mix: out of memory\n");
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
		fprintf(stderr, "mix: can't open ascii input-file\n");
		exit(1);
	}

	fprintf(fp,"Input File -- \t%s\n",argv[1]);
	fprintf(fp,"Filetype  -- \t.inp\n");

/*  ctime() inserts a newline character so we don't need to do it */
	fprintf(fp,"Datestamp -- \t%s", ctime(&datestamp));

	fprintf(fp, "Status -- \t%s\n", status);
	fprintf(fp,"Number of atoms (natoms) -- \t%d\n", natoms-(nCH3OH-nM)*3);
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
/* loop over methanol */
n = 0;
for (i=0; i<nCH3OH*3; i += 3) {
   if (pos[i].fz < z1 || pos[i].fz > z2){/* keep this CH3OH*/
	fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",n*3,atom[i].type,pos[i].fx,pos[i].fy,pos[i].fz,atom[i].flags,n*3,atom[i].param1,atom[i].param2);
	fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",n*3+1,atom[i+1].type,pos[i+1].fx,pos[i+1].fy,pos[i+1].fz,atom[i+1].flags,n*3,atom[i+1].param1,atom[i+1].param2);
	fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",n*3+2,atom[i+2].type,pos[i+2].fx,pos[i+2].fy,pos[i+2].fz,atom[i+2].flags,n*3,atom[i+2].param1,atom[i+2].param2);
	n++;
    }
}
if (n != nM){
fprintf(stderr,"error selecting nM. nM must match the actual number of CH3OH molecules between z1 and z2 in the input file.\n");
exit(1);
}

for (i=nCH3OH*3; i<natoms; i++) {
	j = i-(nCH3OH-nM)*3;
                fprintf(fp,
"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
j,atom[i].type,pos[i].fx,pos[i].fy,pos[i].fz,atom[i].flags,atom[i].parent+j-i,atom[i].param1,atom[i].param2);
        }

/*
 *  Close the ascii file, and deallocate the memory
 */
	fclose(fp);
	cfree(atom);
	cfree(pos);
}
