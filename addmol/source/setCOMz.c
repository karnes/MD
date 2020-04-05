/*
 * set z center of mass to user defined value
 */

#include	<stdio.h>
#include	<typedefs.h>
#include	<globals.h>
#include	<time.h>
#include        "atomtypes.h"

main (argc, argv)
	int	argc;
	char	*argv[];
{
	int	i;
	double atof(),Z,zcom,newZ,totMass,*mass;
	FILE	*fp;
	float getmass();

/*
 *  argv[1] and argv[2] contain the input and output file names
 */
	if (argc != 4) {
		fprintf(stderr, "usage: %s <binary_file> <ascii_file>> <new COM z>\n",
			argv[0]);
		exit(0);
	}

//N = atoi(argv[3]);
Z = atof(argv[3]);
//RFlag = atoi(argv[5]);
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
//	fprintf(stderr,"read in old file...\n");
	fprintf(stderr,"file: %s...\n",argv[1]);

/* find current center of mass */
	if ((mass = (double *) calloc(natoms,sizeof(double))) == NULL) 
	{
	   fprintf(stderr,"out of core\n");
	   exit(0);
	}
	for(i=0;i<natoms;i++)
	{
	   mass[i] = getmass(atom[i].type);
	}

	zcom = 0.0;
	totMass = 0.0;
	for(i=0;i<natoms;i++)
	{
//	   zcom += pos[i].fz*mass[i];
	   zcom += pos[i].fz*getmass(atom[i].type);
//	   totMass += mass[i];
	   totMass += getmass(atom[i].type);
	}
	zcom /= totMass;

//	fprintf(stderr,"mass O  = %5.2f, mass CHx = %5.2f, total mass = %9.2f\n",mass[0],mass[4000],totMass);
	fprintf(stderr,"...old z-COM = %5.2f, new z-COM = %5.2f\n",zcom,Z);
//	fprintf(stderr,"atom type 0 = %d, atom type 4000 = %d\n",atom[0].type,atom[4000].type);


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

	fprintf(fp,"   #  type  x position    y position    z position  flags parent param1 param2\n");
	for (i=0; i<natoms; i++) 
	{
	   newZ = pos[i].fz-zcom+Z;
	   if(newZ > zwall-5.0 || newZ < -zwall+5.0)
	   {
	      fprintf(stderr,"atom outside z PBC after CoM shift. quitting...\n");
	      exit(0);
	   }
	   fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
i,atom[i].type,pos[i].fx,pos[i].fy,newZ,atom[i].flags,atom[i].parent,atom[i].param1,atom[i].param2);
	}

/*
 *  Close the ascii file, and deallocate the memory
 */
	fclose(fp);
	cfree(atom);
	cfree(pos);
}
