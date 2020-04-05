/*
 *	addW.c:	This program converts binary input-files into ascii
 *			text files where a layer of water molecules of
 *			thickness L has been added. The initial file
 *                       is a liquid/liquid interface where the water
 *                      is in the Z<0 side and there is a water/vapor interface.
 *			the file must be converted to binary using writebin and equilibrated.
 *
 *	To call:	addW <binary_file> <ascii_file> <nw> <L> <G>
 *                     nw is the original number of water molecules
 *                     G is the approximate location of the Gibbs surface
 */

#include	<stdio.h>
#include	<typedefs.h>
#include	<globals.h>
#include	<time.h>

main (argc, argv)
	int	argc;
	char	*argv[];
{
	int	i,nw,n,k;
	double L,G,WS,Wmid, atof();
	tripd npos[2000];
	FILE	*fp;

/*
 *  argv[1] and argv[2] contain the input and output file names
 */
	if (argc != 6) {
		fprintf(stderr, "usage: %s <binary_file> <ascii_file> <nw> <L> <G>\n",
			argv[0]);
		exit(0);
	}
nw = atoi(argv[3]);
L = atof(argv[4]);
G = atof(argv[5]);
/*
 *  open the binary input file
 */
	if ((fp = fopen(argv[1], "r")) == NULL) {
		fprintf(stderr, "addW: can't open binary input-file\n");
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
		fprintf(stderr, "addW: out of memory\n");
		exit(1);
	}

	fread(atom, sizeof(parts), natoms, fp);
	fread(&xtrInQ, sizeof(xtrInQ), 1, fp);
	fread(pos, sizeof(tripd), natoms, fp);

/*
 *  Close the binary file
 */
	fclose(fp);

WS = nw/(0.0334*xwall*ywall*4);/* approximate thickness of the water slab*/
Wmid = G-WS/2;/*midpoint location of the water slab*/

n = 0;
for (i=0;i<nw;i++){
/*shift by -L all the water that are at Z < Wmid, creating a void of thickness L*/
	if (pos[i*3].fz < Wmid){
		pos[i*3].fz -= L;
		pos[i*3+1].fz -= L;
		pos[i*3+2].fz -= L;
	}
/*Copy the positions of the water molecules that are between Wmid and Wmid+L and move them to the void*/
	if (pos[i*3].fz > Wmid && pos[i*3].fz < Wmid+L){
		for (k=0;k<3;k++){
			npos[n].fx = pos[i*3+k].fx;
			npos[n].fy = pos[i*3+k].fy;
			npos[n].fz = pos[i*3+k].fz-L;
			n++;
		}
	}
}
fprintf(stderr,"%d water atoms have been added\n",n);
/*new number of atoms in the system is natoms+n*/

/*
 *  Write the output to the ascii input-file in argv[2]
 */
	if ((fp = fopen(argv[2], "w")) == NULL) {
		fprintf(stderr, "addW: can't open ascii input-file\n");
		exit(1);
	}

	fprintf(fp,"Input File -- \t%s\n",argv[1]);
	fprintf(fp,"Filetype  -- \t.inp\n");

/*  ctime() inserts a newline character so we don't need to do it */
	fprintf(fp,"Datestamp -- \t%s", ctime(&datestamp));

	fprintf(fp, "Status -- \t%s\n", status);
	fprintf(fp,"Number of atoms (natoms) -- \t%d\n", natoms+n);
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
/*print the original water*/
	for (i=0; i<nw*3; i++) {
		fprintf(fp,
"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
i,atom[i].type,pos[i].fx,pos[i].fy,pos[i].fz,atom[i].flags,atom[i].parent,atom[i].param1,atom[i].param2);
	}
/*print the additional water*/
	for (i=0; i < n; i++) {
		fprintf(fp,
"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
nw*3+i,atom[i%3].type,npos[i].fx,npos[i].fy,npos[i].fz,atom[i%3].flags,nw*3+3*(i/3),atom[i%3].param1,atom[i%3].param2);
	}
/*print the rest of the system atoms*/
	for (i=nw*3; i < natoms; i++) {
		fprintf(fp,
"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
n+i,atom[i].type,pos[i].fx,pos[i].fy,pos[i].fz,atom[i].flags,atom[i].parent+n,atom[i].param1,atom[i].param2);
	}

/*
 *  Close the ascii file, and deallocate the memory
 */
	fclose(fp);
	cfree(atom);
	cfree(pos);
}
