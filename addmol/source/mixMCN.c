/*
 *	mixMCN.c:	Replace all the CH3OH molecules between z1 and z2
 *		by the same  number of CH3CN molecules. 
 *	To call:	mixMCN <binary_file> <ascii_file> <initial nCH3OH> <z1> <z2>
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
	int	i,j,nCH3OH,nM,nA;
	int C_type, N_type, M_type;
	FILE	*fp;
	double z1,z2,atof();
	tripd pA[3], posA[3000];
/*
 *  argv[1] and argv[2] contain the input and output file names
 */
	if (argc != 6) {
		fprintf(stderr, "usage: %s <binary_file> <ascii_file> <nCH3OH> <z1> <z2>\n",
			argv[0]);
		exit(0);
	}
nCH3OH = atoi(argv[3]);
z1 = atof(argv[4]);
z2 = atof(argv[5]);
C_type = 18;
N_type = 49;
M_type = 10;
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
	fprintf(fp,"Number of atoms (natoms) -- \t%d\n", natoms);/*no change in th enumber of atoms*/
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
nM = 0;
nA = 0;
for (i=0; i<nCH3OH*3; i += 3) {
   if (pos[i].fz > z1 && pos[i].fz < z2){/* replace this CH3OH*/
	selectCH3CN(&pos[i],pA);
	for (j = 0; j<3; j++){
		posA[nA*3+j].fx = pA[j].fx;
		posA[nA*3+j].fy = pA[j].fy;
		posA[nA*3+j].fz = pA[j].fz;
	}
	nA++;
    }
    else {/* print this CH3OH*/
	fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",nM*3,atom[i].type,pos[i].fx,pos[i].fy,pos[i].fz,atom[i].flags,nM*3,atom[i].param1,atom[i].param2);
	fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",nM*3+1,atom[i+1].type,pos[i+1].fx,pos[i+1].fy,pos[i+1].fz,atom[i+1].flags,nM*3,atom[i+1].param1,atom[i+1].param2);
	fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",nM*3+2,atom[i+2].type,pos[i+2].fx,pos[i+2].fy,pos[i+2].fz,atom[i+2].flags,nM*3,atom[i+2].param1,atom[i+2].param2);
	nM++;
    }
}
printf("nA = %d, nM = %d\n",nA, nM);
for (i=0; i<nA; i++) {
	fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",nM*3+i*3,C_type,posA[i*3].fx,posA[i*3].fy,posA[i*3].fz,atom[i*3].flags,nM*3+i*3,atom[i*3].param1,atom[i*3].param2);
	fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",nM*3+i*3+1,N_type,posA[i*3+1].fx,posA[i*3+1].fy,posA[i*3+1].fz,atom[i*3+1].flags,nM*3+i*3,atom[i*3+1].param1,atom[i*3+1].param2);
	fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",nM*3+i*3+2,M_type,posA[i*3+2].fx,posA[i*3+2].fy,posA[i*3+2].fz,atom[i*3+2].flags,nM*3+i*3,atom[i*3+2].param1,atom[i*3+2].param2);
}

for (i=nCH3OH*3; i<natoms; i++) {
                fprintf(fp,
"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
i,atom[i].type,pos[i].fx,pos[i].fy,pos[i].fz,atom[i].flags,atom[i].parent,atom[i].param1,atom[i].param2);
        }

/*
 *  Close the ascii file, and deallocate the memory
 */
	fclose(fp);
	cfree(atom);
	cfree(pos);
}
/*need to change the routine below*/
#define rMC 1.46
#define rCN 1.17 
selectCH3CN(pm,pa)
tripd *pm;
tripd *pa;
{
tripd v1,v2;
double b;
pa[0].fx = pm[0].fx; pa[0].fy = pm[0].fy; pa[0].fz = pm[0].fz; /*C takes O's position*/
/* v1 and v2 are unit vectors along the two bonds in CH3OH*/
v1.fx = pm[1].fx - pm[0].fx; v1.fy = pm[1].fy - pm[0].fy; v1.fz = pm[1].fz - pm[0].fz; 
v2.fx = pm[2].fx - pm[0].fx; v2.fy = pm[2].fy - pm[0].fy; v2.fz = pm[2].fz - pm[0].fz;
b = sqrt (v1.fx*v1.fx + v1.fy*v1.fy + v1.fz*v1.fz);
v1.fx /= b; v1.fy /= b; v1.fz /= b;
b = sqrt (v2.fx*v2.fx + v2.fy*v2.fy + v2.fz*v2.fz);
v2.fx /= b; v2.fy /= b; v2.fz /= b;
/* The N and the CH3 groups are placed along the two old bonds at the old bond length*/
pa[1].fx = pa[0].fx + v1.fx *rCN;
pa[1].fy = pa[0].fy + v1.fy *rCN;
pa[1].fz = pa[0].fz + v1.fz *rCN;
pa[2].fx = pa[0].fx + v2.fx *rMC;
pa[2].fy = pa[0].fy + v2.fy *rMC;
pa[2].fz = pa[0].fz + v2.fz *rMC;
}

