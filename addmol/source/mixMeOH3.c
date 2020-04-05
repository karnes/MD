/*
 *	mixMeOH3.c:	Replace a number of CH3CN molecules between z1 and z2
 *		by the same  number of CH3OH molecules. Limits replaced MeOH to "center" of box.
 *	To call:	mixMeOH <binary_file> <ascii_file> <nCH3OH> <initial nCH3CN> <z1> <z2>
 *
 */

#include	<stdio.h>
#include	<typedefs.h>
#include	<globals.h>
#include	<time.h>
#include        <math.h>
main (argc, argv)
	int	argc;
	char	*argv[];
{
	int	i,k,nCH3OH,nCH3CN,actualM,nA,list[1023];
	int O_type, H_type, M_type;
	FILE	*fp;
	double z1,z2,atof();
	tripd pM[3];
/*
 *  argv[1] and argv[2] contain the input and output file names
 */
	if (argc != 7) {
		fprintf(stderr, "usage: %s <binary_file> <ascii_file> <nCH3OH> <nCH3CN> <z1> <z2>\n",
			argv[0]);
		exit(0);
	}
nCH3OH = atoi(argv[3]);
nCH3CN = atoi(argv[4]);
z1 = atof(argv[5]);
z2 = atof(argv[6]);
O_type = 45;
H_type = 7;
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
/* loop over acetonitrile */
nA = 0;
actualM = 0;
for (i=0; i<nCH3CN*3; i += 3) {
   if (pos[i].fz > z1 && pos[i].fz < z2 && (pos[i].fx*pos[i].fx+pos[i].fy*pos[i].fy) < 25.0 && actualM < nCH3OH){/* replace this CH3CN*/
	selectCH3OH(&pos[i],pM);
	fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",actualM*3,O_type,pM[0].fx,pM[0].fy,pM[0].fz,atom[i].flags,actualM*3,atom[i].param1,atom[i].param2);
	fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",actualM*3+1,H_type,pM[1].fx,pM[1].fy,pM[1].fz,atom[i+1].flags,actualM*3,atom[i+1].param1,atom[i+1].param2);
	fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",actualM*3+2,M_type,pM[2].fx,pM[2].fy,pM[2].fz,atom[i+2].flags,actualM*3,atom[i+2].param1,atom[i+2].param2);
	actualM++;
    }
    else {/* keep this CH3CN*/
	list[nA] = i;
	nA++;
    }
}
if (actualM < nCH3OH){
	fprintf(stderr,"incomplete replacemnt of CH3CN. nCH3OH = %d\n",actualM);
	fprintf(stderr," z1 = %f, z2 = %f, must increase z2-z1\n",z1,z2);
}
printf("nA = %d\n",nA);
for (i=0; i<nA; i++) {
	k = list[i];
	fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",actualM*3+i*3,atom[k].type,pos[k].fx,pos[k].fy,pos[k].fz,atom[k].flags,actualM*3+i*3,atom[k].param1,atom[k].param2);
	fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",actualM*3+i*3+1,atom[k+1].type,pos[k+1].fx,pos[k+1].fy,pos[k+1].fz,atom[k+1].flags,actualM*3+i*3,atom[k+1].param1,atom[k+1].param2);
	fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",actualM*3+i*3+2,atom[k+2].type,pos[k+2].fx,pos[k+2].fy,pos[k+2].fz,atom[k+2].flags,actualM*3+i*3,atom[k+2].param1,atom[k+2].param2);
}

for (i=nCH3CN*3; i<natoms; i++) {
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
#define rMO 1.43
#define rOH 0.945
selectCH3OH(pA,pm)
tripd *pA;
tripd *pm;
{
tripd v1,v2;
double b;
pm[0].fx = pA[0].fx; pm[0].fy = pA[0].fy; pm[0].fz = pA[0].fz; /*O takes C's position*/
/* v1 and v2 are unit vectors along the two bonds in CH3CN*/
v1.fx = pA[1].fx - pA[0].fx; v1.fy = pA[1].fy - pA[0].fy; v1.fz = pA[1].fz - pA[0].fz; 
v2.fx = pA[2].fx - pA[0].fx; v2.fy = pA[2].fy - pA[0].fy; v2.fz = pA[2].fz - pA[0].fz;
b = sqrt (v1.fx*v1.fx + v1.fy*v1.fy + v1.fz*v1.fz);
v1.fx /= b; v1.fy /= b; v1.fz /= b;
b = sqrt (v2.fx*v2.fx + v2.fy*v2.fy + v2.fz*v2.fz);
v2.fx /= b; v2.fy /= b; v2.fz /= b;
/* The H and the CH3 groups are placed along the two old bonds at the right bond length*/
pm[1].fx = pm[0].fx + v1.fx *rOH;
pm[1].fy = pm[0].fy + v1.fy *rOH;
pm[1].fz = pm[0].fz + v1.fz *rOH;
pm[2].fx = pm[0].fx + v2.fx *rMO;
pm[2].fy = pm[0].fy + v2.fy *rMO;
pm[2].fz = pm[0].fz + v2.fz *rMO;
}

