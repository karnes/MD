/*
 *
 *	findHBA.c:	This program reads a binary CH3CN/Silica input file,
 *			find a cenral silica site with nM = one CH3CN molecule 
 *			and writes an ascii file that includes only the
 *			one CH3CN molecule and the silica surface.
 *
 *	To call:	findHBM <binary_file> <ascii_file> <nM>
 *
 */

#include	<stdio.h>
#include	<stdlib.h>
#include	<typedefs.h>
#include	<globals.h>
#include	<math.h>
#include	<time.h>

#define abs(x)          ((x)>0?(x):-(x))
#define nCH3CN 1023
#define NOdist2 (3.5*3.5)
#define cosNOH 0.866025
main (argc, argv)
	int	argc;
	char	*argv[];
{
	int	i, j, k, l, n, ix, iy, nM, nHb;
	int pM[2];
	tripd d;
	double r2;
	FILE	*fp;

/*
 *  argv[1] and argv[2] contain the input and output file names
 */
	if (argc != 4) {
		fprintf(stderr, "usage: %s <binary_file> <ascii_file> <nM>\n",
			argv[0]);
		exit(0);
	}

nM = atoi(argv[3]);/* 1 (Si donates) or -1 (MeOH donates) or  2 (both)*/

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
/*loop over all the central Si units */
for (iy = 4; iy < 6; iy++){
	for (ix = 3; ix < 6; ix++){
		l = nCH3CN*3+(iy*9+ix)*3;/*This is the Oxygen of the current SiOH unit*/
		nHb = 0;/*counts how many CH3CN are bonded to this SiOH unit*/
		for (i=0; i< nCH3CN; i++){/*loop over all acetonitrile molecules*/
			n = 0;
			k = 3*i;
			d.fx = pos[k+1].fx-pos[l].fx;
			d.fy = pos[k+1].fy-pos[l].fy;
			d.fz = pos[k+1].fz-pos[l].fz;
			r2 = d.fx*d.fx+d.fy*d.fy+d.fz*d.fz;/*we do not need to image this vector*/
			if(r2 < NOdist2)
				n = getHb(k,l,d.fx,d.fy,d.fz);
			if (n == 1 && nM == 1){
				pM[nHb] = i;/*save this ACN*/
				nHb++;
			 	goto done;
			}
			if (n == -1 && nM == -1){
				pM[nHb] = i;/*save this ACN*/
				nHb++;
				goto done;
			}
			if (nM == 2 && n != 0) {
printf("nHb = %d i = %d\n",nHb,i*3);
				pM[nHb] = i;/*save this ACN*/
				nHb++;
			}
			if ( nM == 2 && nHb == 2) goto done;
		}
	}
}
if (nHb < abs(nM)){
	fprintf(stderr, " nHb = %d could not find any site\n",nHb);
	exit(0);
}
done:
				
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
	fprintf(fp,"Number of atoms (natoms) -- \t%d\n", natoms-nCH3CN*3+abs(nM)*3);
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
	for (i=0; i<nHb; i++) {
		j= 3*pM[i];
		for (k = 0; k <3; k++)
		fprintf(fp,
"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
i*3+k,atom[j+k].type,pos[j+k].fx,pos[j+k].fy,pos[j+k].fz,atom[j+k].flags,i*3,atom[j+k].param1,atom[j+k].param2);
	}
	for (i=nCH3CN*3; i<natoms; i++) {
		j = nHb*3;
		fprintf(fp,
"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
i-nCH3CN*3+j,atom[i].type,pos[i].fx,pos[i].fy,pos[i].fz,atom[i].flags,atom[i].parent-nCH3CN*3+j,atom[i].param1,atom[i].param2);
	}

/*
 *  Close the ascii file, and deallocate the memory
 */
	fclose(fp);
	cfree(atom);
	cfree(pos);
}
getHb(k,l,Ox,Oy,Oz)
int k,l; /* l=SiOH's O  k=ACN's C  */
double Ox,Oy,Oz;  /* vector from O(Si) to N(ACN) */
{
int n = 0;
tripd SiOHvec;
double rON, rOH, dot;
/* calculate OHO angles. either is < OHO ang, then H bonded */
	/* get  O(silica) -> H(silica) vector */
		SiOHvec.fx = pos[l+2].fx - pos[l].fx;
		SiOHvec.fy = pos[l+2].fy - pos[l].fy;
		SiOHvec.fz = pos[l+2].fz - pos[l].fz;
	/* O(Si) -> N(ACN) vector is (Ox,Oy,Oz) */
		rON = sqrt(Ox*Ox + Oy*Oy + Oz*Oz);
		rOH = sqrt(SiOHvec.fx*SiOHvec.fx + SiOHvec.fy*SiOHvec.fy + SiOHvec.fz*SiOHvec.fz);
		dot = Ox*SiOHvec.fx + Oy*SiOHvec.fy + Oz*SiOHvec.fz;
		if(dot/(rON*rOH) > cosNOH) /* H-bond donated by Silica*/
			n = 1;
	return(n);
}
