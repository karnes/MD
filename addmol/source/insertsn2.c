/*
 *	insertsn2.c:	reads xyz positions of 3 solute species
 *                      and insert in center of liquid-liquid box.
 *
 *	To call:	inssn2 <l/l binary> <solute binary> <nwater (atoms)> <nBCD> <ascii_file>
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
#define NSOL 3
#define BrOs 9
#define BCDs 147

void quicksort(float[][2], int, int);

main (argc, argv)
	int	argc;
	char	*argv[];
{
	int i,j,m,k,l,nw,nBCD,nBrO;
	int aflag,  pa1, pa2, ctype[NSOL];
	double sqrt(), atof(), r; // z1, dz, rdisk;
	double cx, cy, cz, totm, ms, BCDz;
	tripd pc[NSOL];
        char    sbuf[256];
	FILE	*fp;
	FILE	*fpo;
	char filetype2[9];
	time_t datestamp2;
	char status2[9];
	int natoms2;
	int nsolute2;
	double xwall2;
	double ywall2;
	double zwall2;
	double EqTemp2;
	double DEqTemp2;
	double EqPress2;
	double DEqPress2;
	parts *atom2;
	int xtrInQ2;
	tripd *pos2;

/*
 *  argv[1] and argv[2] contain the two input files names
 */
	if (argc != 6) {
		fprintf(stderr, "usage: %s <liq-liq binary> <solute binary> <nwater (atoms)> <nBCD> <ascii_file>\n",
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
/*
 *  open the  BCD binary input file
 */
	if ((fp = fopen(argv[2], "r")) == NULL) {
		fprintf(stderr, "can't open BCD binary input-file\n");
		exit(1);
	}

/*
 *  Read the header data, allocate space for the arrays, and read them too.
 */
	fread(filetype2, sizeof(char), 5, fp);
	fread(&datestamp2, sizeof(datestamp2), 1, fp);
	fread(status2, sizeof(char), 4, fp);
	fread(&natoms2, sizeof(natoms2), 1, fp);
	fread(&nsolute2, sizeof(nsolute2), 1, fp);
	fread(&xwall2, sizeof(xwall2), 1, fp);
	fread(&ywall2, sizeof(ywall2), 1, fp);
	fread(&zwall2, sizeof(zwall2), 1, fp);
	fread(&EqTemp2, sizeof(EqTemp2), 1, fp);
	fread(&DEqTemp2, sizeof(DEqTemp2), 1, fp);
	fread(&EqPress2, sizeof(EqPress2), 1, fp);
	fread(&DEqPress2, sizeof(DEqPress2), 1, fp);

	if ((atom2 = (parts *) calloc(natoms2, sizeof(atom2[0]))) == NULL ||
	    (pos2  = (tripd *) calloc(natoms2, sizeof(pos2[0]))) == NULL)
	 {
		fprintf(stderr, "out of memory\n");
		exit(1);
	}

	fread(atom2, sizeof(parts), natoms2, fp);
	fread(&xtrInQ2, sizeof(xtrInQ2), 1, fp);
	fread(pos2, sizeof(tripd), natoms2, fp);

/*
 *  Close the binary file
 */
	fclose(fp);
if (nsolute2 != 3){
	fprintf(stderr,"nsolute2 = %d must be equal to 3\n",nsolute);
	exit(1);
}
//fprintf(stderr,"natoms2 = %d\n",natoms2);

nBCD = atoi(argv[4]);

/*************************************	
 *	Center solutes at 0,0,0
 ************************************/

cx = cy = cz = totm = 0.0;
for(i=natoms2-nsolute2;i<natoms2;i++){
	if(i==natoms2-nsolute2)
	   ms = CH3_MASS;
	else
	   ms = CL_MASS;
	cx += pos2[i].fx * ms;
	cy += pos2[i].fy * ms;
	cz += pos2[i].fz * ms;
	totm += ms;
}
cx/=totm;
cy/=totm;
cz/=totm;

for(i=natoms2-nsolute2;i<natoms2;i++){
	pos2[i].fx -= cx;
	pos2[i].fy -= cy;
	pos2[i].fz -= cz;
}

nw = atoi(argv[3]);
nBrO = (natoms - nw*3 - nBCD*BCDs)/BrOs;

fprintf(stderr,"nwater = %d\n",nw);
fprintf(stderr,"natoms = %d\n",natoms);
nBrO = (natoms - nw*3 - nBCD*BCDs)/BrOs;		

natoms = natoms + 3;
/*
 *  Write the output to the ascii output-file in argv[3]
 */
	if ((fpo = fopen(argv[5], "w")) == NULL) {
		fprintf(stderr, "can't open ascii output-file %s\n",argv[5]);
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
/*print the water*/
for(j=0;j<nw;j++){/*loop over all the water molecules*/
   l = 3*j;
   fprintf(fpo,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
	l,atom[l].type,pos[l].fx,pos[l].fy,pos[l].fz,atom[l].flags,l,atom[l].param1,atom[l].param2);
   fprintf(fpo,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
	l+1,atom[l+1].type,pos[l+1].fx,pos[l+1].fy,pos[l+1].fz,atom[l+1].flags,l,atom[l+1].param1,atom[l+1].param2);
   fprintf(fpo,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
	l+2,atom[l+2].type,pos[l+2].fx,pos[l+2].fy,pos[l+2].fz,atom[l+2].flags,l,atom[l+2].param1,atom[l+2].param2);
}
/*print the Br-octane except those that were removed*/
for(j=0;j<nBrO;j++){/*loop over all the Br-oct molecules*/
   l = nw*3+BrOs*j;
   i = 0;
   if(j==0) 
      i = 80;
   for(k=0;k<BrOs;k++){
      fprintf(fpo,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
	l+k,atom[l+k].type,pos[l+k].fx,pos[l+k].fy,pos[l+k].fz + (float)(i),atom[l+k].flags,l,atom[l+k].param1,atom[l+k].param2);
   }
}
/*print the BCD, if there was on ein the input file */
if(nBCD==1){
   l = nw*3+nBrO*BrOs;
   for(k=0;k<BCDs;k++){
      fprintf(fpo,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
	l+k,atom[l+k].type,pos[l+k].fx,pos[l+k].fy,pos[l+k].fz+1.0,atom[l+k].flags,atom[l+k].parent,atom[l+k].param1,atom[l+k].param2);
   }
}

/* print the solute */
i = nw*3+nBrO*BrOs+nBCD*BCDs;
for(k=natoms2-nsolute2;k<natoms2;k++){
   fprintf(fpo,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
	i,atom2[k].type,pos2[k].fx,pos2[k].fy,pos2[k].fz,atom2[k].flags,i,atom2[k].param1,atom2[k].param2);
   i++;
}
/*
 *  Close the file, and deallocate the memory
 */
	fclose(fpo);
/*
 *  deallocate the memory
 */
	cfree(atom);
	cfree(pos);
	cfree(atom2);
	cfree(pos2);

}
