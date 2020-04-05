/*
 *	insertBCDw.c:	reads xyz positions of BCD
 *                      and insert in center of water box
 *
 *	To call:	insBCDw <water binary> <BCD binary> <ascii_file>
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
#define NSOL 147
#define rmax 2.0
main (argc, argv)
	int	argc;
	char	*argv[];
{
	int	i,j,k,l,nw,nWR, inC[5100];
	int aflag,  pa1, pa2, ctype[NSOL];
	double sqrt(), atof(), r; // z1, dz, rdisk;
	double cx, cy, cz, totm, ms;
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
	if (argc != 4) {
		fprintf(stderr, "usage: %s <water binary> <BCD binary> <ascii_file>\n",
			argv[0]);
		exit(0);
	}

//rdisk = atof(argv[4]);/* cavity radius*/
//dz = atof(argv[5]);/* half thickness of the cylindrecial cavity)*/
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
if (nsolute != 0){
	fprintf(stderr,"nsolute = %d must be equal to 0\n",nsolute);
	exit(1);
}
//fprintf(stderr,"natoms2 = %d\n",natoms2);
/*	
 *	Center BCD at 0,0,0
 */
/*
cx = cy = cz = totm = 0.0;
for(i=0;i<natoms2;i++){
	if(atom2[i].type==5||atom2[i].type==7)
		ms=1.0;
	else if(atom2[i].type==45||atom2[i].type==0)
		ms=16.0;
	else if(atom2[i].type==20)
		ms=12.0;
	else{
		fprintf(stderr,"unknown atom type in BCD... exiting.\n");
		exit(0);
	}
	
	cx += pos2[i].fx * ms;
	cy += pos2[i].fy * ms;
	cz += pos2[i].fz * ms;
	totm += ms;
}
cx/=totm;
cy/=totm;
cz/=totm;
for(i=0;i<natoms2;i++){
	pos2[i].fx -= cx;
	pos2[i].fy -= cy;
	pos2[i].fz -= cz;
}
*/
nw = natoms/3;/* number of water molecules in the original file*/
nWR = 0;
for(i=0;i<nw;i++){
   inC[i] = 0;/* first assume water is outside BCD radius*/
   for(j=0;j<natoms2;j++){
      for(k=0;k<1;k++){ // only check water oxygens	   
	r = sqrt(sq(pos[i*3+k].fx - pos2[j].fx) + sq(pos[i*3+k].fy - pos2[j].fy) + sq(pos[i*3+k].fz - pos2[j].fz));
//if(j==0)fprintf(stderr,"r = %f\n",r);
	if (r < rmax){
		if(inC[i]==0){
		   nWR++;
		   inC[i] = 1;/* this water will be removed*/
		   j = natoms2;
		   k=3;
		}
		else{
		   fprintf(stderr,"error... exiting. \n");
		   exit(0);
		}
	}
      }
   }
}	

fprintf(stderr," %d water molecules will be removed\n",nWR);

natoms = natoms - nWR*3  + NSOL;
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
	fprintf(fpo,"Number of solute atoms (nsolute) -- \t%d\n", 1/*NSOL*/);
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
/*print the water except those that were removed*/
i = 0;/*index run over atoms of water that are kept*/
for (j=0; j<nw; j++) {/*loop over all the water molecules*/
	if (inC[j] == 0){/*if this water is outside, print it*/
		l = 3*j;
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
//		fprintf(fpo,
//"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
//i,atom[l+3].type,pos[l+3].fx,pos[l+3].fy,pos[l+3].fz,atom[l+3].flags,i-1,atom[l+3].param1,atom[l+3].param2);
//i++;
	}
}
//j = natoms - nWR * 3;
for (k=0; k<natoms2; k++) {/*loop over all the ETH molecules*/
		fprintf(fpo,
"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
k+i,atom2[k].type,pos2[k].fx,pos2[k].fy,pos2[k].fz,atom2[k].flags,atom2[k].parent + i,atom2[k].param1,atom2[k].param2);
}
if(xwall==ywall && ywall==zwall)
   fprintf(stderr,"(OCT) density = %f atoms / Å^3\n",(float)(natoms)/(4.0*xwall*ywall*zwall));
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
