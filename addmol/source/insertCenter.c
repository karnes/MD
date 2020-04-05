/*
 *	insertCenter.c:	reads xyz positions of solute
 *                      and insert in center (0,0,0) of liquid-liquid box.
 *                      May delete atoms closest to center. Customize for use.
 *
 *	To call:	inscent <solvent binary> <solute binary> <ascii_file>
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
#define NSOL 78 //atoms in solute
#define solvSites 14 //atoms in solvent
#define rmSolv 5  //solvent molecules to remove

void quicksort(float[][2], int, int);

main (argc, argv)
	int	argc;
	char	*argv[];
{
	int i,j,k,l,inC[9999];
	float solvPos[9999][2];
//	int ctype[NSOL];
	int nsolv;
	double sqrt(), atof(), r;
	double cx, cy, cz, totm, ms;
//	tripd pc[NSOL];
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
		fprintf(stderr, "usage: %s <solvent binary> <solute binary> <ascii_file>\n",
			argv[0]);
		exit(0);
	}


/*
 *  open the solvent binary input file
 */
	if ((fp = fopen(argv[1], "r")) == NULL) {
		fprintf(stderr, "can't open solvent binary input-file\n");
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
 *  open the  solute binary input file
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
	fprintf(stderr,"nsolute = %d must be equal to 0. Exiting...\n",nsolute);
	exit(1);
}
if (natoms2 != NSOL){
	fprintf(stderr,"natoms2 = %d... must be equal to NSOL (%d). Exiting...\n",natoms2,NSOL);
	exit(1);
}
//fprintf(stderr,"natoms2 = %d\n",natoms2);

/*************************************	
 *	Center solute at 0,0,0
 ************************************/

cx = cy = cz = totm = 0.0;
for(i=0;i<natoms2;i++){
	if(atom2[i].type==5||atom2[i].type==7)
		ms=1.0;
	else if(atom2[i].type==45||atom2[i].type==0)
		ms=16.0;
	else if(atom2[i].type==9)
		ms=12.0;
	else if(atom2[i].type==94)
		ms=79.9;
	else if(atom2[i].type==55)
		ms=14.1;
	else{
		fprintf(stderr,"unknown atom type in BCD. Exiting...\n");
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
// shift solute molecule
for(i=0;i<natoms2;i++){
	pos2[i].fx -= cx;
	pos2[i].fy -= cy;
	pos2[i].fz -= cz;
}

if(natoms%solvSites!=0){
   fprintf(stderr,"natoms (%d) not a multiple of solvSites (%d). Exiting...\n",natoms,solvSites);
   exit(1);
}
   
nsolv = natoms/solvSites;

/* get solvent center atom distance from 0,0,0; put indices in a dummy array */
for(i=0;i<nsolv; i++)
{
   j=i*solvSites;
   r=sqrt(pos[j].fx*pos[j].fx+pos[j].fy*pos[j].fy+pos[j].fz*pos[j].fz);
   solvPos[i][0] = (float)i;
   solvPos[i][1] = r;
} 
/* give high values to indices > nsolv */
for (i= nsolv; i < 9999; i++)
{
   solvPos[i][0] = (float)i;
   solvPos[i][1] = 8888.8;
  }
/* flag lowest z waters for removal */ 
quicksort(solvPos,0,9999-1);
fprintf(stderr,"nsolv-rmSolv = %d\n",nsolv-rmSolv);
for(i=0;i<rmSolv;i++){
   inC[(int)(solvPos[i][0])]=0;
}
for(i=rmSolv;i<nsolv;i++){
   inC[(int)(solvPos[i][0])]=1;
}
/*
fprintf(stderr,"test: sorted water 0: index = %d, %d, z = %f\n",0,(int)(Wpos[0][0]),Wpos[0][1]);
fprintf(stderr,"test: sorted water 1: index = %d, %d, z = %f\n",1,(int)(Wpos[1][0]),Wpos[1][1]);
fprintf(stderr,"test: sorted water 500: index = %d, %d, z = %f\n",500,(int)(Wpos[500][0]),Wpos[500][1]);
fprintf(stderr,"test: sorted Br-oct 0: index = %d, %d, z = %f\n",0,(int)(Bpos[0][0]),Bpos[0][1]);
fprintf(stderr,"test: sorted Br-oct 1: index = %d, %d, z = %f\n",1,(int)(Bpos[1][0]),Bpos[1][1]);
fprintf(stderr,"test: sorted Br-oct 500: index = %d, %d, z = %f\n",500,(int)(Bpos[500][0]),Bpos[500][1]);
*/	  
/*
if((nWR+rmxW)!=rmW || (nBR+rmxB)!=rmB){
   fprintf(stderr,"nWR+rmxW = %d, nBR+rmxB = %d\n",nWR+rmxW,nBR+rmxB);
   fprintf(stderr,"wrong number of molecules removed... exiting.\n");
   exit(0);
}
*/
natoms = natoms + NSOL - solvSites*rmSolv;
/*
 *  Write the output to the ascii output-file in argv[3]
 */
if((fpo = fopen(argv[3], "w")) == NULL) {
   fprintf(stderr, "can't open ascii output-file %s\n",argv[3]);
   exit(1);
}

fprintf(fpo,"Input File -- \t%s\n",argv[1]);
fprintf(fpo,"Filetype  -- \t.inp\n");

/*  ctime() inserts a newline character so we don't need to do it */
fprintf(fpo,"Datestamp -- \t%s", ctime(&datestamp));

fprintf(fpo, "Status -- \t%s\n", status);
fprintf(fpo,"Number of atoms (natoms) -- \t%d\n", natoms);
fprintf(fpo,"Number of solute atoms (nsolute) -- \t%d\n", 0);
fprintf(fpo,"X box size (xwall) -- \t%f\n", xwall);
fprintf(fpo,"Y box size (ywall) -- \t%f\n", ywall);
fprintf(fpo,"Z box size (zwall) -- \t%f\n", zwall);
fprintf(fpo, "Equilibration Temperature (EqTemp) -- \t%f\n",EqTemp);
fprintf(fpo, "St. Dev. of Equil. Temp. (DEqTemp) -- \t%f\n",DEqTemp);
fprintf(fpo, "Equilibration Pressure (EqTemp) -- \t%f\n",EqPress);
fprintf(fpo, "St. Dev. of Equil. Press. (DEqPress) --\t%f\n",DEqPress);
fprintf(fpo, "Extra file flag (xtrInQ) -- \t%d\n", xtrInQ);
/*print the solvent -- no change*/
fprintf(fpo,
"   #  type  x position    y position    z position  flags parent param1 param2\n");
/*print the solvent except those that were removed*/
i = 0;/*index run over atoms of solvent that are kept*/
for(j=0;j<nsolv;j++){/*loop over all the solvent molecules*/
   if(inC[j] == 1){/*print desired solvents*/
      for(k=0;k<solvSites;k++){
         l = solvSites*j+k;
         fprintf(fpo,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
i,atom[l].type,pos[l].fx,pos[l].fy,pos[l].fz,atom[l].flags,(int)(i/solvSites)*solvSites,atom[l].param1,atom[l].param2);
         i++;
      }
   }
}
// print solute
for(k=0;k<natoms2;k++){/*loop over all the ETH molecules*/
   fprintf(fpo,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
k+i,atom2[k].type,pos2[k].fx,pos2[k].fy,pos2[k].fz,atom2[k].flags,atom2[k].parent + i,atom2[k].param1,atom2[k].param2);
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

void quicksort(float x[][2], int first, int last)
{
	int pivot,j,i;
	float temp0,temp1;

	if(first < last)
	{
		pivot = first;
		i = first;
		j = last;

		while(i < j)
		{
			while(x[i][1]<=x[pivot][1] && i<last)
				i++;
			while(x[j][1]>x[pivot][1])
				j--;
			if(i<j)
			{
				temp0=x[i][0];
				temp1=x[i][1];
				x[i][0]=x[j][0];
				x[i][1]=x[j][1];
				x[j][0]=temp0;
				x[j][1]=temp1;
			}
		}

		temp0=x[pivot][0];
		temp1=x[pivot][1];
		x[pivot][0]=x[j][0];
		x[pivot][1]=x[j][1];
		x[j][0]=temp0;
		x[j][1]=temp1;
		
		quicksort(x,first,j-1);
		quicksort(x,j+1,last);
	}
}
