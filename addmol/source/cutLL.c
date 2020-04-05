/*
 *	cutLL.c:reads in water/DDC/Er binary file
 *			chops xwall and ywall to user-defined lengths
 *
 *	To call:	cutll <WDE binary> <output .as> <x,y-box edge (Å)> <nDDC> <nEr>
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
#define DDCs 12
#define Ers 1
//#define rmW 0//1460 // number of total water to remove
//#define rmB 0//370
main (argc, argv)
	int	argc;
	char	*argv[];
{
	int	i,j,k,l,nw,nWR,nDR,inW[12000],inD[12000];
	double sqrt(), atof(), r; // z1, dz, rdisk;
	double cx, cy, cz, totm, ms;
        char    sbuf[256];
	FILE	*fp;
	FILE	*fpo;
	double edge;
	int nDDC,nEr;
//	int rmxW,rmxB;
//	float Wpos[9999][2],Bpos[9999][2];

/*
 *  argv[1] and argv[2] contain the two input files names
 */
	if (argc != 6) {
		fprintf(stderr, "usage: %s <l/l binary> <output .as> <x,y-box edge (x/ywall*2) (Å)> <nDDC in bin input> <nEr in bin input>\n",
			argv[0]);
		exit(0);
	}
	
	edge = atof(argv[3]);
	nDDC = atoi(argv[4]);
	nEr = atoi(argv[5]);

/*
 *  open the l/l binary input file
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


nw = (natoms - nDDC*DDCs - nEr*Ers)/3;/* number of water molecules in the original file*/
/* 
 * remove all water molecules with image center outside +/- desired x/y wall
 */
nWR = 0;
for(i=0;i<nw;i++){
   k = i*3;
   inW[i]=0; // assume in box
   if(fabs(pos[k].fx) > (edge/2.0) || fabs(pos[k].fy) > (edge/2.0)){
      nWR++;
      inW[i] = 1; // remove this water
   }
}	

fprintf(stderr," %d water molecules outside of new box\n",nWR);

/* sort waters, remove enough so that total removed = rmW */	
/* assign water z pos, indices to a dummy array */
/*
if(nw > 0){
   rmxW = 0;
   for (i=0;i< nw; i++)
   {
	Wpos[i][0] = (float)i;
	if(inW[i]==0)
	   Wpos[i][1] = pos[i*3].fz;
	else // if already removed, give a high value.
	   Wpos[i][1] = 888.0;
   } 
   // give high values to indices > nw/3 
   for (i= nw; i < 9999; i++)
   {
	Wpos[i][0] = (float)i;
	Wpos[i][1] = 888.0;
   }
   // flag lowest z waters for removal  
   quicksort(Wpos,0,9999-1);
fprintf(stderr,"rmW-nWR = %d\n",rmW-nWR);
   i=0;
   while(rmxW<(rmW-nWR)){
      if(inW[(int)(Wpos[i][0])]==0){
         inW[(int)(Wpos[i][0])] = 1;
         rmxW++;
      }
      i++;
   }
}
else{
   rmxW = 0;
}
*/
/* 
 * remove all BrO molecules with image center outside +/- 22.5 Å
 */
nDR = 0;
for(i=0;i<nDDC;i++){
   k = nw*3+i*DDCs;
   inD[i] = 0; //assume in box
   if(fabs(pos[k].fx)>(edge/2.0) || fabs(pos[k].fy)>(edge/2.0)){
      nDR++;
      inD[i] = 1; // remove this Br-Octane
   }
}  


fprintf(stderr," %d DDC molecules outside of new box\n",nDR);
/*
rmxB = 0;
if(nBrO > 0){
   for (i=0;i<nBrO;i++)
   {
      Bpos[i][0] = (float)i;
      Bpos[i][1] = -(pos[nw*3 + i*9].fz);
   }
   for (i= nBrO; i < 9999; i++)
   {
	Bpos[i][0] = (float)i;
	Bpos[i][1] = 888.0;
   }
   quicksort(Bpos,0,9999-1);
fprintf(stderr,"rmB-nBR = %d\n",rmB-nBR);
   i=0;
   while(rmxB<(rmB-nBR)){
      if(inB[(int)(Bpos[i][0])]==0){
         inB[(int)(Bpos[i][0])] = 1;
         rmxB++;
      }
      i++;
   }
}

natoms = natoms - (nWR+rmxW)*3 - (nBR+rmxB)*BrOs;
fprintf(stderr,"total removed: water: %d, Br-oct: %d\n",nWR+rmxW,nBR+rmxB);
*/
natoms = natoms - (nWR)*3 - (nDR)*DDCs;
/*
 *  Write the output to the ascii output-file in argv[3]
 */
	if ((fpo = fopen(argv[2], "w")) == NULL) {
		fprintf(stderr, "can't open ascii output-file %s\n",argv[2]);
		exit(1);
	}

	fprintf(fpo,"Input File -- \t%s\n",argv[1]);
	fprintf(fpo,"Filetype  -- \t.inp\n");

/*  ctime() inserts a newline character so we don't need to do it */
	fprintf(fpo,"Datestamp -- \t%s", ctime(&datestamp));

	fprintf(fpo, "Status -- \t%s\n", status);
	fprintf(fpo,"Number of atoms (natoms) -- \t%d\n", natoms);
	fprintf(fpo,"Number of solute atoms (nsolute) -- \t%d\n", 0);
	fprintf(fpo,"X box size (xwall) -- \t%f\n", edge/2.0);
	fprintf(fpo,"Y box size (ywall) -- \t%f\n", edge/2.0);
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
/*print the waters, except those that were removed*/
i = 0;/*index run over atoms of water that are kept*/
for (j=0; j<nw; j++) {/*loop over all the water molecules*/
   if (inW[j] == 0){/*if this water is in the box; print it*/
	l = j*3;
	for(k=0;k<3;k++){
	   fprintf(fpo,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
		i,atom[l+k].type,pos[l+k].fx,pos[l+k].fy,pos[l+k].fz,atom[l+k].flags,(i/3)*3,atom[l+k].param1,atom[l+k].param2);
	   i++;
	}
   }
}
/*print the DDC except those that were removed*/
//i = 0; /*maintain same counter...*/
for (j=0; j<nDDC; j++) {/*loop over all the Br-oct molecules*/
   if (inD[j] == 0){/*if this Br-oct is in the box; print it*/
	l = nw*3+DDCs*j;
	for(k=0;k<DDCs;k++){
	   fprintf(fpo,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
		i,atom[l+k].type,pos[l+k].fx,pos[l+k].fy,pos[l+k].fz,atom[l+k].flags,3*(nw-nWR)+((i-(3*(nw-nWR)))/DDCs)*DDCs,atom[l+k].param1,atom[l+k].param2);
	   i++;
	}
   }
}
/*print the Er*/
//i = 0; /*maintain same counter...*/
for (j=0; j<nEr; j++) {/*loop over BCD molecule*/
	l = nw*3+DDCs*nDDC;
	for(k=0;k<Ers;k++){
	   fprintf(fpo,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
		i,atom[l+k].type,pos[l+k].fx,pos[l+k].fy,pos[l+k].fz,atom[l+k].flags,atom[l+k].parent - 3*(nWR) - DDCs*(nDR),atom[l+k].param1,atom[l+k].param2);
	   i++;
	}
}
//fprintf(stderr,"new natoms = %d, density = %f.\n",i,(float)(i)/(edge*edge*edge));
/*
 *  Close the file, and deallocate the memory
 */
	fclose(fpo);
/*
 *  deallocate the memory
 */
	cfree(atom);
	cfree(pos);

}
/*
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
*/
