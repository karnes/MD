/*
 *	insert3sol.c:	reads xyz positions of 3 solute species
 *                      and insert in center of liquid-liquid box.
 *
 *	To call:	ins3sol <l/l binary> <solute binary> <nwater (atoms)> <nBCD> <ascii_file>
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
#define rmaxw 3.0
#define rmaxb 3.0
#define rmB 0 //20
#define rmW 0 //84

void quicksort(float[][2], int, int);

main (argc, argv)
	int	argc;
	char	*argv[];
{
	int i,j,m,k,l,nw,nBCD,nBrO,nWR,nBR,inC[9999], inB[9999];
	float Wpos[9999][2],Bpos[9999][2];
	int aflag,  pa1, pa2, ctype[NSOL];
	int rmxW, rmxB,totRem;
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
//	pos2[i].fz += BCDz;
}

nw = atoi(argv[3]);
nBrO = (natoms - nw*3 - nBCD*BCDs)/BrOs;
nw/=3;
fprintf(stderr,"nwater = %d\n",nw);
fprintf(stderr,"natoms = %d\n",natoms);
nWR = 0;
for(i=0;i<nw;i++){
   inC[i] = 0;/* first assume water is outside solute radius*/
   for(j=natoms2-nsolute2;j<natoms2;j++){
      for(k=0;k<3;k++){	   
	r = sqrt(sq(pos[i*3+k].fx - pos2[j].fx) + sq(pos[i*3+k].fy - pos2[j].fy) + sq(pos[i*3+k].fz - pos2[j].fz));
//if(j==0)fprintf(stderr,"r = %f\n",r);
	if (r < rmaxw){
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
/* sort waters, remove enough so that total removed = rmW */	
/* assign water z pos, indices to a dummy array */
if(nw > 0){
   rmxW = 0;
   for (i=0;i< nw; i++)
   {
	Wpos[i][0] = (float)i;
	Wpos[i][1] = pos[i*3].fz;
   } 
   /* give high values to indices > nw/3 */
   for (i= nw; i < 9999; i++)
   {
	Wpos[i][0] = (float)i;
	Wpos[i][1] = 888.0;
   }
   /* flag lowest z waters for removal */ 
   quicksort(Wpos,0,9999-1);
fprintf(stderr,"rmW-nWR = %d\n",rmW-nWR);
   i=0;
   while(rmxW<(rmW-nWR)){
      if(inC[(int)(Wpos[i][0])]==0){
         inC[(int)(Wpos[i][0])] = 1;
         rmxW++;
      }
      i++;
   }
}
else{
   rmxW = 0;
}
nBR = 0;
for(i=0;i<nBrO;i++){
   m = nw*3 + i*BrOs;
   inB[i] = 0; /* assume BrO is outside solute radius */
   for(j=natoms2-nsolute2;j<natoms2;j++){
      for(k=0;k<BrOs;k++){
	 r = sqrt(sq(pos[m+k].fx - pos2[j].fx) + sq(pos[m+k].fy - pos2[j].fy) + sq(pos[m+k].fz - pos2[j].fz));
	 if(r<rmaxb){
	    if(inB[i]==0){
		nBR++;
		inB[i] = 1; /* remove this BrOct */
		j = natoms2;
		k = BrOs;
	    }
	    else{
		fprintf(stderr, "error.... exiting. \n");
		exit(0);
	    }
	 }
      }
   }
}
if((natoms - nw*3 - nBCD*BCDs)%BrOs!=0){
   fprintf(stderr,"error in number of atoms? exiting...\n");
   exit(0);
}
nBrO = (natoms - nw*3 - nBCD*BCDs)/BrOs;		
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
/*
fprintf(stderr,"test: sorted water 0: index = %d, %d, z = %f\n",0,(int)(Wpos[0][0]),Wpos[0][1]);
fprintf(stderr,"test: sorted water 1: index = %d, %d, z = %f\n",1,(int)(Wpos[1][0]),Wpos[1][1]);
fprintf(stderr,"test: sorted water 500: index = %d, %d, z = %f\n",500,(int)(Wpos[500][0]),Wpos[500][1]);
fprintf(stderr,"test: sorted Br-oct 0: index = %d, %d, z = %f\n",0,(int)(Bpos[0][0]),Bpos[0][1]);
fprintf(stderr,"test: sorted Br-oct 1: index = %d, %d, z = %f\n",1,(int)(Bpos[1][0]),Bpos[1][1]);
fprintf(stderr,"test: sorted Br-oct 500: index = %d, %d, z = %f\n",500,(int)(Bpos[500][0]),Bpos[500][1]);
*/	  
fprintf(stderr," %d water molecules will be removed to make room for the solute. %d total to remove.\n",nWR,rmW);
fprintf(stderr," %d Br-octane molecules will be removed to make room for the solute. %d total to remove.\n",nBR,rmB);
/*
if((nWR+rmxW)!=rmW || (nBR+rmxB)!=rmB){
   fprintf(stderr,"nWR+rmxW = %d, nBR+rmxB = %d\n",nWR+rmxW,nBR+rmxB);
   fprintf(stderr,"wrong number of molecules removed... exiting.\n");
   exit(0);
}
*/
natoms = natoms - (nWR+rmxW)*3 - (nBR+rmxB)*BrOs + NSOL;
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
	}
}
/*print the Br-octane except those that were removed*/
/* i is index run over atoms of Br-octane that are kept, continued from water */
for (j=0; j<nBrO;j++) {/*loop over all the Br-oct molecules*/
	if (inB[j] == 0){/*if this Br-oct is outside, print it*/
	   l = nw*3+BrOs*j;
	      for(k=0;k<BrOs;k++){
		 fprintf(fpo,
"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
i,atom[l+k].type,pos[l+k].fx,pos[l+k].fy,pos[l+k].fz,atom[l+k].flags,i-k,atom[l+k].param1,atom[l+k].param2);
i++;
              }
	}
}
/*print the BCD, if there was on ein the input file */
/* i is index run over atoms in the original solvent input file, continued from water & Br-oct */
if(nBCD==1){
   l = nw*3+nBrO*BrOs;
   totRem = 3*nWR + nBR*BrOs;
      for(k=0;k<BCDs;k++){
	 fprintf(fpo,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
i,atom[l+k].type,pos[l+k].fx,pos[l+k].fy,pos[l+k].fz,atom[l+k].flags,atom[l+k].parent-totRem,atom[l+k].param1,atom[l+k].param2);
i++;
              }
	}

/* print the solute */
for (k=natoms2-nsolute2; k<natoms2; k++) {
		fprintf(fpo,
"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
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
