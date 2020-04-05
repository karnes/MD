/*
 *	chopWatBrO.c:	reads in water/BrO binary file
 *			deletes water, shifts BrO molecules so CoM ~ 0 Z
 *			and chops into a cube of user-defined edge
 *
 *	To call:	cutWatBox <water binary> <output .as> <box edge (Å)>
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
#define BrOs 9
main (argc, argv)
	int	argc;
	char	*argv[];
{
	int	i,j,k,l,nw,nWR,nBR,inW[5100],inB[5100];
	double sqrt(), atof(), r; // z1, dz, rdisk;
	double cx, cy, cz, totm, ms;
        char    sbuf[256];
	FILE	*fp;
	FILE	*fpo;
	double edge;
	int nBrO;

/*
 *  argv[1] and argv[2] contain the two input files names
 */
	if (argc != 4) {
		fprintf(stderr, "usage: %s <l/l binary> <output .as> <box edge (Å)>\n",
			argv[0]);
		exit(0);
	}
	
	edge = atof(argv[3]);

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


nw = 3944;/* number of water molecules in the original file*/
nWR = 0;
for(i=0;i<nw;i++){
   nWR++;
   inW[i] = 1; // remover all water
}	

fprintf(stderr," %d water molecules will be removed\n",nWR);


/*
 * shift all atoms -35 Z *
 */
for(i=0;i<natoms;i++){
   pos[i].fz-=35.0;
}

/* 
 * remove all BrO molecules with image center outside +/- 25 Å
 */
nBrO=1000;
nBR = 0;
for(i=0;i<nBrO;i++){
   k = nw*3+i*BrOs;
   inB[i] = 0; //assume in box
   if(fabs(pos[k].fx)>(edge/2.0) || fabs(pos[k].fy)>(edge/2.0) || fabs(pos[k].fz)>(edge/2.0)){
      nBR++;
      inB[i] = 1; // remover this Br-Octane
   }
}  


natoms = natoms - nWR*3 - nBR*BrOs;

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
	fprintf(fpo,"Z box size (zwall) -- \t%f\n", edge/2.0);
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
/*print the BrOctane except those that were removed*/
i = 0;/*index run over atoms of Br-oct that are kept*/
for (j=0; j<nBrO; j++) {/*loop over all the Br-oct molecules*/
   if (inB[j] == 0){/*if this Br-oct is in the box; print it*/
	l = nw*3+BrOs*j;
	for(k=0;k<BrOs;k++){
	   fprintf(fpo,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
		i,atom[l+k].type,pos[l+k].fx,pos[l+k].fy,pos[l+k].fz,atom[l+k].flags,(i/BrOs)*BrOs,atom[l+k].param1,atom[l+k].param2);
	   i++;
	}
   }
}
/*
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
//i++;*/
fprintf(stderr,"new natoms = %d, density = %f.\n",i,(float)(i)/(edge*edge*edge));
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
