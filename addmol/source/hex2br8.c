/*
 *	hex2br8.c:	This program reads a water- hexane box 
 *                      and mutates each hexane to a bromo-octane
 *
 *	To call:	hex2br8 <ascii_file>
 *
 */
#include "typedefs.h"
#include "atomtypes.h"
#include "globals.h"
#include <stdio.h>
#include <time.h>
#define abs(x)          ((x)>0.?(x):-(x))
#define sq(x)	((x)*(x))
#define hsgn(x)         ((x)<0.?-.5:.5)
#define NumWater 986
#define NumHEX 250
#define BrC (1.966 - 0.2)
#define CC (1.535 - 0.2)

main (argc, argv)
	int	argc;
	char	*argv[];
{
	int	i,j,k,l,nw,nhex;
	int aflag,  pa1, pa2;
	double sqrt(),atof(),r2,r;
        char    sbuf[256];
	FILE	*fp;
	FILE	*fpo;
	tripd thf, fts; /* carbon [index=]3-to-5 and 4-to-6 vectors */
	tripd Br, C7, C8;

	nw = NumWater;
	nhex = NumHEX;

	if (argc != 2) {
		fprintf(stderr, "usage: %s <output ascii_file>\n",
			argv[0]);
		exit(0);
	}

/*
 *  open the  water-hexane binary input file
 */
	if ((fp = fopen("wathex000.bin", "r")) == NULL) {
		fprintf(stderr, "can't open water-hexane binary input-file\n");
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

natoms = natoms + nhex*3;
/*
 *  Write the output to the ascii output-file in argv[3]
 */
	if ((fpo = fopen(argv[1], "w")) == NULL) {
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
	fprintf(fpo,"   #  type  x position    y position    z position  flags parent param1 param2\n");
	for (i=0;i<nw*3;i++) {/*loop over all the water molecules*/
		fprintf(fpo,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
		i,atom[i].type,pos[i].fx,pos[i].fy,pos[i].fz,atom[i].flags,atom[i].parent,atom[i].param1,atom[i].param2);
	}
	for(i=0;i<nhex;i++){
	   k=nw*3+i*6; //index of hexane center united atom
	   thf.fx = pos[k+4].fx - pos[k+2].fx;
	   thf.fy = pos[k+4].fy - pos[k+2].fy;
	   thf.fz = pos[k+4].fz - pos[k+2].fz;
	   r2 = thf.fx*thf.fx + thf.fy*thf.fy + thf.fz*thf.fz;
	   r = sqrt(r2);
	   thf.fx /= (r/BrC);
	   thf.fy /= (r/BrC);
	   thf.fz /= (r/BrC);
	   Br.fx = pos[k+4].fx + thf.fx;
	   Br.fy = pos[k+4].fy + thf.fy;
	   Br.fz = pos[k+4].fz + thf.fz;

	   fts.fx = pos[k+5].fx - pos[k+3].fx;
	   fts.fy = pos[k+5].fy - pos[k+3].fy;
	   fts.fz = pos[k+5].fz - pos[k+3].fz;
	   r2 = fts.fx*fts.fx + fts.fy*fts.fy + fts.fz*fts.fz;
	   r = sqrt(r2);
	   fts.fx /= (r/CC);
	   fts.fy /= (r/CC);
	   fts.fz /= (r/CC);
	   C7.fx = pos[k+5].fx + fts.fx;
	   C7.fy = pos[k+5].fy + fts.fy;
	   C7.fz = pos[k+5].fz + fts.fz;
	   C8.fx = C7.fx + fts.fx;
	   C8.fy = C7.fy + fts.fy;
	   C8.fz = C7.fz + fts.fz;
	   
	   l = nw*3+i*9;

	fprintf(fpo,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
	    l,atom[k+2].type,pos[k+4].fx + 0.1,pos[k+4].fy,pos[k+4].fz,atom[k].flags,l,8,9);
	fprintf(fpo,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
	    l+1,94,Br.fx,Br.fy,Br.fz,atom[k+1].flags,l,atom[k+1].param1,atom[k+1].param2);
	fprintf(fpo,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
	    l+2,atom[k+2].type,pos[k+2].fx,pos[k+2].fy,pos[k+2].fz,atom[k+1].flags,l,atom[k+1].param1,atom[k+1].param2);
	fprintf(fpo,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
	    l+3,atom[k+2].type,pos[k].fx,pos[k].fy,pos[k].fz,atom[k+1].flags,l,atom[k+1].param1,atom[k+1].param2);
	fprintf(fpo,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
	    l+4,atom[k+2].type,pos[k+1].fx,pos[k+1].fy,pos[k+1].fz,atom[k+1].flags,l,atom[k+1].param1,atom[k+1].param2);
	fprintf(fpo,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
	    l+5,atom[k+2].type,pos[k+3].fx,pos[k+3].fy,pos[k+3].fz,atom[k+1].flags,l,atom[k+1].param1,atom[k+1].param2);
	fprintf(fpo,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
	    l+6,atom[k+2].type,pos[k+5].fx,pos[k+5].fy,pos[k+5].fz,atom[k+1].flags,l,atom[k+1].param1,atom[k+1].param2);
	fprintf(fpo,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
	    l+7,atom[k+2].type,C7.fx + 0.1,C7.fy + 0.1,C7.fz,atom[k+1].flags,l,atom[k+1].param1,atom[k+1].param2);
	fprintf(fpo,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
	    l+8,atom[k+4].type,C8.fx - 0.1,C8.fy - 0.1,C8.fz,atom[k+1].flags,l,atom[k+1].param1,atom[k+1].param2);
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
}
