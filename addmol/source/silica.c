/*
 *
 *	silica.c:	This program generates an ascii file
 *			that contains coordinates and other info
 *			for a silica slab, based upon
 *			Lee and Rossky's 1994 JCP paper.
 *
 *	To call:	silica <#surface Si in x dimension> 
 *				<#surface Si in y dimension> 
 *				<zwall (Angstroms)> <ascii_file_name>
 *
 */

/*	updates: 	switch while loops to for loops, ensure that all O atoms
 *		 	are on same side of periodic boundary as their Si-O-H 
 *			silicon atom.
 */

#include	<stdio.h>
#include	<stdlib.h>
#include	<time.h>
#include	"../include/atomtypes.h"

main (int argc, char *argv[] )
{
	if (argc != 5) {
		fprintf(stderr, "usage: %s <surface Si (x)> <surface Si (y)> <zwall (Angstroms)> <ascii_file_name>\n",
			argv[0]);
		exit(0);
	}
	if (atoi(argv[2])%2 != 0) {
		fprintf(stderr, "<surface Si (y)> must be an even integer\n",
			argv[0]);
		exit(0);
	}
	
	FILE	*fp;
	double xSi = atof(argv[1]);
	double ySi = atof(argv[2]);
	double zwall = atof(argv[3]);
	double xdist = 5.0;
	double ydist = 4.3301;
	double xwall = xdist * xSi / 2;
	double ywall = ydist * ySi / 2;

	int max_atoms = 5000;
	struct atom
	{
		int type;
		double posx;
		double posy;
		double posz;
		int flags;
		int parent;
		int param1;
		int param2;
	};
	struct atom slab[max_atoms];

	int N = 0;
	int nsolute = 0;
	double EqTemp = 300.000000;
	double DEqTemp = 0.000000;
	double EqPress = 0.000000;
	double DEqPress = 0.000000;
	int xtrInQ = 0;
	int flag = A_FIXED | A_MAJOR | A_DRAW | A_MINOR;
	int parm1 = 0;
	int parm2 = 1;
	double x,y,z;
	int alt=0;
	int i,ix,iy;
	double xshift = 1.25 + 0.22857 + 0.0001;
	double yshift = 0.72120;

	//populate E plane, a' plane, H plane
	x = xdist / 2.0; y = 1.4424; z = 1.0206;
	alt = -1;
	for(iy = 0; iy < ySi; iy++)
	{
		for(ix = 0; ix < xSi; ix++)
		{
			slab[N].type = 0;
			slab[N].posx = x - xwall;
			slab[N].posy = y - ywall;
			slab[N].posz = z + 1.53093;
			slab[N].flags = A_MAJOR | A_DRAW | A_MINOR;
			slab[N].parent = N;
			slab[N].param1 = 2;
			slab[N].param2 = 3;
			N++;
			slab[N].type = 113;
			slab[N].posx = x - xwall;
			slab[N].posy = y - ywall;
			slab[N].posz = z;
			slab[N].flags = A_FIXED | A_DRAW;
			slab[N].parent = N-1;
			slab[N].param1 = 0;
			slab[N].param2 = 0;
			N++;
			slab[N].type = 1;
			slab[N].posx = x - xwall + 0.01;
			slab[N].posy = y - ywall + 0.01;
			slab[N].posz = z + 1.53093 + 0.9572; //TIP4P bond length
			slab[N].flags = A_DRAW;
			slab[N].parent = N-2;
			slab[N].param1 = 0;
			slab[N].param2 = 0;
			N++;
			x += xdist;
		}
		y += ydist; alt++;
		x = (double)(alt%2) * (xdist / 2.0);
	}

	//now step through the atom types and populate the slab
	//populate BCD plane
	x = 0; y = 0; z = 0;
	alt = 0;
	for(iy = 0; iy < ySi; iy++)
	{
		for(ix = 0; ix < xSi; ix++)
		{
			slab[N].type = 113;
			slab[N].posx = x - xwall;
			slab[N].posy = y - ywall;
			slab[N].posz = z;
			slab[N].flags = flag;
			slab[N].parent = N;
			slab[N].param1 = parm1;
			slab[N].param2 = parm2;
			x += xdist; N++;
		}
		y += ydist; alt++;
		x = (double)(alt%2) * (xdist / 2.0);
	}

	//populate b  plane
	x = 3.521225; y = 0.588130; z = 0.136355;
	alt = 0;
	for(iy = 0; iy < ySi; iy++)
	{
		for(ix = 0; ix < xSi; ix++)
		{
			slab[N].type = 0;
			slab[N].posx = x - xwall;
			slab[N].posy = y - ywall;
			slab[N].posz = z;
			slab[N].flags = flag;
			slab[N].parent = N;
			slab[N].param1 = parm1;
			slab[N].param2 = parm2;
			x += xdist; N++;
		}
		y += ydist; alt++;
		x = 3.521225 - (double)(alt%2) * (xdist / 2.0);
	}
	
	//populate c  plane
	x = 2.730153; y = 3.017922; z = 0.8842658;
	alt = 0;
	for(iy = 0; iy < ySi; iy++)
	{
		for(ix = 0; ix < xSi; ix++)
		{
			slab[N].type = 0;
			slab[N].posx = x - xwall;
			slab[N].posy = y - ywall;
			slab[N].posz = z;
			slab[N].flags = flag;
			slab[N].parent = N;
			slab[N].param1 = parm1;
			slab[N].param2 = parm2;
			x += xdist; N++;
		}
		y += ydist; alt++;
		x = 2.730153 - (double)(alt%2) * (xdist / 2.0);
	}

	//populate d  plane
	x = 1.021427; y = 1.117937; z = 0.510310;
	alt = 0;
	for(iy = 0; iy < ySi; iy++)
	{
		for(ix = 0; ix < xSi; ix++)
		{
			slab[N].type = 0;
			slab[N].posx = x - xwall;
			slab[N].posy = y - ywall;
			slab[N].posz = z;
			slab[N].flags = flag;
			slab[N].parent = N;
			slab[N].param1 = parm1;
			slab[N].param2 = parm2;
			x += xdist; N++;
		}
		y += ydist; alt++;
		x =  1.021427 - (double)(alt%2) * (xdist / 2.0);
	}

	if ((fp = fopen(argv[4], "w")) == NULL) {
		fprintf(stderr, "readbin: can't open ascii input-file\n");
		exit(1);
	}

	fprintf(fp,"Input File -- \t%s.c\n",argv[0]);
	fprintf(fp,"Filetype  -- \t.inp\n");

/*  ctime() inserts a newline character so we don't need to do it */
	fprintf(fp,"Datestamp -- \tdatestamp\n");   //%s", ctime(&datestamp));

	fprintf(fp, "Status -- \tRST\n");    //%s\n", status);
	fprintf(fp,"Number of atoms (natoms) -- \t%d\n", N);
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
	for (i=0; i<N; i++) {
		fprintf(fp,
"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
i,slab[i].type,slab[i].posx + xshift,slab[i].posy + yshift,slab[i].posz,slab[i].flags,slab[i].parent,slab[i].param1,slab[i].param2);
	}

/*
 *  Close the ascii file
 */
	fclose(fp);
}
