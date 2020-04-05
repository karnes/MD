/*
 *	makeSlab111:	This program creates an ascii text file for a silica
 *			[beta cristobalite] (111) slab surface.
 *
 *			Data obtained from Lee and Rossky, JCP, 1994. 
 *			Rotated unit cell with PBC so that 'low density 111 
 *			crystal face' contacts liquid phase.
 *
 *	To call:	****** not implemented yet*****  <#of unit cells> <ascii_file>
 *
 */

#include	<stdio.h>
#include	<time.h>
#include	<math.h>

main ()
{
	FILE	*fp;
	int layers = 4;
	int cellsx = 12;
	int cellsy = 12;

/*	if (argc != 3) {
		fprintf(stderr, "usage: %s <binary_file> <ascii_file>\n",
			argv[0]);
		exit(0);
	}

*/
/*	Set dummy data for header   */
	
	double unitCellEdge = 5.0 / sqrt(2.0);
	int natoms; 
	//layers * atoms in unit cell * cells in x direction 
	//      * cells in y direction + top layer of hydrogens
	
	//status = "RST";
	int nsolute = 0;
	double xwall = cellsx * unitCellEdge;  //45.60;
	double ywall = cellsy * unitCellEdge;  //43.88;
	double zwall = layers * unitCellEdge;;
	double EqTemp = 300.000000;
	double DEqTemp = 0.000000;
	double EqPress = 0.000000;
	double DEqPress = 0.000000;
	int xtrInQ = 0;
	
	//unit cell data
	int i,j;
	double edge = 5.0 / sqrt(2.0);
	natoms = (layers * (5+4) * cellsx * cellsy);// + cellsx * cellsy; 
	//layers * atoms in unit cell * cells in x direction 
	//      * cells in y direction + top layer of hydrogens
	
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
	
	struct atom uCell[10]; //for top layer, use oxygen a'[9] instead of a[5]
	
	struct atom A = {113, edge*0.788675, edge*0.788675, edge*0.15470054, 20, 2, 0, 1};//A
	uCell[0] = A;
	struct atom B = {113, edge*0.0, edge*1.0, edge*0.0, 20, 2, 0, 1};//B
	uCell[1] = B;
	struct atom C = {113, edge*1.0, edge*0.0, edge*0.0, 20, 2, 0, 1};//C
	uCell[2] = C;
	struct atom D = {113, edge*0.3660254, 0.3660254, 0, 20, 2, 0, 1};//D
	uCell[3] = D;
	struct atom E = {113, edge*0.78867513, edge*0.78867513, edge*0.28867513, 20, 2, 0, 1};//E
	uCell[4] = E;
	
	double ax = 0.78867513;
	double ay = 0.78867513;
	
	struct atom a = {0, edge*ax, edge*ay, edge*0.72168784, 20, 2, 0, 1};//a or a'
	uCell[5] = a;
	struct atom b = {0, edge*0.44722285, edge*0.94722285, edge*0.038567, 20, 2, 0, 1};//b
	uCell[6] = b;
	struct atom c = {0, edge*0.84145228, edge*0.34145228, edge*0.25010814, 20, 2, 0, 1};//c
	uCell[7] = c;
	struct atom d = {0, edge*0.16895027, edge*0.98575027, edge*0.14433757, 20, 2, 0, 1};//d
	uCell[8] = d;
	
	ax = 0.69707513;
	ay = 0.8827513;

	struct atom slab[natoms];
	i=0;

	//build silica layers, starting at top layer
	int xaxis, yaxis, zaxis;
	for(zaxis = layers; zaxis > 0; zaxis--)
	{
		for(xaxis = 0; xaxis <  cellsx; xaxis++)
		{	
			for(yaxis = 0; yaxis <  cellsy; yaxis++)
			{	
		
				for(j=0; j<9; j++)
				{
					slab[i].type = uCell[j].type;
					slab[i].posx = (edge * xaxis) + uCell[j].posx;
					slab[i].posy = (edge * yaxis) + uCell[j].posy;
					slab[i].posz = (edge * (zaxis - 1)) + uCell[j].posz;
					slab[i].flags = uCell[j].flags;
					slab[i].parent = i;
					slab[i].param1 = uCell[j].param1;
					slab[i].param2 = uCell[j].param2;
					i++;
				}
			}
		}

		uCell[5].posx = edge * ax;
		uCell[5].posy = edge * ay;
		
	}

/*
 *  Write the output to the ascii input-file in argv[2]
 */
	if ((fp = fopen("slab.as", "w")) == NULL) {
		fprintf(stderr, "readbin: can't open ascii input-file\n");
		/*exit(1);*/
	}

	fprintf(fp,"Input File -- \tRST\n");  //  %s\n",argv[1]);
	fprintf(fp,"Filetype  -- \t.inp\n");

/*  ctime() inserts a newline character so we don't need to do it */
	fprintf(fp,"Datestamp -- \tdatestamp\n"); //  %s", ctime(&datestamp));

	fprintf(fp, "Status -- \tRST\n");  //  %s\n", status);
	fprintf(fp,"Number of atoms (natoms) -- \t%d\n", natoms);
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
	for (i=0; i<natoms; i++) {
		fprintf(fp,
"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
i,slab[i].type,slab[i].posx,slab[i].posy,slab[i].posz,slab[i].flags,slab[i].parent,slab[i].param1,slab[i].param2);
	}

/*
 *  Close the ascii file, and deallocate the memory
 */
	fclose(fp);
//	cfree(atom);
//	cfree(pos);
}
