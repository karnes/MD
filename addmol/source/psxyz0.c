/*
 *	psxyz.c:	This program prints an xyz file from a binary
 *			input file. If there is a solute, it prints
 *			only the solute and all the solvents inside
 *			a sphere of radius R
 *
 *	To call:	psxyz <binary_file> <xyz file> <R>
 *
 */

#include	<stdio.h>
#include	<stdlib.h>
#include	<typedefs.h>
#include	<atomtypes.h>
#include	<globals.h>
#include	<time.h>
#define BIG  1000
#define abs(x)          ((x)>0.?(x):-(x))
#define sq(x)   ((x)*(x))
#define hsgn(x)         ((x)<0.?-.5:.5)
char *symb[] = {"O","H","H","H","H","H","H","H","H","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","O","O","O","O","O","O","N","N","N","N","N","N","N","N","N","N","S","S","P","Li","Na","K","Cs","F","Cl","F","Cl","Cl","O","C","He","Ar","Kr","Xe","N","D","I","C","H","I","C","H","Cl","C","H","Cl","C","H","Br","C","H","Br","Cl","C","Cl","Ne","X","X","X","X","X","X","X","X","X","X","X","X","X","X","Si","X"};

main (argc, argv)
	int	argc;
	char	*argv[];
{
	int	i,j,n;
	double R, r2, atof();
	FILE	*fp;

/*
 *  argv[1] and argv[2] contain the input and output file names
 */
	if (argc < 3 || argc > 4) {
		fprintf(stderr, "usage: %s <binary_file> <ascii_file> <opt R>\n",
			argv[0]);
		exit(0);
	}

/*
 *  open the binary input file
 */
	if ((fp = fopen(argv[1], "r")) == NULL) {
		fprintf(stderr, "psxyz: can't open binary input-file\n");
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
		fprintf(stderr, "psxyz: out of memory\n");
		exit(1);
	}

	fread(atom, sizeof(parts), natoms, fp);
	fread(&xtrInQ, sizeof(xtrInQ), 1, fp);
	fread(pos, sizeof(tripd), natoms, fp);

/*
 *  Close the binary file
 */
	fclose(fp);
R = BIG;
/*
 *  Write the output to the ascii input-file in argv[2]
 */
	if ((fp = fopen(argv[2], "w")) == NULL) {
		fprintf(stderr, "psxyz: can't open ascii input-file\n");
		exit(1);
	}
	fprintf(fp,"                       \n");

//if (nsolute == 1 && argc == 4){
//	 R = atof(argv[3]);
//	offset(natoms-1);/*center the box on the ion*/
//}
//if (nsolute != 1 && argc == 4) R = atof(argv[3]);

n = 0;
offset(n);

for	(i = 0; i < natoms; i++){
//printf("flag = %d\n",atom[i].flags);
//	if	((atom[i].flags & A_MAJOR) && atom[i].param1){
//		r2 = sq(pos[i].fx)+sq(pos[i].fy)+sq(pos[i].fz);
//		if (1){//r2 < R*R){
//			for 	(j=i; j <= i + atom[i].param1; j++){
//				n++;
				fprintf(fp,"%s %f %f %f\n", symb[atom[i].type],
					pos[i].fx,pos[i].fy,pos[i].fz);
//			}
//		}
//	}
//	else if	(atom[i].flags & A_MAJOR){
//		r2 = sq(pos[i].fx)+sq(pos[i].fy)+sq(pos[i].fz);
//		if (1){//r2 < R*R){
//			n++;
//			fprintf(fp,"%s %f %f %f\n", symb[atom[i].type],
//				pos[i].fx,pos[i].fy,pos[i].fz);
//		}
//	}
}

//fprintf(fp,"move this line to the top and remove the text %d\n",n);
/*
 *  Close the ascii file, and deallocate the memory
 */
	fclose(fp);
	if ((fp = fopen(argv[2], "r+")) == NULL) {
		fprintf(stderr, "psxyz: can't open ascii input-file\n");
		exit(1);
	}
	fprintf(fp,"%d\n",natoms);
	char str[60]; fgets(str, 60,fp);
	fclose(fp);
	cfree(atom);
	cfree(pos);
}
/*
 *	offset:
 *	
 *	This routine offsets all atoms by an appropriate amount so as to bring
 *	atom n to the center of the box (0,0,0).
 */

offset(n)
	int	n;	/* atom number from which to offset positions */
{
	int	i, j;
	tripd	image, off;
//	if (n > natoms-1 || n < 0) return;
	off.fx = pos[n].fx;
	off.fy = pos[n].fy;
	off.fz = 0;//pos[n].fz;

	for (i = 0; i < natoms; i++) {
		pos[i].fx -= off.fx;
		pos[i].fy -= off.fy;
//		pos[i].fz -= off.fz;

		if ((atom[i].flags & A_MAJOR) && atom[i].param1) {
			image.fx = -pos[i].fx;
			image.fy = -pos[i].fy;
//			image.fz = -pos[i].fz;
			mvimage(&pos[i]);
			image.fx += pos[i].fx;
			image.fy += pos[i].fy;
//			image.fz += pos[i].fz;
			for (j=i+1; j <= i + atom[i].param1; j++) {
				pos[j].fx += image.fx;
				pos[j].fy += image.fy;
//				pos[j].fz += image.fz;
			}
		}
		else if (atom[i].flags & A_MAJOR)
			mvimage(&pos[i]);
	}
}

/*
 *	This routine does the periodic imaging of any given
 *	vector for the parallelopiped or the truncated octahedral
 *      boundary conditions.
 */

mvimage(del)
	tripdouble	*del;
{
	double	fx, fy, fz;
	tripd period, iperiod;
	if (natoms == nsolute) /* no solvent */
	     return;
	period.fx = 2.0 * xwall;	/* box period (boundary conditions) */
	period.fy = 2.0 * ywall;	/* box period (boundary conditions) */
	period.fz = 2.0 * zwall;	/* box period (boundary conditions) */
	iperiod.fx = 1.0 / period.fx;		/* inverse period */
	iperiod.fy = 1.0 / period.fy;		/* inverse period */
	iperiod.fz = 1.0 / period.fz;		/* inverse period */
	fx = del->fx * iperiod.fx;  /* Normalizing the positions */
	fy = del->fy * iperiod.fy;  /* |fx|, |fy|, |fz| <= 0.5 for the  */
	fz = del->fz * iperiod.fz;  /* atom to be in the box */
	fx -= (int) fx;
	fy -= (int) fy;   /* bringing the particle to the nearest box */
	fz -= (int) fz;
	fx -= (int) (2 * fx);
	fy -= (int) (2 * fy);  /* bringing the particle to the box */
	fz -= (int) (2 * fz);
	/* truncated octahedral boundaries */
		if	((abs(fx) + abs(fy) + abs(fz)) > 0.75) {
			fx -= hsgn(fx);
			fy -= hsgn(fy);
			fz -= hsgn(fz);
		}

	del->fx = fx * period.fx;
	del->fy = fy * period.fy;
	del->fz = fz * period.fz;
}
