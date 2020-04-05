/*
 *	This program replace one chloroform molecule near the sn2 reactants
 *  	by one water molecule and write the output as an ascii file
 *
 *	To call:	insertH2O <binary_file> <ascii_file>
 *
 */

#include "typedefs.h"
#include "atomtypes.h"
#include "globals.h"
#include <stdio.h>
#include <sys/time.h>
#define abs(x)          ((x)>0.?(x):-(x))
#define sq(x)   ((x)*(x))
#define hsgn(x)         ((x)<0.?-.5:.5)

main (argc, argv)
	int	argc;
	char	*argv[];
{
	int	i,k,j,m,nCF,clst[260], atype, aflag, pa, p1, p2;
	double r2[260], rC, sqrt();
	tripd r, posO, posH[2];
	FILE	*fp;

/*
 *  argv[1] and argv[2] contain the input and output file names
 */
	if (argc != 3) {
		fprintf(stderr, "usage: %s <binary_file> <ascii_file>\n",
			argv[0]);
		exit(0);
	}

/*
 *  open the binary input file which contains one sn2 system and chloroform molecules
 */
	if ((fp = fopen(argv[1], "r")) == NULL) {
		fprintf(stderr, "can't open binary input-file\n");
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
		fprintf(stderr, "readbin: out of memory\n");
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
 * Find the CHCl3 molecule nearest the sn2 system based on the C-CH3 disatnce
 */

if (nsolute != 3){
	fprintf(stderr,"Wrong input file . Sn2 must have nsolute = 3\n");
	exit(1);
}
offset(natoms-3);/*center the box on the CH3 group*/
nCF = (natoms-3)/5; /* number of chloroform molecules*/
for (i=0;i<nCF;i++){
        r.fx = pos[i*5].fx-pos[natoms-3].fx;
        r.fy = pos[i*5].fy-pos[natoms-3].fy;
        r.fz = pos[i*5].fz-pos[natoms-3].fz;
        r2[i] = r.fx*r.fx + r.fy*r.fy + r.fz*r.fz;/*array of CH3-C distances*/
}               
sort_r(nCF,&r2,&clst);/*clst[0] is the index of the nearest chloroform*/
/*position of water oxygen is the same is the C atom of the removed chloroform*/
posO.fx = pos[clst[0]*5].fx;
posO.fy = pos[clst[0]*5].fy;
posO.fz = pos[clst[0]*5].fz;
/*position of the first water H atoms will be the same as the H of chloroform*/
posH[0].fx = pos[clst[0]*5+4].fx;
posH[0].fy = pos[clst[0]*5+4].fy;
posH[0].fz = pos[clst[0]*5+4].fz;
/*position of the second water H atoms will be 1A along the C-Cl vector */
r.fx = pos[clst[0]*5+1].fx-posO.fx;
r.fy = pos[clst[0]*5+1].fy-posO.fy;
r.fz = pos[clst[0]*5+1].fz-posO.fz;
rC = sqrt(r.fx*r.fx + r.fy*r.fy + r.fz*r.fz);/*C-Cl bond length*/
posH[1].fx = posO.fx + r.fx/rC;
posH[1].fy = posO.fy + r.fy/rC;
posH[1].fz = posO.fz + r.fz/rC;

/*
 *  Write the output to the ascii input-file in argv[2].
 *  Note that the number of atoms is 2 less than the original file
 */
	if ((fp = fopen(argv[2], "w")) == NULL) {
		fprintf(stderr, "readbin: can't open ascii input-file\n");
		exit(1);
	}

	fprintf(fp,"Input File -- \t%s\n",argv[1]);
	fprintf(fp,"Filetype  -- \t.inp\n");

/*  ctime() inserts a newline character so we don't need to do it */
	fprintf(fp,"Datestamp -- \t%s", ctime(&datestamp));

	fprintf(fp, "Status -- \t%s\n", status);
	fprintf(fp,"Number of atoms (natoms) -- \t%d\n", natoms-2);
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
/*print the water oxygen*/
i =0; atype = 0 ; aflag = 0026; pa = 0; p1 = 2; p2 = 3;
fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
i,atype,posO.fx,posO.fy,posO.fz,aflag,pa,p1,p2);
/*print the water hydrogen*/
i++; atype = 1 ; aflag = 0004; pa = 0; p1 = 0; p2 = 0;
fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
i,atype,posH[0].fx,posH[0].fy,posH[0].fz,aflag,pa,p1,p2);
i++;
fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
i,atype,posH[1].fx,posH[1].fy,posH[1].fz,aflag,pa,p1,p2);
/*print all the chloroform except for the one replaced by water*/
for (j=0; j<nCF; j++) {
    if (j != clst[0]){
	for (k=0;k<5;k++){
  	    i++;
     	    m = 5*j+k;
	    fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
		i,atom[m].type,pos[m].fx,pos[m].fy,pos[m].fz,atom[m].flags,i-k,atom[m].param1,atom[m].param2);
	}
    }
}
/*print the solute atoms*/
for (k=0;k<3;k++){
	i++;
	m = natoms - 3 +k;
	fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
	i,atom[m].type,pos[m].fx,pos[m].fy,pos[m].fz,atom[m].flags,i,atom[m].param1,atom[m].param2);
}
/*
 *  Close the ascii file, and deallocate the memory
 */
	fclose(fp);
	cfree(atom);
	cfree(pos);
}
sort_r(n,r2,clst)
int n;
double *r2;
int *clst;
{
int i,j,b;
double a;
	for (j=0;j<n;j++) clst[j] = j;
	for (j=1;j<n;j++) {
		a = r2[j];
		b = clst[j];
		i=j-1;
		while (i >= 0 && r2[i] > a) {
			r2[i+1] = r2[i];
			clst[i+1]= clst[i];
			i--;
		}
		r2[i+1]=a;
		clst[i+1]=b;
	}
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
	if (n > natoms-1 || n < 0) return;
	off.fx = pos[n].fx;
	off.fy = pos[n].fy;
	off.fz = pos[n].fz;

	for (i = 0; i < natoms; i++) {
		pos[i].fx -= off.fx;
		pos[i].fy -= off.fy;
		pos[i].fz -= off.fz;

		if ((atom[i].flags & A_MAJOR) && atom[i].param1) {
			image.fx = -pos[i].fx;
			image.fy = -pos[i].fy;
			image.fz = -pos[i].fz;
			mvimage(&pos[i]);
			image.fx += pos[i].fx;
			image.fy += pos[i].fy;
			image.fz += pos[i].fz;
			for (j=i+1; j <= i + atom[i].param1; j++) {
				pos[j].fx += image.fx;
				pos[j].fy += image.fy;
				pos[j].fz += image.fz;
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

