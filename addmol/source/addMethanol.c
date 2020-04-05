/*
 * AddMethanol merge 4 identical  LxLxL cubic boxes of methanol
 * and remove all  methanols molecules outside
 * to form 1 box with given size Lx,Ly < 2L
 */
#include	<stdio.h>
#include	<stdlib.h>
#include	<typedefs.h>
#include	<globals.h>
#include	<time.h>
#define	abs(x)		((x)>0.?(x):-(x))

main(argc,argv)
int argc;
char *argv[];
{
int new_natoms;		/* new number of atoms in the big box	*/
int i,j,k,l,n,nn;	/* dummy indeces			*/
int	anum, type, flags, parent, param1, param2;
double Lx, Ly, Zshift;
double atof();
tripd tryp;
tripd *new_pos;		/* pointer to new positions array	*/
parts *new_atom;	/* pointer to new atom parts array	*/
char filename[80];	/* input file name			*/
char	sbuf[256];
FILE	*fp;

if (argc != 6) {
	fprintf(stderr, "usage: %s <ascii_file> <binary_file> <Lx> <Ly> <Zshift>\n",
		argv[0]);
	exit(0);
}
Lx = atof(argv[3]);	
Ly = atof(argv[4]);	
Zshift = atof(argv[5]);
/*
 *  open the ascii input file and read the header input
 */
	if ((fp = fopen(argv[1], "r")) == NULL) {
		fprintf(stderr, "writebin: can't open ascii input file\n");
	}

	fgets(sbuf, 256, fp);		/* Input file name -- ignore */

	fgets(sbuf, 256, fp);		/* File type */
	sscanf(sbuf, "%*s %*s %s", filetype);
	if (strcmp(filetype, ".inp") != 0) {
		fprintf(stderr, "writebin: file not of type .inp\n");
		exit(1);
	}

	fgets(sbuf, 256, fp);		/* datestamp -- ignore, will recalc. */

	fgets(sbuf, 256, fp);		/* status -- NEW, MIN, or EQU */
	sscanf(sbuf, "%*s %*s %s", status);

	fgets(sbuf, 256, fp);		/* natoms */
	sscanf(sbuf, "%*s %*s %*s %*s %*s %d", &natoms);

	fgets(sbuf, 256, fp);		/* nsolute */
	sscanf(sbuf, "%*s %*s %*s %*s %*s %*s %d", &nsolute);

	fgets(sbuf, 256, fp);		/* xwall */
	sscanf(sbuf, "%*s %*s %*s %*s %*s %lf", &xwall);

	fgets(sbuf, 256, fp);		/* ywall */
	sscanf(sbuf, "%*s %*s %*s %*s %*s %lf", &ywall);

	fgets(sbuf, 256, fp);		/* zwall */
	sscanf(sbuf, "%*s %*s %*s %*s %*s %lf", &zwall);

	fgets(sbuf, 256, fp);		/* Equil. Temp. */
	sscanf(sbuf, "%*s %*s %*s %*s %lf", &EqTemp);

	fgets(sbuf, 256, fp);		/* St. Dev. of equil. temp. */
	sscanf(sbuf, "%*s %*s %*s %*s %*s %*s %*s %lf", &DEqTemp);

	fgets(sbuf, 256, fp);		/* Equil. Press */
	sscanf(sbuf, "%*s %*s %*s %*s %lf", &EqPress);

	fgets(sbuf, 256, fp);		/* St. Dev. of equil. press. */
	sscanf(sbuf, "%*s %*s %*s %*s %*s %*s %*s %lf", &DEqPress);

	fgets(sbuf, 256, fp);		/* Extra file flag */
	sscanf(sbuf, "%*s %*s %*s %*s %*s %d", &xtrInQ);

/*
 *  allocate space for the arrays and then read them
 */
	if ((atom = (parts *) calloc(natoms, sizeof(parts))) == NULL ||
	    (pos  = (tripd *) calloc(natoms, sizeof(tripd))) == NULL) {
		fprintf(stderr, "writebin: out of memory\n");
		exit(1);
	}

	fgets(sbuf, 256, fp);		/* header line -- ignore */
	for (i=0; i<natoms; i++) {
		fgets(sbuf, 256, fp);	/* pos and atom info */
		sscanf(sbuf, "%d %d %lf %lf %lf %o %d %d %d", &anum, &type,
			&pos[i].fx, &pos[i].fy, &pos[i].fz,
			&flags, &parent, &param1, &param2);
		atom[i].type   = type;
		atom[i].flags  = flags;
		atom[i].parent = parent;
		atom[i].param1 = param1;
		atom[i].param2 = param2;
	}
	fclose(fp);

new_natoms = 4*natoms;
if (  (new_atom	= (parts  *) calloc(new_natoms,sizeof(parts))) == NULL
   || (new_pos	= (tripd  *) calloc(new_natoms,sizeof(tripd))) == NULL)
   {
	fprintf(stderr,"Out of core\n");
	exit(0);
   }
nn = 0;

for (l=0; l < natoms; l = l +3){
for (i=0; i<2; i++){
    for (j=0; j<2; j++){
	tryp.fx =  pos[l].fx + (2*i-1)*xwall;
	tryp.fy =  pos[l].fy + (2*j-1)*ywall;
	if (abs(tryp.fx) < Lx && abs(tryp.fy) < Ly){/*inside the new box*/
		new_pos[nn].fx = tryp.fx;
		new_pos[nn].fy = tryp.fy;
		new_pos[nn].fz = pos[l].fz +zwall+Zshift;
		new_atom[nn].type = atom[l].type;
		new_atom[nn].flags = atom[l].flags;
		new_atom[nn].param1 = atom[l].param1;
		new_atom[nn].param2 = atom[l].param2;
		new_atom[nn].parent = nn;
		nn++;
		for (k = 1; k <= atom[l].param1; k++){
			new_pos[nn].fx = pos[l+k].fx + (2*i-1)*xwall;
			new_pos[nn].fy = pos[l+k].fy + (2*j-1)*ywall;
			new_pos[nn].fz = pos[l+k].fz + zwall+Zshift;
			new_atom[nn].type = atom[l+k].type;
			new_atom[nn].flags = atom[l+k].flags;
			new_atom[nn].param1 = atom[l+k].param1;
			new_atom[nn].param2 = atom[l+k].param2;
			new_atom[nn].parent = nn-k;
			nn++;
		}
	}
     }
}		
}
natoms = nn;
xwall = Lx;
ywall = Ly;
zwall = 100;
/*
 *  Write output to binary input file
 */
	if ((fp = fopen(argv[2], "w")) == NULL) {
		fprintf(stderr, "writebin: can't open binary inp-file\n");
		exit(1);
	}

	fwrite(filetype, sizeof(char), 5, fp);
	time(&datestamp);
	fwrite(&datestamp, sizeof(datestamp), 1, fp);
	fwrite(status, sizeof(char), 4, fp);
	fwrite(&natoms, sizeof(natoms), 1, fp);
	fwrite(&nsolute, sizeof(nsolute), 1, fp);
	fwrite(&xwall, sizeof(xwall), 1, fp);
	fwrite(&ywall, sizeof(ywall), 1, fp);
	fwrite(&zwall, sizeof(zwall), 1, fp);
	fwrite(&EqTemp, sizeof(EqTemp), 1, fp);
	fwrite(&DEqTemp, sizeof(DEqTemp), 1, fp);
	fwrite(&EqPress, sizeof(EqPress), 1, fp);
	fwrite(&DEqPress, sizeof(DEqPress), 1, fp);
	fwrite(new_atom, sizeof(parts), natoms, fp);
	fwrite(&xtrInQ, sizeof(xtrInQ), 1, fp);
	fwrite(new_pos, sizeof(tripd), natoms, fp);

/*
 *  close the binary input file
 */
	fclose(fp);
}
