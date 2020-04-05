/*
 * cub_to_oct  merge 8 cubic copies of solvent to form 1 box with
 * twice the linear size and then remove solvent molecules to keep a
 * truncated octahedral  solvent box with a given size.
 */
#include        <stdio.h>
#include        <typedefs.h>
#include        <globals.h>
#include        <time.h>
#include "atomtypes.h"
#define abs(x)          ((x)>0.?(x):-(x))
double octS;	/*new xwall, ywall and zwall		*/

main(argc,argv)
int argc;
char *argv[];
{
int new_natoms;		/* new number of atoms in the big box	*/
int i,j,k,l,n,nn,ipr;	/* dummy indeces			*/
tripd *new_pos, tpos;		/* pointer to new positions array	*/
parts *new_atom;	/* pointer to new atom parts array	*/
char filename[80];	/* input file name			*/
double atof();
FILE *fp;

if	(argc != 4) {
	fprintf(stderr,"usage: %s  <input binary file> <oct box size> <ascii file>\n", argv[0]);
	exit(0);
}
octS = atof(argv[2]);

/*
 *  open the binary input file
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
	if (nsolute > 0) {
		fprintf(stderr, "input file should be neat solvent\n");
		exit(1);
	}
	if ((atom = (parts *) calloc(natoms, sizeof(atom[0]))) == NULL ||
	    (pos  = (tripd *) calloc(natoms, sizeof(pos[0]))) == NULL) {
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

	
new_natoms = 8*natoms;
if (  (new_atom	= (parts  *) calloc(new_natoms,sizeof(parts))) == NULL
   || (new_pos	= (tripd  *) calloc(new_natoms,sizeof(tripd))) == NULL)
   {
	fprintf(stderr,"Out of core\n");
	exit(0);
   }
nn = 0;

for (i=0; i<2; i++)
    {
    for (j=0; j<2; j++)
    	{
	for (k=0; k<2; k++)
	    {
		for (l=0; l < natoms; l++)
		    {
			if	((atom[l].flags & A_MAJOR) && atom[l].param1){
				tpos.fx = pos[l].fx + (2*i-1)*xwall;
				tpos.fy = pos[l].fy + (2*j-1)*ywall;
				tpos.fz = pos[l].fz + (2*k-1)*zwall;
			   if (inbox(&tpos)){
				new_pos[nn].fx = tpos.fx;
				new_pos[nn].fy = tpos.fy;
				new_pos[nn].fz = tpos.fz;
				new_atom[nn].type = atom[l].type;
				new_atom[nn].flags = atom[l].flags;
				new_atom[nn].param1 = atom[l].param1;
				new_atom[nn].param2 = atom[l].param2;
				new_atom[nn].parent = ipr = nn;
				nn++;
				for (n=l+1; n <= l + atom[l].param1; n++){
				   new_pos[nn].fx = pos[n].fx + (2*i-1)*xwall;
				   new_pos[nn].fy = pos[n].fy + (2*j-1)*ywall;
				   new_pos[nn].fz = pos[n].fz + (2*k-1)*ywall;
				   new_atom[nn].type = atom[n].type;
				   new_atom[nn].flags = atom[n].flags;
				   new_atom[nn].param1 = atom[n].param1;
				   new_atom[nn].param2 = atom[n].param2;
				   new_atom[nn].parent = ipr;
/*if ( nn == ipr+3) new_atom[nn].parent = ipr+2;*/ /*DCE patch*/
				   nn++;
				}
			   }
			}
			else if	(atom[l].flags & A_MAJOR){
				tpos.fx = pos[l].fx + (2*i-1)*xwall;
				tpos.fy = pos[l].fy + (2*j-1)*ywall;
				tpos.fz = pos[l].fz + (2*k-1)*zwall;
			   if (inbox(&tpos)){
				new_pos[nn].fx = tpos.fx;
				new_pos[nn].fy = tpos.fy;
				new_pos[nn].fz = tpos.fz;
				new_atom[nn].type = atom[l].type;
				new_atom[nn].flags = atom[l].flags;
				new_atom[nn].param1 = atom[l].param1;
				new_atom[nn].param2 = atom[l].param2;
				new_atom[nn].parent = nn;
				nn++;
			   }
			}
		    }
            }
	}		
    }
    if (nn > new_natoms){
	fprintf(stderr,"nn = %d > new_natoms = %d\n",nn, new_natoms);
	exit(0);
   }
/*
 *  Write the output to the ascii input-file in argv[3]
 */
	if ((fp = fopen(argv[3], "w")) == NULL) {
		fprintf(stderr, "can't open ascii input-file\n");
		exit(1);
	}

	fprintf(fp,"Input File -- \t%s\n",argv[1]);
	fprintf(fp,"Filetype  -- \t.inp\n");

/*  ctime() inserts a newline character so we don't need to do it */
	fprintf(fp,"Datestamp -- \t%s", ctime(&datestamp));

	fprintf(fp, "Status -- \t%s\n", status);
	fprintf(fp,"Number of atoms (natoms) -- \t%d\n", nn);
	fprintf(fp,"Number of solute atoms (nsolute) -- \t%d\n", nsolute);
	fprintf(fp,"X box size (xwall) -- \t%f\n", octS);
	fprintf(fp,"Y box size (ywall) -- \t%f\n", octS);
	fprintf(fp,"Z box size (zwall) -- \t%f\n", octS);
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
	for (i=0; i<nn; i++) {
		fprintf(fp,
"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
i,new_atom[i].type,new_pos[i].fx,new_pos[i].fy,new_pos[i].fz,new_atom[i].flags,new_atom[i].parent,new_atom[i].param1,new_atom[i].param2);
	}
fprintf(stderr," natoms = %d density = %f\n",nn, nn/(4.0*octS*octS*octS));
/*
 *  Close the ascii file, and deallocate the memory
 */
	fclose(fp);
	cfree(atom);
	cfree(pos);
}
int   inbox(p)
tripdouble *p;
{
double	fx, fy, fz, oct;
int inside;

inside = 1;
  /* Normalizing the positions */
	fx = abs (p->fx / (octS*2));
	fy = abs (p->fy / (octS*2));  
	fz = abs (p->fz / (octS*2));  
	oct = fx + fy + fz ;
	if (fx >= 0.5 || fy >= 0.5 || fz >= 0.5 || oct >= 0.75) inside = 0;
	return(inside);
}
