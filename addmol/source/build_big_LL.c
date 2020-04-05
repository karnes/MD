/*
 * Build_big_LL merge 4 boxes of molecular water/liquid interface to form
 * 1 box with four times the cross section, but the same Zwall
 * It is assumed that the input files are identical except for the positions of
 * the atoms. The name of the 4 input files must end with 00n.equ, n=0,1...,3.
 * in the command line do not supply the 00n.equ part
 */

#include "typedefs.h"
#include "atomtypes.h"
#include "globals.h"
#include <stdio.h>
#include <time.h>

main(argc,argv)
int argc;
char *argv[];
{
int new_natoms;		/* new number of atoms in the big box	*/
int nw;		/* the number of water atoms in the small box	*/
int i,j,k,l,n,nn,no,pnum;	/* dummy indeces			*/
tripd *new_pos;		/* pointer to new positions array	*/
parts *new_atom;	/* pointer to new atom parts array	*/
char filename[80];	/* input file name			*/

if	(argc != 4) {
	fprintf(stderr,"usage: %s  <natoms> <# water molecules> <input root file name>\n", argv[0]);
	exit(0);
}
	
new_natoms = 4*atoi(argv[1]);
nw = 3*atoi(argv[2]);
if (  (new_atom	= (parts  *) calloc(new_natoms,sizeof(parts))) == NULL
   || (new_pos	= (tripd  *) calloc(new_natoms,sizeof(tripd))) == NULL)
   {
	fprintf(stderr,"Out of core\n");
	exit(0);
   }
nn = 0;
no = nw*4;

for (i=0; i<2; i++)
    {
    for (j=0; j<2; j++)
    	{
		n = 2*i + j;
		sprintf(filename, "%s%03d.equ", argv[3],n);
		readpos(filename);
		for (l=0; l < nw; l++)/* read water and place in new pos array*/
		    {
			new_pos[nn].fx = pos[l].fx + (2*i-1)*xwall;
			new_pos[nn].fy = pos[l].fy + (2*j-1)*ywall;
			new_pos[nn].fz = pos[l].fz ;
			new_atom[nn].type = atom[l].type;
			new_atom[nn].flags = atom[l].flags;
			new_atom[nn].param1 = atom[l].param1;
			new_atom[nn].param2 = atom[l].param2;
			pnum = nn/3;			
			new_atom[nn].parent = 3*pnum;
			nn++;
		    }
		/*read organic liquid and place in the new_pos array*/
		for (l=nw; l < natoms; l++)
		    {
			new_pos[no].fx = pos[l].fx + (2*i-1)*xwall;
			new_pos[no].fy = pos[l].fy + (2*j-1)*ywall;
			new_pos[no].fz = pos[l].fz ;
			new_atom[no].type = atom[l].type;
			new_atom[no].flags = atom[l].flags;
			new_atom[no].param1 = atom[l].param1;
			new_atom[no].param2 = atom[l].param2;
			pnum = (no-4*nw)/atom[nw].param2;			
			new_atom[no].parent = 4*nw+atom[nw].param2*pnum;
			no++;
		    }
            }
	}		
if (no != new_natoms)
   {
	fprintf(stderr,"no = %d different than new_natoms = %d \n",no,new_natoms);
	exit(0);
   }
sprintf(status, "NEW");
sprintf(filename, "%s%03d.new.big", argv[3],0);
writepos(filename,new_atom, new_pos);
}

/* write the new positions */

writepos(file, new_atom, new_pos)
char	*file;
parts *new_atom;
tripd *new_pos;
{
	char	filetype[5];
	FILE	*fp;

/*	Get the time of day to stamp file ...				*/

	time(&datestamp);

/*	Open file and give error message if there is a problem ...	*/

	if	((fp = fopen(file,"w")) == NULL) {
		fprintf(stderr,"writepos:  can't open %s\n",file);
		exit(1);
	}

/*	Write the pertinent stuff ...					*/

	sprintf(filetype, ".inp");
	fwrite(filetype, sizeof(char), 5, fp);
	fwrite(&datestamp, sizeof(datestamp), 1, fp);
	fwrite(status, sizeof(char), 4, fp);	/* must define status */
	natoms *= 4;
	fwrite(&natoms, sizeof(natoms), 1, fp);
	fwrite(&nsolute, sizeof(nsolute), 1, fp); /* nsolute = 0 */
	xwall *= 2;
	ywall *= 2;
	fwrite(&xwall, sizeof(xwall), 1, fp);
	fwrite(&ywall, sizeof(ywall), 1, fp);
	fwrite(&zwall, sizeof(zwall), 1, fp);
	EqTemp = DEqTemp = EqPress = DEqPress = 0.0;
	fwrite(&EqTemp, sizeof(EqTemp), 1, fp);
	fwrite(&DEqTemp, sizeof(DEqTemp), 1, fp);
	fwrite(&EqPress, sizeof(EqPress), 1, fp);
	fwrite(&DEqPress, sizeof(DEqPress), 1, fp);
	fwrite(new_atom,sizeof new_atom[0],natoms,fp);
	xtrInQ = 0;
	fwrite(&xtrInQ, sizeof(xtrInQ), 1, fp);
	fwrite(new_pos, sizeof new_pos[0], natoms, fp);

/*	All done ...							*/

	fclose(fp);
}			

/* read the old positions */

readpos(filename)
char filename[80];
{
static int	allocate = 1;
FILE		*fp;

/*	Open input file ...						*/

	if	((fp = fopen(filename,"r")) == NULL) {
		fprintf(stderr,"Can't open %s\n",filename);
		exit(1);
	}

/*	Read filetype string and warn if not ".inp" ...		*/

	fread(filetype, sizeof(char), 5, fp);
	if (strcmp(filetype, ".inp") != 0)
		fprintf(stderr, "rinp:  file is not a .inp file");

/*	Read header information ...					*/

	fread(&datestamp, sizeof(datestamp), 1, fp);
	fread(status, sizeof(char), 4, fp);	/* status = NEW, MIN, EQU */
	fread(&natoms, sizeof(natoms), 1, fp);
	fread(&nsolute, sizeof(nsolute), 1, fp);
	if (nsolute != 0)
	   fprintf(stderr,"Warning: %s contains solute\n", filename);
	fread(&xwall, sizeof(xwall), 1, fp);
	fread(&ywall, sizeof(ywall), 1, fp);
	fread(&zwall, sizeof(zwall), 1, fp);
	fread(&EqTemp, sizeof(EqTemp), 1, fp);	/* if status = EQU */
	fread(&DEqTemp, sizeof(DEqTemp), 1, fp); 	/* '' */
	fread(&EqPress, sizeof(EqPress), 1, fp);
	fread(&DEqPress, sizeof(DEqPress), 1, fp);
		
/*	Allocate space for arrays, if this is the first time through ...*/

	if (allocate &&
	(  (atom	= (parts  *) calloc(natoms,sizeof(parts))) == NULL
	|| (pos		= (tripd  *) calloc(natoms,sizeof(tripd))) == NULL
	)) {
		fprintf(stderr,"Out of core\n");
		exit(0);
	}
	allocate = 0;	/* so we won't allocate space the next time */

/*	Read the atom specifications and the positions ...		*/

	fread(atom,sizeof atom[0],natoms,fp);
	fread(&xtrInQ, sizeof(xtrInQ), 1, fp);
	fread(pos, sizeof pos[0], natoms, fp);

	fclose(fp);
}

