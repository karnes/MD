/*
 *	addmol:	This program sets up the initial condition "input files" to
 *		start the MD procedure, usually beginning with minimization.
 *		The program runs interactively and the following information
 *		is prompted for:
 *			type of molecule
 *			number of molecules
 *			density
 *			file name for output
 *		The output file can be either a new file, or you can append
 *		more atoms to the end of another file.
 *
 *		The output from addmol consists of a binary file containing
 *		the following information:
 *			filetype (".inp")
 *			datestamp
 *			status of file ("NEW" -- to indicate just generated)
 *			number of atoms (not molecules!)
 *			number of solute atoms
 *			dimension of the periodic box in x
 *			dimension of the periodic box in y
 *			dimension of the periodic box in z
 *			equilibration temp (set to zero here)
 *			s.d. of equil. temp. (set to zero here)
 *			equilibration pressure (set to zero here)
 *			s.d. of equil. press. (set to zero here)
 *			for each atom:
 *				parameters describing the atom
 *			flag determining if there is another extra input file
 *			for each atom:
 *				the x,y,z coordinates of the atom
 *
 *		To call addmol, simply type "addmol" and the program will
 *		prompt you for the appropriate input.
 *
 *		To build the usual case of a few atom solute in the midst of
 *		a many atom solvent, do the process in two steps.  First, call
 *		addmol to create a new input file containing the solvent atoms.
 *		Then call addmol again, appending to the above file the solute
 *		atoms.
 *
 *		To check that the binary files you have just created are not
 *		grossly in error, you may check them using readbin.  This
 *		program translates the binary file into ascii format, where it
 *		can be read and edited.  Call writebin to convert back into
 *		binary, if you have made any changes in the ascii file.
 */



#include	<typedefs.h>	/* global type definitions */
#include	<atomtypes.h>	/* defines for the atom types and masses */
#include	<globals.h>	/* extern declarations for global variables */

#include	<stdio.h>
#include	<math.h>
#include	<sys/time.h>

#define	RMAX	2147483648.0	/* 2^31 -- random returns in range [0,RMAX-1] */
#define	CUBE_WALL  0.20757085	/* 1 / (Avagadro's # * 8e-24 cubic cm/A)
				 * factor of 8 accounts for cube volume
				 */
#define	PTO_WALL   0.41514171	/* 1 / (Avagadro's # * 4e-24 cubic cm/A)
				 * factor of 4 accounts for t.o. volume
				 */
#define	CUBE_RT_2	1.58739	/* cube root of 2 (SWC im_space correction) */
#define	IMSPACE_CONST	0.7937	 /* intermolecular spacing constant, is
				 * actually cube root of 0.5, corresponding
				 * to a calculation of the volume of n cubes
				 * containing the inscribed sphere of radius
				 * "im_space"
				 */
#define	MAX_FACTOR	1000000	/* # of attempts allowed per molecule added */
#define	MAX_MOLNAME	16
#define sign(x)		(x < 0.0 ? -1.0 : 1.0)

#define FUNITS		1024	/* "globals.c" units per angstrom */

#define	SWC	1
#define	PC	2
#define	PTO	3
#define	POR	4		/* Periodic ORthorhombus */

static char	swc_str[] = "soft walled cubic";
static char	pc_str[] = "periodic cubic";
static char	pto_str[] = "periodic truncated octahedral";
static char	por_str[] = "periodic orthorhombus";

typedef struct	{
	int	bt_cond;	/* Boundary Type - CONDitions */
	char *	bt_str;		/* Boundary Type - STRing */
	float	bt_wall;	/* Boundary Type - WALL constant */
	}	bc_type;

typedef struct	{
	char *	bw_mand;	/* Boundary Word - MANDatory part */
	char *	bw_opt;		/* Boundary Word - OPTional part */
	bc_type	*bw_ptr;	/* Boundary Word - PoinTeR to bc info */
	}	bc_word;

static bc_type	swc	=	{SWC,	(char *) swc_str,	CUBE_WALL};
static bc_type	pc	=	{PC,	(char *) pc_str,	CUBE_WALL};
static bc_type	pto	=	{PTO,	(char *) pto_str,	PTO_WALL};
static bc_type	por	=	{POR,	(char *) por_str,	CUBE_WALL};
	
static bc_word	bound_cond[] = {
	"s",	"wc",	(bc_type *) &swc,
	"pc",	0,	(bc_type *) &pc,
	"pt",	"o",	(bc_type *) &pto,
	"por",	0,	(bc_type *) &por,
	0,	0,	0};

/******************************************************************************/

main(argc,argv)
	int	argc;
	char	*argv[];
{
	float		getmass();
	FILE		*fp;		/* file pointer */
	int		i, j, k;	/* misc indices */
	char		inpname[60];	/* stores input-file name */
	tripfloat	im_space;	/* intermolecular spacing in angstroms
					 * one each for x,y,z */
	float		imspace_const;	/* spacing constant for bc's used */
	float		tot_mass;	/* total mass of system */
	float		addmass;	/* tot mass of to-be-added molecules */
	float		density;	/* desired density of system */
	int		nmoladd;	/* number of molecules to add */
	int		nadd;		/* number of atoms to add */
	int		nput;		/* number of molecules actually placed*/
	int		ninmol;		/* number of atoms in molecule */
	int		attempts;	/* number of attempts while placing */
	float		wall_const;	/* density constant according to
						boundary conditions used */

	struct timeval	timep;		/* used to initialize random() */
	struct timezone	timezp;

	int		newfile;
	tripfloat	new;
	bc_type		*bc_ptr;
	bc_type		*bc_match();
	aword	 	*mwdptr;	/* molecule word pointer */
	aword		*mol_match();
	description	*molptr;	/* molecule description pointer */
	description	*mp;
	char		molname[15];
	char		answer[50];
	float		max;

/*  Initialize some variables */
	newfile = nput = attempts = 0;
	tot_mass = 0.0;

/*  Seed the random number generator and warm it up a bit */
	gettimeofday(&timep, &timezp);
	srandom((int)(timep.tv_sec + timep.tv_usec));
	for (i=0; i<100; i++)
		random();

/*  Get the molecule information from the user */
	printf("Molecule to be added:  ");
	while ((mwdptr = mol_match(gets(molname),molecules)) == NULL) {
		printf("Valid responses are:\n\t");
		for (i = 0; molecules[i].w_mandatory != NULL; i++) {
			sprintf(molname,"%s%s",molecules[i].w_mandatory,
				molecules[i].w_optional);
			printf("%s",molname);
			j = strlen(molname);
			if ((i+1) % 4)
				for (k = MAX_MOLNAME - j; k > 0; k -= 8)
					putchar('\t');
			else {
				putchar('\n');
				putchar('\t');
			}	/* end else */
		}	/* end for */
		printf("\nChoice:  ");
	}	/* end while */
	strcpy(molname,mwdptr->w_mandatory);
	strcat(molname,mwdptr->w_optional);
	molptr = mwdptr->w_pnt;

/*  Calculate the total mass of 1 molecule */
	ninmol = molptr->d_atom.param1 + 1;
	for	(mp = molptr; mp->d_atom.type != NOTANATOM; mp++)
		tot_mass += getmass(mp->d_atom.type);

/*
 *  Get the number of molecules from the user and find total system mass and
 *  the number of atoms
 */
	do {
		printf("Number of %ss to add:  ", molname);
		gets(answer);
		sscanf(answer,"%d",&nmoladd);
	} while	(nmoladd <= 0);
	tot_mass *= nmoladd;
	addmass = tot_mass;
	nadd = nmoladd * ninmol;

/*
 *  Get the boundary conditions, and from that determine the wall and inter-
 *  molecular spacing constants
 */
	printf("Boundary conditions:  ");
	while ((bc_ptr = bc_match(gets(answer),bound_cond)) == NULL) {
		printf("Valid responses are:\n");
		for	(i = 0; bound_cond[i].bw_mand != NULL; i++) {
			printf("\t%s",bound_cond[i].bw_mand);
			if	(bound_cond[i].bw_opt != NULL)
				printf("%s",bound_cond[i].bw_opt);
			printf("\t(%s)\n",bound_cond[i].bw_ptr->bt_str);
		}	/* end for */
		printf("Choice:\t");
	}	/* end while */
	wall_const = bc_ptr->bt_wall;
	switch	(bc_ptr->bt_cond) {
		case	PTO:
			imspace_const = IMSPACE_CONST * 1.26;
			break;
		case	SWC:
			imspace_const = IMSPACE_CONST * CUBE_RT_2;
			break;
		case	PC:
			imspace_const = IMSPACE_CONST * CUBE_RT_2 * 0.95;
			break;
		case	POR:
			imspace_const = IMSPACE_CONST * CUBE_RT_2 * 0.95;
			goto cont;
			break;
	}

/*  Find the density */
	do {
		printf("Density (g/cc):  ");
		gets(answer);
		sscanf(answer,"%f",&density);
	} while	(density <= 0.0);

/*  Prompt user whether to add results to previous file or create new one */
cont:	do {
		printf("Add them to a previous file?  ");
		gets(answer);
	} while	(answer[0] != 'y' && answer[0] != 'n');

/*
 *  Adding results to previous file.  Get file name, and read the information
 *  already there
 */
	if (answer[0] == 'y') {
		do {
			printf("Previous file:  ");
			gets(answer);
			sscanf(answer,"%s",inpname);
		} while	((fp = fopen(inpname,"r")) == NULL);

		fread(filetype, sizeof(char), 5, fp);
		if (strcmp(filetype, ".inp") != 0) {
			fprintf(stderr, "addmol: file not of type .inp\n");
			exit(1);
		}
		fread(&datestamp, sizeof(datestamp), 1, fp);
		fread(status, sizeof(char), 4, fp);
		if (strcmp(status, "NEW") != 0)
			fprintf(stderr, "addmol: file not of status NEW\n");
		fread(&natoms, sizeof(natoms), 1, fp);
		fread(&nsolute, sizeof(nsolute), 1, fp);
		fread(&xwall, sizeof(xwall), 1, fp);
		fread(&ywall, sizeof(ywall), 1, fp);
		fread(&zwall, sizeof(zwall), 1, fp);
		fread(&EqTemp, sizeof(EqTemp), 1, fp);
		fread(&DEqTemp, sizeof(DEqTemp), 1, fp);
		fread(&EqPress, sizeof(EqPress), 1, fp);
		fread(&DEqPress, sizeof(DEqPress), 1, fp);

		if ((pos = (tripd *) calloc(natoms+nadd, sizeof(pos[0])))
			== NULL ||
		   (atom = (parts *) calloc(natoms+nadd, sizeof(atom[0])))
			== NULL) {
			fprintf(stderr, "addmol: out of core\n");
			exit(1);
		}
		fread(atom, sizeof(atom[0]), natoms, fp);
		fread(&xtrInQ, sizeof(xtrInQ), 1, fp);
		fread(pos, sizeof(pos[0]), natoms, fp);

/*
 *  Determine the total number of molecules in the old data (stored in the
 *  variable k in the following), and add the old masses to the previously
 *  determined new masses to get the overall total mass
 */
		for (k = 0,j = 0; j < natoms; j++) {
			if (A_MAJOR & atom[j].flags)
				k++;
			tot_mass += getmass(atom[j].type);
		}

/*
 *  Determine the new wall sizes and density and tell the user.  Calculate the
 *  new intermolecular space constant (for POR and non-POR bc's)
 */
		if (bc_ptr->bt_cond != POR) {
			printf("Old wall was %8.5f angstroms\n",xwall);
			if ((xwall = pow((wall_const * tot_mass / density),
					1.0/3.0)) < 5.0) {
				xwall = 5.0;
				printf("Using a wall of 5.0 A\n");
				printf("New density is %8.5f g/cc\n",
				(density = tot_mass * wall_const * 0.008));
			}
			else
				printf("New wall is %8.5f angstroms\n",xwall);
			ywall = zwall = xwall;
			im_space.fx = imspace_const * xwall / pow((float) (nmoladd + k), 1.0/3.0);
			printf("Intermolecular spacing is %8.5f A\n",im_space.fx);
			im_space.fy = im_space.fz = im_space.fx;
		}
		else {			/* POR */
			printf("Old xwall, ywall, zwall were %f %f %f A\n",
				xwall, ywall, zwall);
			do {
				printf("Do you want to change them?");
				gets(answer);
			} while	(answer[0] != 'y' && answer[0] != 'n');
			if (answer[0] == 'y') {
				do {
					printf("Enter new xwall:\n");
					gets(answer);
					sscanf(answer,"%lf",&xwall);
					printf("Enter new ywall:\n");
					gets(answer);
					sscanf(answer,"%lf",&ywall);
					printf("Enter new zwall:\n");
					gets(answer);
					sscanf(answer,"%lf",&zwall);
				} while (xwall < 5.0 || ywall < 5.0 ||
						zwall < 5.0);
			}
			printf("New density is %8.5f g/cc\n",
			(density = tot_mass * wall_const /(xwall*ywall*zwall)));
			im_space.fx = imspace_const * xwall / pow((float) (nmoladd + k), 1.0/3.0);
			im_space.fy = imspace_const * ywall / pow((float) (nmoladd + k), 1.0/3.0);
			im_space.fz = imspace_const * zwall / pow((float) (nmoladd + k), 1.0/3.0);
			printf("Intermolecular spacing is %8.5f A\n",(max = (im_space.fx + im_space.fy + im_space.fz) / 3.0));	/* average of them */
		}

/*  Close the old input-file */
		fclose(fp);
	}


/*  Else, if we are starting a new input-file . . .
 *  Calculate the wall sizes and im spacings
 */
	else {
		natoms = 0;
		newfile = 1;
		if (bc_ptr->bt_cond != POR) {
			if ((xwall = pow((wall_const * tot_mass / density),
					1./3.)) < 5.0) {
				xwall = 5.0;
				printf("Using a wall of 5.0 A\n");
				printf("New density is %8.5f g/cc\n",
				(density = tot_mass * wall_const * 0.008));
			}
			else
				printf("Wall is %8.5f A\n",xwall);
			im_space.fx = imspace_const * xwall / 
					pow((float) nmoladd, 1.0/3.0);
			im_space.fy = im_space.fz = im_space.fx;
			printf("Intermolecular spacing is %8.5f A\n",
				im_space.fx);
			ywall = zwall = xwall;
		}
		else {		/* POR */
			do {
				printf("Enter new xwall, ywall, zwall:\n");
				gets(answer);
				sscanf(answer,"%lf",&xwall);
				gets(answer);
				sscanf(answer,"%lf",&ywall);
				gets(answer);
				sscanf(answer,"%lf",&zwall);
			} while	(xwall < 5.0 || ywall < 5.0 || zwall < 5.0);
			printf("New density is %8.5f g/cc\n",
			(density = tot_mass * wall_const /(xwall*ywall*zwall)));

			im_space.fx = imspace_const * xwall /
					pow((float) nmoladd, 1.0/3.0);
			im_space.fy = imspace_const * ywall /
					pow((float) nmoladd, 1.0/3.0);
			im_space.fz = imspace_const * zwall /
					pow((float) nmoladd, 1.0/3.0);
			printf("Intermolecular spacing is %8.5f A\n",(max = (im_space.fx + im_space.fy + im_space.fz) / 3.0));	/* average of them */
		}

/*  Allocate space for the position and atom arrays */
		if ((pos = (tripd *) calloc(nadd, sizeof(pos[0])))
				== NULL ||
		   (atom = (parts *) calloc(nadd, sizeof(atom[0])))
				== NULL) {
			fprintf(stderr,"%s: out of core\n",argv[0]);
			exit(1);
		}
	}

/*
 *  For both new input-files and revisions of old ones . . .
 *  Prompt for output file name
 */
	do {
		printf("Output file name:  ");
		gets(answer);
		sscanf(answer,"%s",inpname);
	} while	((fp = fopen(inpname,"w")) == NULL) ;

/*
 *  Start putting atoms in the box.  If only one, put at origin.
 *  Otherwise, call findpt to find a suitable spot, then call addmol to place
 *  it there
 */
	if	(newfile && nmoladd == 1) {
		new.fx = new.fy = new.fz = 0.0;
		i = addmol(molptr,&new,0) + 1;
		nput = attempts = 1;
	}

	else {
		for (i = 0; i < nadd; i++) {
			int	ntrys;
	 
			if ((ntrys = find_pt(&new,&im_space,i,bc_ptr->bt_cond))
					< 0) {
				printf("Failed to find suitable point in %d attempts\n",
					MAX_FACTOR * 10);
				attempts = ntrys = MAX_FACTOR * 10;
				break;
			}

			i = addmol(molptr,&new,i);
			++nput;
			if ((attempts += ntrys) > MAX_FACTOR * nmoladd) {
				printf("Exceeded total maximum allowed trys.\n");
				i++;
				break;
			}
		}
	}

/*  Give statistics on the success of the building process */
	printf("Attempts:\t%d\nSuccesses:\t%d\nSuccess rate\t%4.1f%%\n",
		attempts,nput,100. * (float) nput / (float) attempts);
	if (nput != nmoladd && nmoladd)
		printf("New density = %8.5f g/cc\n",
			density * ((addmass * (float) nput / (float) nmoladd) +
			tot_mass - addmass) / tot_mass);

/*  Update natoms to the number actually entered in the box */
	natoms += i;

/*  Write the data to the input file, and close it up */
	time(&datestamp);	/* update the datestamp */
	EqTemp = DEqTemp = 0.0;	/* zero these fields */
	EqPress = DEqPress = 0.0; /* ditto */
	xtrInQ = 0;		/* no extra file yet */
	nsolute = 0;		/* has not been implemented yet */
	sprintf(status, "NEW");	/* set status to NEW */
	wfile(fp);
	fclose(fp);
}

/*****************************************************************************/
/*
 *  addmol:  Add the molecule, rotated to a random orientation, at the given
 *	     point
 */

addmol(molecule,point,next)
	description	*molecule;	/* pointer to a molecule description */
	tripfloat	*point;	/* pt. about which the molecule is to be built*/
	int		next;	/* offset to next opening in pos/atom arrays */
{
	tripfloat		*pptr;
	description		*dp;
	tripfloat		*molpos;
	int			index;
	
	if ((molpos = (tripfloat *)calloc(molecule->d_atom.param1+1,
				sizeof(tripfloat))) == NULL) {
		fprintf(stderr,"function addmol: out of core\n");
		exit(1);
		}
	for (dp = molecule, pptr = molpos; dp->d_atom.type != NOTANATOM;
				dp++, pptr++) {
		pptr->fx = (float) dp->d_triplet.x;
		pptr->fy = (float) dp->d_triplet.y;
		pptr->fz = (float) dp->d_triplet.z;
	}
	rotcoord(molpos,molecule->d_atom.param1+1);
	dp = molecule;
	pptr = molpos;
	index = natoms + next;
	do {
		pos[index].fx = point->fx + pptr->fx / FUNITS;
		pos[index].fy = point->fy + pptr->fy / FUNITS;
		pos[index].fz = point->fz + pptr->fz / FUNITS;
		atom[index].type = dp->d_atom.type;
		atom[index].flags = dp->d_atom.flags;
		if	(atom[index].flags & A_RING)
			atom[index].param1 = dp->d_atom.param1 + index;
		else
			atom[index].param1 = dp->d_atom.param1;
		atom[index].param2 = dp->d_atom.param2;
		atom[index].parent = index + dp->d_atom.parent;
		dp++;
		pptr++;
		index++;
	} while	(dp->d_atom.type != NOTANATOM);
	cfree(molpos);
	return(index - natoms - 1);
}


/*****************************************************************************/

aword	*
mol_match(acp,wp)
	char *		acp;
	register aword	*wp;
{
	register char *		cp;
	register char *		dp;

	if (acp == NULL)
		return(NULL);
	for ( ; (dp = wp->w_mandatory) != 0; ++wp) {
		cp = acp;
		while (*cp++ == *dp++)
			if (*dp == '\0') {	/* mandatory part matched */
				dp = wp->w_optional;
				do	if	(*cp == '\0')
						return(wp);
					while	(*cp++ == *dp++);
				break;	/* if optional part not matched */
			}
	}
	return(NULL);
}



/*****************************************************************************/

bc_type	*
bc_match(acp,wp)
	char 			*acp;
	register bc_word	*wp;
{
	register char *		cp;
	register char *		dp;

	if (acp == NULL)
		return(NULL);
	for ( ; (dp = wp->bw_mand) != 0; ++wp) {
		cp = acp;
		while (*cp++ == *dp++)
			if (*dp == '\0') {	/* mandatory part matched */
				dp = wp->bw_opt;
				do	if	(*cp == '\0')
						return(wp->bw_ptr);
					while	(*cp++ == *dp++);
				break;	/* if optional part not matched */
			}
	}
	return(NULL);
}



/*  findpt:  find an appropriate place to put the next molecule.  We insist
 *	     that the distance to the nearest neighbor be greater than the
 *	     intermolecular spacing.
 */

/*
 *	IF I WERE TO USE AN imspace FOR EACH OF x,y,z  THEN I WOULD USE
 *	spacing AS A TRIPFLOAT RATHER THAN A FLOAT. THEN IF !POR, 
 *	WE'D HAVE THE SAME VALUE IN EACH OF x,y,z.
 */

find_pt(trial,spacing,naddsofar,bc_type)
	tripfloat	*trial;
	tripfloat	*spacing;	/* one each for x,y,z */
	int		naddsofar;
	int		bc_type;
{
	int		ntrys;
	tripd		*tp;
	tripfloat	edge_const;	/* spacing factor for pbc's so atoms
					don't get too close to a boundary */
	float		edge;

	ntrys = 0;
/*
 *	ASSUMING THE SAME FOR POR AS FOR PC IN TERMS OF edge_const
 *	EXCEPT THAT WE ARE USING DIFFERENT ONES FOR EACH WALL.
 */
	edge = bc_type == SWC ? 0.0 : 0.5;
	edge_const.fx = (bc_type == PTO) ?
			(2. * xwall - spacing->fx * edge) :
			(xwall - spacing->fx * edge);
	edge_const.fy = (bc_type == PTO) ?
			(2. * ywall - spacing->fy * edge) :
			(ywall - spacing->fy * edge);
	edge_const.fz = (bc_type == PTO) ?
			(2. * zwall - spacing->fz * edge) :
			(zwall - spacing->fz * edge);
again:
	if (++ntrys > MAX_FACTOR * 10)
		return(-1);

/*
 *  Generate random numbers in the range [0,1].  This routine will probably
 *  need to be changed if you move to a different computer.
 */
	trial->fx = random()/RMAX;
	trial->fy = random()/RMAX;
	trial->fz = random()/RMAX;

	if (bc_type == PTO) {
		trial->fx -= 0.5;	/* in range [-.5,.5] */
		trial->fy -= 0.5;
		trial->fz -= 0.5;
		if ((fabs(trial->fx) + fabs(trial->fy) + fabs(trial->fz))
			>= 0.75) {
			trial->fx -= 0.5 * sign(trial->fx);
			trial->fy -= 0.5 * sign(trial->fy);
			trial->fz -= 0.5 * sign(trial->fz);
		}
	}
	trial->fx *= edge_const.fx;
	trial->fy *= edge_const.fy;
	trial->fz *= edge_const.fz;
	for (tp = pos; tp < pos + natoms + naddsofar; tp++) {
		if (fabs(tp->fx - trial->fx) < spacing->fx
			&& fabs(tp->fy - trial->fy) <  spacing->fy
			&& fabs(tp->fz - trial->fz) < spacing->fz)
			goto again;
	}
	return(ntrys);
}
