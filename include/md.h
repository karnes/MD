#include	<stdio.h>
#include	<typedefs.h>
#include	<sys/types.h>	/* to declare time_t for datestamp */
#include	<units.h>
#include	<atom.h>

/*
 *	Global arrays for positions, velocities, etc..
 */
tripdouble	*pos;			/* pointer to positions		*/
tripdouble	*vel;			/* pointer to velocities	*/
tripdouble	*force;			/* pointer to forces		*/
double		*mass;			/* pointer to masses		*/
double		*imass;			/* pointer to inverse masses	*/
parts		*atom;			/* pointer to the atom parts */
tripdouble	*tvel;			/* pointer to velocities(verlet)*/
int		*itable;		/* index table for partial dump(pd)*/
parts		*satom;			/* pointer to the atom parts for pd */
tripdouble	*spos;			/* pointer to positions for pd	*/
tripdouble	*pTensor;		/* the pressure tensor 		*/

/*
 *	These macros are for general purpose comparisons, etc.
 */
#define	abs(x)		((x)>0.?(x):-(x))
#define max(x,y)	((x)>(y)?(x):(y))
#define min(x,y)	((x)<(y)?(x):(y))
#define	isgn(x)		((x)>0?1:-1)
#define	isign(x,y)	(abs(x)*isgn(y))
#define	sgn(x)		((x)>0.?1.:-1.)
#define	sign(x,y)	(abs(x)*sgn(y))
#define	hsgn(x)		((x)<0.?-.5:.5)
#define sq(x)           ((x)*(x))

/*
 *	This macro is for printing errors (x must be of the form:
 *	x == (stderr, "...\n", arg...) so that the parentheses are included.),
 *	then returning or exiting (y == return or y == exit).
 */
#define	ERROR(x,y)	{fprintf x ; y (1);}

/*
 *	Values in input files.
 */

char	action[5];	/* type of MD execution. e.g .RUN. */
char	filetype[5];	/* filetype of file (.inp, .pos, .dat, ... ) */
char	pbcType[8];	/* Type of periodic boundary conditions		*
			 **** O - Truncated Octahedron 		     ****
			 **** C - Cubic				     ****	
			 */
time_t	datestamp;	/* datestamp for stamping output files */
char	status[4];	/* indicates status of .inp file (ADD, MIN, EQU) */
int	natoms;		/* number of atoms */
int	nsolute;	/* number of solute atoms */
int	nshel;		/* number of solvent atoms + solute in solvation shell*/
double	EqTemp;		/* if status == EQU, indicates the equil. temp. */
double	DEqTemp;	/* if status == EQU, indicates sd of the eq. temp. */
double	EqPress;	/* indicates the equil. press. */
double	DEqPress;	/* indicates the sd of the eq. press. */
double	xwall;		/* x wall of box in angstroms */
double	ywall;		/* y wall of box in angstroms */
double	zwall;		/* z wall of box in angstroms */
int	xtrInQ;		/* extra input file flag */


/*
 *	Other global variables
 */

int	g;		/* number of degrees of freedom */
double	swr2min;	/* swapping region parameters */
double	swr2max;
double	h;		/* time step */
tripd	period;		/* period of boundary conditions */
tripd	iperiod;	/* inverse period */
double	etime;		/* elapsed time */
double	temp;		/* temperature (instantaneous) */
double	press;		/* pressure (instantaneous) */
double  pslabSize;	/* slab size for calculating pressure tensor	*/
int  nslabs;		/* number of slabs for pressure tensor	*/
double	h1;		/* the previous time step */
double	h2;		/* the latest time step */
double 	Teq;		/* equilibrium temperature */
double	s;		/* scaling variable */
double	vs;		/* scaling variable's equivalent to "vel" */
double	Q;		/* scaling variable's "mass" */
double	E; 		/* total energy */
double PREV_E;		/* for debuging unexpected energy jumps*/
double	H;		/* Hamiltonian */
double  rdf[350];	/* radial distribution function			*/
double	swcoef[4];	/* switching parameters for group cut-offs */
double	pswcoef[3];	/* derivates of above */

tripd	P;		/* total momentum				*/
double	V;		/* total potential energy			*/
double	VLIQ;		/* liquid potential energy			*/
double	VSYS;		/* system potential energy			*/
double	VINT;		/* interaction potential energy			*/
double  KLIQ;           /* liquid kinetic energy                        */
double	K;		/* total kinetic energy				*/

ljcon	lj[3][3];	/* LJ constants array				*/

/***	globals for the general constraints package		***/

char cons_method[5];	/* method of constraints UNC, ANA, ITR		*/
char constrain_file[256];/* constraints file name                       */
char cons_info[2];	/* one character flag for type of constraints	*
			 * a - diatomic liquid, l = list of bond, s -special*/
constr	*cons_list;	/* the list of constraints-atom numbers and distances*/
int cons_num;		/* number of constraints			*/
double cons_bond;	/* the bond distance for a diatomic liquid	*/
tripd *prevpos;
tripd *dsigma;

int numDP_x_dataRate;
int numDPx;
int dataRatex;
int tc; /*time counter for correlation functions. see also data/wdat.c*/
