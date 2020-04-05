/*
 *	Atom description structure.
 */
typedef struct	{
 	int	type;		/* the atom type number */
 	int	flags;		/* flags bits for each atom (see below) */
 	int	param1;		/* misc. atom parameters; currently the
				 * low mantissa is used for either the # of
				 * atoms in this molecule or the atom ring
				 * parameter.  the other fields are not used.
				 */
 	int	param2;		/* holds the number of atoms in the group if
				 * this atom is a minor atom.
				 */
 	int	parent;		/* this atom's parent number */
}	parts;		/* the atom parts (above) */

typedef struct	{
	short	x;
	short	y;
	short	z;
}	tripint;	/* a short integer triplet */

typedef struct	{
	float	fx;
	float	fy;
	float	fz;
}	tripfloat;	/* a float triplet */

typedef	tripfloat	tripf;

typedef struct	{
	double	fx;
	double	fy;
	double	fz;
	}	tripdouble;	/* a double triplet */

typedef	tripdouble	tripd;
// jjk - added here
typedef struct	{
	float	fx;
	float	fy;
	float	fz;
	float  zpos;
	int    inp;       /* is x in b-CD pore? */
	}	tripfzp;	/* a tripf with z-position holder and in-pore boolean */

typedef struct	{
	float	fx;
	float	fy;
	float	fz;
	float  zpos;
	}	tripfz;	/* a tripf with z-position holder */

typedef struct 	{
	double	q;		/* effective squared charge */
	double	a;		/* term for i/r^n */
	double	b;		/* term for i/r^6 */
	int	n;		/* power of first term */
} ljcon;

/*
 *	constants for two and three body interactions
 */
typedef struct {
	double	f3_theta;	/* theta 0 */
	double	f3_ktheta;	/* .5 * Ktheta */
	double	f3_ltheta;	/* L theta */
	double	f3_D1;		/* D */
	double	f3_beq1;	/* equilibrium bond distance */
	double	f3_alpha1;	/* alpha */
	double	f3_D2;		/* D */
	double	f3_beq2;	/* equilibrium bond distance */
	double	f3_alpha2;	/* alpha */
	double	f3_Krarb;	/* bond-bond cross force constant */
	double	f3_Kratheta;	/* bond-angle cross force constant */
	double	f3_Krbtheta;	/* bond-angle cross force constant */
	}	con3list;

/*
 *	structure for listing four body interactions.
 */
typedef struct	{
	int	b4_1;		/* "lesser" outside atom type */
	int	b4_2;		/* "lesser" inside atom type */
	int	b4_3;		/* "greater" inside atom type */
	int	b4_4;		/* "greater" outside atom type */
	int	b4_seq;		/* index to force constant (kphi) */
	}	bod4list;

/*
 *      structure for the vector-vector correlation functions
 */
typedef struct {
	double f1;
	double f2;
	double f3;
	double f4;
	double f5;
	double f6;
}          sixd;

/*
 *	structure for printing out data files
 */
typedef	struct {
		int dat; /* print .dat file   */ 
		int pos; /* print .pos file   */
		int vel; /* print .vel file   */
		int frc; /* print .frc file   */
		int xtr; /* print .xtr file   */
	} outflags;
/*
 * 	structure for constraints information; constrain the distance
 *	between atom A and atom B to be dis.
 */

typedef struct	{
	int	atomA;
	int	atomB;
	double	dis;
	}	constr;
