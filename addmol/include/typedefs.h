/*
 *	Type definitions needed for the addmol/readfill/writefill programs.
 *	#include this in .c files needing access to these types.
 */


typedef struct {
	int	type;
	int	flags;
	int	param1;
	int	param2;
	int	parent;
}	parts;

typedef struct {
	char	*w_mandatory;
	char	*w_optional;
	int	w_token;
}	word;

typedef struct {
	short	x;
	short	y;
	short	z;
}	tripint;

typedef struct {
	float	fx;
	float	fy;
	float	fz;
}	tripfloat;

typedef struct	{
	double	fx;
	double	fy;
	double	fz;
}	tripdouble;

typedef tripdouble tripd;

typedef struct {
	parts	d_atom;
	tripint	d_triplet;
}	description;

typedef	struct {
	char		*w_mandatory;
	char		*w_optional;
	description	*w_pnt;
}	aword;
