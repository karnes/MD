/*
 *	extern declarations for the global variables in globals.c.
 *	#include this in each file needing access to these globals.
 */
#include	<sys/types.h>

/*
 *	The basic input-file variables . . .
 */

extern char	filetype[];
extern time_t	datestamp;
extern char	status[];
extern int	natoms;
extern int	nsolute;
extern double	xwall;
extern double	ywall;
extern double	zwall;
extern double	EqTemp;
extern double	DEqTemp;
extern double	EqPress;
extern double	DEqPress;
extern parts	*atom;
extern int	xtrInQ;
extern tripd	*pos;

/*
 *	The molecule descriptions ...
 */

extern description	xraytest[];
extern description	cl2[];
extern description	i2[];
extern description	chloride[];
extern description	lithium[];
extern description	potassium[];
extern description	helium[];
extern description	fluoride[];
extern description	sodium[];
extern description	cesium[];
extern description	crystal[];
extern description	argon[];
extern description	xenon[];
extern description	krypton[];
extern description	hcl[];
extern description	co[];
extern description	o2[];
extern description	n2[];
extern description	watermol[];
extern description	dmso[];
extern description	methane[];
extern description	ch3d[];
extern description	ch2d2[];
extern description	chd3[];
extern description	cd4[];
extern description	ch2i2[];
extern description	ccl4[];
extern description	chcl3[];
extern description	ccl4_zrp[];
extern description	n2_zrp[];
extern description	cl2_zrp[];
extern description	ethane[];
extern description	hexane[];
extern description	hexadecane[];
extern description	ethanol[];
extern description	glycol[];
extern description	chexane[];
extern description	clch3cl[];
extern description	dabco[];
extern description	dclethane[];
extern description	acetonitrile[];
extern aword		molecules[];
