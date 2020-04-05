/*
 *	This file contains the global variable definitions and
 *	initializations used in the addmol, readfill, and writefill
 *	programs.  The corresponding header file "globals.h" gives
 *	extern declarations to be included in whatever .c files want
 *	to access the globals.  If you want to add a new molecule type,
 *	add the name to the list here, add an extern declaration in 
 *	"globals.h", and put an atom type code and mass code in the
 *	file "atomtypes.h".
 */


#include	"atomtypes.h"
#include	"typedefs.h"
#include	<globals.h>
#include	<sys/types.h>

/*
 *	The input-file variables . . .
 */
	char	filetype[] = ".inp";	/* Label the file with .inp */
	time_t	datestamp;	/* Used to stamp the date on the file */
	char	status[4];	/* Label to indicate file status */
	int	natoms;		/* Number of atoms (not mols!) in system */
	int	nsolute;	/* Number of solute atoms in system */
	double	xwall;		/* box size in x-dimension */
	double	ywall;		/* box size in y-dimension */
	double	zwall;		/* box size in z-dimension */
	double	EqTemp;		/* equil. temp. (zero here) */
	double	DEqTemp;	/* sd of eq. temp. (zero here) */
	double	EqPress;	/* equil. press. (zero here) */
	double	DEqPress;	/* sd of eq. press. (zero here) */
	parts	*atom;		/* pointer to atom info array */
	int	xtrInQ;		/* boolean: non-zero if there is more input */
	tripd	*pos;		/* pointer to atom position array */

/*
 *	The various molecule types available to addmol . . .
 */
#define	HNACL	1444

description	xraytest[] = 	{
O_HOL,	A_MAJOR|A_DRAW|A_FIXED|A_MINOR,	3,	4,	0,	-1024,	1024,	-2048,
C_PRIM,		A_DRAW | A_FIXED,	0,	0,	-1,	1024,	2048,	-1024,
C_PRIM,		A_DRAW | A_FIXED,	0,	0,	-1,	0,	0,	1024,
H_ALK,		A_DRAW | A_FIXED,	0,	0,	-2,	2048,	-2048,	1024,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	cl2[] = 	{
CHLORINE,	A_MAJOR|A_DRAW|A_MINOR,	1,	2,	0,	1018,	0,	0,
CHLORINE,	A_DRAW,			0,	0,	-1,	-1018,	0,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	i2[] = 	{
IODINE,		A_MAJOR|A_DRAW|A_MINOR,	1,	2,	0,	1365,	0,	0,
IODINE,		A_DRAW,			0,	0,	-1,	-1365,	0,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	chloride[] = {
CHLORIDE,	A_MAJOR|A_DRAW|A_MINOR,	0,	1,	0,	0,	0,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	lithium[] = {
LITHIUM,	A_MAJOR|A_DRAW|A_MINOR,	0,	1,	0,	0,	0,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	potassium[] = {
POTASSIUM,	A_MAJOR|A_DRAW|A_MINOR,	0,	1,	0,	0,	0,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	helium[] = {
HELIUM,		A_MAJOR|A_DRAW|A_MINOR,	0,	1,	0,	0,	0,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	fluoride[] = {
FLUORIDE,	A_MAJOR|A_DRAW|A_MINOR,	0,	1,	0,	0,	0,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	sodium[] = {
SODIUM,		A_MAJOR|A_DRAW|A_MINOR,	0,	1,	0,	0,	0,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	cesium[] = {
CESIUM,		A_MAJOR|A_DRAW|A_MINOR,	0,	1,	0,	0,	0,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	crystal[] = {
SODIUM,		A_MAJOR|A_DRAW|A_MINOR,	0,	1,	0,	HNACL,	HNACL,	HNACL,
CHLORIDE,	A_MAJOR|A_DRAW|A_MINOR,	0,	1,	0,	-HNACL,	HNACL,	HNACL,
CHLORIDE,	A_MAJOR|A_DRAW|A_MINOR,	0,	1,	0,	HNACL,	-HNACL,	HNACL,
SODIUM,		A_MAJOR|A_DRAW|A_MINOR,	0,	1,	0,	-HNACL,	-HNACL,	HNACL,
CHLORIDE,	A_MAJOR|A_DRAW|A_MINOR,	0,	1,	0,	HNACL,	HNACL,	-HNACL,
SODIUM,		A_MAJOR|A_DRAW|A_MINOR,	0,	1,	0,	-HNACL,	HNACL,	-HNACL,
SODIUM,		A_MAJOR|A_DRAW|A_MINOR,	0,	1,	0,	HNACL,	-HNACL,	-HNACL,
CHLORIDE,	A_MAJOR|A_DRAW|A_MINOR,	0,	1,	0,	-HNACL,	-HNACL,	-HNACL,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	argon[] = {
ARGON,		A_MAJOR|A_DRAW|A_MINOR,	0,	1,	0,	0,	0,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	xenon[] = {
XENON,		A_MAJOR|A_DRAW|A_MINOR,	0,	1,	0,	0,	0,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	krypton[] = {
KRYPTON,	A_MAJOR|A_DRAW|A_MINOR,	0,	1,	0,	0,	0,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	hcl[] = 	{
CHLORIDE,	A_MAJOR|A_DRAW|A_MINOR,	1,	2,	0,	653,	0,	0,
HION,		A_DRAW,			0,	0,	-1,	-653,	0,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	co[] = 	{
COXYGEN,	A_MAJOR|A_DRAW|A_MINOR,	1,	2,	0,	-578,	0,	0,
COCARB,		A_DRAW,			0,	0,	-1,	578,	0,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	o2[] = 	{
OXYGEN,		A_MAJOR|A_DRAW|A_MINOR,	1,	2,	0,	-619,	0,	0,
OXYGEN,		A_DRAW,			0,	0,	-1,	619,	0,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	n2[] = 	{
NITROGEN,	A_MAJOR|A_DRAW|A_MINOR,	1,	2,	0,	-562,	0,	0,
NITROGEN,	A_DRAW,			0,	0,	-1,	562,	0,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	watermol[] = 	{
OXYGEN,		A_MAJOR|A_DRAW|A_MINOR,	2,	3,	0,	0,	67,	0,
HYDROGEN,	A_DRAW,			0,	0,	-1,	-755,	-535,	0,
HYDROGEN,	A_DRAW,			0,	0,	-2,	755,	-535,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	dmso[] = 	{
S_CYS,		A_MAJOR|A_DRAW|A_MINOR,	3,	4,	0,	0,	0,	0,
O_XYL,		A_DRAW,			0,	0,	-1,	0,	0,	-1567,
CH3_ME,		A_DRAW,			0,	0,	-2,	1765,	0,	531,
CH3_ME,		A_DRAW,			0,	0,	-3,	-408,	1717,	531,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	methane[] = 	{
ENDC,		A_MAJOR|A_DRAW|A_MINOR,	4,	5,	0,	0,	0,	0,
ALKH,		A_DRAW,			0,	0,	-1,	-372,	-1052,	0,
ALKH,		A_DRAW,			0,	0,	-2,	-372,	527,	-911,
ALKH,		A_DRAW,			0,	0,	-3,	-372,	527,	911,
ALKH,		A_DRAW,			0,	0,	-4,	1116,	0,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	ch3d[] = 	{
ENDC,		A_MAJOR|A_DRAW|A_MINOR,	4,	5,	0,	0,	0,	0,
ALKH,		A_DRAW,			0,	0,	-1,	-372,	-1052,	0,
ALKH,		A_DRAW,			0,	0,	-2,	-372,	527,	-911,
ALKH,		A_DRAW,			0,	0,	-3,	-372,	527,	911,
DEUTERIUM,	A_DRAW,			0,	0,	-4,	1116,	0,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	ch2d2[] = 	{
ENDC,		A_MAJOR|A_DRAW|A_MINOR,	4,	5,	0,	0,	0,	0,
ALKH,		A_DRAW,			0,	0,	-1,	-372,	-1052,	0,
ALKH,		A_DRAW,			0,	0,	-2,	-372,	527,	-911,
DEUTERIUM,	A_DRAW,			0,	0,	-3,	-372,	527,	911,
DEUTERIUM,	A_DRAW,			0,	0,	-4,	1116,	0,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	chd3[] = 	{
ENDC,		A_MAJOR|A_DRAW|A_MINOR,	4,	5,	0,	0,	0,	0,
ALKH,		A_DRAW,			0,	0,	-1,	-372,	-1052,	0,
DEUTERIUM,	A_DRAW,			0,	0,	-2,	-372,	527,	-911,
DEUTERIUM,	A_DRAW,			0,	0,	-3,	-372,	527,	911,
DEUTERIUM,	A_DRAW,			0,	0,	-4,	1116,	0,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	cd4[] = 	{
ENDC,		A_MAJOR|A_DRAW|A_MINOR,	4,	5,	0,	0,	0,	0,
DEUTERIUM,	A_DRAW,			0,	0,	-1,	-372,	-1052,	0,
DEUTERIUM,	A_DRAW,			0,	0,	-2,	-372,	527,	-911,
DEUTERIUM,	A_DRAW,			0,	0,	-3,	-372,	527,	911,
DEUTERIUM,	A_DRAW,			0,	0,	-4,	1116,	0,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	ch2i2[] =	{
C_CH2I2,	A_MAJOR|A_DRAW|A_MINOR,	4,	5,	0,	0,	0,	0,
H_CH2I2,	A_DRAW,			0,	0,	-1,	-528,	914,	-373,
H_CH2I2,	A_DRAW,			0,	0,	-2,	1055,	0,	-373,
I_CH2I2,	A_DRAW,			0,	0,	-3,	0,	0,	2171,
I_CH2I2,	A_DRAW,			0,	0,	-4,	-1023,	-1773,	-724,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	ccl4[] = 	{
ENDC,		A_MAJOR|A_DRAW|A_MINOR,	4,	5,	0,	0,	0,	0,
ALKCL,		A_DRAW,			0,	0,	-1,	-596,	-1686,	0,
ALKCL,		A_DRAW,			0,	0,	-2,	-596,	845,	-1460,
ALKCL,		A_DRAW,			0,	0,	-3,	-596,	845,	1460,
ALKCL,		A_DRAW,			0,	0,	-4,	1788,	0,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	chcl3[] = 	{
ENDC,		A_MAJOR|A_DRAW|A_MINOR,	4,	5,	0,	0,	0,	0,
ALKCL,		A_DRAW,			0,	0,	-1,	-596,	-1686,	0,
ALKCL,		A_DRAW,			0,	0,	-2,	-596,	845,	-1460,
ALKCL,		A_DRAW,			0,	0,	-3,	-596,	845,	1460,
ALKH,		A_DRAW,			0,	0,	-4,	1096,	0,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

/*
 *	The following types are for the ZRP system
 */
description	ccl4_zrp[] =	{
CCL4,		A_MAJOR|A_DRAW|A_MINOR,	0,	1,	0,	0,	0,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	n2_zrp[] = {
NITROGEN,	A_MAJOR|A_DRAW|A_MINOR,	1,	2,	0,	-662,	0,	0,
NITROGEN,	A_DRAW,			0,	0,	-1,	662,	0,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	cl2_zrp[] = {
CHLORINE,	A_MAJOR|A_DRAW|A_MINOR,	1,	2,	0,	-968,	0,	0,
CHLORINE,	A_DRAW,			0,	0,	-1,	968,	0,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

/*	End ZRP entries. */

/*	SN2 system	*/

description	clch3cl[] = {
CTCL1,		A_MAJOR|A_DRAW|A_MINOR,	2,	3,	0,	0,	0,	0,
CTCH3,		A_DRAW,			0,	0,	-1,	2440,	0,	0,
CTCL2,		A_DRAW,			0,	0,	-2,	4880,	0,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	ethane[] = 	{
C_PRIM,		A_MAJOR|A_DRAW|A_MINOR,	7,	4,	0,	-788,	0,	0,
H_ALK,		A_DRAW,			0,	0,	-1,	-1150,	-1052,	0,
H_ALK,		A_DRAW,			0,	0,	-2,	-1150,	527,	-911,
H_ALK,		A_DRAW,			0,	0,	-3,	-1150,	527,	911,
C_PRIM,		A_DRAW|A_MINOR,		0,	4,	-4,	778,	0,	0,
H_ALK,		A_DRAW,			0,	0,	-1,	1171,	-526,	911,
H_ALK,		A_DRAW,			0,	0,	-2,	1171,	-526,	-911,
H_ALK,		A_DRAW,			0,	0,	-3,	1171,	1052,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	hexane[] = 	{
C_PRIM,		A_MAJOR|A_DRAW|A_MINOR,	19,	4,	0,	-7672,	-5924,	0,
H_ALK,		A_DRAW,			0,	0,	-1,	-7289,	-6450,	911,
H_ALK,		A_DRAW,			0,	0,	-2,	-7289,	-6450,	-911,
H_ALK,		A_DRAW,			0,	0,	-3,	-8795,	-5924,	0,
C_SECOND,	A_DRAW|A_MINOR,		0,	3,	-4,	-7133,	-4443,	0,
H_ALK,		A_DRAW,			0,	0,	-1,	-7495,	-3916,	-911,
H_ALK,		A_DRAW,			0,	0,	-2,	-7495,	-3916,	911,
C_SECOND,	A_DRAW|A_MINOR,		0,	3,	-3,	-5557,	-4443,	0,
H_ALK,		A_DRAW,			0,	0,	-1,	-5174,	-4969,	911,
H_ALK,		A_DRAW,			0,	0,	-2,	-5174,	-4969,	-911,
C_SECOND,	A_DRAW|A_MINOR,		0,	3,	-3,	-5018,	-2962,	0,
H_ALK,		A_DRAW,			0,	0,	-1,	-5380,	-2435,	-911,
H_ALK,		A_DRAW,			0,	0,	-2,	-5380,	-2435,	911,
C_SECOND,	A_DRAW|A_MINOR,		0,	3,	-3,	-3442,	-2962,	0,
H_ALK,		A_DRAW,			0,	0,	-1,	-3059,	-3488,	911,
H_ALK,		A_DRAW,			0,	0,	-2,	-3059,	-3488,	-911,
C_PRIM,		A_DRAW|A_MINOR,		0,	4,	-3,	-2903,	-1481,	0,
H_ALK,		A_DRAW,			0,	0,	-1,	-3265,	-954,	-911,
H_ALK,		A_DRAW,			0,	0,	-2,	-3265,	-954,	911,
H_ALK,		A_DRAW,			0,	0,	-3,	-1327,	-1481,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	hexadecane[] = 	{
C_PRIM,		A_MAJOR|A_DRAW|A_MINOR,	49,	4,	0,	-7672,	-5924,	0,
H_ALK,		A_DRAW,			0,	0,	-1,	-7289,	-6450,	911,
H_ALK,		A_DRAW,			0,	0,	-2,	-7289,	-6450,	-911,
H_ALK,		A_DRAW,			0,	0,	-3,	-8795,	-5924,	0,
C_SECOND,	A_DRAW|A_MINOR,		0,	3,	-4,	-7133,	-4443,	0,
H_ALK,		A_DRAW,			0,	0,	-1,	-7495,	-3916,	-911,
H_ALK,		A_DRAW,			0,	0,	-2,	-7495,	-3916,	911,
C_SECOND,	A_DRAW|A_MINOR,		0,	3,	-3,	-5557,	-4443,	0,
H_ALK,		A_DRAW,			0,	0,	-1,	-5174,	-4969,	911,
H_ALK,		A_DRAW,			0,	0,	-2,	-5174,	-4969,	-911,
C_SECOND,	A_DRAW|A_MINOR,		0,	3,	-3,	-5018,	-2962,	0,
H_ALK,		A_DRAW,			0,	0,	-1,	-5380,	-2435,	-911,
H_ALK,		A_DRAW,			0,	0,	-2,	-5380,	-2435,	911,
C_SECOND,	A_DRAW|A_MINOR,		0,	3,	-3,	-3442,	-2962,	0,
H_ALK,		A_DRAW,			0,	0,	-1,	-3059,	-3488,	911,
H_ALK,		A_DRAW,			0,	0,	-2,	-3059,	-3488,	-911,
C_SECOND,	A_DRAW|A_MINOR,		0,	3,	-3,	-2903,	-1481,	0,
H_ALK,		A_DRAW,			0,	0,	-1,	-3265,	-954,	-911,
H_ALK,		A_DRAW,			0,	0,	-2,	-3265,	-954,	911,
C_SECOND,	A_DRAW|A_MINOR,		0,	3,	-3,	-1327,	-1481,	0,
H_ALK,		A_DRAW,			0,	0,	-1,	-944,	-2007,	911,
H_ALK,		A_DRAW,			0,	0,	-2,	-944,	-2007,	-911,
C_SECOND,	A_DRAW|A_MINOR,		0,	3,	-3,	-788,	0,	0,
H_ALK,		A_DRAW,			0,	0,	-1,	-1150,	527,	-911,
H_ALK,		A_DRAW,			0,	0,	-2,	-1150,	527,	911,
C_SECOND,	A_DRAW|A_MINOR,		0,	3,	-3,	778,	0,	0,
H_ALK,		A_DRAW,			0,	0,	-1,	1171,	-526,	911,
H_ALK,		A_DRAW,			0,	0,	-2,	1171,	-526,	-911,
C_SECOND,	A_DRAW|A_MINOR,		0,	3,	-3,	1327,	1481,	0,
H_ALK,		A_DRAW,			0,	0,	-1,	965,	2008,	-911,
H_ALK,		A_DRAW,			0,	0,	-2,	965,	2008,	911,
C_SECOND,	A_DRAW|A_MINOR,		0,	3,	-3,	2903,	1481,	0,
H_ALK,		A_DRAW,			0,	0,	-1,	3286,	955,	911,
H_ALK,		A_DRAW,			0,	0,	-2,	3286,	955,	-911,
C_SECOND,	A_DRAW|A_MINOR,		0,	3,	-3,	3442,	2962,	0,
H_ALK,		A_DRAW,			0,	0,	-1,	3080,	3489,	-911,
H_ALK,		A_DRAW,			0,	0,	-2,	3080,	3489,	911,
C_SECOND,	A_DRAW|A_MINOR,		0,	3,	-3,	5018,	2962,	0,
H_ALK,		A_DRAW,			0,	0,	-1,	5401,	2436,	911,
H_ALK,		A_DRAW,			0,	0,	-2,	5401,	2436,	-911,
C_SECOND,	A_DRAW|A_MINOR,		0,	3,	-3,	5557,	4443,	0,
H_ALK,		A_DRAW,			0,	0,	-1,	5195,	4970,	-911,
H_ALK,		A_DRAW,			0,	0,	-2,	5195,	4970,	911,
C_SECOND,	A_DRAW|A_MINOR,		0,	3,	-3,	7133,	4443,	0,
H_ALK,		A_DRAW,			0,	0,	-1,	7516,	3917,	911,
H_ALK,		A_DRAW,			0,	0,	-2,	7516,	3917,	-911,
C_PRIM,		A_DRAW|A_MINOR,		0,	4,	-3,	7672,	5924,	0,
H_ALK,		A_DRAW,			0,	0,	-1,	7310,	6451,	-911,
H_ALK,		A_DRAW,			0,	0,	-2,	7310,	6451,	911,
H_ALK,		A_DRAW,			0,	0,	-3,	8795,	5924,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	ethanol[] =	{
C_PRIM,		A_MAJOR|A_DRAW|A_MINOR,	8,	4,	0,	-788,	0,	0,
H_ALK,		A_DRAW,			0,	0,	-1,	-1150,	-1052,	0,
H_ALK ,		A_DRAW,			0,	0,	-2,	-1150,	527,	-911, 
H_ALK,		A_DRAW,			0,	0,	-3,	-1150,	527,	911,
C_SECOND,	A_DRAW|A_MINOR,		0,	3,	-4,	778,	0,	0,
H_ALK,		A_DRAW,			0,	0,	-1,	1171,	-526,	911,
H_ALK,		A_DRAW,			0,	0,	-2,	1171,	-526,	-911,
O_HOL,		A_DRAW|A_MINOR,		0,	2,	-3,	1302,	1487,	0,
H_HOL,		A_DRAW,			0,	0,	-1,	2326,	1487,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	glycol[] =	{
C_SECOND,	A_MAJOR|A_DRAW|A_MINOR,	9,	3,	0,	-788,	0,	0,
H_ALK,		A_DRAW,			0,	0,	-1,	-1150,	527,	-911, 
H_ALK,		A_DRAW,			0,	0,	-2,	-1150,	527,	911,
O_HOL,		A_DRAW|A_MINOR,		0,	2,	-3,	-1302,	-1487,	0,
H_HOL,		A_DRAW,			0,	0,	-1,	-2326,	-1487,	0,
C_SECOND,	A_DRAW|A_MINOR,		0,	3,	-5,	778,	0,	0,
H_ALK,		A_DRAW,			0,	0,	-1,	1171,	-526,	911,
H_ALK,		A_DRAW,			0,	0,	-2,	1171,	-526,	-911,
O_HOL,		A_DRAW|A_MINOR,		0,	2,	-3,	1302,	1487,	0,
H_HOL,		A_DRAW,			0,	0,	-1,	2326,	1487,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	chexane[] = {
H_ALK,		A_DRAW|A_MAJOR|A_MINOR,	17,	3,	0,	-1161,	2199,	-526,
C_SECOND,	A_RING | A_DRAW,	14,	0,	-1,	-788,	1287,	0,
H_ALK,		A_DRAW,			0,	0,	-1,	-1161,	1286,	1052,
C_SECOND,	A_DRAW|A_MINOR,		0,	3,	-2,	788,	1287,	0,
H_ALK,		A_DRAW,			0,	0,	-1,	1161,	1286,	-1052,
H_ALK,		A_DRAW,			0,	0,	-2,	1161,	2199,	526,
C_SECOND,	A_DRAW|A_MINOR,		0,	3,	-3,	1314,	0,	744,
H_ALK,		A_DRAW,			0,	0,	-1,	2430,	0,	745,
H_ALK,		A_DRAW,			0,	0,	-2,	941,	0,	1796,
C_SECOND,	A_DRAW|A_MINOR,		0,	3,	-3,	788,	-1287,	0,
H_ALK,		A_DRAW,			0,	0,	-1,	1161,	-1286,	-1052,
H_ALK,		A_DRAW,			0,	0,	-2,	1161,	-2199,	526,
C_SECOND,	A_DRAW|A_MINOR,		0,	3,	-3,	-788,	-1287,	0,
H_ALK,		A_DRAW,			0,	0,	-1,	-1161,	-2199,	-526,
H_ALK,		A_DRAW,			0,	0,	-2,	-1161,	-1286,	1052,
C_SECOND,	A_RING|A_DRAW|A_MINOR,	-14,	3,	-3,	-1314,	0,	-744,
H_ALK,		A_DRAW,			0,	0,	-1,	-941,	0,	-1796,
H_ALK,		A_DRAW,			0,	0,	-2,	-2430,	0,	-745,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	dabco[] = {
N_1,	A_RING|A_DRAW|A_MAJOR|A_MINOR,	7, 8,	0,	0,	0,	1292,
C_SECOND,	A_DRAW|A_MINOR,	0,	0,	0,	1421,	0,	790,
C_SECOND,	A_DRAW|A_MINOR,	0,	0,	0,	-710,	1230,	790,
C_SECOND,	A_DRAW|A_MINOR,	0,	0,	0,	-710,	-1230,	790,
N_1,	A_RING|A_DRAW|A_MINOR,	0,	0,	0,	0,	0,	-1292,
C_SECOND,	A_DRAW|A_MINOR,	0,	0,	0,	1421,	0,	-790,
C_SECOND,	A_DRAW|A_MINOR,	0,	0,	0,	-710,	1230,	-790,
C_SECOND,	A_DRAW|A_MINOR,	0,	0,	0,	-710,	-1230,	-790,
NOTANATOM,	0,		0,	0,	0,	0,	0,	0};

description	dpolar[] = {
P_PLUS,	A_MAJOR|A_DRAW|A_MINOR,	1,	0,	0,	-1024,	0,	0,
P_MINUS,	A_MAJOR|A_DRAW|A_MINOR,	1,	0,	0,	 1024,	0,	0,
NOTANATOM,	0,		0,	0,	0,	0,	0,	0};

description	dnpolar[] = {
NON_P,	A_MAJOR|A_DRAW|A_MINOR,	1,	0,	0,	-1024,	0,	0,
NON_P,	A_MAJOR|A_DRAW|A_MINOR,	1,	0,	0,	 1024,	0,	0,
NOTANATOM,	0,		0,	0,	0,	0,	0,	0};

description	dclethane[] = 	{
C_CH2CL2,	A_MAJOR|A_DRAW|A_MINOR,	3,	2,	0,    205,  540,   771,
CL_CH2CL2,	A_DRAW,			0,	0,     -1,   1820,  379,   931,
C_CH2CL2,	A_DRAW|A_MINOR,		0,	2,     -2,     58, -634,  -376,
CL_CH2CL2,	A_DRAW,			0,	0,     -1,  -1909, -347, -1065,
NOTANATOM,	0,			0,	0,	0,	0,	0,   0};

description	acetonitrile[] = 	{
CARBON,		A_MAJOR|A_DRAW|A_MINOR,	2,	3,	0,	139,	0,	0,
N_PEP,		A_DRAW,			0,	0,	-1,	1337,	0,	0,
CH3_ME,		A_DRAW,			0,	0,	-2,	-1356,	0,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

description	united_ethanol[] =	{
C_CH2CL2,	A_MAJOR|A_DRAW|A_MINOR,	3,	2,	0,	-788,	0,	0,
CH3_ME,		A_DRAW,			0,	0,	-1,	778,	0,	0,
O_HOL,		A_DRAW|A_MINOR,		0,	2,	-2,	1302,	1487,	0,
H_HOL,		A_DRAW,			0,	0,	-1,	2326,	1487,	0,
NOTANATOM,	0,			0,	0,	0,	0,	0,	0};

aword molecules[] = {
	"xray",		"test",		xraytest,
	"chorin",	"e",		cl2,
	"f",		"luoride",	fluoride,
	"chlorid",	"e",		chloride,
	"i",		"odine",	i2,
	"l",		"ithium",	lithium,
	"s",		"odium",	sodium,
	"p",		"otassium",	potassium,
	"ce",		"sium",		cesium,
	"cr",		"ystal",	crystal,
	"a",		"rgon",		argon,
	"kr",		"ypton",	krypton,
	"xe",		"non",		xenon,
	"H",		"Cl",		hcl,
	"C",		"O",		co,
	"o",		"xygen",	o2,
	"n",		"itrogen",	n2,
	"w",		"ater",		watermol,
	"d",		"mso",		dmso,
	"m",		"ethane",	methane,
	"ch3",		"d",		ch3d,
	"ch2",		"d2",		ch2d2,
	"chd",		"3",		chd3,
	"chcl",		"3",		chcl3,
	"cd",		"4",		cd4,
	"ch2i",		"2",		ch2i2,
	"cy",		"clohexane",	chexane,
	"he",		"lium",		helium,
	"e",		"thane",	ethane,
	"ethano",	"l",		ethanol,
	"gly",		"col",		glycol,
	"hex",		"adecane",	hexadecane,
	"acetonit",	"rile",		acetonitrile,
	"ca",		"rbontet",	ccl4,
	"hexan",	"e",		hexane,
	"ccl4_u",	"nit",		ccl4_zrp,
	"n2_z",		"rp",		n2_zrp,
	"cl2_z",	"rp",		cl2_zrp,
	"clch3",	"cl",		clch3cl,
	"dab",		"co",		dabco,
	"dpo",		"lar",		dpolar,
	"dnp",		"olar",		dnpolar,
	"dcleth",	"ane",		dclethane,
	"united",	"_ethanol",	united_ethanol,
	0,		0,		0};
