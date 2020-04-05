/*
 * N.B.:
 *	THE UNITS INTERNAL TO THIS PROGRAM ARE THE FOLLOWING:
 *		mass:	gm / mole	(amu)
 *		length:	angstrom	(A)
 *		time:	femtosecond	(fs)
 *
 *	Charges  (qi*qj) should be defined in terms of electrons,
 *	then DIVIDED BY "E2", e.g., QH*QC/E2.
 *
 *	This implies that energies are in the units: amu A**2 / fs**2.
 *
 *	The following definitions are for conversions FROM internal
 *	units (AMUAFS) TO those defined.
 */
#define	AMUAFS		1.
#define	ERG		1.e14		/* --> erg/mol */
#define	JOULE		1.e7		/* --> J/mol */
#define	KJOULE		1.e4		/* --> kJ/mol */
#define	CAL		(1.e7/4.184)	/* --> cal/mol */
#define	KCAL		(1.e4/4.184)	/* --> kcal/mol */
#define	EVOLT		(JOULE/9.65177e4)	/* --> eV/molecule */
#define	CM1		(EVOLT*8.06573e3)	/* wavenumbers */
#define	SEC		1.e-15		/* --> sec*/
#define HARTREE		(EVOLT/27.211396) /* = 3.8075202*/

/*
 *	The following constants are not energy conversion factors,
 *	but are, none the less, very useful.
 */
#define	PI		3.14159265358979323846
#define	C		2.99792458e-5	/* [cm/fs]*/
#define	AVOGADRO	6.022045e23
#define	ELECTRON	1.6021892e-19	/* charge of proton [Coul] */
#define	EPS0_AMU	1.11265e-13	/* for E2 [Coul^2 / AMUAFS-A] */
#define	E2		(EPS0_AMU/(ELECTRON*ELECTRON*AVOGADRO)) /*7.197587876*/
#define	KBOLTZ		8.31441e-7	/* gas constant [AMUAFS/degK] */
#define	TWOPIC		(2.*PI*C)	/* 2 PI * (speed of light) [ cm/fs] */
#define	HBAR		(1.0545887e-34*AVOGADRO/JOULE/SEC)
					/*Planck's constant [AMUAFS * fs]*/
/*	To get the pressure in atm. multiply the one calculated in the 
 *	simulation units (amu/(A*fs*fs)) by ATMOS
 */
#define ATMOS		(82.056/(KBOLTZ*AVOGADRO*1.0e-24)) /* --> atm */
/*reduced pressure = atm*sigma**3/epsilon/REDPRES*/
#define REDPRES		(KBOLTZ*ATMOS)	/* approx = 136.26*/
#define BOHR	0.52917725
