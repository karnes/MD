/*
 *	water variables
 */
//#define	TBLSZ	9000
char waterP[8];		/* Type of water potential			*/
double	VWATTS;		/* watts non-bond water energy			*/
double	WATERV;		/* intramolecular water energy			*/
double	qabOH[6];	/* Lennard-Jones Coulomb parameters for O and H	*/
//double	ootbl[TBLSZ];	/* polynomial table for watts O-O interaction	*/
//double	ohtbl[TBLSZ];	/* polynomial table for watts O-H interaction	*/
//double	hhtbl[TBLSZ];	/* polynomial table for watts H-H interaction	*/
//unsigned short	oonext;	/* actual table usage for O-O interaction	*/
//unsigned short	ohnext;	/* actual table usage for O-H interaction	*/
//unsigned short	hhnext;	/* actual table usage for H-H interaction	*/
double waterRdf[500];
//int hbmat[1000][1000];
//int rmat[1000][1000];
