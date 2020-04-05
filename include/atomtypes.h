/*
 *	Atom Type Codes
 */
#define NTYPES		115	/* number of atom types; cannot exceed 128 */

#define	NOTANATOM	(-1)
#define	OXYGEN		0	/* water oxygen */
#define	HYDROGEN	1	/* water hydrogen */
#define HION		2	/* hydrogen ion & h - cl */
#define ALKH		3	/* alkane hydrogen */

/*
 *	The following types are borrowed from Arnie Hagler's protein simulation
 *	VFF (Valence Force Field) programs.  The two letter code given in the
 *	comments is the VFF equivalent designation for that atom type.
 *	The comments describe in general how a given atom type is used
 *	and is not intended as a comprehensive description of the atom type.
 *	The precise meaning of any atom type must be determined by cross
 *	referencing it in the equivalence array for the criteria of interest
 *	and consulting the non-bonded and/or bonded parameter selection
 *	routines in tables.c and bprep.c respectively.
 *	Also very helpful is the book containing the drawings of all the
 *	amino acid and other molecular structures with the atom types labeled.
 */
#define H_AMIDE 	4	/* "hn" - H bonded to N */
#define H_ALK		5	/* "h " - H bonded to alkyl C */
#define H_AROM		6	/* "hp" - H bonded to aromatic C */
#define H_HOL		7	/* "ho" - H bonded to O (alcohol) */
#define H_SULF		8	/* "hs" - H bonded to S */
/*
 *	The next 9 Carbons include their "Hydrogens" implicitly
 */
#define CH3_ALK		9	/* "c " - alkyl C (usually a methyl) */
#define CH3_ME		10	/* "me" - methyl (CH3) */
#define CH_ALK		11	/* "cm" - alkyl C (1 or 2 H's) */
#define CH_AROM		12	/* "cp" - aromatic CH */
#define CH_PRO		13	/* "cn" - side chain CH2 bonded to N in pro */
#define CH_PYRROLE	14	/* "ch" - pyrrole CH in histidine */
#define CH_5INDOLE	15	/* "c*" - CH in the 5 member indole ring (trp)*/
#define	CH_THIOL	16	/* "cc" - thiol or disulfide CH2 (in cys) */
#define	CH_THIOETH	17	/* "cs" - thioether CH2 (in met) */

#define CARBON		18	/* Non-Protein Carbon (hydrocarbon) */
#define ENDC		19	/* Non-Protein Carbon (hydrocarbon) */

#define	C_PRIM		20	/* "c3" - primary C (bonded to 3 H's) */
#define	C_SECOND	21	/* "c2" - secondary C (bonded to 2 H's) */
#define	C_TERT		22	/* "c1" - tertiary C (bonded to 1 H's) */
#define	C_ALPHA		23	/* "ca" - alpha C */
#define	C_GLY		24	/* "cg" - alpha C of glycine */
#define	C_MTHIOETH	25	/* "cu" - thioether methyl C (-CH3) (met) */
#define	C_ETHIOETH	26	/* "cx" - thioether ethyl C (-CH2-) (met) */
#define	C_THIOL		27	/* "cv" - thiol or disulfide C in cys */
#define	C_PYRROLE	28	/* "cf" - pyrrole C in histidine */
#define	C_PRO		29	/* "cj" - side chain C bonded to N in proline */
#define	C_RIBOSE	30	/* "ce" - ribose C for nucleosides */
#define	C_AROM		31	/* "cy" - aromatic C (6 membered rings) */
#define	C_PYRIM		32	/* "c+" - C in pryimidines */
#define	C_6PURINE	33	/* "cq" - C in the 6 membered ring of purines */
#define	C_5PURINE	34	/* "ck" - C in the 5 membered ring of purines */
#define	C_5INDOLE	35	/* "c5" - C in the 5 membered ring of indoles */
#define	C_NYL		36	/* "c'" - carbonyl C */
#define	C_XYL		37	/* "c-" - carboxyl C */
#define	C_ESTER		38	/* "co" - alkyl C bonded to esterified O */
#define C_METHN		39	/* "ci" - methyl C bonded to N (methotrexate) */
#define	C_N3		40	/* "cr" - charged C (+) bonded to 3 N's (arg) */
#define	C_THIOL2	41	/* "cw" - methyl C bonded to SH or S-S */
#define C_T		42	/* "ct" - unused  */

#define O_XYL		43	/* "o-" - carboxyl O */
#define O_NYL		44	/* "o'" - carbonyl O */
#define O_HOL		45	/* "oh" - O bonded to H (alcohol) */
#define O_HOL2		46	/* "o " - O bonded to H (alcohol) */
#define O_PESTER	47	/* "op" - esterified or neutral phosphate O */
#define O_PO4		48	/* "o"" - charged phosphate O */

#define	N_PEP		49	/* "n " - N in peptide bond (or general N) */
#define	N_1		50	/* "n1" - unused */
#define	N_PYRROLE	51	/* "nh" - N in pyrrole rings */
#define	N_P		52	/* "np" - unused */
#define	N_H2		53	/* "n2" - N bonded to C and 2 H's */
#define	N_LYS		54	/* "nl" - terminal N in lysine */
#define	N_H3		55	/* "n3" - N bonded to C and 3 H's (charged) */
#define	N_PYRIM		56	/* "nm" - N peculiar to pyrimidines */
#define	N_PLUS		57	/* "n+" - unused */
#define N_T		58	/* "nt" - unused */

#define S_CYS		59	/* "s1" - S in cysteine (-SH) */
#define S_CYSCYS	60	/* "s " - S in cystine (-S-S-) */

#define P_PHATE		61	/* "p " - phosphate P in nucleotides */
/*
 *	END of "protein" atom types
 */

#define LITHIUM		62	/* lithium ion */
#define	SODIUM		63	/* sodium ion */
#define	POTASSIUM	64	/* potassium ion */
#define	CESIUM		65	/* cesium ion (+1) */
#define	FLUORIDE	66	/* fluoride ion */
#define	CHLORIDE	67	/* chloride ion */
/*
 *	It is imperative that the preceeding types not be unorganized
 *	since "findtblist()" (in tables.c) does a binary search on the
 *	tablist structure and will be messed up if the order of the
 *	tablist structure does not jive with this order.
 */
#define	FLUORINE	68	/* fluorine */

#define ALKCL		69	/* chlorine (as in CCl4) */
#define CHLORINE	70	/* chlorine (as in cl2) */
#define COXYGEN		71	/* co oxygen */
#define COCARB		72	/* co carbon */
#define HELIUM		73	/* helium */
#define ARGON		74	/* argon */
#define KRYPTON		75	/* krypton: superman's poison */
#define XENON		76	/* xenon */
#define NITROGEN	77	/* diatomic nitrogen */
#define DEUTERIUM	78	/* deuterium (for ch3d, ch2d2, chd3, cd4) */
#define IODINE		79	/* iodine */
#define	C_CH2I2		80	/* carbon for ch2i2 */
#define	H_CH2I2		81	/* hydrogen for ch2i2 */
#define	I_CH2I2		82	/* iodine for ch2i2 */
#define	C_CHCL3		83	/* carbon for chcl3 */
#define	H_CHCL3		84	/* hydrogen for chcl3 */
#define	CL_CHCL3	85	/* chlorine for chcl3 */
#define	C_CH2CL2	86	/* carbon for ch2cl2 */
#define	H_CH2CL2	87	/* hydrogen for ch2cl2 */
#define	CL_CH2CL2	88	/* chlorine for ch2cl2 */
#define	C_CHBR3		89	/* carbon for chbr3 */
#define	H_CHBR3		90	/* hydrogen for chbr3 */
#define	BR_CHBR3	91	/* bromine for chbr3 */
#define	C_CH2BR2	92	/* carbon for ch2br2 */
#define	H_CH2BR2	93	/* hydrogen for ch2br2 */
#define	BR_CH2BR2	94	/* bromine for ch2br2 */
/*
 *	The following 3 atoms are for the ABC Charge Transfer
 *	reaction system intended to model an Sn2 reaction.
 */
#define	CTCL1		95	/* 1st Cl atom for (charge xfer) ABC system */
#define	CTCH3		96	/* CH3  "atom" for (charge xfer) ABC system */
#define	CTCL2		97	/* 2nd Cl atom for (charge xfer) ABC system */
/*
 *	We now have a NEON atom type
 */
#define	NEON		98	/* Neon */

/*
 *	For the zrp system, define unified CCl4 type . . .
 */
#define	CCL4		99

/*	Define atom types for polar and nonpolar diatomic liquid	*/

#define P_PLUS		101	/* the positive end of the diatomic 	*/ 	
#define P_MINUS		102	/* the negative end of the diatomic 	*/ 	
#define NON_P		103	/* atom for nonpolar diatomic liquid	*/
#define P_PLUSL		104	/* positive ion with half the mass	*/ 	
#define PLATINUM	105	/* Pt					*/
#define SE		106	/* solvated electron polymer chain	*/
#define SELEC		109	/* classical electron		*/
#define ethyl		110	/* C2H5		*/
#define pentyl		111	/* C3H7		*/
#define butyl		112	/* C4H9		*/
#define Erbium		113	/* Erbium ion +3		*/

/*
 *	Masses for various atoms
 */
#define	SE_MASS	0.1
#define	CSE_MASS	0.00054859
#define	H_MASS	1.00797
#define	D_MASS	2.0141
#define	HE_MASS	4.00260
#define	LI_MASS	6.941
#define	C_MASS	12.01115
#define	CH3_MASS 15.035
#define	CH2_MASS 14.027
#define	CH_MASS 13.019
#define	N_MASS	14.0067
#define	O_MASS	15.9994
#define	F_MASS	18.99840
#define NE_MASS	20.183
#define	NA_MASS	22.98977
#define	P_MASS	30.97376
#define	S_MASS	32.0
#define	CL_MASS	35.453
#define	AR_MASS	39.948
#define	K_MASS	39.098
#define	KR_MASS	83.80
#define	I_MASS	126.9045
#define	XE_MASS	131.30
#define	CS_MASS 132.9054
#define	BR_MASS	78.9183
#define CCL4_MASS 153.82	/* unified CCl4 for zrp system */
#define DP_MASS		40.00	/* polar and non-polar diatoms	*/
#define Pt_MASS		195.09	/* Pt mass			*/
#define	C2H5_MASS 29.062
#define	C3H7_MASS 43.089
#define	C4H9_MASS 57.116
#define	ER_MASS 167.259
