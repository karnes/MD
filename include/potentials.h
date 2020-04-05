/*
 *	The routine that loads this file must #include md.h before this.
 */
/*
 *	From:  A. Warshel and S. Lifson, J. Chem. Phys., 53, 582 (1970).
 *
 * Nonbonded	0.5*r	 epsilon**.5	e
 *		     eq			 eff
 *
 *   H---H	1.774	   0.0508	0.11
 *   C---C	1.808	   0.4297
 *
 *
 */
#define	QH	.11
#define	QC	(-2.*QH)
    /*	q		a		b 
	AMUAFS*A	AMUAFS*A**9	AMUAFS*A**6 */

ljcon	lj_CBUT[3] = {			/* LJ6-9 */
	QH*QH/E2,	459.84336/KCAL,	15.443664/KCAL,	9,	/* H-H */
	QH*QC/E2,	4238.2765/KCAL,	138.32599/KCAL,	9,	/* H-C */
	QC*QC/E2,	39031.613/KCAL,	1238.2903/KCAL,	9};	/* C-C */

    /*	q		a		b 
	AMUAFS*A	AMUAFS*A**12	AMUAFS*A**6 */

/*
 *	Lennard-Jones Parameters used for the a+bc project
 */
ljcon	Cl_He	= {0.,	5.3764e1,	1.0447e-1,	12};	/* He-Cl */
ljcon   H_He    = {0.,	5.317139e2,	3.124527e0,	12};	/* H-He */
ljcon   He_He   = {0.,	2.981100e0,	1.034703e-2,	12};	/* He-He */
ljcon	Cl_Ar	= {0.,	9.6065e2,	8.0626e-1,	12};	/* Ar-Cl */
ljcon   H_Ar    = {0.,	5.805436e1,	1.287752e-1,	12};	/* H-Ar */
ljcon   Ar_Ar   = {0.,	9.521637e2,	6.163651e-1,	12}; 	/* Ar-Ar */
ljcon	Cl_Xe   = {0.,	2.9137e3,	1.6866e0,	12};	/* Xe-Cl */
ljcon   H_Xe    = {0.,	2.195699e5,	3.054028e2,	12};	/* H-Xe */
ljcon   Xe_Xe   = {0.,	8.757112e3,	2.697424e0,	12};	/* Xe-Xe */

ljcon	C_He    = {0.,	4.226988e3,	7.276628e0,	12};	/* C-He */
ljcon	C_Ar    = {0.,	4.801197e2,	9.658471e-2,	12};	/* C-Ar */
ljcon	C_Xe    = {0.,	1.076100e6,	5.584476e2,	12};	/* C-Xe */

/*
 *	Lennard-Jones Parameters for ZRP interactions
 *	Format:	      Q          A        B     n
 *	Units:  [length] = angstrom, [mass] = amu, [time] = fs
 */
ljcon	Ar_N_z      = {0.0,	432.5,	0.3880,	12};	/* Ar-N        */
ljcon	Ar_Cl_z     = {0.0,	1098.,	0.8688,	12};	/* Ar-Cl       */
ljcon	Ar_Ar_z     = {0.0,  	967.4,  0.6207, 12};	/* Ar-Ar       */
ljcon	CCl4_N_z    = {0.0,     3.544e4,   44.98,	12};	/* CCl_4 - N  */
ljcon	CCl4_CCl4_z = {0.0,     1.861e6,   4.515, 12};	/* CCl_4-CCl_4 */


/*
 *	Lennard-Jones (6-9 for N2, 6-12 for rest) Parameters q, a, b

ljcon	N2_N2 = {0.0,	1.42680e-1,	1.61818e-1,	9};
ljcon	N2_Ar = {0.0,	3.006328e1,	2.32322e-2,	12};

 */

/*
 *	Lennard-Jones Parameters for ICN - rare gas interactions
 *      Taken from "Computer Simulation of Liquids, by M.P. Allen
 *      and D.J. Tildesley, Clarendon Press 1987, Chap. 1".
 */
ljcon He_He_i=  {0.,   1.599e+03/KCAL,  1.139e+01/KCAL, 12};   
ljcon I_He_i  =  {0.,   4.337e+05/KCAL,  3.680e+02/KCAL, 12};
ljcon C_He_i  =  {0.,   4.496e+04/KCAL,  9.036e+01/KCAL, 12};   
ljcon N_He_i  =  {0.,   3.523e+04/KCAL,  7.389e+01/KCAL, 12};   
ljcon Ne_Ne_i =  {0.,   6.125e+04/KCAL,  1.512e+02/KCAL, 12};   
ljcon I_Ne_i  =  {0.,   1.114e+06/KCAL,  9.114e+02/KCAL, 12};   
ljcon C_Ne_i  =  {0.,   2.381e+05/KCAL,  3.046e+02/KCAL, 12};   
ljcon N_Ne_i  =  {0.,   1.877e+05/KCAL,  2.499e+02/KCAL, 12};   
ljcon Ar_Ar_i =  {0.,   2.353e+06/KCAL,  1.497e+03/KCAL, 12};   
ljcon I_Ar_i  =  {0.,   6.599e+06/KCAL,  2.982e+03/KCAL, 12};         
ljcon C_Ar_i  =  {0.,   1.384e+06/KCAL,  9.279e+02/KCAL, 12};   
ljcon N_Ar_i  =  {0.,   1.100e+06/KCAL,  7.643e+02/KCAL, 12};   
ljcon Kr_Kr_i =  {0.,   1.298e+07/KCAL,  4.113e+03/KCAL, 12};   
ljcon I_Kr_i  =  {0.,   1.522e+07/KCAL,  4.897e+03/KCAL, 12};     
ljcon C_Kr_i  =  {0.,   3.337e+06/KCAL,  1.559e+03/KCAL, 12};   
ljcon N_Kr_i  =  {0.,   2.663e+06/KCAL,  1.287e+03/KCAL, 12};   
ljcon Xe_Xe_i =  {0.,   2.093e+07/KCAL,  6.447e+03/KCAL, 12};   
ljcon I_Xe_i  =  {0.,   1.932e+07/KCAL,  6.132e+03/KCAL, 12};   
ljcon C_Xe_i  =  {0.,   4.245e+06/KCAL,  1.953e+03/KCAL, 12};   
ljcon N_Xe_i  =  {0.,   3.389e+06/KCAL,  1.612e+03/KCAL, 12};   

/*
 * LJ parameters for the interactions between Helium and the carbon and
 * nitrogen of dabco calculated from the data in Allen&Tildesley page 21.
 */

ljcon C_He_d  =  {0.,   4.496e+04/KCAL,  9.036e+01/KCAL, 12};   
ljcon N_He_d  =  {0.,   3.523e+04/KCAL,  7.389e+01/KCAL, 12};   

/*
 * LJ parameters for the diatomic liquid. taken from Ciccotti et. al.
 * the 3D debye dipole moment case
 */

ljcon lj_PP_PP = {0.097534/E2, 5.372e+06/KCAL,	 2.922e+03/KCAL, 12};
ljcon lj_PP_PM = {-0.097534/E2, 5.372e+06/KCAL, 2.922e+03/KCAL, 12};
ljcon lj_PP_NP = {0.,	 5.372e+06/KCAL, 2.922e+03/KCAL, 12};
ljcon lj_PM_PM = {0.097534/E2, 5.372e+06/KCAL,	 2.922e+03/KCAL, 12};
ljcon lj_PM_NP = {0.,	 5.372e+06/KCAL, 2.922e+03/KCAL, 12};
ljcon lj_NP_NP = {0.,	 5.372e+06/KCAL, 2.922e+03/KCAL, 12};

/* same as above for the 2Debye case */

ljcon lj_PP_PP2 = {0.043348/E2, 5.372e+06/KCAL,	 2.922e+03/KCAL, 12};
ljcon lj_PP_PM2 = {-0.043348/E2, 5.372e+06/KCAL, 2.922e+03/KCAL, 12};
ljcon lj_PM_PM2 = {0.043348/E2, 5.372e+06/KCAL,	 2.922e+03/KCAL, 12};

/* LJ-Coulomb interaction for water chlorine ion			*/
/*first choice
ljcon O_Cl = {0.6666667/E2,	305.5,	0.2954,	12};
ljcon H_Cl = {-0.3333333/E2,	87.40,	0.0966,	12};
**
**second choice
ljcon O_Cl = {0.6666667/E2,	2413.7,	1.206,	12};
ljcon H_Cl = {-0.3333333/E2,	5.558,	0.05787, 12};
**
second choice with spc water (Pettitt & Rossky, JCP 84 5836 (1986)*/
ljcon O_Cl = {0.82/E2,	2413.7,	1.206,	12};
ljcon H_Cl = {-0.41/E2,	5.558,	0.05787, 12};

/* LJ-Coulomb interaction for water sodium ion		****
 * Pettitt & Rossky, JCP 84 5836 (1986)				****/
ljcon O_Na = {-0.82/E2,	3.674e+01, 9.072e-02, 12};
ljcon H_Na = {0.41/E2,	5.722e-03, 1.132e-03, 12};
