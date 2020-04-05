/*
 *	getmass:
 *
 *	this routine is called in rinp.c to find the masses of the
 *	various atoms in the system to fill the arrays mass and imass.
 */

#include	<stdio.h>
#include	<atomtypes.h>

float
getmass(atm)
	register int	atm;
{
	switch	(atm) {
		case O_XYL:	case O_NYL:	case O_HOL:	
		case O_HOL2:	case O_PESTER:	case O_PO4:	
		case OXYGEN:	case COXYGEN:
			return(O_MASS);
		case HYDROGEN:	case ALKH:	case HION:
		case H_AMIDE:	case H_ALK:	case H_AROM:
		case H_HOL:	case H_SULF:	case H_CH2I2:
		case H_CHCL3:	case H_CH2CL2:	case H_CHBR3:
		case H_CH2BR2:
#ifndef	DEUTER
			return(H_MASS);
#endif
		case DEUTERIUM:
			return(D_MASS);
		case LITHIUM:
			return(LI_MASS);
		case SODIUM:
#ifdef	ARION
			return(AR_MASS);
#else
			return(NA_MASS);
#endif
		case POTASSIUM:
			return(K_MASS);
		case CHLORIDE:	case CHLORINE:	case ALKCL:
		case CL_CHCL3:	case CL_CH2CL2: case CTCL1: case CTCL2:
			return(CL_MASS);
		case CH3_ALK:
			return(C_MASS);	/* for CTLJCE */
//	jjk 2017-11-15		return(CL_MASS);	/* for CTLJCE */
		case CH3_ME: case CTCH3:
			return(CH3_MASS);
		case C_CH2CL2: case C_PRIM:
			return(CH2_MASS);
		case C_SECOND:
			return(CH_MASS);
		case CH_ALK:	case CH_AROM:	case CH_PRO:
		case CH_THIOL:	case CH_PYRROLE:	
		case CARBON:	case CH_THIOETH:
		case ENDC:	case CH_5INDOLE:
		case C_TERT:	case C_MTHIOETH:
		case C_GLY:	case C_ALPHA:	case C_ETHIOETH:	
		case C_THIOL:	case C_PYRROLE:	case C_PRO:	
		case C_RIBOSE:	case C_AROM:	case C_PYRIM:	
		case C_6PURINE:	case C_5PURINE:	case C_5INDOLE:	
		case C_NYL:	case C_XYL:	case C_ESTER:	
		case C_METHN:	case C_N3:	case C_THIOL2:	
		case C_T:	case COCARB:	case C_CH2I2:
		case C_CHCL3:	case C_CHBR3:	case C_CH2BR2:
			return(C_MASS);
		case NEON:
			return(NE_MASS);
		case ARGON:
			return(AR_MASS);
		case XENON:
			return(XE_MASS);
		case KRYPTON:
			return(KR_MASS);
		case N_PEP:	case N_1:	case N_PYRROLE:	
		case N_P:	case N_H2:	case N_LYS:	
		case N_H3:	case N_PYRIM:	case N_PLUS:	
		case N_T:	case NITROGEN:
			return(N_MASS);
		case S_CYS:	case S_CYSCYS:	
			return(S_MASS);
		case HELIUM:
			return(HE_MASS);
		case IODINE:	case I_CH2I2:
			return(I_MASS);
		case FLUORINE:	case FLUORIDE:
			return(F_MASS);
		case P_PHATE:	
			return(P_MASS);
		case CESIUM:
			return(CS_MASS);
		case BR_CHBR3:	case BR_CH2BR2:
			return(BR_MASS);
		case CCL4:
			return(CCL4_MASS);
		case P_PLUS:	case P_MINUS:	case NON_P:
			return(DP_MASS);
		case P_PLUSL:
			return(DP_MASS/2.);
		case PLATINUM:
			return(Pt_MASS);
		case SELEC:
			return(CSE_MASS);
		case SE:
			return(SE_MASS);
		case ethyl:
			return(C2H5_MASS);
		case pentyl:
			return(C3H7_MASS);
		case butyl:
			return(C4H9_MASS);
		case Erbium:
			return(ER_MASS);
		default:
			fprintf(stderr,"getmass: can't mass of type %d\n",atm);
			exit(1);
	}
}
