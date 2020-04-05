#include	<md.h>
#include	<system.h>

/* calculate the electric forces on all the atoms due to a constant electric
 * field */

double ElecPot(pos, force)
	tripd	*pos;
	tripd	*force;
{
	int	i, index,nw;
	double pe;

	pe = 0.;
	EVNIT = EVH2O = EVION[0] = EVION[1] = 0.;
	nw = natoms-nsolute-14*nNIT;/* number of water atoms*/
	if (EFieldQ[0] > 0){
		for (i = 0;i < nw;i++){
			index = i%3+14;/* index = 14 for O and 15,16 for H*/
			force[i].fz += EForce[index];
			pe -= EForce[index]*pos[i].fz;
			EVH2O -= EForce[index]*pos[i].fz;
		}
	}
	if (EFieldQ[1] > 0){
		for (i = nw;i < natoms-nsolute;i++){
			index = (i-nw)%14;
			force[i].fz += EForce[index];
			pe -= EForce[index]*pos[i].fz;
			EVNIT -= EForce[index]*pos[i].fz;
		}
	}
	if (EFieldQ[2] > 0){
		for (i=natoms-nsolute; i<natoms;i++){
			index = i-natoms+nsolute+17;
			force[i].fz += EForce[index];
			pe -= EForce[index]*pos[i].fz;
			EVION[i-natoms+nsolute] = -EForce[index]*pos[i].fz;
		}
	}
	return(pe);
}
