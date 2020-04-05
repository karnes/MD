#include	<md.h>

pimage(del)
tripdouble	*del;
	{
	double	fx, fy, fz;
	if (natoms == nsolute) return;
	fx = del->fx * iperiod.fx;
	fy = del->fy * iperiod.fy;
	fz = del->fz * iperiod.fz;
	fx -= (int) (2 * fx);
	fy -= (int) (2 * fy);
	fz -= (int) (2 * fz);
	if (pbcType[0] == 'O'){
	/* truncated octahedral boundaries */
		if	((abs(fx) + abs(fy) + abs(fz)) > 0.75) {
			fx -= hsgn(fx);
			fy -= hsgn(fy);
			fz -= hsgn(fz);
		}
	}
	else if (pbcType[0] != 'C')
		ERROR((stderr,"pimage is set for either OCT or CUBE\n"), exit);

	del->fx = fx * period.fx;
	del->fy = fy * period.fy;
	del->fz = fz * period.fz;
	}
