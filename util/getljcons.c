#include	<md.h>

	/*	ljcon	lj[3][3] is declared in md.h.
	 *      lj[0][0] contains the solvent-solvent parameters.
	 *      lj[1][j] contains the interaction between the solvent and
	 *      solute atom j
	 */
	 
getljcons()
{
	int	liquid_type, system_type, j;
	ljcon	*(lookup());

	if (natoms > nsolute)
	   {
		if (nsolute > 3)
	   	  ERROR((stderr,"getljcons: ljcons defined for 3 solutes\n"), exit);
		liquid_type = atom[0].type;
		lj[0][0] = *(lookup(liquid_type, liquid_type));
		for (j = 0; j < nsolute; j++)
		    {
		     system_type = atom[natoms-nsolute+j].type;
		     lj[1][j] = *(lookup(liquid_type,system_type));
                    }
	   }
}
