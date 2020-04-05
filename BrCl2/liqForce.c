#include	<md.h>
#include        <water.h>
#include	<system.h>

liqForce(pos, force)
	tripd	*pos;
	tripd	*force;
{
   int i;
	
   tc++;
   
   VLIQ = Cl2BondE = INTER_X = X_C = 0.0;
   
   if(nCl2==1){
      Cl2Force();
      VLIQ+=Cl2BondE+INTER_X;
   }  
//fprintf(stderr,"liqForce.c complete\n");
}
