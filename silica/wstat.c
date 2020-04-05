/*
 *	This routine will write a status line to the file pointed to by fp,
 */

#include	<md.h>
#include        <system.h>
#include        <water.h>
wstat(fp)
	FILE	*fp;
{
 fprintf(fp, "%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s\n",
 "TIME","TEMP","ETOT","VLIQ","INTRAV_S","VNB_S","INTRAV_M","VNB_M","VMS","INTRAV_A","VNB_A","VAS","VMA","HB_StoM","HB_MtoS","HB_StoA");
	fprintf(fp,"%-9.1f%-9.2f%-9.2f%-9.2f%-9.2f%-9.2f%-9.2f%-9.2f%-9.2f%-9.2f%-9.2f%-9.2f%-9.2f%-9d%-9d%-9d\n"
    ,etime,temp,E*KCAL,VLIQ*KCAL,INTRAV_S*KCAL,VNB_S*KCAL,INTRAV_M*KCAL,VNB_M*KCAL,VMS*KCAL,INTRAV_A*KCAL,VNB_A*KCAL,VAS*KCAL,VMA*KCAL,HBSitoM,HBMtoSi,HBSitoA);

// print distance stats if there's only one methanol
/*if(nCH3OH == 1){
 fprintf(fp, "%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s%-9s\n",
 "TIME","MHtoSiO","SiOtoMOd","cosMd","iSiOMH","SiAtoA","MOtoSiH","SiOtoMOa","cosMa","iSiHMO","SiDtoA","MeOHd","durMeOHd","MeOHa","durMeOHa","grepFlag");
	fprintf(fp,"%-9.1f%-9.4f%-9.4f%-9.5f%-9d%-9.4f%-9.4f%-9.4f%-9.5f%-9d%-9.4f%-9d%-9.3f%-9d%-9.3f%-9d\n"
    ,etime,MHtoSiO,SiOtoMOd,cosMd,SiOMH,SiAtoA,SiHtoMO,SiOtoMOa,cosMa,SiHMO,SiDtoA,MeOHd,durMeOHd,MeOHa,durMeOHa,88888);
}
*/
}
/*
VINTSA, VINTSM, VINTSS solute interactions with silica and solvents
INTRAV_S, VNB_S	silica intramolecular and intermolecular energies
INTRAV_M, VNB_M	methanol intramolecular and intermolecular energies
INTRAV_A, VNB_A	acetonitrile intramolecular and intermolecular energies
VMA		acetonitrile-methanol interaction
VMS		SiOH-methanol interaction
VAS		SiOH-acetonitrile interaction
HB_StoM		H bonds: Si donor, methanol acceptor
HB_MtoS		H bonds: methanol donor, Si acceptor
HB_StoA		H bonds: Si donor, acetonitrile acceptor
*/
