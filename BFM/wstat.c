/*
 *	This routine will write a status line to the file pointed to by fp,
 */

#include	<md.h>
#include	<system.h>

wstat(fp)
FILE	*fp;
{
fprintf(fp,"%-7s%-10s%-10s%-8s%-8s\n","TIME","TEMP","ETOT","INTRAV","CFMNB");
fprintf(fp,"%-8.1f%-9.2f%-9.2f%-10.2f%-10.2f\n",etime,temp,E*KCAL,INTRAV*KCAL,CFMNB*KCAL);
}
