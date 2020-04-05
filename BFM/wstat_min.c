/*
 *	This routine will write a status line to the file pointed to by fp,
 *	for use by minimize().
 */

#include	<md.h>
#include	<system.h>

wstat_min(fp)
FILE	*fp;
{
fprintf(fp, "%-10s%-10s%-10s%-10s%-10s\n",
    "TIME", "POT",  "VLIQ", "CFMNB", "INTRAV");
fprintf(fp, "%-8.2f%-14.2f%-14.2f%-14.2f%-14.2f\n"
    ,etime,  V*KCAL, VLIQ*KCAL, CFMNB*KCAL, INTRAV*KCAL);
}
