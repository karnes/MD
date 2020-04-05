/*
 *	This routine will write a status line to the file pointed to by fp,
 *	for use by minimize().
 */

#include	<md.h>
#include	<system.h>

wstat_min(fp)
FILE	*fp;
{
fprintf(fp, "%-10s%-12s%-12s%-10s%-14s\n",
    "TIME", "POT",  "VLIQ", "BCDNB", "INTRABCD");
fprintf(fp, "%-8.2f%-14.2f%-14.2f%-14.2f%-14.2f\n"
    ,etime,  V*KCAL, VLIQ*KCAL, 0*KCAL, INTRABCD*KCAL);
}
