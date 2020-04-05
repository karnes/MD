/*
 *	This routine will write a status line to the file pointed to by fp,
 */

#include	<md.h>
#include	<system.h>
#include	<water.h>

wstat(fp)
    FILE	*fp;
{
    // adding solvation shell populations
    fprintf(fp, "%-19s%-12s%-12s%-11s%-12s%-12s%-12s%-12s%-12s%-12s%-12s%-12s%-10s%-10s%-10s%-10s%-9s","TIME", "TEMP", "ETOT", "INTRADDC","INTRAH2O","Er-DDC","Er-H2O","H2O-DDC","H2O-NB","DDC-NB","Er-z","sysCoM","wShl1","wShl2","V_bias","V_win","w_cor");
    fprintf(fp, "\n");
    fprintf(fp,"%-18.1f%-12.2f%-12.2f%-12.2f%-12.2f%-12.2f%-12.2f%-12.2f%-12.2f%-12.2f%-12.3f%-12.3f%-10d%-10d%-10.3f%-10.3f%-9.3f",etime,temp,E*KCAL,INTRADDC*KCAL,WATERV*KCAL,V_ED*KCAL,V_EW*KCAL,V_WD*KCAL,VWATTS*KCAL,DDCNB*KCAL,Erz,osysCmz,shellO_1,shellO_2,V_bias*KCAL,V_w*KCAL,wcor);
    fprintf(fp, "\n");

    if(nEr==1){
	offsetEr(natoms-1);
    }

}


offsetEr(n)
    int	n;	/* atom number from which to offset positions */
{
    int	i, j;
    tripd	image, off;
    if (n > natoms-1 || n < 0) return;
    off.fx = pos[n].fx;
    off.fy = pos[n].fy;
    if (pbcType[0] == 'O'){
	off.fz = pos[n].fz;
    }
    for (i = 0; i < natoms; i++) {
	pos[i].fx -= off.fx;
	pos[i].fy -= off.fy;
	if (pbcType[0] == 'O'){
	    pos[i].fz -= off.fz;
	}

	if ((atom[i].flags & A_MAJOR) && atom[i].param1) {
	    image.fx = -pos[i].fx;
	    image.fy = -pos[i].fy;
	    if (pbcType[0] == 'O'){
		image.fz = -pos[i].fz;
	    }
	    mvimage(&pos[i]);
	    image.fx += pos[i].fx;
	    image.fy += pos[i].fy;
	    if (pbcType[0] == 'O'){
		image.fz += pos[i].fz;
	    }
	    for (j=i+1; j <= i + atom[i].param1; j++) {
		pos[j].fx += image.fx;
		pos[j].fy += image.fy;
		if (pbcType[0] == 'O'){
		    pos[j].fz += image.fz;
		}
	    }
	}
	else if (atom[i].flags & A_MAJOR)
	    mvimage(&pos[i]);
    }
}
