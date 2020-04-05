/*
 *	wdat:
 *
 *	This routine writes the .dat file, which we currently using to
 *	collect ...
 */

#include	<md.h>
#include        <math.h>
#include	<system.h>
#include	<string.h>

wdat(fp,initQ,pFreq)
	FILE	*fp;
	int initQ,pFreq;

{

char h2o[3] = {'O','H','H'};
char nitb[14] = {'C','C','C','C','C','C','N','O','O','H','H','H','H','H'};
int i,j,k;
int nw;

if (tc % pFreq != 0) return;
//if (tc > 1) return;
fprintf(stderr," max snap = %5.2f\n",center_w+9.8 + width_w/2.0 - 1.1);
if(snap > (center_w + 9.8) + ((width_w/2.0) - 1.1))
   return;
else if(fabs(pos[natoms-1].fz - snap) < snap_w/5.0){ 
/* then increment snap and print a frame */

snap += snap_w;

if(nsolute==1){
        offsetWatx(natoms-1);
}


nw = natoms - nsolute - 14*nNIT;

fprintf(fp,"%d\n",natoms);
fprintf(fp,"Cl z %6.3f\n",pos[natoms-1].fz);
for(i=0;i<nw/3;i++){
   for(j=0;j<3;j++){
      k=3*i + j;
      fprintf(fp,"%c %f %f %f\n",h2o[j],pos[k].fx,pos[k].fy,pos[k].fz);
   }
}
for(i=0;i<nNIT;i++){
   for(j=0;j<14;j++){
      k=nw + i*14 + j;
      fprintf(fp,"%c %f %f %f\n",nitb[j],pos[k].fx,pos[k].fy,pos[k].fz);
   }
}

fprintf(fp,"Cl %f %f %f\n",pos[natoms-1].fx,pos[natoms-1].fy,pos[natoms-1].fz);
}

}

offsetWatx(n)
	int	n;	/* atom number from which to offset positions */
{
	int	i, j;
	tripd	image, off;
	if (n > natoms-1 || n < 0) return;
	off.fx = pos[n].fx;
	off.fy = pos[n].fy;
//	off.fz = pos[n].fz;

	for (i = 0; i < natoms; i++) {
		pos[i].fx -= off.fx;
		pos[i].fy -= off.fy;
//		pos[i].fz -= off.fz;

		if ((atom[i].flags & A_MAJOR) && atom[i].param1) {
			image.fx = -pos[i].fx;
			image.fy = -pos[i].fy;
//			image.fz = -pos[i].fz;
			mvimage(&pos[i]);
			image.fx += pos[i].fx;
			image.fy += pos[i].fy;
//			image.fz += pos[i].fz;
			for (j=i+1; j <= i + atom[i].param1; j++) {
				pos[j].fx += image.fx;
				pos[j].fy += image.fy;
//				pos[j].fz += image.fz;
			}
		}
		else if (atom[i].flags & A_MAJOR)
			mvimage(&pos[i]);
	}
}
