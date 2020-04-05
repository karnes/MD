/*
 *	wdat:
 *
 *	This routine writes the .dat file, which we currently using to
 *	collect the positions of nCH3OH and nCH3CN molecules within 5A from each of the O(Si) atoms
 */

#include	<md.h>
#include	<system.h>
#include	<string.h>
tripd		shpos[3699];
char MeOH[3] = {'O','H','C'};
char ACN[3] = {'C','N','C'};
char SiOH[3][3] = {"O","Si","H"};


wdat(fp,initQ,pFreq)
	FILE	*fp;
	int initQ,pFreq;
{
//printf("in wdat.c \n");
int i,j,nshelM,	/*number of inshell CH3OH molecules*/
	nshelA,	/* number of inshel CH3CN molecules*/
	nshel;  /* total number of atoms including surface OHs*/
	if (initQ) {
//		fprintf(fp,"natoms = %d,  xwall = %f, ywall = %f, zwall = %f\n",natoms, xwall, ywall, zwall);
	}
//printf("wdat.c -- after if (initQ)\n");
	if (tc % pFreq != 0) return;
//	getShel(5.,&nshelM,&nshelA);
//printf("wdat.c -- aftergetShel\n");
	nshel = 3*nCH3OH+3*nCH3CN+3*nSi;//nshelM*3+nshelA*3+nSi*3;
	fprintf(fp,"%d\n\n",nshel);
//	printf(" tc = %d, nshelM = %d, nshelA = %d, nSI = %d\n",tc,nshelM, nshelA,nSi);
	for (i=0;i< nCH3OH;i++)
	    for (j=0;j<3; j++)
		fprintf(fp,"%c %f %f %f\n",MeOH[j],pos[i*3+j].fx,pos[i*3+j].fy,pos[i*3+j].fz);
	for (i=0;i< nCH3CN; i++)
	    for (j=0;j<3; j++)
		fprintf(fp,"%c %f %f %f\n",ACN[j],pos[nCH3OH*3+i*3+j].fx,pos[nCH3OH*3+i*3+j].fy,pos[nCH3OH*3+i*3+j].fz);
	for (i=0;i<nSi; i++)
	    for (j=0;j<3; j++)
		fprintf(fp,"%s %f %f %f\n",SiOH[j],pos[nCH3OH*3+nCH3CN*3 +i*3+j].fx,pos[nCH3OH*3+nCH3CN*3+i*3+j].fy,pos[nCH3OH*3+nCH3CN*3+i*3+j].fz);
//printf("leave wdat.c\n");
}/*
getShel(rad,nM,nA)
double rad;
int *nM,*nA;
{
int i,j,k,l,m,n;
tripd  dis;
double rad2;

rad2 = rad*rad;
*nM = 0;
*nA = 0;i
*/
/*Inshell CH3OH	*/
//for (j=0;j<nSi;j++){
/* 
 j=40;
   l = nCH3OH*3+nCH3CN*3+3*j;
   for 	(i=0; i< nCH3OH; i++) {
	k = i*3;
	dis.fx = pos[k].fx - pos[l].fx;
	dis.fy = pos[k].fy - pos[l].fy;
	dis.fz = pos[k].fz - pos[l].fz;
	mvimage(&dis);*//* current CH3OH-O(Si) atom vector*/
//printf("wdat.c -- in getShel MeOH before if\n");
/*	if 	(dis.fx*dis.fx + dis.fy*dis.fy + dis.fz*dis.fz < rad2) {
		for (m=0;m<3;m++){
			shpos[3*(*nM)+m].fx = pos[k+m].fx;
			shpos[3*(*nM)+m].fy = pos[k+m].fy;
			shpos[3*(*nM)+m].fz = pos[k+m].fz;
			
//printf("wdat.c -- in getShel MeOH if\n");
		}
	*nM += 1;
	}
//printf("wdat.c -- in getShel MeOH after if\n");
   }
//}
*/
/*In shell CH3CN*/
//for (j=0;j<nSi;j++){
/*   l = nCH3OH*3+nCH3CN*3+3*j;
   for	(i = 0; i < nCH3CN; i++){
	k = nCH3OH*3+i*3;
	dis.fx = pos[l].fx - pos[k].fx;
	dis.fy = pos[l].fy - pos[k].fy;
	dis.fz = pos[l].fz - pos[k].fz;
	mvimage(&dis);*//* current CH3CN-O(Si) atom vector*/
/*	if 	(dis.fx*dis.fx + dis.fy*dis.fy + dis.fz*dis.fz < rad2) {
		for (m=0;m<3;m++){
			shpos[3*(*nM)+3*(*nA)+m].fx = pos[m+k].fx;
			shpos[3*(*nM)+3*(*nA)+m].fy = pos[m+k].fy;
			shpos[3*(*nM)+3*(*nA)+m].fz = pos[m+k].fz;
			
		}
		*nA += 1;
	}
   }
//}
for (j=0;j<nSi;j++){
   	l = nCH3OH*3+nCH3CN*3+3*j;
	n = 3*(*nM)+3*(*nA)+j*3;
	shpos[n].fx = pos[l].fx;
	shpos[n].fy = pos[l].fy;
	shpos[n].fz = pos[l].fz;
	shpos[n+1].fx = pos[l+1].fx;
	shpos[n+1].fy = pos[l+1].fy;
	shpos[n+1].fz = pos[l+1].fz;
	shpos[n+2].fx = pos[l+2].fx;
	shpos[n+2].fy = pos[l+2].fy;
	shpos[n+2].fz = pos[l+2].fz;
}
}*/
