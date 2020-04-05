/*
 *	wdat:
 *
 *	This routine writes the .dat file, which we currently using to
 *	collect the positions of nCH3OH and nCH3CN molecules within 5A from each of the O(Si) atoms
 */

#include	<md.h>
#include	<system.h>
#include	<string.h>
float	ACNpos[1023][2];
float   Sipos[90][2];
float	MeOHpos[1023][2];
char MeOH[3] = {'O','H','C'};
char ACN[3] = {'C','N','C'};
char SiOH[3][3] = {"O","Si","H"};
void quicksort(float[][2], int, int);
int numACN = 20;
int numMeOH = 1;
int numSi = 12;
int rad2 = 49; // square of radius of cylinder
double dx,dy;

wdat(fp,initQ,pFreq)
	FILE	*fp;
	int initQ,pFreq;
{
//printf("in wdat.c \n");
int i,j,k,l,nshel;  /* total number of atoms including surface OHs*/
	if (initQ) {
//		fprintf(fp,"natoms = %d,  xwall = %f, ywall = %f, zwall = %f\n",natoms, xwall, ywall, zwall);
	}
//printf("wdat.c -- after if (initQ)\n");
	if (tc % pFreq != 0) return;
//	getShel(5.,&nshelM,&nshelA);
//printf("wdat.c -- aftergetShel\n");
//	nshel = 3*numMeOH+3*numACN+3*nSi;//nshelM*3+nshelA*3+nSi*3;
	nshel = 3*numMeOH+3*numACN+3*numSi;//nshelM*3+nshelA*3+nSi*3;
	fprintf(fp,"%d\n\n",nshel);
//	printf(" tc = %d, nshelM = %d, nshelA = %d, nSI = %d\n",tc,nshelM, nshelA,nSi);
	/* assign MeOH z pos, indices to a dummy array */
	if(nCH3OH > 0){
	for (i=0;i< nCH3OH; i++)
	{
		MeOHpos[i][0] = (float)i;
		MeOHpos[i][1] = pos[i*3].fz;
	}
	for (i= nCH3OH; i < 1023; i++)
	{
		MeOHpos[i][0] = (float)i;
		MeOHpos[i][1] = 88.0;
	}
	quicksort(MeOHpos,0,1023-1);

	for (i=0;i< numMeOH;i++)
	    for (j=0;j<3; j++)
		fprintf(fp,"%c %f %f %f\n",MeOH[j],pos[(int)(MeOHpos[i][0])*3+j].fx,pos[(int)(MeOHpos[i][0])*3+j].fy,pos[(int)(MeOHpos[i][0])*3+j].fz);
	}
	/* assign ACN z pos, indices to a dummy array */
	if(nCH3CN > 0){
	for (i=0;i< nCH3CN; i++)
	{
		ACNpos[i][0] = (float)i;
		ACNpos[i][1] = pos[nCH3OH*3 + i*3].fz;
	}
	for (i= nCH3CN; i < 1023; i++)
	{
		ACNpos[i][0] = (float)i;
		ACNpos[i][1] = 88.0;
	}
	quicksort(ACNpos,0,1023-1);

//	for (i=0;i< numACN; i++){
	i=0;
	while(k < numACN && i < nCH3CN){
		l= nCH3OH*3 + (int)(ACNpos[i][0])*3; // loop through ACNs by z 
		dx = pos[0].fx - pos[l].fx;
		dy = pos[0].fy - pos[l].fy;
		if(dx*dx + dy*dy < rad2){   // then print the ACN
		   k++;
		   for (j=0;j<3; j++)
		      fprintf(fp,"%c %f %f %f\n",ACN[j],pos[nCH3OH*3+(int)(ACNpos[i][0])*3+j].fx,pos[nCH3OH*3+(int)(ACNpos[i][0])*3+j].fy,pos[nCH3OH*3+(int)(ACNpos[i][0])*3+j].fz);
		}
		i++; // increment the ACN
	}
	}
	//similar approach for SiOH-- except sorted by radius from MeOH in x-y plane
	for(i=0;i<nSi;i++)
	{
		Sipos[i][0] = (float)i;
		dx = pos[0].fx - pos[3*nCH3OH+3*nCH3CN+3*i].fx;
		dy = pos[0].fy - pos[3*nCH3OH+3*nCH3CN+3*i].fy;
		Sipos[i][1] = dx*dx+dy*dy;
	}
	quicksort(Sipos,0,90-1);

	for (i=0;i<numSi; i++)
	    for (j=0;j<3; j++)
		fprintf(fp,"%s %f %f %f\n",SiOH[j],pos[nCH3OH*3+nCH3CN*3 +(int)(Sipos[i][0])*3+j].fx,pos[nCH3OH*3+nCH3CN*3+(int)(Sipos[i][0])*3+j].fy,pos[nCH3OH*3+nCH3CN*3+(int)(Sipos[i][0])*3+j].fz);
//printf("leave wdat.c\n");
}

void quicksort(float x[][2], int first, int last)
{
	int pivot,j,i;
	float temp0,temp1;

	if(first < last)
	{
		pivot = first;
		i = first;
		j = last;

		while(i < j)
		{
			while(x[i][1]<=x[pivot][1] && i<last)
				i++;
			while(x[j][1]>x[pivot][1])
				j--;
			if(i<j)
			{
				temp0=x[i][0];
				temp1=x[i][1];
				x[i][0]=x[j][0];
				x[i][1]=x[j][1];
				x[j][0]=temp0;
				x[j][1]=temp1;
			}
		}

		temp0=x[pivot][0];
		temp1=x[pivot][1];
		x[pivot][0]=x[j][0];
		x[pivot][1]=x[j][1];
		x[j][0]=temp0;
		x[j][1]=temp1;
		
		quicksort(x,first,j-1);
		quicksort(x,j+1,last);
	}
}
