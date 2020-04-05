/*
 *	wdatBCD.c:
 *	This routine writes the .dat file, writes BCD molecule xyz file.
 */
#include	<md.h>
#include	<system.h>
wdat(fp,initQ,pFreq)
	FILE	*fp;
	int initQ,pFreq;
{
void quicksort(float[][2], int, int);
int i,j,k;
int nw,nb,n2p,nsolvent;
tripd center,sdl;
double r, totmass;
char symw[3] = {'O','H','H'};
char sym[21] = {'O',
		'C','H','C','H',
		'C','H','C','H',
		'C','H','C','H',
		'H','O','H','O',
		'O','H','O','H'};
nw = natoms - nBCD*BCDs - nBrO*BrOs - nsolute;
nb = natoms - nBCD*BCDs - nw - nsolute;
nsolvent = nw + nb;

fprintf(fp,"%d\n\n",natoms);
for(i=0;i<nw/3;i++){
  for(j=0;j<3;j++){
    k=3*i+j; 
    fprintf(fp,"%c %f %f %f\n",symw[j],pos[k].fx,pos[k].fy,pos[k].fz);
  }
}

for(i=0;i<nb/BrOs;i++){
  for(j=0;j<BrOs;j++){
    k=nw+BrOs*i+j;
    if(j==1)
	fprintf(fp,"Br %f %f %f\n",pos[k].fx,pos[k].fy,pos[k].fz);
    else
	fprintf(fp,"C %f %f %f\n",pos[k].fx,pos[k].fy,pos[k].fz);
  }
}

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
