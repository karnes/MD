/*
 *	wdat:
 *
 *	This routine writes the .dat file, which we currently using to
 *	collect the positions of nCH3OH and nCH3CN molecules within 5A from each of the O(Si) atoms
 */

#include	<md.h>
#include        <math.h>
#include	<system.h>
#include	<string.h>
char h2o[3] = {'O','H','H'};
char nitb[14] = {'C','C','C','C','C','C','N','O','O','H','H','H','H','H'};

float	Wpos[1023][2];
float   NITpos[1023][2];
int numW = 200;
int numNIT = 100;

void quicksort(float[][2], int, int);

wdat(fp,initQ,pFreq)
	FILE	*fp;
	int initQ,pFreq;
{
   if (tc % pFreq != 0) return;
   int i,j; 
   int nw;
   nw = natoms - nsolute - 14*nNIT;
  fprintf(fp,"%d\n\n",numW*3+numNIT*14);
 /* assign water z pos, indices to a dummy array */
   if(nw > 0){
      for (i=0;i< nw/3;i++)
      {
         Wpos[i][0] = (float)i;
	 Wpos[i][1] = -pos[i*3].fz; //we want highest water
      }
      for(i= nw/3;i < 1023;i++)
      {
         Wpos[i][0] = (float)i;
	 Wpos[i][1] = 888.0;
      }
      quicksort(Wpos,0,1023-1);
/*
if(tc==0){
   fname++;
   fprintf(stderr,"wdat: Evap %d: highest water O is water %d, z = %f\n",fname+1000,(int)(Wpos[0][0]),-Wpos[0][1]);
}
*/
      for(i=0;i<numW;i++)
         for(j=0;j<3; j++)
//		fprintf(fp,"%c %f %f %f\n",h2o[j],pos[(int)(Wpos[i][0])*3+j].fx,pos[(int)(Wpos[i][0])*3+j].fy,pos[(int)(Wpos[i][0])*3+j].fz);
		fprintf(fp,"%d %c %f %f %f %f %f %f\n",(int)(Wpos[i][0])*3+j,h2o[j],pos[(int)(Wpos[i][0])*3+j].fx,pos[(int)(Wpos[i][0])*3+j].fy,pos[(int)(Wpos[i][0])*3+j].fz,tvel[(int)(Wpos[i][0])*3+j].fx,tvel[(int)(Wpos[i][0])*3+j].fy,tvel[(int)(Wpos[i][0])*3+j].fz);
   }
   /* assign nitrobenzene z pos, indices to a dummy array */
   if(nNIT > 0){
      for (i=0;i< nNIT; i++)
      {
         NITpos[i][0] = (float)i;
	 NITpos[i][1] = pos[nw + i*14].fz; //we want lowest NIT's
      }
      for (i= nNIT; i < 1023; i++)
      {
         NITpos[i][0] = (float)i;
	 NITpos[i][1] = 888.0;
      }
      quicksort(NITpos,0,1023-1);
//fprintf(stderr,"i think the lowest NIT is NIT %d, z = %f\n",(int)(NITpos[0][0]),NITpos[0][1]);

      for(i=0;i<numNIT; i++)
         for(j=0;j<14; j++)
	 {
//		fprintf(fp,"%c %f %f %f\n",nitb[j],pos[nw+(int)(NITpos[i][0])*14+j].fx,pos[nw+(int)(NITpos[i][0])*14+j].fy,pos[nw+(int)(NITpos[i][0])*14+j].fz);
		fprintf(fp,"%d %c %f %f %f %f %f %f\n",nw+(int)(NITpos[i][0])*14+j,nitb[j],pos[nw+(int)(NITpos[i][0])*14+j].fx,pos[nw+(int)(NITpos[i][0])*14+j].fy,pos[nw+(int)(NITpos[i][0])*14+j].fz,tvel[nw+(int)(NITpos[i][0])*14+j].fx,tvel[nw+(int)(NITpos[i][0])*14+j].fy,tvel[nw+(int)(NITpos[i][0])*14+j].fz);
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
