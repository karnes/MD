/* wnFC.c		read rsp files, determine if water finger broken, 
 * 			Cl- ion windowed runs at equilibrium. 
 * 									*/
#include	<stdio.h>
#include	<stdlib.h>
#include	<typedefs.h>
#include	<atomtypes.h>
#include	<globals.h>
#include	<time.h>
//#define BIG  100
#define abs(x)          ((x)>0.?(x):-(x))
#define sq(x)   ((x)*(x))
#define hsgn(x)         ((x)<0.?-.5:.5)
#include	<math.h>
#include	<string.h>
#define numWindows 21 // num of 3Å windows to calculate probability distributions
#define watMinZ -1.5
#define maxWat 1200
#define stWin 10

main(){

int i,j,k,p;
char inFile[80],buf[360];
FILE *fp;
int allocFlag = 0;
float Wpos[maxWat][2];
int ii,jj;
int nw,nwD;
double Zgap,maxZgap;
int ionflag;
tripd sdist;
double wcor;
double maxZ[10*(numWindows-stWin)][10]={0.0};
/*
*/
char Cl[22] = {'H','G','F','E','D','C','B','A',
	'0','1','2','3','4','5','6','7','8','9',
	'R','S','T','U'};

for(i=stWin;i<numWindows;i++){
   for(j=0;j<10;j++){
	   for(k=0;k<10;k++){
		sprintf(inFile,"/home/jkarnes/MD/WNIT/run/data/Cl%c/Cl%c000%d%d.rsp",Cl[i],Cl[i],j,k);
		fprintf(stderr,"open file %s\n",inFile);
		fp = fopen(inFile,"r");

/*
 *  Read the header data, allocate space for the arrays, and read them too.
 */
		fread(filetype, sizeof(char), 5, fp);
		fread(&datestamp, sizeof(datestamp), 1, fp);
		fread(status, sizeof(char), 4, fp);
		fread(&natoms, sizeof(natoms), 1, fp);
		fread(&nsolute, sizeof(nsolute), 1, fp);
		fread(&xwall, sizeof(xwall), 1, fp);
		fread(&ywall, sizeof(ywall), 1, fp);
		fread(&zwall, sizeof(zwall), 1, fp);
		fread(&EqTemp, sizeof(EqTemp), 1, fp);
		fread(&DEqTemp, sizeof(DEqTemp), 1, fp);
		fread(&EqPress, sizeof(EqPress), 1, fp);
		fread(&DEqPress, sizeof(DEqPress), 1, fp);
		
		if(allocFlag==0){
			if ((atom = (parts *) calloc(natoms, sizeof(atom[0]))) == NULL ||
		    	   (pos  = (tripd *) calloc(natoms, sizeof(pos[0]))) == NULL) {
				fprintf(stderr, "psxyz: out of memory\n");
				exit(1);
			}
			allocFlag=1;
		}

		fread(atom, sizeof(parts), natoms, fp);
		fread(&xtrInQ, sizeof(xtrInQ), 1, fp);
		fread(pos, sizeof(tripd), natoms, fp);
	
	        fclose(fp);
		
		
 		//get waters of interest
		nw = 2958;//natoms-14*nNIT-nsolute;
	nwD = 0;
	wcor = 0.0;
	//sort waters by oxygen z position
	for(ii=0;ii<nw/3;ii++)
	{
   	   Wpos[ii][0]=(float)(ii*3);
   	   Wpos[ii][1]=-pos[ii*3].fz;
     	  if(pos[ii*3].fz > watMinZ){
        	  nwD++;
       	  }
        }
	for(ii=nw/3;ii<maxWat;ii++){
	   Wpos[ii][0] = (float)(ii*3);
	   Wpos[ii][1] = 888.0;
	}
	qSort(Wpos,0,maxWat-1);
	// first attempt: max z separation -> wcor
	// step through waters by decreasing z.
	ii=0;
	//find solute
	while(pos[natoms-1].fz < pos[(int)(Wpos[ii][0])].fz){
	   ii++;
	}
	maxZgap = pos[natoms-1].fz - pos[(int)(Wpos[ii][0])].fz;
//	fprintf(stderr,"ORIG maxZgap = %f\n",maxZgap);
	ionflag=8888;
	for(jj=ii;jj<nwD;jj++) // perhaps replace with "while z > 0.0" and step through
	// although the "move to closest until z < 0.0" sounds pretty good
	{
	   Zgap=pos[(int)(Wpos[jj][0])].fz - pos[(int)(Wpos[jj+1][0])].fz;
	   if(Zgap>maxZgap){
		   maxZgap=Zgap;
//fprintf(stderr,"NEW maxZgap = %f\n",maxZgap);
		   ionflag=jj;
	   }
        }
	if(ionflag==8888)
	{
	     jj=natoms-1;
	     ii=(int)(Wpos[ii][0]);
	}
	else
	{
	   jj=(int)(Wpos[ionflag][0]);
	   ii=(int)(Wpos[ionflag+1][0]);
   	}
	sdist.fx = pos[ii].fx - pos[jj].fx;
	sdist.fy = pos[ii].fy - pos[jj].fy;
	sdist.fz = pos[ii].fz - pos[jj].fz;
//	mvimage(&sdist);
	wcor = sqrt(sdist.fx*sdist.fx +  sdist.fy*sdist.fy +  sdist.fz*sdist.fz);
//	fprintf(stderr,"wcor = %f\n",wcor);
        maxZ[10*(i-stWin)+j][k]=maxZgap;

    }
  }
}

sprintf(inFile,"wnERunFinger.txt");
fp = fopen(inFile,"w");

fprintf(fp,"run\t");
for(i=stWin;i<numWindows;i++){
	for(j=0;j<10;j++){
		fprintf(fp,"Cl%c_%d\t",Cl[i],j);
	}
}
fprintf(fp,"\n");
for(k=0;k<10;k++){
	fprintf(fp,"%d\t",k);
	for(i=stWin;i<numWindows;i++){
		for(j=0;j<10;j++){
			fprintf(fp,"%f\t",maxZ[10*(i-stWin)+j][k]);
		}
	}
	fprintf(fp,"\n");
}
fclose(fp);

}


void qSort(float x[][2], int first, int last)
{
	int pivot, i, j;
	float temp0, temp1;
	if(first < last)
	{
		pivot = first;
		i = first;
		j = last;

		while(i < j)
		{
			while(x[i][1] <= x[pivot][1] && i < last)
				i++;
			while(x[j][1] > x[pivot][1])
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

		qSort(x,first,j-1);
		qSort(x,j+1,last);
	}
}
     
/*
 *	This routine does the periodic imaging of any given
 *	vector for the parallelopiped or the truncated octahedral
 *      boundary conditions.
 */

mvimage(del)
	tripdouble	*del;
{
	double	fx, fy, fz;
	tripd period, iperiod;
	if (natoms == nsolute) /* no solvent */
	     return;
	period.fx = 2.0 * xwall;	/* box period (boundary conditions) */
	period.fy = 2.0 * ywall;	/* box period (boundary conditions) */
	period.fz = 2.0 * zwall;	/* box period (boundary conditions) */
	iperiod.fx = 1.0 / period.fx;		/* inverse period */
	iperiod.fy = 1.0 / period.fy;		/* inverse period */
	iperiod.fz = 1.0 / period.fz;		/* inverse period */
	fx = del->fx * iperiod.fx;  /* Normalizing the positions */
	fy = del->fy * iperiod.fy;  /* |fx|, |fy|, |fz| <= 0.5 for the  */
	fz = del->fz * iperiod.fz;  /* atom to be in the box */
	fx -= (int) fx;
	fy -= (int) fy;   /* bringing the particle to the nearest box */
	fz -= (int) fz;
	fx -= (int) (2 * fx);
	fy -= (int) (2 * fy);  /* bringing the particle to the box */
	fz -= (int) (2 * fz);
	/* truncated octahedral boundaries */
		if	((abs(fx) + abs(fy) + abs(fz)) > 0.75) {
			fx -= hsgn(fx);
			fy -= hsgn(fy);
			fz -= hsgn(fz);
		}

	del->fx = fx * period.fx;
	del->fy = fy * period.fy;
	del->fz = fz * period.fz;
}
