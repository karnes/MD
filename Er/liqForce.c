#include	<md.h>
#include        <water.h>
#include	<system.h>

liqForce(pos, force)
    tripd	*pos;
    tripd	*force;
{
    int i,nw;
    double ER, Reflect();

    nw = natoms - nEr - nDDC*DDCsites;
    tc++;
    VLIQ = V_EW = V_ED = DDCNB = V_WD = 0.0;
    EW_C = ED_C = DD_C = WD_C = VINT = 0.0;
    INTRADDC = VWATTS = WATERV = H2OC = 0.0;
    shellO_1 = shellO_2 = 0;

    if(nEr==1){
	ErForce();
	VLIQ += V_EW + V_ED;
	Erz = pos[natoms-1].fz;
    }
    if(nDDC>0){
	DDCForce();
	VLIQ += DDCNB + INTRADDC;
    }
    if(nw>0){
	waterForce(pos,force,nw);
	VLIQ += WATERV + VWATTS;
	if(nDDC>0){
	    waterDDC(pos,force);
	    VLIQ += V_WD;
	}
    }
// soft reflecting wall and water finger coordinate aded by john k 2018-08-27
	if(nw > 9 && nDDC > 0){
	      ER=0.0;
              ER=Reflect();
	      VLIQ += ER;
	      if(ER > 0.0) 
		fprintf(stderr," ER = %f\n",ER*KCAL);
	}
	if(dataRatex > 0){ // (EQU runs crash if % 0 tested)
		if(tc % dataRatex == 0){
	    		wcor = 0.0;
	    		if(nEr == 1 && nw > 9){
		   	    getWcoord();
	        	}
		}
	}

}



// breadth first search routine to get wcor
#define maxWat 5000
void getWcoord(void)
{
   //get waters of interest
   double watMinZ;
   int left,right,mid,irwf;
   double Wpos[maxWat][2]; // 0 = index, 1 = z 
   double **dist;
   double *distList;
   int **graph;
   int i,j,rr,ff,all,flag;
   int nw,nwD;
   int ele,numAbove;
   int *queue, *visited;
   int fbulk; // first index of bulk water
   tripd sdist;
   nw = natoms-DDCsites*nDDC-nEr;
   nwD = 0;
   tripd ivec, jvec;
   tripd image, sdl;
   double r2,delfx,delfy,delfz;
   //sort waters by oxygen z position
   if(pos[natoms-1].fz > 0.0){
	   watMinZ = -6.0;
   }
   else{
	   watMinZ=pos[natoms-1].fz - 3.4;
   }

   for(i=0;i<nw/3;i++)
   {
      Wpos[i][0]=(float)(i*3);
      Wpos[i][1]=-pos[i*3].fz;
      if(pos[i*3].fz > watMinZ){
	  nwD++;
      }
   }
   for(i=nw/3;i<maxWat;i++){
      Wpos[i][0] = (float)(i*3);
      Wpos[i][1] = 888.0;
   }
   qSort(Wpos,0,maxWat-1);
   //find solute. this loop eliminates water above the ion.
   numAbove=0;
   while(pos[natoms-1].fz < pos[(int)(Wpos[numAbove][0])].fz){
	   numAbove++;
   }
   // allocate distance matrix,graph,etc.
   ele = nwD - numAbove + nEr;
//fprintf(stderr,"numAbove = %d, nwD = %d, ele = %d\n",numAbove, nwD, ele);
   if((dist  = calloc(ele, sizeof(double *))) == NULL ||
      (distList = calloc(ele*ele, sizeof(double))) == NULL ||
      (graph = calloc(ele, sizeof(int *))) == NULL || 
      (queue = calloc(ele, sizeof(int))) == NULL ||
      (visited = calloc(ele, sizeof(int))) == NULL)
	  ERROR((stderr,"liqForce.c: out of core. wcor array initilization.\n"),exit);
   for(i=0;i<ele;i++){
   	if(( dist[i] = calloc(ele, sizeof(double))) == NULL ||
   	   (graph[i] = calloc(ele, sizeof(int))) == NULL)
	   ERROR((stderr,"liqForce.c: out of core. wcor array initilization.\n"),exit);
   }
   //populate dist matrix, initialize others
   fbulk=9999;
   for(i=0;i<ele;i++){
	   if(i==ele-1){
	         ivec.fx=pos[natoms-1].fx;
	         ivec.fy=pos[natoms-1].fy;
	         ivec.fz=pos[natoms-1].fz;
	   }
	   else{
		 ivec.fx=pos[(int)Wpos[numAbove+i][0]].fx;
		 ivec.fy=pos[(int)Wpos[numAbove+i][0]].fy;
		 ivec.fz=pos[(int)Wpos[numAbove+i][0]].fz;
		 if(pos[(int)Wpos[numAbove+i][0]].fz<=0 && fbulk==9999){
			fbulk=i;
		 }
	   }
	   for(j=0;j<ele;j++){
		   if(j==ele-1){
		         jvec.fx=pos[natoms-1].fx;
		         jvec.fy=pos[natoms-1].fy;
		         jvec.fz=pos[natoms-1].fz;
		   }
		   else{
			 jvec.fx=pos[(int)Wpos[numAbove+j][0]].fx;
			 jvec.fy=pos[(int)Wpos[numAbove+j][0]].fy;
			 jvec.fz=pos[(int)Wpos[numAbove+j][0]].fz;
		   }
                   /***	Determine image vector for i -  j ***/
		   image.fx = -(sdl.fx = ivec.fx - jvec.fx);
		   image.fy = -(sdl.fy = ivec.fy - jvec.fy);
		   image.fz = -(sdl.fz = ivec.fz - jvec.fz);
		   mvimage(&sdl);
		   image.fx += (delfx =sdl.fx);
		   image.fy += (delfy =sdl.fy);
		   image.fz += (delfz =sdl.fz);
		   r2 = delfx*delfx + delfy*delfy + delfz*delfz;
		   if(i==ele-1 || j==ele-1){
		   	distList[i*ele+j]=dist[i][j]=/*-0.4*/0.5+sqrt(r2); // adjust for ion-water g(r) here
		   }
	 	   else{
		   	distList[i*ele+j]=dist[i][j]=sqrt(r2);
		   }
		   graph[i][j]=0;
	   }
   }
   qSort2(distList,0,(ele*ele)-1); 
   // finished setup
   left = 0;
   right = ele*ele;
   while(left < right){
	   mid = (left+right)/2;
	   makeGraph(ele,dist,graph,distList[mid]);
    	   for(i=0;i<ele;i++){
	      queue[i] = visited[i] = 0;  //reset list
   	   }
   	   ff=0;rr=-1;
	   visited[ele-1]=1;
   	   bfs(ele,ele-1,graph,queue,visited,rr,ff);
   	   all=1;		
    	   for(i=0;i<ele;i++){
	      if(visited[i]==0){
	          all=0;
	      }
	      if(i >= fbulk && i != (ele-1)){
		  if(visited[i]==1){
	 	      all=1;
		      i=ele+1;
		  }
	      }
   	   }
	   if(all==1){
		irwf=mid;
		right=mid-1;
	   }
	   else{
		left=mid+1;
	        irwf=left;
	   }
   }
   wcor = distList[irwf];
       
   for(i=0;i<ele;i++){
	   free(dist[i]);
	   free(graph[i]);
   }
   free(dist);
   free(graph);
   free(distList);
   free(queue);
   free(visited);
}

void makeGraph(int ele,double **dist,int **graph,double maxDist)
{
	int i,j;
	for(i=0;i<ele;i++){
	   for(j=0;j<ele;j++){
	      if(i==j){
		   graph[i][j]=0;
	      }
	      else if(dist[i][j] < maxDist)
	      {
		   graph[i][j]=1;
	      }
	      else
	      {
		   graph[i][j]=0;
	      }
	   }
	}
}

void bfs(int ele, int v, int **a, int *q, int *visited, int r, int f) 
{
	int i;
	for(i=0;i<ele;i++){
		if(a[v][i]==1 && visited[i]==0){
			q[++r]=i;
		}
	        if(f<=r){
		   visited[q[f]]=1;
		   bfs(ele,q[f++],a,q,visited,r,f);
		}
	}
}

void qSort2(double *x, int first, int last)
{
	int pivot, i, j;
	float temp;
	if(first < last)
	{
		pivot = first;
		i = first;
		j = last;

		while(i < j)
		{
			while(x[i] <= x[pivot] && i < last)
				i++;
			while(x[j] > x[pivot])
				j--;
			if(i<j)
			{
				temp=x[i];
				x[i]=x[j];
				x[j]=temp;
			}
		}
		temp=x[pivot];
		x[pivot]=x[j];
		x[j]=temp;

		qSort2(x,first,j-1);
		qSort2(x,j+1,last);
	}
}
void qSort(double x[][2], int first, int last)
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
     

// soft reflecting wall added by john k 2015-06-19
#define ExCon 0.5
double Reflect(){
	int i,j,nw;
	double cmposZ,ztop,zbot;
	double RFIELD;
	RFIELD = 0.0;

	nw = natoms - DDCsites*nDDC - nEr;

	ztop = zwall - 5; // reflecting wall 5A from top of box
	zbot = -zwall+ 5; // reflecting wall 5A from bottom of box

	for(i=0;i<nw;i=i+3){
		cmposZ = 0.0;
		for(j=0;j<3;j++){
			cmposZ += mass[i+j]*pos[i+j].fz;
		}
		cmposZ /= WMass;
		if(cmposZ > ztop){
			fprintf(stderr,"water %d above. cmposZ = %f\n",i/3,cmposZ);
			RFIELD += ExCon*(cmposZ-ztop)*(cmposZ-ztop);
			for(j=0;j<3;j++){ //force > 0
				force[i+j].fz -= 2*ExCon*mass[i+j]*(cmposZ-ztop)/WMass;
			}
		}
		else if(cmposZ < zbot){
			fprintf(stderr,"water %d below. cmposZ = %f\n",i/3,cmposZ);
			RFIELD += ExCon*(cmposZ-zbot)*(cmposZ-zbot);
			for(j=0;j<3;j++){ //force < 0
				force[i+j].fz -= 2*ExCon*mass[i+j]*(cmposZ-zbot)/WMass;
			}
		}
	}
	for(i=0;i< nDDC*DDCsites; i=i+DDCsites){
		cmposZ = 0.0;
		for(j=0;j<DDCsites;j++){
				cmposZ += mass[nw+i+j]*pos[nw+i+j].fz;
		}
		cmposZ /= DDCmass;
		if(cmposZ > ztop){
			fprintf(stderr,"DDC %d above. cmposZ = %f\n",i/DDCsites,cmposZ);
			RFIELD += ExCon*(cmposZ-ztop)*(cmposZ-ztop);
			for(j=0;j<DDCsites;j++){
				force[nw+i+j].fz -= 2*ExCon*mass[nw+i+j]*(cmposZ-ztop)/DDCmass;
			}
		}
		else if(cmposZ < zbot){
			fprintf(stderr,"DDC %d below. cmposZ = %f\n",i/DDCsites,cmposZ);
			RFIELD += ExCon*(cmposZ-zbot)*(cmposZ-zbot);
			for(j=0;j<DDCsites;j++){
				force[nw+i+j].fz -= 2*ExCon*mass[nw+i+j]*(cmposZ-zbot)/DDCmass;
			}
		}
	}
	if(nEr==1){
		cmposZ = pos[natoms-1].fz;
		if(cmposZ > ztop){
			fprintf(stderr,"ion above. cmposZ = %f\n",cmposZ);
			RFIELD += ExCon*(cmposZ-ztop)*(cmposZ-ztop);
			force[natoms-1].fz -= 2*ExCon*mass[natoms-1]*(cmposZ-ztop)/ErMass;
		}
		else if(cmposZ < zbot){
			fprintf(stderr,"ion below. cmposZ = %f\n",cmposZ);
			RFIELD += ExCon*(cmposZ-zbot)*(cmposZ-zbot);
			force[natoms-1].fz -= 2*ExCon*mass[natoms-1]*(cmposZ-zbot)/ErMass;
		}
	}
	return(RFIELD);
}


