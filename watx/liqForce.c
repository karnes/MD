#include	<md.h>
#include        <water.h>
#include	<system.h>

void qSort(double[][2],int,int);
void qSort2(double*,int,int);

liqForce(pos, force)
tripd	*pos;
tripd	*force;
{
	
int i,j,nw,endRun;
double ER,r;
double Reflect();
tripd dist;
// CENTER OF MASS OF WATER TO BE TRANSFERRED        
watx.fx = (pos[0].fx*mass[0] + pos[1].fx*mass[1] + pos[2].fx*mass[2])/WMass;
watx.fy = (pos[0].fy*mass[0] + pos[1].fy*mass[1] + pos[2].fy*mass[2])/WMass;
watx.fz = (pos[0].fz*mass[0] + pos[1].fz*mass[1] + pos[2].fz*mass[2])/WMass;
/*
if (natoms-nsolute == 0)
   return;
*/
nw = natoms-nsolute-14*nNIT;

tc++; 
static int liqAlloc = 1;
if(tc==0){
   snap = center_w + 9.8 - (width_w/2.0);
//   fprintf(stderr,"liqForce.c: 'snapshot start' z = %5.2f\n",snap);
   initDataMatrices(liqAlloc);
   if(liqAlloc){
      fname=0;
      getTeth(nw);
   }
   else{	   
      for(i=0;i<fixShel;i++){
         fprintf(stderr,"liqForce.c - keep water %d tethered to ion\n",fixw[i]);
      }
   }
}
liqAlloc = 0;
cMer=0;
Umer = Umer0 = Umer1 = 0.0;
// killfinger energy
kfnit = kfwat = 0.0;
//last h-bond information
rOO = ADmin = 999.9;
AorD = 2;
// NITROBENZENE
VINTx = VINT_NIT = VC_NIT = 0.0;
NITV = NITNB = NITC = 0.0;
NITROForce();
VLIQ = NITV + NITNB;
// WATER
H2ONITV = VWATTS = WATERV = H2OC = WNITC = 0.0;
VCH2O[0] = 0.;/* first shell water eletrostatic*/
VCH2O[1] = 0.;/* all other water electrostatic*/
VINTH2O = 0.;/*total interaction with all water*/
VINTH2O_s = 0.;/*total interaction with 1st shell water*/
V_teth = 0.0; /*total energy tethering fixShel waters */
Nshel = Nshel2 = 0; /*number of water molecules in first shell*/
if(nw > 0){
   waterForceXWat(pos, force, nw);
   VLIQ += VWATTS + WATERV;
   waterNITRO(pos,force);
   VLIQ += H2ONITV;
   xWatForce(pos, force);
}
if(KillFinger==1){
   extForce(nw,nNIT*14);
   VLIQ += EXFIELD;
}
if(nw > 9 && nNIT > 0){
   ER=0.0;
   ER=Reflect();
   VLIQ += ER;
   if(ER > 0.0) 
      fprintf(stderr,"liqForce.c: ER = %f\n",ER*KCAL);
}

VINTx  = VINT_NIT + VINTH2O;

if(dataRatex > 0){ // (EQU runs crash if % 0 tested)
   if(k10flag==1)
   {
      endRun = 0;
      for(i=0;i<nw;i+=3)
      {
	 if(pos[i].fz >= 10.0)
	 {
	    endRun = 1; 
	    for(j=0;j<nw;j+=3)
	    {
	       dist.fx = pos[i].fx - pos[j].fx; 
	       dist.fy = pos[i].fy - pos[j].fy; 
	       dist.fz = pos[i].fz - pos[j].fz;
	       mvimage(&dist);
	       r = sqrt(dist.fx*dist.fx+dist.fy*dist.fy+dist.fz*dist.fz);
	       if(r<6.0 && i!=j)
		  endRun = 0;
	    }
	    if(endRun==1){
		fprintf(stderr,"liqForce.c: atom %d (water O) z = %f, rOO min > 6.0 Å ... terminating run.\n",i,pos[i].fz);
		exit(0);
	    }
	 }
      }
   }
   else if(tc % dataRatex == 0){
      wcor = 0.0;
      if(nsolute == 0 && nw > 9){
         getWcoord();
      }
   }
}
	
}

// breadth first search routine to get wcor
#define maxWat 1200
void getWcoord(void)
{
   //get waters of interest
//   double watMinZ;
   int left,right,mid,irwf;
   double Wpos[maxWat][2]; // 0 = index, 1 = z 
   double **dist;
   double *distList;
//   double distList2[987*987][2]; //jjk 
   int **graph;
   int i,j,rr,ff,all,flag;
//   int ii,jj;
   int nw,nwD;
   int ele,numAbove;
   int *queue, *visited;
   int fbulk; // first index of bulk water
   tripd sdist;
   nw = natoms-14*nNIT-nsolute;
   nwD = 1;
   tripd ivec, jvec;
   tripd image, sdl;
   double r2,delfx,delfy,delfz;
   // For water transfer, use xwat's oxygen atom as reference 
   // so we can refer to hydrogen bond distances
   
   //sort waters by oxygen z position
/*   if(pos[0].fz > 0.0){
	   watMinZ=-Rshel;
   }
   else{
	   watMinZ= pos[0].fz - Rshel;
   }
*/
   for(i=0;i<maxWat;i++){
      Wpos[i][0] = (float)(i*3);
      Wpos[i][1] = 888.0;
   }
   for(i=1;i<nw/3;i++)
   {
      Wpos[i][0]=(float)(i*3);
      Wpos[i][1]=-pos[i*3].fz;
//      if(pos[i*3].fz > watMinZ){
	  nwD++;
//      }
   }
   /*
   for(i=nw/3;i<maxWat;i++){
      Wpos[i][0] = (float)(i*3);
      Wpos[i][1] = 888.0;
   }
   */
   qSort(Wpos,0,maxWat-1);
   //find solute. this loop eliminates water above the ion.
   numAbove=0;
/*   while(pos[0].fz < pos[(int)(Wpos[numAbove][0])].fz){
	   numAbove++;
   }
*/
   // allocate distance matrix,graph,etc.
   ele = nwD - numAbove + 1; // '1' was nsolute
   //fprintf(stderr,"numAbove = %d, nwD = %d, ele = %d\n",numAbove, nwD, ele);
   if((dist  = calloc(ele, sizeof(double *))) == NULL ||
      (distList = calloc(ele*ele, sizeof(double))) == NULL ||
//      (distList2 = calloc(ele*ele, sizeof(double [2]))) == NULL || //jjk
      (graph = calloc(ele, sizeof(int *))) == NULL || 
      (queue = calloc(ele, sizeof(int))) == NULL ||
      (visited = calloc(ele, sizeof(int))) == NULL)
	  ERROR((stderr,"liqForce.c: out of core. wcor array initilization.\n"),exit);
   for(i=0;i<ele;i++){
   	if(( dist[i] = calloc(ele, sizeof(double))) == NULL ||
   	   (graph[i] = calloc(ele, sizeof(int))) == NULL)
	   ERROR((stderr,"liqForce.c: out of core. wcor array initilization.\n"),exit);
   }
   //whole loop is for getting rid of nearest water
//   for(i=0;i<ele*ele;i++){
//      if((distList2[i] = calloc(2,sizeof(double))) == NULL)
//	   ERROR((stderr,"liqForce.c: out of core. wcor array initilization.\n"),exit);
//   }
//fprintf(stderr,"initiatlized, ele = %d\n",ele);
   //populate dist matrix, initialize others
   fbulk=9999;
   for(i=0;i<ele;i++){
      if(i==ele-1){
         ivec.fx=pos[0].fx;
         ivec.fy=pos[0].fy;
         ivec.fz=pos[0].fz;
//	 ii=986;//jjk
      }
      else{
	 ivec.fx=pos[(int)Wpos[numAbove+i][0]].fx;
	 ivec.fy=pos[(int)Wpos[numAbove+i][0]].fy;
	 ivec.fz=pos[(int)Wpos[numAbove+i][0]].fz;
//	 ii=((int)(Wpos[numAbove+i][0]))/3;//jjk
	 if(pos[(int)Wpos[numAbove+i][0]].fz<=-5.0 && fbulk==9999){
		fbulk=i;
	 }
      }
      for(j=0;j<ele;j++){
         if(j==ele-1){
            jvec.fx=pos[0].fx;
            jvec.fy=pos[0].fy;
            jvec.fz=pos[0].fz;
//	    jj = 986;//jjk
	 }
         else{
	    jvec.fx=pos[(int)Wpos[numAbove+j][0]].fx;
	    jvec.fy=pos[(int)Wpos[numAbove+j][0]].fy;
	    jvec.fz=pos[(int)Wpos[numAbove+j][0]].fz;
//	    jj=((int)Wpos[numAbove+j][0])/3;
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
	 if(r2 < 0.1){
   	    distList[i*ele+j]=dist[i][j]=99.;//sqrt(r2);//99.;//got rid of self-compare
//	    distList2[i*ele+j][1]=99.;//jjk
	 }
/*	 if(i==ele-1 || j==ele-1){
   	    distList[i*ele+j]=dist[i][j]=sqrt(r2);//got rid of "-0.4"
	 }*/
	 else{
   	    distList[i*ele+j]=dist[i][j]=sqrt(r2);
//	    distList2[i*ele+j][1]=dist[i][j];//jjk
	 }
//	 distList2[i*ele+j][0]=(double)((int)(ii*ele+jj));//jjk
         graph[i][j]=0;
      }
   }
//fprintf(stderr,"fbulk = %d\n", fbulk);
   qSort2(distList,0,(ele*ele)-1); 
//   qSort(distList2,0,(ele*ele)-1); //jjk
//fprintf(stderr,"sorted\n");
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
//fprintf(stderr,"liqForce.c: %d is wcor position in distList... %d mod %d = %d.distList[irwf] =  %f ... find me r = %f, water1 = %d water2 = %d\n",irwf,irwf,ele,irwf%ele,distList[irwf],distList2[irwf][1],(int)(distList2[irwf][0])%ele, ((int)(distList2[irwf][0])/ele));       
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
	
int i,j,k,nw;
double cmposZ,ztop,zbot;
double RFIELD;
RFIELD = 0.0;
double totMass;

nw = natoms - 14*nNIT - nsolute;

ztop = zwall - 5; // reflecting wall 5A from top of box
zbot = -zwall+ 5; // reflecting wall 5A from bottom of box

for(i=0;i<nw;i=i+3){
   cmposZ = 0.0;
   for(j=0;j<3;j++){
      cmposZ += mass[i+j]*pos[i+j].fz;
   }
   cmposZ /= WMass;
   if(cmposZ > ztop){
      fprintf(stderr,"liqForce.c: water %d above. cmposZ = %f\n",i/3,cmposZ);
      RFIELD += ExCon*(cmposZ-ztop)*(cmposZ-ztop);
      for(j=0;j<3;j++){ //force > 0
	 force[i+j].fz -= 2*ExCon*mass[i+j]*(cmposZ-ztop)/WMass;
      }
      // compensating force added 11-July-2017 by john k
      totMass = 0.0;
      for(k=0;k<natoms;k++)
         totMass+=mass[k];
      totMass-=WMass;
      for(k=0;k<natoms;k++)
      {
	 if(k<i || k>=i+3)
	 {
	    force[k].fz += 2*ExCon*mass[k]*(cmposZ-ztop)/totMass;
	 }
      }
   }
   else if(cmposZ < zbot){
      fprintf(stderr,"liqForce.c: water %d below. cmposZ = %f\n",i/3,cmposZ);
      RFIELD += ExCon*(cmposZ-zbot)*(cmposZ-zbot);
      for(j=0;j<3;j++){ //force < 0
	 force[i+j].fz -= 2*ExCon*mass[i+j]*(cmposZ-zbot)/WMass;
      }
      // compensating force added 10-July-2017 by john k
      totMass = 0.0;
      for(k=0;k<natoms;k++)
         totMass+=mass[k];
      totMass-=WMass;
      for(k=0;k<natoms;k++)
      {
	 if(k<i || k>=i+3)
	 {
	    force[k].fz += 2*ExCon*mass[k]*(cmposZ-zbot)/totMass;
	 }
      }
   }
}
for(i=0;i< nNIT*14; i=i+14){
   cmposZ = 0.0;
   for(j=0;j<14;j++){
      cmposZ += mass[nw+i+j]*pos[nw+i+j].fz;
   }
   cmposZ /= NITMass;
   if(cmposZ > ztop){
      fprintf(stderr,"liqForce.c: nitrobenzene %d above. cmposZ = %f\n",i/14,cmposZ);
      RFIELD += ExCon*(cmposZ-ztop)*(cmposZ-ztop);
      for(j=0;j<14;j++){
	 force[nw+i+j].fz -= 2*ExCon*mass[nw+i+j]*(cmposZ-ztop)/NITMass;
      }
      // compensating force added 10-July-2017 by john k
      totMass = 0.0;
      for(k=0;k<natoms;k++)
         totMass+=mass[k];
      totMass-=NITMass;
      for(k=0;k<natoms;k++)
      {
	 if(k<nw+i || k>=nw+i+14)
	 {
	    force[k].fz += 2*ExCon*mass[k]*(cmposZ-ztop)/totMass;
	 }
      }
   }
   else if(cmposZ < zbot){
      fprintf(stderr,"liqForce.c: nitrobenzene %d below. cmposZ = %f\n",i/14,cmposZ);
      RFIELD += ExCon*(cmposZ-zbot)*(cmposZ-zbot);
      for(j=0;j<14;j++){
	 force[nw+i+j].fz -= 2*ExCon*mass[nw+i+j]*(cmposZ-zbot)/NITMass;
      }
      // compensating force added 10-July-2017 by john k
      totMass = 0.0;
      for(k=0;k<natoms;k++)
         totMass+=mass[k];
      totMass-=NITMass;
      for(k=0;k<natoms;k++)
      {
	 if(k<nw+i || k>=nw+i+14)
	 {
	    force[k].fz += 2*ExCon*mass[k]*(cmposZ-zbot)/totMass;
	 }
      }
   }
}
if(nsolute==1){
   cmposZ = pos[natoms-1].fz;
   if(cmposZ > ztop){
      fprintf(stderr,"liqForce.c: ion above. cmposZ = %f\n",cmposZ);
      RFIELD += ExCon*(cmposZ-ztop)*(cmposZ-ztop);
      force[natoms-1].fz -= 2*ExCon*mass[natoms-1]*(cmposZ-ztop)/ionMass;
   }
   else if(cmposZ < zbot){
      fprintf(stderr,"liqForce.c: ion below. cmposZ = %f\n",cmposZ);
      RFIELD += ExCon*(cmposZ-zbot)*(cmposZ-zbot);
      force[natoms-1].fz -= 2*ExCon*mass[natoms-1]*(cmposZ-zbot)/ionMass;
   }
}

return(RFIELD);
}


/* calculate the external forces and energy that remove surface roughness *
 * We assume that the water is in the z<0 region and the NIT  		 *
 * in the z > 0 region.							 */
extForce(nw,nnb)
	int nw;/* number of water atoms*/
	int nnb;/* number of nitrobenzene atoms*/
{
int i,j;
double cmposZ,waterM,nitM;
EXFIELD = 0.;
//
if(fixShel>0 && tc==0) return;//because fixShel waters not assigned yet
//
waterM = mass[0]+2.*mass[1];
nitM = 0.0;
for(i=0;i<14;i++){
	nitM+=mass[nw+i];
}
for(i=3;i<nw;i=i+3){ // first water is immune
    // ignore selected water molecules if 'fixShel' is active
    if(i/3 != fixw[0] && i/3 != fixw[1] && i/3 != fixw[2]){
	cmposZ = 0.;
	for (j=0;j<3;j++)
		cmposZ += mass[i+j]*pos[i+j].fz;
	if (cmposZ > 0) {/* water in the z>0 side, force it to the z<0 side */
		cmposZ /= waterM;
		EXFIELD += ExCon*cmposZ*cmposZ;
		kfwat += ExCon*cmposZ*cmposZ;
		for (j=0;j<3;j++)
			force[i+j].fz -= 2*ExCon*mass[i+j]*cmposZ/waterM;
	}
    }
}
for(i=0;i<nnb;i=i+14){
	cmposZ = 0.;
	for (j=0;j<14;j++)
		cmposZ += mass[nw+i+j]*pos[nw+i+j].fz;
	if (cmposZ < 0) {/* NIT in the z<0 side, force it to the z>0 side */
		cmposZ /= nitM;
		EXFIELD += ExCon*cmposZ*cmposZ;
		kfnit += ExCon*cmposZ*cmposZ;
		for (j=0;j<14;j++)
			force[nw+i+j].fz -= 2*ExCon*mass[nw+i+j]*cmposZ/nitM;
	}
}
}

// initialize matrices that need initialization
void initDataMatrices(int liqAlloc)
{
   if(liqAlloc){
      int i;
      for(i=0;i<10;i++){
         fixw[i]=-99; //initialize indices of  waters tethered to ion as unreal values
      }
      fixed = 0;
    }
}

//assign ids of water molecules to tether to the ion
void getTeth(int nw){
tripd image, sdist;
int i,j;
double rOI, r2;
double watr[5000][2];
for(i=0;i<5000;i++){
  watr[i][0]=(float)i;
  watr[i][1]=999.9;
}

for(i=3;i<nw;i+=3)/* use nw=number water molecules *3 */{
  image.fx = -(sdist.fx = pos[i].fx - pos[0].fx);
  image.fy = -(sdist.fy = pos[i].fy - pos[0].fy);
  image.fz = -(sdist.fz = pos[i].fz - pos[0].fz);
  mvimage(&sdist);
  image.fx += sdist.fx;
  image.fy += sdist.fy;
  image.fz += sdist.fz;
  r2 = sdist.fx*sdist.fx +  sdist.fy*sdist.fy +  sdist.fz*sdist.fz;
  rOI = sqrt(r2);
  watr[i][0]=(float)(i/3);
  watr[i][1]=rOI;
}

qSort(watr,0,4999);

// here we can override if specially prepared input files
//if(fixShel == 3 && KillFinger == 0){
//	watr[0][0] = 1.;
//	watr[1][0] = 2.;
//	watr[2][0] = 3.;
//} 
//watr[0][0] = 502.;
//watr[1][0] = 415.;
//watr[2][0] = 520.;

for(i=0;i<fixShel;i++){
  fixw[i]=(int)(watr[i][0]);
  fixed++;
  fprintf(stderr,"liqForce.c - tether water %d to ion (r=%4.3f). (%d/%d)\n",fixw[i],fixed,watr[i][1],fixShel);
}

}
