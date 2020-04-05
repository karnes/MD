#include	<md.h>
#include	<system.h>
#include	 <math.h>

BrForce()
{
int	i, j, k, l, m;
double	totMass, r, r2, dedr, delfx, delfy, delfz;
double	eg, wnc, s, sp;
tripd	com, frc, image, sdl, f[THAsites+1]; //f[GLYsites+THAsites];
double	bgljq();
double	btljq();

//return(0);

for(i=0;i<nGLY;i++){
      k = i*GLYsites;
      l = nGLY*GLYsites + nTHA*THAsites; //index of Br- ion
      Xz = pos[l].fz;
      eg = wnc = 0.0;

/*****	Determine image vector, GLY parent - Br-.  ******/
      image.fx = -(sdl.fx = pos[k].fx - pos[l].fx);
      image.fy = -(sdl.fy = pos[k].fy - pos[l].fy);
      image.fz = -(sdl.fz = pos[k].fz - pos[l].fz);
      mvimage(&sdl);
      image.fx += (delfx =sdl.fx);
      image.fy += (delfy =sdl.fy);
      image.fz += (delfz =sdl.fz);
      r2 = delfx*delfx + delfy*delfy + delfz*delfz;
      if(r2 >= swr2max)
         continue;
      if(r2 <= swr2min)
      {
         s = 1.0;
	 sp = 0.0;
      }
      else
	 swtch(r2-swr2min,&s,&sp);
      for(m=0;m<(THAsites+1); m++) //clear full array (same for GLY/Br and THA/Br) 
         f[m].fx = f[m].fy = f[m].fz = 0.;
      for(m=0;m<GLYsites;m++){ /*Loop over atoms in GLY molecules*/
/*****	Determine image vector ******/
	 delfx = pos[k+m].fx - pos[l].fx + image.fx;
	 delfy = pos[k+m].fy - pos[l].fy + image.fy;
	 delfz = pos[k+m].fz - pos[l].fz + image.fz;
	 r2 = delfx*delfx + delfy*delfy + delfz*delfz;
/*****	Get (1/r dV/dr) ******/
	 dedr = bgljq(r2,m,GLYsites,&eg,k,&wnc,l);
/*****	Resolve forces on atoms.  ******/
	 f[GLYsites].fx += (delfx *= dedr);
	 f[GLYsites].fy += (delfy *= dedr);
	 f[GLYsites].fz += (delfz *= dedr);
	 f[m].fx -= delfx;
	 f[m].fy -= delfy;
	 f[m].fz -= delfz;
      }
      if(sp != 0.0){
	 for(m=0;m<GLYsites+1;m++){
	       f[m].fx *= s;
	       f[m].fy *= s;
	       f[m].fz *= s;
	 }
	 f[0].fx -= (frc.fx = sp*eg*sdl.fx);
	 f[0].fy -= (frc.fy = sp*eg*sdl.fy);
	 f[0].fz -= (frc.fz = sp*eg*sdl.fz);
	 f[GLYsites].fx += frc.fx;
	 f[GLYsites].fy += frc.fy;
	 f[GLYsites].fz += frc.fz;
	 eg *= s;
	 wnc *= s;
      }
      k=i*GLYsites;
      l = nGLY*GLYsites+nTHA*THAsites;
	/* GLY */
	force[k].fx += f[0].fx;
	force[k].fy += f[0].fy;
	force[k].fz += f[0].fz;
	force[++k].fx += f[1].fx;
	force[k].fy += f[1].fy;
	force[k].fz += f[1].fz;
	force[++k].fx += f[2].fx;
	force[k].fy += f[2].fy;
	force[k].fz += f[2].fz;
	force[++k].fx += f[3].fx;
	force[k].fy += f[3].fy;
	force[k].fz += f[3].fz;
	force[++k].fx += f[4].fx;
	force[k].fy += f[4].fy;
	force[k].fz += f[4].fz;
	force[++k].fx += f[5].fx;
	force[k].fy += f[5].fy;
	force[k].fz += f[5].fz;
	force[++k].fx += f[6].fx;
	force[k].fy += f[6].fy;
	force[k].fz += f[6].fz;
	force[++k].fx += f[7].fx;
	force[k].fy += f[7].fy;
	force[k].fz += f[7].fz;
	force[++k].fx += f[8].fx;
	force[k].fy += f[8].fy;
	force[k].fz += f[8].fz;
	force[++k].fx += f[9].fx;
	force[k].fy += f[9].fy;
	force[k].fz += f[9].fz;
	force[++k].fx += f[10].fx;
	force[k].fy += f[10].fy;
	force[k].fz += f[10].fz;
	force[++k].fx += f[11].fx;
	force[k].fy += f[11].fy;
	force[k].fz += f[11].fz;
	force[++k].fx += f[12].fx;
	force[k].fy += f[12].fy;
	force[k].fz += f[12].fz;
	force[++k].fx += f[13].fx;
	force[k].fy += f[13].fy;
	force[k].fz += f[13].fz;

	/* Br */
	force[l].fx += f[GLYsites].fx;
	force[l].fy += f[GLYsites].fy;
	force[l].fz += f[GLYsites].fz;
	/*   
	*/
	XGLYV += eg;
	XGLYC += wnc;
}

if(nTHA==1){
  // find cos of angle formed by +z and Br-->THAcom
  com.fx = com.fy = com.fz = totMass = 0.0;
  l = nGLY*GLYsites + nTHA*THAsites; //index of Br- ion
  for(i=0;i<THAsites;i++){
     k = nGLY*GLYsites+i;
     totMass+=mass[k];
     com.fx+=pos[k].fx*mass[k];
     com.fy+=pos[k].fy*mass[k];
     com.fz+=pos[k].fz*mass[k];
  }
  com.fx/=totMass;
  com.fy/=totMass;
  com.fz/=totMass;
  
  sdl.fx = com.fx - pos[l].fx;
  sdl.fy = com.fy - pos[l].fy;
  sdl.fz = com.fz - pos[l].fz;
  mvimage(&sdl);
  r2 = sdl.fx*sdl.fx+sdl.fy*sdl.fy+sdl.fz*sdl.fz;
  r=sqrt(r2);
  cosTHAX = sdl.fz/r;

  for(i=0;i<THAsites;i++){
      k = nGLY*GLYsites+i;
      l = nGLY*GLYsites + nTHA*THAsites; //index of Br- ion

      eg = wnc = 0.0;

/*****	Determine vector, THA atom - Br-.  ******/
      sdl.fx = pos[k].fx - pos[l].fx;
      sdl.fy = pos[k].fy - pos[l].fy;
      sdl.fz = pos[k].fz - pos[l].fz;
      mvimage(&sdl);
      delfx =sdl.fx;
      delfy =sdl.fy;
      delfz =sdl.fz;
      r2 = delfx*delfx + delfy*delfy + delfz*delfz;
      if(r2 >= swr2max)
         continue;
      if(r2 <= swr2min)
      {
         s = 1.0;
	 sp = 0.0;
      }
      else
	 swtch(r2-swr2min,&s,&sp);
      for(m=0;m<(THAsites+1); m++){ //clear full array (same for GLY/Br and THA/Br) 
         f[m].fx = f[m].fy = f[m].fz = 0.;
      }
//      for(m=0;m<GLYsites;m++){ /*Loop over atoms in GLY molecules*/
/*****	Determine image vector ******/
//      delfx = pos[k].fx - pos[l].fx + image.fx;
//      delfy = pos[k].fy - pos[l].fy + image.fy;
//      delfz = pos[k].fz - pos[l].fz + image.fz;
//      r2 = delfx*delfx + delfy*delfy + delfz*delfz;
/*****	Get (1/r dV/dr) ******/
      dedr = btljq(r2,i,THAsites,&eg,k,&wnc,l);
/*****	Resolve forces on atoms.  ******/
      f[1].fx += (delfx *= dedr);//Br-
      f[1].fy += (delfy *= dedr);
      f[1].fz += (delfz *= dedr);
      f[0].fx -= delfx;//THA atom
      f[0].fy -= delfy;
      f[0].fz -= delfz;
//    }
      if(sp != 0.0){
	 for(m=0;m<2;m++){
	       f[m].fx *= s;
	       f[m].fy *= s;
	       f[m].fz *= s;
	 }
	 f[1].fx -= (frc.fx = sp*eg*sdl.fx);
	 f[1].fy -= (frc.fy = sp*eg*sdl.fy);
	 f[1].fz -= (frc.fz = sp*eg*sdl.fz);
	 f[0].fx += frc.fx;
	 f[0].fy += frc.fy;
	 f[0].fz += frc.fz;
	 eg *= s;
	 wnc *= s;
      }
      k = nGLY*GLYsites+i; //THA atom
      l = nGLY*GLYsites+nTHA*THAsites;//Br- ion
	/* THA atom */
	force[k].fx += f[0].fx;
	force[k].fy += f[0].fy;
	force[k].fz += f[0].fz;

	/* Br */
	force[l].fx += f[1].fx;
	force[l].fy += f[1].fy;
	force[l].fz += f[1].fz;
	/*   
	*/
	XTHAV += eg;
	XTHAC += wnc;
//	fprintf(stderr,"%f\n",KCAL*(eg));
   }

}

}

double
bgljq(r2,m,n,eg,nN,wnc,l)
double r2, *eg,*wnc;
int m, /* m is the GLY atom index */
    n, /* n is the Br- (GLYsites) atom index */
    nN,/* the number of the nitrobenzene molecule for the rdf calculations*/
    l; /* index of water O*/
{
double der,r,ir,ir6, ec,a,b,q,aa,bb;
int bin;
double en1;
int index;
	ir = 1. / r2;
	ir6 = ir * ir * ir;
	a = GLYlj[m][n].a;
	b = GLYlj[m][n].b;
	q = GLYlj[m][n].q;
	aa =12*a;
	bb = 6*b;
	r = sqrt(r2);
	ec = q/r;
	*wnc += ec;
	en1 =  ( a * ir6 - b ) * ir6 + ec;
	*eg +=  ( a * ir6 - b ) * ir6 + ec;
	der = (bb - aa * ir6) * ir6 * ir - ec/r2 ;
        if(tc%grdt==0){
	   index = (int)(r/binRDF);
	   if(index<300){// g(r)'s
	      if(m==0||m==4||m==9){//C-Br
	         grGLY[7][index]+=nGLY*1.0/3.0;
		 if(r<gGCXmax[0]){
		    gGCX[0]++;
		    XGLYs[0]+=en1;
		 }
	      }
	      else if(m==1||m==5||m==10){//O-Br
	         grGLY[8][index]+=nGLY*1.0/3.0;
		 if(r<gGOXmax[0]){
		    gGOX[0]++;
		    XGLYs[1]+=en1;
		 }
	      }
	      else if(m==2||m==6||m==11){//H(oh)-Br
	         grGLY[9][index]+=nGLY*1.0/3.0;
		 if(r<gGHXmax[0]){
		    gGHX[0]++;
		    XGLYs[2]+=en1;
		 }
	      }
	      else if(r<gGCXmax[0]){
		 XGLYs[3]+=en1;
	      }
	   }
        }
	return(der);
}
double
btljq(r2,m,n,eg,nN,wnc,l)
double r2, *eg,*wnc;
int m, /* m is the THA atom index */
    n, /* n is the Br- (THAsites) index */
    nN,/* the number of the nitrobenzene molecule for the rdf calculations*/
    l; /* index of water O*/
{
double der,r,ir,ir6, ec,a,b,q,aa,bb;
int bin;
int rindex;
	ir = 1. / r2;
	ir6 = ir * ir * ir;
	a = THAlj[m][n].a;
	b = THAlj[m][n].b;
	q = THAlj[m][n].q;
	aa =12*a;
	bb = 6*b;
	r = sqrt(r2);
        if(m==0)
	   rXN = r;
	ec = q/r;
	*wnc += ec;
	*eg +=  ( a * ir6 - b ) * ir6 + ec;
	der = (bb - aa * ir6) * ir6 * ir - ec/r2 ;
//	fprintf(stderr,"%f\n",KCAL*(*eg));
        return(der);
}
