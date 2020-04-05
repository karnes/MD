#include	<md.h>
#include	<system.h>
#include	 <math.h>

TSForce()
{
int	i, j, k, l,m;
double	totMass, r, r2, dedr, delfx, delfy, delfz;
double	eg, wnc, s, sp, deriv;
tripd	d,com, frc, image, sdl, f[THAsites+3]; //f[GLYsites+THAsites];
double  gTS();
double  tTS();

TSbondE = TSbendE = 0.0;

if(nTS!=1){
   return(0);
}

intraTS();
//get z, cos vs. z of TS
l = nGLY*GLYsites+nTHA*THAsites; // center Cl in TS
Xz = (pos[l].fz*mass[l]+pos[l+1].fz*mass[l+1]+pos[l+2].fz*mass[l+2])/(mass[l]+mass[l+1]+mass[l+2]);
r = sqrt(sq(pos[l+1].fx-pos[l+2].fx)+sq(pos[l+1].fy-pos[l+2].fy)+sq(pos[l+1].fz-pos[l+2].fz));
cosXz = (pos[l+1].fz-pos[l+2].fz)/r;

for(i=0;i<nGLY;i++){
      k = i*GLYsites;
      l = nGLY*GLYsites+nTHA*THAsites; // center Cl in TS
// *****	Determine image vector, GLY parent - Cl center of TS.  ******
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
   //GLY - TS interactions
   for(m=0;m<GLYsites;m++)
      gflag[m] = 0;
   for(m=0;m<(THAsites+3); m++) //clear full array (same for GLY/Br and THA/Br) 
      f[m].fx = f[m].fy = f[m].fz = 0.;
   eg = wnc = 0.0;
   for(j=0;j<3;j++){
      l = nGLY*GLYsites + nTHA*THAsites + j; //index of TS atom
      k = i*GLYsites;
      for(m=0;m<GLYsites;m++){ //Loop over atoms in GLY molecules
// *****	Determine image vector ******
	 delfx = pos[k+m].fx - pos[l].fx + image.fx;
	 delfy = pos[k+m].fy - pos[l].fy + image.fy;
	 delfz = pos[k+m].fz - pos[l].fz + image.fz;
	 r2 = delfx*delfx + delfy*delfy + delfz*delfz;
// *****	Get (1/r dV/dr) ******
	 dedr = gTS(r2,m,3+j,&eg,&wnc);
// *****	Resolve forces on atoms.  ******
	 f[GLYsites+j].fx += (delfx *= dedr);
	 f[GLYsites+j].fy += (delfy *= dedr);
	 f[GLYsites+j].fz += (delfz *= dedr);
	 f[m].fx -= delfx;
	 f[m].fy -= delfy;
	 f[m].fz -= delfz;
      }
   }
   if(sp != 0.0){
	 for(m=0;m<GLYsites+3;m++){
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
	// * GLY *
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
        k = i* GLYsites;
	// TS (LJq) 
	force[l].fx += f[GLYsites].fx;
	force[l].fy += f[GLYsites].fy;
	force[l].fz += f[GLYsites].fz;
	force[++l].fx += f[GLYsites+1].fx;
	force[l].fy += f[GLYsites+1].fy;
	force[l].fz += f[GLYsites+1].fz;
	force[++l].fx += f[GLYsites+2].fx;
	force[l].fy += f[GLYsites+2].fy;
	force[l].fz += f[GLYsites+2].fz;
	//   
	//
	XGLYV += eg;
	XGLYC += wnc;
}


if(nTHA==1){

   for(i=0;i<THAsites;i++){
      k = nGLY*GLYsites+i;
      l = nGLY*GLYsites+nTHA*THAsites; // center Cl in TS
// *****	Determine image vector, GLY parent - Cl center of TS.  ******
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
   //GLY - TS interactions
   for(m=0;m<(THAsites+3); m++) //clear full array (same for GLY/Br and THA/Br) 
      f[m].fx = f[m].fy = f[m].fz = 0.;
   eg = wnc = 0.0;
   for(j=0;j<3;j++){
      l = nGLY*GLYsites + nTHA*THAsites + j; //index of TS atom
      k = nGLY*GLYsites + i;
// *****	Determine image vector ******
	 delfx = pos[k].fx - pos[l].fx + image.fx;
	 delfy = pos[k].fy - pos[l].fy + image.fy;
	 delfz = pos[k].fz - pos[l].fz + image.fz;
	 r2 = delfx*delfx + delfy*delfy + delfz*delfz;
// *****	Get (1/r dV/dr) ******
	 dedr = tTS(r2,i,3+j,&eg,&wnc);
// *****	Resolve forces on atoms.  ******
	 f[1+j].fx += (delfx *= dedr);//TS atom center
	 f[1+j].fy += (delfy *= dedr);
	 f[1+j].fz += (delfz *= dedr);
	 f[0].fx -= delfx;//THA atom
	 f[0].fy -= delfy;
	 f[0].fz -= delfz;
   }
   if(sp != 0.0){
	 for(m=0;m<4;m++){
	       f[m].fx *= s;
	       f[m].fy *= s;
	       f[m].fz *= s;
	 }
	 f[0].fx -= (frc.fx = sp*eg*sdl.fx);//THA
	 f[0].fy -= (frc.fy = sp*eg*sdl.fy);
	 f[0].fz -= (frc.fz = sp*eg*sdl.fz);
	 f[1].fx += frc.fx;//TS center
	 f[1].fy += frc.fy;
	 f[1].fz += frc.fz;
	 eg *= s;
	 wnc *= s;
   }
   k=nGLY*GLYsites+i;
   l = nGLY*GLYsites+nTHA*THAsites;
	// * THA  *
	force[k].fx += f[0].fx;
	force[k].fy += f[0].fy;
	force[k].fz += f[0].fz;
	// TS (LJq) 
	force[l].fx += f[1].fx;
	force[l].fy += f[1].fy;
	force[l].fz += f[1].fz;
	force[++l].fx += f[2].fx;
	force[l].fy += f[2].fy;
	force[l].fz += f[2].fz;
	force[++l].fx += f[3].fx;
	force[l].fy += f[3].fy;
	force[l].fz += f[3].fz;
	//   
	XTHAV += eg;
	XTHAC += wnc;
   }

}

}

double
gTS(r2,m,n,eg,wnc)
double r2,*eg,*wnc;
int m, /* m is the THA atom index */
    n; /* n is the TS atom index (3-5) */
{
double der,r,ir,ir6, ec,a,b,q,aa,bb;
double en1;
int index;
	ir = 1. / r2;
	ir6 = ir * ir * ir;
	a = solGLYlj[m][n].a;
	b = solGLYlj[m][n].b;
	q = solGLYlj[m][n].q;
	aa =12*a;
	bb = 6*b;
	r = sqrt(r2);
	ec = q/r;
	*wnc += ec;
	en1 =  ( a * ir6 - b ) * ir6 + ec;
	*eg +=  ( a * ir6 - b ) * ir6 + ec;
        if(tc%grdt==0){
	   index = (int)(r/binRDF);
	   if(index<300){// g(r)'s
	      if(m==0||m==4||m==9){//C-Br
	         grGLY[(n-3)*3+7][index]+=nGLY*1.0/3.0;
		 if(gflag[m]==0 && r<gGCXmax[n-3]){
		    gGCX[n-3]++;
		    gflag[m] = 1;
		    XGLYs[0] += en1;
		 }
	      }
	      else if(m==1||m==5||m==10){//O-Br
	         grGLY[(n-3)*3+8][index]+=nGLY*1.0/3.0;
		 if(gflag[m]==0 && r<gGOXmax[n-3]){
		    gGOX[n-3]++;
		    gflag[m] = 1;
		    XGLYs[1] += en1;
		 }
	      }
	      else if(m==2||m==6||m==11){//H(oh)-Br
	         grGLY[(n-3)*3+9][index]+=nGLY*1.0/3.0;
		 if(gflag[m]==0 && r<gGHXmax[n-3]){
		    gGHX[n-3]++;
		    gflag[m] = 1;
		    XGLYs[2] += en1;
		 }
	      }
	      else if(gflag[m]==0 && r<gGCXmax[n-3]){
		 XGLYs[3] += en1;
		 gflag[m] = 1;
	      }
	   }
        }
	der = (bb - aa * ir6) * ir6 * ir - ec/r2 ;
//	fprintf(stderr,"%f\n",KCAL*(*eg));
        return(der);
}

int intraTS()
{
   int i,j,k,i0,i1,i2;
   double r[3],q1,q2,dist,deriv,pe,dedr[3];
   tripd grad[3],d,r0,r1;
   double r1r2,n1n2,da,dedd;

//Cl(j)-Cl(i)-Br(k)
i=nGLY*GLYsites+nTHA*THAsites;
j=i+1;
k=i+2;
//
d.fx = pos[i].fx - pos[j].fx;
d.fy = pos[i].fy - pos[j].fy;
d.fz = pos[i].fz - pos[j].fz;
dist = sqrt(d.fx*d.fx+d.fy*d.fy+d.fz*d.fz);
deriv = TSkstr[0]*(dist-TSreq[0])/dist;
pe = 0.5*TSkstr[0]*sq(dist-TSreq[0]);
force[i].fx -= (d.fx *= deriv);
force[i].fy -= (d.fy *= deriv);
force[i].fz -= (d.fz *= deriv);
force[j].fx += d.fx;
force[j].fy += d.fy;
force[j].fz += d.fz;
TSbondE += pe;
//if(tc%20==0)
//   fprintf(stderr,"bond 1: %f,",dist);
//
d.fx = pos[i].fx - pos[k].fx;
d.fy = pos[i].fy - pos[k].fy;
d.fz = pos[i].fz - pos[k].fz;
dist = sqrt(d.fx*d.fx+d.fy*d.fy+d.fz*d.fz);
deriv = TSkstr[1]*(dist-TSreq[1])/dist;
pe = 0.5*TSkstr[1]*sq(dist-TSreq[1]);
force[i].fx -= (d.fx *= deriv);
force[i].fy -= (d.fy *= deriv);
force[i].fz -= (d.fz *= deriv);
force[k].fx += d.fx;
force[k].fy += d.fy;
force[k].fz += d.fz;
TSbondE += pe;
//if(tc%20==0)
//   fprintf(stderr," bond 2: %f\n",dist);
//fprintf(stderr,"check bond eq's: 1 = %f, 2 = %f\n",TSreq[0],TSreq[1]);
// calcbend

//Cl(1)-Cl(0)-Br(2)
i0 = i; //center atom
i1 = j;
i2 = k;

r0.fx = pos[i1].fx - pos[i0].fx;
r0.fy = pos[i1].fy - pos[i0].fy;
r0.fz = pos[i1].fz - pos[i0].fz;
r1.fx = pos[i2].fx - pos[i0].fx;
r1.fy = pos[i2].fy - pos[i0].fy;
r1.fz = pos[i2].fz - pos[i0].fz;

r[0] = sqrt(r0.fx*r0.fx + r0.fy*r0.fy + r0.fz*r0.fz);
r[1] = sqrt(r1.fx*r1.fx + r1.fy*r1.fy + r1.fz*r1.fz);
r[2] = r0.fx*r1.fx+r0.fy*r1.fy+r0.fz*r1.fz;

n1n2 = r[2] /(r[0]*r[1]);

// *	for angles very close to 180 degrees, n1n2 = -1 + eps,
// *	where eps is very small and can be both positive and
// *	negative. If its negative acos will blow up. If its
// *	positive but very small, da will be very small for a linear
// *	equilibrium state and dedr will be the ratio between two
// *	very small quantities. We therefore use a switch	
if(n1n2+1 < 1.0e-8) 
   n1n2 = -1;
da = acos(n1n2) - TSeqbend;
q1 = r[0]-TSreq[0];//eqCN;
q2 = r[1]-TSreq[1];//eqCC;
TSbendE += 0.5*(TSkstr[0]*q1*q1+TSkstr[1]*q2*q2+TSkbend*da*da);
dedr[0] = TSkstr[0]*q1/r[0];
dedr[1] = TSkstr[1]*q2/r[1];
if(n1n2+1 < 1.0e-8)
   dedr[2] = -TSkbend/(r[0]*r[1]);
else
   dedr[2] = TSkbend*da/ (sqrt(1.-n1n2*n1n2)*r[0]*r[1]);
force[i1].fx -= dedr[0]*r0.fx;
force[i1].fy -= dedr[0]*r0.fy;
force[i1].fz -= dedr[0]*r0.fz;
force[i2].fx -= dedr[1]*r1.fx;
force[i2].fy -= dedr[1]*r1.fy;
force[i2].fz -= dedr[1]*r1.fz;
force[i0].fx += dedr[0]*r0.fx + dedr[1]*r1.fx;
force[i0].fy += dedr[0]*r0.fy + dedr[1]*r1.fy;
force[i0].fz += dedr[0]*r0.fz + dedr[1]*r1.fz;
force[i1].fx += (grad[0].fx = dedr[2]*(r1.fx-r[2]*r0.fx/(r[0]*r[0])));
force[i1].fy += (grad[0].fy = dedr[2]*(r1.fy-r[2]*r0.fy/(r[0]*r[0])));
force[i1].fz += (grad[0].fz = dedr[2]*(r1.fz-r[2]*r0.fz/(r[0]*r[0])));
force[i2].fx += (grad[1].fx = dedr[2]*(r0.fx-r[2]*r1.fx/(r[1]*r[1])));
force[i2].fy += (grad[1].fy = dedr[2]*(r0.fy-r[2]*r1.fy/(r[1]*r[1])));
force[i2].fz += (grad[1].fz = dedr[2]*(r0.fz-r[2]*r1.fz/(r[1]*r[1])));
force[i0].fx -= grad[1].fx + grad[0].fx;
force[i0].fy -= grad[1].fy + grad[0].fy;
force[i0].fz -= grad[1].fz + grad[0].fz;
}

double
tTS(r2,m,n,eg,wnc)
double r2,*eg,*wnc;
int m, /* m is the THA atom index */
    n; /* n is the Cl (THAsites+1) index */
{
double der,r,ir,ir6, ec,a,b,q,aa,bb;
int bin;
int rindex;
	ir = 1. / r2;
	ir6 = ir * ir * ir;
	a = solTHAlj[m][n].a;
	b = solTHAlj[m][n].b;
	q = solTHAlj[m][n].q;
	aa =12*a;
	bb = 6*b;
	r = sqrt(r2);
	ec = q/r;
	*wnc += ec;
	*eg +=  ( a * ir6 - b ) * ir6 + ec;
	der = (bb - aa * ir6) * ir6 * ir - ec/r2 ;
//	fprintf(stderr,"%f\n",KCAL*(*eg));
        return(der);
}
