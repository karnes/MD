#include	<md.h>
#include	<system.h>
#include	 <math.h>
#define ldist 0.8
Cl2Force()
{

int	i, j, k, l,C1,C2,m;
double	totMass, r, r2, dedr, delfx, delfy, delfz;
double	eg, wnc, s, sp, deriv, r_12, ir12;
tripd	d,com, frc, image, sdl, f[3]; 
tripd	r_a,r_b,r_c,r_1,r_2,r1mr2,r2mr1; 
double	gCl(),tCl();
double  gM0(),tM0(),gMp();
double sqrt();

/////////////////
// ATOM INDICES
// Br = 0
// Cl(1) = 1
// Cl(2) = 2
////////////////

//M0e = Me = Cle = 0.0;
Cl2BondE=0.0;

if(nCl2!=1){
   return(0);
}

// intra Cl2
i=1;
j=i+1;
//r1mr2 is vector Cl(2)-->Cl(1)
r2mr1.fx = -(r1mr2.fx = pos[i].fx - pos[j].fx);
r2mr1.fy = -(r1mr2.fy = pos[i].fy - pos[j].fy);
r2mr1.fz = -(r1mr2.fz = pos[i].fz - pos[j].fz);
r_12 = sqrt(r1mr2.fx*r1mr2.fx+r1mr2.fy*r1mr2.fy+r1mr2.fz*r1mr2.fz);
ir12 = ldist/(r_12*r_12*r_12);

//stretch is the only intra Cl2 force
deriv = Clkstr*(r_12-Clreq)/r_12;
eg = 0.5*Clkstr*sq(r_12-Clreq);
force[i].fx -= (d.fx = deriv*r1mr2.fx);
force[i].fy -= (d.fy = deriv*r1mr2.fy);
force[i].fz -= (d.fz = deriv*r1mr2.fz);
force[j].fx += d.fx;
force[j].fy += d.fy;
force[j].fz += d.fz;
Cl2BondE += eg;

//
// a--1--b--2--c
//
//Calculate r_a, r_b, r_c
r_b.fx = (pos[i].fx + pos[j].fx)/2.0;
r_b.fy = (pos[i].fy + pos[j].fy)/2.0;
r_b.fz = (pos[i].fz + pos[j].fz)/2.0;

r_a.fx = pos[i].fx + ldist*(r1mr2.fx / r_12);
r_a.fy = pos[i].fy + ldist*(r1mr2.fy / r_12);
r_a.fz = pos[i].fz + ldist*(r1mr2.fz / r_12);
r_c.fx = pos[j].fx + ldist*(r2mr1.fx / r_12);
r_c.fy = pos[j].fy + ldist*(r2mr1.fy / r_12);
r_c.fz = pos[j].fz + ldist*(r2mr1.fz / r_12);

//test 5-site definition
//if(tc==0){
//   fprintf(stderr,"5\n");
//   fprintf(stderr,"M Cl M0 Cl' M'\n");
//   fprintf(stderr,"H %f %f %f\n",M[1].fx,M[1].fy,M[1].fz);
//   fprintf(stderr,"C %f %f %f\n",pos[i].fx,pos[i].fy,pos[i].fz);
//   fprintf(stderr,"H %f %f %f\n",M[0].fx,M[0].fy,M[0].fz);
//   fprintf(stderr,"C %f %f %f\n",pos[i+1].fx,pos[i+1].fy,pos[i+1].fz);
//   fprintf(stderr,"H %f %f %f\n",M[2].fx,M[2].fy,M[2].fz);
//}
//
k = 0; 
// *****	Determine image vector, Br atom - b (center of Cl2).  ******
image.fx = -(sdl.fx = r_b.fx - pos[k].fx);
image.fy = -(sdl.fy = r_b.fy - pos[k].fy);
image.fz = -(sdl.fz = r_b.fz - pos[k].fz);
mvimage(&sdl);
image.fx += (delfx =sdl.fx);
image.fy += (delfy =sdl.fy);
image.fz += (delfz =sdl.fz);
r2 = delfx*delfx + delfy*delfy + delfz*delfz;
XXdist = sqrt(r2);
//if(r2 >= swr2max)
//   continue;
if(r2 <= swr2min)
{
   s = 1.0;
   sp = 0.0;
}
else
   swtch(r2-swr2min,&s,&sp);
////
//// BEGIN PROBLEMATIC PART
////


// Br - a interactions -- project forces onto Cl atoms

// Interaction between Br and Cl virtual site 'a'

eg = wnc = 0.0;
for(m=0;m<3;m++){ //clear full array 
   f[m].fx = f[m].fy = f[m].fz = 0.;
}
// *****	Determine image vector ******
delfx = r_a.fx - pos[k].fx + image.fx;
delfy = r_a.fy - pos[k].fy + image.fy;
delfz = r_a.fz - pos[k].fz + image.fz;
r2 = delfx*delfx + delfy*delfy + delfz*delfz;
// *****   Get (1/r dV/dr) ******
dedr = gMp(r2,0,&eg);
// *****   Resolve forces on atoms.  ******
f[1].fx -= delfx*dedr*(1.0 + ldist/r_12 - ir12*r1mr2.fx*r1mr2.fx);
f[1].fy -= delfy*dedr*(1.0 + ldist/r_12 - ir12*r1mr2.fy*r1mr2.fy);
f[1].fz -= delfz*dedr*(1.0 + ldist/r_12 - ir12*r1mr2.fz*r1mr2.fz);

f[2].fx -= delfx*dedr*( -ldist/r_12 + ir12*r1mr2.fx*r1mr2.fx);
f[2].fy -= delfy*dedr*( -ldist/r_12 + ir12*r1mr2.fy*r1mr2.fy);
f[2].fz -= delfz*dedr*( -ldist/r_12 + ir12*r1mr2.fz*r1mr2.fz);

f[0].fx += delfx*dedr;
f[0].fy += delfy*dedr;
f[0].fz += delfz*dedr;

// cutoff stuff here (deleted for debugging)

force[0].fx += f[0].fx;
force[0].fy += f[0].fy;
force[0].fz += f[0].fz;

force[1].fx += f[1].fx;
force[1].fy += f[1].fy;
force[1].fz += f[1].fz;
force[2].fx += f[2].fx;
force[2].fy += f[2].fy;
force[2].fz += f[2].fz;

INTER_X += eg;

// Interaction between Br and Cl virtual site 'c'

eg = wnc = 0.0;
for(m=0;m<3;m++){ //clear full array 
   f[m].fx = f[m].fy = f[m].fz = 0.;
}
// *****	Determine image vector ******
delfx = r_c.fx - pos[k].fx + image.fx;
delfy = r_c.fy - pos[k].fy + image.fy;
delfz = r_c.fz - pos[k].fz + image.fz;
r2 = delfx*delfx + delfy*delfy + delfz*delfz;
// *****   Get (1/r dV/dr) ******
dedr = gMp(r2,0,&eg);
// *****   Resolve forces on atoms.  ******
f[2].fx -= delfx*dedr*(1.0 + ldist/r_12 - ir12*r1mr2.fx*r1mr2.fx);
f[2].fy -= delfy*dedr*(1.0 + ldist/r_12 - ir12*r1mr2.fy*r1mr2.fy);
f[2].fz -= delfz*dedr*(1.0 + ldist/r_12 - ir12*r1mr2.fz*r1mr2.fz);

f[1].fx -= delfx*dedr*( -ldist/r_12 + ir12*r1mr2.fx*r1mr2.fx);
f[1].fy -= delfy*dedr*( -ldist/r_12 + ir12*r1mr2.fy*r1mr2.fy);
f[1].fz -= delfz*dedr*( -ldist/r_12 + ir12*r1mr2.fz*r1mr2.fz);

f[0].fx += delfx*dedr;
f[0].fy += delfy*dedr;
f[0].fz += delfz*dedr;

// cutoff stuff here (deleted for debugging)

force[0].fx += f[0].fx;
force[0].fy += f[0].fy;
force[0].fz += f[0].fz;

force[2].fx += f[2].fx;
force[2].fy += f[2].fy;
force[2].fz += f[2].fz;
force[1].fx += f[1].fx;
force[1].fy += f[1].fy;
force[1].fz += f[1].fz;

INTER_X += eg;

////
//// END PROBLEMATIC PART
////

// b - GLY interactions
eg = wnc = 0.0;
for(m=0;m<3; m++) //clear full array (same for GLY/Br and THA/Br) 
   f[m].fx = f[m].fy = f[m].fz = 0.;

// ***** Determine image vector ******
delfx = r_b.fx - pos[0].fx + image.fx;
delfy = r_b.fy - pos[0].fy + image.fy;
delfz = r_b.fz - pos[0].fz + image.fz;
r2 = delfx*delfx + delfy*delfy + delfz*delfz;
// *****	Get (1/r dV/dr) ******
dedr = gM0(r2,1,&eg);
// *****	Resolve forces on atoms.  ******
f[1].fx -= (delfx *= dedr);
f[1].fy -= (delfy *= dedr);
f[1].fz -= (delfz *= dedr);
f[0].fx += delfx;
f[0].fy += delfy;
f[0].fz += delfz;
if(sp != 0.0){
      for(m=0;m<3;m++){
	 f[0].fx *= s;
	 f[0].fy *= s;
	 f[0].fz *= s;
      }
   f[0].fx -= (frc.fx = sp*eg*sdl.fx);
   f[0].fy -= (frc.fy = sp*eg*sdl.fy);
   f[0].fz -= (frc.fz = sp*eg*sdl.fz);
   f[1].fx += frc.fx;
   f[1].fy += frc.fy;
   f[1].fz += frc.fz;
      
   eg *= s;
   wnc *= s;
}
// * Br *
force[0].fx += f[0].fx;
force[0].fy += f[0].fy;
force[0].fz += f[0].fz;
// * Cl2 M[0]q: half force on each Cl *
force[1].fx += f[1].fx/2.0;
force[1].fy += f[1].fy/2.0;
force[1].fz += f[1].fz/2.0;
force[2].fx += f[1].fx/2.0;
force[2].fy += f[1].fy/2.0;
force[2].fz += f[1].fz/2.0;
//   
INTER_X += eg;
X_C += eg;
 
// Br - Cl interactions
eg = wnc = 0.0;
// to ensure no double counting in solvation shell
for(m=0;m<3; m++) //clear full array (same for GLY/Br and THA/Br) 
   f[m].fx = f[m].fy = f[m].fz = 0.;
for(j=0;j<2;j++){
      l = 1 + j; //index of Cl
      k = 0;
// *****	Determine image vector ******
	 delfx = pos[l].fx - pos[k].fx + image.fx;
	 delfy = pos[l].fy - pos[k].fy + image.fy;
	 delfz = pos[l].fz - pos[k].fz + image.fz;
	 r2 = delfx*delfx + delfy*delfy + delfz*delfz;
// *****	Get (1/r dV/dr) ******
	 dedr = gCl(r2,2,&eg,&wnc);
// *****	Resolve forces on atoms.  ******
	 f[l].fx -= (delfx *= dedr);
	 f[l].fy -= (delfy *= dedr);
	 f[l].fz -= (delfz *= dedr);
	 f[k].fx += delfx;
	 f[k].fy += delfy;
	 f[k].fz += delfz;
}   
if(sp != 0.0){
   for(m=0;m<3;m++){
         f[k].fx *= s;
	 f[k].fy *= s;
	 f[k].fz *= s;
   }
   f[0].fx += (frc.fx = sp*eg*sdl.fx);
   f[0].fy += (frc.fy = sp*eg*sdl.fy);
   f[0].fz += (frc.fz = sp*eg*sdl.fz);
   f[1].fx -= frc.fx/2.0;
   f[1].fy -= frc.fy/2.0;
   f[1].fz -= frc.fz/2.0;
   f[2].fx -= frc.fx/2.0;
   f[2].fy -= frc.fy/2.0;
   f[2].fz -= frc.fz/2.0;
   eg *= s;
   wnc *= s;
}
// * GLY *
force[k].fx += f[0].fx;
force[k].fy += f[0].fy;
force[k].fz += f[0].fz;
// Cl2 (LJq) 
force[1].fx += f[1].fx;
force[1].fy += f[1].fy;
force[1].fz += f[1].fz;
force[2].fx += f[2].fx;
force[2].fy += f[2].fy;
force[2].fz += f[2].fz;
//   
INTER_X += eg;
X_C += wnc;

}

double
gM0(r2,m,eg)
double r2, *eg;
int m; /* m is the Cl site index */
{
   double der,r,ec,q, sqrt();
   int index;
   q = clj[m].q;
   r = sqrt(r2);
   ec = q/r;
   *eg += ec;
   der = -ec/r2 ;
   return(der);
}

double
gMp(r2,m,eg)
double r2, *eg;
int m; /* m is the Cl site index */
{
   double der,r,ec,q;
   int index;
   q = clj[m].q;
   r = sqrt(r2);
   ec = q/r;
   *eg += ec;
   der = -ec/r2 ;
   return(der);
}

double
gCl(r2,m,eg,wnc)
double r2,*eg,*wnc;
int m; /* m is the Cl site index */
{
   double der,r,ir,ir6, ec,a,b,q,aa,bb;
   int bin;
   int index;
   ir = 1. / r2;
   ir6 = ir * ir * ir;
   a = clj[m].a;
   b = clj[m].b;
   q = clj[m].q;
   aa =12*a;
   bb = 6*b;
   r = sqrt(r2);
   ec = q/r;
   *wnc += ec;
   *eg +=  ( a * ir6 - b ) * ir6 + ec;
   der = (bb - aa * ir6) * ir6 * ir - ec/r2 ;
   return(der);
}
