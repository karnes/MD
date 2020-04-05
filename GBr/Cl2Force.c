#include	<md.h>
#include	<system.h>
#include	 <math.h>
#define ldist 0.8
Cl2Force()
{
int	i, j, k, l,C1,C2,m;
double	totMass, r, r2, dedr, delfx, delfy, delfz;
double	eg, wnc, s, sp, deriv, r_12, ir12;
tripd	d,com, frc, image, sdl, f[THAsites+2]; 
tripd	r_a,r_b,r_c,r_1,r_2,r1mr2,r2mr1; 
double	gCl(),tCl();
double  gM0(),tM0(),gMp();
double sqrt();

//M0e = Me = Cle = 0.0;
Cl2bondE=0.0;

if(nCl2!=1){
   return(0);
}

// intra Cl2
i=nGLY*GLYsites+nTHA*THAsites;
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
Cl2bondE += eg;
//
cosXz = r1mr2.fz/r_12;
Xz = (pos[i].fz+pos[j].fz)/2.0;
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
for(i=0;i<nGLY;i++){
   k = i*GLYsites;
   // *****	Determine image vector, GLY parent atom - b (center of Cl2).  ******
   image.fx = -(sdl.fx = r_b.fx - pos[k].fx);
   image.fy = -(sdl.fy = r_b.fy - pos[k].fy);
   image.fz = -(sdl.fz = r_b.fz - pos[k].fz);
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
////
//// BEGIN PROBLEMATIC PART
////


// GLY - a interactions -- project forces onto Cl atoms
   eg = wnc = 0.0;
   for(m=0;m<THAsites;m++){ //clear full array 
      f[m].fx = f[m].fy = f[m].fz = 0.;
   }
   l = nGLY*GLYsites + nTHA*THAsites; // first Cl atom
   k = i*GLYsites; //center atom of GLY i
// *****	Determine image vector ******
   for(m=0;m<GLYsites;m++){ //
	 delfx = r_a.fx - pos[k+m].fx + image.fx;
	 delfy = r_a.fy - pos[k+m].fy + image.fy;
	 delfz = r_a.fz - pos[k+m].fz + image.fz;
	 r2 = delfx*delfx + delfy*delfy + delfz*delfz;
// *****   Get (1/r dV/dr) ******
	 dedr = gM0(r2,m,0,&eg);
// *****   Resolve forces on atoms.  ******
	 f[GLYsites].fx -= delfx*dedr*(1.0 + ldist/r_12 - ir12*r1mr2.fx*r1mr2.fx);
	 f[GLYsites].fy -= delfy*dedr*(1.0 + ldist/r_12 - ir12*r1mr2.fy*r1mr2.fy);
	 f[GLYsites].fz -= delfz*dedr*(1.0 + ldist/r_12 - ir12*r1mr2.fz*r1mr2.fz);

	 f[GLYsites+1].fx -= delfx*dedr*( -ldist/r_12 + ir12*r1mr2.fx*r1mr2.fx);
	 f[GLYsites+1].fy -= delfy*dedr*( -ldist/r_12 + ir12*r1mr2.fy*r1mr2.fy);
	 f[GLYsites+1].fz -= delfz*dedr*( -ldist/r_12 + ir12*r1mr2.fz*r1mr2.fz);

	 f[m].fx += delfx*dedr;
	 f[m].fy += delfy*dedr;
	 f[m].fz += delfz*dedr;
   // cutoff stuff here (deleted for debugging)
   }//end m loop
   for(m=0;m<GLYsites;m++){
      force[k+m].fx += f[m].fx;
      force[k+m].fy += f[m].fy;
      force[k+m].fz += f[m].fz;
   }

   force[l].fx += f[GLYsites].fx;
   force[l].fy += f[GLYsites].fy;
   force[l].fz += f[GLYsites].fz;
   force[l+1].fx += f[GLYsites+1].fx;
   force[l+1].fy += f[GLYsites+1].fy;
   force[l+1].fz += f[GLYsites+1].fz;

   XGLYV += eg;


// GLY - c interactions -- project forces onto Cl atoms
   eg = wnc = 0.0;
   for(m=0;m<THAsites;m++){ //clear full array 
      f[m].fx = f[m].fy = f[m].fz = 0.;
   }
   l = nGLY*GLYsites + nTHA*THAsites; // second Cl atom
   k = i*GLYsites; //center atom of GLY i
// *****	Determine image vector ******
   for(m=0;m<GLYsites;m++){ //
	 delfx = r_c.fx - pos[k+m].fx + image.fx;
	 delfy = r_c.fy - pos[k+m].fy + image.fy;
	 delfz = r_c.fz - pos[k+m].fz + image.fz;
	 r2 = delfx*delfx + delfy*delfy + delfz*delfz;
// *****   Get (1/r dV/dr) ******
	 dedr = gM0(r2,m,0,&eg);
// *****   Resolve forces on atoms.  ******
	 f[GLYsites+1].fx -= delfx*dedr*(1.0 + ldist/r_12 - ir12*r1mr2.fx*r1mr2.fx);
	 f[GLYsites+1].fy -= delfy*dedr*(1.0 + ldist/r_12 - ir12*r1mr2.fy*r1mr2.fy);
	 f[GLYsites+1].fz -= delfz*dedr*(1.0 + ldist/r_12 - ir12*r1mr2.fz*r1mr2.fz);

	 f[GLYsites].fx -= delfx*dedr*( -ldist/r_12 + ir12*r1mr2.fx*r1mr2.fx);
	 f[GLYsites].fy -= delfy*dedr*( -ldist/r_12 + ir12*r1mr2.fy*r1mr2.fy);
	 f[GLYsites].fz -= delfz*dedr*( -ldist/r_12 + ir12*r1mr2.fz*r1mr2.fz);

	 f[m].fx += delfx*dedr;
	 f[m].fy += delfy*dedr;
	 f[m].fz += delfz*dedr;
   // cutoff stuff here (deleted for debugging)
   }//end m loop
   for(m=0;m<GLYsites;m++){
      force[k+m].fx += f[m].fx;
      force[k+m].fy += f[m].fy;
      force[k+m].fz += f[m].fz;
   }

   force[l+1].fx += f[GLYsites+1].fx;
   force[l+1].fy += f[GLYsites+1].fy;
   force[l+1].fz += f[GLYsites+1].fz;
   force[l].fx += f[GLYsites].fx;
   force[l].fy += f[GLYsites].fy;
   force[l].fz += f[GLYsites].fz;

   XGLYV += eg;

////
//// END PROBLEMATIC PART
////
   // b - GLY interactions
   eg = wnc = 0.0;
   l = nGLY*GLYsites + nTHA*THAsites; //index of first Cl
   for(m=0;m<(THAsites+2); m++) //clear full array (same for GLY/Br and THA/Br) 
      f[m].fx = f[m].fy = f[m].fz = 0.;
   for(m=0;m<GLYsites;m++){ // *Loop over atoms in GLY molecules*
// ***** Determine image vector ******
      delfx = r_b.fx - pos[k+m].fx + image.fx;
      delfy = r_b.fy - pos[k+m].fy + image.fy;
      delfz = r_b.fz - pos[k+m].fz + image.fz;
      r2 = delfx*delfx + delfy*delfy + delfz*delfz;
      // *****	Get (1/r dV/dr) ******
      dedr = gM0(r2,m,1,&eg);
      // *****	Resolve forces on atoms.  ******
      f[GLYsites].fx -= (delfx *= dedr);
      f[GLYsites].fy -= (delfy *= dedr);
      f[GLYsites].fz -= (delfz *= dedr);
      f[m].fx += delfx;
      f[m].fy += delfy;
      f[m].fz += delfz;
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
   l=nGLY*GLYsites+nTHA*THAsites;
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
   k = i*GLYsites;
   // * Cl2 M[0]q: half force on each Cl *
   force[l].fx += f[GLYsites].fx/2.0;
   force[l].fy += f[GLYsites].fy/2.0;
   force[l].fz += f[GLYsites].fz/2.0;
   force[l+1].fx += f[GLYsites].fx/2.0;
   force[l+1].fy += f[GLYsites].fy/2.0;
   force[l+1].fz += f[GLYsites].fz/2.0;
   //   
   XGLYV += eg;
   XGLYC += eg;
   

 
// GLY - Cl interactions
   eg = wnc = 0.0;
// to ensure no double counting in solvation shell
   for(m=0;m<GLYsites;m++)
      gflag[m] = 0;
   for(m=0;m<THAsites+2; m++) //clear full array (same for GLY/Br and THA/Br) 
      f[m].fx = f[m].fy = f[m].fz = 0.;
   for(j=0;j<2;j++){
      l = nGLY*GLYsites + nTHA*THAsites + j; //index of Cl
      k = i*GLYsites;
      for(m=0;m<GLYsites;m++){ //Loop over atoms in GLY molecules
// *****	Determine image vector ******
	 delfx = pos[l].fx - pos[k+m].fx + image.fx;
	 delfy = pos[l].fy - pos[k+m].fy + image.fy;
	 delfz = pos[l].fz - pos[k+m].fz + image.fz;
	 r2 = delfx*delfx + delfy*delfy + delfz*delfz;
// *****	Get (1/r dV/dr) ******
	 dedr = gCl(r2,m,2,&eg,&wnc);
// *****	Resolve forces on atoms.  ******
	 f[GLYsites+j].fx -= (delfx *= dedr);
	 f[GLYsites+j].fy -= (delfy *= dedr);
	 f[GLYsites+j].fz -= (delfz *= dedr);
	 f[m].fx += delfx;
	 f[m].fy += delfy;
	 f[m].fz += delfz;
      }
   }   
   if(sp != 0.0){
      for(m=0;m<GLYsites+2;m++){
         f[m].fx *= s;
	 f[m].fy *= s;
	 f[m].fz *= s;
      }
      f[0].fx += (frc.fx = sp*eg*sdl.fx);
      f[0].fy += (frc.fy = sp*eg*sdl.fy);
      f[0].fz += (frc.fz = sp*eg*sdl.fz);
      f[GLYsites].fx -= frc.fx/2.0;
      f[GLYsites].fy -= frc.fy/2.0;
      f[GLYsites].fz -= frc.fz/2.0;
      f[GLYsites+1].fx -= frc.fx/2.0;
      f[GLYsites+1].fy -= frc.fy/2.0;
      f[GLYsites+1].fz -= frc.fz/2.0;
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
   // Cl2 (LJq) 
   force[l].fx += f[GLYsites].fx;
   force[l].fy += f[GLYsites].fy;
   force[l].fz += f[GLYsites].fz;
   force[l+1].fx += f[GLYsites+1].fx;
   force[l+1].fy += f[GLYsites+1].fy;
   force[l+1].fz += f[GLYsites+1].fz;
   //   
   XGLYV += eg;
   XGLYC += wnc;

} // end GLY (k) loop
/*
if(nTHA==1){
   for(i=0;i<THAsites;i++){
      k = nGLY*GLYsites + i;
// *****	Determine image vector, THA atom - M[0] (center of Cl2).  ******
      image.fx = -(sdl.fx = pos[k].fx - M[0].fx);
      image.fy = -(sdl.fy = pos[k].fy - M[0].fy);
      image.fz = -(sdl.fz = pos[k].fz - M[0].fz);
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
     
     // M[0] - THA interactions
     for(m=0;m<(THAsites+2); m++) //clear full array (same for GLY/X and THA/X) 
        f[m].fx = f[m].fy = f[m].fz = 0.;
     eg = wnc = 0.0;
     l = nGLY*GLYsites + nTHA*THAsites; //index of first Cl
// *****	Get (1/r dV/dr) ******
     dedr = tM0(r2,i,1,&eg);
// *****	Resolve forces on atoms.  ******
      f[1].fx += (delfx * dedr);// M0
      f[1].fy += (delfy * dedr);
      f[1].fz += (delfz * dedr);
      f[0].fx -= (delfx * dedr);// THA atom
      f[0].fy -= (delfy * dedr);
      f[0].fz -= (delfz * dedr);
      if(sp != 0.0){
	 for(m=0;m<2;m++){
	       f[m].fx *= s;
	       f[m].fy *= s;
	       f[m].fz *= s;
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
      k=nGLY*GLYsites+i;
      l=nGLY*GLYsites+nTHA*THAsites;
      // * THA *
      force[k].fx += f[0].fx;
      force[k].fy += f[0].fy;
      force[k].fz += f[0].fz;
      k = i*GLYsites;
      // * Cl2 M[0]q: half force on each Cl *
      force[l].fx += f[1].fx/2.0;
      force[l].fy += f[1].fy/2.0;
      force[l].fz += f[1].fz/2.0;
      force[l+1].fx += f[1].fx/2.0;
      force[l+1].fy += f[1].fy/2.0;
      force[l+1].fz += f[1].fz/2.0;
      //   
      //
      XTHAV += eg;
      XTHAC += wnc;
//	M0e+=eg;
      eg = wnc = 0.0;
      for(m=0;m<(THAsites+2); m++) //clear full array 
         f[m].fx = f[m].fy = f[m].fz = 0.;
// THA - M,M' interactions -- project forces onto Cl atoms
      for(j=0;j<2;j++){ // loop for each M...
         k = nGLY*GLYsites+i;
         if(j==0){
            C1=0;
	    C2=1;
            Mn=1;
         }
         else if(j==1){
	    C1=1;
	    C2=0;
	    Mn=2;
         }
//	 Get (1/r dV/dr)
//	 Determine image vector
	 delfx = pos[k].fx - M[Mn].fx + image.fx;
	 delfy = pos[k].fy - M[Mn].fy + image.fy;
	 delfz = pos[k].fz - M[Mn].fz + image.fz;
	 r2 = delfx*delfx + delfy*delfy + delfz*delfz;
	 //get force on THA atom
	 dedr = tM0(r2,i,0,&eg);
	 //calculate forces on Cl(a),Cl(b)
	 //get force on Cl close to M[Mn]
	 fCl.fx = (delfx*dedr)*(1+0.8/rClCl-0.8*(Clvec.fx*Clvec.fx)/(rClCl*rClCl*rClCl));
	 fCl.fy = (delfy*dedr)*(1+0.8/rClCl-0.8*(Clvec.fy*Clvec.fy)/(rClCl*rClCl*rClCl));
	 fCl.fz = (delfz*dedr)*(1+0.8/rClCl-0.8*(Clvec.fz*Clvec.fz)/(rClCl*rClCl*rClCl));
	 //...and force on Cl further from M[Mn]
	 fClp.fx = (delfx*dedr)*(-0.8/rClCl+0.8*(Clvec.fx*Clvec.fx)/(rClCl*rClCl*rClCl));
	 fClp.fy = (delfy*dedr)*(-0.8/rClCl+0.8*(Clvec.fy*Clvec.fy)/(rClCl*rClCl*rClCl));
	 fClp.fz = (delfz*dedr)*(-0.8/rClCl+0.8*(Clvec.fz*Clvec.fz)/(rClCl*rClCl*rClCl));
	 
//	 Resolve forces on atoms. 
	 f[1+C1].fx += fCl.fx;
	 f[1+C1].fy += fCl.fy;
	 f[1+C1].fz += fCl.fz;
	 f[1+C2].fx += fClp.fx;
	 f[1+C2].fy += fClp.fy;
	 f[1+C2].fz += fClp.fz;
	 f[0].fx -= delfx*dedr;
	 f[0].fy -= delfy*dedr;
	 f[0].fz -= delfz*dedr;
      }
      if(sp != 0.0){
	 for(m=0;m<3;m++){
	       f[m].fx *= s;
	       f[m].fy *= s;
	       f[m].fz *= s;
	 }
	 f[0].fx -= (frc.fx = sp*eg*sdl.fx);
	 f[0].fy -= (frc.fy = sp*eg*sdl.fy);
	 f[0].fz -= (frc.fz = sp*eg*sdl.fz);
	 f[1].fx += frc.fx*0.5;
	 f[1].fy += frc.fy*0.5;
	 f[1].fz += frc.fz*0.5;
	 f[2].fx += frc.fx*0.5;
	 f[2].fy += frc.fy*0.5;
	 f[2].fz += frc.fz*0.5;
	 eg *= s;
	 wnc *= s;
      }
      k=nGLY*GLYsites+i;
      l = nGLY*GLYsites+nTHA*THAsites;
      // THA 
      force[k].fx += f[0].fx;
      force[k].fy += f[0].fy;
      force[k].fz += f[0].fz;
      k=nGLY*GLYsites+i;
      // M projected on C and Cl'(q) 
      force[l].fx += f[1].fx;
      force[l].fy += f[1].fy;
      force[l].fz += f[1].fz;
      force[l+1].fx += f[2].fx;
      force[l+1].fy += f[2].fy;
      force[l+1].fz += f[2].fz;
      //   
      //
      XTHAV += eg;
      XTHAC += wnc;
      
//GLY - Cl interactions
      eg = wnc = 0.0;
      for(m=0;m<(THAsites+2); m++) //clear full array (same for GLY/Br and THA/Br) 
         f[m].fx = f[m].fy = f[m].fz = 0.;
      for(j=0;j<2;j++){
         l = nGLY*GLYsites + nTHA*THAsites + j; //index of Cl
         k = nGLY*GLYsites+i;
// *****	Determine image vector ******
	 delfx = pos[k].fx - pos[l].fx + image.fx;
	 delfy = pos[k].fy - pos[l].fy + image.fy;
	 delfz = pos[k].fz - pos[l].fz + image.fz;
	 r2 = delfx*delfx + delfy*delfy + delfz*delfz;
// *****	Get (1/r dV/dr) ******
	 dedr = tCl(r2,i,2,&eg,&wnc);
// *****	Resolve forces on atoms.  ******
	 f[1+j].fx += (delfx * dedr);
	 f[1+j].fy += (delfy * dedr);
	 f[1+j].fz += (delfz * dedr);
	 f[0].fx -= delfx * dedr; // THA atom
	 f[0].fy -= delfy * dedr;
	 f[0].fz -= delfz * dedr;
      }
      if(sp != 0.0){
	 for(m=0;m<3;m++){
	    f[m].fx *= s;
	    f[m].fy *= s;
	    f[m].fz *= s;
         }
	 f[0].fx -= (frc.fx = sp*eg*sdl.fx);
	 f[0].fy -= (frc.fy = sp*eg*sdl.fy);
	 f[0].fz -= (frc.fz = sp*eg*sdl.fz);
	 f[1].fx += frc.fx/2.0;
	 f[1].fy += frc.fy/2.0;
	 f[1].fz += frc.fz/2.0;
	 f[2].fx += frc.fx/2.0;
	 f[2].fy += frc.fy/2.0;
	 f[2].fz += frc.fz/2.0;
	 eg *= s;
	 wnc *= s;
      }
      k=nGLY*GLYsites+i;
      l = nGLY*GLYsites+nTHA*THAsites;	// * THA *
      force[k].fx += f[0].fx;
      force[k].fy += f[0].fy;
      force[k].fz += f[0].fz;
      k = nGLY* GLYsites;
      // Cl2 (LJq) 
      force[l].fx += f[1].fx;
      force[l].fy += f[1].fy;
      force[l].fz += f[1].fz;
      force[l+1].fx += f[2].fx;
      force[l+1].fy += f[2].fy;
      force[l+1].fz += f[2].fz;
      //   
      //
      XTHAV += eg;
      XTHAC += wnc;
      

   } // end THAsites (i) loop

} // end if THA loop
*/
}

double
gM0(r2,m,n,eg)
double r2, *eg;
int m, /* m is the GLY atom index */
    n; /* n is the Cl (GLYsites+1) atom index */
{
   double der,r,ec,q, sqrt();
   int index;
   q = solGLYlj[m][n].q;
   r = sqrt(r2);
   if(tc%grdt==0){
      index = (int)(r/binRDF);
      if(index<300){// g(r)'s
         if(m==0||m==4||m==9){//C-M0
            grGLY[13][index]+=nGLY*1.0/3.0;
	    if(r<gGCXmax[2]){
	       gGCX[2]++;
	    }
         }
         else if(m==1||m==5||m==10){//O-M0
            grGLY[14][index]+=nGLY*1.0/3.0;
	    if(r<gGOXmax[2]){
	       gGOX[2]++;
	    }
         }
         else if(m==2||m==6||m==11){//H(oh)-M0
            grGLY[15][index]+=nGLY*1.0/3.0;
	    if(r<gGHXmax[2]){
	       gGHX[2]++;
	    }
         }
      }
   }
   ec = q/r;
   *eg += ec;
   der = -ec/r2 ;
   return(der);
}

double
gMp(r2,m,n,eg)
double r2, *eg;
int m, /* m is the GLY atom index */
    n; /* n is the a or c virtual site index */
{
   double der,r,ec,q;
   int index;
   q = solGLYlj[m][n].q;
   r = sqrt(r2);
   ec = q/r;
   *eg += ec;
   der = -ec/r2 ;
   return(der);
}

double
gCl(r2,m,n,eg,wnc)
double r2,*eg,*wnc;
int m, /* m is the THA atom index */
    n; /* n is the Cl (THAsites+1) index */
{
   double der,r,ir,ir6, ec,a,b,q,aa,bb;
   int bin;
   int index;
   ir = 1. / r2;
   ir6 = ir * ir * ir;
   a = solGLYlj[m][n].a;
   b = solGLYlj[m][n].b;
   q = solGLYlj[m][n].q;
   aa =12*a;
   bb = 6*b;
   r = sqrt(r2);
   if(tc%grdt==0){
      index = (int)(r/binRDF);
      if(index<300){// g(r)'s
         if(m==0||m==4||m==9){//C-Cl
            grGLY[7][index]+=nGLY*1.0/(2.0*3.0);
	    
	    if(gflag[m]==0 && r<gGCXmax[0]){
	       gGCX[0]++;
	       gflag[m] = 1;
	    }
         }
         else if(m==1||m==5||m==10){//O-Cl
            grGLY[8][index]+=nGLY*1.0/(2.0*3.0);
	    if(gflag[m]==0  && r<gGOXmax[0]){
	       gGOX[0]++;
	       gflag[m] = 1;
	    }
         }
         else if(m==2||m==6||m==11){//H(oh)-Cl
            grGLY[9][index]+=nGLY*1.0/(2.0*3.0);
	    if(gflag[m]==0 && r<gGHXmax[0]){
	       gGHX[0]++;
	       gflag[m] = 1;
	    }
         }
      }
   }
   ec = q/r;
   *wnc += ec;
   *eg +=  ( a * ir6 - b ) * ir6 + ec;
   der = (bb - aa * ir6) * ir6 * ir - ec/r2 ;
//	fprintf(stderr,"%f\n",KCAL*(*eg));
   return(der);
}

double
tM0(r2,m,n,eg)
double r2, *eg;
int m, /* m is the THA atom index */
    n; /* n is the b virtual atom index */
{
   double der,r,ec,q;
   q = solTHAlj[m][n].q;
   r = sqrt(r2);
   ec = q/r;
   *eg += ec;
   der = -ec/r2 ;
   return(der);
}

double
tCl(r2,m,n,eg,wnc)
double r2,*eg,*wnc;
int m, /* m is the THA atom index */
    n; /* n is the Cl index */
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


