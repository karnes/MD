#include	<md.h>
#include	<system.h>
#include	<math.h>
#include        <water.h>

// Halogen biasing potential along z
intForce()
{
    VINT = V_bias = V_w = 0.0;
    if(nBr!=1)
	return;
    if(windowOn==1)
	window();
}

window()
{
    double totMass,zeta,deriv,pow(),fabs(),zpos;
    tripd sysCm;
    int i;
    double bias_com = 0.0;
    static int printCOM = 1;
    V_w = 0;
    sysCm.fx = sysCm.fy = sysCm.fz = totMass = 0.;
    for(i = 0; i < natoms; i++) {
	sysCm.fx += mass[i] * pos[i].fx;
	sysCm.fy += mass[i] * pos[i].fy;
	sysCm.fz += mass[i] * pos[i].fz;
	totMass += mass[i];
    }
    sysCm.fx /= totMass;
    sysCm.fy /= totMass;
    sysCm.fz /= totMass;
    if (fabs(sysCm.fz - osysCmz) > 0.001 || printCOM){
	fprintf(stderr,"system center of mass = %f \n",sysCm.fz);
	printCOM = 0;
    }
    osysCmz = sysCm.fz;
    //    *** Add biasing potential for ion. Form of potential is tanh or cubic *** 
    if(biasingOn>0){
	bias_com = osysCmz;
	zpos = pos[natoms-1].fz - bias_com /*+ BIAS_G +*/; //bias_com is system COM, BIAS_G is an offset value for the equation.
	if(biasingOn==1){ //hyperbolic tangent
	    //     *** V_bias = BIAS_W (1 + tanh(BIAS_A*x) )/2             ***
	    //     ***   where x = pos[ion].fz - bias_com + BIAS_G  ***
	    V_bias = - ( BIAS_W*(1.0+tanh(BIAS_A*zpos)) / 2.0 ) / KCAL ;
	    deriv = - ( (BIAS_W*BIAS_A/2.0) / 
		    ( cosh(BIAS_A*zpos) * cosh(BIAS_A*zpos) ) ) / KCAL;
	}
	else if(biasingOn==2){  // cubic form here
	    V_bias = - (BIAS_A * zpos * zpos * zpos +  BIAS_B * zpos * zpos 
		    + BIAS_C * zpos + BIAS_D ) / KCAL;
	    deriv = - (3.0 * BIAS_A * zpos * zpos + 2.0 *BIAS_B * zpos
		    + BIAS_C ) / KCAL;
	}
	force[natoms-1].fz -= deriv;
	for(i=0;i<natoms;i++){
	    force[i].fz+=deriv*mass[i]/totMass;
	}
	VINT+=V_bias;
    }
    //** Add the window potential if solute outside the window ***
    zpos = pos[natoms-1].fz - sysCm.fz - center_w;
    zeta = fabs(zpos)-width_w/2.;
    if(zeta > 0.){// outside window, calculate energy 
	if (zeta > 1.5/**width_w*/){
	    fprintf(stderr,"window: too much outside the window\n");
	    exit(1);
	}
	else { // 0 < zeta < width_w
	    V_w = pot_w*pow(zeta,power_w-1);
	    deriv = V_w*power_w*sgn(zpos);
	    V_w *= zeta;
	    force[natoms-1].fz -= deriv;
            for (i=0;i<natoms;i++){
	        force[i].fz += (deriv * mass[i]/totMass);
	    }
	    VINT += V_w;
	}
    }
}

/*
#include	<md.h>
#include	<system.h>
#include	<math.h>

// R_THA-X windowing (and biasing) potentials
intForce()
{

   int i,j,k,nX;
   tripd Rn,fTR,fXR,fRnet,Rvec,Xcom,Tcom;
   double Xtotmass,Ttotmass,fRt;
//   XwindowOn = 1;
//   Xmin = 15.0;
   if(nBr==1){
      nX = 1;
   }
   else if(nCl2==1){
      nX = 2;
   }
   else if(nTS==1){
      nX = 3;
   }
   else{
      return;
   }
   if(nTHA != 1){ 
      return;
   }

   V_window = VINT = V_bias = V_Xwindow = 0.0;
*/
   /***	Restrict the distance between the THA CoM and Br- ***/
/*   if(window_w > zwall){ // no window 
      return;
   }
   else if(window_w < 0.0){ // if the size of the region is negative
      fprintf(stderr,"intForce.c: window < 0 ? error. exiting...\n");
      exit(0);	
   }
   else{
      window(nX);
   }
*/
   /*** restrict halogen to liquid/vapor interface ***/
/*   Xwindow(nX);
   VINT += V_window + V_bias + V_Xwindow;
}
*/
/* radial windowing potential */
/*
window(nX)
int nX;
{
   int i,j,k;
   double r2,r,deriv;
   double dx,dy,dz,zeta;
   tripd sdist,Tcom,Xcom;
   double Ttotmass,Xtotmass;

   if(windowOn!=1)
      return;
   // Find THA,X center of mass
   Ttotmass = Tcom.fx = Tcom.fy = Tcom.fz = 0.0;
   Xtotmass = Xcom.fx = Xcom.fy = Xcom.fz = 0.0;
   for(i=0;i<THAsites;i++){
      k = nGLY*GLYsites + i;
      Tcom.fx+=pos[k].fx*mass[k];
      Tcom.fy+=pos[k].fy*mass[k];
      Tcom.fz+=pos[k].fz*mass[k];
      Ttotmass+=mass[k];
   }
   Tcom.fx/=Ttotmass;
   Tcom.fy/=Ttotmass;
   Tcom.fz/=Ttotmass;

   for(i=0;i<nX;i++){
      k = nGLY*GLYsites+nTHA*THAsites+i;
      Xcom.fx+=pos[k].fx*mass[k];
      Xcom.fy+=pos[k].fy*mass[k];
      Xcom.fz+=pos[k].fz*mass[k];
      Xtotmass+=mass[k];
   }
   Xcom.fx/=Xtotmass;
   Xcom.fy/=Xtotmass;
   Xcom.fz/=Xtotmass;
*/
   /***	Determine THA com - Br- distance    ***/
/*
   sdist.fx = Tcom.fx-Xcom.fx;
   sdist.fy = Tcom.fy-Xcom.fy;
   sdist.fz = Tcom.fz-Xcom.fz;
   mvimage(&sdist);
   r2 = sdist.fx*sdist.fx +  sdist.fy*sdist.fy +  sdist.fz*sdist.fz;
   r = sqrt(r2);

   zeta = fabs(r-center_r) - window_w/2.0;
   #define TethCon 1.0 
   if(zeta > 0.0){
      V_window += TethCon*zeta*zeta*zeta;
      deriv = -3.0*TethCon*zeta*zeta*sgn(r-center_r);
      dx = sdist.fx/r;
      dy = sdist.fy/r;
      dz = sdist.fz/r;
      for(i=0;i<THAsites;i++){
         k=nGLY*GLYsites+i;
         force[k].fx+=deriv*dx*mass[k]/Ttotmass;
         force[k].fy+=deriv*dy*mass[k]/Ttotmass;
         force[k].fz+=deriv*dz*mass[k]/Ttotmass;
      }
      for(i=0;i<nX;i++){
         k=nGLY*GLYsites+nTHA*THAsites+i;
         force[k].fx-=deriv*dx*mass[k]/Xtotmass;
         force[k].fy-=deriv*dy*mass[k]/Xtotmass;
         force[k].fz-=deriv*dz*mass[k]/Xtotmass;
      }
   }   

   if(biasOn==1){
      V_bias = -(a_bias*(r-o_bias)*(r-o_bias) + b_bias*(r-o_bias))/KCAL;
      deriv = (a_bias*(r-o_bias) + b_bias)/KCAL;
      dx = sdist.fx/r;
      dy = sdist.fy/r;
      dz = sdist.fz/r;
      for(i=0;i<THAsites;i++){
         k=nGLY*GLYsites+i;
         force[k].fx+=deriv*dx*mass[k]/Ttotmass;
         force[k].fy+=deriv*dy*mass[k]/Ttotmass;
         force[k].fz+=deriv*dz*mass[k]/Ttotmass;
      }
      for(i=0;i<nX;i++){
         k=nGLY*GLYsites+nTHA*THAsites+i;
         force[k].fx-=deriv*dx*mass[k]/Xtotmass;
         force[k].fy-=deriv*dy*mass[k]/Xtotmass;
         force[k].fz-=deriv*dz*mass[k]/Xtotmass;
      }
   }
}
*/
/* halogen (X) windowing potential */
/*Xwindow(nX)
int nX;
{
   if(XWindowOn!=1)
      return;
   int i,j,k;
   double deriv;
   double dz,zeta;
   tripd sdist,Tcom,Xcom;
   double sysMass,Xtotmass;

   if(windowOn!=1)
      return;
   // Find X z-axis center of mass
   sysMass = Xtotmass = Xz = 0.0;
   for(i=0;i<nX;i++){
      k = nGLY*GLYsites+nTHA*THAsites+i;
      Xz+=pos[k].fz*mass[k];
      Xtotmass+=mass[k];
   }
   Xz/=Xtotmass;
   for(i=0;i<natoms-nX;i++)
      sysMass+=mass[i];

   zeta = Xmin - (Xz - syscomz);
   #define WinCon 1.0 
   if(zeta > 0.0){
      V_window += WinCon*zeta*zeta*zeta;
      deriv = -3.0*WinCon*zeta*zeta;
      for(i=0;i<nX;i++){
         k=nGLY*GLYsites+nTHA*THAsites+i;
         force[k].fz-=deriv*mass[k]/Xtotmass;
      }
      for(i=0;i<natoms-nX;i++){
         force[i].fz+=deriv*mass[i]/sysMass;
      }
   }   
}
i*/
