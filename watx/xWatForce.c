#include	<md.h>
#include	<system.h>
#include	<water.h>
#include	<math.h>
/* forces pertaining to the transferring water molecule: first water in input files.*/

xWatForce(pos, force)
	tripd	*pos;
	tripd	*force;
{
	int	i,j,k,nw;
	double c, pe, ElecPot();
//	tripd wcom;

	k = natoms - nsolute;
	nw = k - 14*nNIT;
	V_teth = 0.0; /*total energy tethering fixShel waters */
	for(i=0;i<nw;i+=3){
	   if(fixShel>0){
/*
	      wcom.fx = wcom.fy = wcom.fz = 0.0;
	      for(j=0;j<3;j++){
		 wcom.fx += pos[j].fx*mass[j];
		 wcom.fy += pos[j].fy*mass[j];
		 wcom.fz += pos[j].fz*mass[j];
	      }
	      wcom.fx/=WMass;
	      wcom.fy/=WMass;
	      wcom.fz/=WMass;
*/
	      if(i/3==fixw[0] || i/3==fixw[1] || i/3==fixw[2]){
	         teth(i); // pass actual index of tethered water
	      }
	   }
	}
	VINT += V_teth;
	if(EFieldOn) 
	   VINT += ElecPot(pos,force);
/***	Restrict the center of mass of the solute to a region.***/
	if(width_w > zwall) /* no window */
	   return;
	if(width_w < 0.) /* if the size of the region is negative*/
	   freez_cm(); /* freez center of mass	*/
	else /* use windowing potential centered at the closest integer
	      * and width width_w outside which there is a power law
	      * potential of power power_w and coeficients pot_w
	      */
	   window();
	if(tc==0){
	   if(fixShel>0){
	      if(fixed!=fixShel){
		fprintf(stderr, "intForce.c: ion lost its fixed hydration shell. terminating run.\n");
		exit(1);
	      }
	   }
	}

}
#define	EXBIAS	0.5
//#define BIAS_W 11.0
//#define BIAS_A 0.18
//#define BIAS_M 3.0
//#define BIAS_B 5.0
//#define BIAS_a 2.9
//#define BIAS_b -40.8
//#define BIAS_c 189.4
//#define BIAS_d -267
//#define Gibbs_surf 7
//#define BIAS_MOD 1.0 // multiply by this factor to scale the biasing force
window()
{
double totMass,zeta,deriv,pow(),fabs(),zpos;
tripd sysCm; //,scom;
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
	fprintf(stderr,"C.O.M = %f \n",sysCm.fz);
	printCOM = 0;
}
osysCmz = sysCm.fz;
if(nsolute == 0){
	/*** Add biasing potential for water x. Form of potential is *** 
	 *** V_bias = bias_a (1 + tanh(bias_b*x) )/2             ***
	 ***   where x = pos[water x].fz - bias_gs - bias_com  ***/
	if(biasingOn==1){
	    bias_com = osysCmz;
	    zpos = watx.fz - bias_gs + (double)(1-KillFinger)*(bias_com - 9.8);
//  hyperbolic tangent form here
	    V_bias = - ( bias_a*(1+tanh(bias_b*zpos)) / 2.0 ) / KCAL;
	    deriv = - ( (bias_a*bias_b/2.0) / 
		( cosh(bias_b*zpos) * cosh(bias_b*zpos) ) ) / KCAL;

//  linear form here
//	    V_bias = - (bias_a * zpos) + bias_b  / KCAL;
//	    deriv = - bias_a / KCAL;
	    
//  quadratic form here
/*	    V_bias = - (bias_a *(zpos - bias_b) * (zpos - bias_b) + bias_c ) / KCAL;
	    deriv = - 2 * bias_a * (zpos - bias_b) / KCAL;
*/
/*  cubic form here
	    V_bias = - (bias_a * zpos * zpos * zpos +  bias_b * zpos * zpos 
		        + bias_c * zpos + bias_d ) / KCAL;
	    deriv = - (3.0 * bias_a * zpos * zpos + 2.0 *bias_b * zpos
			+ bias_c ) / KCAL;
*/
	    force[0].fz -= deriv*mass[0]/WMass;
	    force[1].fz -= deriv*mass[1]/WMass;
	    force[2].fz -= deriv*mass[2]/WMass;
	    for(i=3;i<natoms;i++){
		force[i].fz+=deriv*mass[i]/(totMass-WMass);
	    }
	    VINT+=V_bias;
	}
	/*** Add the window potential if solute outside the window ***/
	if(KillFinger==0){
	   zpos = watx.fz - sysCm.fz - center_w;
	}
	else{
	   zpos = watx.fz - 9.8 - center_w;
	}
	zeta = fabs(zpos)-width_w/2.;
	if(zeta > 0.){/* outside window, calculate energy */
	   if (zeta > 1.5*width_w){
	      fprintf(stderr,"window: too much outside the window\n");
	      exit(1);
	   }
	   else { /* 0 < zeta < width_w*/
	      V_w = pot_w*pow(zeta,power_w-1);
	      deriv = V_w*power_w*sgn(zpos);
	      V_w *= zeta;
	      force[0].fz -= deriv*mass[0]/WMass;
	      force[1].fz -= deriv*mass[1]/WMass;
	      force[2].fz -= deriv*mass[2]/WMass;
	      if(KillFinger==0){
		 for(i=3;i<natoms;i++){
		    force[i].fz += deriv * mass[i]/(totMass-WMass);
		 }
		 VINT += V_w;
	      }
	      else{
		 VINT += V_w;
	      }
	   }
	}
}
else{
	if(nsolute == 2){
	    int solA, solB;
	    double zpos, muA, muB;
	    solA = natoms-2;
	    solB = natoms-1;
	    muA = mass[solA]/(mass[solA] + mass[solB]);
	    muB = 1. - muA;
	    zpos = muA * pos[solA].fz + muB * pos[solB].fz - center_w - sysCm.fz;
	    zeta = fabs(zpos)-width_w/2.;
	    if (zeta < 0.) /*inside window, no interactions*/
		return;
	    else {/*zeta > 0.; outside window, calculate energy */
		if (zeta > width_w){
	   	    fprintf(stderr,"window:too much outside the window\n");
		    exit(1);
		}
		else { /* 0 < zeta < width_w*/
		    V_w = pot_w*pow(zeta,power_w-1);
		    deriv = V_w*power_w*sgn(zpos);
		    V_w *= zeta;
		    force[solA].fz -= muA*deriv;
		    force[solB].fz -= muB*deriv;
		    for (i=0;i<natoms;i++)
			force[i].fz += (deriv * mass[i]/totMass);
		    VINT += V_w;
		}
	    }
	}
	else{
	    fprintf(stderr,"No window setup for polyatomic solute\n");
	    exit(1);
	}
}
}

freez_cm()
{
   fprintf(stderr,"ERROR: freez_cm() not implemented. exiting...\n");
   exit(0);
/*
if (nsolute == 1){
	force[natoms-1].fz = 0.;
}
else	{
	if (nsolute == 2){
		int solA, solB;
		double zCMfrc, muA, muB, zCMvel;
		solA = natoms-2;
		solB = natoms -1;
		muA = mass[solA]/(mass[solA] + mass[solB]);
		muB = 1. - muA;
		zCMfrc = force[solA].fz + force[solB].fz;
		force[solA].fz -= muA*zCMfrc;
		force[solB].fz -= muB*zCMfrc;
}
}*/
}

/* tethering forces */
teth(tw)
//tripd *ri;/* ri[0], ri[1] and ri[2] are position of the O H1 and H2 of water*/
int tw;
//tripd wcom;
//tripd *rj;/* rj[0] is the position of the 'solute' water O atom*/
//tripd *fi;/* fi[0], fi[1] and fi[2] are forces on the O H1 and H2 of water*/
//tripd *fj;/* force on solute atom*/
{
int i;
double sqrt(),r2,rOI,deriv;
double dx,dy,dz;
tripd tcom,image,sdist;

/***	Determine Oxygen-solute image vector			***/
tcom.fx = tcom.fy = tcom.fz = 0.0;
for(i=0;i<3;i++){
   tcom.fx += pos[tw+i].fx*mass[tw+i];
   tcom.fy += pos[tw+i].fy*mass[tw+i];
   tcom.fz += pos[tw+i].fz*mass[tw+i];
}
tcom.fx/=WMass;
tcom.fy/=WMass;
tcom.fz/=WMass;

image.fx = -(sdist.fx = tcom.fx - watx.fx);
image.fy = -(sdist.fy = tcom.fy - watx.fy);
image.fz = -(sdist.fz = tcom.fz - watx.fz);
mvimage(&sdist);
image.fx += sdist.fx;
image.fy += sdist.fy;
image.fz += sdist.fz;
r2 = sdist.fx*sdist.fx +  sdist.fy*sdist.fy +  sdist.fz*sdist.fz;
rOI = sqrt(r2);
/*
image.fx = -(sdist.fx = ri[0].fx - rj[0].fx);
image.fy = -(sdist.fy = ri[0].fy - rj[0].fy);
image.fz = -(sdist.fz = ri[0].fz - rj[0].fz);
mvimage(&sdist);
image.fx += sdist.fx;
image.fy += sdist.fy;
image.fz += sdist.fz;
r2 = sdist.fx*sdist.fx +  sdist.fy*sdist.fy +  sdist.fz*sdist.fz;
rOI = sqrt(r2);
 */
//fixShel tethering forces
#define TethCon 1.0 
if(rOI > Rshel){
   V_teth += TethCon*(rOI-Rshel)*(rOI-Rshel);
   deriv = -2.0*TethCon*(rOI-Rshel);
   dx=(tcom.fx-watx.fx+image.fx)/rOI;
   dy=(tcom.fy-watx.fy+image.fy)/rOI;
   dz=(tcom.fz-watx.fz+image.fz)/rOI;
   for(i=0;i<3;i++){
      force[tw+i].fx += deriv*dx*mass[i]/WMass;
      force[tw+i].fy += deriv*dy*mass[i]/WMass;
      force[tw+i].fz += deriv*dz*mass[i]/WMass;
      force[i].fx -= deriv*dx*mass[i]/WMass;
      force[i].fy -= deriv*dy*mass[i]/WMass;
      force[i].fz -= deriv*dz*mass[i]/WMass;
   }
   /*
   fi[0].fx += deriv*dx;
   fi[0].fy += deriv*dy;
   fi[0].fz += deriv*dz;
   fj[0].fx -= deriv*dx;
   fj[0].fy -= deriv*dy;
   fj[0].fz -= deriv*dz;
   */
}

}
