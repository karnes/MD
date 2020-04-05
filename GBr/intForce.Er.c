#include	<md.h>
#include	<system.h>
#include	<math.h>
#include        <water.h>

// R_THA-X windowing (and biasing) potentials
intForce()
{
    VINT = V_bias = V_w = 0.0;
    if(nEr!=1)
	return;
    if(windowOn==1)
	window();
}


//#define Gibbs_surf 0.0
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
    //    *** Add biasing potential for ion. Form of potential is *** 
    if(biasingOn>0){
	bias_com = osysCmz;
	zpos = pos[natoms-1].fz - BIAS_G + bias_com;
	if(biasingOn==1){ //hyperbolic tangent
	    //     *** V_bias = BIAS_W (1 + tanh(BIAS_A*x) )/2             ***
	    //     ***   where x = pos[ion].fz - BIAS_G + bias_com  ***
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
    if(KillFinger==0){
	zpos = pos[natoms-1].fz - sysCm.fz - center_w;
    }
    else{
	zpos = pos[natoms-1].fz - 9.8 - center_w;
    }
    zeta = fabs(zpos)-width_w/2.;
    if(zeta > 0.){// outside window, calculate energy 
	if (zeta > 1.5*width_w){
	    fprintf(stderr,"window: too much outside the window\n");
	    exit(1);
	}
	else { // 0 < zeta < width_w
	    V_w = pot_w*pow(zeta,power_w-1);
	    deriv = V_w*power_w*sgn(zpos);
	    V_w *= zeta;
	    force[natoms-1].fz -= deriv;
	    if(KillFinger==0){
		for (i=0;i<natoms;i++){
		    force[i].fz += (deriv * mass[i]/totMass);
		}
		VINT += V_w;
	    }
	    else{
		VINT += V_w;
	    }
	}
    }
}


