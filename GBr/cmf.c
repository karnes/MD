#include	<md.h>
#include	<system.h>
#include	<math.h>

cmf()
{

    // calculates the force on THA and X alonf R_THA-X
    //
    int i,j,k,nX;
    tripd Rn,Rvec,Xcom,Tcom;
    double Xtotmass,Ttotmass;

    if(nBr==1){
	nX = 1;
    }
    else if(nCl2==1){
	nX = 2;
    }
    else if(nTS==1){
	nX = 3;
    }

    /***	find THA CoM and X CoM   ***/

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
	k = nGLY*GLYsites+nTHA*THAsites + i;
	Xcom.fx+=pos[k].fx*mass[k];
	Xcom.fy+=pos[k].fy*mass[k];
	Xcom.fz+=pos[k].fz*mass[k];
	Xtotmass+=mass[k];
    }
    Xcom.fx/=Xtotmass;
    Xcom.fy/=Xtotmass;
    Xcom.fz/=Xtotmass;   
    // R is X pointing to THA
    Rvec.fx = Tcom.fx-Xcom.fx;
    Rvec.fy = Tcom.fy-Xcom.fy;
    Rvec.fz = Tcom.fz-Xcom.fz;

    rXCom = sqrt(sq(Rvec.fx)+sq(Rvec.fy)+sq(Rvec.fz));
    // Rn is R unit vector
    Rn.fx = Rvec.fx/rXCom;
    Rn.fy = Rvec.fy/rXCom;
    Rn.fz = Rvec.fz/rXCom;

    // calculate net THA and X forces projected on R //
    fTR = fXR = 0.0;
    for(i=0;i<THAsites;i++){
	k = nGLY*GLYsites + i;
	fTR += force[k].fx*Rn.fx;
	fTR += force[k].fy*Rn.fy;
	fTR += force[k].fz*Rn.fz;
    }
    for(i=0;i<nX;i++){
	k = nGLY*GLYsites + nTHA*THAsites + i;
	fXR += force[k].fx*Rn.fx;
	fXR += force[k].fy*Rn.fy;
	fXR += force[k].fz*Rn.fz;
    }
    //fprintf(stderr,"fXR = %f, fTR = %f. fXR + fTR = %f, Ucon = %f\n",fXRt*KCAL,fTRt*KCAL,fXRt*KCAL+fTRt*KCAL,Ucon*KCAL);

}
