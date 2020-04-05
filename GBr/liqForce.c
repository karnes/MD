#include	<md.h>
#include        <water.h>
#include	<system.h>

liqForce(pos, force)
    tripd	*pos;
    tripd	*force;
{
    int i,nw;
    double totMass;
    double com_zVel;

    if (natoms-nsolute == 0)
	return;
    tc++;
    for(i=0;i<3;i++){
	gGCX[i] = gGOX[i] = gGHX[i] = 0;
    } 
    // Glycerols
    VLIQ = INTRA_GLY = INTER_GLY = GLYC = 0.0;
    if(nGLY>0)
	GLYForce();
    VLIQ += INTER_GLY + INTRA_GLY;
    // THA
    INTRA_THA = INTER_THA = THAC = VkapWin = 0.0;
    if(nTHA>0)
	THAForce();
    VLIQ += INTER_THA + INTRA_THA + VkapWin;
    // THA-Glycerols
    GLYTHAV = GLYTHAC = 0.0;
    if(nGLY>0 && nTHA>0)
	GLYTHA();
    VLIQ += GLYTHAV;
    // halogen solutes
    Xz = cosXz = cosTHAX = rXN = XGLYV = XGLYC = XTHAV = XTHAC = 0.0;
    for(i=0;i<4;i++){
	XGLYs[i] = 0.0;
    }
    if(nBr==1){
	BrForce();
    }
    else if(nCl2==1){
	Cl2Force();
	VLIQ+=Cl2bondE;
    } 
    else if(nTS==1){
	TSForce();
	VLIQ+=TSbondE+TSbendE;
    }  
    VLIQ += XGLYV + XTHAV;
    VWATTS = WATERV = 0.0;
    nw = natoms - nTHA*THAsites - nGLY*GLYsites - nBr - nCl2 - nTS;
    if(nw>0){
	waterForce(pos,force,nw);
	VLIQ += VWATTS + WATERV;
	if(nTHA==1){
	    waterTHA();
	    VLIQ += GLYTHAV;
	}
    }
//fprintf(stderr,"nw = %d, waterV = %f, VWatts = %f, waterTHA = %f\n",nw,KCAL*WATERV,KCAL*VWATTS,KCAL*GLYTHAV);
    getSysCom();
    // record forces along R_THA-X, if apropriate
    if(nTHA == 1 && nBr+nTS+nCl2 == 1)
	cmf();
    // FORCE INTEGRATION ON Br-
    vel[natoms-1].fz = 0.0;
    BrZForce = force[natoms-1].fz;
    force[natoms-1].fz = 0.0;
    totMass =com_zVel =  0.0;
    for(i=0;i<natoms-1;i++){
	totMass+=mass[i];
	com_zVel+=vel[i].fz*mass[i];
    }
    com_zVel/=totMass;
    for(i=0;i<natoms-1;i++){
	vel[i].fz-=com_zVel;
    }

}

int getSysCom(){
    int i,j;
    double sysMass = 0.0;
    syscomz = 0.0;
    for(i=0;i<natoms;i++){
	sysMass += mass[i];
	syscomz += pos[i].fz*mass[i];
    }
    syscomz/=sysMass;
}
