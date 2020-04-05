#include	<md.h>
#include	<system.h>

liqForce(pos, force)
	tripd	*pos;
	tripd	*force;
{
	tc++;
	int initMat(void);
	static int init = 1;
//	fprintf(stderr,"init = %d\n",init);
	if(tc==0){
	   if(init){
		init = initMat();
	   }
	   zeroMat();
	}
	INTRAV = CFMNB = 0.;
	CFMForce();
	VLIQ = INTRAV + CFMNB;
}

int initMat(void)
{
	int i,j,k;
	for(i=0;i<4;i++){
		if((stkLife[i] = (int **)calloc(nCFM,sizeof(int *))) == NULL)
			ERROR((stderr,"wxtr: out of core, polar stacking.\n"),exit);
	}
	for(i=0;i<4;i++){
	   for(j=0;j<nCFM;j++){
		if((stkLife[i][j] = (int *)calloc((int)(numDP_x_dataRate/TCFdt)+1,sizeof(int))) == NULL)
		     ERROR((stderr,"wxtr: out of core, polar stacking.\n"),exit);
	   }
	}
	if((dipVec = (tripd **)calloc(nCFM, sizeof(tripd *))) == NULL
	|| (c3Vec  = (tripd **)calloc(nCFM, sizeof(tripd *))) == NULL
	|| (chVec  = (tripd **)calloc(nCFM, sizeof(tripd *))) == NULL
	|| (pstack  = (int **)calloc(nCFM, sizeof(int *))) == NULL
	|| (pstack2  = (int **)calloc(nCFM, sizeof(int *))) == NULL)
		ERROR((stderr,"liqForce.c: out of core, TCF vector initializations.\n"),exit);
	
	for(i=0;i< nCFM;i++){
		if((dipVec[i] = (tripd *)calloc((int)(numDP_x_dataRate / TCFdt)+1, sizeof(tripd))) == NULL
		|| (c3Vec[i]  = (tripd *)calloc((int)(numDP_x_dataRate / TCFdt)+1, sizeof(tripd))) == NULL
		|| (chVec[i]  = (tripd *)calloc((int)(numDP_x_dataRate / TCFdt)+1, sizeof(tripd))) == NULL
		|| (pstack[i]  = (int *)calloc((int)(numDP_x_dataRate / TCFdt)+1, sizeof(int))) == NULL
		|| (pstack2[i]  = (int *)calloc((int)(numDP_x_dataRate / TCFdt)+1, sizeof(int))) == NULL)
		  ERROR((stderr,"liqForce.c: out of core, TCF vector initializations.\n"),exit);
	}
	return 0;
}

void zeroMat(void)
{
	int i,j,k;
	avgDipDip = avgDipDip2 = avgDipDipn = avgDipDip2n = totInHCC = 0.0;
	planCos = fabs(cos((planAng+90.0)*PI/180.0));
//	CHPolarCount = 0;
//	CClPolarCount = 0;
//	for(i=0;i<5;i++){
//	   dipPolarCount[i] = 0;
//	   for(j=0;j<MAXPOLAR;j++){
//	      dipPolar[i][j][0] = dipPolar[i][j][1] = 0.0;
//	      CHPolar[j][0] = CHPolar[j][1] = 0.0;
//	      CClPolar[j][0] = CClPolar[j][1] = 0.0;
//	   }
//	   for(j=0;j<360;j++){
//	      cont1D[i][j] = 0;
//	   }
//	   for(j=0;j<80;j++){
//	      cont1Dcos[i][j] = 0;
//	      for(k=0;k<20;k++){
//		contCos[k][i][j] = 0;
//	      }
//	   }
//	   for(j=0;j<20;j++){
//	      contCosn[j][i] = 0;
//	   }
//	   cont1Dn[i] = 0;
//	}
	for(i=0;i<(2*(int)(contRad/contBin)+1);i++){
	   for(j=0;j<(2*(int)(contRad/contBin)+1);j++){
	      CCmap[i][j] = 0;
	      CHmap[i][j] = 0;
	      CClmap[i][j] = 0;
	      for(k=0;k<5;k++){
		 dipMap[k][i][j]=0.0;
		 dipMapsm[k][i][j]=0.0;
	      }
	   }
	}
	for(i=0;i<nCFM;i++){
		for(j=0;j<(int)(numDP_x_dataRate / TCFdt)+1;j++){
			dipVec[i][j].fx = dipVec[i][j].fy = dipVec[i][j].fz = 0.0;
			c3Vec[i][j].fx = c3Vec[i][j].fy = c3Vec[i][j].fz = 0.0;
			chVec[i][j].fx = chVec[i][j].fy = chVec[i][j].fz = 0.0;
			pstack[i][j] = -1;
			pstack2[i][j] = -1;
		}
	}
	for(i=0;i<(int)(gKrad/gKbin);i++){
	   gK[i] = 0.0;
	   gKn[i] = 0.0;
	}

//	fprintf(stderr,"init'd TCF in liqforce.c to zeroes \n");
	return;
}
