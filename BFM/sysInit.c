#include	<md.h>
#include	<system.h>
#include	<math.h>

sysInit(confile)
	char	*confile;
{
	FILE	*fp;
	char	sbuf[256];
	int	i, j;
	double 	sig,eps,s6,sqrt(),sin(), cos(), asin(),
		R_CCl,R_CH,gama,gamH,k_CCl,k_CH,k_ClCCl,k_ClCH,convF,
		grkCI,
	 	sigC1,epsC1,qC1,
	 	sigCl,epsCl,qCl,
	 	sigH,epsH,qH,
	 	sigI,epsI,Qgr,Qex;
//fprintf(stderr,"bin -0.2 = %f, 0.2 = %f, -0.1 = %f,0.1 = %f\n",round((-.2+contRad)/contBin),round((.2+contRad)/contBin),round((-0.1+contRad)/contBin),round((0.1+contRad)/contBin));
	if ((fp = fopen(confile,"r")) == NULL) {
		fprintf(stderr, "sysinit: cannot open %s\n", confile);
		exit(1);
	}
	nCFM = (natoms - nsolute)/5;/*number of Chloroform molecules*/

/*	Get the appropriate input paramters for the Chloroform intra- forces */

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf%lf",&R_CCl,&R_CH,&gamH,&gama) != 4){
		fprintf(stderr, "sysInit: error reading line 1\n");
		exit(1);
	}

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf%lf",&k_CCl,&k_CH,&k_ClCCl,&k_ClCH) != 4){
		fprintf(stderr, "sysInit: error reading line 2\n");
		exit(1);
	}

	/* Convert force constant from dyneX1e5/cm to internal units*/
//	convF = AVOGADRO*1.0e5/(ERG*1.0e16);/*0.06022045*/
	ks[1] = ks[2] = ks[3] = 2.0*k_CCl/KCAL;//*convF;
	ks[4] = 2.0*k_CH/KCAL;//*convF;
	kb[1][2] = kb[1][3] = kb[2][3] = 2.0*k_ClCCl/KCAL;//*R_CCl*R_CCl*convF;
	kb[1][4] = kb[2][4] = kb[3][4] = 2.0*k_ClCH/KCAL;//*R_CH*R_CH*convF;
	/* calculte the pair equilibrium distances matrix*/
	gamH = gamH*PI/180.0;//PI- asin(sin(gama*PI/360.0)*2.0/sqrt(3.0));/*The H-C-Cl angle*/
	Req[1] = Req[2] = Req[3] = R_CCl;
	Req[4] = R_CH;
	Beq[1][2] = Beq[1][3] = Beq[2][3] = gama*PI/180.0;
	Beq[1][4] = Beq[2][4] = Beq[3][4] = gamH;

/*	Get the appropriate input paramters for the Chloroform inter- forces */

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sigC1, &epsC1, &qC1) != 3) {
		fprintf(stderr, "sysInit: error reading line 3\n");
		exit(1);
	}
	s6 = sigC1*sigC1*sigC1;
	s6 = s6*s6;
	CFMlj[0][0].b = 4*epsC1*s6/KCAL;
	CFMlj[0][0].a = CFMlj[0][0].b * s6;
	CFMlj[0][0].q = qC1*qC1/E2;

//	QCFM[0] = qC1;/*need for dipole calculations*/

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sigCl, &epsCl, &qCl) != 3) {
		fprintf(stderr, "sysInit: error reading line 4\n");
		exit(1);
	}
	s6 = sigCl*sigCl*sigCl;
	s6 = s6*s6;
	CFMlj[1][1].b = CFMlj[1][2].b = CFMlj[1][3].b = 4*epsCl*s6/KCAL;
	CFMlj[2][1].b = CFMlj[2][2].b = CFMlj[2][3].b = 4*epsCl*s6/KCAL;
	CFMlj[3][1].b = CFMlj[3][2].b = CFMlj[3][3].b = 4*epsCl*s6/KCAL;
	CFMlj[1][1].a = CFMlj[1][2].a = CFMlj[1][3].a = CFMlj[1][1].b * s6;
	CFMlj[2][1].a = CFMlj[2][2].a = CFMlj[2][3].a = CFMlj[1][1].b * s6;
	CFMlj[3][1].a = CFMlj[3][2].a = CFMlj[3][3].a = CFMlj[1][1].b * s6;
	CFMlj[1][1].q = CFMlj[1][2].q = CFMlj[1][3].q = qCl*qCl/E2;
	CFMlj[2][1].q = CFMlj[2][2].q = CFMlj[2][3].q = qCl*qCl/E2;
	CFMlj[3][1].q = CFMlj[3][2].q = CFMlj[3][3].q = qCl*qCl/E2;

//	QCFM[1] = QCFM[2] = QCFM[3] = qCl;/*need for dipole calculations*/

	sig = (sigC1+sigCl)/2.;
	eps = sqrt(epsC1*epsCl);
	s6 = sig*sig*sig;
	s6 = s6*s6;
	CFMlj[0][1].b = CFMlj[1][0].b = 4*eps*s6/KCAL;
	CFMlj[0][2].b = CFMlj[2][0].b = 4*eps*s6/KCAL;
	CFMlj[0][3].b = CFMlj[3][0].b = 4*eps*s6/KCAL;
	CFMlj[0][1].a = CFMlj[1][0].a = CFMlj[0][1].b * s6;
	CFMlj[0][2].a = CFMlj[2][0].a = CFMlj[0][1].b * s6;
	CFMlj[0][3].a = CFMlj[3][0].a = CFMlj[0][1].b * s6;
	CFMlj[0][1].q = CFMlj[1][0].q = qC1*qCl/E2;
	CFMlj[0][2].q = CFMlj[2][0].q = qC1*qCl/E2;
	CFMlj[0][3].q = CFMlj[3][0].q = qC1*qCl/E2;

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sigH, &epsH, &qH) != 3) {
		fprintf(stderr, "sysInit: error reading line 5\n");
		exit(1);
	}

	s6 = sigH*sigH*sigH;
	s6 = s6*s6;
	CFMlj[4][4].b = 4*epsH*s6/KCAL;
	CFMlj[4][4].a = CFMlj[4][4].b * s6;
	CFMlj[4][4].q = qH*qH/E2;

//	QCFM[4] = qH;/*need for dipole calculations*/

// this is the switch for CFM / CTC identity-dependent parameters	
// Hard-code QCFM for dipole calculations.
QCFM[0] = 0.32;
QCFM[1] = QCFM[2] = QCFM[3] = -0.14;
QCFM[4] = 0.10;

//	if(QCFM[4]==QCFM[3]){ // CTC
	// Switch density and cutoffs on sigH
	if(sigH == 3.5){ //then CTC
	   CFMden = 0.006225;
//	   QCFM[0] = 0.32;
//	   QCFM[4] = 0.10;
	   CHmax2 = 5.0*5.0;
	   CHmin2 = 2.0*2.0;
	   rCCmax = 7.9;
	   rCCmin = 3.0;
           CFMMass = 153.82;
	   fprintf(stderr,"sysInit: CCl4 run detected. Density = %f, rCH max = %f\n", CFMden,sqrt(CHmax2));
	}
	else if(sigCl == 3.4){  // then a CFM
	   CFMden = 0.007466;
	   CHmax2 = 5.0*5.0;
	   CHmin2 = 2.0*2.0;
	   rCCmax = 7.6;
	   rCCmin = 3.0;
           CFMMass = 119.38;
	   fprintf(stderr,"sysInit: CHCl3 run detected. Density = %f, rCH max = %f\n", CFMden,sqrt(CHmax2));
	}
	else if(sigCl == 3.47){  // then a BFM
	   CFMden = 0.00688685;
	   CHmax2 = 5.0*5.0;
	   CHmin2 = 2.0*2.0;
	   rCCmax = 7.6;
	   rCCmin = 3.0;
           CFMMass = 252.731;
	   fprintf(stderr,"sysInit: CHBr3 run detected. Density = %f, rCH max = %f\n", CFMden,sqrt(CHmax2));
	}
	else{ // unknown thing...
	   fprintf(stderr,"ssyInit.c: Problem with input file. unidentified sigCl. sigCl = %f\n",sigH);
	   exit(0);
	}
	

	sig = (sigC1+sigH)/2.;
	eps = sqrt(epsC1*epsH);
	s6 = sig*sig*sig;
	s6 = s6*s6;
	CFMlj[0][4].b = CFMlj[4][0].b = 4*eps*s6/KCAL;
	CFMlj[0][4].a = CFMlj[4][0].a = CFMlj[0][4].b * s6;
	CFMlj[0][4].q = CFMlj[4][0].q = qC1*qH/E2;

	sig = (sigCl+sigH)/2.;
	eps = sqrt(epsCl*epsH);
	s6 = sig*sig*sig;
	s6 = s6*s6;
	CFMlj[1][4].b = CFMlj[4][1].b = 4*eps*s6/KCAL;
	CFMlj[2][4].b = CFMlj[4][2].b = 4*eps*s6/KCAL;
	CFMlj[3][4].b = CFMlj[4][3].b = 4*eps*s6/KCAL;
	CFMlj[1][4].a = CFMlj[4][1].a = CFMlj[1][4].b * s6;
	CFMlj[2][4].a = CFMlj[4][2].a = CFMlj[1][4].b * s6;
	CFMlj[3][4].a = CFMlj[4][3].a = CFMlj[1][4].b * s6;
	CFMlj[1][4].q = CFMlj[4][1].q = qCl*qH/E2;
	CFMlj[2][4].q = CFMlj[4][2].q = qCl*qH/E2;
	CFMlj[3][4].q = CFMlj[4][3].q = qCl*qH/E2;

fclose(fp);
}
