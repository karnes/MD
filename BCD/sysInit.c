#include	<md.h>
#include	<system.h>
#include	<math.h>
#include        <water.h>

typedef struct lj {
	double S, E, Q;
	} LJ;

sysInit(confile)
	char	*confile;
{
	FILE	*fp;
	char	sbuf[256];
	int	i, j;
        int p1, p2;
	double 	sig,eps,q,s6,kr,kst,eqb,kb;
	double ptor[4] = {0.0};
	double f0,f1,f2,f3,f4;
	double BCDq[21]={0.0};
	double BCDsig[21]={0.0};
	double BCDeps[21]={0.0};
	LJ lj[BrOs+3];
	double swmax, swmin;
     	double dz, dz3, dz4, dz5, bc_factor;
	int chSw;

	if ((fp = fopen("/home/jkarnes/MD/BCD/run/inpfiles/LJq15pairs.txt","r")) == NULL) {
		fprintf(stderr, "sysinit: cannot open b-CD non-bonded pairs\n");
		exit(1);
	}
	fgets(sbuf,256,fp);// not using first line
	fgets(sbuf,256,fp);// not using second line
	for(i=0;i<9849;i++){
		fgets(sbuf,256,fp);
		if(sscanf(sbuf,"%d%d",&p1,&p2) != 2) {
		   fprintf(stderr, "sysInit: error reading BCD pairs\n");
		   exit(1);
		}
		pairNB[i][0]=p1;
		pairNB[i][1]=p2;
	}
//fprintf(stderr,"first pairNB: %d, %d. last pairNB: %d, %d.\n",pairNB[0][0],pairNB[0][1],pairNB[9848][0],pairNB[9848][1]);

	close(fp);

	if ((fp = fopen(confile,"r")) == NULL) {
		fprintf(stderr, "sysinit: cannot open %s\n", confile);
		exit(1);
	}

/* Begin reading confile */
	fgets(sbuf,256,fp);// first line
	if (sscanf(sbuf,"%d%d%s",&nBrO, &nBCD, waterP) != 3) {
		fprintf(stderr, "sysInit: error reading line 1\n");
		exit(1);
	}
	else{
	   if (sscanf(sbuf,"%d%d%s%lf%lf",&nBrO, &nBCD, waterP,&zGH,&wGH) != 5) {
		zGH = 0.0;
	        wGH = 999.9;
	   }
	   else{
	      fprintf(stderr,"sysInit: Host-guest com dist along BCDz = %f, window width = %f\n",zGH,wGH);
	   }
	}
//	fprintf(stderr,"sysInit check: Host-guest r_com = %f, window width = %f\n",zGH,wGH);
//	nsolute = 0; // if using evb, will redefine later
//fprintf(stderr,"nBrO=%d, nBCD=%d, waterP = %s\n",nBrO, nBCD, waterP);	
/* Read in Br-octane parameters */
	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading Br-octane Br LJq\n");
		exit(1);
	}
	lj[4].S = sig;
        lj[4].E = eps;
        lj[4].Q = q;
//fprintf(stderr, "sysInit: sig:%f=sig, %f=eps, %f=q\n",sig,eps,q);
	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading Br-octane CH2(Br) LJq\n");
		exit(1);
	}
	lj[3].S = sig;
        lj[3].E = eps;
        lj[3].Q = q;
//fprintf(stderr, "sysInit: sig:%f=sig, %f=eps, %f=q\n",sig,eps,q);
	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading Br-octane CH2 LJq\n");
		exit(1);
	}
	for(i=5;i<11;i++){
 	   lj[i].S = sig;
           lj[i].E = eps;
           lj[i].Q = q;
	}
//fprintf(stderr, "sysInit: sig:%f=sig, %f=eps, %f=q\n",sig,eps,q);
	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading Br-octane CH3 LJq\n");
		exit(1);
	}
	lj[11].S = sig;
        lj[11].E = eps;
        lj[11].Q = q;
	/* read Br-oct stretches */
//fprintf(stderr, "sysInit: sig:%f=sig, %f=eps, %f=q\n",sig,eps,q);
	for(i=0;i<4;i++){
		fgets(sbuf,256,fp);
		if(sscanf(sbuf,"%lf%lf",&kr,&kst) != 2){
		   fprintf(stderr, "sysInit: error reading Br-oct bond stretches\n");
		   exit(1);
		}
		BeqBond[i]=kr;
		BkStr[i]=kst / KCAL;
//fprintf(stderr, "sysInit: BrO str[%d]:%4.3f=r, %7.3f=kStr\n",i,BeqBond[i],BkStr[i]*KCAL);
	}
	for(i=0;i<4;i++){
		fgets(sbuf,256,fp);
		if(sscanf(sbuf,"%lf%lf",&eqb,&kb) != 2){
		   fprintf(stderr, "sysInit: error reading Br-oct bond bends\n");
		   exit(1);
		}
		BeqBend[i]=eqb*(PI/180.0);
		BkBend[i]=kb / KCAL;
//fprintf(stderr, "sysInit: BrO bend[%d]:%5.3f=ang, %7.3f=kBend\n",i,BeqBend[i]/(PI/180),BkBend[i]*KCAL);
	}
	for(i=0;i<4;i++){
		fgets(sbuf,256,fp);
		if(sscanf(sbuf,"%lf%lf%lf%lf",&ptor[0],&ptor[1],&ptor[2]) != 3){
		   fprintf(stderr, "sysInit: error reading BrO torsions\n");
		   exit(1);
		}
		for(j=0;j<3;j++){
		   BpTors[i][j]=ptor[j]/KCAL;
		}
	}
//		if(i>=0){ /* overwrite with con_hexC torsions */
//		   BpTors[i][0]=0.7055/KCAL;
//		   BpTors[i][1]=-0.1355/KCAL;
//		   BpTors[i][2]=1.50725/KCAL;	
//		   BpTors[i][3]=0.0;
//		}
	for(i=0;i<4;i++){
//	   fprintf(stderr, "sysInit: BpTors[%d]:%5.3f = 0, %5.3f = 1, %5.3f= 2 \n",i,BpTors[i][0]*KCAL,BpTors[i][1]*KCAL,BpTors[i][2]*KCAL);
	}
/* Read in water parameters */
	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading water O LJq\n");
		exit(1);
	}
	lj[0].S = sig;
        lj[0].E = eps;
        lj[0].Q = q;
//	OQ=q;

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading water H LJq\n");
		exit(1);
	}
	lj[1].S = lj[2].S = sig;
        lj[1].E = lj[2].E = eps;
        lj[1].Q = lj[2].Q = q;
//	HQ=q;
//fprintf(stderr,"water H LJq: %f, %f, %f\n",sig, eps, q);	
/* Read in BCD parameters */
	for(i=0;i<21;i++){
		fgets(sbuf,256,fp);
		if(sscanf(sbuf,"%lf%lf%lf",&sig,&eps,&q) != 3) {
		   fprintf(stderr, "sysInit: error reading BCD charges\n");
		   exit(1);
		}
		BCDq[i] = q;
		BCDsig[i]=sig;
		BCDeps[i]=eps;
//fprintf(stderr,"BCD LJq: %d:  %f, %f, %f\n",i,sig, eps, q);	
	}
	/* make mixed LJq terms */
	for(i=0;i<21;i++){
	   /* make intra BCD mixed terms */	
	   for(j=0;j<21;j++){
		eps = sqrt(BCDeps[i]*BCDeps[j]);
		sig = (BCDsig[i]+BCDsig[j])/2.0;
		s6 = sig*sig;
		s6 = s6*s6*s6;
		BCDlj[i][j].b = 4 * eps * s6 / KCAL;
		BCDlj[i][j].a = BCDlj[i][j].b * s6;
		BCDlj[i][j].q = BCDq[i] * BCDq[j] / E2;
//fprintf(stderr,"BCDlj[%d][%d] sig = %f, eps = %f\n",i,j,sig,eps);
//fprintf(stderr,"BCDlj[%d][%d].b = %f, a = %f, q = %f\n",i,j,BCDlj[i][j].b*KCAL, BCDlj[i][j].a*KCAL,BCDlj[i][j].q*E2);
	   }
	   /* make water - BCD mixed terms */
	   for(j=0;j<3;j++){
		eps = sqrt(BCDeps[i]*lj[j].E);
		sig = (BCDsig[i]+lj[j].S)/2.0;
		s6 = sig*sig;
		s6 = s6*s6*s6;
		wBCDlj[i][j].b = 4 * eps * s6 / KCAL;
		wBCDlj[i][j].a = wBCDlj[i][j].b * s6;
		wBCDlj[i][j].q = BCDq[i] * lj[j].Q / E2;
//fprintf(stderr,"wBCDlj[%d][%d] sig = %f, eps = %f\n",i,j,sig,eps);
//if(i==0&&j==0)fprintf(stderr,"wBCDlj[%d][%d].b = %f, a = %f, q = %f\n",i,j,wBCDlj[i][j].b*KCAL, wBCDlj[i][j].a*KCAL,wBCDlj[i][j].q*E2);
	   }
	   /* make Br-octane - BCD mixed terms */
	   for(j=0;j<BrOs;j++){
		eps = sqrt(BCDeps[i]*lj[j+3].E);
		sig = (BCDsig[i]+lj[j+3].S)/2.0;
		s6 = sig*sig;
		s6 = s6*s6*s6;
		bBCDlj[i][j].b = 4 * eps * s6 / KCAL;
		bBCDlj[i][j].a = bBCDlj[i][j].b * s6;
		bBCDlj[i][j].q = BCDq[i] * lj[j+3].Q / E2;
	   }
	}
//i=0;j=0;
//fprintf(stderr,"wBCDlj[%d][%d].b = %f, a = %f, q = %f\n",i,j,wBCDlj[i][j].b*KCAL, wBCDlj[i][j].a*KCAL,wBCDlj[i][j].q*E2);
	for(i=0;i<BrOs+3;i++){
	   for(j=0;j<BrOs+3;j++){
		eps = sqrt(lj[i].E*lj[j].E);
		sig = (lj[i].S + lj[j].S)/2.0;
		s6 = sig*sig;
		s6 = s6*s6*s6;
		slj[i][j].b = 4 * eps *s6 / KCAL;
		slj[i][j].a = slj[i][j].b * s6;
		slj[i][j].q = lj[i].Q * lj[j].Q / E2;
	   }
	}
//fprintf(stderr,"LJq:BCD7,BCD0:b = %f, a = %f, q = %f\n",BCDlj[7][0].b, BCDlj[7][0].a, BCDlj[7][0].q);
//fprintf(stderr,"LJq:BCD7,H2O0:b = %f, a = %f, q = %f\n",wBCDlj[7][0].b, wBCDlj[7][0].a, wBCDlj[7][0].q);
//i=0;j=0;
//fprintf(stderr,"wBCDlj[%d][%d].b = %f, a = %f, q = %f\n",i,j,wBCDlj[i][j].b*KCAL, wBCDlj[i][j].a*KCAL,wBCDlj[i][j].q*E2);
	/* continue reading confile */
	for(i=0;i<6;i++){
		fgets(sbuf,256,fp);
		if(sscanf(sbuf,"%lf%lf",&kr,&kst) != 2){
		   fprintf(stderr, "sysInit: error reading BCD bond stretches\n");
		   exit(1);
		}
		eqBond[i]=kr;
		kStr[i]=kst / KCAL;
//fprintf(stderr, "sysInit: str[%d]:%4.3f=r, %7.3f=kStr\n",i,eqBond[i],kStr[i]*KCAL);
	}
	for(i=0;i<11;i++){
		fgets(sbuf,256,fp);
		if(sscanf(sbuf,"%lf%lf",&eqb,&kb) != 2){
		   fprintf(stderr, "sysInit: error reading BCD bond bends\n");
		   exit(1);
		}
		eqBend[i]=eqb*(PI/180.0);
		kBend[i]=kb / KCAL;
//fprintf(stderr, "sysInit: bend[%d]:%5.3f=ang, %7.3f=kBend\n",i,eqBend[i]/(PI/180),kBend[i]*KCAL);
	}
	for(i=0;i<9;i++){
		fgets(sbuf,256,fp);
		if(sscanf(sbuf,"%lf%lf%lf",&ptor[0],&ptor[1],&ptor[2]) != 3){
		   fprintf(stderr, "sysInit: error reading BCD torsions\n");
		   exit(1);
		}
		/* convert Ryckaert-Bellemans to OPLS */
	/*	f4 = 0.0;
		f3 = ptor[3] / (-2.0);
		f2 = (ptor[2] - (4.0*f4)) * -1.0;
		f1 = ((2.0*ptor[1])-(3.0*f3)) * -1.0;
	*/
		/* assign and convert kJ to kcal */
	/*	pTors[i][0]=f1/(2.0*4.184*KCAL);
		pTors[i][1]=f2/(2.0*4.184*KCAL);
		pTors[i][2]=f3/(2.0*4.184*KCAL);	
	*/	/* convert Ryckaert-Bellemans to OPLS */
/*		f4 = 0.0;
		f3 = ptor[3] / (-2.0);
		f2 = ptor[2] / (-1.0);
		f1 = ((2.0*ptor[1])-(3.0*f3))*-1.0;
		f0 = 0;	
*/
		/* assign and convert kcal */
		pTors[i][0]=ptor[0]/KCAL;
		pTors[i][1]=ptor[1]/KCAL;
		pTors[i][2]=ptor[2]/KCAL;	
fprintf(stderr, "sysInit: BCD torsions (Vn/2): pTors[%d]:%5.3f = 0, %5.3f = 1, %5.3f= 2 \n",i,pTors[i][0]*(KCAL),pTors[i][1]*(KCAL),pTors[i][2]*(KCAL));
//		pTors[i][0]=0.0;
//		pTors[i][1]=0.0;
//		pTors[i][2]=0.0;	
	}
fprintf(stderr, "sysInit: nsolute = %d\n",nsolute);
	/* get the 63 1-4 L-J pairs */
//i=0;j=0;
//fprintf(stderr,"wBCDlj[%d][%d].b = %f, a = %f, q = %f\n",i,j,wBCDlj[i][j].b*KCAL, wBCDlj[i][j].a*KCAL,wBCDlj[i][j].q*E2);
	for(i=0;i<63;i++){
		fgets(sbuf,256,fp);
		if(sscanf(sbuf,"%d%d",&p1,&p2) != 2){
		   fprintf(stderr, "sysInit: error reading BCD 1-4 LJ pairs\n");
		   exit(1);
		}
		pair[i][0] = p1;
		pair[i][1] = p2;
	}
//	fprintf(stderr, "sysinit: pair[62][0] = %d, pair[62][1] = %d\n",pair[62][0],pair[62][1]);
//i=0;j=0;
//fprintf(stderr,"wBCDlj[%d][%d].b = %f, a = %f, q = %f\n",i,j,wBCDlj[i][j].b*KCAL, wBCDlj[i][j].a*KCAL,wBCDlj[i][j].q*E2);
	waterInit();
	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%d%lf%lf%lf%lf",
		&chSw,&swmax,&swmin,&swSoMax,&swSoMin) != 5) {
		fprintf(stderr, "sysinit: error reading switching parameters. set to default.\n");
		chSw = 0;
	//	exit(1);
	}
//	fprintf(stderr, "sysinit: chSw = %d, swmax = %f, swmin = %f\n",chSw,swmax,swmin);
	if (chSw==1){/* Change switching parameters */
	        fprintf(stderr, "sysinit: manual switching parameters: swmax = %f, swmin = %f, swSoMax = %f, swSoMin = %f\n",swmax,swmin,swSoMin,swSoMax);
		if (pbcType[0] == 'C'){
			bc_factor = 1.; /* for cubic */
		}
		else if (pbcType[0] == 'O'){
			bc_factor = 0.75; /* for truncated octahedron */
		}
		else
			ERROR((stderr,"sysinit:periodic boundaries not defined\n"), exit);

        	if (xwall > ywall || xwall > zwall)
			ERROR((stderr,"sysinit: xwall must be smaller than y or z\n"), exit);
		swr2max = bc_factor * xwall * xwall;
		/* This is the square of 1/2 the distance between the hexagonal
		 * faces of the pto box, or the sq of 1/2 the small length of
		 * the cube */
	 
		swr2min = bc_factor * (xwall - 1.0) * (xwall - 1.0);
		if (swmax*swmax > swr2max)	
			ERROR((stderr,"sysinit: swmax larger than maximum\n"), exit);
		if (swmin*swmin > swr2min)	
			ERROR((stderr,"sysinit: swmin larger than maximum\n"), exit);
		swSoMax = swSoMax*swSoMax;
		swSoMin = swSoMin*swSoMin;
		if (swSoMax > swr2max)	
			ERROR((stderr,"sysinit: swSoMax larger than maximum\n"), exit);
		if (swSoMin > swr2min)	
			ERROR((stderr,"sysinit: swSoMin larger than maximum\n"), exit);
		swr2max = swmax*swmax;
		swr2min = swmin*swmin;
		dz = swr2max - swr2min; /*r**2 (min max) for switching region */
		dz3 = (dz*dz*dz);
		dz4 = dz3*dz;
		dz5 = dz4*dz;
		swcoef[0] =  1.0; /* switching parameters for group cut-offs */
		swcoef[1] =  -10.0 / dz3;
		swcoef[2] =   15.0 / dz4;
		swcoef[3] =   -6.0 / dz5;
		pswcoef[0] = -60.0 / dz3;	/* 2 times swcoef's derivs */
		pswcoef[1] = 120.0 / dz4;
		pswcoef[2] = -60.0 / dz5;
		/* for solute-solvent*/
		dz = swSoMax - swSoMin;/* r**2 (min max) for switching region */
		dz3 = (dz*dz*dz);
		dz4 = dz3*dz;
		dz5 = dz4*dz;
		swSocoef[0] =1.0; /* switching parameters for solvent-solute*/
		swSocoef[1] =  -10.0 / dz3;
		swSocoef[2] =   15.0 / dz4;
		swSocoef[3] =   -6.0 / dz5;
		pswSocoef[0] = -60.0 / dz3;	/* 2 times swSocoef's derivs */
		pswSocoef[1] = 120.0 / dz4;
		pswSocoef[2] = -60.0 / dz5;

	}
	else {/*set solute-solvent to defaut values*/
	        fprintf(stderr, "sysinit: use default switching parameters\n");
		swSoMax = swr2max;
		swSoMin = swr2min;
		swSocoef[0] = swcoef[0];
		swSocoef[1] = swcoef[1];
		swSocoef[2] = swcoef[2];
		swSocoef[3] = swcoef[3];
		pswSocoef[0] = pswcoef[0];
		pswSocoef[1] = pswcoef[1];
		pswSocoef[2] = pswcoef[2];
	}
	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%d",
		&KillFinger) != 1) {
		KillFinger = 0;
		fprintf(stderr, "sysinit: error reading KillFinger flag. KillFinger off. Guest Biasing off.\n");
	}
	else if (sscanf(sbuf,"%d%d%lf%lf%lf",&KillFinger,&guestBiasOn,&gbA,&gbB,&gbO) != 5) {
		fprintf(stderr, "sysinit: error reading Guest Biasing. Guest Biasing off.\n");
		guestBiasOn = 0;
		gbA = gbB = gbO = 0.0;
	}
	fprintf(stderr,"sysinit: KillFinger = %d, guestBiasing = %d.\n",KillFinger,guestBiasOn);
	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%d%lf%lf%lf%lf",
		&windowOn,&center_w,&width_w,&power_w,&pot_w) != 5) {
		windowOn = 0;
		center_w = width_w = power_w = pot_w = 0.0;
		fprintf(stderr, "sysinit: error reading BCD z-windowing.\n");
	}
	fprintf(stderr,"sysinit: windowOn = %d, center = %4.2f, width = %4.2f, pow = %3.1f, pot = %3.2f.\n",windowOn,center_w,width_w,power_w,pot_w);

	int	jns;
	double	sigS,epsS,qS;
/*read EVB potential parameters*/
    if(!fgets(sbuf, 256, fp)){
		fprintf(stderr, "sysInit: no EVB initialization line. Skipping EVB initialization...\n");
        	evbFlag = 0; nsolute = 0,solBCD_c=999.9,solBCD_w=0.0;
    }
    else if(sscanf(sbuf, "%d%d%lf%lf", &evbFlag, &nsolute, &solBCD_c, &solBCD_w) != 4){
		fprintf(stderr, "sysInit: error reading EVB initialization line. Skipping EVB initialization...\n");
        	evbFlag = 0; nsolute = 0,solBCD_c=999.9,solBCD_w=0.0,parabA=0.0,parabB=0.0;
    }
    if(sscanf(sbuf, "%d%d%lf%lf%lf%lf", &evbFlag, &nsolute, &solBCD_c, &solBCD_w, &parabA, &parabB) != 6){
		fprintf(stderr, "sysInit: no SN2 - BCD biasing along p-vector. setting to 0.0...\n");
        	parabA=0.0,parabB=0.0;
    }
    fprintf(stderr,"sysinit: test: evbFlag = %d, nsolute = %d, solBCD_c = %5.3f, solBCD_w = %5.3f, parabA = %5.2f, parabB = %5.2f\n",evbFlag,nsolute,solBCD_c,solBCD_w,parabA,parabB);
    parabA /= KCAL;
    if(evbFlag==1){
        fprintf(stderr, "sysInit: Begin EVB initialization. nsolute = %d...",nsolute);
	/* allocate space for EVB data */
	evb_space();
//	fprintf(stderr,"done with evb_space\n");
	/* start reading in EVB parameters */
	fgets(sbuf, 256, fp);
	if (sscanf(sbuf, "%lf %lf %lf %lf", &Qa, &Qb, &Qc, &Qd) != 4){
		fprintf(stderr, "sysInit: error reading EVB line 2\n");
        	exit(1);
	}
	Qd /= KCAL;
	fgets(sbuf, 256, fp);
	if (sscanf(sbuf, "%lf %lf %lf", &Dmors, &amors, &reqmors) != 3){
		fprintf(stderr, "sysInit: error reading EVB line 3\n");
        	exit(1);
	}
	Dmors /= KCAL;
	fgets(sbuf, 256, fp);
	if (sscanf(sbuf, "%lf %lf %lf", &Zion, &rstar, &EA) != 3){
		fprintf(stderr, "sysInit: error reading EVB line 4\n");
        	exit(1);
	}
	Zion /= KCAL;
	EA /= KCAL;
	fgets(sbuf, 132, fp);
	if (sscanf(sbuf, "%lf %lf %lf %lf", &bid, &epsid, &sigid, &nid) != 4){
		fprintf(stderr, "sysInit: error reading line 5\n");
        	exit(1);
	}
	epsid /= KCAL;
	fgets(sbuf, 132, fp);
	if (sscanf(sbuf, "%lf %lf %lf", &mua, &mub, &Qbeta) != 3){
		fprintf(stderr, "sysInit: error reading line 5\n");
        	exit(1);
	}
	Qbeta /= KCAL;

if	(nsolute > 0){
/*
 *	Read information about modified switching parameters
 */
//	fgets(sbuf,256,fp);
//	if (sscanf(sbuf,"%d%lf%lf%lf%lf",
//		&chSw,&swmax,&swmin,&swSoMax,&swSoMin) != 5) {
//		fprintf(stderr, "sysinit: error reading line 12\n");
//		exit(1);
//	}
//	if (chSw){/* Change switching parameters */
//		if (pbcType[0] == 'C'){
//			bc_factor = 1.; /* for cubic */
//		}
//		else if (pbcType[0] == 'O'){
//			bc_factor = 0.75; /* for truncated octahedron */
//		}
//		else
//			ERROR((stderr,"sysinit:periodic boundaries not defined\n"), exit);
//
  //      	if (xwall > ywall || xwall > zwall)
//			ERROR((stderr,"sysinit: xwall must be smaller than y or z\n"), exit);
//		swr2max = bc_factor * xwall * xwall;
		/* This is the square of 1/2 the distance between the hexagonal
		 * faces of the pto box, or the sq of 1/2 the small length of
		 * the cube */
	 
//		swr2min = bc_factor * (xwall - 1.0) * (xwall - 1.0);
//		if (swmax*swmax > swr2max)	
//			ERROR((stderr,"sysinit: swmax larger than maximum\n"), exit);
//		if (swmin*swmin > swr2min)	
//			ERROR((stderr,"sysinit: swmin larger than maximum\n"), exit);
//		swSoMax = swSoMax*swSoMax;
//		swSoMin = swSoMin*swSoMin;
//		if (swSoMax > swr2max)	
//			ERROR((stderr,"sysinit: swSoMax larger than maximum\n"), exit);
//		if (swSoMin > swr2min)	
//			ERROR((stderr,"sysinit: swSoMin larger than maximum\n"), exit);
//		swr2max = swmax*swmax;
//		swr2min = swmin*swmin;
//		dz = swr2max - swr2min; /*r**2 (min max) for switching region */
//		dz3 = (dz*dz*dz);
//		dz4 = dz3*dz;
//		dz5 = dz4*dz;
//		swcoef[0] =  1.0; /* switching parameters for group cut-offs */
//		swcoef[1] =  -10.0 / dz3;
//		swcoef[2] =   15.0 / dz4;
//		swcoef[3] =   -6.0 / dz5;
//		pswcoef[0] = -60.0 / dz3;	/* 2 times swcoef's derivs */
//		pswcoef[1] = 120.0 / dz4;
//		pswcoef[2] = -60.0 / dz5;
		/* for solute-solvent*/
//		dz = swSoMax - swSoMin;/* r**2 (min max) for switching region */
//		dz3 = (dz*dz*dz);
//		dz4 = dz3*dz;
//		dz5 = dz4*dz;
//		swSocoef[0] =1.0; /* switching parameters for solvent-solute*/
//		swSocoef[1] =  -10.0 / dz3;
//		swSocoef[2] =   15.0 / dz4;
//		swSocoef[3] =   -6.0 / dz5;
//		pswSocoef[0] = -60.0 / dz3;	/* 2 times swSocoef's derivs */
//		pswSocoef[1] = 120.0 / dz4;
//		pswSocoef[2] = -60.0 / dz5;
//
//	}
//	else {/*set solute-solvent to defaut values*/
//		swSoMax = swr2max;
//		swSoMin = swr2min;
//		swSocoef[0] = swcoef[0];
//		swSocoef[1] = swcoef[1];
//		swSocoef[2] = swcoef[2];
//		swSocoef[3] = swcoef[3];
//		pswSocoef[0] = pswcoef[0];
//		pswSocoef[1] = pswcoef[1];
//		pswSocoef[2] = pswcoef[2];
//	}
	/* Read solute lj parameters	*/

	for (jns =0;jns<nsolute;jns++){
		/*loop over all solute atoms*/
		fgets(sbuf,256,fp);
		if (sscanf(sbuf,"%lf%lf%lf",&sigS,&epsS,&qS) != 3) {
			fprintf(stderr, "sysInit: error reading solute LJ\n");
			exit(1);
		}
		for(i=0;i<3;i++){
		   sig = (lj[i].S+sigS)/2.;
		   eps = sqrt(lj[i].E*epsS);
		   s6 = sig*sig*sig;
		   s6 = s6*s6;
		   wslj[i][jns].b = 4*eps*s6/KCAL;
		   wslj[i][jns].a = wslj[i][jns].b * s6;
		   wslj[i][jns].q = qS*lj[i].Q/E2;
		}
		
		for(i=3;i<3+BrOs;i++){
		   sig = (lj[i].S+sigS)/2.;
		   eps = sqrt(lj[i].E*epsS);
		   s6 = sig*sig*sig;
		   s6 = s6*s6;
		   bslj[i-3][jns].b = 4*eps*s6/KCAL;
		   bslj[i-3][jns].a = bslj[i-3][jns].b * s6;
		   bslj[i-3][jns].q = qS*lj[i].Q/E2;
		}
		for(i=0;i<21;i++){
		   sig = (BCDsig[i]+sigS)/2.;
		   eps = sqrt(BCDeps[i]*epsS);
		   s6 = sig*sig*sig;
		   s6 = s6*s6;
		   cdslj[i][jns].b = 4*eps*s6/KCAL;
		   cdslj[i][jns].a = cdslj[i][jns].b * s6;
		   cdslj[i][jns].q = qS*BCDq[i]/E2;
		}
	}
	
        fgets(sbuf,256,fp);
        if (sscanf(sbuf,"%lf%lf%lf%lf",&Ecenter_w,&Ewidth_w,&Epower_w,&Epot_w) != 4) {
          fprintf(stderr, "sysinit: error reading EVB window potential parameters\n");
            exit(1);
        }
        fgets(sbuf,256,fp);
        if (sscanf(sbuf,"%lf%lf%lf",&x1RC, &x2RC, &kRC) != 3) {
          fprintf(stderr, "sysinit: error reading reaction coordinate constraints parameters \n");
            exit(1);
        }
        fgets(sbuf,256,fp);
        if (sscanf(sbuf,"%lf%lf%lf%lf",&Abias, &Bbias, &alpha_bias, &beta_bias) != 4) {
          fprintf(stderr, "sysinit: error reading PMF biasing pot parameters \n");
            exit(1);
        }
        fprintf(stderr, "sysinit: EVB biasing: A = %f, B = %f, alpha = %f, beta = %f\n",Abias,Bbias,alpha_bias,beta_bias);
	Abias /= KCAL;
	Bbias /= KCAL;

fprintf(stderr, "initialization complete.\n");
}/*endif nsolute > 0 */
}/*endif evbFlag == 1 */

/*	fgets(sbuf,256,fp);
        if (sscanf(sbuf,"%d %lf %d",&KillFinger,&ExCon,&ReflectFQ) != 3) {
             fprintf(stderr,"sysinit:error reading KillFinger,ExCon,ReflectFQ\n");

             exit(1);
        }
        if (KillFinger)
             fprintf(stderr,"Flat interface option on\n");
	     */
fclose(fp);
fprintf(stderr, "sysInit: nsolute = %d\n",nsolute);
}

evb_space()
{
static int	allocate_evb = 1;
int i;	
/*	Allocate space for the EVB forces arrays if first time through ...*/

	if(allocate_evb && ((fevb = (tripd  *) calloc(natoms,sizeof(tripd))) == NULL)) {
		ERROR((stderr,"evb_space: out of core\n"), exit);
	}
	for(i=0;i<natoms;i++){
	   fevb[i].fx = fevb[i].fy = fevb[i].fz = 0.0;
	}
//fprintf(stderr, "evbi space allocated.\n");
	allocate_evb = 0;	/* so we won't allocate space the next time */

}
