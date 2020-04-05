#include	<md.h>
#include	<system.h>
#include	<water.h>

typedef struct lj {
    double S, E, Q;
} LJ;
sysInit(confile)
    char	*confile;
{
    FILE	*fp;
    char	sbuf[256];
    int i,j,k;
    double 	sig,eps,q,s6;
    double kstr,req,kbend,theq,t0,t1,t2,f14;
    double c1sig,c1eps,c1q; //2 unique C terms
    double c2sig,c2eps,c2q; //2 unique C terms
    double osig,oeps,oq; 
    double esig,eeps,eq; 
    double hsig,heps,hq; 
    double swmax, swmin;
    double dz, dz3, dz4, dz5, bc_factor;
    int chSw;
    double dsig[DDCsites];
    double deps[DDCsites];

    if ((fp = fopen(confile,"r")) == NULL) {
	fprintf(stderr, "sysinit: cannot open %s\n", confile);
	exit(1);
    }

    fgets(sbuf,256,fp);
    if (sscanf(sbuf,"%d%d%s", &nDDC,&nEr,waterP) != 3) {
	fprintf(stderr, "sysInit: error reading line 1\n");
	exit(1);
    }
    // DDC LJq
    fgets(sbuf,256,fp);
    if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
	fprintf(stderr, "sysInit: error reading CH3 LJq\n");
	exit(1);
    }
    c1sig = sig;
    c1eps = eps;
    c1q = q;
    fgets(sbuf,256,fp);
    if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
	fprintf(stderr, "sysInit: error reading CH2 LJq\n");
	exit(1);
    }
    c2sig = sig;
    c2eps = eps;
    c2q = q;
    // DDC stretch, bend	
    fgets(sbuf,256,fp);
    if (sscanf(sbuf,"%lf%lf%lf%lf",&req,&kstr,&theq,&kbend) != 4) {
	fprintf(stderr, "sysInit: error reading DDC stretch/bend\n");
	exit(1);
    }
    DDCkstr = kstr/KCAL;
    DDCreq = req;
    DDCkbend = kbend/KCAL;
    DDCtheq = theq*PI/180.0;
    // DDC tors, 1-4 factor
    fgets(sbuf,256,fp);
    if (sscanf(sbuf,"%lf%lf%lf%lf",&t0, &t1, &t2, &f14 ) != 4) {
	fprintf(stderr, "sysInit: error reading DDC tors, 1-4 factor\n");
	exit(1);
    }
    pTors[0]=t0/=KCAL;
    pTors[1]=t1/=KCAL;
    pTors[2]=t2/=KCAL;
    factor14=f14;
    // water O,H sig, eps, q
    fgets(sbuf,256,fp);
    if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
	fprintf(stderr, "sysInit: error reading water O LJq\n");
	exit(1);
    }
    osig=sig;
    oeps=eps;
    oq=q;
    fgets(sbuf,256,fp);
    if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
	fprintf(stderr, "sysInit: error reading water H LJq\n");
	exit(1);
    }
    hsig=sig;
    heps=eps;
    hq=q;
    // Er3+ sig, eps, q
    fgets(sbuf,256,fp);
    if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
	fprintf(stderr, "sysInit: error reading Er3+ LJq\n");
	exit(1);
    }
    esig=sig;
    eeps=eps;
    eq=q;

    for(i=0;i<DDCsites;i++){
	if(i==0 || i==DDCsites-1){
	    sig = (esig+c1sig)/2.0;
	    eps = sqrt(eeps*c1eps);
	    EDlj[i].q = eq*c1q/E2;
	}
	else{
	    sig = (esig+c2sig)/2.0;
	    eps = sqrt(eeps*c2eps);
	    EDlj[i].q = eq*c2q/E2;
	}
	s6 = sig*sig*sig;
	s6 = s6*s6;
	EDlj[i].b=4.0*eps*s6/KCAL;
	EDlj[i].a=EDlj[i].b*s6;
    }
    i=0;
    for(j=0;j<DDCsites;j++){
	if(j==0 || j==DDCsites-1){
	    sig = (c1sig+osig)/2.0;
	    eps = sqrt(c1eps*oeps);
	    WDlj[i][j].q=0.0;//oq*c1q/E2;
	}
	else{
	    sig = (c2sig+osig)/2.0;
	    eps = sqrt(c2eps*oeps);
	    WDlj[i][j].q=0.0;//oq*c2q/E2;
	}
	s6 = sig*sig*sig;
	s6 = s6*s6;
	WDlj[i][j].b=4.0*eps*s6/KCAL;
	WDlj[i][j].a=WDlj[i][j].b*s6;
    }
    for(i=1;i<3;i++){
	for(j=0;j<DDCsites;j++){
	    if(j==0 || j==11){
		sig = (c1sig+hsig)/2.0;
		eps = sqrt(c1eps*heps);
		WDlj[i][j].q=0.0;//hq*c1q/E2;
	    }
	    else{
		sig = (c2sig+hsig)/2.0;
		eps = sqrt(c2eps*heps);
		WDlj[i][j].q=0.0;//hq*c2q/E2;
	    }
	    s6 = sig*sig*sig;
	    s6 = s6*s6;
	    WDlj[i][j].b=4.0*eps*s6/KCAL;
	    WDlj[i][j].a=WDlj[i][j].b*s6;
	}
    }
    i=0;
    sig = (osig+esig)/2.0;
    eps = sqrt(oeps*eeps);
    EWlj[i].q=oq*eq/E2;
    s6 = sig*sig*sig;
    s6 = s6*s6;
    EWlj[i].b=4.0*eps*s6/KCAL;
    EWlj[i].a=EWlj[i].b*s6;

    for(i=1;i<3;i++){
	sig = (hsig+esig)/2.0;
	eps = sqrt(heps*eeps);
	EWlj[i].q=hq*eq/E2;
	s6 = sig*sig*sig;
	s6 = s6*s6;
	EWlj[i].b=4.0*eps*s6/KCAL;
	EWlj[i].a=EWlj[i].b*s6;
    }
    dsig[0] = dsig[DDCsites-1] = c1sig;
    deps[0] = deps[DDCsites-1] = c1eps;
    for(i=1;i<DDCsites-1;i++){
	dsig[i] = c2sig;
	deps[i] = c2eps;
    }
    for(i=0;i<DDCsites;i++){
	for(j=0;j<DDCsites;j++){
	    sig=(dsig[i]+dsig[j])/2.0;
	    eps=sqrt(deps[i]*deps[j]);
	    DDlj[i][j].q = 0.0;
	    s6 = sig*sig*sig;
	    s6 = s6*s6;
	    DDlj[i][j].b=4.0*eps*s6/KCAL;
	    DDlj[i][j].a=DDlj[i][j].b*s6;
	}
    }
    // windowing parmeters
    fgets(sbuf,256,fp);
    if (sscanf(sbuf,"%d%lf%lf%lf%lf",&windowOn, &center_w, &width_w, &power_w,&pot_w) != 5) {
	fprintf(stderr, "sysInit: error reading windowing potential\n");
	fprintf(stderr, "sysInit: Windowing potential disabled.\n");
	windowOn = 0;
	center_w = width_w = power_w  = pot_w = 0.0;
    }
    // biasing parmeters
    fgets(sbuf,256,fp);
    if (sscanf(sbuf,"%d%lf%lf%lf%lf%lf%lf",&biasingOn, &BIAS_W, &BIAS_A, &BIAS_B, &BIAS_C, &BIAS_D, &BIAS_G ) != 7) {
	fprintf(stderr, "sysInit: error reading biasing potential\n");
	fprintf(stderr, "sysInit: Biasing potential disabled.\n");
	biasingOn = 0;
	BIAS_W = BIAS_A = BIAS_B = BIAS_C = BIAS_D = BIAS_G = 0.0;
    }
    // flat interface
    fgets(sbuf,256,fp);
    if (sscanf(sbuf,"%d",&KillFinger) != 1) {
	fprintf(stderr, "sysInit: error reading KillFinger\n");
	fprintf(stderr, "sysInit: flat interface not enabled.\n");
	KillFinger = 0;
    }
    fgets(sbuf,256,fp);
    if (sscanf(sbuf,"%d%lf%lf%lf%lf",
		&chSw,&swmax,&swmin,&swSoMax,&swSoMin) != 5) {
	fprintf(stderr, "sysinit: error reading switching parameters. set to default.\n");
	chSw = 0;
    }
    if (chSw==1){// Change switching parameters 
	fprintf(stderr, "sysinit: manual switching parameters: swmax = %f, swmin = %f, swSoMax = %f, swSoMin = %f\n",swmax,swmin,swSoMax,swSoMin);
	if (pbcType[0] == 'C'){
	    bc_factor = 1.; // for cubic 
	}
	else if (pbcType[0] == 'O'){
	    bc_factor = 0.75; // for truncated octahedron 
	}
	else
	    ERROR((stderr,"sysinit:periodic boundaries not defined\n"), exit);

	if (xwall > ywall || xwall > zwall)
	    ERROR((stderr,"sysinit: xwall must be smaller than y or z\n"), exit);
	swr2max = bc_factor * xwall * xwall;
	// This is the square of 1/2 the distance between the hexagonal
	// faces of the pto box, or the sq of 1/2 the small length of
	// the cube 

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
	dz = swr2max - swr2min; // r**2 (min max) for switching region 
	dz3 = (dz*dz*dz);
	dz4 = dz3*dz;
	dz5 = dz4*dz;
	swcoef[0] =  1.0; // switching parameters for group cut-offs 
	swcoef[1] =  -10.0 / dz3;
	swcoef[2] =   15.0 / dz4;
	swcoef[3] =   -6.0 / dz5;
	pswcoef[0] = -60.0 / dz3;	// 2 times swcoef's derivs 
	pswcoef[1] = 120.0 / dz4;
	pswcoef[2] = -60.0 / dz5;
	// for solute-solvent
	dz = swSoMax - swSoMin;// r**2 (min max) for switching region 
	dz3 = (dz*dz*dz);
	dz4 = dz3*dz;
	dz5 = dz4*dz;
	swSocoef[0] =1.0; // switching parameters for solvent-solute
	swSocoef[1] =  -10.0 / dz3;
	swSocoef[2] =   15.0 / dz4;
	swSocoef[3] =   -6.0 / dz5;
	pswSocoef[0] = -60.0 / dz3;	// 2 times swSocoef's derivs 
	pswSocoef[1] = 120.0 / dz4;
	pswSocoef[2] = -60.0 / dz5;
    }
    else {//set solute-solvent to defaut values
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
    fclose(fp);
}
