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
	int	ns, chSw,jns, i, j;
	LJ	lj[19];
	double swmax, swmin, EField;
	double 	sig,eps,q, coul,s6,sqrt();
	 	
     	double	dz, dz3, dz4, dz5, bc_factor;

	if ((fp = fopen(confile,"r")) == NULL) {
		fprintf(stderr, "sysinit: cannot open %s\n", confile);
		exit(1);
	}

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%d%d%s",&nNIT, &ns, waterP) != 3) {
		fprintf(stderr, "sysInit: error reading line 1\n");
		exit(1);
	}
/*read nitrobenzen intramolecular parameters*/
	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&pTors[0][0], &pTors[0][1], &pTors[0][2]) != 3) {
		fprintf(stderr, "sysInit: error reading H-CA-CA-H\n");
		exit(1);
	}

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&pTors[1][0], &pTors[1][1], &pTors[1][2]) != 3) {
		fprintf(stderr, "sysInit: error reading N-CA-CA-H\n");
		exit(1);
	}

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&pTors[2][0], &pTors[2][1], &pTors[2][2]) != 3) {
		fprintf(stderr, "sysInit: error reading N-CA-CA-CA\n");
		exit(1);
	}

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&pTors[3][0], &pTors[3][1], &pTors[3][2]) != 3) {
		fprintf(stderr, "sysInit: error reading O-N-CA-CA\n");
		exit(1);
	}

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&pTors[4][0], &pTors[4][1], &pTors[4][2]) != 3) {
		fprintf(stderr, "sysInit: error reading ring improper\n");
		exit(1);
	}

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&pTors[5][0], &pTors[5][1], &pTors[5][2]) != 3) {
		fprintf(stderr, "sysInit: error reading ring-hydrogen improper\n");
		exit(1);
	}

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&pTors[6][0], &pTors[6][1], &pTors[6][2]) != 3) {
		fprintf(stderr, "sysInit: error reading ring-N improper\n");
		exit(1);
	}

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&pTors[7][0], &pTors[7][1], &pTors[7][2]) != 3) {
		fprintf(stderr, "sysInit: error reading ring-NOT improper\n");
		exit(1);
	}
	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf",&kBend[0], &eqBend[0]) != 2) {
		fprintf(stderr, "sysInit: error reading CA-CA-CA\n");
		exit(1);
	}

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf",&kBend[1], &eqBend[1]) != 2) {
		fprintf(stderr, "sysInit: error reading line H-CA-CA\n");
		exit(1);
	}

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf",&kBend[2], &eqBend[2]) != 2) {
		fprintf(stderr, "sysInit: error reading line N-CA-CA\n");
		exit(1);
	}

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf",&kBend[3], &eqBend[3]) != 2) {
		fprintf(stderr, "sysInit: error reading line O-N-CA\n");
		exit(1);
	}

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf",&kBend[4], &eqBend[4]) != 2) {
		fprintf(stderr, "sysInit: error reading line O-N-O\n");
		exit(1);
	}

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf",&kStr[0], &eqBond[0]) != 2) {
		fprintf(stderr, "sysInit: error reading line CA-CA\n");
		exit(1);
	}

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf",&kStr[1], &eqBond[1]) != 2) {
		fprintf(stderr, "sysInit: error reading line H-CA\n");
		exit(1);
	}

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf",&kStr[2], &eqBond[2]) != 2) {
		fprintf(stderr, "sysInit: error reading line N-CA\n");
		exit(1);
	}

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf",&kStr[3], &eqBond[3]) != 2) {
		fprintf(stderr, "sysInit: error reading line N-O\n");
		exit(1);
	}

	for	(i = 0; i < 3; i++)
		for	(j = 0; j < 8; j++)
			pTors[j][i] /= KCAL;
	
	for	(i = 0; i < 5; i++){
		kStr[i] /= KCAL;
		kBend[i] /= KCAL;
		eqBend[i] *= (PI/180.);
	}

	for ( i = 0; i < 19; i++ )
		lj[i].S = lj[i].E = lj[i].Q = 0.0; 

/*	Get the appropriate input paramters for the NIT inter- forces */
/* first read in CA carbon parameters */
	fgets(sbuf,256,fp);
        if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q) != 3) {
                fprintf(stderr, "sysInit: error reading line 2\n");
                exit(1);
        }
	lj[0].S = sig; 
	lj[0].E = eps;
	lj[0].Q = q;

	fgets(sbuf,256,fp);
        if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q) != 3) {
                fprintf(stderr, "sysInit: error reading line 3\n");
                exit(1);
        }
        lj[1].S = lj[5].S = sig; 
        lj[1].E = lj[5].E = eps;
        lj[1].Q = lj[5].Q = q;

	fgets(sbuf,256,fp);
        if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q) != 3) {
                fprintf(stderr, "sysInit: error reading line 4\n");
                exit(1);
        }
        lj[2].S = lj[4].S = sig;
        lj[2].E = lj[4].E = eps;
        lj[2].Q = lj[4].Q = q;

	fgets(sbuf,256,fp);
        if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q) != 3) {
                fprintf(stderr, "sysInit: error reading line 5\n");
                exit(1);
        }
        lj[3].S = sig; 
        lj[3].E = eps;
        lj[3].Q = q;


/* Now read in the N2 parameters */
	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q) != 3) {
		fprintf(stderr, "sysInit: error reading line 6\n");
		exit(1);
	}
	lj[6].S = sig; 
        lj[6].E = eps;
        lj[6].Q = q;


/* Now read in the O2 parameters */
	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q) != 3) {
		fprintf(stderr, "sysInit: error reading line 7\n");
		exit(1);
	}
	for ( i = 7; i < 9; i++ )
                {
                lj[i].S = sig;
                lj[i].E = eps;
                lj[i].Q = q;
                }


/* Now read in the benzene hydrogen parameters */
	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q) != 3) {
		fprintf(stderr, "sysInit: error reading line 8\n");
		exit(1);
	}
        lj[9].S = lj[13].S = sig;
        lj[9].E = lj[13].E = eps;
        lj[9].Q = lj[13].Q = q;

	fgets(sbuf,256,fp);
        if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q) != 3) {
                fprintf(stderr, "sysInit: error reading line 9\n");
                exit(1);
        }
        lj[10].S = lj[12].S = sig;
        lj[10].E = lj[12].E = eps;
        lj[10].Q = lj[12].Q = q;


	fgets(sbuf,256,fp);
        if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q) != 3) {
                fprintf(stderr, "sysInit: error reading line 10\n");
                exit(1);
        }
        lj[11].S = sig;
        lj[11].E = eps;
        lj[11].Q = q;

/* Read in water parameters */
	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading line 11\n");
		exit(1);
	}
	lj[14].S = sig;
        lj[14].E = eps;
        lj[14].Q = q;
	OQ=q;


	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading line 12\n");
		exit(1);
	}
	for ( i = 15; i < 17; i++ )
                {
                lj[i].S = sig;
                lj[i].E = eps;
                lj[i].Q = q;
                }
	HQ=q;

/*  read in solute parameters here */
if (ns == 1){
	fgets(sbuf,256,fp);
        if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q) != 3) {
            fprintf(stderr, "sysInit: error reading line 13\n");
                 exit(1);
	}
	lj[17].S = sig;
        lj[17].E = eps;
     	lj[17].Q = q;
}
		
/* Now apply mixing rules  and fill in LJNIT global array */
	for ( i  = 0; i < 17+ns; i++ )
		for ( j = 0; j < 17+ns; j++ )
		{
		sig = (lj[i].S+lj[j].S)/2.0; 
		eps = sqrt( lj[i].E * lj[j].E );
		s6 = sig*sig*sig;
		s6 = s6*s6;
		coul = lj[i].Q * lj[j].Q / E2;
		NITlj[i][j].b = NITlj[j][i].b = 4*eps*s6/KCAL;
		NITlj[i][j].a = NITlj[j][i].a = NITlj[i][j].b * s6;
		NITlj[i][j].q = NITlj[j][i].q = coul;
		}

	waterInit();

/*
 *      Read information about modified switching parameters
 */

if      (ns > 0){

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%d%lf%lf%lf%lf",
		&chSw,&swmax,&swmin,&swSoMax,&swSoMin) != 5) {
		fprintf(stderr, "sysinit: error reading line 14\n");
		exit(1);
	}
	if (chSw){/* Change switching parameters */
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

/* read in neqFlags */	
	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%d",&neqFlag) != 1){
		fprintf(stderr, "sysinit: error reading line 15\n");
		exit(1);
	}
/* read in window parameters */	
	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf%lf",&center_w,&width_w,&power_w,&pot_w)!=4) {
                fprintf(stderr, "sysinit: error reading line 16\n");
		exit(1);
	}
        if (width_w > 0. && width_w < 10.){
                fprintf(stderr,"solute constrain to window between\n");
                fprintf(stderr,"%f and %f\n",center_w-width_w/2.,center_w+width_w/2.);
                fprintf(stderr,"relative to system center of mass\n");
        }
	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%d",&EFieldOn) != 1) {
		fprintf(stderr,"sysinit:error reading EFieldOn\n");
		exit(1);
	}
	if (EFieldOn){/* read electric potential (Volts/A) and flags*/
		fgets(sbuf,256,fp);
		if (sscanf(sbuf,"%lf %d %d %d",&EField,&EFieldQ[0],&EFieldQ[1],&EFieldQ[2]) !=4){
			fprintf(stderr,"sysinit:error reading EField parameters\n");
			exit(1);
		}
		for (i=0;i<17+ns;i++)
		EForce[i] = EField*lj[i].Q/EVOLT;
	}
/* read info about hydration shell radius*/
	fgets(sbuf,256,fp);
        if (sscanf(sbuf,"%lf",&Rshel) !=1){
                fprintf(stderr,"sysinit:error reading Rshel\n");
                exit(1);
        }
/* read potential biasing, KillFinger flags */  
	fgets(sbuf,256,fp);
        if (sscanf(sbuf,"%d %d %d",&biasingOn, &KillFinger, &fixShel) !=3){
		biasingOn = 0;
		KillFinger = 0;
		fixShel = 0;
		fprintf(stderr,"sysinit: error reading bias,KilFin, fixShel flags. Set biasingOn=%d, KillFinger=%d fixShel=%d\n",biasingOn,KillFinger,fixShel);
        }
	else fprintf(stderr,"sysInit: biasingOn = %d, KillFinger = %d fixShel = %d\n",biasingOn,KillFinger,fixShel);

}    /* end if ns>0 */
/*
for (i=0; i < 17+ns; i++){
	for (j=0; j < 17+ns; j++){
		printf("%3.3f ",NITlj[i][j].q);
		printf("\n");
	}
}
for (i=0; i < 17+ns; i++){
	for (j=0; j < 17+ns; j++){
		printf("%3.3f ",NITlj[i][j].a);
		printf("\n");
	}
}
for (i=0; i < 17+ns; i++){
	for (j=0; j < 17+ns; j++){
		printf("%3.3f ",NITlj[i][j].b);
		printf("\n");
	}
}
exit(1);
*/
fclose(fp);
}
