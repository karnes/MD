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
	int ii,jj,kk,ll,id;
	double 	sig,eps,q,s6;
	double kr,req,kth,theq,p1,p2,p3;
	double csig[3],ceps[3],cq[3];//3 unique Cl terms
	double bsig,beps,bq; 
	double swmax, swmin;
     	double dz, dz3, dz4, dz5, bc_factor;
	int chSw;

	if ((fp = fopen(confile,"r")) == NULL) {
		fprintf(stderr, "sysinit: cannot open %s\n", confile);
		exit(1);
	}

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%d%d", &nBr,&nCl2) != 2) {
		fprintf(stderr, "sysInit: error reading line 1\n");
		exit(1);
	}

        // Br- LJq
	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading Br- LJq\n");
		exit(1);
	}
	bsig = sig;
	beps = eps;
	bq = q;
        // Cl2 LJq
	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading Cl2 LJq\n");
		exit(1);
	}
	csig[2]=sig;
	ceps[2]=eps;
	cq[2]=q;

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading M LJq\n");
		exit(1);
	}
	csig[0]=sig;
	ceps[0]=eps;
	cq[0]=q;

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading M0 LJq\n");
		exit(1);
	}
	csig[1]=sig;
	ceps[1]=eps;
	cq[1]=q;
	//Cl2 stretch	
        fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf",&kr, &req ) != 2) {
		fprintf(stderr, "sysInit: error reading Cl2 stretch\n");
		exit(1);
	}
	Clkstr = kr/KCAL;
	Clreq = req;

// calculate Br-Cl mixed interactions
	for(i=0;i<3;i++){
	      sig = (csig[i]+bsig)/2.0;
	      eps = sqrt(ceps[i]*beps);
	      s6 = sig*sig*sig;
	      s6 = s6*s6;
	      clj[i].b = 4.0*eps*s6/KCAL;
	      clj[i].a = clj[i].b * s6;
	      clj[i].q = cq[i]*bq/E2;
	}

//	fprintf(stderr, "sysinit: chSw = %d, swmax = %f, swmin = %f\n",chSw,swmax,swmin);
//	fprintf(stderr, "sysinit: chSw = %d, swmax = %f, swmin = %f\n",chSw,swmax,swmin);
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

fclose(fp);

//fprintf(stderr,"sysInit.c complete\n");

}
