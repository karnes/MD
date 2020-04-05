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
	double 	sig,eps,q,q2,s6;
	double kr,req,kth,theq,p1,p2,p3;
	double gsig[GLYsites+2],geps[GLYsites+2],gq[GLYsites+2];//+1 for Br- terms +1 for Cl terms
	double tsig[THAsites+2],teps[THAsites+2],tq[THAsites+2];
	double wsig[3],weps[3],wq[3];
	double ssig[6],seps[6],sq[6]; 
	double swmax, swmin;
     	double dz, dz3, dz4, dz5, bc_factor;
	int chSw;

	*waterP = 'S';

	if ((fp = fopen("/home/jkarnes/MD/GLY/run/inpfiles/THAbonds.txt","r")) == NULL){
	   fprintf(stderr, "sysinit: cannot open THA bonding info\n");
	   exit(1);
	}
	// read in THA intra-interactions
	fgets(sbuf,256,fp);
	fgets(sbuf,256,fp);
	fgets(sbuf,256,fp);
	for(i=0;i<nTHAstr;i++){
		fgets(sbuf,256,fp);
		if(sscanf(sbuf,"%d%d%d",&ii,&jj,&id) != 3) {
		   fprintf(stderr, "sysInit: error reading THA stretches\n");
		   exit(1);
		}
		THAstr[i][0]=ii;
		THAstr[i][1]=jj;
		THAstr[i][2]=id;
	}
	fgets(sbuf,256,fp);
	fgets(sbuf,256,fp);
	fgets(sbuf,256,fp);
	fgets(sbuf,256,fp);
	for(i=0;i<nTHAbend1;i++){
		fgets(sbuf,256,fp);
		if(sscanf(sbuf,"%d%d%d%d",&ii,&jj,&kk,&id) != 4) {
		   fprintf(stderr, "sysInit: error reading THA bends (intra-unit)\n");
		   exit(1);
		}
		THAbend1[i][0]=ii;
		THAbend1[i][1]=jj;
		THAbend1[i][2]=kk;
		THAbend1[i][3]=id;
	}
	fgets(sbuf,256,fp);
	fgets(sbuf,256,fp);
	fgets(sbuf,256,fp);
	fgets(sbuf,256,fp);
	for(i=0;i<nTHAbend2;i++){
		fgets(sbuf,256,fp);
		if(sscanf(sbuf,"%d%d%d%d",&ii,&jj,&kk,&id) != 4) {
		   fprintf(stderr, "sysInit: error reading THA bends (inter-unit)\n");
		   exit(1);
		}
		THAbend2[i][0]=ii;
		THAbend2[i][1]=jj;
		THAbend2[i][2]=kk;
		THAbend2[i][3]=id;
	}
	fgets(sbuf,256,fp);
	fgets(sbuf,256,fp);
	fgets(sbuf,256,fp);
	fgets(sbuf,256,fp);
	for(i=0;i<nTHAtors1;i++){
		fgets(sbuf,256,fp);
		if(sscanf(sbuf,"%d%d%d%d%d",&ii,&jj,&kk,&ll,&id) != 5) {
		   fprintf(stderr, "sysInit: error reading THA torsions (intra-unit); i=%d\n",i);
		   exit(1);
		}
		THAtors1[i][0]=ii;
		THAtors1[i][1]=jj;
		THAtors1[i][2]=kk;
		THAtors1[i][3]=ll;
		THAtors1[i][4]=id;
	}
	fgets(sbuf,256,fp);
	fgets(sbuf,256,fp);
	fgets(sbuf,256,fp);
	fgets(sbuf,256,fp);
	for(i=0;i<nTHAtors2;i++){
		fgets(sbuf,256,fp);
		if(sscanf(sbuf,"%d%d%d%d%d",&ii,&jj,&kk,&ll,&id) != 5) {
		   fprintf(stderr, "sysInit: error reading THA torsions (inter-unit); i=%d\n",i);
		   exit(1);
		}
		THAtors2[i][0]=ii;
		THAtors2[i][1]=jj;
		THAtors2[i][2]=kk;
		THAtors2[i][3]=ll;
		THAtors2[i][4]=id;
	}
	fgets(sbuf,256,fp);
	fgets(sbuf,256,fp);
	fgets(sbuf,256,fp);
	fgets(sbuf,256,fp);
	for(i=0;i<nTHA15;i++){
		fgets(sbuf,256,fp);
		if(sscanf(sbuf,"%d%d",&ii,&jj) != 2) {
		   fprintf(stderr, "sysInit: error reading THA 1-5 pairs\n");
		   exit(1);
		}
		THA15[i][0]=ii;
		THA15[i][1]=jj;
// fprintf(stderr,"15a = %d, 15b = %d\n",GLY15[i][0],GLY15[i][1]);
	}
	fgets(sbuf,256,fp);
	fgets(sbuf,256,fp);
//	fprintf(stderr,"this should say 'end': %s\n",sbuf);
	fgets(sbuf,256,fp);

	close(fp);

	if ((fp = fopen("/home/jkarnes/MD/GLY/run/inpfiles/GLYbonds.txt","r")) == NULL){
	   fprintf(stderr, "sysinit: cannot open b-CD non-bonded pairs\n");
	   exit(1);
	}
	// read in GLY intra-interactions
	fgets(sbuf,256,fp);
	fgets(sbuf,256,fp);
	fgets(sbuf,256,fp);
	for(i=0;i<nGLYstr;i++){
		fgets(sbuf,256,fp);
		if(sscanf(sbuf,"%d%d%d",&ii,&jj,&id) != 3) {
		   fprintf(stderr, "sysInit: error reading GLY stretches\n");
		   exit(1);
		}
		GLYstr[i][0]=ii;
		GLYstr[i][1]=jj;
		GLYstr[i][2]=id;
	}
	fgets(sbuf,256,fp);
	fgets(sbuf,256,fp);
	fgets(sbuf,256,fp);
	fgets(sbuf,256,fp);
	for(i=0;i<nGLYbend;i++){
		fgets(sbuf,256,fp);
		if(sscanf(sbuf,"%d%d%d%d",&ii,&jj,&kk,&id) != 4) {
		   fprintf(stderr, "sysInit: error reading GLY bends\n");
		   exit(1);
		}
		GLYbend[i][0]=ii;
		GLYbend[i][1]=jj;
		GLYbend[i][2]=kk;
		GLYbend[i][3]=id;
	}
	fgets(sbuf,256,fp);
	fgets(sbuf,256,fp);
	fgets(sbuf,256,fp);
	fgets(sbuf,256,fp);
	for(i=0;i<nGLYtors;i++){
		fgets(sbuf,256,fp);
		if(sscanf(sbuf,"%d%d%d%d%d",&ii,&jj,&kk,&ll,&id) != 5) {
		   fprintf(stderr, "sysInit: error reading GLY torsions; i=%d\n",i);
		   exit(1);
		}
		GLYtors[i][0]=ii;
		GLYtors[i][1]=jj;
		GLYtors[i][2]=kk;
		GLYtors[i][3]=ll;
		GLYtors[i][4]=id;
	}
	fgets(sbuf,256,fp);
	fgets(sbuf,256,fp);
	fgets(sbuf,256,fp);
	fgets(sbuf,256,fp);
	for(i=0;i<nGLY15;i++){
		fgets(sbuf,256,fp);
		if(sscanf(sbuf,"%d%d",&ii,&jj) != 2) {
		   fprintf(stderr, "sysInit: error reading GLY 1-5 pairs\n");
		   exit(1);
		}
		GLY15[i][0]=ii;
		GLY15[i][1]=jj;
// fprintf(stderr,"15a = %d, 15b = %d\n",GLY15[i][0],GLY15[i][1]);
	}
	fgets(sbuf,256,fp);
	fgets(sbuf,256,fp);
//	fprintf(stderr,"this should say 'end': %s\n",sbuf);
	fgets(sbuf,256,fp);

	close(fp);

	if ((fp = fopen(confile,"r")) == NULL) {
		fprintf(stderr, "sysinit: cannot open %s\n", confile);
		exit(1);
	}

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%d%d%d%d%d", &nGLY,&nTHA,&nBr,&nCl2,&nTS) != 5) {
		fprintf(stderr, "sysInit: error reading line 1\n");
		exit(1);
	}
	setGmax();

/* Read in GLY parameters */
	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading GLY cC LJq\n");
		exit(1);
	}
	gsig[0] = sig;
	geps[0] = eps;
	gq[0] = q;

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading GLY C LJq\n");
		exit(1);
	}
	gsig[4] = gsig[9] = sig;
	geps[4] = geps[9] = eps;
	gq[4] = gq[9] = q;

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading GLY O LJq\n");
		exit(1);
	}
	gsig[1] = gsig[5] = gsig[10] = sig;
	geps[1] = geps[5] = geps[10] = eps;
	gq[1] = gq[5] = gq[10] = q;

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading GLY H(O) LJq\n");
		exit(1);
	}
	gsig[2] = gsig[6] = gsig[11] = sig;
	geps[2] = geps[6] = geps[11] = eps;
	gq[2] = gq[6] = gq[11] = q;

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf%lf",&sig, &eps, &q, &q2 ) != 4) {
		fprintf(stderr, "sysInit: error reading GLY H(C) LJq\n");
		exit(1);
	}
	gsig[7] = gsig[8] = gsig[3] = gsig[12] = gsig[13] = sig;
	geps[7] = geps[8] = geps[3] = geps[12] = geps[13] = eps;
	gq[7] = gq[8] = gq[3] = gq[12] = gq[13] = q;
	gq[3] = q2;

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading Br- LJq\n");
		exit(1);
	}
	gsig[GLYsites] = tsig[THAsites] = sig;
	geps[GLYsites] = teps[THAsites] = eps;
	gq[GLYsites] = tq[THAsites] = q;

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading Cl2 LJq\n");
		exit(1);
	}
	ssig[2]=sig;
	seps[2]=eps;
	sq[2]=q;

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading M LJq\n");
		exit(1);
	}
	ssig[0]=sig;
	seps[0]=eps;
	sq[0]=q;

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading M0 LJq\n");
		exit(1);
	}
	ssig[1]=sig;
	seps[1]=eps;
	sq[1]=q;

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading Cl(TS1) LJq\n");
		exit(1);
	}
	ssig[4]=sig;
	seps[4]=eps;
	sq[4]=q;

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading Cl(TSc) LJq\n");
		exit(1);
	}
	ssig[3]=sig;
	seps[3]=eps;
	sq[3]=q;

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading Br(TS) LJq\n");
		exit(1);
	}
	ssig[5]=sig;
	seps[5]=eps;
	sq[5]=q;
	for(i=0;i<4;i++){
	   fgets(sbuf,256,fp);
	   if (sscanf(sbuf,"%lf%lf",&kr, &req ) != 2) {
		fprintf(stderr, "sysInit: error reading GLY stretches\n");
		exit(1);
	   }
	   gkstr[i] = kr/KCAL;
	   greq[i] = req;
	}

	for(i=0;i<6;i++){
	   fgets(sbuf,256,fp);
	   if (sscanf(sbuf,"%lf%lf",&kth, &theq ) != 2) {
		fprintf(stderr, "sysInit: error reading GLY bends\n");
		exit(1);
	   }
	   gkbend[i] = kth/KCAL;
	   geqbend[i] = theq*PI/180.0;
	}

	for(i=0;i<7;i++){
	   fgets(sbuf,256,fp);
	   if (sscanf(sbuf,"%lf%lf%lf",&p1, &p2, &p3 ) != 3) {
		fprintf(stderr, "sysInit: error reading GLY torsions.. i=%d\n",i);
		exit(1);
	   }
	   gtors[i][0] = p1/KCAL;
	   gtors[i][1] = p2/KCAL;
	   gtors[i][2] = p3/KCAL;
	}
// calculate GLY-GLY mixed interactions
	for(i=0;i<GLYsites+1;i++){
	   for(j=0;j<GLYsites+1;j++){
	      sig = (gsig[i]+gsig[j])/2.0;
	      eps = sqrt(geps[i]*geps[j]);
//	      if(i==j)
//	        fprintf(stderr,"GLY:i=j=%d, sig = %f, eps = %f\n",i,sig,eps);
	      s6 = sig*sig*sig;
	      s6 = s6*s6;
	      GLYlj[i][j].b = GLYlj[j][i].b = 4.0*eps*s6/KCAL;
	      GLYlj[i][j].a = GLYlj[j][i].a = GLYlj[i][j].b * s6;
	      GLYlj[i][j].q = GLYlj[j][i].q = gq[i]*gq[j]/E2;
	   }
	}
/* Read in THA parameters */
	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading THA N LJq\n");
		exit(1);
	}
	tsig[0] = sig;
	teps[0] = eps;
	tq[0] = q;

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading THA C1 LJq\n");
		exit(1);
	}
	tsig[1] = sig;
	teps[1] = eps;
	tq[1] = q;

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading THA H2-4 LJq\n");
		exit(1);
	}
	tsig[2] = tsig[3] = tsig[4] = sig;
	teps[2] = teps[3] = teps[4] = eps;
	tq[2] = tq[3] = tq[4] = q;

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading THA C5 LJq\n");
		exit(1);
	}
	tsig[5] = sig;
	teps[5] = eps;
	tq[5] = q;

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading THA H6-7 LJq\n");
		exit(1);
	}
	tsig[6] = tsig[7] = sig;
	teps[6] = teps[7] = eps;
	tq[6] = tq[7] = q;

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading THA C8 LJq\n");
		exit(1);
	}
	tsig[8] = sig;
	teps[8] = eps;
	tq[8] = q;

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading THA H9-10 LJq\n");
		exit(1);
	}
	tsig[9] = tsig[10] = sig;
	teps[9] = teps[10] = eps;
	tq[9] = tq[10] = q;

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading THA C11 LJq\n");
		exit(1);
	}
	tsig[11] = sig;
	teps[11] = eps;
	tq[11] = q;

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading THA H12-13 LJq\n");
		exit(1);
	}
	tsig[12] = tsig[13] = sig;
	teps[12] = teps[13] = eps;
	tq[12] = tq[13] = q;

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading THA C14 LJq\n");
		exit(1);
	}
	tsig[14] = sig;
	teps[14] = eps;
	tq[14] = q;

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading THA H15-16 LJq\n");
		exit(1);
	}
	tsig[15] = tsig[16] = sig;
	teps[15] = teps[16] = eps;
	tq[15] = tq[16] = q;

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading THA C17 LJq\n");
		exit(1);
	}
	tsig[17] = sig;
	teps[17] = eps;
	tq[17] = q;

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf%lf",&sig, &eps, &q ) != 3) {
		fprintf(stderr, "sysInit: error reading THA H18-19 LJq\n");
		exit(1);
	}
	tsig[18] = tsig[19] = sig;
	teps[18] = teps[19] = eps;
	tq[18] = tq[19] = q;

	for(i=0;i<3;i++){
	   fgets(sbuf,256,fp);
	   if (sscanf(sbuf,"%lf%lf",&kr, &req ) != 2) {
		fprintf(stderr, "sysInit: error reading THA stretches\n");
		exit(1);
	   }
	   tkstr[i] = kr/KCAL;
	   treq[i] = req;
	}

	for(i=0;i<6;i++){
	   fgets(sbuf,256,fp);
	   if (sscanf(sbuf,"%lf%lf",&kth, &theq ) != 2) {
		fprintf(stderr, "sysInit: error reading THA bends\n");
		exit(1);
	   }
	   tkbend[i] = kth/KCAL;
	   teqbend[i] = theq*PI/180.0;
	}

	for(i=0;i<7;i++){
	   fgets(sbuf,256,fp);
	   if (sscanf(sbuf,"%lf%lf%lf",&p1, &p2, &p3 ) != 3) {
		fprintf(stderr, "sysInit: error reading THA torsions\n");
		exit(1);
	   }
	   ttors[i][0] = p1/KCAL;
	   ttors[i][1] = p2/KCAL;
	   ttors[i][2] = p3/KCAL;
	}
	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf",&kr, &req ) != 2) {
		fprintf(stderr, "sysInit: error reading Cl2 stretch\n");
		exit(1);
	}
	Clkstr = kr/KCAL;
	Clreq = req;

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf",&kr, &req ) != 2) {
		fprintf(stderr, "sysInit: error reading Cl-Cl(TS) stretch\n");
		exit(1);
	}
	TSkstr[0] = kr/KCAL;
	TSreq[0] = req;

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf",&kr, &req ) != 2) {
		fprintf(stderr, "sysInit: error reading Cl-Br(TS) stretch\n");
		exit(1);
	}
	TSkstr[1] = kr/KCAL;
	TSreq[1] = req;

	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%lf%lf",&kth, &theq ) != 2) {
		fprintf(stderr, "sysInit: error reading THA bends\n");
		exit(1);
	}
	TSkbend = kth/KCAL;
	TSeqbend = theq*PI/180.0;


// populate rest of sig eps q
	for(j=1;j<4;j++){
	   for(i=0;i<19;i++){
	      tsig[j*19+i+1] = tsig[i];
	      teps[j*19+i+1] = teps[i];
	      tq[j*19+i+1] = tq[i];
	   }
	}
// calculate THA-THA mixed interactions
	for(i=0;i<THAsites+1;i++){
	   for(j=0;j<THAsites+1;j++){
	      sig = (tsig[i]+tsig[j])/2.0;
	      eps = sqrt(teps[i]*teps[j]);
//	      if(i==j)
//	        fprintf(stderr,"THA:i=j=%d, sig = %f, eps = %f\n",i,sig,eps);
	      s6 = sig*sig*sig;
	      s6 = s6*s6;
	      THAlj[i][j].b = THAlj[j][i].b = 4.0*eps*s6/KCAL;
	      THAlj[i][j].a = THAlj[j][i].a = THAlj[i][j].b * s6;
	      THAlj[i][j].q = THAlj[j][i].q = tq[i]*tq[j]/E2;
	   }
	}
// calculate GLY-THA mixed interactions
	for(i=0;i<GLYsites;i++){
	   for(j=0;j<THAsites;j++){
	      sig = (gsig[i]+tsig[j])/2.0;
	      eps = sqrt(geps[i]+teps[j]);
	      s6 = sig*sig*sig;
	      s6 = s6*s6;
	      GTlj[i][j].b = 4.0*eps*s6/KCAL;
	      GTlj[i][j].a = GTlj[i][j].b*s6;
	      GTlj[i][j].q = gq[i]*tq[j]/E2;
	   }
	}
wsig[0] = 3.16554;
wsig[1] = wsig[2] = 0.0;
weps[0] = 0.1554;
weps[1] = weps[2] = 0.0;
wq[0] = -0.82;
wq[1] = wq[2] = 0.41;
// calculate water-THA mixed interactions
	for(i=0;i<3;i++){
	   for(j=0;j<THAsites;j++){
	      sig = (wsig[i]+tsig[j])/2.0;
	      eps = sqrt(weps[i]+teps[j]);
	      s6 = sig*sig*sig;
	      s6 = s6*s6;
	      WTlj[i][j].b = 4.0*eps*s6/KCAL;
	      WTlj[i][j].a = WTlj[i][j].b*s6;
	      WTlj[i][j].q = wq[i]*tq[j]/E2;
	   }
	}
// calculate GLY-solute mixed interactions
	for(i=0;i<GLYsites;i++){
	   for(j=0;j<6;j++){
	      sig = (gsig[i]+ssig[j])/2.0;
	      eps = sqrt(geps[i]+seps[j]);
	      s6 = sig*sig*sig;
	      s6 = s6*s6;
	      solGLYlj[i][j].b = 4.0*eps*s6/KCAL;
	      solGLYlj[i][j].a = solGLYlj[i][j].b*s6;
	      solGLYlj[i][j].q = gq[i]*sq[j]/E2;
	   }
	}
// calculate THA-solute mixed interactions
	for(i=0;i<THAsites;i++){
	   for(j=0;j<6;j++){
	      sig = (tsig[i]+ssig[j])/2.0;
	      eps = sqrt(teps[i]+seps[j]);
	      s6 = sig*sig*sig;
	      s6 = s6*s6;
	      solTHAlj[i][j].b = 4.0*eps*s6/KCAL;
	      solTHAlj[i][j].a = solTHAlj[i][j].b*s6;
	      solTHAlj[i][j].q = tq[i]*sq[j]/E2;
	   }
	}
	fgets(sbuf,256,fp);
	if(sscanf(sbuf,"%d%lf%lf%lf%lf",
	   &chSw,&swmax,&swmin,&swSoMax,&swSoMin) != 5) {
	   fprintf(stderr, "sysinit: error reading switching parameters. set to default.\n");
	   chSw = 0;
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

/*	fgets(sbuf,256,fp);
	if(sscanf(sbuf,"%d%lf%lf",
	   &windowOn,&center_r,&window_w) != 3) {
	   fprintf(stderr, "sysinit: error reading windowing potential. windowing potential off.\n");
	   windowOn = 0;
	   center_r = 0.0;
	   window_w = 999.9;
	}
	else if(windowOn==1){
	   fprintf(stderr, "sysinit: windowing potential on: center = %4.2f, width = %4.2f\n",center_r,window_w);
	}
	fgets(sbuf,256,fp);
	if(sscanf(sbuf,"%d%lf%lf%lf",
	   &biasOn,&a_bias,&b_bias,&o_bias) != 4) {
	   fprintf(stderr, "sysinit: error reading biasing. biasing off.\n");
	   biasOn = 0;
	   a_bias = b_bias = o_bias = 0.0;
	}
	else if(biasOn==1){
	   fprintf(stderr, "sysinit: biasin on: Ub = a(x-o)^2 + b(x-o) : a = %4.2f, b = %4.2f = o = %4.2f\n",a_bias,b_bias,o_bias);
	}
	fgets(sbuf,256,fp);
	if(sscanf(sbuf,"%d%lf",
	   &XWindowOn,&Xmin) != 2) {
	   fprintf(stderr, "sysinit: No constraint on THA z-position\n");
	   XWindowOn = 0;
	   Xmin = -999.0;
	}
	else if(XWindowOn==1){
	   fprintf(stderr, "sysinit: Halogen CoM confined to z > %5.2f\n",Xmin);
	}
	else{
	   fprintf(stderr, "sysinit: No constraint on halogen z-position\n");
	}
	fgets(sbuf,256,fp);
	if(sscanf(sbuf,"%d%lf%lf",
	   &kappaWindowOn,&kap_c,&kap_w) != 3) {
	   fprintf(stderr, "sysinit: error reading kappa windowing. kappa windowing off.\n");
	   kappaWindowOn = 0;
	   kap_c = 0.0;
	   kap_w = 99.9;
	}
	else{
		fprintf(stderr,"sysInit: kappaWindowOn = %d, window center = %4.3f Å, window width = %4.3f Å.\n",kappaWindowOn,kap_c,kap_w);
	}
*/
	if(chSw==1){/* Change switching parameters */
	   fprintf(stderr, "sysinit: manual switching parameters: swmax = %f, swmin = %f, swSoMax = %f, swSoMin = %f\n",swmax,swmin,swSoMax,swSoMin);
	   if(pbcType[0] == 'C'){
	      bc_factor = 1.; /* for cubic */
	   }
	   else if(pbcType[0] == 'O'){
	      	bc_factor = 0.75; /* for truncated octahedron */
	   }
	   else
	      ERROR((stderr,"sysinit:periodic boundaries not defined\n"), exit);

           if(xwall > ywall || xwall > zwall)
	      ERROR((stderr,"sysinit: xwall must be smaller than y or z\n"), exit);
	      swr2max = bc_factor * xwall * xwall;
	      /* This is the square of 1/2 the distance between the hexagonal
	       * faces of the pto box, or the sq of 1/2 the small length of
	       * the cube */
	 
	      swr2min = bc_factor * (xwall - 1.0) * (xwall - 1.0);
	      if(swmax*swmax > swr2max)	
	         ERROR((stderr,"sysinit: swmax larger than maximum\n"), exit);
	      if(swmin*swmin > swr2min)	
	         ERROR((stderr,"sysinit: swmin larger than maximum\n"), exit);
	      swSoMax = swSoMax*swSoMax;
	      swSoMin = swSoMin*swSoMin;
	      if(swSoMax > swr2max)	
	      	 ERROR((stderr,"sysinit: swSoMax larger than maximum\n"), exit);
	      if(swSoMin > swr2min)	
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
	else
        {  /*set solute-solvent to defaut values*/
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
/*
// generate list of 1-5+ pairs
// can't be in a torsion
	int totTors[4*nTHAtors1+nTHAtors2][4];
	int prnt;
	for(i=0;i<nTHAtors1;i++){
	   for(j=0;j<4;j++){
	      for(k=0;k<4;k++){
	         totTors[i+j*nTHAtors1][k]=(THAtors1[i][k]==0)?0:THAtors1[i][k]+j*19;
	      }
	   }
	}
	for(i=0;i<nTHAtors2;i++){
	   for(k=0;k<4;k++){
	      totTors[i+4*nTHAtors1][k]=THAtors2[i][k];
	   }
	}
//can't be in a bend
	int totBend[4*nTHAbend1][3];
	for(i=0;i<nTHAbend1;i++){
	   for(j=0;j<4;j++){
              for(k=0;k<3;k++){
		 totBend[i+j*nTHAbend1][k]=THAbend1[i][k]==0?0:j*19+THAbend1[i][k];
	      }
	   }
	}
	for(i=0;i<THAsites;i++){
	   for(j=i+1;j<THAsites;j++){
	      prnt = 1;
	      for(k=0;k<4*nTHAtors1+nTHAtors2;k++){
		 if(i==totTors[k][0]||i==totTors[k][1]||i==totTors[k][2]||i==totTors[k][3]){
		    if(j==totTors[k][0]||j==totTors[k][1]||j==totTors[k][2]||j==totTors[k][3]){
		       prnt = 0;
	               if(i==70&&j==76)
		          fprintf(stderr,"PANIC: %d %d %d %d, k=%d\n",totTors[k][0],totTors[k][1],totTors[k][2],totTors[k][3],k);
		    }
		 }
	      }
	      //also check bends
	      for(k=0;k<4*nTHAbend1;k++){
		 if(i==totBend[k][0]||i==totBend[k][1]||i==totBend[k][2]){
		    if(j==totBend[k][0]||j==totBend[k][1]||j==totBend[k][2]){
		       prnt = 0;
	            //   fprintf(stderr,"FOUND A 1-4~\n");
		    }
		 }
	      }
	      if(prnt==1)
	         fprintf(stderr,"%d %d\n",i,j);
	   }
	}
	
*/
fclose(fp);


}

int setGmax()
{
	int i;
	   if(nBr==1){
	   gGOXmax[0] = 4.40;
	   gGCXmax[0] = 5.05;
	   gGHXmax[0] = 3.60;
	}
	else if(nCl2==1){
	   gGOXmax[0] = 4.35;
	   gGCXmax[0] = 4.95;
	   gGHXmax[0] = 4.70;
	   gGOXmax[2] = 4.90;
	   gGCXmax[2] = 7.00;
	   gGHXmax[2] = 5.30;
	}
	else if(nTS==1){
	   gGOXmax[0] = 4.25;
	   gGCXmax[0] = 4.90;
	   gGHXmax[0] = 3.75;
	   gGOXmax[1] = 4.00;
	   gGCXmax[1] = 5.05;
	   gGHXmax[1] = 3.60;
	   gGOXmax[2] = 4.60;
	   gGCXmax[2] = 5.35;
	   gGHXmax[2] = 4.10;
	}
for(i=0;i<3;i++){
//   fprintf(stderr,"gGCXmax[ %d ] = %4.2f , gGOXmax[ %d ] = %4.2f, gGHXmax[ %d ] = %4.2f \n", i, gGCXmax[i], i, gGOXmax[i], i, gGHXmax[i]);
}

	return(0);
}
	   
	   




