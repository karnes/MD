/*
 *	collect:
 *
 *	This routine runs the dynamics for a single trajectory.
 */
#define point_to_skip 0 /* calculate averages using numDP -point_to_skip*/
#include	<md.h>
#include	<system.h>

collectS(numDP, dataRate, outRoot, prFlags, prStatQ, fixTQ) /*ilan*/
	int	numDP;		/* # of data points */
	int	dataRate;	/* # time steps per data point */
	char	outRoot[];	/* root name of output file(s) */
	outflags prFlags;	/* Flags to print output files
				 * if prFlags.dat == 1 print .dat file
				 * with frequency = prFlags
				 * if prFlags.pos == 1 - 4 print .pos file
				 * if prFlags.vel == 1 print .vel file
				 * if prFlags.frc == 1 print .frc file
				 * if prFlags.xtr > 0 print .xtr file
				 * every prFlags.xtr*100 configurations.
				 */
	int	prStatQ;	/* boolean: print status line */
	int	fixTQ;		/* boolean: run constant temperature dynamics */
{
//printf("start collect.c\n");
	double	H2_avg, H_avg, T2_avg, T_avg;
	double P_avg, P2_avg, tension, convFactor ;
	double sqrt(), log();
	int	i, restart_freq;
	FILE	*dfp, *pfp, *vfp, *ffp, *xfp, *ifp;
	char	dfile[80], pfile[80], vfile[80], infofile[80];
	char	ffile[80], xfile[80], rspfile[80], rsvfile[80];
	
#ifdef RDF
	double uu(), factor, pmf; /* john k -- renamed function (u to uu)
* to avoid compilation error with 'u' function in collect.c */
	
	fprintf(ifp,"-Calculate radial distribution function\n");
	fprintf(ifp,"	which will be printed in this file\n");
	for (i =0; i<350; i++)
		rdf[i] = 0.;
#endif RDF
	restart_freq = numDP/10;
		
/*	Generate the desired output file names and open them		*/
	sprintf(rspfile, "%s.rsp", outRoot);
	/*sprintf(rsvfile, "%s.rsv", outRoot);*/
	sprintf(infofile, "%s.inf", outRoot);
	if ((ifp = fopen(infofile, "w")) == NULL) {
		fprintf(stderr, "collect: cannot open infofile\n");
		exit(1);
	}
	if (prFlags.dat) {
		sprintf(dfile, "%s.dat", outRoot);
		if ((dfp = fopen(dfile, "w")) == NULL) {
			fprintf(stderr, "collect: cannot open dfile\n");
			exit(1);
		}
	}
	
	if (prFlags.pos) {
		sprintf(pfile, "%s.pos", outRoot);
		if ((pfp = fopen(pfile, "w")) == NULL) {
			fprintf(stderr, "collect: cannot open pfile\n");
			exit(1);
		}
	}
	
	if (prFlags.vel) {
		sprintf(vfile, "%s.vel", outRoot);
		if ((vfp = fopen(vfile, "w")) == NULL) {
			fprintf(stderr, "collect: cannot open vfile\n");
			exit(1);
		}
	}
	
	if (prFlags.frc) {
		sprintf(ffile, "%s.frc", outRoot);
		if ((ffp = fopen(ffile, "w")) == NULL) {
			fprintf(stderr, "collect: cannot open ffile\n");
			exit(1);
		}
	}
	
	if (prFlags.xtr) {
		sprintf(xfile, "%s.xtr", outRoot);
		if ((xfp = fopen(xfile, "w")) == NULL) {
			fprintf(stderr, "collect: cannot open xfile\n");
			exit(1);
		}
	}

/*	Get the date stamp ...						*/

	time(&datestamp);

/*	Do appropriate initializations ...				*/

//printf("before collect/integ\n");
	integ(1, fixTQ, 1);
	if (fixTQ)
		fprintf(ifp,"-Constant temperature run\n");
	else
		fprintf(ifp,"-Constant energy run\n");
	fprintf(ifp,"-Time step = %f fs\n",h);
	fprintf(ifp,"-Number of configurations = %d\n",numDP+1);
	fprintf(ifp,"	integration step per configuration = %d\n",dataRate);
	fprintf(ifp,"-Temperature = %f\n",Teq);
	if (prFlags.pos == 1){
		fprintf(ifp,"-Writing .pos file every configuration\n");
		wpos(pfp, 1);
	}
	if (prFlags.pos == 2){
		fprintf(ifp,"-Writing shell .pos file every configuration\n");
		fprintf(ifp,"Writing a shell of 8A, but need to fix it\n");
		fprintf(ifp,"(It does not update the shell)\n");
		wshel(pfp, 1);
	}
	if (prFlags.pos == 3){
		fprintf(ifp,"-Writing .pos file every 10th configuration\n");
		wpos(pfp, 1);
	}
	if (prFlags.pos == 4 && nsolute > 0){
		fprintf(ifp,"-Writing .pos file for solute only every configuration\n");
		wpsolute(pfp, 1);
	}
	if (prFlags.vel){
		fprintf(ifp,"-Writing .vel file every configuration\n");
		wvel(vfp, 1);
	}
	if (prFlags.frc == 1){
		fprintf(ifp,"-Writing .frc file every configuration\n");
		wfrc(ffp, 1);
	}
	if (prFlags.frc == 4 && nsolute > 0){
		fprintf(ifp,"-Writing .frc file for solute only every configuration\n");
		wfsolute(ffp, 1);
	}
	if (prFlags.xtr){
		fprintf(ifp,"-Writing .xtr file every %d configurations\n",prFlags.xtr*100);
		wxtr(xfp, 1, prFlags.xtr*100);
	}
	if (prFlags.dat){
		fprintf(ifp,"-Writing .dat file every %d configuration\n",prFlags);
		wdat(dfp,1,prFlags.dat);
	}
if	(point_to_skip == 0)
	{
	H2_avg = H * H;
	H_avg = H;
	T2_avg = temp * temp;
	T_avg = temp;
	}
	else{
		fprintf(ifp,"-Skip %d points in calculating\n",point_to_skip);
		fprintf(ifp,"	average energy and temperature\n");
		H2_avg = H_avg = T2_avg = T_avg = 0.;
	}
	if (prStatQ)
		wstat(stderr);
	fflush(ifp);
/*	Integrate and collect data ...					*/
	i = 0; /*ilan*/
	while (i < numDP /*&& !MeQ*/) { /*ilan*/
#ifdef DEBJ
PREV_E = E;
#endif
		integ(dataRate, fixTQ, 1);
#ifndef DEBJ
	/* print restart files every tenth of the run if the run is long
	 * otherwise print it at the end*/
		if (restart_freq > 10) /* long run - numDP > 100 */
			if (i % restart_freq == 0 && i > 0) {
				sprintf(status, "RST");
				winp(rspfile);
				/*wrsv(rsvfile);*/
			}
#endif
		
		if (prFlags.pos == 1)
			wpos(pfp, 0);
		if (prFlags.pos == 2)
			wshel(pfp, 0);
		if (prFlags.pos == 3 && i%10 == 0 && i > 0)
			wpos(pfp, 0);
		if (prFlags.pos == 4 && nsolute > 0)
			wpsolute(pfp, 0);
		if (prFlags.vel)
			wvel(vfp, 0);
		if (prFlags.frc == 1)
			wfrc(ffp, 0);
		if (prFlags.frc == 4 && nsolute > 0)
			wfsolute(ffp, 0);
		if (prFlags.xtr)
			wxtr(xfp, 0, prFlags.xtr*100);
		if (prFlags.dat)
			wdat(dfp, 0, prFlags.dat);
		if (prStatQ)
			wstat(stderr);
	if (i >= point_to_skip-1)
		{
		H2_avg += H * H;
		H_avg += H;
		T2_avg += temp * temp;
		T_avg += temp;
		}
#ifdef DEBJ
if (abs(E-PREV_E)*KCAL > 10){
fprintf(stderr,"E = %f PREV_E = %f\n",E*KCAL,PREV_E*KCAL);
sprintf(status, "RST");
winp(rspfile);
exit(1);
}
#endif
	i++; /*ilan*/
	}
//	MeQ = 0; // reset flag
	fprintf(stderr,"BEGIN NEW TRAJECTORY\n");
#ifndef DEBJ
	if	(etime > 1000.){ /* print restart files at the end for all
				    runs that are longer than 1 ps. This
				    replaces the previous situation where
				    for long runs (numDP > 100) the final
				    rsp file was not the same as the final
				    position in the .pos file. Begining
				    July 4, 1991 all trajectories will have
				    as a final rsp file a file which is the 
				    same as the final configuration in the
				    .pos file */
		sprintf(status, "RST");
		winp(rspfile);
		/*wrsv(rsvfile);*/
	}
#endif
/*	All done.  Close output files, calculate averages, and return... */

	if (prFlags.pos)
		fclose(pfp);
	if (prFlags.vel)
		fclose(vfp);
	if (prFlags.frc)
		fclose(ffp);
	if (prFlags.xtr)
		fclose(xfp);
	if (prFlags.dat)
		fclose(dfp);

#ifdef RDF
	factor = 2*xwall*ywall*zwall*20*20*20;
	factor /= PI*natoms*natoms*(numDP*dataRate+1);
	fprintf(ifp,"-Solvent part of p.m.f:\n");
	for (i=0; i< 350; i++)
	    {
	    pmf =(rdf[i]/(i*(i+1)+1.))*factor;	
	    if (pmf > 0.)
   	         fprintf(ifp,"%f %f\n",(i+0.5)/20.,-log(pmf) -uu((i+0.5)/20.)/(KBOLTZ*108.));
            }
#endif RDF

	numDP++;
	numDP -=point_to_skip;
	fprintf(ifp,"numDP = %d\n",numDP);
	H_avg /= (double)numDP;
	H2_avg = H2_avg / (double)numDP - H_avg * H_avg;
	T_avg /= (double)numDP;
	T2_avg = T2_avg / (double)numDP - T_avg * T_avg;
	fprintf(ifp, "File = %s\n", outRoot);
	fprintf(ifp, "H_avg = %g	std dev = %g	ratio = %g\n",
		KCAL*H_avg, KCAL*sqrt(H2_avg), sqrt(H2_avg)/H_avg);
	fprintf(ifp, "T_avg = %g	std dev = %g	ratio = %g\n",
		T_avg, sqrt(T2_avg), sqrt(T2_avg)/T_avg);
#ifdef TENSION
/***	pressure tensor calculation		***/
	P_avg = P2_avg = tension = 0.;
	nslabs = 2*zwall/pslabSize;
	convFactor = 3.4*3.4/(KBOLTZ*120.); /* specific for Ar*/
fprintf(ifp,"-----------------------------------\n");
fprintf(ifp,"components of pTensor:\n");
	for (i = 0; i < nslabs; i++) {
		pTensor[i].fy += pTensor[i].fx * KBOLTZ *T_avg ;
		pTensor[i].fy /= ((numDP*dataRate+1)*4*xwall*ywall*pslabSize);
		pTensor[i].fz += pTensor[i].fx * KBOLTZ *T_avg ;
		pTensor[i].fz /= ((numDP*dataRate+1)*4*xwall*ywall*pslabSize);
		pTensor[i].fx = pslabSize*(pTensor[i].fz - pTensor[i].fy ); 
fprintf(ifp,"%f %f %f\n",pTensor[i].fx *convFactor, pTensor[i].fy*convFactor*3.4, pTensor[i].fz*convFactor*3.4);
		tension += pTensor[i].fx/2.;
		P_avg += (2*pTensor[i].fy + pTensor[i].fz)/3.;
		P2_avg += sq(2*pTensor[i].fy + pTensor[i].fz)/9.;
	}
fprintf(ifp,"-----------------------------------\n");
	P_avg /= (double)nslabs;
	P2_avg = P2_avg / (double)nslabs - P_avg * P_avg;
	fprintf(ifp, "Contribution to pressure from LJ forces only:\n");
	fprintf(ifp, "P_avg = %g	std dev = %g	ratio = %g\n",
		P_avg*convFactor*3.4, sqrt(P2_avg)*convFactor*3.4, sqrt(P2_avg)/P_avg);
	fprintf(ifp,"surface tension = %f\n",tension*convFactor);
#endif TENSION
	fclose(ifp);
//printf("end collect.c\n");
}
/* calculate the LJ interaction at the distance r	*/
double uu(r)
double r;
{
	double	ir, ir6, a, aa, b, bb;
	double	r2, pe, f;
	double	s, ur;
	a = lj[0][0].a;
	aa = 12.*a;
	b = lj[0][0].b;
	bb = 6.*b;
        r2 = r*r;
	ir = 1. / r2;
	ir6 = ir * ir * ir;
	if	(r2 <= swr2min)
		ur = ( a * ir6 - b ) * ir6;
	else {
		pe = r2-swr2min;
		f = pe * pe;
		s = swcoef[0] + pe * f *
		(swcoef[1]+pe*swcoef[2]+f*swcoef[3]);
		pe = ( a * ir6 - b ) * ir6;
		ur = s * pe;
		}

	   return(ur);
}
