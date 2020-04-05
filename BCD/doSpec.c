/*
 *	doSpec
 *
 *	This routine read input rsp files where the sn2 system is near the TS
 *	it brings it to the TS and equilibrate than let it react and calculate
 *	the reactive flux. Repeat for different velocities and requilibrate 
 *	and do again.
 */

#include	<md.h>
#include	<system.h>

void doSpec(fp)
	FILE	*fp;
{
	char	buf[132],	/* input lines are read into this string*/
		endTag[10],	/* end marker				*/
		inFile[80],	/* name of input file			*/
		outFile[80],	/* name of output file			*/
		outRoot[80],	/* root name of output file		*/
		inRoot[80],	/* root name of input file		*/
		conFile[80];	/* system's constants parameters file	*/
	int	numIn,		/* number of input files		*/
		ofsetQ,		/* set last atom to 0,0,0		*/
		numEq,		/* number of equil files per input file	*/
		numTraj,	/*number of dynamic trajectories per equ*/
		numRand,	/* number of data points to collect	*/
		numSteps,	/* number of time steps per data point	*/
		numDP,		/*data points per trajectory  		*/
		dataRate,	/*integration steps per data point	*/
		inVers,		/* initial input file number		*/
		outVers,	/* initial output file number		*/
		i,j,k,jj;	/* dummy indexes			*/
	double	kelvin,		/* equilibriation temperature		*/
		tsR, tsD;	/* integration time step		*/
	int	
		fixTQ,		/* run constant temperature dynamics	*/
		prStatQ,	/* print status line			*/
		vScaleQ; 	/* scale initial velocity		*/
	outflags prFlags;
	double denom;

	fgets(buf, 132, fp);
	if (sscanf(buf, "%d %d %d", &numIn, &numEq, &numTraj) != 3) {
		fprintf(stderr, "doSpec: error reading numIn, numEq and numTraj\n");
		exit(1);
	}
	if (numIn*numEq*numTraj > 1000){	
		fprintf(stderr, "doSpec: too many trajectories\n");
		exit(1);
	}
	fgets(buf, 132, fp);
	if (sscanf(buf, "%d %d", &inVers, &outVers) != 2) {
		fprintf(stderr, "doSpec: error reading inVers and outVers\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%s", inRoot) != 1) {
		fprintf(stderr, "doSpec: error reading inRoot\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%s", outRoot) != 1) {
		fprintf(stderr, "doSpec: error reading outRoot\n");
		exit(1);
	}
	fgets(buf, 132, fp);
	if (sscanf(buf, "%lf", &kelvin) != 1) {
		fprintf(stderr, "doSpec: error reading kelvin\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);/*for the equilibration part*/
	if (sscanf(buf, "%d %d", &numRand, &numSteps) != 2) {
		fprintf(stderr, "doSpec: error reading numRand, numSteps\n");
		exit(1);
	}

	fgets(buf, 132, fp);/*for the dynamics*/
        if (sscanf(buf, "%d %d", &numDP, &dataRate) != 2) {
                fprintf(stderr, "doSpec: error reading numDP, dataRate\n");
                exit(1);
        }

	fgets(buf, 132, fp);
	if (sscanf(buf, "%d %d %d %d", &fixTQ, &prStatQ, &vScaleQ, &ofsetQ) != 4) {
		fprintf(stderr,"doSpec:error reading fixTQ,prStatQ, vScaleQ and ofsetQ\n");
		exit(1);
	}

	fgets(buf, 132, fp);
        if (sscanf(buf, "%d %d %d %d %d", &prFlags.dat, &prFlags.pos,
                &prFlags.vel, &prFlags.frc, &prFlags.xtr) != 5) {
                fprintf(stderr, "doSpec: error reading output file options\n");
                exit(1);
        }
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%lf %lf", &tsR, &tsD) != 2) {
		fprintf(stderr, "doSpec: error reading ts for equ and dynamics\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%s", conFile) != 1) {
		fprintf(stderr, "doSpec: error reading conFile\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%s", constrain_file) != 1) {
		fprintf(stderr,"doSpec: error reading constraints file name\n");
		exit(1);
	}

	fgets(buf, 132, fp);
	if (sscanf(buf, "%s", pbcType) != 1){
		fprintf(stderr, "doSpec: error reading pbcType\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%s", endTag) != 1 ||
		strcmp(endTag, ".ENDSPL.") != 0) {
		fprintf(stderr, "doSpec: error reading end marker\n");
		exit(1);
	}
	fluxQ = 0;/*flag to instruct flux calculations in sysForce*/
endC = 0.0;
for (i=0;i<numDP*dataRate+1;i++){
	corflux[i] = 0;
	totc1sq[i] = 0.0;
	totc[i] = 0.0;
}
for (j = 0; j < numIn; j++) {
	sprintf(inFile, "%s%04d.rsp", inRoot, inVers);
	rinp(inFile);
	init(conFile);
	settime(0.01);
	x1RC= x2RC = 0.;
	kRC = 1000.;
/*short equilibration to bring the system to the TS*/
	fprintf(stderr,"initial  equilibration %d\n",j);
	equil(50, 10, kelvin, fixTQ,1,vScaleQ);
	for (k = 0; k < numEq; k++) {
		settime(tsR);
		x1RC= x2RC = 0.;
		kRC = 1000.;
/*long equilibration to prepare the next independet configuration*/
		fprintf(stderr,"long  equilibration %d %d\n",j,k);
		sprintf(outFile,"%s%03d%02db.equ",outRoot,outVers,k);
		equil(numRand, numSteps, kelvin, fixTQ,prStatQ,vScaleQ);
		if (abs(xwall - zwall) < 0.1 && ofsetQ) /*bulk*/
			offset(natoms-3);
		winp(outFile);
		settime(0.01);
		kRC = 10000.;
/*short equilibration to bring system to the TS exactly*/

		fprintf(stderr,"short final equilibration %d %d\n",j,k);
		equil(10, 10, kelvin, fixTQ,1,vScaleQ);
		settime(tsD);
		save(1);
		kRC = 0.;
		fprintf(stderr,"begin flux calculations\n");
		fluxQ = 1;
	//	plat = 9.0; // release BCDz constraint
		for (i=0;i<numTraj;i++){
			for(jj=0;jj<numDP*dataRate+1;jj++)
			   cprod[jj]=0.0;
			fprintf(stderr,"trajectory = %d %d %d\n",j,k,i);
			sprintf(outFile, "%s%04d", outRoot, outVers);
			ranvel(kelvin, 1, vScaleQ, 1);
			collect(numDP, dataRate, outFile,prFlags,1,fixTQ);
			outVers++;
			/* make cprod go to products */
			for(jj=0;jj<numDP*dataRate+1;jj++){
			   if(cprod[numDP*dataRate] < 0.5){
			      totc[jj] += 1.0 - cprod[jj];
			      endC += 2.0;
			   } 
			   else{
			      totc[jj] += cprod[jj];
			      endC += 1.0;
			   }
			}
			restore(1);
		}
		fluxQ = 0;
	//	plat = 0.1;
	}
	inVers++;
}
/*print normalized flux correlation*/
denom = numTraj*numIn*numEq*1.0;
for (i=1;i<numDP*dataRate+1;i++)
	fprintf(stderr,"%f\t%f\t%f\t%f\t%f\n",i*tsD,corflux[i]/denom,totc1sq[i]/denom,totc[i]/denom,endC/(denom*(float)(numDP*dataRate+1)));
}
