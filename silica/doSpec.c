/*
 *	doSpec:
 *
 *	This routine controls the running of trajectories that can be stopped
 *      if a certain condition applies
 */

#include	<md.h>
#include        <system.h>

doSpec(fp)
	FILE	*fp;
{
	char	buf[132],	/* input lines are read into this string*/
		endTag[10],	/* end marker				*/
		inFile[80],	/* name of input file			*/
		outFile[80],	/* name of output file			*/
		outRoot[80],	/* root name of output file		*/
		inRoot[80],	/* root name of input file		*/
		conFile[80];	/* system's constants parameters file	*/
	int	numIn,		/* number of initial condition files    */
		numTraj,	/* number of trajectories per init. cond*/
		reset,		/* flag to restore the original input file*/
		numDP,		/* number of data points to collect	*/
		dataRate,	/* number of time steps per data point	*/
		inVers,		/* initial input file number		*/
		outVers,	/* initial output file number		*/
		j,k;		/* dummy indexes			*/
	double	kelvin,		/* equilibriation temperature		*/
		ts;		/* integration time step		*/
	int	rcmQ,		/* remove center mass velocity		*/
		fixTQ,		/* run constant temperature dynamics	*/
		prStatQ,	/* print status line			*/
		vScaleQ;	/* scale initial velocity 		*/
	outflags prFlags;	/* Flags to print output files
				 * if prFlags.dat != 0 print .dat file
				 * if prFlags.pos != 0 print .pos file
				 * if prFlags.vel != 0 print .vel file
				 * if prFlags.frc != 0 print .frc file
				 * if prFlags.xtr != 0 print .xtr file
				 */
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%d %d %d", &numIn, &numTraj, &reset) != 3) {
		fprintf(stderr, "doSpec:error reading 1st line in setup.run\n");
		exit(1);
	}
	reset = 1; /*force the subsequent trajectories to start from the same input file*/
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
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%d %d", &numDP, &dataRate) != 2) {
		fprintf(stderr, "doSpec: error reading numDP, dataRate\n");
		exit(1);
	}
	numDP_x_dataRate = numDP*dataRate;	
	numDPx = numDP;
	dataRatex = dataRate;
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%lf", &ts) != 1) {
		fprintf(stderr, "doSpec: error reading ts\n");
		exit(1);
	}
	settime(ts);

	fgets(buf, 132, fp);
	if (sscanf(buf, "%d %d %d %d", &rcmQ, &fixTQ, &prStatQ, &vScaleQ) != 4) {
		fprintf(stderr," doSpec: error reading rcmQ, fixTQ, prStatQ, vScaleQ\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%d %d %d %d %d", &prFlags.dat, &prFlags.pos, 
		&prFlags.vel, &prFlags.frc, &prFlags.xtr) != 5) {
		fprintf(stderr, "doSpec: error reading output file options\n");
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
	for (j = 0; j < numIn; j++) {
		sprintf(inFile, "%s%03d.equ", inRoot, inVers);
		rinp(inFile);
		init(conFile);
		if (reset == 1)
			save(1);
//		MeQ = 0;/*global flag changed in MSiForce*/
		for (k = 0; k < numTraj; k++) {
			sprintf(outFile, "%s%04d", outRoot, outVers);
			ranvel(kelvin, rcmQ, vScaleQ, 1);
			collectS(numDP, dataRate, outFile,prFlags,prStatQ,fixTQ);
			outVers++;
			if (reset == 1)
			/* start the next trajectory from the same input file*
			 * but different velocities*/
				restore(1);
		}
		inVers++;
	}
}
