/*
 *	doTraj:
 *
 *	This routine controls the running of all requested trajectories.
 *	It reads the data in the setup file, reads the input files, does
 *	the needed initializations, calls ranvel, and calls collect to
 *	run each trajectory.
 */

#include	<md.h>

doTraj(fp)
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
		fprintf(stderr, "doTraj:error reading 1st line in setup.run\n");
		exit(1);
	}
	/* check that the new version of the setup file is used*/
	if (reset != 1 && reset != 0){
		fprintf(stderr, "doTraj:error reading 1st line in setup.run\n");
		fprintf(stderr, "check that reset is equal to 0 or 1\n");
		fprintf(stderr, "(new version of setup.run must be used)\n");
		exit(1);
	}
	fgets(buf, 132, fp);
	if (sscanf(buf, "%d %d", &inVers, &outVers) != 2) {
		fprintf(stderr, "doTraj: error reading inVers and outVers\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%s", inRoot) != 1) {
		fprintf(stderr, "doTraj: error reading inRoot\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%s", outRoot) != 1) {
		fprintf(stderr, "doTraj: error reading outRoot\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%lf", &kelvin) != 1) {
		fprintf(stderr, "doTraj: error reading kelvin\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%d %d", &numDP, &dataRate) != 2) {
		fprintf(stderr, "doTraj: error reading numDP, dataRate\n");
		exit(1);
	}
	numDP_x_dataRate = numDP*dataRate;	
	numDPx = numDP;
	dataRatex = dataRate;
	
//	printf("doTraj.c - numDP = %d, dataRate = %d\n",numDP,dataRate);
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%lf", &ts) != 1) {
		fprintf(stderr, "doTraj: error reading ts\n");
		exit(1);
	}
	settime(ts);

	fgets(buf, 132, fp);
	if (sscanf(buf, "%d %d %d %d", &rcmQ, &fixTQ, &prStatQ, &vScaleQ) != 4) {
		fprintf(stderr," doTraj: error reading rcmQ, fixTQ, prStatQ, vScaleQ\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%d %d %d %d %d", &prFlags.dat, &prFlags.pos, 
		&prFlags.vel, &prFlags.frc, &prFlags.xtr) != 5) {
		fprintf(stderr, "doTraj: error reading output file options\n");
		exit(1);
	}
	fgets(buf, 132, fp);
	if (sscanf(buf, "%s", conFile) != 1) {
		fprintf(stderr, "doTraj: error reading conFile\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%s", constrain_file) != 1) {
		fprintf(stderr,"doTraj: error reading constraints file name\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%s", pbcType) != 1){
		fprintf(stderr, "doTraj: error reading pbcType\n");
		exit(1);
	}

	fgets(buf, 132, fp);
	if (sscanf(buf, "%s", endTag) != 1 ||
		strcmp(endTag, ".ENDRUN.") != 0) {
		fprintf(stderr, "doTraj: error reading end marker\n");
		exit(1);
	}
	for (j = 0; j < numIn; j++) {
		sprintf(inFile, "%s%03d.equ", inRoot, inVers);
		rinp(inFile);
		init(conFile);
		if (reset == 1)
			save(1);
		for (k = 0; k < numTraj; k++) {
			sprintf(outFile, "%s%04d", outRoot, outVers);
			ranvel(kelvin, rcmQ, vScaleQ, 1);
			collect(numDP, dataRate, outFile,prFlags,prStatQ,fixTQ);
			outVers++;
			if (reset == 1)
			/* start the next trajectory from the same input file*
			 * but different velocities*/
				restore(1);
		}
		inVers++;
	}
}
