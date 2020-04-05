/*
 *	reStart:
 *
 *	This routine controls the running of trajectories which need to be
 *	restarted. It reads the data in the setup file, reads the input 
 *	position (.rsp file) and velocities (.rsv file), does the needed
 *	initializations, call ranvel (optional) and calls collect to run each
 *	trajectory.
 *	
 */

#include	<md.h>

reStart(fp)
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
				 * if prFlags.dat == 1 print .dat file
				 * if prFlags.pos == 1 print .pos file
				 * if prFlags.vel == 1 print .vel file
				 * if prFlags.frc == 1 print .frc file
				 * if prFlags.xtr == 1 print .xtr file
				 */
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%d %d", &numIn, &numTraj) != 2) {
		fprintf(stderr, "reStart: error reading numIn and numTraj\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%d %d", &inVers, &outVers) != 2) {
		fprintf(stderr, "reStart: error reading inVers and outVers\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%s", inRoot) != 1) {
		fprintf(stderr, "reStart: error reading inRoot\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%s", outRoot) != 1) {
		fprintf(stderr, "reStart: error reading outRoot\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%lf", &kelvin) != 1) {
		fprintf(stderr, "reStart: error reading kelvin\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%d %d", &numDP, &dataRate) != 2) {
		fprintf(stderr, "reStart: error reading numDP, dataRate\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%lf", &ts) != 1) {
		fprintf(stderr, "reStart: error reading ts\n");
		exit(1);
	}
	settime(ts);

	fgets(buf, 132, fp);
	if (sscanf(buf, "%d %d %d %d", &rcmQ, &fixTQ, &prStatQ, &vScaleQ) != 4) {
		fprintf(stderr," reStart: error reading rcmQ, fixTQ, prStatQ, vScaleQ\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%d %d %d %d %d", &prFlags.dat, &prFlags.pos, 
		&prFlags.vel, &prFlags.frc, &prFlags.xtr) != 5) {
		fprintf(stderr, "reStart: error reading output file options\n");
		exit(1);
	}
	fgets(buf, 132, fp);
	if (sscanf(buf, "%s", conFile) != 1) {
		fprintf(stderr, "reStart: error reading conFile\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%s", constrain_file) != 1) {
		fprintf(stderr,"reStart: error reading constraints file name\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%s", pbcType) != 1){
		fprintf(stderr, "doTraj: error reading pbcType\n");
		exit(1);
	}

	fgets(buf, 132, fp);
	if (sscanf(buf, "%s", endTag) != 1 ||
		strcmp(endTag, ".ENDSTR.") != 0) {
		fprintf(stderr, "reStart: error reading end marker\n");
		exit(1);
	}
	for (j = 0; j < numIn; j++) {
		sprintf(inFile, "%s%03d.rsp", inRoot, inVers);
		rinp(inFile);
		init(conFile);
		sprintf(inFile, "%s%03d.rsv", inRoot, inVers);
		rvel(inFile,1);
		save(1);
		for (k = 0; k < numTraj; k++) {
			sprintf(outFile, "%s%04d", outRoot, outVers);
			ranvel(kelvin, rcmQ, vScaleQ, k);
			collect(numDP, dataRate, outFile,prFlags,prStatQ,fixTQ);
			outVers++;
			restore(1);
		}
		inVers++;
	}
}
