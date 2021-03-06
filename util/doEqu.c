/*
 *	doEqu:
 *
 *	This routine handles a requested equilibration run.  It is called
 *	by main.  The routine gets the necessary input data from the setup
 *	file, and makes the necessary calls to init, offset, and equil 
 *	to get things going.
 */

#include	<md.h>

void doEqu(fp)
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
		numEq,		/* number of equil files per input file	*/
		numRand,	/* number of data points to collect	*/
		numSteps,	/* number of time steps per data point	*/
		inVers,		/* initial input file number		*/
		outVers,	/* initial output file number		*/
		j,k;		/* dummy indexes			*/
	double	kelvin,		/* equilibriation temperature		*/
		ts;		/* integration time step		*/
	int	rcmQ,		/* remove center mass velocity		*/
		ofsetQ,		/* set atom number natoms - ofsetQ	*
				 * to position 0,0,0			*/
		fixTQ,		/* run constant temperature dynamics	*/
		prStatQ,	/* print status line			*/
		vScaleQ; 	/* scale initial velocity		*/

	fgets(buf, 132, fp);
	if (sscanf(buf, "%d %d", &numIn, &numEq) != 2) {
		fprintf(stderr, "doEqu: error reading numIn and numEq\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%d %d", &inVers, &outVers) != 2) {
		fprintf(stderr, "doEqu: error reading inVers and outVers\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%s", inRoot) != 1) {
		fprintf(stderr, "doEqu: error reading inRoot\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%s", outRoot) != 1) {
		fprintf(stderr, "doEqu: error reading outRoot\n");
		exit(1);
	}
	fgets(buf, 132, fp);
	if (sscanf(buf, "%lf", &kelvin) != 1) {
		fprintf(stderr, "doEqu: error reading kelvin\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%d %d", &numRand, &numSteps) != 2) {
		fprintf(stderr, "doEqu: error reading numRand, numSteps\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%d %d %d %d", &fixTQ, &prStatQ, &vScaleQ, &ofsetQ) != 4) {
		fprintf(stderr,"doEqu:error reading fixTQ,prStatQ, vScaleQ and ofsetQ\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%lf", &ts) != 1) {
		fprintf(stderr, "doEqu: error reading ts\n");
		exit(1);
	}
	settime(ts);
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%s", conFile) != 1) {
		fprintf(stderr, "doEqu: error reading conFile\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%s", constrain_file) != 1) {
		fprintf(stderr,"doEqu: error reading constraints file name\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%s", pbcType) != 1){
		fprintf(stderr, "doTraj: error reading pbcType\n");
		exit(1);
	}

	fgets(buf, 132, fp);
	if (sscanf(buf, "%s", endTag) != 1 ||
		strcmp(endTag, ".ENDEQU.") != 0) {
		fprintf(stderr, "doEqu: error reading end marker\n");
		exit(1);
	}
	for (j = 0; j < numIn; j++) {
		sprintf(inFile, "%s%03d.min", inRoot, inVers);
		rinp(inFile);
		init(conFile);

		for (k = 0; k < numEq; k++) {
			sprintf(outFile, "%s%03d.equ", outRoot, outVers);
			equil(numRand, numSteps, kelvin, fixTQ,prStatQ,vScaleQ);
			EqTemp = kelvin;
			DEqTemp = 0.0;
			sprintf(status, "EQU");
			if (ofsetQ) offset(natoms-ofsetQ);
			winp(outFile);
			outVers++;
		}
		
		inVers++;
	}
}
