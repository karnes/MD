/*
 *	doMin:
 *
 *	This routine handles all the necessary tasks to run a requested 
 *	minimization.  It gets the necessary input data from the setup
 *	file and input file, and then calls init, offset, and minimize.
 *	When done, it calls winp to rewrite the input file.
 */

#include	<md.h>

doMin(fp)
	FILE	*fp;		/* pointer to setup file		*/
{

	char	buf[132],	/* input lines are read into this string*/
		endTag[10],	/* end marker				*/
		inFile[80],	/* name of input file			*/
		outFile[80],	/* name of output file			*/
		outRoot[80],	/* root name of output file		*/
		inRoot[80],	/* root name of input file		*/
		conFile[80];	/* system's constants parameters file	*/
	int	numMin,		/* number of initial files to minimize	*/
		numIter_cg,	/* number of conj. grad. mini. steps	*/
		numIter_sd,	/* number of steep. dec. mini. steps	*/
		niter,		/* number of actual iterations		*/
		vers,		/* initial input file number		*/
		j;		/* dummy index				*/
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%d %d %d", &numMin, &numIter_cg, &numIter_sd) != 3) {
		fprintf(stderr, "doMin: error reading numMin and numIter\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%d", &vers) != 1) {
		fprintf(stderr, "doMin: error reading vers\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%s", inRoot) != 1) {
		fprintf(stderr, "doMin: error reading inRoot\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%s", outRoot) != 1) {
		fprintf(stderr, "doMin: error reading outRoot\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%s", conFile) != 1) {
		fprintf(stderr, "doMin: error reading conFile\n");
		exit(1);
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%s", pbcType) != 1){
		fprintf(stderr, "doTraj: error reading pbcType\n");
		exit(1);
	}

	fgets(buf, 132, fp);
	if (sscanf(buf, "%s", endTag) != 1 ||
		strcmp(endTag, ".ENDMIN.") != 0) {
		fprintf(stderr, "doMin: error reading end marker\n");
		exit(1);
	}
	
	for (j = 0; j < numMin; j++) {
		sprintf(inFile, "%s%03d.new", inRoot, vers);
		sprintf(outFile, "%s%03d.min", outRoot, vers);
		rinp(inFile);
		init(conFile);
/*		offset(natoms-1); */

		niter = numIter_sd;
		niter -= minimize_sd(1.0e-13, niter);
		printf("V = %.8f\n", V*KCAL);
		printf("Number of steepest-descent iterations = %d\n", niter);
		
		niter = numIter_cg;
		niter -= minimize_cg(1.0e-13, niter);
		printf("V = %.8f\n", V*KCAL);
		printf("Number of conjugate-gradient iterations = %d\n", niter);
		
		EqTemp = DEqTemp = 0.0;
		EqPress = DEqPress = 0.0;
		sprintf(status, "MIN");
		winp(outFile);
		
		vers++;
	}
}
