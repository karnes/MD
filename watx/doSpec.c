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

//	void checkRand(void); /* verify that randomly generated velocities move reaction
//				 coordinate toward product */

//	void prinRanVel(void); // debugging routine: print out selected results from ranvel.c
			       // and post checkRand velocities */
	fgets(buf, 132, fp);
	if (sscanf(buf, "%d %d %d", &numIn, &numTraj, &reset) != 3) {
		fprintf(stderr, "doSpec:error reading 1st line in setup.run\n");
		exit(1);
	}
//	reset = 1; /*force the subsequent trajectories to start from the same input file*/
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
	//for creating TS  
	strcpy(fileHead,outRoot);  // for TS files
	fileVers = outVers;    // for TS files
//fprintf(stderr,"filehead = %s\n",fileHead);	
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
//			printRanVel();
//			checkRand();
//			printRanVel();
//			collectS(numDP, dataRate, outFile,prFlags,prStatQ,fixTQ);
			outVers++;
			fileVers = outVers;  // for TS files
			if (reset == 1)
			/* start the next trajectory from the same input file*
			 * but different velocities*/
				restore(1);
		}
		inVers++;
	}
}

void checkRand(void)
{
/*
	int i,j;
	int mo,sio;
	int MO,SiO,N; // indices of relevant atoms:Si as accept/donor
	double dist;
	double rab,rbc;
	double ab,bc;
	double d,dxdt;
	tripd Ra,Rb,Rc,Va,Vb,Vc;

	//set index modifiers based on TS of inpfile
	if(inpTS==0)
		return;
	else if(inpTS==1)//Xa; Si acceptor
	{
		mo=1; //MeOH's H
		sio=0; //SiOH's O
	}
	else if(inpTS==2)//Xd; Si donor
	{
		mo=0; //MeOH's O
		sio=2; //SiOH's H
	}
	else
	{
		fprintf(stderr,"sys error: doSpec.c--don't recognize inpTS. No ranvel modification.\n");
		return;
	}
	MO = 0;
// find shortest MeOH(H/O)-SiOH(O/H) distances
	d = dist = 999.0;
	for(i=0;i<nSi;i++)
	{
		j=3*nCH3OH+3*nCH3CN+3*i; // SiOH (O)
		d = (pos[MO+mo].fx-pos[j+sio].fx)*(pos[MO+mo].fx-pos[j+sio].fx)+
		   (pos[MO+mo].fy-pos[j+sio].fy)*(pos[MO+mo].fy-pos[j+sio].fy)+
		   (pos[MO+mo].fz-pos[j+sio].fz)*(pos[MO+mo].fz-pos[j+sio].fz);
		if(d<dist){
			dist=d;
			SiO=j;
		}
	}
//fprintf(stderr,"relevant Si(O) = %d, distance = %f, x(SiO) = %f, y(SiO)= %f, z(SiO) = %f  \n",(SiO-3*nCH3OH-3*nCH3CN)/3,sqrt(dist),pos[SiO].fx,pos[SiO].fy,pos[SiO].fz);
	d = dist = 999.0;
// SiO-=2; //for debug
// find shortest SiOH(O)-ACN(N) distances
	for(i=0;i<nCH3CN;i++)
	{	
		j=3*nCH3OH+3*i+1;  // ACN(N)
		d = (pos[SiO].fx-pos[j].fx)*(pos[SiO].fx-pos[j].fx)+
			(pos[SiO].fy-pos[j].fy)*(pos[SiO].fy-pos[j].fy)+
			(pos[SiO].fz-pos[j].fz)*(pos[SiO].fz-pos[j].fz);
		if(d<dist){
			dist=d;
			N=j;
		}
	}
//fprintf(stderr,"relevant ACN(N) = %d, distance = %f, x(N) = %f, y(N) = %f, z(N) = %f\n",(N-1-nCH3OH*3)/3,sqrt(dist),pos[N].fx,pos[N].fy,pos[N].fz);
// calculate dx/dt for rxn coordinate
	Ra=pos[MO];
	Rb=pos[SiO];
	Rc=pos[N];
	Va=vel[MO];
	Vb=vel[SiO];
	Vb=vel[N];
	rab=sqrt((Rb.fx-Ra.fx)*(Rb.fx-Ra.fx)+
			(Rb.fy-Ra.fy)*(Rb.fy-Ra.fy)+
			(Rb.fz-Ra.fz)*(Rb.fz-Ra.fz));
	rbc=sqrt((Rc.fx-Rb.fx)*(Rc.fx-Rb.fx)+
			(Rc.fy-Rb.fy)*(Rc.fy-Rb.fy)+
			(Rc.fz-Rb.fz)*(Rc.fz-Rb.fz));
	ab = (Rb.fx-Ra.fx)*(Vb.fx-Va.fx)+
		(Rb.fy-Ra.fy)*(Vb.fy-Va.fy)+
		(Rb.fz-Ra.fz)*(Vb.fz-Va.fz);
	bc = (Rc.fx-Rb.fx)*(Vc.fx-Vb.fx)+
		(Rc.fy-Rb.fy)*(Vc.fy-Vb.fy)+
		(Rc.fz-Rb.fz)*(Vc.fz-Vb.fz);
	dxdt = (ab/rab)-(bc/rbc);
// if dxdt > 0 then flip all velocities
	if(dxdt > 0.0)
	{
		for(i=0;i<natoms;i++)
		{
			vel[i].fx *= -1.0;
			vel[i].fy *= -1.0;
			vel[i].fz *= -1.0;
		}
		fprintf(stderr,"system dxdt = %f -- reversed ranvel velocities.\n",dxdt);
	}
	else
	{
		fprintf(stderr,"system dxdt = %f -- accepted ranvel velocities.\n",dxdt);
	}
*/
}

void printRanVel(void){
/*
  	int i,j,k;
	for(i=0;i<nCH3OH;i++)
	{
	  fprintf(stderr,"MeOH(O) v(x)=%-5.4f,v(y)=%-5.4f,v(z)=%-5.4f\n",vel[i+0].fx,vel[i+0].fy,vel[i+0].fz);
	  fprintf(stderr,"MeOH(H) v(x)=%-5.4f,v(y)=%-5.4f,v(z)=%-5.4f\n",vel[i+1].fx,vel[i+1].fy,vel[i+1].fz);
	  fprintf(stderr,"MeOH(M) v(x)=%-5.4f,v(y)=%-5.4f,v(z)=%-5.4f\n",vel[i+2].fx,vel[i+2].fy,vel[i+2].fz);
	}
	for(i=0;i<nCH3CN;i=i+100)
	{
	  fprintf(stderr,"ACN (C) v(x)=%-5.4f,v(y)=%-5.4f,v(z)=%-5.4f\n",vel[i+0].fx,vel[i+0].fy,vel[i+0].fz);
	  fprintf(stderr,"ACN (N) v(x)=%-5.4f,v(y)=%-5.4f,v(z)=%-5.4f\n",vel[i+1].fx,vel[i+1].fy,vel[i+1].fz);
	  fprintf(stderr,"ACN (M) v(x)=%-5.4f,v(y)=%-5.4f,v(z)=%-5.4f\n",vel[i+2].fx,vel[i+2].fy,vel[i+2].fz);
	}
	for(i=0;i<nSi;i=i+10)
	{
	  fprintf(stderr,"SiOH(O) v(x)=%-5.4f,v(y)=%-5.4f,v(z)=%-5.4f\n",vel[i+0].fx,vel[i+0].fy,vel[i+0].fz);
	  fprintf(stderr,"SiOS(S) v(x)=%-5.4f,v(y)=%-5.4f,v(z)=%-5.4f\n",vel[i+1].fx,vel[i+1].fy,vel[i+1].fz);
	  fprintf(stderr,"SiOH(H) v(x)=%-5.4f,v(y)=%-5.4f,v(z)=%-5.4f\n",vel[i+2].fx,vel[i+2].fy,vel[i+2].fz);
	}
*/
}

