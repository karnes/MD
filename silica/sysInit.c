#include	<md.h>
#include	<system.h>

sysInit(confile)
	char	*confile;
{

FILE	*fp;
char	sbuf[256];
int	i,j, ns, jns;
int	jx, jy, jz;
double 	sig,eps,s6,sqrt(),
 	sigS[3],epsS[3],qS[3],/*O, Si, H*/
 	sigM[3],epsM[3],qM[3],/*O, H, CH3*/
 	sigA[3],epsA[3],qA[3];/*C, N, CH3*/
double	sigI,epsI,qI;/*solute atoms*/

if ((fp = fopen(confile,"r")) == NULL) {
	fprintf(stderr, "sysinit: cannot open %s\n", confile);
	exit(1);
}
/*initialization of the SiOH and SIOO intramolecular parameters*/
kStr[0] = (545.0/KCAL);/* OH force constant*/
eqBond[0] =  0.96; /* OH equilibrium bond length*/
kStr[1] = (428.0/KCAL);/* SiO_h force constant */
eqBond[1] =  1.61; /* SiO_h equilibrium bond length*/
kStr[2] = (885.10/KCAL);/* SiO_b force constant*/
eqBond[2] =  1.61; /* SiO_b equilibrium bond length*/
kBend[0] =(57.5/KCAL) ;/*SiOH bending force constant*/
eqBend[0] = (106.0*PI/180.0);/* SiOH equilibrium bond angle*/
kBend[1] =  (153.26/KCAL);/* SiOO bending force constant*/
eqBend[1] =  (111.09*PI/180.0);/* SiOO equilibrium bond angle*/
	fgets(sbuf,256,fp);
	if (sscanf(sbuf,"%d%d%d%d",&nSi, &nCH3OH, &nCH3CN, &ns) != 4) {
		fprintf(stderr, "sysInit: error reading line 1\n");
		exit(1);
	}
	if (ns > 5){
	  fprintf(stderr, "sysInit: increase nuber of solute atom types in sollj structure array\n");
	  exit(1);
	}

/*	Get the Lennard-Jones paramters for the silica intermolecular forces */
	for (i = 0; i<3; i++){
		fgets(sbuf,256,fp);
		if (sscanf(sbuf,"%lf%lf%lf",&sigS[i], &epsS[i], &qS[i]) != 3) {
			fprintf(stderr, "sysInit: error reading line i+2\n");
			exit(1);
		}
		s6 = sigS[i]*sigS[i]*sigS[i];
		s6 = s6*s6;
		Silj[i+3][i+3].b = 4*epsS[i]*s6/KCAL;
		Silj[i+3][i+3].a = Silj[i+3][i+3].b * s6;
		Silj[i+3][i+3].q = qS[i]*qS[i]/E2;
		gridQ2[i+6] = qS[i]/E2;
	}
/*	calculate mixed interaction parameters*/
	for (i = 0; i<3; i++){
	    for (j = i+1; j<3; j++){
		sig = (sigS[i]+sigS[j])/2.;
		eps = sqrt(epsS[i]*epsS[j]);
		s6 = sig*sig*sig;
		s6 = s6*s6;
		Silj[i+3][j+3].b = Silj[j+3][i+3].b = 4*eps*s6/KCAL;
		Silj[i+3][j+3].a = Silj[j+3][i+3].a = Silj[i+3][j+3].b * s6;
		Silj[i+3][j+3].q = Silj[j+3][i+3].q = qS[i]*qS[j]/E2;
	    }
	}


/*	Get the Lennard-Jones paramters for the CH3OH intermolecular forces */
	for (i = 0; i<3; i++){
		fgets(sbuf,256,fp);
		if (sscanf(sbuf,"%lf%lf%lf",&sigM[i], &epsM[i], &qM[i]) != 3) {
			fprintf(stderr, "sysInit: error reading line i+5\n");
			exit(1);
		}
		s6 = sigM[i]*sigM[i]*sigM[i];
		s6 = s6*s6;
		liqlj[i][i].b = 4*epsM[i]*s6/KCAL;
		liqlj[i][i].a = liqlj[i][i].b * s6;
		liqlj[i][i].q = qM[i]*qM[i]/E2;
		gridQ2[i] = qM[i]/E2;
	}
/*	calculate the CH3OH-CH3OH mixed interaction parameters*/
	for (i = 0; i<3; i++){
	    for (j = i+1; j<3; j++){
		sig = (sigM[i]+sigM[j])/2.;
		eps = sqrt(epsM[i]*epsM[j]);
		s6 = sig*sig*sig;
		s6 = s6*s6;
		liqlj[i][j].b = liqlj[j][i].b = 4*eps*s6/KCAL;
		liqlj[i][j].a = liqlj[j][i].a = liqlj[i][j].b * s6;
		liqlj[i][j].q = liqlj[j][i].q = qM[i]*qM[j]/E2;
	    }
	}
/*	Calculate the Silica-CH3OH mixed interaction parameters*/
	for (i = 0; i<3; i++){
	    for (j = 0; j<3; j++){
		sig = (sigM[i]+sigS[j])/2.;
		eps = sqrt(epsM[i]*epsS[j]);
		s6 = sig*sig*sig;
		s6 = s6*s6;
		Silj[i][j].b = 4*eps*s6/KCAL;
		Silj[i][j].a = Silj[i][j].b * s6;
		Silj[i][j].q = qM[i]*qS[j]/E2;
	    }
	}

/*	Get the Lennard-Jones paramters for the CH3CN intermolecular forces */
	for (i = 0; i<3; i++){
		fgets(sbuf,256,fp);
		if (sscanf(sbuf,"%lf%lf%lf",&sigA[i], &epsA[i], &qA[i]) != 3) {
			fprintf(stderr, "sysInit: error reading line i+8\n");
			exit(1);
		}
		s6 = sigA[i]*sigA[i]*sigA[i];
		s6 = s6*s6;
		liqlj[i+3][i+3].b = 4*epsA[i]*s6/KCAL;
		liqlj[i+3][i+3].a = liqlj[i+3][i+3].b * s6;
		liqlj[i+3][i+3].q = qA[i]*qA[i]/E2;
		gridQ2[i+3] = qA[i]/E2;
	}
/*	calculate the CH3CN-CH3CN mixed interaction parameters*/
	for (i = 0; i<3; i++){
	    for (j = i+1; j<3; j++){
		sig = (sigA[i]+sigA[j])/2.;
		eps = sqrt(epsA[i]*epsA[j]);
		s6 = sig*sig*sig;
		s6 = s6*s6;
		liqlj[i+3][j+3].b = liqlj[j+3][i+3].b = 4*eps*s6/KCAL;
		liqlj[i+3][j+3].a = liqlj[j+3][i+3].a = liqlj[i+3][j+3].b * s6;
		liqlj[i+3][j+3].q = liqlj[j+3][i+3].q = qA[i]*qA[j]/E2;
	    }
	}
/*	calculate the CH3OH-CH3CN mixed interaction parameters*/
	for (i = 0; i<3; i++){
	    for (j = 0; j<3; j++){
		sig = (sigM[i]+sigA[j])/2.;
		eps = sqrt(epsM[i]*epsA[j]);
		s6 = sig*sig*sig;
		s6 = s6*s6;
		liqlj[i][j+3].b = 4*eps*s6/KCAL;
		liqlj[i][j+3].a = liqlj[i][j+3].b * s6;
		liqlj[i][j+3].q = qM[i]*qA[j]/E2;
	    }
	}
/*	Calculate the Silica-CH3CN mixed interaction parameters*/
	for (i = 0; i<3; i++){
	    for (j = 0; j<3; j++){
		sig = (sigA[i]+sigS[j])/2.;
		eps = sqrt(epsA[i]*epsS[j]);
		s6 = sig*sig*sig;
		s6 = s6*s6;
		Silj[i+3][j].b = 4*eps*s6/KCAL;
		Silj[i+3][j].a = Silj[i+3][j].b * s6;
		Silj[i+3][j].q = qA[i]*qS[j]/E2;
	    }
	}

	fgets(sbuf,256,fp);/*read dt values, added by jjk 9/10/14. See system.h*/
	if (sscanf(sbuf,"%d%d%d",&TCFdt, &HBdt, &EPdt) != 3) {
		if(action[0] != 'E' || action[1] != 'E'){
			fprintf(stderr, "sysInit: error reading dt parameters. ");
			fprintf(stderr, "Setting dt's to 10.\n");
		}
			HBdt = TCFdt = EPdt = 10;
	//		printf("%c%c%c%c%c\n", action[0], action[1], action[2], action[3],action[4]);
	}
	else{
		fprintf(stderr,"sysInit: TCFdt = %d, HBdt = %d, EPdt = %d from confile.\n",TCFdt,HBdt,EPdt);
	}
	/* if status = EQU, set dt's to large numbers */
	if(action[0] == 'E' || action[1] == 'E')
	{
	//	fprintf(stderr,"%c%c%c%c\n", action[0], action[1], action[2], action[3]);
		HBdt = TCFdt = EPdt = 1000000;
		fprintf(stderr,"sysInit.c - .EQU. made all dt's = %d\n",TCFdt);
	}
	fgets(sbuf,256,fp);/*read window parameters, added by ilan 9/9/14. See system.h*/
	if (sscanf(sbuf,"%d%lf%lf",&setWindow, &W1, &W2) != 3) {
		fprintf(stderr, "sysInit: error reading window parameters. ");
		fprintf(stderr, "Setting setWindow to zero. Continue without window.\n");
		setWindow = 0;
		W1 = W2 = 0.0;
	}
	else{
		fprintf(stderr, "sysInit: setWindow = %d, W1 = %f, W2 = %f from confile.\n",setWindow,W1,W2);
	}
	fgets(sbuf,256,fp);/* read descriptor line, added by jjk 23-Jan-2015 */
/*	fgets(sbuf,256,fp);*//* read MeOHtau (min time to countwindow parameters, added by ilan 9/9/14. See system.h*/
/*	if(action[0] == 'S' || action[1] == 'S'){   // if special run (non equilib, 1 MeOH in ACN)
		if (sscanf(sbuf,"%lf",&MeOHtau) != 1) {
			fprintf(stderr, "sysInit: error reading MeOHtau. ");
			fprintf(stderr, "Setting MeOHtau to 2.0 ps.\n");
			MeOHtau = 1.99;
		}
		fprintf(stderr,"Special run:  MeOHtau = %f\n",MeOHtau);
	}*/
if	(ns > 0){
	/* Read solute lj parameters	*/
	for (jns =0;jns<ns;jns++){/*loop over all solute atoms types*/
		fgets(sbuf,256,fp);
		if (sscanf(sbuf,"%lf%lf%lf",&sigI, &epsI, &qI) != 3) {
			fprintf(stderr, "sysInit: error reading solute atom %d parameters\n",jns);
			exit(1);
		}
	/*solute-CH3OH*/
	    for (i=0;i<3;i++){
		sig = (sigM[i]+sigI)/2.;
		eps = sqrt(epsM[i]*epsI);
		s6 = sig*sig*sig;
		s6 = s6*s6;
		sollj[jns][i].b = 4*eps*s6/KCAL;
		sollj[jns][i].a = sollj[jns][i].b * s6;
		sollj[jns][i].q = qM[i]*qI/E2;
	    }
	/*solute-CH3CN*/
	    for (i=0;i<3;i++){
		sig = (sigA[i]+sigI)/2.;
		eps = sqrt(epsA[i]*epsI);
		s6 = sig*sig*sig;
		s6 = s6*s6;
		sollj[jns][i+3].b = 4*eps*s6/KCAL;
		sollj[jns][i+3].a = sollj[jns][i+3].b * s6;
		sollj[jns][i+3].q = qA[i]*qI/E2;
	    }
	/*solute-SiOH*/
	    for (i=0;i<3;i++){
		sig = (sigS[i]+sigI)/2.;
		eps = sqrt(epsS[i]*epsI);
		s6 = sig*sig*sig;
		s6 = s6*s6;
		sollj[jns][i+6].b = 4*eps*s6/KCAL;
		sollj[jns][i+6].a = sollj[jns][i+6].b * s6;
		sollj[jns][i+6].q = qS[i]*qI/E2;
	    }
	}
}/*endif ns > 0*/
/*setup the grid for the electrostatic potential calculation*/
for     (jz = 0; jz < NZ; jz++) {
        for (jx = 0; jx < NX; jx++) {
            for (jy = 0; jy < NY; jy++) {
                j = NX*NY*jz+NY*jx+jy;
//		potGrid[j].fx = -xwall+(jx+0.5)/(NX*1.0);
		potGrid[j].fx = 0.999*(-xwall+(2.0*xwall)*((double)jx)/((double)NX));
//		potGrid[j].fy = -ywall+(jy+0.5)/(NY*1.0);
		potGrid[j].fy = 0.999*(-ywall+(2.0*ywall)*((double)jy)/((double)NY));
		potGrid[j].fz = 2.6+jz*0.25;
    	    }
	}
}
//for(j=0;j<NX*NY*NZ;j++)  fprintf(stderr,"%f\t%f\t%f\n",potGrid[j].fx,potGrid[j].fy,potGrid[j].fz);
fclose(fp);
}
