/*
 * takes a binary input file that begins with N MeOH molecules
 * and mutates these N MeOH into EtOH molecules and creates an
 * input file with E EtOH molecules.
 */

#include	<stdio.h>
#include	<typedefs.h>
#include	<globals.h>
#include	<time.h>

main (argc, argv)
	int	argc;
	char	*argv[];
{
	int	i, k, N, E;
//, RFlag;
//	double atof(),Z;
	FILE	*fp;

/*
 *  argv[1] and argv[2] contain the input and output file names
 */
	if (argc != 5) {
		fprintf(stderr, "usage: %s <binary_file> <ascii_file> <N> <E>\n",
			argv[0]);
		exit(0);
	}

N = atoi(argv[3]);
E = atoi(argv[4]);

	if (E> N) {
		fprintf(stderr, "usage: %s <bin_file> <as_file> <N> <E>... N >= E.\n",
			argv[0]);
		exit(0);
	}
//Z = atof(argv[4]);
//RFlag = atoi(argv[5]);
tripd HO, HCH3, EtCH3;
double normHO, normHCH3, normEtCH3;
double wHO, wHCH3, dist, offset;
/*
 *  open the binary input file
 */
	if ((fp = fopen(argv[1], "r")) == NULL) {
		fprintf(stderr, "readbin: can't open binary input-file\n");
		exit(1);
	}

/*
 *  Read the header data, allocate space for the arrays, and read them too.
 */
	fread(filetype, sizeof(char), 5, fp);
	fread(&datestamp, sizeof(datestamp), 1, fp);
	fread(status, sizeof(char), 4, fp);
	fread(&natoms, sizeof(natoms), 1, fp);
	fread(&nsolute, sizeof(nsolute), 1, fp);
	fread(&xwall, sizeof(xwall), 1, fp);
	fread(&ywall, sizeof(ywall), 1, fp);
	fread(&zwall, sizeof(zwall), 1, fp);
	fread(&EqTemp, sizeof(EqTemp), 1, fp);
	fread(&DEqTemp, sizeof(DEqTemp), 1, fp);
	fread(&EqPress, sizeof(EqPress), 1, fp);
	fread(&DEqPress, sizeof(DEqPress), 1, fp);

	if ((atom = (parts *) calloc(natoms, sizeof(atom[0]))) == NULL ||
	    (pos  = (tripd *) calloc(natoms, sizeof(pos[0]))) == NULL) {
		fprintf(stderr, "readbin: out of memory\n");
		exit(1);
	}

	fread(atom, sizeof(parts), natoms, fp);
	fread(&xtrInQ, sizeof(xtrInQ), 1, fp);
	fread(pos, sizeof(tripd), natoms, fp);

/*
 *  Close the binary file
 */
	fclose(fp);

/*
 *  Write the output to the ascii input-file in argv[2]
 */
	if ((fp = fopen(argv[2], "w")) == NULL) {
		fprintf(stderr, "readbin: can't open ascii input-file\n");
		exit(1);
	}

	fprintf(fp,"Input File -- \t%s\n",argv[2]);
	fprintf(fp,"Filetype  -- \t.inp\n");

/*  ctime() inserts a newline character so we don't need to do it */
	fprintf(fp,"Datestamp -- \t%s", ctime(&datestamp));

	fprintf(fp, "Status -- \t%s\n", status);
	fprintf(fp,"Number of atoms (natoms) -- \t%d\n", natoms - (3*N) + (4*E)); /* 1 additional atom per MeOH->EtOH */
	fprintf(fp,"Number of solute atoms (nsolute) -- \t%d\n", nsolute);
	fprintf(fp,"X box size (xwall) -- \t%f\n", xwall);
	fprintf(fp,"Y box size (ywall) -- \t%f\n", ywall);
	fprintf(fp,"Z box size (zwall) -- \t%f\n", zwall);
	fprintf(fp, "Equilibration Temperature (EqTemp) -- \t%f\n",
		EqTemp);
	fprintf(fp, "St. Dev. of Equil. Temp. (DEqTemp) -- \t%f\n",
		DEqTemp);
	fprintf(fp, "Equilibration Pressure (EqTemp) -- \t%f\n",
		EqPress);
	fprintf(fp, "St. Dev. of Equil. Press. (DEqPress) --\t%f\n",
		DEqPress);
	fprintf(fp, "Extra file flag (xtrInQ) -- \t%d\n", xtrInQ);

	fprintf(fp,
"   #  type  x position    y position    z position  flags parent param1 param2\n");
/*for (i=0; i<natoms; i++) {
	fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
i,atom[i].type,pos[i].fx,pos[i].fy,pos[i].fz,atom[i].flags,atom[i].parent,atom[i].param1,atom[i].param2);
	}*/

	/** loop over MeOH atoms **/
	for(i=0; i < E; i++)
	{
		k = i * 3;
		/** get MeOH H->CH3 vector, normalize **/
		HCH3.fx = pos[k+2].fx - pos[k+1].fx;
		HCH3.fy = pos[k+2].fy - pos[k+1].fy;
		HCH3.fz = pos[k+2].fz - pos[k+1].fz;
		normHCH3 = sqrt(HCH3.fx*HCH3.fx + HCH3.fy*HCH3.fy + HCH3.fz*HCH3.fz);
		HCH3.fx /= normHCH3;
		HCH3.fy /= normHCH3;
		HCH3.fz /= normHCH3;

		/** get MeOH H->O vector, normalize **/
		HO.fx = pos[k].fx - pos[k+1].fx;
		HO.fy = pos[k].fy - pos[k+1].fy;
		HO.fz = pos[k].fz - pos[k+1].fz;
		normHO = sqrt(HO.fx*HO.fx + HO.fy*HO.fy + HO.fz*HO.fz);
		HO.fx /= normHO;
		HO.fy /= normHO;
		HO.fz /= normHO;

		/** 'bisect' the angle, normalize **/
		wHCH3 = 44.5;
		wHO = 55.5;
		EtCH3.fx = wHCH3 * HCH3.fx + wHO * HO.fx;
		EtCH3.fy = wHCH3 * HCH3.fy + wHO * HO.fy;
		EtCH3.fz = wHCH3 * HCH3.fz + wHO * HO.fz;
		normEtCH3 = sqrt(EtCH3.fx*EtCH3.fx + EtCH3.fy*EtCH3.fy + EtCH3.fz*EtCH3.fz);
		EtCH3.fx /= normEtCH3;
		EtCH3.fy /= normEtCH3;
		EtCH3.fz /= normEtCH3;

		/**  lengthen, add to H position **/
		dist = 3.22442;
		EtCH3.fx *= dist;
		EtCH3.fy *= dist;
		EtCH3.fz *= dist;
		EtCH3.fx += pos[k+1].fx;
		EtCH3.fy += pos[k+1].fy;
		EtCH3.fz += pos[k+1].fz;

		/** re-order atoms: 0=CH2, 1=CH3, 2=O, 3=H **/
		/* CH2 is MeOH's CH3 */
		offset = 0; /*pos[k].fz * (46/32) - pos[k].fz;*/
		fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
			4*i,20/*CH2*/,pos[k+2].fx,pos[k+2].fy,pos[k+2].fz + offset,atom[k].flags,4*i,3,4);
		/* CH3 is the new atom */
		fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
			4*i+1,atom[k+2].type,EtCH3.fx,EtCH3.fy,EtCH3.fz + offset,atom[k+1].flags,4*i,atom[k+1].param1,atom[k+1].param2);
		/* O */
		fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
			4*i+2,atom[k].type,pos[k].fx,pos[k].fy,pos[k].fz + offset,atom[k+1].flags,4*i,atom[k+1].param1,atom[k+1].param2);
		/* H */
		fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
			4*i+3,atom[k+1].type,pos[k+1].fx,pos[k+1].fy,pos[k+1].fz + offset,atom[k+1].flags,4*i,
			atom[k+1].param1,atom[k+1].param2);
	}

for (i=N*3; i< natoms; i++) {
	fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",
	i-(3*N)+(4*E),atom[i].type,pos[i].fx,pos[i].fy,pos[i].fz,atom[i].flags,atom[i].parent - (3*N) + (4*E),atom[i].param1,atom[i].param2);
	}

/*
 *  Close the ascii file, and deallocate the memory
 */
	fclose(fp);
	cfree(atom);
	cfree(pos);
}
