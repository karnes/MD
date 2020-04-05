/*
 *	Converts a simple xyz file to an ascii file (without header)
 *
 *	To call:	xyz2as <xyz_file> <ascii_file>
 *
 */

#include <stdio.h>
#define abs(x)          ((x)>0.?(x):-(x))
#define sq(x)   ((x)*(x))
#define hsgn(x)         ((x)<0.?-.5:.5)
#define atoms 77
#define molec 1
#define per 19
#define numUp 8946
#define zUp 24.0

main (argc, argv)
	int argc;
	char *argv[];
{
	int i,j,k;
	FILE	*fp;
	char	sbuf[256],xx;
	double ii,jj,kk;
	double fx[atoms],fy[atoms],fz[atoms];
/*	char gly[14] = {'C','O','H','H',
			'C','O','H','H','H',
			'C','O','H','H','H'};
	int glyID[14]={20,45,7,5,
			20,45,7,5,5,
			20,45,7,5,5};

*/
	char THA[19] = {'C','H','H','H',
			'C','H','H',
			'C','H','H',
			'C','H','H',
			'C','H','H',
			'C','H','H'};
	int THAID[19] = {9,5,5,5,
			 9,5,5,
			 9,5,5,
			 9,5,5,
			 9,5,5,
			 9,5,5};
/*
 *  argv[1] and argv[2] contain the input and output file names
 */
	if (argc != 3) {
		fprintf(stderr, "usage: %s <xyz_file> <ascii_file>\n",
			argv[0]);
		exit(0);
	}

/*
 *  open the binary input file which contains one sn2 system and chloroform molecules
 */
	if ((fp = fopen(argv[1], "r")) == NULL) {
		fprintf(stderr, "can't open xyz input-file\n");
		exit(1);
	}

/*
 *  Read the header and xyz data
 */

	fgets(sbuf,256,fp);
	fgets(sbuf,256,fp);
	for(i=0;i<atoms;i++){
		fgets(sbuf,256,fp);
		if(sscanf(sbuf,"%s%lf%lf%lf",&xx,&ii,&jj,&kk) != 4) {
		   fprintf(stderr, "sysInit: error reading xyz data... i=%d\n",i);
		   exit(1);
		}
		fx[i]=ii;
		fy[i]=jj;
		fz[i]=kk;
	}
/*
 *  Close the binary file
 */
	fclose(fp);

fprintf(stderr,"read xyz...\n");
/*
 *  Write the output to the ascii input-file in argv[2].
 */
	if ((fp = fopen(argv[2], "w")) == NULL) {
		fprintf(stderr, "readbin: can't open ascii input-file\n");
		exit(1);
	}

fprintf(stderr,"opened file...\n");
/*print the atoms*/
// parent, family, children
	i=0;
	fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",i+numUp,55,fx[i],fy[i],fz[i]+zUp,0026,per*(int)(i/per)+numUp,77,76);
	for (i=1; i<atoms; i++){
	    fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",i+numUp,THAID[(i-1)%per],fx[i],fy[i],fz[i]+zUp,0004,0+numUp,0,0);
	}
/*
 *  Close the ascii file, and deallocate the memory
 */
	fclose(fp);
}
