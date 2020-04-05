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
#define atoms (12*125)
#define waters 986
#define molec 4
#define per 19 
#define zUp 0.0

main (argc, argv)
    int argc;
    char *argv[];
{
    int i,j,k;
    FILE	*fp;
    char	sbuf[256],xx;
    double ii,jj,kk;
    double fx[atoms],fy[atoms],fz[atoms];
    double ffx[waters*3],ffy[waters*3],ffz[waters*3];
    double Ex,Ey,Ez;
    int THAID[19]={9,5,5,5,
	9,5,5,
	9,5,5,
	9,5,5,
	9,5,5,
	9,5,5};

    int H2OID[3]={0,1,1};
    int NID = 55; 

    /*
     *  argv[1] and argv[2] contain the input and output file names
     */
    if (argc != 4) {
	fprintf(stderr, "usage: %s <THAxyz_file> <ErWxyz_file> <ascii_file>\n",
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

    /*
     *  open the binary input file which contains one sn2 system and chloroform molecules
     */
    if ((fp = fopen(argv[2], "r")) == NULL) {
	fprintf(stderr, "can't open xyz input-file\n");
	exit(1);
    }

    fgets(sbuf,256,fp);
    fgets(sbuf,256,fp);
    for(i=0;i<waters*3;i++){
	fgets(sbuf,256,fp);
	if(sscanf(sbuf,"%s%lf%lf%lf",&xx,&ii,&jj,&kk) != 4) {
	    fprintf(stderr, "sysInit: error reading xyz data... i=%d\n",i);
	    exit(1);
	}
	ffx[i]=ii;
	ffy[i]=jj;
	ffz[i]=kk;
    }
    fgets(sbuf,256,fp);
    if(sscanf(sbuf,"%s%lf%lf%lf",&xx,&ii,&jj,&kk) != 4) {
	fprintf(stderr, "sysInit: error reading xyz data... i=%d\n",i);
	exit(1);
    }
    Ex=ii;
    Ey=jj;
    Ez=kk;
    /*
     *  Close the binary file
     */
    fclose(fp);

    fprintf(stderr,"read xyz...\n");
    /*
     *  Write the output to the ascii input-file in argv[2].
     */
    if ((fp = fopen(argv[3], "w")) == NULL) {
	fprintf(stderr, "readbin: can't open ascii input-file\n");
	exit(1);
    }

    fprintf(stderr,"opened file...\n");
    /*print the atoms*/
    // parent, family, children
    for(k=0;k<waters;k++){
	for(j=0;j<3;j++){
	    i=k*3+j;
	    fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",i,H2OID[j],ffx[i],ffy[i],ffz[i],(j==0?0026:0004),k*3,(j==0?2:0),(j==0?3:0));
	}
    }
//print nitrogen parent
    fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",waters*3,NID,Nx,Ny,Nz,0026,waters*3,76,77);
    for(k=0;k<molec;k++){
	for(j=0;j<per;j++){
	    i=per*k+j+1;

	    fprintf(fp,"%4d %3d  %12.5f  %12.5f  %12.5f   %04o  %4d %4d %4d\n",i+waters*3,THAID[j],fx[i],fy[i],fz[i],0004,waters*3,0,0);
	}
    }
    /*
     *  Close the ascii file, and deallocate the memory
     */
    fclose(fp);
}
