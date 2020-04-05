/*
 *
 *	tsAngs2.c:	Reads a1ne (1-ACN) transition state files and outputs 
 *			an orientational distribution of the C-N vectors
 *
 *	To call:	tsang2 //<binary_file> <ascii_file>
 *
 */
#include	<stdio.h>
#include	<stdlib.h>
#include	<typedefs.h>
#include	<globals.h>
#include	<time.h>

#define cosBins 10

main (/*argc, argv*/)
/*	int	argc;
	char	*argv[]; */
{
	int	i,j,k;
	int flag;
	char fileName[80];
	FILE	*fp;
	int ODs[10][cosBins] = { 0 };
	int nODs[10] = { 0 };
	double dx,dy,dz,len;

	const char *head[5];
	head[0] = "hb";
	head[1] = "Xd";
	head[2] = "Ma";
	head[3] = "Md";
	head[4] = "Mt";
	const char *head2[5];
	head2[0] = "ACN";
	head2[1] = "ts";
	head2[2] = "hb";
	head2[3] = "hb";
	head2[4] = "hb";


for(k=0;k<1;k++){  
	j=0;
	flag=1;
	sprintf(fileName,"../../NEsilica2/run/inpfiles/%sStart/%s%s%03d.equ",head[k],head[k],head2[k],j);
	if ((fp = fopen(fileName, "r")) == NULL) {
		flag=0;
		fprintf(stderr,"Failed to open %s\n",fileName);
	}
	else
		fprintf(stderr,"Opened %s\n",fileName);
	while(flag==1){
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
//fprintf(stderr,"0=%d, 1=%d, 2=%d, 3=%d, 3063=%d, 3066=%d, 3067=%d, 3068=%d, 3069=%d\n", atom[0].type, atom[1].type, atom[2].type, atom[3].type, atom[3063].type, atom[3066].type,atom[3067].type,atom[3068].type,atom[3069].type);
	//tabulate vectors
	//C->N vector
	dx=pos[3067].fx-pos[3066].fx;
	dy=pos[3067].fy-pos[3066].fy;
	dz=pos[3067].fz-pos[3066].fz;
	len=sqrt(dx*dx+dy*dy+dz*dz);
	
	ODs[k][(int)(((double)cosBins)/2.0*dz/len+((double)cosBins)/2.0)]++;
	nODs[k]++;

	//O->C vector
//	dx=pos[1].fx-pos[0].fx;
//	dy=pos[1].fy-pos[0].fy;
//	dz=pos[1].fz-pos[0].fz;
//	len=sqrt(dx*dx+dy*dy+dz*dz);
//	
//	ODs[k+5][(int)(((double)cosBins)/2.0*dz/len+((double)cosBins)/2.0)]++;
//	nODs[k+5]++;

	//open next file
	j++;
	sprintf(fileName,"../../NEsilica2/run/inpfiles/%sStart/%s%s%03d.equ",head[k],head[k],head2[k],j);
	if ((fp = fopen(fileName, "r")) == NULL) {
		flag=0;
	}
	else
		fprintf(stderr,"Opened %s\n",fileName);
	}
}
	fp=fopen("ACNhbTS_ODs.txt","w");
	fprintf(fp,"cos\t");
	for(i=0;i<1;i++)
		fprintf(fp,"%s_CN_%d\t",head[i],nODs[i]);
//	for(i=0;i<5;i++)
//		fprintf(fp,"%s_OH_%d\t",head[i],nODs[i+5]);
	fprintf(fp,"\n");

	for(i=0;i<cosBins;i++){
		fprintf(fp,"%f\t",((2.0*(double)(i))/((double)cosBins))-1.0);
		for(j=0;j<1;j++)
			fprintf(fp,"%f\t",(double)ODs[j][i]/(double)nODs[j]);
		fprintf(fp,"\n");
	}
	fclose(fp);

}
