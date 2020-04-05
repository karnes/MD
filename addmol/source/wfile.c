#include	<stdio.h>
#include	<math.h>
#include	"typedefs.h"
#include	"globals.h"

wfile(fp)
	FILE *fp;
{
	static double fieldtemp = 0.0;

	fwrite(filetype, sizeof(char), 5, fp);
	fwrite(&datestamp, sizeof(datestamp), 1, fp);
	fwrite(status, sizeof(char), 4, fp);
	fwrite(&natoms, sizeof(natoms), 1, fp);
	fwrite(&nsolute, sizeof(nsolute), 1, fp);
	fwrite(&xwall, sizeof(xwall), 1, fp);
	fwrite(&ywall, sizeof(ywall), 1, fp);
	fwrite(&zwall, sizeof(zwall), 1, fp);
	fwrite(&EqTemp, sizeof(EqTemp), 1, fp);
	fwrite(&DEqTemp, sizeof(DEqTemp), 1, fp);
	fwrite(&EqPress, sizeof(EqPress), 1, fp);
	fwrite(&DEqPress, sizeof(DEqPress), 1, fp);
	fwrite(atom, sizeof(atom[0]), natoms, fp);
	fwrite(&xtrInQ, sizeof(xtrInQ), 1, fp);
	fwrite(pos, sizeof(pos[0]), natoms, fp);

/* write 6 fields with zeroes in them so that .inp file has same format
as one time step of .pos file */

	fwrite(&fieldtemp, sizeof(fieldtemp), 1, fp);
	fwrite(&fieldtemp, sizeof(fieldtemp), 1, fp);
	fwrite(&fieldtemp, sizeof(fieldtemp), 1, fp);
	fwrite(&fieldtemp, sizeof(fieldtemp), 1, fp);
	fwrite(&fieldtemp, sizeof(fieldtemp), 1, fp);
	fwrite(&fieldtemp, sizeof(fieldtemp), 1, fp);
}
