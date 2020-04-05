#include	<md.h>
#include	<math.h>

/*
 * fix the bond length (real or virtual) between two atoms, or between
 * all pair of a diatomic liquid, depending if the flag contsrain equals
 * 2 or 1 respectively. Do not fix if constrain = 0.
 */

fix(constrain, type, atom1, atom2, d)
int constrain, atom1, atom2;
double d;
char type; 	/* type = r or v for position or velocity adjustments*/
{
int i;
if (! constrain)
	return;
if (constrain == 1)
	/* constrain the bond in all diatomic molecules to equal d.
	 * It is assumed that we have even number of atoms and that the
	 * atoms 2i and 2i+1 are bonded, i = 0,1,...natoms/2 -1.
	 */
	 {
	 if (type == 'r')
		for (i = 0; i < natoms/2; i++)
			adjust_r(2*i,2*i+1,d);
	 if (type == 'v')
		for (i = 0; i < natoms/2; i++)
			adjust_v(2*i,2*i+1,d);
	 }
if (constrain == 2)
	/* use the indeces atom1 and atom2 to constrain only the bond between
	 * these two atoms to be of length d.
	 */
	 {
	 if (type == 'r')
		adjust_r(atom1,atom2,d);
	 if (type == 'v')
		adjust_v(atom1,atom2,d);
	 }
}

adjust_r(i,j,d)
int i,j;
double d;
{ 
tripd rij, vij;
double v2, rv;
double factor;
rij.fx = pos[i].fx - pos[j].fx;
rij.fy = pos[i].fy - pos[j].fy;
rij.fz = pos[i].fz - pos[j].fz;
mvimage(&rij);
vij.fx = vel[i].fx - vel[j].fx;
vij.fy = vel[i].fy - vel[j].fy;
vij.fz = vel[i].fz - vel[j].fz;
v2 = sq(vij.fx) + sq(vij.fy) + sq(vij.fz);
rv = rij.fx * vij.fx + rij.fy * vij.fy + rij.fz * vij.fz;
factor =  h1*h1*rv*rv + d*d*(d*d - h1*h1*v2);
factor = ((sqrt(factor) - h1*rv)/(d*d) -1)/(mass[i]+mass[j]);
pos[i].fx += rij.fx * factor * mass[j];
pos[j].fx -= rij.fx * factor * mass[i];
pos[i].fy += rij.fy * factor * mass[j];
pos[j].fy -= rij.fy * factor * mass[i];
pos[i].fz += rij.fz * factor * mass[j];
pos[j].fz -= rij.fz * factor * mass[i];
}


adjust_v(i,j,d)
int i,j;
double d;
{ 
tripd rij, vij;
double rv, factor;
rij.fx = pos[i].fx - pos[j].fx;
rij.fy = pos[i].fy - pos[j].fy;
rij.fz = pos[i].fz - pos[j].fz;
mvimage(&rij);
vij.fx = vel[i].fx - vel[j].fx;
vij.fy = vel[i].fy - vel[j].fy;
vij.fz = vel[i].fz - vel[j].fz;
rv = rij.fx * vij.fx + rij.fy * vij.fy + rij.fz * vij.fz;
factor = -rv/(d*d*(mass[i] + mass[j]));
vel[i].fx += rij.fx * factor * mass[j];
vel[j].fx -= rij.fx * factor * mass[i];
vel[i].fy += rij.fy * factor * mass[j];
vel[j].fy -= rij.fy * factor * mass[i];
vel[i].fz += rij.fz * factor * mass[j];
vel[j].fz -= rij.fz * factor * mass[i];
}


