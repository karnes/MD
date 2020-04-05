#include	<md.h>
#include	<math.h>

/*
 * fix the bond length (real or virtual) between two atoms, or between
 * all pairs of a diatomic liquid, depending if the flag contsrain equals
 * 2 or 1 respectively. Do not fix if constrain = 0.
 */

fix(constrain, type, atom1, atom2, d)
int constrain, atom1, atom2;
double d;
char type; 	/* type = r, v or t for pos, vel or KE adjustments*/
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
	 if (type == 't')
		for (i = 0; i < natoms/2; i++)
			adjust_t(2*i,2*i+1,d);
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
	 if (type == 't')
		adjust_t(atom1,atom2,d); /* total KE adjusment	*/
	 }
}

adjust_r(i,j,d)
int i,j;
double d;
{ 
tripd rij, old_rij, oldposi, oldposj, vij;
double r2, rold_r;
double factor;
/*	The unconstrained new positions				*/
rij.fx = pos[i].fx - pos[j].fx;
rij.fy = pos[i].fy - pos[j].fy;
rij.fz = pos[i].fz - pos[j].fz;
mvimage(&rij);
/*	The correct preveous positions				*/
oldposi.fx = pos[i].fx - h1*vel[i].fx;
oldposi.fy = pos[i].fy - h1*vel[i].fy;
oldposi.fz = pos[i].fz - h1*vel[i].fz;
mvimage(&oldposi);
oldposj.fx = pos[j].fx - h1*vel[j].fx;
oldposj.fy = pos[j].fy - h1*vel[j].fy;
oldposj.fz = pos[j].fz - h1*vel[j].fz;
mvimage(&oldposj);
old_rij.fx = oldposi.fx - oldposj.fx;
old_rij.fy = oldposi.fy - oldposj.fy;
old_rij.fz = oldposi.fz - oldposj.fz;
mvimage(&old_rij);

r2 = sq(rij.fx) + sq(rij.fy) + sq(rij.fz);
rold_r = rij.fx * old_rij.fx + rij.fy * old_rij.fy + rij.fz * old_rij.fz;
factor =  sq(rold_r) - d*d*(r2-d*d);
factor = (sqrt(factor) - rold_r)/(d*d*(mass[i]+mass[j]));
pos[i].fx += old_rij.fx * factor * mass[j];
pos[j].fx -= old_rij.fx * factor * mass[i];
pos[i].fy += old_rij.fy * factor * mass[j];
pos[j].fy -= old_rij.fy * factor * mass[i];
pos[i].fz += old_rij.fz * factor * mass[j];
pos[j].fz -= old_rij.fz * factor * mass[i];
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


adjust_t(i,j,d)
int i,j;
double d;
{ 
tripd rij, vij;
double rv, deltaK;
rij.fx = pos[i].fx - pos[j].fx;
rij.fy = pos[i].fy - pos[j].fy;
rij.fz = pos[i].fz - pos[j].fz;
mvimage(&rij);
vij.fx = tvel[i].fx - tvel[j].fx;
vij.fy = tvel[i].fy - tvel[j].fy;
vij.fz = tvel[i].fz - tvel[j].fz;
rv = rij.fx * vij.fx + rij.fy * vij.fy + rij.fz * vij.fz;
deltaK = -rv*rv*mass[i]*mass[j]/(d*d*(mass[i] + mass[j]));
K += deltaK; /* twice the kinetic energy of the system	*/
if	(i < natoms-nsolute && j < natoms-nsolute)
	/* atoms i and j are solvent atoms. must adjust liquid KE */ 
		KLIQ += deltaK/2.;
}


