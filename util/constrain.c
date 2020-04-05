#include	<md.h>
#include	<math.h>

/*
 * This is a general constraint package.
 * It first read the information specific to the given system (what bond
 * to constrain, what length and what method), which is located
 * in the run directory for that system under the name indicated in
 * the setup file, and stored in the global variable constrain_file,
 * then does the position velocities or KE adjustments according
 * to the flag "type". There is also a provision for more specific constraints.
 */

constrain(type)
char type; 	/* type = r, v or t for pos, vel or KE adjustments*/
{
FILE *fp;
char buf[132];
int i;
static first_call = 1;
if	(first_call) {
	
/*	read information about the system	*/
	if	((fp = fopen(constrain_file,"r")) == NULL){
		 fprintf(stderr,"cannot open file %s\n", constrain_file);
		 exit(1);
	}

	fgets(buf, 132, fp);
	if (sscanf(buf, "%s", cons_method) != 1) {
		fprintf(stderr, "constrain: error reading cons_method\n");
		exit(1);
	}
	if (strcmp(cons_method, ".UNC.") == 0) {
		fclose(fp);
		first_call = 0;
		return; /* no constrains */
	}
	
	fgets(buf, 132, fp);
	if (sscanf(buf, "%s", cons_info) != 1) {
		fprintf(stderr, "constrain: error reading cons_info\n");
		exit(1);
	}

		/* read number of bond constrains */
	fgets(buf, 132, fp);
	if (sscanf(buf, "%d", &cons_num) != 1) { 
		fprintf(stderr, "constrain: error reading cons_num\n");
		exit(1);
	}
	if	(cons_info[0] == 'l') {
	/* normal list of pair of atom numbers and distance to constrain*/
		/* allocate space */

	if ( (cons_list   = (constr *) calloc(cons_num,sizeof(constr))) == NULL)
		ERROR((stderr,"constrain: out of core\n"), exit);
	
		/* read atom indices and the bond length */

		for (i = 0; i < cons_num; i++){
			fgets(buf, 132, fp);
			if (sscanf(buf, "%d %d %lf", &cons_list[i].atomA, &cons_list[i].atomB, &cons_list[i].dis) != 3) {
				fprintf(stderr, "constrain: error reading \n");
				exit(1);
			}
		}
	}
	if	(cons_info[0] == 'a') {
	/* assume that we have a diatomic liquids and all pairs are const. */
		fgets(buf, 132, fp);
		if (sscanf(buf, "%lf", &cons_bond) != 1) {
			fprintf(stderr, "constrain: error reading cons_bond\n");
			exit(1);
		}
	
        }
	if	(cons_info[0] == 's') {
	/* special list of pair of atom numbers and distance to constrain*/
		/* allocate space */

	if ( (cons_list   = (constr *) calloc(cons_num,sizeof(constr))) == NULL)
		ERROR((stderr,"constrain: out of core\n"), exit);
	
		/* read atom indices and the bond length */

		for (i = 0; i < cons_num; i++){
			fgets(buf, 132, fp);
			if (sscanf(buf, "%d %d %lf", &cons_list[i].atomA, &cons_list[i].atomB, &cons_list[i].dis) != 3) {
				fprintf(stderr, "constrain: error reading \n");
				exit(1);
			}
		}
	}
	fclose(fp);
	first_call = 0;
}
if (strcmp(cons_method, ".UNC.") == 0) {
	return; /* no constrains */
}
if (strcmp(cons_method, ".ANA.") == 0) {
	/* analytical shake method. applicable to diatomic molecules only */
	shake_a(type);
}

if (strcmp(cons_method, ".ITR.") == 0) {
	/* General iteration over constraints method */
	shake_i(type);
}
}
shake_a(type)
char type;
{
int i;
if (cons_info[0] == 'a')
	/* constrain the bond in all diatomic molecules to equal cons_bond.
	 * It is assumed that the 2*num_cons atoms are in the beginning of the
	 * list and that atoms 2i and 2i+1 are bonded, i = 0,1,...num_cons -1.
	 */
	 {
	 if (type == 'r')
		for (i = 0; i < cons_num; i++)
			adjust_r(2*i,2*i+1,cons_bond);
	 if (type == 'v')
		for (i = 0; i < cons_num; i++)
			adjust_v(2*i,2*i+1,cons_bond);
	 if (type == 't')
		for (i = 0; i < cons_num; i++)
			adjust_t(2*i,2*i+1,cons_bond);
	 }
if (cons_info[0] == 'l')
	/* use the global array cons_list to constrain only the bond between
	 * atoms cons_list[i].atomA, cons_list[i].atomB to be cons_list[i].dis,
	 * where i = 0 - num_cons-1. It is assumed that no atom belongs to more
	 * than one constraint.
	 */
	 {
	 if (type == 'r')
		for (i = 0; i < cons_num; i++)
			adjust_r(cons_list[i].atomA, cons_list[i].atomB, cons_list[i].dis);
	 if (type == 'v')
		for (i = 0; i < cons_num; i++)
			adjust_v(cons_list[i].atomA, cons_list[i].atomB, cons_list[i].dis);
	 if (type == 't')
		for (i = 0; i < cons_num; i++)
			adjust_t(cons_list[i].atomA, cons_list[i].atomB, cons_list[i].dis);
	 }
	if (cons_info[0] == 's') {
		fprintf(stderr,"can't do special constraints using\n");
		fprintf(stderr,"analytical method. Must use ITR\n");
		exit(1);
	}
}

shake_i(type)
char type;
{
int i;
	if (cons_info[0] == 'a') {
		fprintf(stderr,"For a diatomic liquids, use only\n");
		fprintf(stderr,"analytical method.\n");
		exit(1);
	}
	if (cons_info[0] == 'l') { /* simple list of bond constraints*/
	/*do normal shake using iteration method over constraints*/
	 if (type == 'r')
		iter_r(90000,1.e-7);
	 if (type == 'v')
		iter_v(90000,1.e-7);
	 if (type == 't')
		iter_t(90000,1.e-7);
	}
	if (cons_info[0] == 's') { /* special list of constraints*/
	/*do normal shake using iteration method on special constraints*/
	/* call a special routine which is system specific */
		fprintf(stderr,"Iteration method not implemented yet.\n");
		exit(1);
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
/*	The correct previous positions				*/
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
/*	Use the iteration over constraints method to correct the 
 *	unconstrained positions
 */
iter_r(iter,tol)
int iter;
double tol;	
{ 
tripd rAB;
double discr, denom, oldSigmaj, gamaj;
int A,B,i,j;
static int	allocate = 1;
/**************************************************
 *** ALLOCATE SPACE IF FIRST TIME THROUGH	***
 **************************************************/
if (allocate &&
(  (prevpos	= (tripd  *) calloc(natoms,sizeof(tripd))) == NULL
|| (dsigma	= (tripd *) calloc(cons_num,sizeof(tripd))) == NULL
)) {
	ERROR((stderr,"constrain:iter_r: out of core\n"), exit);
	}
allocate = 0;	/* so we won't allocate space the next time */
/**************************************************
 ***	           INITIALIZATION		***
 **************************************************/
/*	find previous-step positions off all atoms */
for	(i = 0; i < natoms; i++)
	{
	prevpos[i].fx = pos[i].fx - h1*vel[i].fx;
	prevpos[i].fy = pos[i].fy - h1*vel[i].fy;
	prevpos[i].fz = pos[i].fz - h1*vel[i].fz;
	mvimage(&prevpos[i]);
	}
/*	calculate the derivatives of the constraints with respect to atom pos *
 *	Here, since we consider bond distances only, each constraint involve  *
 *	only two atoms. We put in array sigma the position vector between     *
 *	these two atoms divided by twice the mass of each atom		      */
for	(j = 0; j < cons_num; j++)
	{
	dsigma[j].fx = prevpos[cons_list[j].atomA].fx
	             - prevpos[cons_list[j].atomB].fx;
	dsigma[j].fy = prevpos[cons_list[j].atomA].fy
	             - prevpos[cons_list[j].atomB].fy;
	dsigma[j].fz = prevpos[cons_list[j].atomA].fz
	             - prevpos[cons_list[j].atomB].fz;
	mvimage(&dsigma[j]);
	}
/**************************************************
 ***	MAIN PART OF PROGRAM			***
 **************************************************/
/*	Outermost loop over iterations		*/
discr = 1.;
while	(iter-- && discr > 0.)
	{
	discr = 0.;
/*	Loop over all constraints		*/
	for	(j = 0; j < cons_num; j++)
		{
		A = cons_list[j].atomA;
		B = cons_list[j].atomB;
		rAB.fx = pos[A].fx-pos[B].fx;
		rAB.fy = pos[A].fy-pos[B].fy;
		rAB.fz = pos[A].fz-pos[B].fz;
		mvimage(&rAB);
		/*	calculate the constraint using old positions	*/
		oldSigmaj = rAB.fx*rAB.fx + rAB.fy*rAB.fy + rAB.fz * rAB.fz -
				cons_list[j].dis * cons_list[j].dis ;
		if (abs(oldSigmaj) >= tol) 
		/*	if convergence criterion is not satisfied for	*
		 *	this constraint update maximum deviation and	*
	 	 *	update positions of the two atoms.		*/
			{
			if (discr < abs(oldSigmaj)) discr = abs(oldSigmaj);
			denom = 1./mass[A] + 1./mass[B];
			denom *= (rAB.fx*dsigma[j].fx + rAB.fy*dsigma[j].fy +
				 rAB.fz*dsigma[j].fz) ;
			gamaj = - oldSigmaj/denom;
			pos[A].fx += gamaj*dsigma[j].fx/(2*mass[A]);
			pos[A].fy += gamaj*dsigma[j].fy/(2*mass[A]);
			pos[A].fz += gamaj*dsigma[j].fz/(2*mass[A]);
			pos[B].fx -= gamaj*dsigma[j].fx/(2*mass[B]);
			pos[B].fy -= gamaj*dsigma[j].fy/(2*mass[B]);
			pos[B].fz -= gamaj*dsigma[j].fz/(2*mass[B]);
			}
		}
	}/*	End while.	*/
if(discr > 0.)
	{
	fprintf(stderr,"convergence problem, largest deviation of square ");
	fprintf(stderr,"bond length is %f\n",discr);
	exit(1);
	}
}
/*	Use the iteration over constraints method to correct the 
 *	unconstrained velocities
 */
iter_v(iter,tol)
int iter;
double tol;	
{ 
tripd vAB;
double discr, denom, VdotR, gamaj;
int A,B,j;
static int	allocate = 1;
/**************************************************
 *** ALLOCATE SPACE IF FIRST TIME THROUGH	***
 **************************************************/
if (allocate &&
   ((dsigma	= (tripd *) calloc(cons_num,sizeof(tripd))) == NULL
)) {
	ERROR((stderr,"constrain:iter_v: out of core\n"), exit);
	}
allocate = 0;	/* so we won't allocate space the next time */
/**************************************************
 ***	           INITIALIZATION		***
 **************************************************/
/*	calculate the derivatives of the constraints with respect to atom pos *
 *	Here, since we consider bond distances only, each constraint involve  *
 *	only two atoms. We put in array sigma the position vector between     *
 *	these two atoms 						      */
for	(j = 0; j < cons_num; j++)
	{
	dsigma[j].fx = pos[cons_list[j].atomA].fx
	             - pos[cons_list[j].atomB].fx;
	dsigma[j].fy = pos[cons_list[j].atomA].fy
	             - pos[cons_list[j].atomB].fy;
	dsigma[j].fz = pos[cons_list[j].atomA].fz
	             - pos[cons_list[j].atomB].fz;
	mvimage(&dsigma[j]);
	}
/**************************************************
 ***	MAIN PART OF PROGRAM			***
 **************************************************/
/*	Outermost loop over iterations		*/
discr = 1.;
while	(iter-- && discr > 0.)
	{
	discr = 0.;
/*	Loop over all constraints		*/
	for	(j = 0; j < cons_num; j++)
		{
		A = cons_list[j].atomA;
		B = cons_list[j].atomB;
		vAB.fx = vel[A].fx-vel[B].fx;
		vAB.fy = vel[A].fy-vel[B].fy;
		vAB.fz = vel[A].fz-vel[B].fz;
		
		/*	calculate the component of vAB along rAB	*/
		VdotR = vAB.fx*dsigma[j].fx + vAB.fy*dsigma[j].fy +
				vAB.fz*dsigma[j].fz;
		if (abs(VdotR) >= tol) 
		/*	if convergence criterion is not satisfied for	*
		 *	this constraint update maximum deviation and	*
	 	 *	update velocities of the two atoms.		*/
			{
			if (discr < abs(VdotR)) discr = abs(VdotR);
			denom = (1./mass[A] + 1./mass[B]) * cons_list[j].dis *
				 cons_list[j].dis;
			gamaj =  VdotR/denom;
			vel[A].fx -= gamaj*dsigma[j].fx/mass[A];
			vel[A].fy -= gamaj*dsigma[j].fy/mass[A];
			vel[A].fz -= gamaj*dsigma[j].fz/mass[A];
			vel[B].fx += gamaj*dsigma[j].fx/mass[B];
			vel[B].fy += gamaj*dsigma[j].fy/mass[B];
			vel[B].fz += gamaj*dsigma[j].fz/mass[B];
			}
		}
	}/*	End while.	*/
if(discr > 0.)
	{
	fprintf(stderr,"convergence problem, largest component of velocity ");
	fprintf(stderr,"along bond constraint is %f\n",discr);
	exit(1);
	}
}
/*	Use the iteration over constraints method to correct the 
 *	unconstrained kinetic energies.
 */
iter_t(iter,tol)
int iter;
double tol;	
{ 
tripd vAB;
double discr, denom, VdotR, gamaj;
int A,B,i,j;
static int	allocate = 1;
/**************************************************
 *** ALLOCATE SPACE IF FIRST TIME THROUGH	***
 **************************************************/
if (allocate &&
   ((dsigma	= (tripd *) calloc(cons_num,sizeof(tripd))) == NULL
   )) {
	ERROR((stderr,"constrain:iter_v: out of core\n"), exit);
	}
allocate = 0;	/* so we won't allocate space the next time */
/**************************************************
 ***	           INITIALIZATION		***
 **************************************************/
/*	calculate the derivatives of the constraints with respect to atom pos *
 *	Here, since we consider bond distances only, each constraint involve  *
 *	only two atoms. We put in array sigma the position vector between     *
 *	these two atoms 						      */
for	(j = 0; j < cons_num; j++)
	{
	dsigma[j].fx = pos[cons_list[j].atomA].fx
	             - pos[cons_list[j].atomB].fx;
	dsigma[j].fy = pos[cons_list[j].atomA].fy
	             - pos[cons_list[j].atomB].fy;
	dsigma[j].fz = pos[cons_list[j].atomA].fz
	             - pos[cons_list[j].atomB].fz;
	mvimage(&dsigma[j]);
	}
/**************************************************
 ***	MAIN PART OF PROGRAM			***
 **************************************************/
/*	Outermost loop over iterations		*/
discr = 1.;
while	(iter-- && discr > 0.)
	{
	discr = 0.;
/*	Loop over all constraints		*/
	for	(j = 0; j < cons_num; j++)
		{
		A = cons_list[j].atomA;
		B = cons_list[j].atomB;
		vAB.fx = tvel[A].fx-tvel[B].fx;
		vAB.fy = tvel[A].fy-tvel[B].fy;
		vAB.fz = tvel[A].fz-tvel[B].fz;
		
		/*	calculate the component of vAB along rAB	*/
		VdotR = vAB.fx*dsigma[j].fx + vAB.fy*dsigma[j].fy +
				vAB.fz*dsigma[j].fz;
		if (abs(VdotR) >= tol) 
		/*	if convergence criterion is not satisfied for	*
		 *	this constraint update maximum deviation and	*
	 	 *	update velocities of the two atoms.		*/
			{
			if (discr < abs(VdotR)) discr = abs(VdotR);
			denom = (1./mass[A] + 1./mass[B]) * cons_list[j].dis *
				 cons_list[j].dis;
			gamaj =  VdotR/denom;
			tvel[A].fx -= gamaj*dsigma[j].fx/mass[A];
			tvel[A].fy -= gamaj*dsigma[j].fy/mass[A];
			tvel[A].fz -= gamaj*dsigma[j].fz/mass[A];
			tvel[B].fx += gamaj*dsigma[j].fx/mass[B];
			tvel[B].fy += gamaj*dsigma[j].fy/mass[B];
			tvel[B].fz += gamaj*dsigma[j].fz/mass[B];
			}
		}
	}/*	End while.	*/
if(discr > 0.)
	{
	fprintf(stderr,"convergence problem, largest component of velocity ");
	fprintf(stderr,"along bond constraint is %f\n",discr);
	exit(1);
	}
/* at this point we can correct the kinetic energy	*/
K = 0;
for	(i = 0; i < natoms; i++)
	{
	K += mass[i] * (tvel[i].fx * tvel[i].fx +
			tvel[i].fy * tvel[i].fy +
			tvel[i].fz * tvel[i].fz);
	if (i == natoms - nsolute -1)
	   /* we just finish the accumulation of
	    * the kinetic energy of the solvent atoms*/
	    KLIQ = K/2.;
	}
}
