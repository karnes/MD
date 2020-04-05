#include	<md.h>
#include	<system.h>
#include	 <math.h>

CH3OHForce()
{
int	i, j, k, l, m, n, index;
double	r2, dedr, delfx, delfy, delfz;
double		eg, s, sp;
tripd		frc, image, sdl, f[6];
double		ljq1();
//printf("CH3OHForce.c\n");
for	(i = 0; i < nCH3OH; i++){
	k = i*3;

/*****	Get intramolecular forces for i'th CH3OH.  ******/
	intraCH3OH(&pos[k],&force[k]);
/*****	Get cos of angles: O-C & O-H, versus +z vector normal to silica.
	Store normalized vectors for TCF */
	if(tc % TCFdt == 0) getMAng(k,i);
//	printf("CH3OHForce.c - just did getMAng\n");
/*****	Get intermolecular forces for i - j CH3OH-CH3OH interactions.  ******/
	for	(j = i+1; j < nCH3OH; j++,k = i*3){
		l = j*3;
		eg = 0.0;

/*****	Determine image vector for O i - O j.  ******/
		image.fx = -(sdl.fx = pos[k].fx - pos[l].fx);
		image.fy = -(sdl.fy = pos[k].fy - pos[l].fy);
		image.fz = -(sdl.fz = pos[k].fz - pos[l].fz);
		mvimage(&sdl);
		image.fx += (delfx =sdl.fx);
		image.fy += (delfy =sdl.fy);
		image.fz += (delfz =sdl.fz);
		r2 = delfx*delfx + delfy*delfy + delfz*delfz;
		if	(r2 >= swr2max)
			continue;
		if	(r2 <= swr2min)
			{
			s = 1.0;
			sp = 0.0;
			}
		else
			swtch(r2-swr2min,&s,&sp);

		for	(m = 0; m < 6; m++)
			f[m].fx = f[m].fy = f[m].fz = 0.;

/*****	Routine for determining H-bonds  *******/
//		printf("checking for hbonds in CH3OHForce.c\n");
		if(r2 < OOdist2 && tc % HBdt == 0){
//			printf("CH3OHForce.c -- heading to getHbonds...\n");
			getHbonds(k,l,delfx,delfy,delfz);
//			printf("CH3OHForce.c - just did getHbonds\n");
			
		}
/*****	Loop over atoms in CH3OH molecules ******/
		for	(m = 0; m < 3; m++)
			for	(n = 0; n < 3; n++)
				{
/*****	Determine image vector ******/
				delfx = pos[k+m].fx - pos[l+n].fx + image.fx;
				delfy = pos[k+m].fy - pos[l+n].fy + image.fy;
				delfz = pos[k+m].fz - pos[l+n].fz + image.fz;
				r2 = delfx*delfx + delfy*delfy + delfz*delfz;
/*****  CH3-CH3 rdf calculations *****/
				if(m==2 && n==2){
					if((pos[k+m].fz < 6.0 && pos[l+n].fz > 6.0) ||
					   (pos[k+m].fz > 6.0 && pos[l+n].fz < 6.0)){
						index=(sqrt(r2))/binRDF;
						if(index<400){
							ccRDF[index]+=1.0;
						}
					}
				}

/*****	Get (1/r dV/dr) ******/
				dedr = ljq1(r2,m,n,&eg);

/*****	Resolve forces on atoms.  ******/
				f[n+3].fx += (delfx *= dedr);
				f[n+3].fy += (delfy *= dedr);
				f[n+3].fz += (delfz *= dedr);
				f[m].fx -= delfx;
				f[m].fy -= delfy;
				f[m].fz -= delfz;
				}
		if	(sp != 0.0)
			{
			for	(m = 0; m < 6; m++)
				{
				f[m].fx *= s;
				f[m].fy *= s;
				f[m].fz *= s;
				}
			f[0].fx -= (frc.fx = sp*eg*sdl.fx);
			f[0].fy -= (frc.fy = sp*eg*sdl.fy);
			f[0].fz -= (frc.fz = sp*eg*sdl.fz);
			f[3].fx += frc.fx;
			f[3].fy += frc.fy;
			f[3].fz += frc.fz;
			eg *= s;
			}
		k = i*3;
		l = j*3;
		force[k].fx += f[0].fx;
		force[k].fy += f[0].fy;
		force[k].fz += f[0].fz;
		force[++k].fx += f[1].fx;
		force[k].fy += f[1].fy;
		force[k].fz += f[1].fz;
		force[++k].fx += f[2].fx;
		force[k].fy += f[2].fy;
		force[k].fz += f[2].fz;

		force[l].fx += f[3].fx;
		force[l].fy += f[3].fy;
		force[l].fz += f[3].fz;
		force[++l].fx += f[4].fx;
		force[l].fy += f[4].fy;
		force[l].fz += f[4].fz;
		force[++l].fx += f[5].fx;
		force[l].fy += f[5].fy;
		force[l].fz += f[5].fz;
		VNB_M += eg;
		}
	}
}

double
ljq1(r2,m,n,eg)
double r2, *eg;
int m,n;
{
double der,r,ir,ir6, ec,a,b,q,aa,bb;

	ir = 1. / r2;
	ir6 = ir * ir * ir;
	a = liqlj[m][n].a;
	b = liqlj[m][n].b;
	q = liqlj[m][n].q;
	aa =12*a;
	bb = 6*b;
	r = sqrt(r2);
	ec = q/r;
	*eg +=  ( a * ir6 - b ) * ir6 +ec;
	der = (bb - aa * ir6) * ir6 * ir -ec/r2;
	return(der);
}


/* intramolecular potential for Methanol*/
/* pos[0] = O, pos[1] = H, pos[2] = CH3*/
intraCH3OH(pos,force)
tripd *pos,*force;
{
double		dedr[3],r[3];
tripdouble	r0,r1,r2,grad[2];
/*OH bond*/
r0.fx = pos[1].fx - pos[0].fx;
r0.fy = pos[1].fy - pos[0].fy;
r0.fz = pos[1].fz - pos[0].fz;
/*CO bond*/
r1.fx = pos[2].fx - pos[0].fx;
r1.fy = pos[2].fy - pos[0].fy;
r1.fz = pos[2].fz - pos[0].fz;

r[0] = sqrt(r0.fx*r0.fx + r0.fy*r0.fy + r0.fz*r0.fz);
r[1] = sqrt(r1.fx*r1.fx + r1.fy*r1.fy + r1.fz*r1.fz);
r[2] = r0.fx*r1.fx+r0.fy*r1.fy+r0.fz*r1.fz;
harmfd(r,dedr);
force[1].fx -= dedr[0]*r0.fx;
force[1].fy -= dedr[0]*r0.fy;
force[1].fz -= dedr[0]*r0.fz;
force[2].fx -= dedr[1]*r1.fx;
force[2].fy -= dedr[1]*r1.fy;
force[2].fz -= dedr[1]*r1.fz;
force[0].fx += dedr[0]*r0.fx + dedr[1]*r1.fx;
force[0].fy += dedr[0]*r0.fy + dedr[1]*r1.fy;
force[0].fz += dedr[0]*r0.fz + dedr[1]*r1.fz;
force[1].fx += (grad[0].fx = dedr[2]*(r1.fx-r[2]*r0.fx/(r[0]*r[0])));
force[1].fy += (grad[0].fy = dedr[2]*(r1.fy-r[2]*r0.fy/(r[0]*r[0])));
force[1].fz += (grad[0].fz = dedr[2]*(r1.fz-r[2]*r0.fz/(r[0]*r[0])));
force[2].fx += (grad[1].fx = dedr[2]*(r0.fx-r[2]*r1.fx/(r[1]*r[1])));
force[2].fy += (grad[1].fy = dedr[2]*(r0.fy-r[2]*r1.fy/(r[1]*r[1])));
force[2].fz += (grad[1].fz = dedr[2]*(r0.fz-r[2]*r1.fz/(r[1]*r[1])));
force[0].fx -= grad[1].fx + grad[0].fx;
force[0].fy -= grad[1].fy + grad[0].fy;
force[0].fz -= grad[1].fz + grad[0].fz;
}

/*
 *	**** function to calculate a triatom harmonic potential		****
 *	**** with respect to the internal coordinates r.      		****
 *      **** r[0] and r[1] are the bond length and r[2]/(r[0]*r[1]) is 	****
 *	**** cos the bond angle						****
 */
#define eqCOH (108.9*PI/180.0)
#define eqOH 0.945
#define eqCO 1.430
#define kCOH (110.0/KCAL)
#define kOH  (1106.0/KCAL)
#define kCO (640.0/KCAL)
harmfd(r,dedr)
double	*r,*dedr;
{
double n1n2, da, q1,q2;
	n1n2 = r[2] /(r[0]*r[1]);
/*	for angles very close to 180 degrees, n1n2 = -1 + eps,
 *	where eps is very small and can be both positive and
 *	negative. If its negative acos will blow up. If its
 *	positive but very small, da will be very small for a linear
 *	equilibrium state and dedr will be the ratio between two
 *	very small quantities. We therefore use a switch	*/
	if (n1n2+1 < 1.0e-8) n1n2 = -1;
	da = acos(n1n2) - eqCOH;
	q1 = r[0]-eqOH;
	q2 = r[1]-eqCO;
	INTRAV_M += 0.5*(kOH*q1*q1+kCO*q2*q2+kCOH*da*da);
	dedr[0] = kOH*q1/r[0];
	dedr[1] = kCO*q2/r[1];
	if (n1n2+1 < 1.0e-8)
		dedr[2] = -kCOH/(r[0]*r[1]);
	else
		dedr[2] = kCOH*da/ (sqrt(1.-n1n2*n1n2)*r[0]*r[1]);
}

void getHbonds(k,l,Ox,Oy,Oz)
int k,l;
double Ox,Oy,Oz;  /* vector from O to O  */
{
	tripd OHk, OHl;
	double rOO, rOHk, rOHl, dotOHOk, dotOHOl;
	/* get O->H (k) vector */
	OHk.fx = pos[k+1].fx - pos[k].fx;
	OHk.fy = pos[k+1].fy - pos[k].fy;
	OHk.fz = pos[k+1].fz - pos[k].fz;  
	/* get O->H (l) vector */
	OHl.fx = pos[l+1].fx - pos[l].fx;
	OHl.fy = pos[l+1].fy - pos[l].fy;
	OHl.fz = pos[l+1].fz - pos[l].fz; 
	/* O(l)->O(k) vector is (Ox,Oy,Oz)  */
	/* where methanol k is H donor */
	dotOHOk = -Ox*OHk.fx - Oy*OHk.fy - Oz*OHk.fz;
	rOHk = sqrt(OHk.fx*OHk.fx + OHk.fy*OHk.fy + OHk.fz*OHk.fz);
	dotOHOl = Ox*OHl.fx + Oy*OHl.fy + Oz*OHl.fz;
	rOHl = sqrt(OHl.fx*OHl.fx + OHl.fy*OHl.fy + OHl.fz*OHl.fz);
	rOO = sqrt(Ox*Ox + Oy*Oy + Oz*Oz);
	if(dotOHOk/(rOO*rOHk) > cosHOO)
	{ 
		/*H-bond donated by methanol k*/
		HBcount[5][(int)((pos[k+1].fz + pos[l].fz)/(2*HBZBinSize))]++;  
	//	printf("assigned Hbond in CH3OHForce\n"); 
	//	printf("donor: %d accept: %d tc: %d z: %f\n",nSi + k/3,nSi + l/3, tc, (pos[k+1].fz+pos[l].fz)/2) ; 
	}
	if(dotOHOl/(rOO*rOHl) > cosHOO)
	{
		/*H-bond donated by methanol l*/
		HBcount[5][(int)((pos[l+1].fz + pos[k].fz)/(2*HBZBinSize))]++;
	//	printf("assigned Hbond in CH3OHForce\n");
	//	printf("donor: %d accept: %d tc: %d z: %f\n",nSi + l/3,nSi + k/3, tc, (pos[l+1].fz+pos[k].fz)/2); 
	}
}
	
void getMAng(int k, int i)
{
//	printf("in getMAng\n");
	double	r, r2, delfx, delfy, delfz;
	double 	cosOC, cosOH;
	
/**Determine cos of angle of O->C vector (versus +z vector normal to silica)
	and store normalized vector for TCF */
	delfx = pos[k+2].fx - pos[k].fx;
	delfy = pos[k+2].fy - pos[k].fy;
	delfz = pos[k+2].fz - pos[k].fz;
	r2 = delfx*delfx + delfy*delfy + delfz*delfz;
	r = sqrt(r2);
	/* normalized O->C vector */
	MOCvec[i][(int)(tc/TCFdt)].fx = (float)(delfx/r);
	MOCvec[i][(int)(tc/TCFdt)].fy = (float)(delfy/r);
	MOCvec[i][(int)(tc/TCFdt)].fz = (float)(delfz/r);//jjk-store, calc TCFs in wxtr.c
	MOCvec[i][(int)(tc/TCFdt)].zpos = (float)(pos[k].fz);
	cosOC = delfz/r;
//	printf("got MOCvec\n");
//	printf("(int)(pos[k].fz/ZBinSize) = %d, (int)((cosOC + 1.0) / CosBinSize) = %d\n", (int)(pos[k].fz/ZBinSize),(int)((cosOC + 1.0) / CosBinSize));
	odOC[(int)(pos[k].fz/ZBinSize)][(int)((cosOC + 1.0) / CosBinSize)]++;
//	printf("got odOC\n");
/**Determine cos of angle of O->H vector (versus +z vector normal to silica)
	and store normalized vector for TCF */
	delfx = pos[k+1].fx - pos[k].fx;
	delfy = pos[k+1].fy - pos[k].fy;
	delfz = pos[k+1].fz - pos[k].fz;
	r2 = delfx*delfx + delfy*delfy + delfz*delfz;
	r = sqrt(r2);
	/* normalized O->H vector */
	MOHvec[i][(int)(tc/TCFdt)].fx = (float)(delfx/r);
	MOHvec[i][(int)(tc/TCFdt)].fy = (float)(delfy/r);
	MOHvec[i][(int)(tc/TCFdt)].fz = (float)(delfz/r); 
	MOHvec[i][(int)(tc/TCFdt)].zpos = (float)(pos[k].fz);//jjk-store, calculate TCFs in wxtr.c
	cosOH = delfz/r;
	odOH[(int)(pos[k].fz/ZBinSize)][(int)((cosOH + 1.0)/CosBinSize)]++;
	int z;
	for(z=0; z < NumSBins; z++)
	{
		if(pos[k].fz > sbins[2 * z] && pos[k].fz < sbins[2 * z + 1])
		{
			sodOC[z][(int)((cosOC + 1.0) / CosBinSize)]++;
			sodOH[z][(int)((cosOH + 1.0) / CosBinSize)]++;
		}
	}
}
