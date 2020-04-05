#include	<md.h>
#include	<system.h>
#include	 <math.h>

SiForce()
{
int	i, j, k, l, m, n, nL;
double	r2, dedr, delfx, delfy, delfz;
double		eg, s, sp;
tripd		frc, image, sdl, f[6];
double		ljqs();
nL = nCH3OH*3+nCH3CN*3;/*total number of liquid atoms*/
for	(i = 0; i < nSi; i++){
	k = nL+i*3;

/*****	Get intramolecular forces for i'th SiOH and SiOO's  ******/
	intraSi(k,nL);
/*****	Get intermolecular forces for i - j SiOH-SiOH interactions.  ******/
	for	(j = i+1; j < nSi; j++,k = nL+i*3){
		l = nL+j*3;
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

/*****	Loop over atoms in SiOH units ******/
		for	(m = 0; m < 3; m++)
			for	(n = 0; n < 3; n++)
				{
/*****	Determine image vector ******/
				delfx = pos[k+m].fx - pos[l+n].fx + image.fx;
				delfy = pos[k+m].fy - pos[l+n].fy + image.fy;
				delfz = pos[k+m].fz - pos[l+n].fz + image.fz;
				r2 = delfx*delfx + delfy*delfy + delfz*delfz;

/*****	Get (1/r dV/dr) ******/
				dedr = ljqs(r2,m,n,&eg,i,j); //added i,j

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
		k = nL + i*3;
		l = nL + j*3;
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
		VNB_S += eg;
		}
	}
}

double
ljqs(r2,m,n,eg,i,j)
double r2, *eg;
int m,n,i,j;
{
double der,r,ir,ir6, ec,a,b,q,aa,bb;

	ir = 1. / r2;
	ir6 = ir * ir * ir;
	a = Silj[m+3][n+3].a;
	b = Silj[m+3][n+3].b;
//	if((i%SiSkip==0 || j%SiSkip==0) && SiSkip < nSi)
//		q=0.0;
//	else
		q = Silj[m+3][n+3].q;
	aa =12*a;
	bb = 6*b;
	r = sqrt(r2);
	ec = q/r;
	*eg +=  ( a * ir6 - b ) * ir6 +ec;
	der = (bb - aa * ir6) * ir6 * ir -ec/r2;
	return(der);
}

intraSi(k,nL)
int k;/*the index of the Oh atom*/
int nL;/*number of liquid atoms in system*/
{
double Si_str, Si_bend, calcstr(), calcbend();
int i;
Si_str = 0;
Si_bend = 0;
/*Oh-H*/
Si_str += calcstr(k,k+2,0);
/*Si-Oh*/
Si_str += calcstr(k,k+1,1);
/*Si-Ob*/
for (i=0;i<3;i++)
	Si_str += calcstr(k+1,nL+(k-nL)/3+(4+i)*nSi,2);
/*H-Oh-Si*/
Si_bend += calcbend(k+2,k,k+1,0);
/*Oh-Si-Ob*/
for (i=0;i<3;i++)
	Si_bend += calcbend(k,k+1,nL+(k-nL)/3+(4+i)*nSi,1);
INTRAV_S += Si_str+Si_bend;
}
double calcstr(i,j,k)
int i,j,k;
{
tripd 	d; 
double	bond, dedr, pe;

d.fx = pos[i].fx - pos[j].fx;
d.fy = pos[i].fy - pos[j].fy;
d.fz = pos[i].fz - pos[j].fz;
bond = sqrt(d.fx*d.fx+d.fy*d.fy+d.fz*d.fz);
dedr = kStr[k]*(bond-eqBond[k])/bond;
pe = 0.5*kStr[k]*sq(bond-eqBond[k]);
force[i].fx -= (d.fx = dedr*d.fx);
force[i].fy -= (d.fy = dedr*d.fy);
force[i].fz -= (d.fz = dedr*d.fz);
force[j].fx += d.fx;
force[j].fy += d.fy;
force[j].fz += d.fz;
//fprintf(stderr,"%d %d %d %f\n",i, j, k, bond);
return(pe);

}
double calcbend(i,j,k,l)
int i,j,k,l;
{
int i0,i1,i2;
tripd r[2], grad[3];
double pe,r1,r2,r1r2,n1n2,da,dedd;

i0 = j;/*center atom*/
i1 = i;
i2 = k;

r[0].fx = pos[i0].fx - pos[i1].fx;
r[0].fy = pos[i0].fy - pos[i1].fy;
r[0].fz = pos[i0].fz - pos[i1].fz;
r[1].fx = pos[i0].fx - pos[i2].fx;
r[1].fy = pos[i0].fy - pos[i2].fy;
r[1].fz = pos[i0].fz - pos[i2].fz;

r1 = sqrt(r[0].fx*r[0].fx + r[0].fy*r[0].fy + r[0].fz*r[0].fz);
r2 = sqrt(r[1].fx*r[1].fx + r[1].fy*r[1].fy + r[1].fz*r[1].fz);
r1r2 =   (r[0].fx*r[1].fx + r[0].fy*r[1].fy + r[0].fz*r[1].fz);
n1n2 = r1r2 / (r1*r2);
da = acos(n1n2) - eqBend[l];
pe = 0.5*kBend[l]*da*da;
dedd = kBend[l]*da/ (sqrt(1.-n1n2*n1n2)*r1*r2);
force[i1].fx -= (grad[1].fx = dedd*(r[1].fx - r1r2*r[0].fx/(r1*r1)));
force[i1].fy -= (grad[1].fy = dedd*(r[1].fy - r1r2*r[0].fy/(r1*r1)));
force[i1].fz -= (grad[1].fz = dedd*(r[1].fz - r1r2*r[0].fz/(r1*r1)));
force[i2].fx -= (grad[2].fx = dedd*(r[0].fx - r1r2*r[1].fx/(r2*r2)));
force[i2].fy -= (grad[2].fy = dedd*(r[0].fy - r1r2*r[1].fy/(r2*r2)));
force[i2].fz -= (grad[2].fz = dedd*(r[0].fz - r1r2*r[1].fz/(r2*r2)));
force[i0].fx += grad[1].fx + grad[2].fx;
force[i0].fy += grad[1].fy + grad[2].fy;
force[i0].fz += grad[1].fz + grad[2].fz;
    return(pe);
}
