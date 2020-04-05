#include	<md.h>
#include	<system.h>
#include	 <math.h>
#include	<water.h>

waterNITRO(pos,force)
	tripd	*pos;
	tripd	*force;
{
int	i, j, k, l, m, n,nw;
double	r2, dedr, delfx, delfy, delfz;
double		eg, wnc, s, sp;
tripd		frc, image, sdl, f[20];
double		wdljq();

nw = (natoms-nsolute-14*nNIT)/3;
if ( nw == 0 )
	return;  

for	(i = 0; i < nNIT; i++){
	k = nw*3+i*14;

	for	(j = 0; j < nw; j++,k = nw*3+i*14){
		l = j*3;
		eg = wnc = 0.0;

/*****	Determine image vector for primeC i - oxygen j.  ******/
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

		for	(m = 0; m < 20; m++)
			f[m].fx = f[m].fy = f[m].fz = 0.;


		for	(m = 0; m < 14; m++)/*Loop over atoms in NIT molecules*/
			for	(n = 0; n < 3; n++)/*Loop over atoms in H2O*/
				{
/*****	Determine image vector ******/
				delfx = pos[k+m].fx - pos[l+n].fx + image.fx;
				delfy = pos[k+m].fy - pos[l+n].fy + image.fy;
				delfz = pos[k+m].fz - pos[l+n].fz + image.fz;
				r2 = delfx*delfx + delfy*delfy + delfz*delfz;

/*****	Get (1/r dV/dr) ******/
				dedr = wdljq(r2,m,n,&eg,k,&wnc,l);

/*****	Resolve forces on atoms.  ******/
				f[n+14].fx += (delfx *= dedr);
				f[n+14].fy += (delfy *= dedr);
				f[n+14].fz += (delfz *= dedr);
				f[m].fx -= delfx;
				f[m].fy -= delfy;
				f[m].fz -= delfz;
				}
		if	(sp != 0.0)
			{
			for	(m = 0; m < 20; m++)
				{
				f[m].fx *= s;
				f[m].fy *= s;
				f[m].fz *= s;
				}
			f[0].fx -= (frc.fx = sp*eg*sdl.fx);
			f[0].fy -= (frc.fy = sp*eg*sdl.fy);
			f[0].fz -= (frc.fz = sp*eg*sdl.fz);
			f[14].fx += frc.fx;
			f[14].fy += frc.fy;
			f[14].fz += frc.fz;
			eg *= s;
			wnc *= s;
			}
		k = nw*3+i*14;
		l = j*3;
		/* NIT */
		force[k].fx += f[0].fx;
		force[k].fy += f[0].fy;
		force[k].fz += f[0].fz;
		force[++k].fx += f[1].fx;
		force[k].fy += f[1].fy;
		force[k].fz += f[1].fz;
		force[++k].fx += f[2].fx;
		force[k].fy += f[2].fy;
		force[k].fz += f[2].fz;
		force[++k].fx += f[3].fx;
		force[k].fy += f[3].fy;
		force[k].fz += f[3].fz;
		force[++k].fx += f[4].fx;
		force[k].fy += f[4].fy;
		force[k].fz += f[4].fz;
		force[++k].fx += f[5].fx;
		force[k].fy += f[5].fy;
		force[k].fz += f[5].fz;
		force[++k].fx += f[6].fx;
		force[k].fy += f[6].fy;
		force[k].fz += f[6].fz;
		force[++k].fx += f[7].fx;
		force[k].fy += f[7].fy;
		force[k].fz += f[7].fz;
		force[++k].fx += f[8].fx;
		force[k].fy += f[8].fy;
		force[k].fz += f[8].fz;
		force[++k].fx += f[9].fx;
		force[k].fy += f[9].fy;
		force[k].fz += f[9].fz;
		force[++k].fx += f[10].fx;
		force[k].fy += f[10].fy;
		force[k].fz += f[10].fz;
		force[++k].fx += f[11].fx;
		force[k].fy += f[11].fy;
		force[k].fz += f[11].fz;
		force[++k].fx += f[12].fx;
		force[k].fy += f[12].fy;
		force[k].fz += f[12].fz;
		force[++k].fx += f[13].fx;
		force[k].fy += f[13].fy;
		force[k].fz += f[13].fz;
		

		/* water */
		force[l].fx += f[14].fx;
		force[l].fy += f[14].fy;
		force[l].fz += f[14].fz;
		force[++l].fx += f[15].fx;
		force[l].fy += f[15].fy;
		force[l].fz += f[15].fz;
		force[++l].fx += f[16].fx;
		force[l].fy += f[16].fy;
		force[l].fz += f[16].fz;
		H2ONITV += eg;
		WNITC += wnc;
		if(j==0){
		   VINT_NIT += eg;
		   VC_NIT += wnc;
		}
/*		if(fixShel>0){
		  if(j==0||j==fixw[0]||j==fixw[1]){
		        Umer+=eg;
			cMer++;
		  }
		}
*/
		}
	}
}

double
wdljq(r2,m,n,eg,nN,wnc,l)
double r2, *eg,*wnc;
int m, /* m is the NIT atom index */
    n, /* n is the water atom index */
    nN,/* the number of the nitrobenzene molecule for the rdf calculations*/
    l; /* index of water O*/
{
double der,r,ir,ir6, ec,a,b,q,aa,bb;
int nn,bin;
int rindex;
	nn = n + 14;   /* water - NIT cross terms begin at 14 */
	ir = 1. / r2;
	ir6 = ir * ir * ir;
	a = NITlj[m][nn].a;
	b = NITlj[m][nn].b;
	q = NITlj[m][nn].q;
	aa =12*a;
	bb = 6*b;
	r = sqrt(r2);
	if(l==0){
	  if(n==0){
	   rindex = (int)(r/binRDF);
	   if(rindex<300){
	      if(m==0) grSOL[2][rindex]+=1.;
	      if(m==1||m==5) grSOL[3][rindex]+=0.5;
	      if(m==2||m==4) grSOL[4][rindex]+=0.5;
	      if(m==3) grSOL[5][rindex]+=1.;
	      if(m==6) grSOL[6][rindex]+=1.;
	      if(m==7||m==8) grSOL[7][rindex]+=0.5;
	   }
	  }
	}

/* water oxygen - non-hydrogen nitrobenzene rdfs for interfacial molecules
	if (pos[nN].fz < 2 && n == 0 && m<9){
		bin = r/binRDF;
		H2ONITRdf[m][bin] += 1.0;
	}
*/
	ec = q/r;
	*wnc += ec;
	*eg +=  ( a * ir6 - b ) * ir6 + ec;
	der = (bb - aa * ir6) * ir6 * ir - ec/r2 ;
	return(der);
}
