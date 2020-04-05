#include	<md.h>
#include	<system.h>
#include	 <math.h>
#include	<water.h>

waterBrO()
{
int	i, j, k, l, m, n, nw, d;
double	r2, dedr, delfx, delfy, delfz;
double	eg, s, sp;
tripd	frc, image, sdl, f[BrOs+3];
double	wblj();
double	wbq();

nw = (natoms-nBCD*BCDs-nBrO*BrOs - nsolute)/3;

for(i = 0; i < nBrO; i++){
   k = nw*3+i*BrOs;

   for(j = 0; j < nw; j++,k = nw*3+i*BrOs){
 	l = j*3;

/*****	L-J interaction, Br-octane - water.  ******/
/*****	Determine image vector for BrO atom center(i) - water O(j)  ******/
	   for(m = 0; m < BrOs; m++)/*Loop over atoms in BrOct*/
	   {
	      k = nw*3+i*BrOs+m;
	      l = j*3;
	      image.fx = -(sdl.fx = pos[k].fx - pos[l].fx);
	      image.fy = -(sdl.fy = pos[k].fy - pos[l].fy);
	      image.fz = -(sdl.fz = pos[k].fz - pos[l].fz);
	      mvimage(&sdl);
	      image.fx += (delfx =sdl.fx);
	      image.fy += (delfy =sdl.fy);
	      image.fz += (delfz =sdl.fz);
	      r2 = delfx*delfx + delfy*delfy + delfz*delfz;
	      if(r2 >= swr2max)
	   	continue;
	      if(r2 <= swr2min)
	      {
	   	s = 1.0;
		sp = 0.0;
	      }
	      else
		swtch(r2-swr2min,&s,&sp);

	      for(d = 0; d < 4; d++)
		f[d].fx = f[d].fy = f[d].fz = 0.;
	      eg = 0.0;

	      for(n = 0; n < 3; n++)/*Loop over atoms in water molecules*/
	      {
/*****	Determine image vector ******/
		delfx = pos[k].fx - pos[l+n].fx + image.fx;
		delfy = pos[k].fy - pos[l+n].fy + image.fy;
		delfz = pos[k].fz - pos[l+n].fz + image.fz;
		r2 = delfx*delfx + delfy*delfy + delfz*delfz;

/*****	Get (1/r dV/dr) ******/
		dedr = wblj(r2,m,n,&eg);

/*****	Resolve L-J forces on BrOct atom center and water atom centers.  ******/
		f[n+1].fx += (delfx *= dedr);
		f[n+1].fy += (delfy *= dedr);
		f[n+1].fz += (delfz *= dedr);
		f[0].fx -= delfx;
		f[0].fy -= delfy;
		f[0].fz -= delfz;
	      }
	      if(sp != 0.0)
	      {
	         for(d = 0; d < 4; d++)
		 {
		   f[d].fx *= s;
		   f[d].fy *= s;
		   f[d].fz *= s;
		 }
		 f[0].fx -= (frc.fx = sp*eg*sdl.fx);
		 f[0].fy -= (frc.fy = sp*eg*sdl.fy);
		 f[0].fz -= (frc.fz = sp*eg*sdl.fz);
		 f[1].fx += frc.fx;
		 f[1].fy += frc.fy;
		 f[1].fz += frc.fz;
		 eg *= s;
	      }
	      k = nw*3+i*BrOs+m;
	      l = j*3;
	      force[k].fx += f[0].fx;
	      force[k].fy += f[0].fy;
	      force[k].fz += f[0].fz;

	      force[l].fx += f[1].fx;
	      force[l].fy += f[1].fy;
	      force[l].fz += f[1].fz;
	      force[++l].fx += f[2].fx;
	      force[l].fy += f[2].fy;
	      force[l].fz += f[2].fz;
	      force[++l].fx += f[3].fx;
	      force[l].fy += f[3].fy;
	      force[l].fz += f[3].fz;
	      H2OBrOV += eg;
	}
/*****	Coulombic interaction: Br-C head group - water  ******/
/*****	Determine image vector for BrO alpha C(i) - water O(j).  ******/
	k = nw*3+i*BrOs;
	l = j*3;
	image.fx = -(sdl.fx = pos[k].fx - pos[l].fx);
	image.fy = -(sdl.fy = pos[k].fy - pos[l].fy);
	image.fz = -(sdl.fz = pos[k].fz - pos[l].fz);
	mvimage(&sdl);
	image.fx += (delfx =sdl.fx);
	image.fy += (delfy =sdl.fy);
	image.fz += (delfz =sdl.fz);
	r2 = delfx*delfx + delfy*delfy + delfz*delfz;
	if(r2 >= swr2max)
		continue;
	if(r2 <= swr2min)
	{
		s = 1.0;
		sp = 0.0;
	}
	else
		swtch(r2-swr2min,&s,&sp);

	for(d = 0; d < 5; d++)
		f[d].fx = f[d].fy = f[d].fz = 0.;
	eg = 0.0;

	for(n = 0; n < 3; n++)/*Loop over atoms in water molecules*/
   	   for(m = 0; m < 2; m++)/*Loop over atoms in BrO Br-C head*/
	   {
		delfx = pos[k+m].fx - pos[l+n].fx + image.fx;
		delfy = pos[k+m].fy - pos[l+n].fy + image.fy;
		delfz = pos[k+m].fz - pos[l+n].fz + image.fz;
		r2 = delfx*delfx + delfy*delfy + delfz*delfz;

/*****	Get (1/r dV/dr) ******/
		dedr = wbq(r2,m,n,&eg);

/*****	Resolve coulombic forces on Br-C headgroup and water atom centers.  ******/
		f[n+2].fx += (delfx *= dedr);
		f[n+2].fy += (delfy *= dedr);
		f[n+2].fz += (delfz *= dedr);
		f[m].fx -= delfx;
		f[m].fy -= delfy;
		f[m].fz -= delfz;
	   }
	   if(sp != 0.0)
	   {
	      for(d = 0; d < 5; d++)
	      {
		f[d].fx *= s;
		f[d].fy *= s;
		f[d].fz *= s;
	      }
	      f[0].fx -= (frc.fx = sp*eg*sdl.fx);
	      f[0].fy -= (frc.fy = sp*eg*sdl.fy);
	      f[0].fz -= (frc.fz = sp*eg*sdl.fz);
	      f[2].fx += frc.fx;
	      f[2].fy += frc.fy;
	      f[2].fz += frc.fz;
	      eg *= s;
	   }
	   k = nw*3+i*BrOs;
	   l = j*3;
	   force[k].fx += f[0].fx;
	   force[k].fy += f[0].fy;
	   force[k].fz += f[0].fz;
	   force[++k].fx += f[1].fx;
	   force[k].fy += f[1].fy;
	   force[k].fz += f[1].fz;
   	   
	   force[l].fx += f[2].fx;
	   force[l].fy += f[2].fy;
	   force[l].fz += f[2].fz;
	   force[++l].fx += f[3].fx;
	   force[l].fy += f[3].fy;
	   force[l].fz += f[3].fz;
	   force[++l].fx += f[4].fx;
	   force[l].fy += f[4].fy;
	   force[l].fz += f[4].fz;
	   H2OBrOV += eg;

    }
}
}

double
wblj(r2,m,n,eg)
double r2;
int m; /* 0<=m<= 8 Br-octane atoms*/
int n; /* 0<=n<= 2 Oxygen and Hydrogen.*/
double *eg;
{
double der,ir,ir6, a,b,aa,bb;
	ir = 1. / r2;
	ir6 = ir * ir * ir;
	a = slj[m+3][n].a;
	b = slj[m+3][n].b;
	aa =12.*a;
	bb = 6.*b;
	*eg +=  ( a * ir6 - b ) * ir6;
	der = (bb - aa * ir6) * ir6 * ir;
	return(der);
}

double
wbq(r2,m,n,eg)
double r2;
int m; /* 0<=m<= 8 Br-octane atoms*/
int n; /* 0<=n<= 2 Oxygen and Hydrogen.*/
double *eg;
{
double der;
double q2,ec,r;
	q2 = slj[m+3][n].q;
	r=sqrt(r2);
	ec = q2/r;
	*eg += ec;
	der = -ec/r2;
	return(der);
}
