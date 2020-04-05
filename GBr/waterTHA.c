#include	<md.h>
#include	<system.h>
#include	 <math.h>

waterTHA()
{
int	i, j, k, l, m;
double	r2, dedr, delfx, delfy, delfz;
double	eg, wnc, s, sp;
tripd	frc, image, sdl, f[4];
double	wtljq();
int nw;

nw = natoms - nBr - nCl2*2 - nTS*3 - nGLY*GLYsites - nTHA*THAsites;

for(i=0;i<nw/3;i++){
   for(j=0;j<THAsites;j++){
      k = i*3;
      l = nw + j; //index of THA atom
      eg = wnc = 0.0;

/*****	Determine image vector, water parent - THA j.  ******/
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
      for(m=0;m<4; m++)
         f[m].fx = f[m].fy = f[m].fz = 0.;
      for(m=0;m<3;m++){ /*Loop over atoms in water molecules*/
/*****	Determine image vector ******/
	 delfx = pos[k+m].fx - pos[l].fx + image.fx;
	 delfy = pos[k+m].fy - pos[l].fy + image.fy;
	 delfz = pos[k+m].fz - pos[l].fz + image.fz;
	 r2 = delfx*delfx + delfy*delfy + delfz*delfz;
/*****	Get (1/r dV/dr) ******/
	 dedr = wtljq(r2,m,j,&eg,&wnc);
/*****	Resolve forces on atoms.  ******/
	 f[3].fx += (delfx *= dedr);//THA atom (l)
	 f[3].fy += (delfy *= dedr);
	 f[3].fz += (delfz *= dedr);
	 f[m].fx -= delfx; //water atom (m or k+m)
	 f[m].fy -= delfy;
	 f[m].fz -= delfz;
      }
      if(sp != 0.0){
	 for(m=0;m<4;m++){
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
	 wnc *= s;
      }
      k=i*3;
      l = nw+j;
	/* water */
	force[k].fx += f[0].fx;
	force[k].fy += f[0].fy;
	force[k].fz += f[0].fz;
	force[++k].fx += f[1].fx;
	force[k].fy += f[1].fy;
	force[k].fz += f[1].fz;
	force[++k].fx += f[2].fx;
	force[k].fy += f[2].fy;
	force[k].fz += f[2].fz;

	/* THA */
	force[l].fx += f[3].fx;
	force[l].fy += f[3].fy;
	force[l].fz += f[3].fz;

	GLYTHAV += eg;
	GLYTHAC += wnc;
   }

}

//fprintf(stderr,"waterTHA_V = %f\n",GLYTHAV);

}

double
wtljq(r2,m,n,eg,wnc)
double r2, *eg,*wnc;
int m, /* m is the water atom index */
    n; /* n is the THA atom index */
{
double der,r,ir,ir6, ec,a,b,q,aa,bb;
int bin;
int rindex;
	ir = 1. / r2;
	ir6 = ir * ir * ir;
	a = WTlj[m][n].a;
	b = WTlj[m][n].b;
	q = WTlj[m][n].q;
	aa =12*a;
	bb = 6*b;
	r = sqrt(r2);
	ec = q/r;
	*wnc += ec;
	*eg +=  ( a * ir6 - b ) * ir6 + ec;
	der = (bb - aa * ir6) * ir6 * ir - ec/r2 ;
	return(der);
}
