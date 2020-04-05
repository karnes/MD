#include	<md.h>
#include	<system.h>
#include	 <math.h>

//GLYTHA(pos,force)
//	tripd	*pos;
//	tripd	*force;
GLYTHA()
{
int	i, j, k, l, m;
double	r2, dedr, delfx, delfy, delfz;
double	eg, wnc, s, sp;
tripd	frc, image, sdl, f[GLYsites+1]; //f[GLYsites+THAsites];
double	gtljq();

//return(0);

for(i=0;i<nGLY;i++){
   for(j=0;j<THAsites;j++){
      k = i*GLYsites;
      l = nGLY*GLYsites + j; //index of THA atom
      eg = wnc = 0.0;

/*****	Determine image vector, GLY parent - THA j.  ******/
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
      for(m=0;m<(GLYsites+1); m++)
         f[m].fx = f[m].fy = f[m].fz = 0.;
      for(m=0;m<GLYsites;m++){ /*Loop over atoms in GLY molecules*/
/*****	Determine image vector ******/
	 delfx = pos[k+m].fx - pos[l].fx + image.fx;
	 delfy = pos[k+m].fy - pos[l].fy + image.fy;
	 delfz = pos[k+m].fz - pos[l].fz + image.fz;
	 r2 = delfx*delfx + delfy*delfy + delfz*delfz;
/*****	Get (1/r dV/dr) ******/
	 dedr = gtljq(r2,m,j,&eg,k,&wnc,l);
/*****	Resolve forces on atoms.  ******/
	 f[GLYsites].fx += (delfx *= dedr);//THA atom (l)
	 f[GLYsites].fy += (delfy *= dedr);
	 f[GLYsites].fz += (delfz *= dedr);
	 f[m].fx -= delfx; //GLY atom (m or k+m)
	 f[m].fy -= delfy;
	 f[m].fz -= delfz;
      }
      if(sp != 0.0){
	 for(m=0;m<GLYsites+1;m++){
	       f[m].fx *= s;
	       f[m].fy *= s;
	       f[m].fz *= s;
	 }
	 f[0].fx -= (frc.fx = sp*eg*sdl.fx);
	 f[0].fy -= (frc.fy = sp*eg*sdl.fy);
	 f[0].fz -= (frc.fz = sp*eg*sdl.fz);
	 f[GLYsites].fx += frc.fx;
	 f[GLYsites].fy += frc.fy;
	 f[GLYsites].fz += frc.fz;
	 eg *= s;
	 wnc *= s;
      }
      k=i*GLYsites;
      l = nGLY*GLYsites+j;
	/* GLY */
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

	/* THA */
	force[l].fx += f[GLYsites].fx;
	force[l].fy += f[GLYsites].fy;
	force[l].fz += f[GLYsites].fz;
	/*   
	force[l].fx += f[14].fx;
	force[l].fy += f[14].fy;
	force[l].fz += f[14].fz;
	force[++l].fx += f[15].fx;
	force[l].fy += f[15].fy;
	force[l].fz += f[15].fz;
	force[++l].fx += f[16].fx;
	force[l].fy += f[16].fy;
	force[l].fz += f[16].fz;
	*/
	GLYTHAV += eg;
	GLYTHAC += wnc;
   }

}

}

double
gtljq(r2,m,n,eg,nN,wnc,l)
double r2, *eg,*wnc;
int m, /* m is the GLY atom index */
    n, /* n is the THA atom index */
    nN,/* the number of the nitrobenzene molecule for the rdf calculations*/
    l; /* index of water O*/
{
double der,r,ir,ir6, ec,a,b,q,aa,bb;
int bin;
int rindex;
//	nn = n + 14;   /* water - NIT cross terms begin at 14 */
	ir = 1. / r2;
	ir6 = ir * ir * ir;
	a = GTlj[m][n].a;
	b = GTlj[m][n].b;
	q = GTlj[m][n].q;
	aa =12*a;
	bb = 6*b;
	r = sqrt(r2);
/*	if(l==0){
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
*/
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
