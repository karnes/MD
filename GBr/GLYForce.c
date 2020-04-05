#include	<md.h>
#include	<system.h>
#include	 <math.h>
#define EXCON 0.1


GLYForce()
{
int i, j, k, l, m, n, index;
double r2, dedr, delfx, delfy, delfz;
double eg, nc, s, sp;
tripd frc, image, sdl, f[GLYsites*2];
double gtraljq();
int GLYgr();

//fprintf(stderr,"swr2min = %f, swr2max = %f\n",swr2min,swr2max);

gbondE = gbendE = gtorsE = g14E = g15E = 0.;
for(i=0;i<nGLY;i++){
   k = i*GLYsites;
   intraGLY(k);

/*****	Get intermolecular forces for i - j GLY-GLY interactions.  ******/
   for(j=i+1;j<nGLY;j++){
      k = i*14;
      l = j*14;
      eg = nc = 0.0;

/*****	Determine image vector for primeC i - primeC j.  ******/
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
      if(r2 <= swr2min){
         s = 1.0;
	 sp = 0.0;
      }
      else
         swtch(r2-swr2min,&s,&sp);

      for(m=0;m<2*GLYsites; m++)
	 f[m].fx = f[m].fy = f[m].fz = 0.;

/*****	Loop over atoms in GLY molecules ******/
      for(m=0;m<GLYsites;m++)
	 for(n=0;n<GLYsites;n++){
/*****	Determine image vector ******/
	    delfx = pos[k+m].fx - pos[l+n].fx + image.fx;
	    delfy = pos[k+m].fy - pos[l+n].fy + image.fy;
	    delfz = pos[k+m].fz - pos[l+n].fz + image.fz;
	    r2 = delfx*delfx + delfy*delfy + delfz*delfz;
/*****	Get (1/r dV/dr) ******/
	    dedr = gtraljq(r2,m,n,&eg,&nc);
	    if(tc%grdt==0){
	       GLYgr(r2,m,n);
	    }
/*****	Resolve forces on atoms.  ******/
	    f[n+14].fx += (delfx *= dedr);
	    f[n+14].fy += (delfy *= dedr);
	    f[n+14].fz += (delfz *= dedr);
	    f[m].fx -= delfx;
	    f[m].fy -= delfy;
	    f[m].fz -= delfz;
	 }
      if(sp != 0.0){
	 for(m=0;m<2*GLYsites;m++){
	    f[m].fx *= s;
	    f[m].fy *= s;
	    f[m].fz *= s;
	 }
/* make sure that f[0] and f[14] are the center of mass atoms */
	 f[0].fx -= (frc.fx = sp*eg*sdl.fx);
	 f[0].fy -= (frc.fy = sp*eg*sdl.fy);
	 f[0].fz -= (frc.fz = sp*eg*sdl.fz);
	 f[14].fx += frc.fx;
	 f[14].fy += frc.fy;
	 f[14].fz += frc.fz;
	 eg *= s;
	 nc *= s;
      }
      k = i*GLYsites;
      l = j*GLYsites;
      
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
		      
      force[l].fx += f[14].fx;
      force[l].fy += f[14].fy;
      force[l].fz += f[14].fz;
      force[++l].fx += f[15].fx;
      force[l].fy += f[15].fy;
      force[l].fz += f[15].fz;
      force[++l].fx += f[16].fx;
      force[l].fy += f[16].fy;
      force[l].fz += f[16].fz;
      force[++l].fx += f[17].fx;
      force[l].fy += f[17].fy;
      force[l].fz += f[17].fz;
      force[++l].fx += f[18].fx;
      force[l].fy += f[18].fy;
      force[l].fz += f[18].fz;
      force[++l].fx += f[19].fx;
      force[l].fy += f[19].fy;
      force[l].fz += f[19].fz;
      force[++l].fx += f[20].fx;
      force[l].fy += f[20].fy;
      force[l].fz += f[20].fz;
      force[++l].fx += f[21].fx;
      force[l].fy += f[21].fy;
      force[l].fz += f[21].fz;
      force[++l].fx += f[22].fx;
      force[l].fy += f[22].fy;
      force[l].fz += f[22].fz;
      force[++l].fx += f[23].fx;
      force[l].fy += f[23].fy;
      force[l].fz += f[23].fz;
      force[++l].fx += f[24].fx;
      force[l].fy += f[24].fy;
      force[l].fz += f[24].fz;
      force[++l].fx += f[25].fx;
      force[l].fy += f[25].fy;
      force[l].fz += f[25].fz;
      force[++l].fx += f[26].fx;
      force[l].fy += f[26].fy;
      force[l].fz += f[26].fz;
      force[++l].fx += f[27].fx;
      force[l].fy += f[27].fy;
      force[l].fz += f[27].fz;
		      
      INTER_GLY += eg;
      GLYC += nc;
      }
   }

}

double
gtraljq(r2,m,n,eg,nc)
double r2, *eg, *nc;
int m,n;
{
double der,r,ir,ir6, ec,a,b,q,aa,bb,sqrt();

	ir = 1. / r2;
	ir6 = ir * ir * ir;
	a = GLYlj[m][n].a;
	b = GLYlj[m][n].b;

	q = GLYlj[m][n].q;
	aa =12*a;
	bb = 6*b;

	r = sqrt(r2);
	ec = q/r;

	*eg +=  ( a * ir6 - b ) * ir6 + ec;
	*nc += ec;
	der = (bb - aa * ir6) * ir6 * ir - ec/r2;
	return(der);
}

intraGLY(k)
int k;
{
int i,n, n1, n2, n3, m;
double dx, dy, dz, r2;
double calcstrch(), calcbend(), calctorq(),gtraljq14(), gtraljq15();  
double gbond_e, gtors_e, gbend_e, g15_e, g14_e;
double e14,c14,e15,c15;
double dedr;
/* loop over all 1-5+ pairs */
g15_e = 0.;
for(i=0;i<nGLY15;i++){
   e15 = c15 = 0.;
   dx=pos[k+GLY15[i][0]].fx-pos[k+GLY15[i][1]].fx;
   dy=pos[k+GLY15[i][0]].fy-pos[k+GLY15[i][1]].fy;
   dz=pos[k+GLY15[i][0]].fz-pos[k+GLY15[i][1]].fz;
   r2=dx*dx + dy*dy + dz*dz;
   dedr = gtraljq15(r2,GLY15[i][0],GLY15[i][1],&e15,&c15);
// ****	Resolve forces on atoms.  ******
   force[k+GLY15[i][1]].fx += (dx *= dedr);
   force[k+GLY15[i][1]].fy += (dy *= dedr);
   force[k+GLY15[i][1]].fz += (dz *= dedr);
   force[k+GLY15[i][0]].fx -= dx;
   force[k+GLY15[i][0]].fy -= dy;
   force[k+GLY15[i][0]].fz -= dz;
//   if(tc==0)
//      fprintf(stderr,"1-5 i=%02d (%d,%d) U = %f\n",i,GLY15[i][0],GLY15[i][1],e15*KCAL);
   g15_e += e15;
}
/* loop over all 1-4 pairs */
g14_e = 0.;
for(i=0;i<nGLYtors;i++){
   dx=pos[k+GLYtors[i][0]].fx-pos[k+GLYtors[i][3]].fx;
   dy=pos[k+GLYtors[i][0]].fy-pos[k+GLYtors[i][3]].fy;
   dz=pos[k+GLYtors[i][0]].fz-pos[k+GLYtors[i][3]].fz;
   r2=dx*dx + dy*dy + dz*dz;
   e14 = c14 = 0.;
   dedr = gtraljq14(r2,GLYtors[i][0],GLYtors[i][3],&e14,&c14);
// ****	Resolve forces on atoms.  ******
   force[k+GLYtors[i][3]].fx += (dx *= dedr);
   force[k+GLYtors[i][3]].fy += (dy *= dedr);
   force[k+GLYtors[i][3]].fz += (dz *= dedr);
   force[k+GLYtors[i][0]].fx -= dx;
   force[k+GLYtors[i][0]].fy -= dy;
   force[k+GLYtors[i][0]].fz -= dz;
//   if(tc==0)
//      fprintf(stderr,"1-4 i=%02d (%d,%d) U = %f\n",i,GLYtors[i][0],GLYtors[i][3],e14*KCAL);
   g14_e += e14;
}
/* loop over all bond stretches*/
gbond_e = 0.;
for(i=0;i<nGLYstr;i++){
   n = k+GLYstr[i][0];
   m = k+GLYstr[i][1];
   gbond_e += calcstrch(n,m,GLYstr[i][2]);
}
/* loop over all torsions*/
gtors_e=0.;
/*i-j-k-l*/
for(i=0;i<nGLYtors;i++){
   n = k + GLYtors[i][0];/*i*/
   n1= k + GLYtors[i][1];/*j*/
   n2= k + GLYtors[i][2];/*k*/
   n3= k + GLYtors[i][3];/*l*/
   gtors_e += calctorq(n,n1,n2,n3,GLYtors[i][4]);
}
/* loop over all angle bends*/
gbend_e = 0.;
/*i-j-k*/
for(i=0;i<nGLYbend;i++){
   n = k + GLYbend[i][0];/*i*/
   n1= k + GLYbend[i][1];/*j*/
   n2= k + GLYbend[i][2];/*k*/
   gbend_e += calcbend(n,n1,n2,GLYbend[i][3]);
}
/* total sum for all molecules*/
gbondE += gbond_e;
gbendE += gbend_e;
gtorsE += gtors_e;
g14E += g14_e;
g15E += g15_e;
/* total intramolecular potential energy*/
INTRA_GLY += gbond_e + gtors_e + gbend_e + g14_e + g15_e;
}

double calcstrch(i,j,k)
int i,j,k;
{
tripd 	d; 
double	bond, dedr, pe;

d.fx = pos[i].fx - pos[j].fx;
d.fy = pos[i].fy - pos[j].fy;
d.fz = pos[i].fz - pos[j].fz;
bond = sqrt(d.fx*d.fx+d.fy*d.fy+d.fz*d.fz);
dedr = gkstr[k]*(bond-greq[k])/bond;
pe = 0.5*gkstr[k]*sq(bond-greq[k]);
force[i].fx -= (d.fx = dedr*d.fx);
force[i].fy -= (d.fy = dedr*d.fy);
force[i].fz -= (d.fz = dedr*d.fz);
force[j].fx += d.fx;
force[j].fy += d.fy;
force[j].fz += d.fz;
return(pe);

}
double calctorq (i1,i2,i3,i4,m)
int i1,i2,i3,i4,m;
{
int i,j,k,index,ats[4];
double dmat[3][3], 
      cmat[3][3],
      grad[4][3], 
      d[3][3], 
      da,cosa, 
      K1,K10,K11,K12,K13,K21,K22,K23, 
      D1201, C1100, 
      dudcos, 
      pe,V1,V2;
ats[0] = i1;
ats[1] = i2;
ats[2] = i3;
ats[3] = i4;
    for (i=0;  i<3;  i++){/* get vectors for atoms */
        d[i][0]=pos[ats[i+1]].fx-pos[ats[i]].fx;
        d[i][1]=pos[ats[i+1]].fy-pos[ats[i]].fy;
        d[i][2]=pos[ats[i+1]].fz-pos[ats[i]].fz;
    }
    for (i=0; i<3;  i++){
        for (j=0;  j<3;  j++){
	    cmat[i][j]=0.0;
            for (k=0;  k<3;  k++){
                cmat[i][j] += d[i][k]*d[j][k];
            }
        }
    }
    for (i=0; i<3;  i++){
        for (j=0;  j<3;  j++){
            dmat[i][j] = cmat[i][i]*cmat[j][j]-cmat[i][j]*cmat[i][j];
        }
    }

    
    K1 = cmat[0][1]*cmat[1][2] - cmat[0][2]*cmat[1][1];
    D1201 = sqrt (dmat[1][2]*dmat[0][1]);

    if (D1201 == 0){
	fprintf(stderr,"No torsion defined\n");
	exit(1);
    }


    /*
    **	    minus gradient at atom1 with respect to the cosine of the angle
    */

    for (i=0;  i<3;  i++){
        grad[0][i]=( ( (K1* (cmat[0][1] * d[1][i]-cmat[1][1] * d[0][i]) )
	/dmat[0][1]) + cmat[1][2]*d[1][i] - cmat[1][1]*d[2][i])/D1201;
    }

    /*
    **	    minus gradient at atom2 with respect to the cosine of the angle
    */

    K11 = cmat[0][1] + cmat[1][1];
    K12 = cmat[1][2] + 2*cmat[0][2];
    K13 = cmat[0][0] + cmat[0][1];
   
    for (i=0;  i<3;  i++){
	V1 = cmat[1][2]*d[2][i] - cmat[2][2]*d[1][i];
	V2 = K11*d[0][i] - K13*d[1][i];
        grad[1][i] = (K11*d[2][i] - K12*d[1][i] + cmat[1][2]*d[0][i] + 
			K1*(V1/dmat[1][2] + V2/dmat[0][1]))/D1201;
    }

    /*
    **	   minus gradient at atom3 with respect to the cosine of the angle
    */
    K21= cmat[0][1] + 2*cmat[0][2];
    K22= cmat[1][2] + cmat[1][1];
    K23= cmat[2][2] + cmat[1][2];
    for (i=0;  i<3;  i++){
        V1 = cmat[0][0]*d[1][i] - cmat[0][1]*d[0][i];
	V2 = K23*d[1][i] - K22*d[2][i];
	grad[2][i] = (K21*d[1][i] - K22*d[0][i] - cmat[0][1]*d[2][i] +
			K1*(V1/dmat[0][1] + V2/dmat[1][2]))/D1201;
    }

    /*
    **	    minus gradient at atom4 with respect to the cosine of the angle
    */
    for (i=0;  i<3;  i++){
        grad[3][i] = (-cmat[0][1]*d[1][i] + cmat[1][1]*d[0][i] + 
	    K1*(cmat[1][1]*d[2][i] - cmat[1][2]*d[1][i])/dmat[1][2])/D1201;
    }
    /*
    **	    calculate cos angle
    */
    cosa = K1/D1201;
    /*IMPA[m] = cosa;*/
/*
    index = acos(cosa)*36/PI;
    gTors[l][index]++;
*/
    /*
    **	    evaluate dU/dcos(a)
    **  U = pTors[l][0]*(1+cos(a))+pTors[l][1]*(1-cos(2*a))+pTors[l][2]*(1+cos(3*a))
    */
//  OPLS-style
//    pe = (1+cosa)*(gtors[m][0]+2*gtors[m][1]*(1-cosa)+gtors[m][2]*sq(2*cosa-1));
//    dudcos = gtors[m][0] - 4*gtors[m][1]*cosa + 3*gtors[m][2]*(4*cosa*cosa-1);
//  AMBER style
    pe = gtors[m][1]*2.0*cosa*cosa + (1+cosa)*(gtors[m][0]+gtors[m][2]*sq(2*cosa-1));
    dudcos = gtors[m][0] + 4*gtors[m][1]*cosa + 3*gtors[m][2]*(4*cosa*cosa-1);
    /*
    **	    add resulting force to atoms for angle
    */
    for (i=0;  i<4;  i++){
        force[ats[i]].fx += dudcos*grad[i][0];
	force[ats[i]].fy += dudcos*grad[i][1];
	force[ats[i]].fz += dudcos*grad[i][2];
    }
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
da = acos(n1n2) - geqbend[l];
pe = 0.5*gkbend[l]*da*da;
dedd = gkbend[l]*da/ (sqrt(1.-n1n2*n1n2)*r1*r2);
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
double
gtraljq15(r2,m,n,eg,nc)
double r2, *eg, *nc;
int m,n;
{
double der,r,ir,ir6, ec,a,b,q,aa,bb,sqrt();

	ir = 1. / r2;
	ir6 = ir * ir * ir;
	a = GLYlj[m][n].a;
	b = GLYlj[m][n].b;

	q = GLYlj[m][n].q;
	aa =12*a;
	bb = 6*b;

	r = sqrt(r2);
	ec = q/r;

	*eg +=  ( a * ir6 - b ) * ir6 + ec;
	*nc += ec;
	der = (bb - aa * ir6) * ir6 * ir - ec/r2;
	return(der);
}


double
gtraljq14(r2,m,n,eg,nc)
double r2, *eg, *nc;
int m,n;
{
double der,r,ir,ir6, ec,a,b,q,aa,bb,sqrt();

	ir = 1. / r2;
	ir6 = ir * ir * ir;
	a = GLYlj[m][n].a;
	b = GLYlj[m][n].b;

	q = GLYlj[m][n].q;
	aa =12*a;
	bb = 6*b;

	r = sqrt(r2);
	ec = q/r;

	*eg += (g14ljscale)*(( a * ir6 - b ) * ir6) + (g14qscale)*ec;
	*nc += (g14qscale)*ec;
	der  = (g14ljscale)*((bb - aa * ir6) * ir6 * ir) - (g14qscale)*ec/r2;
	return(der);
}

GLYgr(r2,m,n)
double r2;
int m,n;
{
 int index;

 index=(int)((sqrt(r2))/binRDF);
 if(index<300){
   /*****  rdf calculations *****/
   // m = CARBON
   //0: C-C, 2: C-O, 3: C-H
   if(m==0 || m==4 || m==9){ //m = C
      if(n==0 || n==4 || n==9){ // C
	 grGLY[0][index]+=(2.0/9.0);
      }
      else if(n==1 || n==5 || n==10){ // O
	 grGLY[2][index]+=(1.0/9.0);
      }
      else{ // n is H
	 grGLY[3][index]+=(1.0/24.0);
      }
   }
   // m = OXYGEN
   //1: O-H(hydrox), 4: O-H (any), 6: O-O
   else if(m==1 || m==5 || m==10){//m = O
      if(n==2 || n==6 || n==11){//n = H(hyd)
         grGLY[1][index]+=(1.0/9.0);
         grGLY[4][index]+=(1.0/24.0);
      }
      else if(n==1 || n==5 || n==10){// n= O
         grGLY[6][index]+=(2.0/9.0);
      }
      else if(n==3 || n==7 || n==8 || n==12 || n==13){//n = H
         grGLY[4][index]+=(1.0/24.0);
      }
      else{ //n = C
	 grGLY[2][index]+=(1.0/9.0);
      }
   }
   // m = HYDROGEN
   //5: H-H, 3: H-C, 4: H(any)-O
   else{
      if(n==2 || n==3 || n==6 || n==7 || n==8 || n==11 || n==12 || n==13){ //n = H
         grGLY[5][index]+=(1.0/32.0);
      }
      else if(n==1 || n==5 || n==10){ // n = O
         grGLY[4][index]+=(1.0/24.0);
      }
      else{ // n = C
         grGLY[3][index]+=(1.0/24.0);
      }
   }
   if(m==2 || m==6 || m==11){ // Hydroxyl H
      if(n==1 || n==5 || n==10){ //O
         grGLY[1][index]+=(1.0/9.0);
      }
   }
 }
}
