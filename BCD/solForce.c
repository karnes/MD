#include	<md.h>
#include	<system.h>
#include	<math.h>
#include	<stdio.h>

/*
 *	this routine loops over the solvent and the solute  atoms and
 *	calculates the non-bond forces and the energies of interaction.
 */

solForce(pos, force)
tripd	*pos;
tripd	*force;
{
int k,i,j,m,n,nw;
double pe,vc;
tripd fw[3],fs[2],fb[BrOs],fcd[1];

//fprintf(stderr,"solForce.c -- entered.\n");
if(nsolute != 3)
	return;
for(i=0;i<natoms;i++){
   fevb[i].fx = fevb[i].fy = fevb[i].fz = 0.0;
}
nw = natoms - BrOs*nBrO - nBCD*BCDs - nsolute;
k = natoms - nsolute;
VINT_WI1 = VINT_WI2 = VINT_BI1 = VINT_BI2 = 0.0;
VINT_WD1 = VINT_WD2 = VINT_BD1 = VINT_BD2 = 0.0;
VINT_CDI1 = VINT_CDI2 = VINT_CDD1 = VINT_CDD2 = 0.0;
VCL_WI1 = VCL_WI2 = VCL_BI1 = VCL_BI2 = 0.0;
VCL_WD1 = VCL_WD2 = VCL_BD1 = VCL_BD2 = 0.0;
VCL_CDI1 = VCL_CDI2 = VCL_CDD1 = VCL_CDD2 = 0.0;
VINT_oCDI1 = VINT_oCDI2 = VINT_oCDD1 = VINT_oCDD2 = 0.0;
VCL_oCDI1 = VCL_oCDI2 = VCL_oCDD1 = VCL_oCDD2 = 0.0;
V_teth = 0.0;
for(i = 0; i < nw; i = i+3){        /* loop over H2O molecules*/
	/* H2O  - ion1 interactions - forces only half*/
	pe = 0.;
	vc = 0.;
	n = 1;
	wion(&pos[i],&pos[k+1],&fw,&fs,&pe,&vc,n);
	VINT_WI1 += pe;
	VCL_WI1 += vc;
/*update global force array*/
	force[i].fx += fw[0].fx;force[i].fy += fw[0].fy;force[i].fz += fw[0].fz;
	force[i+1].fx += fw[1].fx;force[i+1].fy += fw[1].fy;force[i+1].fz += fw[1].fz;
	force[i+2].fx += fw[2].fx;force[i+2].fy += fw[2].fy;force[i+2].fz += fw[2].fz;
	force[k+1].fx += fs[0].fx;force[k+1].fy += fs[0].fy;force[k+1].fz += fs[0].fz;
/*update global EVB force array*/
	fevb[i].fx += fw[0].fx;fevb[i].fy += fw[0].fy;fevb[i].fz += fw[0].fz;
	fevb[i+1].fx += fw[1].fx;fevb[i+1].fy += fw[1].fy;fevb[i+1].fz += fw[1].fz;
	fevb[i+2].fx += fw[2].fx;fevb[i+2].fy += fw[2].fy;fevb[i+2].fz += fw[2].fz;
	fevb[k+1].fx += fs[0].fx;fevb[k+1].fy += fs[0].fy;fevb[k+1].fz += fs[0].fz;
	/* H2O  - ion2 interactions - forces only half*/
	pe = 0.;
	vc = 0.;
	n = 2;
	wion(&pos[i],&pos[k+2],&fw,&fs,&pe,&vc,n);
	VINT_WI2 += pe;
	VCL_WI2 += vc;
/*update global force array*/
	force[i].fx += fw[0].fx;force[i].fy += fw[0].fy;force[i].fz += fw[0].fz;
	force[i+1].fx += fw[1].fx;force[i+1].fy += fw[1].fy;force[i+1].fz += fw[1].fz;
	force[i+2].fx += fw[2].fx;force[i+2].fy += fw[2].fy;force[i+2].fz += fw[2].fz;
	force[k+2].fx += fs[0].fx;force[k+2].fy += fs[0].fy;force[k+2].fz += fs[0].fz;
/*update global EVB force array*/
	fevb[i].fx -= fw[0].fx;fevb[i].fy -= fw[0].fy;fevb[i].fz -= fw[0].fz;
	fevb[i+1].fx -= fw[1].fx;fevb[i+1].fy -= fw[1].fy;fevb[i+1].fz -= fw[1].fz;
	fevb[i+2].fx -= fw[2].fx;fevb[i+2].fy -= fw[2].fy;fevb[i+2].fz -= fw[2].fz;
	fevb[k+2].fx -= fs[0].fx;fevb[k+2].fy -= fs[0].fy;fevb[k+2].fz -= fs[0].fz;
	/* H2O  - dipol1 interactions - forces only half*/
	pe = 0.;
	vc = 0.;
	wdip(&pos[i],&pos[k+2],&fw,&fs,&pe,&vc);
//fprintf(stderr,"solForce.c: after wdip 1\n");
	VINT_WD1 += pe;
	VCL_WD1 += vc;
/*update global force array*/
	force[i].fx += fw[0].fx;force[i].fy += fw[0].fy;force[i].fz += fw[0].fz;
	force[i+1].fx += fw[1].fx;force[i+1].fy += fw[1].fy;force[i+1].fz += fw[1].fz;
	force[i+2].fx += fw[2].fx;force[i+2].fy += fw[2].fy;force[i+2].fz += fw[2].fz;
	force[k].fx += fs[0].fx;force[k].fy += fs[0].fy;force[k].fz += fs[0].fz;
	force[k+2].fx += fs[1].fx;force[k+2].fy += fs[1].fy;force[k+2].fz += fs[1].fz;
/*update global EVB force array*/
	fevb[i].fx += fw[0].fx;fevb[i].fy += fw[0].fy;fevb[i].fz += fw[0].fz;
	fevb[i+1].fx += fw[1].fx;fevb[i+1].fy += fw[1].fy;fevb[i+1].fz += fw[1].fz;
	fevb[i+2].fx += fw[2].fx;fevb[i+2].fy += fw[2].fy;fevb[i+2].fz += fw[2].fz;
	fevb[k].fx += fs[0].fx;fevb[k].fy += fs[0].fy;fevb[k].fz += fs[0].fz;
	fevb[k+2].fx += fs[1].fx;fevb[k+2].fy += fs[1].fy;fevb[k+2].fz += fs[1].fz;
	/* H2O  - dipol2 interactions - forces only half*/
	pe = 0.;
	vc = 0.;
	wdip(&pos[i],&pos[k+1],&fw,&fs,&pe,&vc);
//fprintf(stderr,"solForce.c: after wdip 2\n");
	VINT_WD2 += pe;
	VCL_WD2 += vc;
/*update global force array*/
	force[i].fx += fw[0].fx;force[i].fy += fw[0].fy;force[i].fz += fw[0].fz;
	force[i+1].fx += fw[1].fx;force[i+1].fy += fw[1].fy;force[i+1].fz += fw[1].fz;
	force[i+2].fx += fw[2].fx;force[i+2].fy += fw[2].fy;force[i+2].fz += fw[2].fz;
	force[k].fx += fs[0].fx;force[k].fy += fs[0].fy;force[k].fz += fs[0].fz;
	force[k+1].fx += fs[1].fx;force[k+1].fy += fs[1].fy;force[k+1].fz += fs[1].fz;
/*update global EVB force array*/
	fevb[i].fx -= fw[0].fx;fevb[i].fy -= fw[0].fy;fevb[i].fz -= fw[0].fz;
	fevb[i+1].fx -= fw[1].fx;fevb[i+1].fy -= fw[1].fy;fevb[i+1].fz -= fw[1].fz;
	fevb[i+2].fx -= fw[2].fx;fevb[i+2].fy -= fw[2].fy;fevb[i+2].fz -= fw[2].fz;
	fevb[k].fx -= fs[0].fx;fevb[k].fy -= fs[0].fy;fevb[k].fz -= fs[0].fz;
	fevb[k+1].fx -= fs[1].fx;fevb[k+1].fy -= fs[1].fy;fevb[k+1].fz -= fs[1].fz;
}
//fprintf(stderr,"after water. k = %d\n",k);
VINT_EVB += (VINT_WI1+VINT_WD1+VINT_WI2+VINT_WD2)/2.0;
//fprintf(stderr,"solForce.c: water part complete.\n");
/*Now for Br-octane*/
for(i = nw; i < nw + nBrO*BrOs; i = i+BrOs){    /* loop over Br-octane molecules*/
	/* Br-oct  - ion1 interactions - forces only half*/
	pe = 0.;
	vc = 0.;
	n = 1;
	bion(&pos[i],&pos[k+1],&fb,&fs,&pe,&vc,n);
	VINT_BI1 += pe;
	VCL_BI1 += vc;
/*update global force array*/
	for (j=0;j<BrOs;j++){
	 force[i+j].fx += fb[j].fx;force[i+j].fy += fb[j].fy;force[i+j].fz += fb[j].fz;
	}
	force[k+1].fx += fs[0].fx;force[k+1].fy += fs[0].fy;force[k+1].fz += fs[0].fz;
/*update global EVB force array*/
	for (j=0;j<BrOs;j++){
	 fevb[i+j].fx += fb[j].fx;fevb[i+j].fy += fb[j].fy;fevb[i+j].fz += fb[j].fz;
	}
	fevb[k+1].fx += fs[0].fx;fevb[k+1].fy += fs[0].fy;fevb[k+1].fz += fs[0].fz;
	/* Br-oct  - ion2 interactions - forces only half*/
	pe = 0.;
	vc = 0.;
	n = 2;
	bion(&pos[i],&pos[k+2],&fb,&fs,&pe,&vc,n);
	VINT_BI2 += pe;
	VCL_BI2 += vc;
/*update global force array*/
	for (j=0;j<BrOs;j++){
	 force[i+j].fx += fb[j].fx;force[i+j].fy += fb[j].fy;force[i+j].fz += fb[j].fz;
	}
	force[k+2].fx += fs[0].fx;force[k+2].fy += fs[0].fy;force[k+2].fz += fs[0].fz;
/*update global EVB force array*/
	for (j=0;j<BrOs;j++){
	 fevb[i+j].fx -= fb[j].fx;fevb[i+j].fy -= fb[j].fy;fevb[i+j].fz -= fb[j].fz;
	}
	fevb[k+2].fx -= fs[0].fx;fevb[k+2].fy -= fs[0].fy;fevb[k+2].fz -= fs[0].fz;
	/* Br-oct  - dipol1 interactions - forces only half*/
	pe = 0.;
	vc = 0.;
	bdip(&pos[i],&pos[k+2],&fb,&fs,&pe,&vc);
	VINT_BD1 += pe;
	VCL_BD1 += vc;
/*update global force array*/
	for (j=0;j<BrOs;j++){
	 force[i+j].fx += fb[j].fx;force[i+j].fy += fb[j].fy;force[i+j].fz += fb[j].fz;
	}
	force[k].fx += fs[0].fx;force[k].fy += fs[0].fy;force[k].fz += fs[0].fz;
	force[k+2].fx += fs[1].fx;force[k+2].fy += fs[1].fy;force[k+2].fz += fs[1].fz;
/*update global EVB force array*/
	for (j=0;j<BrOs;j++){
	 fevb[i+j].fx += fb[j].fx;fevb[i+j].fy += fb[j].fy;fevb[i+j].fz += fb[j].fz;
	}
	fevb[k].fx += fs[0].fx;fevb[k].fy += fs[0].fy;fevb[k].fz += fs[0].fz;
	fevb[k+2].fx += fs[1].fx;fevb[k+2].fy += fs[1].fy;fevb[k+2].fz += fs[1].fz;
	/* Br-oct  - dipol2 interactions - forces only half*/
	pe = 0.;
	vc = 0.;
	bdip(&pos[i],&pos[k+1],&fb,&fs,&pe,&vc);
	VINT_BD2 += pe;
	VCL_BD2 += vc;
/*update global force array*/
	for (j=0;j<BrOs;j++){
	 force[i+j].fx += fb[j].fx;force[i+j].fy += fb[j].fy;force[i+j].fz += fb[j].fz;
	}
	force[k].fx += fs[0].fx;force[k].fy += fs[0].fy;force[k].fz += fs[0].fz;
	force[k+1].fx += fs[1].fx;force[k+1].fy += fs[1].fy;force[k+1].fz += fs[1].fz;
/*update global EVB force array*/
	for (j=0;j<BrOs;j++){
	 fevb[i+j].fx -= fb[j].fx;fevb[i+j].fy -= fb[j].fy;fevb[i+j].fz -= fb[j].fz;
	}
	fevb[k].fx -= fs[0].fx;fevb[k].fy -= fs[0].fy;fevb[k].fz -= fs[0].fz;
	fevb[k+1].fx -= fs[1].fx;fevb[k+1].fy -= fs[1].fy;fevb[k+1].fz -= fs[1].fz;
}
VINT_EVB += ( VINT_BI1 + VINT_BD1 + VINT_BI2 + VINT_BD2 ) / 2.0;
//fprintf(stderr,"solForce.c: Br-octane part complete.\n");
/* Finally b-CD - solute forces */
if(nBCD==1){
  // fprintf(stderr,"in BCD. k = %d\n",k);
   for(j=0;j<BCDs;j++){
	m = j%21; // site in glucose unit
	i = nw + nBrO*BrOs + j;
//fprintf(stderr,"solForec... BCD atom i = %d\n",i);
//fprintf(stderr,"in BCD. i = %d\n",i);
	/* bCD - ion1 interactions - forces only half*/
	pe = 0.;
	vc = 0.;
	n = 1;
	cdion(&pos[i],&pos[k+1],&fcd,&fs,&pe,&vc,n,m);
	VINT_CDI1 += pe;
	VCL_CDI1 += vc;
	if(m>13 && m!=16){
	   VINT_oCDI1 += pe;
	   VCL_oCDI1 += vc;
	}
/*	if(tc==0){
	   fprintf(stderr,"m = %d, n = %d, pe = %8.3f, vc = %8.3f\n",m,n,pe*KCAL,vc*KCAL);
	}*/
/*update global force array*/
	force[i].fx += fcd[0].fx;force[i].fy += fcd[0].fy;force[i].fz += fcd[0].fz;
	force[k+1].fx += fs[0].fx;force[k+1].fy += fs[0].fy;force[k+1].fz += fs[0].fz;
/*update global EVB force array*/
	fevb[i].fx += fcd[0].fx;fevb[i].fy += fcd[0].fy;fevb[i].fz += fcd[0].fz;
	fevb[k+1].fx += fs[0].fx;fevb[k+1].fy += fs[0].fy;fevb[k+1].fz += fs[0].fz;
	/* b-CD  - ion2 interactions - forces only half*/
	pe = 0.;
	vc = 0.;
	n = 2;
	cdion(&pos[i],&pos[k+2],&fcd,&fs,&pe,&vc,n,m);
	VINT_CDI2 += pe;
	VCL_CDI2 += vc;
	if(m>13 && m!=16){
	   VINT_oCDI2 += pe;
	   VCL_oCDI2 += vc;
	}
/*	if(tc==0){
	   fprintf(stderr,"m = %d, n = %d, pe = %8.3f, vc = %8.3f\n",m,n,pe*KCAL,vc*KCAL);
	}*/
/*update global force array*/
	force[i].fx += fcd[0].fx;force[i].fy += fcd[0].fy;force[i].fz += fcd[0].fz;
	force[k+2].fx += fs[0].fx;force[k+2].fy += fs[0].fy;force[k+2].fz += fs[0].fz;
/*update global EVB force array*/
	fevb[i].fx -= fcd[0].fx;fevb[i].fy -= fcd[0].fy;fevb[i].fz -= fcd[0].fz;
	fevb[k+2].fx -= fs[0].fx;fevb[k+2].fy -= fs[0].fy;fevb[k+2].fz -= fs[0].fz;
	/* b-CD  - dipol1 interactions - forces only half*/
	pe = 0.;
	vc = 0.;
	cddip(&pos[i],&pos[k+2],&fcd,&fs,&pe,&vc,m);
	VINT_CDD1 += pe;
	VCL_CDD1 += vc;
	if(m>13 && m!=16){
	   VINT_oCDD1 += pe;
	   VCL_oCDD1 += vc;
	}
/*update global force array*/
	force[i].fx += fcd[0].fx;force[i].fy += fcd[0].fy;force[i].fz += fcd[0].fz;
	force[k].fx += fs[0].fx;force[k].fy += fs[0].fy;force[k].fz += fs[0].fz;
	force[k+2].fx += fs[1].fx;force[k+2].fy += fs[1].fy;force[k+2].fz += fs[1].fz;
/*update global EVB force array*/
	fevb[i].fx += fcd[0].fx;fevb[i].fy += fcd[0].fy;fevb[i].fz += fcd[0].fz;
	fevb[k].fx += fs[0].fx;fevb[k].fy += fs[0].fy;fevb[k].fz += fs[0].fz;
	fevb[k+2].fx += fs[1].fx;fevb[k+2].fy += fs[1].fy;fevb[k+2].fz += fs[1].fz;
	/* b-CD  - dipol2 interactions - forces only half*/
	pe = 0.;
	vc = 0.;
	cddip(&pos[i],&pos[k+1],&fcd,&fs,&pe,&vc,m);
	VINT_CDD2 += pe;
	VCL_CDD2 += vc;
	if(m>13 && m!=16){
	   VINT_oCDD2 += pe;
	   VCL_oCDD2 += vc;
	}
/*update global force array*/
	force[i].fx += fcd[0].fx;force[i].fy += fcd[0].fy;force[i].fz += fcd[0].fz;
	force[k].fx += fs[0].fx;force[k].fy += fs[0].fy;force[k].fz += fs[0].fz;
	force[k+1].fx += fs[1].fx;force[k+1].fy += fs[1].fy;force[k+1].fz += fs[1].fz;
/*update global EVB force array*/
	fevb[i].fx -= fcd[0].fx;fevb[i].fy -= fcd[0].fy;fevb[i].fz -= fcd[0].fz;
	fevb[k].fx -= fs[0].fx;fevb[k].fy -= fs[0].fy;fevb[k].fz -= fs[0].fz;
	fevb[k+1].fx -= fs[1].fx;fevb[k+1].fy -= fs[1].fy;fevb[k+1].fz -= fs[1].fz;
   }
}
//fprintf(stderr,"after BCD. k = %d\n",k);
VINT_EVB += ( VINT_CDI1 + VINT_CDD1 + VINT_CDI2 + VINT_CDD2 ) / 2.0;
VSN2_BCD += ( VINT_CDI1 + VINT_CDD1 + VINT_CDI2 + VINT_CDD2 ) / 2.0;
VSN2_BCDoh += ( VINT_oCDI1 + VINT_oCDI2 + VINT_oCDD1 + VINT_oCDD2 ) / 2.0;
//fprintf(stderr,"solForce.c: b-CD part complete.\n");

if(nBCD==1 && nsolute>=3){
   BCDsolTether();
   VINT += V_teth + V_solB;
}
//fprintf(stderr,"Ewidth_w = %f, zwall = %f\n",Ewidth_w,zwall);
if(Ewidth_w > zwall) return; /* no window */
EVBwindow();
}

wion(ri,rj,fi,fj,pe,vc,n)
tripd *ri;  /* ri[0-2], are position of the H2O atoms*/
tripd *rj;  /* rj[0] is the position of the ion */
tripd *fi;  /* fi[0-2], are the forces on H2O's atoms */
tripd *fj;  /* force on solute atom*/
double *pe; /* interaction energy between one water and the solute */
double *vc; /* electrostatic interaction energy between one water and the solute */
int n;/* for the RDF*/
{
double s, sp, f, a, aa, b, bb, q2, r2, r, ir, ir6, EC;
tripd image, rij, sdist;
int i,index;

fj[0].fx = fj[0].fy = fj[0].fz = 0;
fj[1].fx = fj[1].fy = fj[1].fz = 0;/*not used here*/
fi[0].fx = fi[0].fy = fi[0].fz = 0;
fi[1].fx = fi[1].fy = fi[1].fz = 0;
fi[2].fx = fi[2].fy = fi[2].fz = 0;
/***	Determine O-solute image vector			***/

image.fx = -(sdist.fx = ri[0].fx - rj[0].fx);
image.fy = -(sdist.fy = ri[0].fy - rj[0].fy);
image.fz = -(sdist.fz = ri[0].fz - rj[0].fz);
mvimage(&sdist);
image.fx += sdist.fx;
image.fy += sdist.fy;
image.fz += sdist.fz;
r2 = sdist.fx*sdist.fx +  sdist.fy*sdist.fy +  sdist.fz*sdist.fz;
if ( r2 >= swSoMax )
	return;
if ( r2 <= swSoMin ){
	s = 1.0;
	sp = 0.0;
 }
else
	swtchSolu(r2 - swSoMin, &s, &sp);

if(r2< rWatCl*rWatCl){
   Nshel[n-1]++;
}
/***	Loop over atoms in H2O  molecule			***/
for (i=0; i< 3; i++){
	rij.fx = ri[i].fx - rj[0].fx + image.fx;
	rij.fy = ri[i].fy - rj[0].fy + image.fy;
	rij.fz = ri[i].fz - rj[0].fz + image.fz;
	r2 = rij.fx*rij.fx +  rij.fy*rij.fy +  rij.fz*rij.fz ;
	a = wslj[i][0].a;
	aa = 12.*a;
	b = wslj[i][0].b;
	bb = 6.*b;
	q2 = wslj[i][0].q;
	ir = 1. / r2;
	ir6 = ir * ir * ir;
	r = sqrt(r2);
	index = r/binRDF;
	if(i == 0 && index < RDFbins){
	   grCl[n][3][index] += 1.0;/*ion - O RDF*/
	   if(r<rWatCl)
	      Clsol[n][3]++;	
	}
	EC = q2/r;
	*vc += EC*s;
	*pe +=  ( a * ir6 - b ) * ir6 + EC;
	f = s*((aa * ir6 - bb ) * ir6 * ir + EC/r2); /* -dV/dr/r	*/
	f /= 2.;
	fi[i].fx += (rij.fx *= f);
	fi[i].fy += (rij.fy *= f);
	fi[i].fz += (rij.fz *= f);
	fj[0].fx -= rij.fx;
	fj[0].fy -= rij.fy;
	fj[0].fz -= rij.fz;

}
if (sp != 0.0){
	f =sp*(*pe)/2;
	fi[0].fx -= (sdist.fx *= f);
	fi[0].fy -= (sdist.fy *= f);
	fi[0].fz -= (sdist.fz *= f);

	fj[0].fx += sdist.fx;
	fj[0].fy += sdist.fy;
	fj[0].fz += sdist.fz;
	*pe *= s;
	}
}
wdip(ri,rj,fi,fj,pe,vc)
tripd *ri;/* ri[0-2], are position of the H2O's atoms */
tripd *rj;/* rj[0] is the position of the Cl atom */
tripd *fi;/* fi[0-2], are the forces on the H2O's atoms */
tripd *fj;/* force on solute atoms*/
double *pe;/* interaction energy between one water and the solute*/
double *vc;/* electrostatic interaction energy between one water and the solute*/
{
double s, sp, f, a, aa, b, bb, q2, r2, r, ir, ir6, EC;
tripd dip[2],image, rij, sdist,Dcent;
int i,j,index;
fj[0].fx = fj[0].fy = fj[0].fz = 0;
fj[1].fx = fj[1].fy = fj[1].fz = 0;
fi[0].fx = fi[0].fy = fi[0].fz = 0;
fi[1].fx = fi[1].fy = fi[1].fz = 0;
fi[2].fx = fi[2].fy = fi[2].fz = 0;
dip[0].fx = pos[natoms-nsolute].fx;
dip[0].fy = pos[natoms-nsolute].fy;
dip[0].fz = pos[natoms-nsolute].fz;
dip[1].fx = rj[0].fx;
dip[1].fy = rj[0].fy;
dip[1].fz = rj[0].fz;
/***	Determine O-solute image vector			***/
Dcent.fx = 0.5*(dip[0].fx+dip[1].fx);
Dcent.fy = 0.5*(dip[0].fy+dip[1].fy);
Dcent.fz = 0.5*(dip[0].fz+dip[1].fz);

image.fx = -(sdist.fx = ri[0].fx - Dcent.fx);
image.fy = -(sdist.fy = ri[0].fy - Dcent.fy);
image.fz = -(sdist.fz = ri[0].fz - Dcent.fz);
mvimage(&sdist);
image.fx += sdist.fx;
image.fy += sdist.fy;
image.fz += sdist.fz;
r2 = sdist.fx*sdist.fx +  sdist.fy*sdist.fy +  sdist.fz*sdist.fz;

if ( r2 >= swSoMax )
	return;
if ( r2 <= swSoMin ){
	s = 1.0;
	sp = 0.0;
 }
else
	swtchSolu(r2 - swSoMin, &s, &sp);

/***	Loop over atoms in H2O  molecule			***/
for (i=0; i< 3; i++){
   for (j=0; j< 2; j++){
	rij.fx = ri[i].fx - dip[j].fx + image.fx;
	rij.fy = ri[i].fy - dip[j].fy + image.fy;
	rij.fz = ri[i].fz - dip[j].fz + image.fz;
	r2 = rij.fx*rij.fx +  rij.fy*rij.fy +  rij.fz*rij.fz ;
	a = wslj[i][j+1].a;
	aa = 12.*a;
	b = wslj[i][j+1].b;
	bb = 6.*b;
	q2 = wslj[i][j+1].q;
        ir = 1. / r2;
	ir6 = ir * ir * ir;
	r = sqrt(r2);
	index = r/binRDF;
	if(i == 0 && j == 0 && index < RDFbins){
	   grCl[0][3][index] += 0.5;/*CH3 -O RDF*/
	   if(r<rWatCl)
	      Clsol[0][3]++;
	}     
	EC = q2/r;
	*vc += EC*s;
	*pe +=  ( a * ir6 - b ) * ir6 + EC;
	f = s*((aa * ir6 - bb ) * ir6 * ir + EC/r2); /* -dV/dr/r	*/
	f /= 2.;
	fi[i].fx += (rij.fx *= f);
	fi[i].fy += (rij.fy *= f);
	fi[i].fz += (rij.fz *= f);

	fj[j].fx -= rij.fx;
	fj[j].fy -= rij.fy;
	fj[j].fz -= rij.fz;
   }  
}
if (sp != 0.0){
	f =sp*(*pe)/2;
	fi[0].fx -= (sdist.fx *= f);
	fi[0].fy -= (sdist.fy *= f);
	fi[0].fz -= (sdist.fz *= f);

	fj[0].fx += 0.5*sdist.fx;
	fj[0].fy += 0.5*sdist.fy;
	fj[0].fz += 0.5*sdist.fz;
	fj[1].fx += 0.5*sdist.fx;
	fj[1].fy += 0.5*sdist.fy;
	fj[1].fz += 0.5*sdist.fz;
	*pe *= s;
	}
}

bion(ri,rj,fi,fj,pe,vc,n)
tripd *ri;  /* ri[0-8], are position of the Br-octane atoms*/
tripd *rj;  /* rj[0] is the position of the ion */
tripd *fi;  /* fi[0-8], are the forces on the Br-octane atoms */
tripd *fj;  /* force on solute atom*/
double *pe; /* interaction energy between one chloroform and the solute */
double *vc; /* electrostatic interaction energy between one chloroform and the solute */
int n; /* for RDF*/
{
double s, sp, f, a, aa, b, bb, q2, r2, r, ir, ir6, EC;
double pe2,f2;
tripd image, rij, sdist;
int i, index;
fj[0].fx = fj[0].fy = fj[0].fz = 0.;
fj[1].fx = fj[1].fy = fj[1].fz = 0.;
for (i=0;i<BrOs;i++){
	fi[i].fx = fi[i].fy = fi[i].fz = 0.;
}

/***	Loop over all atoms in Br-oct for L-J   ***/
for(i=0;i<BrOs;i++){
	rij.fx = ri[i].fx - rj[0].fx;
	rij.fy = ri[i].fy - rj[0].fy;
	rij.fz = ri[i].fz - rj[0].fz;
	mvimage(&rij);
	r2 = rij.fx*rij.fx + rij.fy*rij.fy + rij.fz*rij.fz;
	
	if(r2 >= swSoMax)
	   continue;
	if(r2 <= swSoMin){
	   s = 1.0;
	   sp = 0.0;
	}
	else
	   swtchSolu(r2 - swSoMin, &s, &sp);
	
	ir = 1. / r2;
	ir6 = ir*ir*ir;
	a = bslj[i][0].a;
	aa = 12.*a;
	b = bslj[i][0].b;
	bb = 6.*b;
	pe2 = (a * ir6 - b) * ir6;
	f2 = s*(aa * ir6 - bb) * ir6 * ir; /* -dV/fr/r */
	f2 /= 2.;
	fi[i].fx += rij.fx * f2;
	fi[i].fy += rij.fy * f2;
	fi[i].fz += rij.fz * f2;
	
	fj[0].fx -= rij.fx * f2;
	fj[0].fy -= rij.fy * f2;
	fj[0].fz -= rij.fz * f2;
	
	if (sp != 0.0){
	   f2 = sp*(pe2)/2.;
	   fi[i].fx -= (rij.fx *= f2);
	   fi[i].fy -= (rij.fy *= f2);
	   fi[i].fz -= (rij.fz *= f2);

	   fj[0].fx += rij.fx;
	   fj[0].fy += rij.fy;
	   fj[0].fz += rij.fz;
	   pe2 *= s;
	}
	*pe += pe2;

}

pe2 = 0.;
/***	Determine alpha-C - solute image vector			***/

image.fx = -(sdist.fx = ri[0].fx - rj[0].fx);
image.fy = -(sdist.fy = ri[0].fy - rj[0].fy);
image.fz = -(sdist.fz = ri[0].fz - rj[0].fz);
mvimage(&sdist);
image.fx += sdist.fx;
image.fy += sdist.fy;
image.fz += sdist.fz;
r2 = sdist.fx*sdist.fx +  sdist.fy*sdist.fy +  sdist.fz*sdist.fz;
if ( r2 >= swSoMax )
	return;
if( r2 <= swSoMin ){
	s = 1.0;
	sp = 0.0;
}
else
	swtchSolu(r2 - swSoMin, &s, &sp);

/***	Loop over charged atoms in Br-oct headgroup	***/
for (i=0; i< 2; i++){
	rij.fx = ri[i].fx - rj[0].fx + image.fx;
	rij.fy = ri[i].fy - rj[0].fy + image.fy;
	rij.fz = ri[i].fz - rj[0].fz + image.fz;
	r2 = rij.fx*rij.fx +  rij.fy*rij.fy +  rij.fz*rij.fz ;
//	a = bslj[i][0].a;
//	aa = 12.*a;
//	b = bslj[i][0].b;
//	bb = 6.*b;
	q2 = bslj[i][0].q;

	ir = 1. / r2;
//	ir6 = ir * ir * ir;
	r = sqrt(r2);
	index = r/binRDF;
	if(index < RDFbins)
	   grCl[n][i][index] += 1.0;/*ion - (aC or Br) RDF*/
	if(i == 1 && r < rBrCl)
	   Clsol[n][i]++;
	if(i == 0 && r < raCCl)
	   Clsol[n][i]++;
	EC = q2/r;
	*vc += EC*s;
	pe2 += /* ( a * ir6 - b ) * ir6 + */ EC;
	f = s*(/*(aa * ir6 - bb ) * ir6 * ir +*/ EC/r2); /* -dV/dr/r	*/

	f /= 2.;

	fi[i].fx += (rij.fx *= f);
	fi[i].fy += (rij.fy *= f);
	fi[i].fz += (rij.fz *= f);
	fj[0].fx -= rij.fx;
	fj[0].fy -= rij.fy;
	fj[0].fz -= rij.fz;

}
if (sp != 0.0){
	f =sp*(pe2)/2.;
	fi[0].fx -= (sdist.fx *= f);
	fi[0].fy -= (sdist.fy *= f);
	fi[0].fz -= (sdist.fz *= f);

	fj[0].fx += sdist.fx;
	fj[0].fy += sdist.fy;
	fj[0].fz += sdist.fz;
	pe2 *= s;
}

*pe += pe2;


}
bdip(ri,rj,fi,fj,pe,vc)
tripd *ri;/* ri[0-8], are position of the Br-octane atoms */
tripd *rj;/* rj0 is the position of the  Cl atom */
tripd *fi;/* fi[0-8], are the forces on the Br-octane atoms */
tripd *fj;/* force on solute atoms*/
double *pe;/* interaction energy between one Br-octane and the solute*/
double *vc;/* electrostatic interaction energy between one Br-octane and the solute*/
{
double s, sp, f, a, aa, b, bb, q2, r2, r, ir, ir6, EC;
double pe2;
tripd image, rij, sdist,Dcent,dip[2];
int i,j,index;
fj[0].fx = fj[0].fy = fj[0].fz = 0.;
fj[1].fx = fj[1].fy = fj[1].fz = 0.;
for (i=0;i<BrOs;i++){
	fi[i].fx = fi[i].fy = fi[i].fz = 0.;
}
dip[0].fx = pos[natoms-nsolute].fx;
dip[0].fy = pos[natoms-nsolute].fy;
dip[0].fz = pos[natoms-nsolute].fz;
dip[1].fx = rj[0].fx;
dip[1].fy = rj[0].fy;
dip[1].fz = rj[0].fz;
/***	Determine O-solute image vector			***/
Dcent.fx = 0.5*(dip[0].fx+dip[1].fx);
Dcent.fy = 0.5*(dip[0].fy+dip[1].fy);
Dcent.fz = 0.5*(dip[0].fz+dip[1].fz);

/***	Loop over atoms in Br-octane molecule for L-J	***/
for(i=0;i<BrOs;i++){
   image.fx = -(sdist.fx = ri[i].fx - Dcent.fx);
   image.fy = -(sdist.fy = ri[i].fy - Dcent.fy);
   image.fz = -(sdist.fz = ri[i].fz - Dcent.fz);
   mvimage(&sdist);
   image.fx += sdist.fx;
   image.fy += sdist.fy;
   image.fz += sdist.fz;
   r2 = sdist.fx*sdist.fx +  sdist.fy*sdist.fy +  sdist.fz*sdist.fz;

   if ( r2 >= swSoMax )
	continue;
   if ( r2 <= swSoMin ){
	s = 1.0;
	sp = 0.0;
   }
   else
	swtchSolu(r2 - swSoMin, &s, &sp);
   pe2 = 0.;
   for(j=0; j< 2; j++){
	rij.fx = ri[i].fx - dip[j].fx + image.fx;
	rij.fy = ri[i].fy - dip[j].fy + image.fy;
	rij.fz = ri[i].fz - dip[j].fz + image.fz;
	r2 = rij.fx*rij.fx + rij.fy*rij.fy + rij.fz*rij.fz;
	a = bslj[i][j+1].a;
	aa = 12.*a;
	b = bslj[i][j+1].b;
	bb = 6.*b;
//	q2 = bslj[i][j+1].q;
	
        ir = 1. / r2;
	ir6 = ir * ir * ir;
//	r = sqrt(r2);
//	index = r/binRDF;
//	if (i == 0 && j == 0) SolRdfCF[0][index] += 0.5;/*CH3 - CF RDF*/
//	EC = q2/r;
//	*vc += EC*s;
	pe2 +=  ( a * ir6 - b ) * ir6 /*+ EC*/;
	f = s*((aa * ir6 - bb ) * ir6 * ir /*+ EC/r2*/); /* -dV/dr/r	*/
	f /= 2.;
	fi[i].fx += (rij.fx *= f);
	fi[i].fy += (rij.fy *= f);
	fi[i].fz += (rij.fz *= f);

	fj[j].fx -= rij.fx;
	fj[j].fy -= rij.fy;
	fj[j].fz -= rij.fz;
   }  
   if(sp != 0.0){
	f =sp*(pe2)/2.;
	fi[0].fx -= (sdist.fx *= f);
	fi[0].fy -= (sdist.fy *= f);
	fi[0].fz -= (sdist.fz *= f);

	fj[0].fx += 0.5*sdist.fx;
	fj[0].fy += 0.5*sdist.fy;
	fj[0].fz += 0.5*sdist.fz;
	fj[1].fx += 0.5*sdist.fx;
	fj[1].fy += 0.5*sdist.fy;
	fj[1].fz += 0.5*sdist.fz;
	pe2 *= s;
   }
   *pe += pe2;
}
/***	Now get coulombic interactions	***/
image.fx = -(sdist.fx = ri[0].fx - Dcent.fx);
image.fy = -(sdist.fy = ri[0].fy - Dcent.fy);
image.fz = -(sdist.fz = ri[0].fz - Dcent.fz);
mvimage(&sdist);
image.fx += sdist.fx;
image.fy += sdist.fy;
image.fz += sdist.fz;
r2 = sdist.fx*sdist.fx +  sdist.fy*sdist.fy +  sdist.fz*sdist.fz;

if ( r2 >= swSoMax )
	return;
if ( r2 <= swSoMin ){
	s = 1.0;
	sp = 0.0;
}
else
	swtchSolu(r2 - swSoMin, &s, &sp);
pe2 = 0.;
for(i=0;i<2;i++){
   for(j=0; j< 2; j++){
	rij.fx = ri[i].fx - dip[j].fx + image.fx;
	rij.fy = ri[i].fy - dip[j].fy + image.fy;
	rij.fz = ri[i].fz - dip[j].fz + image.fz;
	r2 = rij.fx*rij.fx +  rij.fy*rij.fy +  rij.fz*rij.fz ;
//	a = bslj[i][j+1].a;
//	aa = 12.*a;
//	b = bslj[i][j+1].b;
//	bb = 6.*b;
	q2 = bslj[i][j+1].q;
	
        ir = 1. / r2;
//	ir6 = ir * ir * ir;
	r = sqrt(r2);
	if(j==0){
	   index = r/binRDF;
	   if (index < RDFbins)
	      grCl[0][i][index] += 0.5;/*CH3 - (aC or Br) RDF*/
	   if(i==1 && r < rBrCl)
	      Clsol[0][i]++;
	   if(i==0 && r < raCCl)
	      Clsol[0][i]++;
	}
	EC = q2/r;
	*vc += EC*s;
	pe2 +=  /*( a * ir6 - b ) * ir6 +*/ EC;
	f = s*(/*(aa * ir6 - bb ) * ir6 * ir +*/ EC/r2); /* -dV/dr/r	*/
	f /= 2.;
	fi[i].fx += (rij.fx *= f);
	fi[i].fy += (rij.fy *= f);
	fi[i].fz += (rij.fz *= f);

	fj[j].fx -= rij.fx;
	fj[j].fy -= rij.fy;
	fj[j].fz -= rij.fz;
   }
}  
if(sp != 0.0){
	f =sp*(pe2)/2;
	fi[0].fx -= (sdist.fx *= f);
	fi[0].fy -= (sdist.fy *= f);
	fi[0].fz -= (sdist.fz *= f);

	fj[0].fx += 0.5*sdist.fx;
	fj[0].fy += 0.5*sdist.fy;
	fj[0].fz += 0.5*sdist.fz;
	fj[1].fx += 0.5*sdist.fx;
	fj[1].fy += 0.5*sdist.fy;
	fj[1].fz += 0.5*sdist.fz;
	pe2 *= s;
}
*pe += pe2;
}

cdion(ri,rj,fi,fj,pe,vc,n,m)
tripd *ri;  /* ri[0] is the position of the b-CD atom*/
tripd *rj;  /* rj[0] is the position of the ion */
tripd *fi;  /* fi[0] is the force on the b-CD atom */
tripd *fj;  /* force on solute atom*/
double *pe; /* interaction energy between b-CD atom and the solute */
double *vc; /* electrostatic interaction energy between b-CD atom and the solute */
int n;/* for the RDF*/
int m;/* index of b-CD atom (0-20) */
{
double /*s, sp,*/ f, a, aa, b, bb, q2, r2, r, ir, ir6, EC;
tripd image, rij, sdist;
int index;

fj[0].fx = fj[0].fy = fj[0].fz = 0.;
//fj[1].fx = fj[1].fy = fj[1].fz = 0;/*not used here*/
fi[0].fx = fi[0].fy = fi[0].fz = 0.;
//fi[1].fx = fi[1].fy = fi[1].fz = 0;
//fi[2].fx = fi[2].fy = fi[2].fz = 0;
/***	Determine O-solute image vector			***/

//image.fx = -(sdist.fx = ri[0].fx - rj[0].fx);
//image.fy = -(sdist.fy = ri[0].fy - rj[0].fy);
//image.fz = -(sdist.fz = ri[0].fz - rj[0].fz);
//mvimage(&sdist);
//image.fx += sdist.fx;
//image.fy += sdist.fy;
//image.fz += sdist.fz;
//r2 = sdist.fx*sdist.fx +  sdist.fy*sdist.fy +  sdist.fz*sdist.fz;
/*
if ( r2 >= swSoMax )
	return;
if ( r2 <= swSoMin ){
	s = 1.0;
	sp = 0.0;
 }
else
	swtchSolu(r2 - swSoMin, &s, &sp);

*/
rij.fx = ri[0].fx - rj[0].fx;
rij.fy = ri[0].fy - rj[0].fy;
rij.fz = ri[0].fz - rj[0].fz;
mvimage(&rij);
r2 = rij.fx*rij.fx +  rij.fy*rij.fy +  rij.fz*rij.fz;
a = cdslj[m][0].a;
aa = 12.*a;
b = cdslj[m][0].b;
bb = 6.*b;
q2 = cdslj[m][0].q;
ir = 1. / r2;
ir6 = ir * ir * ir;
r = sqrt(r2);
index = r/binRDF;
if((m == 14 || m == 17 || m ==19)  && index < RDFbins){
   grCl[n][2][index] += 1.0;/*ion - bCD hydroxyl O RDF*/
   if(r < rOHCl)
      Clsol[n][2]++;
}
EC = q2/r;
*vc += EC/**s*/;
*pe +=  ( a * ir6 - b ) * ir6 + EC;
f = ((aa * ir6 - bb ) * ir6 * ir + EC/r2); /* -dV/dr/r	*/
f /= 2.;
fi[0].fx += (rij.fx *= f);
fi[0].fy += (rij.fy *= f);
fi[0].fz += (rij.fz *= f);

fj[0].fx -= rij.fx;
fj[0].fy -= rij.fy;
fj[0].fz -= rij.fz;

}
cddip(ri,rj,fi,fj,pe,vc,m)
tripd *ri;/* ri[0] is the position of the b-CD atoms */
tripd *rj;/* rj[0] is the position of the Cl atom */
tripd *fi;/* fi[0] are the forces on the b-CD atoms */
tripd *fj;/* force on solute atoms*/
double *pe;/* interaction energy between b-CD atom and the solute*/
double *vc;/* electrostatic interaction energy between b-CD atom and the solute*/
int m; /* insex of b-CD atom (0-20) */
{
//fprintf(stderr,"in cddip... m = %d\n",m);
double /*s, sp,*/ f, a, aa, b, bb, q2, r2, r, ir, ir6, EC;
tripd dip[2],rij/*,image,sdist,Dcent*/;
int j,index;
fj[0].fx = fj[0].fy = fj[0].fz = 0.;
fj[1].fx = fj[1].fy = fj[1].fz = 0.;
fi[0].fx = fi[0].fy = fi[0].fz = 0.;
//fi[1].fx = fi[1].fy = fi[1].fz = 0;
//fi[2].fx = fi[2].fy = fi[2].fz = 0;
dip[0].fx = pos[natoms-nsolute].fx;
dip[0].fy = pos[natoms-nsolute].fy;
dip[0].fz = pos[natoms-nsolute].fz;
dip[1].fx = rj[0].fx;
dip[1].fy = rj[0].fy;
dip[1].fz = rj[0].fz;
/***	Determine O-solute image vector			***/
/*
Dcent.fx = 0.5*(dip[0].fx+dip[1].fx);
Dcent.fy = 0.5*(dip[0].fy+dip[1].fy);
Dcent.fz = 0.5*(dip[0].fz+dip[1].fz);

image.fx = -(sdist.fx = ri[0].fx - Dcent.fx);
image.fy = -(sdist.fy = ri[0].fy - Dcent.fy);
image.fz = -(sdist.fz = ri[0].fz - Dcent.fz);
mvimage(&sdist);
image.fx += sdist.fx;
image.fy += sdist.fy;
image.fz += sdist.fz;
r2 = sdist.fx*sdist.fx +  sdist.fy*sdist.fy +  sdist.fz*sdist.fz;

if ( r2 >= swSoMax )
	return;
if ( r2 <= swSoMin ){
	s = 1.0;
	sp = 0.0;
 }
else
	swtchSolu(r2 - swSoMin, &s, &sp);
*/
/***	Loop over dipole			***/
//for (i=0; i< 3; i++){
for(j=0;j<2;j++){
   rij.fx = ri[0].fx - dip[j].fx;
   rij.fy = ri[0].fy - dip[j].fy;
   rij.fz = ri[0].fz - dip[j].fz;
   mvimage(&rij);
   r2 = rij.fx*rij.fx + rij.fy*rij.fy + rij.fz*rij.fz;
   a = cdslj[m][j+1].a;
   aa = 12.*a;
   b = cdslj[m][j+1].b;
   bb = 6.*b;
   q2 = cdslj[m][j+1].q;
   ir = 1. / r2;
   ir6 = ir * ir * ir;
   r = sqrt(r2);
   index = r/binRDF;
   if(j==0){
      if((m == 14 || m == 17 || m ==19)  && index < RDFbins){
         grCl[0][2][index] += 0.5;/*CH3 - bCD hydroxyl O RDF*/
         if(r < rOHCl)
            Clsol[0][2]++;
      }
   }
   EC = q2/r;
   *vc += EC/**s*/;
   *pe +=  ( a * ir6 - b ) * ir6 + EC;
   f = ((aa * ir6 - bb ) * ir6 * ir + EC/r2); /* -dV/dr/r	*/
   f /= 2.;
   fi[0].fx += (rij.fx *= f);
   fi[0].fy += (rij.fy *= f);
   fi[0].fz += (rij.fz *= f);

   fj[j].fx -= rij.fx;
   fj[j].fy -= rij.fy;
   fj[j].fz -= rij.fz;
}  
/*
if (sp != 0.0){
	f =sp*(*pe)/2;
	fi[0].fx -= (sdist.fx *= f);
	fi[0].fy -= (sdist.fy *= f);
	fi[0].fz -= (sdist.fz *= f);

	fj[0].fx += 0.5*sdist.fx;
	fj[0].fy += 0.5*sdist.fy;
	fj[0].fz += 0.5*sdist.fz;
	fj[1].fx += 0.5*sdist.fx;
	fj[1].fy += 0.5*sdist.fy;
	fj[1].fz += 0.5*sdist.fz;
	*pe *= s;
	}
*/
}


/*
 *	swtchSolu:
 *
 *	This subroutine generates the switching function and its derivative
 *	given the value of zz.  it is normalized to the size of the box.
 */

swtchSolu(zz,s,sp)
double	zz, *s, *sp;
	{
	double	zz2;
	zz2 = zz * zz;
	*s = swSocoef[0]+zz*zz2*(swSocoef[1]+zz*swSocoef[2]+zz2*swSocoef[3]);
	*sp = zz2 * (pswSocoef[0] + zz*pswSocoef[1] + zz2*pswSocoef[2]);
	}


EVBwindow(){

        int     i;
        double  totMassol,totMassys,zpos,zeta,deriv,fabs(),pow(),relMassol[3];
        tripd   sysCm;
	static int 	init_print = 1;

V_sn2w = totMassol = 0.;
for(i=0;i<nsolute;i++) totMassol += mass[natoms-nsolute+i]; /* total mass of the solute */
for(i=0;i<nsolute;i++) relMassol[i] = mass[natoms-nsolute+i]/totMassol; /* relative mass of solute atoms */

sysCm.fx = sysCm.fy = sysCm.fz = totMassys = 0.;
for(i=0;i<natoms;i++){
        sysCm.fx += mass[i] * pos[i].fx;
        sysCm.fy += mass[i] * pos[i].fy;
        sysCm.fz += mass[i] * pos[i].fz;
        totMassys += mass[i];
        }
sysCm.fx /= totMassys; /* system center of mass */
sysCm.fy /= totMassys;
sysCm.fz /= totMassys;
if (init_print == 1 || fabs(sysCm.fz - osysCmz) > 0.001){
        fprintf(stderr,"system C.O.M = (x=%f y=%f z=%f)\n",sysCm.fx,sysCm.fy,sysCm.fz);
        init_print = 0;
}

osysCmz = sysCm.fz;
zpos = 0.;
for(i=0;i<nsolute;i++) zpos += relMassol[i]*pos[natoms-nsolute+i].fz;
sn2z = zpos;
zpos -= (Ecenter_w + sysCm.fz);
zeta = fabs(zpos) - Ewidth_w / 2.;
//fprintf(stderr,"sysCm.fz = %f, Ecenter_w = %f, zpos = %f, zeta = %f, Epot_w = %f, Epower_w = %f\n",sysCm.fz,Ecenter_w,zpos,zeta,Epot_w,Epower_w);
if(zeta < 0.) return;
else{
        if(zeta > Ewidth_w){
                fprintf(stderr,"window:too much outside the window\n");
                exit(1);
                }
        else{
                V_sn2w = Epot_w*pow(zeta,Epower_w-1);
                deriv = V_sn2w*Epower_w*sgn(zpos);
                V_sn2w *= zeta;
                for(i=0;i<nsolute;i++) force[natoms-nsolute+i].fz -= relMassol[i]*deriv;
                for(i=0;i<natoms;i++) force[i].fz += (deriv * mass[i]/totMassys);
                VINT += V_sn2w;
                }
        }
}

#define TETHCON 1.0
// re-test BCDz tethering!! 
BCDsolTether(){

int i,j,nsolvent;
double sn2mass,MCD,deriv,ccr;
double zpos,zeta,gamma,platMod,solBCDmod,Wc,W_width;
tripd sn2com,cs,cs_hat;

nsolvent = natoms - nBCD*BCDs - nsolute;
Wc = solBCD_c;
W_width = solBCD_w;
V_teth = V_solB = 0.0;
sn2mass = 0.0;
sn2com.fx = sn2com.fy = sn2com.fz = 0.0;
for(i=nsolute;i>0;i--){
   sn2com.fx += pos[natoms-i].fx*mass[natoms-i];
   sn2com.fy += pos[natoms-i].fy*mass[natoms-i];
   sn2com.fz += pos[natoms-i].fz*mass[natoms-i];
   sn2mass += mass[natoms-i];
}
//fprintf(stderr,"sn2 mass = %f\n",sn2mass);
// sn2 CoM
sn2com.fx /= sn2mass;
sn2com.fy /= sn2mass;
sn2com.fz /= sn2mass;

//BCDcom -> sn2com vector
cs.fx = sn2com.fx - BCDcom.fx;
cs.fy = sn2com.fy - BCDcom.fy;
cs.fz = sn2com.fz - BCDcom.fz;
mvimage(&cs);
// gamma calculated before making comcom a unit vector
rBCDsol = gamma = BCDz.fx*cs.fx+BCDz.fy*cs.fy+BCDz.fz*cs.fz;
// ccr is com-com radius
ccr = sqrt(sq(cs.fx)+sq(cs.fy)+sq(cs.fz));
cs_hat.fx/=ccr;
cs_hat.fy/=ccr;
cs_hat.fz/=ccr;

cosccz = BCDz.fx*cs_hat.fx + BCDz.fy*cs_hat.fy + BCDz.fz*cs_hat.fz;

MCD = 0.0; 
for(i=0;i<BCDs;i++){
   MCD += mass[nsolvent+i];
}

// sn2 - BCD: apply biasing 
V_solB =  -parabA*pow(gamma-parabB,2.0);
deriv = 2.0*parabA*(gamma-parabB);

for(i=nsolute;i>0;i--){
   force[natoms-i].fx += deriv*BCDz.fx*mass[natoms-i]/sn2mass;
   force[natoms-i].fy += deriv*BCDz.fy*mass[natoms-i]/sn2mass;
   force[natoms-i].fz += deriv*BCDz.fz*mass[natoms-i]/sn2mass;
}
for(i=0;i<BCDs;i++){
   j = nsolvent + i;
   if(i%21==3 || i%21==5){
      force[j].fx -= deriv*(BCDz.fx*mass[j]/MCD + cs.fx*mass[j]/(capR*Ml) - BCDz.fx*gamma*mass[j]/(capR*Ml));
      force[j].fy -= deriv*(BCDz.fy*mass[j]/MCD + cs.fy*mass[j]/(capR*Ml) - BCDz.fy*gamma*mass[j]/(capR*Ml));
      force[j].fz -= deriv*(BCDz.fz*mass[j]/MCD + cs.fz*mass[j]/(capR*Ml) - BCDz.fz*gamma*mass[j]/(capR*Ml));
   }
   else if(i%21==9 || i%21==16){
      force[j].fx -= deriv*(BCDz.fx*mass[j]/MCD - cs.fx*mass[j]/(capR*Ms) + BCDz.fx*gamma*mass[j]/(capR*Ms));
      force[j].fy -= deriv*(BCDz.fy*mass[j]/MCD - cs.fy*mass[j]/(capR*Ms) + BCDz.fy*gamma*mass[j]/(capR*Ms));
      force[j].fz -= deriv*(BCDz.fz*mass[j]/MCD - cs.fz*mass[j]/(capR*Ms) + BCDz.fz*gamma*mass[j]/(capR*Ms));
   }
   else{
      force[j].fx -= deriv*(BCDz.fx*mass[j]/MCD);
      force[j].fy -= deriv*(BCDz.fy*mass[j]/MCD);
      force[j].fz -= deriv*(BCDz.fz*mass[j]/MCD);
   }
}

// sn2 - BCD window potential
zeta = fabs(gamma - Wc) - (W_width/2.0);
if(zeta > 0.0){
      V_teth += TETHCON*zeta*zeta*zeta;
      deriv = -3.0*TETHCON*zeta*zeta*sgn(gamma-Wc);

      for(i=nsolute;i>0;i--){
	 force[natoms-i].fx += deriv*BCDz.fx*mass[natoms-i]/sn2mass;
	 force[natoms-i].fy += deriv*BCDz.fy*mass[natoms-i]/sn2mass;
	 force[natoms-i].fz += deriv*BCDz.fz*mass[natoms-i]/sn2mass;
      }
/*      MCD = 0.0; 
      for(i=0;i<BCDs;i++){
	 MCD += mass[nsolvent+i];
      }*/
      for(i=0;i<BCDs;i++){
	j = nsolvent + i;
	if(i%21==3 || i%21==5){
	 force[j].fx -= deriv*(BCDz.fx*mass[j]/MCD + cs.fx*mass[j]/(capR*Ml) - BCDz.fx*gamma*mass[j]/(capR*Ml));
	 force[j].fy -= deriv*(BCDz.fy*mass[j]/MCD + cs.fy*mass[j]/(capR*Ml) - BCDz.fy*gamma*mass[j]/(capR*Ml));
	 force[j].fz -= deriv*(BCDz.fz*mass[j]/MCD + cs.fz*mass[j]/(capR*Ml) - BCDz.fz*gamma*mass[j]/(capR*Ml));
	}
	else if(i%21==9 || i%21==16){
	 force[j].fx -= deriv*(BCDz.fx*mass[j]/MCD - cs.fx*mass[j]/(capR*Ms) + BCDz.fx*gamma*mass[j]/(capR*Ms));
	 force[j].fy -= deriv*(BCDz.fy*mass[j]/MCD - cs.fy*mass[j]/(capR*Ms) + BCDz.fy*gamma*mass[j]/(capR*Ms));
	 force[j].fz -= deriv*(BCDz.fz*mass[j]/MCD - cs.fz*mass[j]/(capR*Ms) + BCDz.fz*gamma*mass[j]/(capR*Ms));
	}
	else{
	 force[j].fx -= deriv*(BCDz.fx*mass[j]/MCD);
	 force[j].fy -= deriv*(BCDz.fy*mass[j]/MCD);
	 force[j].fz -= deriv*(BCDz.fz*mass[j]/MCD);
	}
      }
}

}

