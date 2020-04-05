/*methanol interactions with the silica*/
#include	<md.h>
#include	<system.h>
#include	 <math.h>

MSiForce()
{
//printf("in MSiForce.c\n");
int	i, j, k, l, m, n;
double	r2, dedr, delfx, delfy, delfz, OOr2;
double	eg, s, sp;
tripd	frc, image, sdl, f[6];
double	ljqm();

for	(i = 0; i < nCH3OH; i++){
	k = i*3;

	for	(j = 0; j < nSi; j++,k = i*3){
		l = nCH3OH*3+nCH3CN*3+3*j;
		eg = 0.0;

/*****	Determine image vector for O(methanol) i - O (SiOH) j.  ******/
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
/****	Routine for determining H-bonds  *****/
		if(r2 < OOdist2 && tc % HBdt == 0){
		//	printf("--MSiForce.c -- heading to getHb\n");
			getHb(k,l,j,delfx,delfy,delfz);
		//	printf("--MSiForce.c - just got back from getHb\n");
		}
/*****	Loop over atoms ******/
		for	(m = 0; m < 3; m++)
			for	(n = 0; n < 3; n++)
				{
/*****	Determine image vector ******/
				delfx = pos[k+m].fx - pos[l+n].fx + image.fx;
				delfy = pos[k+m].fy - pos[l+n].fy + image.fy;
				delfz = pos[k+m].fz - pos[l+n].fz + image.fz;
				r2 = delfx*delfx + delfy*delfy + delfz*delfz;
/*				if(tc%HBdt==0 && nCH3OH==1){
					if(m==0 && n==0)
						OOr2 = r2;
				
			// loops to determine distances between MeOH and nearest
			// possible SiOH HB partners
					if(m==0 && n==2){
						if(r2 < SiHtoMO*SiHtoMO){
							SiHtoMO = sqrt(r2);
							SiOtoMOa = sqrt(OOr2);
							SiHMO = j;
						}
					}
					if(m==1 && n==0){
						if(r2 < MHtoSiO*MHtoSiO){
							MHtoSiO = sqrt(r2);
							SiOtoMOd = sqrt(OOr2);
							SiOMH = j;
						}
					}
				}
*/
/*****	Get (1/r dV/dr) ******/
				dedr = ljqm(r2,m,n,&eg,j); // added j
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
		l = j*3+nCH3OH*3+nCH3CN*3;
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
		VMS += eg;
		}
	}
/*
if(tc%HBdt == 0 && nCH3OH ==1){
	if(MeOHa == 1)
		durMeOHa += (double)HBdt / 2. / 1000.; //time in ps
	if(MeOHa == 0)
		durMeOHa = 0.0;
	if(MeOHd == 1)
		durMeOHd += (double)HBdt / 2. / 1000.; //time in ps
	if(MeOHd == 0)
		durMeOHd = 0.0;
	if(durMeOHa >= MeOHtau && durMeOHd >= MeOHtau)
		MeQ = 1;
}
*/
}

double
ljqm(r2,m,n,eg,j)
double r2, *eg;
int m,n,j;
{
double der,r,ir,ir6, ec,a,b,q,aa,bb;
int index;
	ir = 1. / r2;
	ir6 = ir * ir * ir;
	a = Silj[m][n].a;
	b = Silj[m][n].b;
//	if(j%SiSkip==0 && SiSkip < nSi)
//		q=0.0;
//	else
		q = Silj[m][n].q;
	aa =12*a;
	bb = 6*b;
	r = sqrt(r2);
	index = r/binRDF;
        if (index < 400)
                SiLiqRdf[m][n][index] += 1.0;
	ec = q/r;
	*eg +=  ( a * ir6 - b ) * ir6 +ec;
	der = (bb - aa * ir6) * ir6 * ir -ec/r2;
	return(der);
}

void getHb(k,l,j,Ox,Oy,Oz)
int k,l,j;
double Ox,Oy,Oz;  /* vector from O(Si) to O(Meth) */
{
tripd MOHvec, SiOHvec;
double rOO, rSiOH, rMOH, dotMOHO, dotSiOHO;
/* calculate OHO angles. either is < OHO ang, then H bonded */
	/* get O(methanol) -> H(methanol) vector */
		MOHvec.fx = pos[k+1].fx - pos[k].fx;
		MOHvec.fy = pos[k+1].fy - pos[k].fy;
		MOHvec.fz = pos[k+1].fz - pos[k].fz;
	/* get  O(silica) -> H(silica) vector */
		SiOHvec.fx = pos[l+2].fx - pos[l].fx;
		SiOHvec.fy = pos[l+2].fy - pos[l].fy;
		SiOHvec.fz = pos[l+2].fz - pos[l].fz;
	/* O(Si) -> O(M) vector is (Ox,Oy,Oz) */
		dotSiOHO = Ox*SiOHvec.fx + Oy*SiOHvec.fy + Oz*SiOHvec.fz;
		rSiOH = sqrt(SiOHvec.fx*SiOHvec.fx + SiOHvec.fy*SiOHvec.fy + SiOHvec.fz*SiOHvec.fz);
		dotMOHO = -Ox*MOHvec.fx - Oy*MOHvec.fy - Oz*MOHvec.fz;
		rMOH = sqrt(MOHvec.fx*MOHvec.fx + MOHvec.fy*MOHvec.fy + MOHvec.fz*MOHvec.fz);
		rOO = sqrt(Ox*Ox + Oy*Oy + Oz*Oz);
//		cosMa = dotSiOHO/(rOO*rSiOH);
		if(dotSiOHO/(rOO*rSiOH) > cosHOO)
		{	
//			printf("a) in MSiForce getHb: donor = %d, accept = %d, tc = %d, z= %f\n", j, nSi + k/3, tc, (pos[k].fz + pos[l+2].fz)/2);
			if(SiHB[j][(int)(tc/HBdt)] != -1)
			{
				doubleSiD++;
				SiHB2[j][(int)(tc/HBdt)] = k/3; /*H-bond donated by Si*/
			}
			else
			{
				SiHB[j][(int)(tc/HBdt)] = k/3; /*H-bond donated by Si*/
			}
//			printf("a) assigned bond-- element=%f\n",SiMHB[j][nSi + k/3][tc]); 
			HBSitoM++;
			HBcount[0][(int)((pos[k].fz + pos[l+2].fz)/(2*HBZBinSize))]++;
			HBcount[2][(int)((pos[k].fz + pos[l+2].fz)/(2*HBZBinSize))]++;
			HBcount[4][(int)((pos[k].fz + pos[l+2].fz)/(2*HBZBinSize))]++;
			HBcount[5][(int)((pos[k].fz + pos[l+2].fz)/(2*HBZBinSize))]++;
//			MeOHa = 1;
		}
//		cosMd = dotMOHO/(rOO*rMOH); 
		if(dotMOHO/(rOO*rMOH) > cosHOO) 
		{
//			printf("b) in MSiForce getHb: donor = %d, accept = %d, tc = %d, z= %f\n", nSi + k/3, j, tc, (pos[k+1].fz + pos[l].fz)/2);
			if(SiHB[j + nSi][(int)(tc/HBdt)] != -1)
			{
				doubleSiA++;
				SiHB2[j + nSi][(int)(tc/HBdt)] = k/3; /*H-bond donated by MeOH*/
			}
			else
			{
				SiHB[j + nSi][(int)(tc/HBdt)] = k/3; /*H-bond donated by MeOH*/
			}
//			printf("b) assigned bond -- element=%f\n",SiMHB[nSi + k/3][j][tc]);
			HBMtoSi++;
			HBcount[1][(int)((pos[k+1].fz + pos[l].fz)/(2*HBZBinSize))]++;
			HBcount[2][(int)((pos[k+1].fz + pos[l].fz)/(2*HBZBinSize))]++;
			HBcount[5][(int)((pos[k+1].fz + pos[l].fz)/(2*HBZBinSize))]++;
//			MeOHd = 1;
		}
}
