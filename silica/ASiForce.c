/*acetonitrile interactions with the silica*/
#include	<md.h>
#include	<system.h>
#include	 <math.h>

ASiForce()
{
int	i, j, k, l, m, n;
double	r2, dedr, delfx, delfy, delfz;
double		eg, s, sp;
tripd		frc, image, sdl, f[6];
double		ljqa();
int 		ACN[90] = { 0 }; /* id of ACN partner in HB */
//void quicksort2(int [nSi], int, int);

/*** reset unique ACN h-bond counter ***/
//uHB = 0;

for	(i = 0; i < nCH3CN; i++){
	k = nCH3OH*3+i*3;

	for	(j = 0; j < nSi; j++,k = nCH3OH*3+i*3){
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

/*****	Loop over atoms ******/
		for	(m = 0; m < 3; m++)
			for	(n = 0; n < 3; n++)
				{
/*****	Determine image vector ******/
				delfx = pos[k+m].fx - pos[l+n].fx + image.fx;
				delfy = pos[k+m].fy - pos[l+n].fy + image.fy;
				delfz = pos[k+m].fz - pos[l+n].fz + image.fz;
				r2 = delfx*delfx + delfy*delfy + delfz*delfz;
				
/****	Routine for determining H-bonds (N--O distance)  *****/
				if(r2 < NOdist2 && m == 1 && n == 0 && tc%HBdt ==0)
				{
					getHbACN(k,l,i,j,delfx,delfy,delfz);
				}
		// record distance from relevant SiOH to nearest ACN N		
/*				if(n==0 && m==1 && tc%HBdt == 0){
					if(j==SiHMO && r2 < SiDtoA*SiDtoA)
					{
						SiDtoA = sqrt(r2);
					}
					if(j==SiOMH && r2 < SiAtoA*SiAtoA)
					{
						SiAtoA = sqrt(r2);
					}
				}
*/
/*****	Get (1/r dV/dr) ******/
				dedr = ljqa(r2,m,n,&eg,j); // added j

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
		k = i*3+nCH3OH*3;
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
		VAS += eg;
		}
	}
/* count unique appearances of ACN molecules in HBs   */
//	if(tc%HBdt == 0){
//	for(j=0;j<nSi;j++)
//	{
//		ACN[j] = SiHB[j][tc/HBdt];
//	}
	/* sort ascending to make counting uniques easy */
//	quicksort2(ACN,0,nSi-1);
	/* count uniques */
//	for(i=0; i< nSi-1;i++)
//	{
//		if(ACN[i] != ACN[i+1]) uHB++;
//	}
//	}

}
/*
void quicksort2(int x[nSi], int first, int last){
	int pivot, j, temp, i;

	if(first<last){
		pivot=first;
		i=first;
		j=last;

		while(i<j) {
			while(x[i]<=x[pivot] && i<last)
				i++;
			while(x[j]>x[pivot])
				j--;
			if(i<j){
				temp=x[i];
				x[i]=x[j];
				x[j]=temp;
			}
		}
	
		temp=x[pivot];
		x[pivot]=x[j];
		x[j]=temp;
		quicksort2(x,first,j-1);
		quicksort2(x,j+1,last);
	}
}

*/

double
ljqa(r2,m,n,eg,j) // added j
double r2, *eg;
int m,n,j;
{
double der,r,ir,ir6, ec,a,b,q,aa,bb;
int index;
	ir = 1. / r2;
	ir6 = ir * ir * ir;
	a = Silj[m+3][n].a;
	b = Silj[m+3][n].b;
//	if(j%SiSkip==0 && SiSkip < nSi)
//		q=0.0;
//	else
		q = Silj[m+3][n].q;
	aa =12*a;
	bb = 6*b;
	r = sqrt(r2);
	index = r/binRDF;
	if (index < 400)
		SiLiqRdf[m+3][n][index] += 1.0;
	ec = q/r;
	*eg +=  ( a * ir6 - b ) * ir6 +ec;
	der = (bb - aa * ir6) * ir6 * ir -ec/r2;
	return(der);
}

void getHbACN(k,l,i,j,Ox,Oy,Oz)
int k,l,i,j;
double Ox,Oy,Oz; /* O(Si)->N(ACN) vector */
{
	tripd SiOHvec;
	double rON, rOH, dot;
	/* get O->H vector */
	SiOHvec.fx = pos[l+2].fx - pos[l].fx;
	SiOHvec.fy = pos[l+2].fy - pos[l].fy;
	SiOHvec.fz = pos[l+2].fz - pos[l].fz;

	rON = sqrt(Ox*Ox + Oy*Oy + Oz*Oz);
	rOH = sqrt(SiOHvec.fx*SiOHvec.fx + SiOHvec.fy*SiOHvec.fy + SiOHvec.fz*SiOHvec.fz);
       	dot = Ox*SiOHvec.fx + Oy*SiOHvec.fy + Oz*SiOHvec.fz;
	if(dot/(rON*rOH) > cosNOH)
	{
		SiHB[j][(int)(tc/HBdt)] = nCH3OH + i;
		HBSitoA++;
		HBcount[3][(int)((pos[l].fz + pos[k].fz)/(2*HBZBinSize))]++;
		HBcount[4][(int)((pos[l].fz + pos[k].fz)/(2*HBZBinSize))]++;
		HBcount[5][(int)((pos[l].fz + pos[k].fz)/(2*HBZBinSize))]++;
	}	
}
