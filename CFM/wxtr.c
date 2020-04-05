#include	<md.h>
#include	<system.h>
#include	<math.h>
//#define CFMden 0.007466

static int rdfN[6] = {1,6,2,9,6,1};
wxtr(fp, initQ, pFreq)
	FILE	*fp;
	int	initQ;	/*initialization flag	*/
	int	pFreq;	/* frequency to print output	*/
{
int flag,i,j,k,l,bin,most,t,tau;
double gfactor1,gfactor2,gfactor3,gKi1,gKi2,gKi3;
double jacob,ir,jr,jr1,gfactor,r,gKi,gKnet,cosi;
double gKnet1, gKnet2,gKnet3;
double ph1,ph2,r1,r2;
tripd d;
//int **stkLife[2];
int *stkTCF[4];
int *stkTCFn[4];
double *dipTCF; // orientational TCF of dipole 
double *c3TCF; //orientational TCF norm to c3 axis
double *chTCF; //orientational TCF, C->H
double *normdip;
double *normc3;
double *normch;
#define stx 7
int totPstx[7];
int *pstx[stx]; //polar stacking analysis array
// 0 = more than 6, 1-6 = num in "stack"
//fprintf(stderr,"in wxtr.c  tdpoint=%d pFreq=%d\n",tdpoint,pFreq);
if (initQ == 1){
/*	if(QCFM[4]==QCFM[3])
	   CFMden = 0.006225;
	else
	   CFMden = 0.007466;*/
	for(k=0;k<(int)(2.0/dStkCos);k++)
	   stackCos[k] = 0;
	apollo = nApollo = 0;
	   	   
	tdpoint = 0;
	for(k=0;k<6;k++)
		for(bin=0;bin<xwall/binRDF;bin++)
			CFMRdf[k][bin] = 0.;
}

tdpoint++;
if ((tdpoint-1) % pFreq == 0 && tdpoint > 1){
/* Print radial distribution functions*/
	gfactor = numDP_x_dataRate*nCFM*2*PI*binRDF*binRDF*binRDF*CFMden;
	fprintf(fp,"r\tCC\tCCl\tCH\tClCl\tClH\tHH\tCoMCoM\n");
	for (i=20;i<500;i++){
		fprintf(fp,"%-5.3f\t",i*binRDF);
		for (k=0;k<6;k++)
		   fprintf(fp,"%-6.4f\t",CFMRdf[k][i]/(rdfN[k]*gfactor*(1.0/3+i*(i+1))));
		   fprintf(fp,"%-6.4f\t",CFMRdf[6][i]/(( (1+(tc/TCFdt))*gfactor/numDP_x_dataRate)*(1.0/3+i*(i+1))));
		fprintf(fp,"\n");
	}
//fprintf(stderr,"wxtr.c - rdf done \n");
/*************************
 *   orientational TCFs
 *************************/
	if((dipTCF = (double *)calloc((int)(numDP_x_dataRate/TCFdt)+1,sizeof(double))) == NULL
	 || (chTCF = (double *)calloc((int)(numDP_x_dataRate/TCFdt)+1,sizeof(double))) == NULL
	 || (c3TCF = (double *)calloc((int)(numDP_x_dataRate/TCFdt)+1,sizeof(double))) == NULL
	 ||(normdip= (double *)calloc((int)(numDP_x_dataRate/TCFdt)+1,sizeof(double))) == NULL
	 ||(normch= (double *)calloc((int)(numDP_x_dataRate/TCFdt)+1,sizeof(double))) == NULL
	 ||(normc3 = (double *)calloc((int)(numDP_x_dataRate/TCFdt)+1,sizeof(double))) == NULL)
		ERROR((stderr,"wxtr: out of core, O-TCFs.\n"),exit);
	for(i=0;i<(int)(numDP_x_dataRate/TCFdt)+1;i++){
		dipTCF[i]=chTCF[i]=c3TCF[i]=normdip[i]=normc3[i]=normch[i]=0.0;
	}
	tcfs(dipVec,dipTCF,normdip);
	tcfs(c3Vec,c3TCF,normc3);
	tcfs(chVec,chTCF,normch);
	// print TCFs
	fprintf(fp,"Orientational TCFs\n");
	fprintf(fp,"t(ps)\tdip\tC3\tn_C3\tnorm_dip\n");
	for(i=0;i<(tc/TCFdt);i+=2){
		fprintf(fp,"%f\t%f\t%f\t%f\t%d\n",(double)(i*TCFdt)*0.5/1000.,dipTCF[i]/normdip[i],chTCF[i]/normch[i],c3TCF[i]/normc3[i],(int)(normdip[i]));
	}

/********************************
 *    polar stacking analysis   *
 *******************************/
	//initialize array
	for(i=0;i<stx;i++){
		if((pstx[i] = (int *)calloc((int)(numDP_x_dataRate/TCFdt)+1,sizeof(int))) == NULL)
			ERROR((stderr,"wxtr: out of core, polar stacking.\n"),exit);
	}
	for(i=0;i<stx;i++){
	   totPstx[i] = 0;
	   for(j=0;j<numDP_x_dataRate/TCFdt+1;j++){
		pstx[i][j] = 0;
	   }
	}
	for(i=0;i<4;i++){
		if((stkLife[i] = (int **)calloc(nCFM,sizeof(int *))) == NULL)
			ERROR((stderr,"wxtr: out of core, polar stacking.\n"),exit);
	}
	for(i=0;i<4;i++){
	   for(j=0;j<nCFM;j++){
		if((stkLife[i][j] = (int *)calloc((int)(numDP_x_dataRate/TCFdt)+1,sizeof(int))) == NULL)
		     ERROR((stderr,"wxtr: out of core, polar stacking.\n"),exit);
	   }
	}
	for(i=0;i<4;i++){
	   if((stkTCF[i] = (int *)calloc((int)(numDP_x_dataRate/TCFdt)+1,sizeof(int))) == NULL)
	     ERROR((stderr,"wxtr: out of core, polar stacking.\n"),exit);
	   if((stkTCFn[i] = (int *)calloc((int)(numDP_x_dataRate/TCFdt)+1,sizeof(int))) == NULL)
	     ERROR((stderr,"wxtr: out of core, polar stacking.\n"),exit)
	}
	for(i=0;i<4;i++)
	   for(k=0;k<(int)(numDP_x_dataRate/TCFdt)+1;k++)
	      stkTCF[i][k] = stkTCFn[i][k] = 0;

	//analyze data
//fprintf(stderr,"wxtr.c - before readStacks.  \n");
	readStacks(pstx);
//fprintf(stderr,"wxtr.c - after readStacks.  \n");
	//print stack data
	fprintf(fp,"Polar stacking\n");
	for(i=0;i<(tc/TCFdt);i++){
//	   fprintf(fp,"%f\t",(double)(i*TCFdt)*0.5/1000.0);
	   for(j=0;j<stx;j++){
	       totPstx[j]+=pstx[j][i];
	   }
	}	
	fprintf(fp,"\nPolar stacking (totals). avgDipDip = %f\n",avgDipDip/avgDipDipn);
	fprintf(fp,">5\t1\t2\t3\t4\t5\tloops\n");
	for(j=0;j<stx;j++){
	    fprintf(fp,"%f\t",(float)(totPstx[j])/((float)(tc)/(float)(TCFdt)));
	}
	fprintf(fp,"\n");
// polar stacking - new algorithm (non-Salzmann)
	for(i=0;i<stx;i++){
	   totPstx[i] = 0;
	   for(j=0;j<numDP_x_dataRate/TCFdt+1;j++){
		pstx[i][j] = 0;
	   }
	}
	readStacks2(pstx);
	fprintf(fp,"Polar stacking - IB group algorithm\n");
	for(i=0;i<(tc/TCFdt);i++){
	   for(j=0;j<stx;j++){
	       totPstx[j]+=pstx[j][i];
	   }
	}	
	fprintf(fp,"\nPolar stacking (IB) (totals) avgDipDip2 = %f, pct.stacked in wedge = %f\n",avgDipDip2/avgDipDip2n,100.0*avgDipDip2n/totInHCC);
	fprintf(fp,">5\t1\t2\t3\t4\t5\tloops\n");
	for(j=0;j<stx;j++){
	    fprintf(fp,"%f\t",(float)(totPstx[j])/((float)(tc)/(float)(TCFdt)));
	}
	fprintf(fp,"\n");
/************************
 *    2-D SDFs          *
 ***********************/
	gfactor = 1./(CFMden * numDP_x_dataRate/TCFdt * nCFM); 
	fprintf(fp,"\nmaps: CH\n");
	for(i=0;i<2*(int)(contRad/contBin)+1;i++){
	   for(j=0;j<(int)(contRad/contBin)+1;j++){
		fprintf(fp,"%f\t",gfactor*CHmap[i][j]/( ( sq((j+1)*contBin) - sq(j*contBin) ) * contBin * PI));
	   }
	   fprintf(fp,"\n");
	}
	fprintf(fp,"\nmaps: CCl\n");
	for(i=0;i<2*(int)(contRad/contBin)+1;i++){
	   for(j=0;j<(int)(contRad/contBin)+1;j++){
		fprintf(fp,"%f\t",(gfactor/3)*CClmap[i][j]/( ( sq((j+1)*contBin) - sq(j*contBin) ) * contBin * PI));
	   }
	   fprintf(fp,"\n");
	}
	fprintf(fp,"\nmaps: CC\n");
	for(i=0;i<2*(int)(contRad/contBin)+1;i++){
	   for(j=0;j<(int)(contRad/contBin)+1;j++){
		fprintf(fp,"%f\t",gfactor*CCmap[i][j]/( ( sq((j+1)*contBin) - sq(j*contBin) ) * contBin * PI));
	   }
	   fprintf(fp,"\n");
	}
	for(k=0;k<5;k++){
	   if(k==0){
		ph1 = 0.0;
		ph2 = dAng*PI/180.0;
		flag = 1;
	   }
	   else if(k==4){
	        ph1 = (180.0-dAng)*PI/180.0;
		ph2 = PI;
		flag = 1;
	   }
	   else{
		ph1 = ((k*45.0)-dAng)*PI/180.0;
		ph2 = ((k*45.0)+dAng)*PI/180.0;
		flag = 2;
	   }
	   fprintf(fp,"\nmaps: dTheta = %3.2f dipTheta%d\n",dAng,45*k);
	   for(i=0;i<2*(int)(contRad/contBin)+1;i++){
	      for(j=0;j<flag*(int)(contRad/contBin)+1;j++){
		ir = fabs(i*contBin - contRad);
		jr = fabs(j*contBin - (flag-1)*contRad);
		r1 = sqrt(sq(ir) + sq(jr));
		r2 = sqrt(sq(ir+contBin) + sq(jr+contBin));
		jacob = /*nCFM*(numDP_x_dataRate/TCFdt)**/fabs((2.0/3.0)*PI*(pow(r2,3.0) - pow(r1,3.0))*(cos(ph1) - cos(ph2)))/*/CFMden*/;
		if(flag==2)
		   fprintf(fp,"%f\t",/*(float)*/(dipMap[k][i][j])/jacob);
		else if (flag==1)
		   fprintf(fp,"%f\t",/*(float)*/(dipMap[k][i][j])/jacob);
		   
	      }
	      fprintf(fp,"\n");
	   }
	}
	
	for(k=0;k<5;k++){
	   if(k==0){
		ph1 = 0.0;
		ph2 = dAngsm*PI/180.0;
		flag = 1;
	   }
	   else if(k==4){
	        ph1 = (180.0-dAngsm)*PI/180.0;
		ph2 = PI;
		flag = 1;
	   }
	   else{
		ph1 = ((k*45.0)-dAngsm)*PI/180.0;
		ph2 = ((k*45.0)+dAngsm)*PI/180.0;
		flag = 2;
	   }
	   fprintf(fp,"\nmaps: dTheta = %3.2f dipTheta%d\n",dAngsm,45*k);
	   for(i=0;i<2*(int)(contRad/contBin)+1;i++){
	      for(j=0;j<flag*(int)(contRad/contBin)+1;j++){
		ir = fabs(i*contBin - contRad);
		jr = fabs(j*contBin - (flag-1)*contRad);
		r1 = sqrt(sq(ir) + sq(jr));
		r2 = sqrt(sq(ir+contBin) + sq(jr+contBin));
		jacob = /*nCFM*(numDP_x_dataRate/TCFdt)**/fabs((2.0/3.0)*PI*(pow(r2,3.0) - pow(r1,3.0))*(cos(ph1) - cos(ph2)))/*/CFMden*/;
		if(flag==2)
		   fprintf(fp,"%f\t",/*(float)*/(dipMapsm[k][i][j])/jacob);
		else if (flag==1)
		   fprintf(fp,"%f\t",/*(float)*/(dipMapsm[k][i][j])/jacob);
		   
	      }
	      fprintf(fp,"\n");
	   }
	}

//	stackLife(stkLife,stkTCF,stkTCFn);
//	int j,i,t,tau;
	for(i=0;i<2;i++)
	{
	   for(j=0;j<nCFM;j++)
	   {
		for(tau=0;tau<(int)(numDP_x_dataRate/TCFdt);tau++)
		{
			for(t=0;t<=tau;t++)
			{
//fprintf(stderr,"i=%d/nCFM=%d/tau=%d/t=%d/stk=%d/",i,j,tau,t,stkLife[i][j][tau]);
				stkTCF[i][t] += stkLife[i][j][tau]*stkLife[i][j][tau-t];
				stkTCFn[i][t] += stkLife[i][j][tau];
			}
		}
	   }
	}
	// calculate lifetime of each discrete stack
	i = 2;
	for(j=0;j<nCFM;j++)
	{
	   for(tau=0;tau<(int)(numDP_x_dataRate/TCFdt);tau++)
	   {
	      for(t=0;t<=tau;t++){
	      //check if in a stack
	         if(pstack[j][tau]!=-1){
		    stkTCFn[i][t]+=1;
	  	    //check if stack intact at time tau-t
		    if(pstack[j][tau]==pstack[j][tau-t]){
		       stkTCF[i][t] += 1;
		    }
		 }
	      }
	   }
	}
//fprintf(stderr,"completed TCF for stacks of n=2\n");
	i = 3;
	for(j=0;j<nCFM;j++)
	{
	   for(tau=0;tau<(int)(numDP_x_dataRate/TCFdt);tau++)
	   {
	      for(t=0;t<=tau;t++){
	      //check if in a stack of 3
	         if(pstack[j][tau]!=-1){
		    if(pstack[pstack[j][tau]][tau]!=-1){
		       stkTCFn[i][t] += 1;
	  	       //check if stack of 3 intact at time tau-t
		       if(pstack[j][tau]==pstack[j][tau-t]){
			 if(pstack[pstack[j][tau]][tau]==pstack[pstack[j][tau-t]][tau-t]){
		            stkTCF[i][t] += 1;
			 }
		       }
		    }
	         }
	      }
	   }
	}
//fprintf(stderr,"completed TCF for stacks of n=3\n");

//fprintf(stderr,"got stack lifetimes OK. tc = %d. \n",tc);
	// print polar stack lifetime CFs
	fprintf(fp,"Polar Stack Lifetime CFs\n");
	fprintf(fp,"t(ps)\tn>=2\tn>=3\tnorm_n>=2\tn>=2\tn>=3\tnorm2\n");
	for(i=0;i<(tc/TCFdt);i+=2){
		fprintf(fp,"%f\t%f\t%f\t%d\t%f\t%f\t%d\n",(double)(i*TCFdt)*0.5/1000.,
			(double)(stkTCF[0][i])/(double)(stkTCFn[0][i]),
			(double)(stkTCF[1][i])/(double)(stkTCFn[1][i]),stkTCFn[0][i],
			(double)(stkTCF[2][i])/(double)(stkTCFn[2][i]),
			(double)(stkTCF[3][i])/(double)(stkTCFn[3][i]),stkTCFn[2][i]);
	}
	// print extra stuff for Apollo versus polar stack stuff
	fprintf(fp,"Apollo vs. polar stacks.\n");
	fprintf(fp,"Fraction of stk in Apollo (cos(dip-dip angle)>%5.4f\n",cosApollo);
	fprintf(fp,"%f\n",(float)(apollo)/(float)(nApollo));
	fprintf(fp,"cos(dipdip) of stacked CFMs\n");
	fprintf(fp,"cos(mu-mu)\tcount\n");
	for(i=0;i<(int)(2.0/dStkCos);i++){
		fprintf(fp,"%4.3f\t%d\n",(double)(i)*dStkCos-1.0,stackCos[i]);
	}

free(dipTCF);
free(chTCF);
free(c3TCF);
free(normdip);
free(normc3);
free(normch);
for(i=0;i<stx;i++){
	free(pstx[i]);
}
for(i=0;i<4;i++){
   free(stkTCF[i]);
   free(stkTCFn[i]);
   free(stkLife[i]);
}

fflush(fp);
}
}
void readStacks(int **pstx)
{
   int i,j,k,m;
   int s[6];
   int bottom,top;
   int stack;
   int current;
//   for(i=0;i<2000;i++){
//      printStx[i] = -1;
//   }
   for(i=0;i<(int)(numDP_x_dataRate/TCFdt)+1;i++){
	for(j=0;j<nCFM;j++){
	    stack=top=0;
	    bottom=1;

	    //is j above any other CFM?
	    for(k=0;k<nCFM;k++){
	       if(pstack[k][i]==j){
		  bottom=0;
		  k=nCFM+1;
	       }
	    }
	    if(bottom==1){
	       s[0]=s[1]=s[2]=s[3]=s[4]=s[5]=-2;
	       current = j;
	       stack++;
	       s[0]=current;
	       while(top==0){
		 //find CFM above
		 if(pstack[current][i]!=-1 && stack < 7){
	            current=pstack[current][i];
		    if(current==s[0]||current==s[1]||current==s[2]||current==s[3]||current==s[4]||current==s[5]){
			    //loop found!
			    top=1;
			    pstx[6][i]++;
		    }
		    else{
		    	    s[stack]=current;
			    stack++;
		    }
		 }
		 else{
		    top=1;
		 }
	       }
	       //populate array for polar stacking lifetime TCF
	       if(stack>=2){
	          for(m=0;m<6;m++){
	             if(s[m]!=-2){
		        stkLife[0][s[m]][i] = 1;
		     }
	          }
	       }
	       if(stack>=3){
	          for(m=0;m<6;m++){
	             if(s[m]!=-2){
		        stkLife[1][s[m]][i] = 1;
//fprintf(stderr,"stackLife: 1,%d,%d/",s[m],i);

		     }
	          }
	       }
	    }
	    if(stack< 7 && stack > 0){
		    pstx[stack][i]++;

	    }
	    else if(stack >= 7){
		    pstx[0][i]++;
	    }
       }
   }
}

void tcfs(tripd **vec, double *tcf, double *norm)
{
	int i,t,tau;
	double dot;
	for(i=0;i<nCFM;i++)
	{
		for(tau=0;tau<(int)(tc/TCFdt)+1;tau++)
		{
			for(t=0;t<tau;t++)
			{
				dot = vec[i][tau].fx * vec[i][tau-t].fx +
				      vec[i][tau].fy * vec[i][tau-t].fy +
				      vec[i][tau].fz * vec[i][tau-t].fz;
				tcf[t] += dot;
				norm[t] += 1.0;
			}
		}
	}
}

void readStacks2(int **pstx)
{
   int i,j,k;
   int s[6];
   int bottom,top;
   int stack;
   int current;
   for(i=0;i<(int)(numDP_x_dataRate/TCFdt)+1;i++){
	for(j=0;j<nCFM;j++){
	    stack=top=0;
	    bottom=1;

	    //is j above any other CFM?
	    for(k=0;k<nCFM;k++){
	       if(pstack2[k][i]==j){
		  bottom=0;
		  k=nCFM+1;
	       }
	    }
	    if(bottom==1){
	       s[0]=s[1]=s[2]=s[3]=s[4]=s[5]=-2;
	       current = j;
	       stack++;
	       s[0]=current;
	       while(top==0){
		 //find CFM above
		 if(pstack2[current][i]!=-1 && stack < 7){
	            current=pstack2[current][i];
		    if(current==s[0]||current==s[1]||current==s[2]||current==s[3]||current==s[4]||current==s[5]){
			    //loop found!
			    top=1;
			    pstx[6][i]++;
		    }
		    else{
		    	    s[stack]=current;
			    stack++;
		    }
		 }
		 else{
		    top=1;
		 }
	       }
	    }
	    if(stack< 7 && stack > 0){
		    pstx[stack][i]++;

/*		    if(stack==3 && i==(int)(numDP_x_dataRate/TCFdt)){
		       printStx[(pstx[3][i]-1)*3] = s[0];
		       printStx[(pstx[3][i]-1)*3+1] = s[1];
		       printStx[(pstx[3][i]-1)*3+2] = s[2];
		    }
*/
	    }
	    else if(stack >= 7){
		    pstx[0][i]++;
	    }
	    if(stack>1){
	//	    fprintf(stderr,"%d %d %d %d %d \n",s[0],s[1],s[2],s[3],s[4],s[5]);
	    }
       }
   }
}

