#include <crunch.h>
/* integrate the electronic action-angle variables for the non-adiabatic ICN
 * project
 */

rk_aav(hq)
double hq;
{
double fn(), fq(), fxi(), fat(), rkn1, rkn2, rkq1, rkq2;
find_canon();/* find the correct non-singular system of action_angle variables*/
switch ( transform)
       {
       case 0: rkn1 = fn(aav[0],aav[1]) * hq;
	       rkq1 = fq(aav[0],aav[1]) * hq;
               rkn2 = fn(aav[0]+rkn1,aav[1]+rkq1) * hq;
               rkq2 = fq(aav[0]+rkn1,aav[1]+rkq1) * hq;
	       aav[0] += (rkn1 + rkn2)/2.;
	       aav[1] += (rkq1 + rkq2)/2.;
	       action = aav[0];
	       break;
       case 1: rkn1 = fxi(aav[0],aav[1]) * hq;
	       rkq1 = fat(aav[0],aav[1]) * hq;
               rkn2 = fxi(aav[0]+rkn1,aav[1]+rkq1) * hq;
               rkq2 = fat(aav[0]+rkn1,aav[1]+rkq1) * hq;
	       aav[0] += (rkn1 + rkn2)/2.;
	       aav[1] += (rkq1 + rkq2)/2.;
	       action = (sq(aav[0]) + sq(aav[1]) - 1.)/2.;
	       break;
       case 2: rkn1 = fxi(-aav[0],-aav[1]) * hq;
	       rkq1 = fat(-aav[0],-aav[1]) * hq;
               rkn2 = fxi(-aav[0]-rkn1,-aav[1]-rkq1) * hq;
               rkq2 = fat(-aav[0]-rkn1,-aav[1]-rkq1) * hq;
	       aav[0] += (rkn1 + rkn2)/2.;
	       aav[1] += (rkq1 + rkq2)/2.;
	       action = (3. - sq(aav[0]) - sq(aav[1]) )/2.;
	       break;
       default: 
               ERROR((stderr,"rk_aav: illegal value transform= %d\n",transform),exit);
       }
}

/* find the correct non-singular system of action_angle variables */
 
find_canon()
{
double ac1, ac2, epsilon, sqrt(), sin(), cos();

epsilon = 0.0001; 
ac1 = sqrt (2 * action + 1.);
ac2 = sqrt (3. - 2 * action);
switch ( transform)
       {
       case 0: if ( action < (-0.5 + epsilon) )
		  {
		    transform = 1;
		    aav[0] = ac1 * sin(aav[1]);
		    aav[1] = ac1 * cos(aav[1]);
                  }
               else 
		    if ( action > (1.5 - epsilon) )
                       {
		         transform = 2;
		         aav[0] = -ac2 * sin(aav[1]);
		         aav[1] = ac2 * cos(aav[1]);
		       }
	       break;
       case 1: if ( action > (1.5 - epsilon) )
		  {
		    transform = 2;
		    aav[0] = -aav[0]*ac2/ac1;
		    aav[1] = aav[1]*ac2/ac1;
                  }
	       break;
       case 2: if ( action < (-0.5 + epsilon) )
		  {
		    transform = 1;
		    aav[0] = -aav[0]*ac1/ac2;
		    aav[1] = aav[1]*ac1/ac2;
                  }
	       break;
       default: 
               ERROR((stderr,"rk_aav: illegal value transform= %d\n",transform),exit);
       }
}

double fn(x,y)
double x, y;
{
double sqrt(), sin(), val;
val =  (2*x+1.)*(3.-2*x) ;
if (val <= 0. )
   ERROR((stderr,"rk_aav:fn: illegal value of action, reduce h\n"),exit);
val = sqrt ( val ) * sin (y) * VEXCOUP;
return(val);
}

double fq(x,y)
double x, y;
{
double sqrt(), cos(), val;
val =  (2*x+1.)*(3.-2*x) ;
if (val <= 0. )
   ERROR((stderr,"rk_aav:fq: illegal value of action, reduce h\n"),exit);
val = sqrt (val);
val = 2*(1.- 2*x)*cos(y)*VEXCOUP / val + VEXLIN - VEXBEN;
return(val);
}

double fxi(x,y)
double x,y;
{
double sqrt(), val;
val = 4. - sq(x) - sq(y) ;
if (val <= 0. )
   ERROR((stderr,"rk_aav:fxi: illegal value of action, reduce h\n"),exit);
val = sqrt( val );
val = (val - sq(y)/val) * VEXCOUP + y * (VEXLIN - VEXBEN) ;
return(val);
}

double fat(x,y)
double x,y;
{
double sqrt(), val;
val = 4. - sq(x) - sq(y) ;
if (val <= 0. )
   ERROR((stderr,"rk_aav:fat: illegal value of action, reduce h\n"),exit);
val = sqrt( val );
val = x * y * VEXCOUP / val + x * (VEXBEN - VEXLIN );
return(val);
}
