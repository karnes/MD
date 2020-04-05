#include	<md.h>
#include	<math.h>
#include	<water.h>
#include	<atomtypes.h>
#define 	en_prec 1.0E-5

/*
 *	this subroutine initializes the polynomials for the water routines.
 */
waterInit()
{
if (waterP[0] == 'W'){ /* watts water */
	oo_poly();
	oh_poly();
	hh_poly();
}
else if (waterP[0] == 'S'){ /* spc water */
	return;
}
else {
		fprintf(stderr, "waterForce: water model undefined\n");
		exit(1);
}
}
#define P_WATT
oo_poly()
{
#ifdef	CTLJCE	/**************************************************************/
if	(natoms <= 6)
	{
	fprintf(stderr, "O-O polynomial was not generated\n");
	return;
	}
#endif	
if	((oonext = waterfit(ootbl,OXYGEN,OXYGEN)) > TBLSZ)
	{
	fprintf(stderr, "error in oo_poly:  out of table memory; ");
	fprintf(stderr, "%d too many words!!\n", oonext - TBLSZ);
	fprintf(stderr, "exiting now!!!\n");
	abort();
	}
else
	{
	/*fprintf(stderr, "O-O polynomial ends at: %d; ", oonext - 1);*/
	/*fprintf(stderr, "length is %d words, with %d words left\n",
		oonext, TBLSZ - oonext);*/
	}
}

oh_poly()
{
#ifdef	CTLJCE	/**************************************************************/
if	(natoms <= 6)
	{
	fprintf(stderr, "O-H polynomial was not generated\n");
	return;
	}
#endif
if	((ohnext = waterfit(ohtbl,OXYGEN,HYDROGEN)) > TBLSZ)
	{
	fprintf(stderr, "error in oh_poly:  out of table memory; ");
	fprintf(stderr, "%d too many words!!\n", ohnext - TBLSZ);
	fprintf(stderr, "exiting now!!!\n");
	abort();
	}
else
	{
	/*fprintf(stderr, "O-H polynomial ends at: %d; ", ohnext - 1);*/
	/*fprintf(stderr, "length is %d words, with %d words left\n",
		ohnext, TBLSZ - ohnext);*/
	}
}

hh_poly()
{
#ifdef	CTLJCE	/**************************************************************/
if	(natoms <= 6)
	{
	fprintf(stderr, "H-H polynomial was not generated\n");
	return;
	}
#endif
if	((hhnext = waterfit(hhtbl,HYDROGEN,HYDROGEN)) > TBLSZ)
	{
	fprintf(stderr, "error in hh_poly:  out of table memory; ");
	fprintf(stderr, "%d too many words!!\n", hhnext - TBLSZ);
	fprintf(stderr, "exiting now!!!\n");
	abort();
	}
else
	{
	/*fprintf(stderr, "H-H polynomial ends at: %d; ", hhnext - 1);*/
	/*fprintf(stderr, "length is %d words, with %d words left\n",
		hhnext, TBLSZ - hhnext);*/
	}
}

/*
 *	Non-bond water potential polynomial fitting subroutine
 *	
 *	This routine is designed to fit general functions with piece
 *	wise polynomials.  By changing "POTCALC" to calculate other
 *	functions, it can be made to fit functions other than the water
 *	potentials in this application.  The polynomial coefficients have been
 *	organized in the AP for optimal usage by ljaafe.aps and aafe.aps
 *	and so may have little resemblance to anything rational.
 */

waterfit(buf,tA,tB)
double	buf[];
int	tA, tB;
{
int	i,j,id,intcnt,is,potopt,nchecks;
int	count;
double	getr(), fval;
double	a,b,c,al,be,ga,fer,eer,chksz,trysz,rmin,rmax;
double	forval,r2,r1,ta,tb,tc,td,te,tf,tg,ti,th,tj,tk,ztmp,dz,en,er,r;
double	idz, invr2, invr1;
double	df, f1, fp1, fpp1, f2, fp2, fpp2;

rmin = getr(tA, tB);
potopt = tA * NTYPES + tB;
switch	(potopt)
	{
	case (OXYGEN * NTYPES + OXYGEN):
		rmax = (int)(sqrt(swr2max) + 1.0);
		break;
	case (HYDROGEN * NTYPES + OXYGEN):
	case (OXYGEN * NTYPES + HYDROGEN):
		rmax = (int)(sqrt(swr2max) + 3.0);
		break;
	case (HYDROGEN * NTYPES + HYDROGEN):
		rmax = (int)(sqrt(swr2max) + 5.0);
		break;
	default:
		fprintf(stderr,"No potential defined for interaction");
		fprintf(stderr," %d-%d\n",tA,tB);
		return;
	}
trysz = 2.;
chksz = trysz/64.;
fer = en_prec / KCAL;		/* from fillfile */
eer = en_prec / KCAL;		/* same for energy and force */
buf[0] = rmax * rmax;
r2 = rmax;
r1 = r2 - trysz;
id = 1;
count = intcnt = 0;
/*
 *	This loop calculates the coefficients for each test interval
 *	If the error at this point in the potential is less than the
 *	allowed limits, then the interval is increased until the error
 *	limits are exceeded.
 */
while	(r2 > rmin)
	{
/*
 *	The first order of business is to parameterize the polynomial fit by
 *	using the values of the function and its first and second derivatives
 *	at the extremes of the interval currently under consideration.
 *	POTCALC provides us with the analytical values.
 */
	dz = r2 * r2 - r1 * r1;
	if	(!count)
		{
		potcalc(r2,&f2,&fp2,&fpp2,potopt);
		invr2 = 1.0 / r2;
		fp2 *= .5 * invr2;
		fpp2 = (fpp2 * .5 - fp2) * invr2 * invr2 * .25;
		count = 1;
		}
	potcalc(r1,&f1,&fp1,&fpp1,potopt);
	invr1 = 1.0 / r1;
	fp1 *= .5 * invr1;
	fpp1 = (fpp1 * .5 - fp1) * invr1 * invr1 * .25;
	idz = 1.0 / dz;
	df = f2 - f1;
	tf = f1;
	tg = fp1;
	th = fpp1;
	ti = -((fpp1*dz + fp1) * dz - df) * idz * idz * idz;
	tj = ((fpp1 *dz + fp2 + 2.*fp1) * dz - 3. * df ) *
		idz * idz * idz * idz;
	tk = (((fpp2 - fpp1)*dz - 3. * (fp2 + fp1))*dz + 6. * df) *
		idz * idz * idz * idz * idz;
	ti += (tk*dz - tj) * dz;
	tj -= 2. * tk * dz;
	ta =  2. * tg;
	tb =  4. * th;
	tc =  6. * ti;
	td =  8. * tj;
	te = 10. * tk;
/*	
 *	Now that we have the coefficients needed, we can proceed to check
 *	the precision of our fit by comparing the analytical values (provided
 *	again by POTCALC) to the polynomial approximation at regular small
 *	intervals ("chksz") along our larger interval ("r2 - r1").
 */
	nchecks = (int) (0.1 + (r2 - r1)/chksz);
	for	(ztmp = r1; ztmp < r2; ztmp += chksz)
		{
		r = ztmp;
		dz = ztmp * ztmp - r1 * r1;
		forval = r * (ta + dz * (tb + dz * (tc + dz * (td + dz * te))));
		en = tf + dz * (tg + dz * (th + dz * (ti + dz * (tj + dz*tk))));
		potcalc(r,&a,&b,&c,potopt);
		if	(fabs(forval-b) > 2.*fer || fabs(en-a) > eer)
/*
 *	The error has exceeded our predetermined limits.  Now we must
 *	determine whether this means we have exceeded the limits of the
 *	current interval (either because we have overextended the interval
 *	or because the check was not fine enough) or whether we can extend
 *	this interval further.  
 */
			break;
		}
	if	(ztmp >= r2)
		{
/*
 *	We succeeded getting through the entire interval unscathed. This may be
 *	the fit we were after.  Save it at this point in case it is.
 */
		buf[id] = r1 * r1;
		buf[id + 1] = ta;
		buf[id + 2] = tb;
		buf[id + 3] = tc;
		buf[id + 4] = td;
		buf[id + 5] = te;
		buf[id + 6] = tf;
		buf[id + 7] = tg;
		buf[id + 8] = th;
		buf[id + 9] = ti;
		buf[id + 10] = tj;
		buf[id + 11] = tk;
		if	(r1 > rmin)
			{
/*
 *	We aren't to the overall maximum interval size yet, so go back
 *	and try our luck with a larger interval.
 */
			if	((r1 -= trysz) < rmin)
				r1 = rmin;
			continue;
			}
		}
	else if	(nchecks == (int)(0.1 + trysz/chksz))
		{
/*
 *	WOOPS:  Not only did the polynomial fit break down here, our
 *	check interval is bigger than the interval itself!  So, back up
 *	and try again with smaller parameters.
 */
		trysz *= 0.5;
		if	((r1 = r2 - trysz) <= rmin)
			{
			r1 = rmin;
			trysz = r2 - r1;
			}
		chksz = 0.03125 * trysz;
		continue;
		}
/*
 *	SUCCESS!!  We have identified a valid polynomial for this interval.
 *	The breakdown of the polynomial fit was only due to overextending the
 *	interval.  Since we saved the previous set of coefficients (which
 *	must have worked), we can just throw these ones away, and start over
 *	with a new interval, after doing some accounting.
 */
	intcnt++;
	count = 0;
	r2 = sqrt(buf[id]);
	if	((r1 = r2 - trysz) <= rmin)
		{
		r1 = rmin;
		trysz = r2 - r1;
		}
	if	((id += 12) > TBLSZ)
		return(TBLSZ+1);
	}
/*
 *	Now, calculate the inner cutoff coefficients which take the function
 *	from the innermost boundary to zero in a smooth continuous manner
 *	such that the forces are zero at the origin.
 */
for	(j = 0; j < 12; j++)
	buf[id + j] = 0.0;
ta = buf[id - 8] * (invr2*invr2*invr2*invr2*invr2*invr2*invr2*invr2);
buf[id + 5] = ta;
buf[id + 8] = buf[id - 4]  +  0.1 * buf[id - 8] * r2 * r2;
buf[id + 11] =  - 0.1 * ta;
r1 = 0.;
intcnt++;
return(12*intcnt+1);
}

#ifdef	P_WATT	/**************************************************************/
#define	ALFHH	3.3392
#define	ALFOH	8.75
#define	ALFOO	4.83
#define	AHH	(306.933/KCAL)
#define	AOH	(2.064/KCAL)
#define	AOO	(1.59202e6/KCAL)
#define	C6	(625.45/KCAL)
#define	C8	(2626.79/KCAL)
#define	C10	(17947.5/KCAL)
#define	QSQ	(36.04/KCAL)
#define	RM	1.75
#define	RS	4.69325
#endif
#ifdef	P_SPC	/**************************************************************/
#define	OO_SIGMA	3.16554
#define	OO_EPSIL	0.1554
#define	Q_CONST		332.053
#define HQ		0.41
#define	OQ		(-0.82)
#endif


/*
 *	POTCALC computes analytically the value ("a") and its first ("b")
 *	and second ("c") derivatives for any defined function at a given value
 *	("r") to be returned to the calling routine for use in deriving a
 *	set of polynomial coefficients and/or for comparison of the derived
 *	polynomial fit approximation to the analytical result.
 *	The version below currently supports either the WATTS or SPC non-bond
 *	water potential though other functions can be found elsewhere
 */
potcalc(r,a,b,c,option)
double	r,*a,*b,*c;
int	option;
{
double	t,t1,t2,t3,t4,t5,aref,ep,ir;

ir = 1./ r;
switch	(option)
	{
#ifdef	P_WATT	/**************************************************************/
/****** oxygen-oxygen ******/ 
	case (OXYGEN * NTYPES + OXYGEN):
		{
		register double	ar, cr, cpr, cppr, sr, spr, sppr;
		ar = AOO*exp(-ALFOO * r);
		t3 = ir * ir * ir;
		t4 = t3 * ir;
		t5 = t4 * ir;
		t3 *= t3 * C6;
		t4 *= t4 * C8;
		t5 *= t5 * C10;
		cr = t3 + t4 + t5;
		cpr = -(6.*t3 + 8.*t4 + 10.*t5) * ir;
		cppr = (42.*t3 + 72.*t4 + 110.*t5) * ir * ir;
		if	(r <= RS)
			{
			sr = exp(-(RS * ir - 1.) * (RS * ir - 1.));
			spr = 2.*(RS*ir-1.)*RS*ir*ir;
			sppr = spr*spr + 2.*(-3.*RS*ir+2.)*RS*ir*ir*ir;
			spr *= sr;
			sppr *= sr;
			}
		else
			{
			sr = 1.;
			spr = sppr = 0.;
			}
		*a = ar             + 4.*QSQ*ir    - sr*cr;
		*b = -ALFOO*ar      - 4.*QSQ*ir*ir - (spr*cr + sr*cpr);
		*c = ALFOO*ALFOO*ar + 8.*QSQ*ir*ir*ir
			- (sppr*cr + 2.*spr*cpr + sr*cppr);
		break;
		}
/****** oxygen-hydrogen	******/ 
	case (HYDROGEN * NTYPES + OXYGEN):
	case (OXYGEN * NTYPES + HYDROGEN):
		t2 = 1. / RM;
		t1 = exp(-ALFOH * (r * t2 - 1.0));
		*a = AOH * t1 * (t1 - 2.0) - 2.0 * QSQ * ir;
		*b = -2.0 * AOH * ALFOH * t1 * (t1 - 1.0) * t2
			+ 2.0 * QSQ * ir * ir;
		*c = 4.0*AOH*ALFOH*ALFOH * t1 * (t1 - 0.5) * t2 * t2
			- 4.0 * QSQ * ir * ir * ir;
		break;
/****** hydrogen-hydrogen ******/ 
	case (HYDROGEN * NTYPES + HYDROGEN):
		t1 = AHH * exp(-ALFHH * r);
		*a = t1 + QSQ * ir;
		*b = -ALFHH * t1 - QSQ * ir * ir;
		*c = ALFHH * ALFHH * t1 + 2. * QSQ * ir * ir * ir;
		break;
#endif
#ifdef	P_SPC	/**************************************************************/
 /*	oxygen-oxygen	*/ 
	case (OXYGEN * NTYPES + OXYGEN):
		t = OO_SIGMA * ir;
		t = t * t * t * t * t * t;
		ep = -4. * OO_EPSIL;
/*  commented out to see if moving the potl curve up and down affects things.
		aref = OO_SIGMA / *rmax;
		aref = aref * aref * aref * aref * aref * aref;
		aref = ep * (aref-aref * aref);
		*a = ep * (t-t * t)-aref;
 */	
		*a = ep * (t-t * t);
		*b = ep * (-6. * t + 12. * t * t) * ir;
		*c = ep * (42. * t-156. * t * t) * ir * ir;
		t = Q_CONST * OQ * OQ;
/*  see above
		*a +=  t * ir - t / *rmax;
 */
		*a +=  t * ir;
		*b += -t * ir * ir;
		*c +=  2.0 * t * ir * ir * ir;
		if	(r <= 2.0)
			close = TRUE;
		break;
 /*	oxygen-hydrogen	*/ 
	case (HYDROGEN * NTYPES + OXYGEN):
	case (OXYGEN * NTYPES + HYDROGEN):
		t = Q_CONST * HQ * OQ;
/*   see above
		*a =  t * ir - t / *rmax;
 */
		*a =  t * ir;
		*b = -t * ir * ir;
		*c =  2.0 * t * ir * ir * ir;
		if	(r <= 1.3)
			close = TRUE;
		break;
 /*	hydrogen-hydrogen	*/ 
	case (HYDROGEN * NTYPES + HYDROGEN):
		t = Q_CONST * HQ * HQ;
/*   see above
		*a =  t * ir - t / *rmax;
 */
		*a =  t * ir;
		*b = -t * ir * ir;
		*c =  2.0 * t * ir * ir * ir;
		if	(r <= 1.0)
			close = TRUE;
		break;
#endif	
	default:
		fprintf(stderr,"waterfit: potential #%d not defined\n", option);
		break;
	}
}
/*
 *	The following was taken from the now defunct getradius() routine.
 *	It is used here instead of ljrad() for convenience, so that the
 *	Lennard-Jones parameters do not need to be fetched.
 *	Radii for various atoms
 */

double
getr(atm1, atm2)
register int	atm1, atm2;
{
#define	OXY_RADIUS	(1.5000)
#define	HYD_RADIUS	(1.2500)
#define	O_H_RADIUS	(1.4000)
switch	(atm1)
	{
case OXYGEN:
	switch	(atm2)
		{
	case OXYGEN:
		return(OXY_RADIUS);
	case HYDROGEN:
		return(O_H_RADIUS);
	default:
		fprintf(stderr, "getr: can't find atom #%d's radius\n", atm2);
		abort();
		}
case HYDROGEN:
	switch	(atm2)
		{
	case OXYGEN:
		return(O_H_RADIUS);
	case HYDROGEN:
		return(HYD_RADIUS);
	default:
		fprintf(stderr, "getr: can't find atom #%d's radius\n", atm2);
		abort();
		}
default:
	fprintf(stderr, "getr: can't find atom #%d's radius\n", atm1);
	abort();
	}
}
