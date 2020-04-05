/*
 * Conjugate Gradient Minimizer adapted from "Numerical Recipes"
 * by Press, Flannery, Teukolsky, and Vetterling,
 * Cambridge, 1986.
 */
#include	<md.h>

#define	TINY	1.e-20 	/*special case constant to prevent division by 0*/
#define	TOL	1.0e-3 		/*tolerance*/
#define	GOLD	1.618034	/* golden mean constant*/
#define	GLIMIT	10.	/*?????*/
#define	EPS	1e-10
#define	SCALE	1.

int	it_com;
double	df_com;
tripd	*ptmp, *ftmp, *pcom, *fcom, *gg, *hh;



/*
 * performs Fletcher-Reeves-Polak-Ribiere minimization of the
 * potential energy function, using its gradient.  The convergence
 * tolerance of the function is given as ftol.  Returned is the minimum
 * value of the function (as fret).  The subroutine "linmin" is called
 * to perform line minimizations.
 */

minimize_cg(ftol, iter)	/*iter = number of iterations*/
double	ftol;
int	iter;
	{
	double	fret, fp, g, dg, gam, fmax, rtol, fmax1;
	double	linmin(), getf(), sqrt();
	int	i, jj;
	if	((gg = (tripd *)malloc(natoms*sizeof(tripd))) == NULL
			||(hh  = (tripd *)malloc(natoms*sizeof(tripd)))==NULL
			||(ptmp =(tripd *)malloc(natoms*sizeof(tripd)))==NULL
			||(ftmp =(tripd *)malloc(natoms*sizeof(tripd)))==NULL
			||(pcom =(tripd *)malloc(natoms*sizeof(tripd)))==NULL
			||(fcom =(tripd *)malloc(natoms*sizeof(tripd)))==NULL)
		ERROR((stderr, "minimize:  out of core\n"),exit);


	it_com = iter;
	rtol = sqrt(ftol);
	ftol *= 0.5;
	/*fp = V*/
	fp = getf(pos, force);
	for	(i = 0; i < natoms; i++)
		{
		hh[i].fx = gg[i].fx = force[i].fx;
		hh[i].fy = gg[i].fy = force[i].fy;
		hh[i].fz = gg[i].fz = force[i].fz;
		}
	/*
	 * loop over iterations
	 */
	do	{
		iter = it_com - 4;
		fret = linmin(pos, force, ftol);
		fmax = getf(pos, force);
		fmax1 = fmax = 0.;
		for	(i = 0; i < natoms; i++)
			{
			fmax1 = max(fmax, abs(force[i].fx));
			fmax1 = max(fmax1, abs(force[i].fy));
			fmax1 = max(fmax1, abs(force[i].fz));
			if	(fmax1 > fmax)
				{
				jj = i;
				fmax = fmax1;
				}
			}
		fprintf(stderr, "i = %d F = %.13g (%d)\n", it_com, fmax, jj);
		if	((abs(fret - fp) <= ftol*(abs(fret)+abs(fp)+EPS)) &&
			fmax < rtol * SCALE)
			{ 	/* Normal return */
			cfree(gg);
			cfree(hh);
			cfree(ptmp);
			cfree(pcom);
			cfree(fcom);
			return(it_com);
			}
		fp = getf(pos, force);
		g = dg = 0.;
		if	(iter <= it_com)
			for	(i = 0; i < natoms; i++)
				{
				hh[i].fx = gg[i].fx = force[i].fx;
				hh[i].fy = gg[i].fy = force[i].fy;
				hh[i].fz = gg[i].fz = force[i].fz;
				}
		else
			{
			for	(i = 0; i < natoms; i ++)
				/*
				 * Polak-Ribiere variant
				 */
				{
				g +=	gg[i].fx * gg[i].fx +
					gg[i].fy * gg[i].fy +
					gg[i].fz * gg[i].fz;
				dg +=	(force[i].fx + gg[i].fx) * force[i].fx +
					(force[i].fy + gg[i].fy) * force[i].fy +
					(force[i].fz + gg[i].fz) * force[i].fz;
				}
		/*
		 * test for unlikely condition that gradient is already 0
		 */
			if	(g == 0.)
				{
				cfree(gg);
				cfree(hh);
				cfree(ptmp);
				cfree(pcom);
				cfree(fcom);
				return(it_com);
				}
			gam = dg/g;
			for	(i = 0; i < natoms; i++)
				{
				gg[i].fx = force[i].fx;
				gg[i].fy = force[i].fy;
				gg[i].fz = force[i].fz;
				/*
				 *	hh[i] = force[i]
				 */
				force[i].fx= hh[i].fx = gg[i].fx + gam*hh[i].fx;
				force[i].fy= hh[i].fy = gg[i].fy + gam*hh[i].fy;
				force[i].fz= hh[i].fz = gg[i].fz + gam*hh[i].fz;
				}
			}

/*
 *	Print status information on standard error:
 */
		wstat_min(stderr);

		} while(--it_com > 0);
	cfree(gg);
	cfree(hh);
	cfree(ptmp);
	cfree(pcom);
	cfree(fcom);
	return(0);
	}

/*
 *	minimizes along the direction force and replaces pcom by the actual
 *	vector displacement.
 */
double
linmin(pos, force, ftol)
tripd	*pos;
tripd	*force;
double	ftol;
	{
	int	j;
	double	f1dim(), dbrent(), fret, xmin;
	double	ax, bx, cx, fa, fb, fc;
	for	(j = 0, bx = 0.; j < natoms; j++)
		{
		/*
		 * xi = force, p = pos
		 */
		pcom[j].fx = pos[j].fx;
		pcom[j].fy = pos[j].fy;
		pcom[j].fz = pos[j].fz;
		fcom[j].fx = force[j].fx;
		fcom[j].fy = force[j].fy;
		fcom[j].fz = force[j].fz;
		bx = max(bx, abs(fcom[j].fx));
		bx = max(bx, abs(fcom[j].fy));
		bx = max(bx, abs(fcom[j].fz));
		}
	/*
	 * initial guess for brackets
	 */
	bx = min(1., (.05/bx));	/* force first step to be less than .05A */
	ax = -bx;
	mnbrak(&ax, &bx, &cx, &fa, &fb, &fc, f1dim);
	fret = dbrent(ax, bx, cx, f1dim, ftol, &xmin);
	/*
	 * Construct the vector results to return
	 */
	for	(j = 0; j < natoms; j++)
		{
		pos[j].fx += (force[j].fx = xmin*force[j].fx);
		pos[j].fy += (force[j].fy = xmin*force[j].fy);
		pos[j].fz += (force[j].fz = xmin*force[j].fz);
		}
	minbox(pos);
	return(fret);
	}



/*
 * given a function "func" and distinct initial points "ax" and "bx",
 * this routine searches in a downhill direction (defined by the function
 * as evaluated at the initial points) and returns new points "ax", "bx",
 * and "cx" which bracket a minimum of the function.  Also returned
 * are the function values at the three point "fa", "fb", and "fc".
 */
mnbrak(ax, bx, cx, fa, fb, fc, func)
double *ax, *bx, *cx, *fa, *fb, *fc, (*func)();
	{
	double	dum, r, q, u, ulim, fu;
	*fa = (*func)(*ax);
	*fb = (*func)(*bx);
	/*
	 * Switch roles of a and b so that we can go downhill in the
	 * direction from a to b
	 */
	if	(*fb > *fa)
		{
		dum = *ax;
		*ax = *bx;
		*bx = dum;
		dum = *fa;
		*fa = *fb;
		*fb = dum;
		}
	/*
	 * first guess for cx
	 */
	*cx = *bx + GOLD * (*bx - *ax); /*1st guess for c*/
	*fc = (*func)(*cx);
	/*
	 * loop will run until minimum is bracketed
	 */
	while	(*fb >= *fc)
		{
		r = (*bx - *ax)*(*fb - *fc);
		q = (*bx - *cx)*(*fb - *fa);
		/*
		 *Compute u by parabolic extrapolation from ax, bx, cx.
		 *TINY is used to prevent division by 0
		 */
		u = *bx - ((*bx - *cx)*q - (*bx - *ax)*r)/(2.*
			sign(max(abs(q - r), TINY), q - r));
		ulim = *bx + GLIMIT*(*cx - *bx);
		/*
		 * Parabolic u is between bx and cx -- try it
		 */
		if	((*bx - u)*(u - *cx) > 0.)
			{
			fu = (*func)(u);
			/*
			 * minimum between bx and cx
			 */
			if	(fu < *fc)
				{
				*ax = *bx;
				*fa = *fb;
				*bx = u;
				*fb = fu;
				continue;
				}
			/*
			 * minimum between ax and u
			 */
			else if (fu >*fb)
				{
				*cx = u;
				*fc = fu;
				continue;
				}
			/*
			 * parabolic not useful -- use default magnification
			 */
			u = *cx + GOLD*(*cx - *bx);
			fu = (*func)(u);
			}
		/*
		 * Parabolic fit is between cx and its allowed limit
		 */
		else if	((*cx - u)*(u - ulim) > 0.)
			{
			fu = (*func)(u);
			if	(fu < *fc)
				{
				*bx = *cx;
				*cx = u;
				u = *cx + GOLD*(*cx - *bx);
				*fb = *fc;
				*fc = fu;
				fu = (*func)(u);
				}
			}	
		/*
		 * limit parabolic u to maximum allowed value
		 */
		else if ((u - ulim)*(ulim - *cx) >= 0.)
			{
			u = ulim;
			fu = (*func)(u);
			}
		/*
		 * reject parabolic u -- use default magnification
		 */
		else
			{
			u = *cx + GOLD*(*cx - *bx);
			fu = (*func)(u);
			}
		/*
		 * eliminate oldest point and continue
		 */
		*ax = *bx;
		*bx = *cx;
		*cx = u;
		*fa = *fb;
		*fb = *fc;
		*fc = fu;
		fprintf(stderr,"m");
		}
	return;
	}



/*
 * Interface between single dimensional line minimizer and the
 * multidimensional minimization
 */
double
f1dim(x)
double	x;
	{
	/*ncom = natoms*/
	int	j;
	double	fp, getf();
	for	(j = 0; j < natoms; j++)
		{
		ptmp[j].fx = pcom[j].fx + x*fcom[j].fx;
		ptmp[j].fy = pcom[j].fy + x*fcom[j].fy;
		ptmp[j].fz = pcom[j].fz + x*fcom[j].fz;
		}
	minbox(ptmp);
	fp = getf(ptmp, ftmp);
	df_com = 0.;
	for	(j = 0; j < natoms; j++)
		{
		df_com -= force[j].fx*fcom[j].fx;
		df_com -= force[j].fy*fcom[j].fy;
		df_com -= force[j].fz*fcom[j].fz;
		}
	return(fp);
	}




/*
 * given a function f and a bracketing triplet of abscissas, this
 * routine isolates the minimum to a fractional precision of about
 * tol using a modification of Brent's method (that uses derivatives).
 * The abscissa of the minimum s returned as "xmin" and the minimum 
 * function value is returned as "dbrent".
 */
double
dbrent(ax, bx, cx, f, tol, xmin)
double	ax, bx, cx;
double	(*f)();
double	tol;
double	*xmin;
	{
	double	a, b, x, w, v, e, fw, fv, fx, dx, dv, dw, iter, xm, d2, d1;
	double	u1, u2, olde, tol1, tol2, d, u, du, fu;
	int	ok1;
	a = min(ax, cx);
	b = max(ax, cx);
	x = w = v = bx;
	e = 0.;
	fw = fv = fx = (*f)(x);
	dw = dv = dx = df_com;
	do	{
		xm = 0.5*(a + b);
		tol1 = tol*abs(x) + EPS;
		tol2 = 2.*tol1;
		if	(abs(x - xm) <= (tol2 - .5*(b - a)))
			{
			*xmin = x;
			return(fx);
			}
		if	(abs(e) > tol1)
			{
			d2 = d1 = 2.*(b - a); /*initialize to an out of bracket
						*value
						*/
			if	(dw != dx)	/*secant method*/
				d1 = (w - x)*dx/(dx - dw);
			if	(dv != dx)	/*secant method with 
						 *stored point
						 */
				d2 = (v - x)*dx/(dx - dv);
			u1 = x + d1;
			u2 = x + d2;
			ok1 = (((a - u1)*(u1 - b) > 0.)&&(dx*d1 <= 0.))?1:0;
			ok1 |= (((a - u2)*(u2 - b) > 0.)&&(dx*d2 <= 0.))?2:0;
			olde = e;
			e = d;
			/*
			 *take smallest value of d
			 */
			switch	(ok1)
				{
			case 1:
				d = d1;
				break;
			case 2:
				d = d2;
				break;
			case 3:
				if	(abs(d1) < abs(d2))
					d = d1;
				else
					d = d2;
				break;
				}
			if	(abs(d) > 0.5*abs(olde))
				{
				/*
				 * decide which segment by sign of
				 * the derivative
				 */
				if	(dx >= 0.)
					e = a - x;
				else
					e = b - x;
				d = 0.5*e; /*bisection*/
				}
			else
				{
				u = x + d;
				if	((u - a) < tol2||(b - u) < tol2)
					d = sign(tol1,xm - x);
				}
			}
		else
			{
			if	(dx >= 0.)
				e = a - x;
			else
				e = b - x;
			d = 0.5*e;	/*bisection*/
			}
		if	(abs(d) >= tol1)
			{
			u = x + d;
			fu = (*f)(u);
			}
		else
			{
			u = x + sign(tol1,d);
			fu = (*f)(u);
			/*
			 * if the minimum step in the downhill direction
			 * goes uphill -- done!
			 */
			if	(fu > fx)
				{
				*xmin = x;
				return(fx);
				}
			}
		du = df_com;
		if	(fu <= fx)
			{
			if	(u >= x)
				a = x;
			else
				b = x;
			/*housekeeping*/
			v = w;
			fv = fw;
			dv = dw;
			w = x;
			fw = fx;
			dw = dx;
			x = u;
			fx = fu;
			dx = du;
			}
		else
			{
			if	(u < x)
				a = u;
			else
				b = u;
			if	(fu <= fw || w == x)
				{
				v = w;
				fv = fw;
				dv = dw;
				w = u;
				fw = fu;
				dw = du;
				}
			else if	(fu <= fv || v == x || v == w)
				{
				v = u;
				fv = fu;
				dv = du;
				}
			}
		fprintf(stderr,"b");
		} while(--it_com > 1);
	*xmin = x;
	return(fx);
	}


/*
 * imaging system to keep molecules within the proper box
 */
minbox(p)
tripd	*p;

	{
	tripd	image;
	int	i, j;

	for	(i = 0; i < natoms; i++)
		{
		if	(atom[i].flags & A_MAJOR)
			{
			if	(atom[i].param1)
				{
				image.fx = -p[i].fx;
				image.fy = -p[i].fy;
				image.fz = -p[i].fz;
				mvimage(&p[i]);
				image.fx += p[i].fx;
				image.fy += p[i].fy;
				image.fz += p[i].fz;
				for	(j = i+1; j <= i+atom[i].param1; j++)
					{
					p[j].fx += image.fx;
					p[j].fy += image.fy;
					p[j].fz += image.fz;
					}
				}
			else
				mvimage(&p[i]);
			}
		}
	}





double
getf(pos, force)
tripd	*pos, *force;

	{
	int	i;
	getforce(pos, force);
	for(i = 0; i < natoms; i++)
		{
		if	(atom[i].flags & A_FIXED)
			force[i].fx = force[i].fy = force[i].fz = 0.;
		}
	return(V);
	}
