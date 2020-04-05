/*
 *	Steepest Descent Algorithm.
 */
#include	<md.h>

#define	EPS	1e-20
#define	FMAX	.01
#define	SCALE	1.
#define	LARGE	1e5

int	it_com;

minimize_sd(ftol, iter)	/*iter = number of iterations*/
double	ftol;
int	iter;
{
register double	fp, fmax, rtol, fmax1, fret;
double	get_f(), sqrt();
register int	i, jj;

it_com = iter;
rtol = sqrt(ftol);
ftol *= 0.5;
fp = LARGE;
/*
 * loop over iterations
 */
do	{
	fret = get_f(pos, force);
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
	fprintf(stderr, "I = %d F = %.13g (%d)\n", it_com, fmax, jj);
	if	((abs(fret - fp) <= ftol*(abs(fret)+abs(fp)+EPS)) &&
		fmax < rtol * SCALE) 	/* Normal return */
		return(it_com);
	if	(fmax > h * FMAX)
		{
		fmax = h * FMAX / fmax;
		for	(i = 0; i < natoms; i++)
			{
			pos[i].fx += force[i].fx * fmax;
			pos[i].fy += force[i].fy * fmax;
			pos[i].fz += force[i].fz * fmax;
			}
		}
	else
		for	(i = 0; i < natoms; i++)
			{
			pos[i].fx += force[i].fx;
			pos[i].fy += force[i].fy;
			pos[i].fz += force[i].fz;
			}
	min_box(pos);
	fp = fret;

/*
 *	Print status information on standard error:
 */
		wstat_min(stderr);

	} while(--it_com > 0);
return(0);
}

/*
 * imaging system to keep molecules within the proper box
 */
min_box(p)
tripd	*p;
{
tripd	image;
register int	i, j;

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
get_f(pos, force)
tripd	*pos, *force;
{
register int	i;
getforce(pos, force);
for(i = 0; i < natoms; i++)
	{
	if	(atom[i].flags & A_FIXED)
		force[i].fx = force[i].fy = force[i].fz = 0.;
	}
return(V);
}
