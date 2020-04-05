#include	<md.h>

#define NBUFS 2

/*
 * The state of the system is defined by the following parameters:
 *	Standard Variables:
 *		h	time increment
 *		h1	previous integration step
 *		h2	previous elapsed time step
 *		etime	elapsed real time
 *		pos	positions
 *		vel	(pos(t) - pos(t-h))/h
 *		"flags"	atom flags
 *
 *	Constant Temperature Algorithm Variables:
 *		g	number of degrees of freedom
 *		Teq	equilibrium temperature
 *		s	scaling variable
 *		vs	scaling variable's equivalent to "vel"
 *		
 */
static int	sr_init[NBUFS] = {0};
static float	sr_h[NBUFS];
static double	sr_h1[NBUFS];
static double	sr_h2[NBUFS];
static float	sr_etime[NBUFS];
static tripd *	sr_pos[NBUFS];
static tripd *	sr_vel[NBUFS];
static int *	sr_flags[NBUFS];
static int	sr_g[NBUFS];
static double	sr_Teq[NBUFS];
static double	sr_s[NBUFS];
static double	sr_vs[NBUFS];

save(buf)
int	buf;
	{
	int	i;

	if	(buf < 1 || buf > NBUFS)
		ERROR((stderr, "save:  illegal buffer %d\n", buf),return);

	if	(!sr_init[--buf])	/* we need to allocate space */
		{
		if	((sr_pos[buf] = (tripd *)malloc(natoms * sizeof(tripd)))
				== NULL ||
			(sr_vel[buf] = (tripd *)malloc(natoms * sizeof(tripd)))
				== NULL ||
			(sr_flags[buf] = (int *)malloc(natoms * sizeof(int)))
				== NULL)
			ERROR((stderr, "save:  out of core\n"), return);
		sr_init[buf] = 1;
		}
	sr_h[buf] = h;
	sr_h1[buf] = h1;
	sr_h2[buf] = h2;
	sr_etime[buf] = etime;
	sr_g[buf] = g;
	sr_Teq[buf] = Teq;
	sr_s[buf] = s;
	sr_vs[buf] = vs;
	for	(i = 0; i < natoms; i++)
		{
		sr_pos[buf][i].fx = pos[i].fx;
		sr_pos[buf][i].fy = pos[i].fy;
		sr_pos[buf][i].fz = pos[i].fz;
		sr_vel[buf][i].fx = vel[i].fx;
		sr_vel[buf][i].fy = vel[i].fy;
		sr_vel[buf][i].fz = vel[i].fz;
		sr_flags[buf][i] = atom[i].flags;
		}
	return(0);
	}

restore(buf)
int	buf;
	{
	int	i;

	if	(buf < 1 || buf > NBUFS)
		ERROR((stderr, "restore:  illegal buffer %d\n", buf),return);

	if	(!sr_init[--buf])	/* nothing is stored yet */
		ERROR((stderr,"restore:  nothing saved in buffer %d\n",++buf),
			return);
	h = sr_h[buf];
	h1 = sr_h1[buf];
	h2 = sr_h2[buf];
	etime = sr_etime[buf];
	g = sr_g[buf];
	Teq = sr_Teq[buf];
	s = sr_s[buf];
	vs = sr_vs[buf];
	for	(i = 0; i < natoms; i++)
		{
		pos[i].fx = sr_pos[buf][i].fx;
		pos[i].fy = sr_pos[buf][i].fy;
		pos[i].fz = sr_pos[buf][i].fz;
		vel[i].fx = sr_vel[buf][i].fx;
		vel[i].fy = sr_vel[buf][i].fy;
		vel[i].fz = sr_vel[buf][i].fz;
		atom[i].flags = sr_flags[buf][i];
		}
	return(0);
	}
