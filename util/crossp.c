#include	<md.h>

crossp(a,b,c)
tripdouble	*a,*b,*c;
	{
	c->fx = a->fy*b->fz - a->fz*b->fy;
	c->fy = a->fz*b->fx - a->fx*b->fz;
	c->fz = a->fx*b->fy - a->fy*b->fx;
	}
