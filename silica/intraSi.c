/* this replaces all the lines in SiForce.c from line 129*/
also the call to intraSi sould be changes from intraSi(&pos[k],&force[k]); to just intraSi(k);*/
intraSi(k)
int k;
{
double Si_str, Si_bend, calcstr(), calcbend();
int i;
Si_str = 0;
Si_bend = 0;
/*Oh-H*/
Si_str += calcstr(k,k+2,0);
/*Si-Oh*/
Si_str += calcstr(k,k+1,1);
/*Si-Ob*/
for (i=0;i<3;i++)
	Si_str += calcstr(k+1,k+(4+i)*nSi,2);
/*H-Oh-Si*/
Si_bend += calcbend(k+2,k,k+1,0);
/*Oh-Si-Ob*/
for (i=0;i<3;i++)
	Si_bend += calcbend(k,k+1,k+(4+i)*nSi,1);
INTRAV_S = Si_str+Si_bend;
}
double calcstr(i,j,k)
int i,j,k;
{
tripd 	d; 
double	bond, dedr, pe;

d.fx = pos[i].fx - pos[j].fx;
d.fy = pos[i].fy - pos[j].fy;
d.fz = pos[i].fz - pos[j].fz;
bond = sqrt(d.fx*d.fx+d.fy*d.fy+d.fz*d.fz);
dedr = kStr[k]*(bond-eqBond[k])/bond;
pe = 0.5*kStr[k]*sq(bond-eqBond[k]);
force[i].fx -= (d.fx = dedr*d.fx);
force[i].fy -= (d.fy = dedr*d.fy);
force[i].fz -= (d.fz = dedr*d.fz);
force[j].fx += d.fx;
force[j].fy += d.fy;
force[j].fz += d.fz;
return(pe);

}
double calcbend(i,j,k,l)
int i,j,k,l;
{
int i0,i1,i2;
tripd r[2], grad[3];
double pe,r1,r2,r1r2,n1n2,da,dedd;

i0 = j;/*center atom*/
i1 = i;
i2 = k;

r[0].fx = pos[i0].fx - pos[i1].fx;
r[0].fy = pos[i0].fy - pos[i1].fy;
r[0].fz = pos[i0].fz - pos[i1].fz;
r[1].fx = pos[i0].fx - pos[i2].fx;
r[1].fy = pos[i0].fy - pos[i2].fy;
r[1].fz = pos[i0].fz - pos[i2].fz;

r1 = sqrt(r[0].fx*r[0].fx + r[0].fy*r[0].fy + r[0].fz*r[0].fz);
r2 = sqrt(r[1].fx*r[1].fx + r[1].fy*r[1].fy + r[1].fz*r[1].fz);
r1r2 =   (r[0].fx*r[1].fx + r[0].fy*r[1].fy + r[0].fz*r[1].fz);
n1n2 = r1r2 / (r1*r2);
da = acos(n1n2) - eqBend[l];
pe = 0.5*kBend[l]*da*da;
dedd = kBend[l]*da/ (sqrt(1.-n1n2*n1n2)*r1*r2);
force[i1].fx -= (grad[1].fx = dedd*(r[1].fx - r1r2*r[0].fx/(r1*r1)));
force[i1].fy -= (grad[1].fy = dedd*(r[1].fy - r1r2*r[0].fy/(r1*r1)));
force[i1].fz -= (grad[1].fz = dedd*(r[1].fz - r1r2*r[0].fz/(r1*r1)));
force[i2].fx -= (grad[2].fx = dedd*(r[0].fx - r1r2*r[1].fx/(r2*r2)));
force[i2].fy -= (grad[2].fy = dedd*(r[0].fy - r1r2*r[1].fy/(r2*r2)));
force[i2].fz -= (grad[2].fz = dedd*(r[0].fz - r1r2*r[1].fz/(r2*r2)));
force[i0].fx += grad[1].fx + grad[2].fx;
force[i0].fy += grad[1].fy + grad[2].fy;
force[i0].fz += grad[1].fz + grad[2].fz;
    return(pe);
}
