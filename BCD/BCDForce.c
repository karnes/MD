#include	<md.h>
#include	<system.h>
#include	<math.h>
#include	<water.h>
//#define EXCON 0.1

/* NOTE: currently verbose. replace i * 21 and ((i + 1)%7)*21 with k,l later */

BCDForce()
{
int	i,j,k,nw,nsolvent, index;
tripd pri,sec,sdl,solc;
double r,totm,dz;
double thx, thy, thz;
double xix, xiy, xiz;
double zex, zey, zez;
bondE = bendE = torsE = nb14E  = nb15E = 0.;

nsolvent = natoms - nBCD*BCDs - nsolute;
nw = nsolvent - BrOs*nBrO;

/* non-PBC algorithm */
/*
BCDcom.fx = BCDcom.fy = BCDcom.fz = totm = 0.0;
for(i=0;i<BCDs;i++){
	j=nsolvent+i;
	BCDcom.fx += mass[j]*pos[j].fx;
	BCDcom.fy += mass[j]*pos[j].fy;
	BCDcom.fz += mass[j]*pos[j].fz;
	totm += mass[j];
}
BCDcom.fx /= totm;
BCDcom.fy /= totm;
BCDcom.fz /= totm;
//fprintf(stderr,"std CoM: x = %5.3f, y = %5.3f, z = %5.3f\n",BCDcom.fx,BCDcom.fy,BCDcom.fz);
*/
/* calculate BCD center of mass */
/* Bai-Breen algorithm */

BCDcom.fx = BCDcom.fy = BCDcom.fz = totm = 0.0;
xix = xiy = xiz = 0.0;
zex = zey = zez = 0.0;
for(i=0;i<BCDs;i++){
	j = nsolvent+i;
	thx = (pos[j].fx + xwall)*PI/xwall;
	xix += mass[j]*cos(thx);
	zex += mass[j]*sin(thx);
	thy = (pos[j].fy + ywall)*PI/ywall;
	xiy += mass[j]*cos(thy);
	zey += mass[j]*sin(thy);
	thz = (pos[j].fz + zwall)*PI/zwall;
	xiz += mass[j]*cos(thz);
	zez += mass[j]*sin(thz);
	totm += mass[j];
}
thx = atan2(-(zex/totm),-(xix/totm)) + PI;
BCDcom.fx = (xwall * thx / PI) - xwall;
thy = atan2(-(zey/totm),-(xiy/totm)) + PI;
BCDcom.fy = (ywall * thy / PI) - ywall;
thz = atan2(-(zez/totm),-(xiz/totm)) + PI;
BCDcom.fz = (zwall * thz / PI) - zwall;
//fprintf(stderr,"std CoM: x = %5.3f, y = %5.3f, z = %5.3f\n",BCDcom.fx,BCDcom.fy,BCDcom.fz);

/* confine BCD to potential window */
if(windowOn==1){
   window();
}
else if(windowOn!=0){
   fprintf(stderr,"BCDForce.c: error-- windowOn should be 1 or 0. value = %d\n",windowOn);
   exit(0);
}

pri.fx = pri.fy = pri.fz = totm = 0.0;
xix = xiy = xiz = 0.0;
zex = zey = zez = 0.0;
BCDz.fx = BCDz.fy = BCDz.fz = 0.0;
// get BCD molecular vector 
for(i=0;i<7;i++){
	j=nsolvent+(i*21);
//9 and 16 are C(attached to -CH3-O-H) and ether O of glucose ring
	thx = (pos[j+9].fx + xwall)*PI/xwall;
	xix += mass[j+9]*cos(thx);
	zex += mass[j+9]*sin(thx);
	thx = (pos[j+16].fx + xwall)*PI/xwall;
	xix += mass[j+16]*cos(thx);
	zex += mass[j+16]*sin(thx);
	thy = (pos[j+9].fy + ywall)*PI/ywall;
	xiy += mass[j+9]*cos(thy);
	zey += mass[j+9]*sin(thy);
	thy = (pos[j+16].fy + ywall)*PI/ywall;
	xiy += mass[j+16]*cos(thy);
	zey += mass[j+16]*sin(thy);
	thz = (pos[j+9].fz + zwall)*PI/zwall;
	xiz += mass[j+9]*cos(thz);
	zez += mass[j+9]*sin(thz);
	thz = (pos[j+16].fz + zwall)*PI/zwall;
	xiz += mass[j+16]*cos(thz);
	zez += mass[j+16]*sin(thz);
	totm += mass[j+9] + mass[j+16];
}
Ms = totm;
thx = atan2(-(zex/totm),-(xix/totm)) + PI;
btop.fx = pri.fx = (xwall * thx / PI) - xwall;
thy = atan2(-(zey/totm),-(xiy/totm)) + PI;
btop.fy = pri.fy = (ywall * thy / PI) - ywall;
thz = atan2(-(zez/totm),-(xiz/totm)) + PI;
btop.fz = pri.fz = (zwall * thz / PI) - zwall;

sec.fx = sec.fy = sec.fz = totm = 0.0;
xix = xiy = xiz = 0.0;
zex = zey = zez = 0.0;
for(i=0;i<7;i++){
	j=nsolvent+(i*21);
// 3 and 5 are carbons on glucose ring connected to OH's
	thx = (pos[j+3].fx + xwall)*PI/xwall;
	xix += mass[j+3]*cos(thx);
	zex += mass[j+3]*sin(thx);
	thx = (pos[j+5].fx + xwall)*PI/xwall;
	xix += mass[j+5]*cos(thx);
	zex += mass[j+5]*sin(thx);
	thy = (pos[j+3].fy + ywall)*PI/ywall;
	xiy += mass[j+3]*cos(thy);
	zey += mass[j+3]*sin(thy);
	thy = (pos[j+5].fy + ywall)*PI/ywall;
	xiy += mass[j+5]*cos(thy);
	zey += mass[j+5]*sin(thy);
	thz = (pos[j+3].fz + zwall)*PI/zwall;
	xiz += mass[j+3]*cos(thz);
	zez += mass[j+3]*sin(thz);
	thz = (pos[j+5].fz + zwall)*PI/zwall;
	xiz += mass[j+5]*cos(thz);
	zez += mass[j+5]*sin(thz);
	totm += mass[j+3] + mass[j+5];
}
Ml = totm;
thx = atan2(-(zex/totm),-(xix/totm)) + PI;
bbot.fx = sec.fx = (xwall * thx / PI) - xwall;
thy = atan2(-(zey/totm),-(xiy/totm)) + PI;
bbot.fy = sec.fy = (ywall * thy / PI) - ywall;
thz = atan2(-(zez/totm),-(xiz/totm)) + PI;
bbot.fz = sec.fz = (zwall * thz / PI) - zwall;

/*
if(tc==0){
// get pore radius 
poreRad = 0.0;
for(i=0;i<7;i++){
   k = nsolvent + i*21;
   sdl.fx = pos[k+9].fx - pri.fx;
   sdl.fy = pos[k+9].fy - pri.fy;
   sdl.fz = pos[k+9].fz - pri.fz;
   mvimage(&sdl);
   poreRad += sqrt(sdl.fx*sdl.fx + sdl.fy*sdl.fy + sdl.fz*sdl.fz);
   sdl.fx = pos[k+16].fx - pri.fx;
   sdl.fy = pos[k+16].fy - pri.fy;
   sdl.fz = pos[k+16].fz - pri.fz;
   mvimage(&sdl);
   poreRad += sqrt(sdl.fx*sdl.fx + sdl.fy*sdl.fy + sdl.fz*sdl.fz);
   sdl.fx = pos[k+3].fx - sec.fx;
   sdl.fy = pos[k+3].fy - sec.fy;
   sdl.fz = pos[k+3].fz - sec.fz;
   mvimage(&sdl);
   poreRad += sqrt(sdl.fx*sdl.fx + sdl.fy*sdl.fy + sdl.fz*sdl.fz);
   sdl.fx = pos[k+5].fx - sec.fx;
   sdl.fy = pos[k+5].fy - sec.fy;
   sdl.fz = pos[k+5].fz - sec.fz;
   mvimage(&sdl);
   poreRad += sqrt(sdl.fx*sdl.fx + sdl.fy*sdl.fy + sdl.fz*sdl.fz);
}

poreRad /= (double)(i*4);
fprintf(stderr,"poreRad = %f\n",poreRad);
}*/
//fprintf(stderr,"pri CoM: x = %5.3f, y = %5.3f, z = %5.3f\n",pri.fx,pri.fy,pri.fz);
//fprintf(stderr,"sec CoM: x = %5.3f, y = %5.3f, z = %5.3f\n",sec.fx,sec.fy,pri.fz);

/* non-PBC algorithm */
/*
pri.fx = pri.fy = pri.fz = totm = 0.0;
for(i=0;i<7;i++){
	j=nsolvent+(i*21);
	pri.fx += mass[j+9]*pos[j+9].fx + mass[j+16]*pos[j+16].fx;
	pri.fy += mass[j+9]*pos[j+9].fy + mass[j+16]*pos[j+16].fy;
	pri.fz += mass[j+9]*pos[j+9].fz + mass[j+16]*pos[j+16].fz;
	totm += mass[j+9] + mass[j+16];
}
pri.fx /= totm;
pri.fy /= totm;
pri.fz /= totm;

sec.fx = sec.fy = sec.fz = totm = 0.0;
for(i=0;i<7;i++){
	j=nsolvent+(i*21);
	sec.fx += mass[j+5]*pos[j+5].fx + mass[j+3]*pos[j+3].fx;
	sec.fy += mass[j+5]*pos[j+5].fy + mass[j+3]*pos[j+3].fy;
	sec.fz += mass[j+5]*pos[j+5].fz + mass[j+3]*pos[j+3].fz;
	totm += mass[j+5] + mass[j+3];
}
sec.fx /= totm;
sec.fy /= totm;
sec.fz /= totm;

//fprintf(stderr,"pri CoM: x = %5.3f, y = %5.3f, z = %5.3f\n",pri.fx,pri.fy,pri.fz);
//fprintf(stderr,"sec CoM: x = %5.3f, y = %5.3f, z = %5.3f\n",sec.fx,sec.fy,pri.fz);
*/
BCDz.fx = pri.fx - sec.fx;
BCDz.fy = pri.fy - sec.fy;
BCDz.fz = pri.fz - sec.fz;

r = sqrt(BCDz.fx*BCDz.fx + BCDz.fy*BCDz.fy + BCDz.fz*BCDz.fz);
capR = r;
BCDz.fx /= r;
BCDz.fy /= r;
BCDz.fz /= r;
//fprintf(stderr,"length BCDz = %f\n",BCDz.fx*BCDz.fx+BCDz.fy*BCDz.fy+BCDz.fz*BCDz.fz);

BCDx.fx = BCDx.fy = BCDx.fz = totm = 0.0;
for(i=0;i<21;i++){
   j = i + nsolvent;
   BCDx.fx += pos[j].fx*mass[j];
   BCDx.fy += pos[j].fy*mass[j];
   BCDx.fz += pos[j].fz*mass[j];
   totm+=mass[j];
}
BCDx.fx /= totm;
BCDx.fy /= totm;
BCDx.fz /= totm;
BCDx.fx-=BCDcom.fx;
BCDx.fy-=BCDcom.fy;
BCDx.fz-=BCDcom.fz;
// make orthogonal to BCDz
// dz = deviation from BCDcom on BCDz
// (project onto BCDz)
dz = BCDx.fx*BCDz.fx + BCDx.fy*BCDz.fy + BCDx.fz*BCDz.fz;
//fprintf(stderr,"dz = %f\n",dz);
BCDx.fx-=dz*BCDz.fx;
BCDx.fy-=dz*BCDz.fy;
BCDx.fz-=dz*BCDz.fz;  
r = sqrt(sq(BCDx.fx)+sq(BCDx.fy)+sq(BCDx.fz));
BCDx.fx /= r;
BCDx.fy /= r;
BCDx.fz /= r;
// orthogonal?
//dz = BCDx.fx*BCDz.fx + BCDx.fy*BCDz.fy + BCDx.fz*BCDz.fz;
//fprintf(stderr,"cos Z*X = %f\n",dz);
//fprintf(stderr,"length BCDx = %f\n",BCDx.fx*BCDx.fx + BCDx.fy*BCDx.fy + BCDx.fz*BCDx.fz);

BCDy.fx = BCDx.fy*BCDz.fz - BCDx.fz*BCDz.fy;
BCDy.fy = BCDx.fz*BCDz.fx - BCDx.fx*BCDz.fz;
BCDy.fz = BCDx.fx*BCDz.fy - BCDx.fy*BCDz.fx;
r = sqrt(sq(BCDy.fx)+sq(BCDy.fy)+sq(BCDy.fz));
BCDy.fx /= r;
BCDy.fy /= r;
BCDy.fz /= r;
//fprintf(stderr,"length BCDy = %f\n",BCDy.fx*BCDy.fx+BCDy.fy*BCDy.fy+BCDy.fz*BCDy.fz);

INTRABCD = 0.0;
//fprintf(stderr,"enter BCDForce.c nsolvent=%d\n",nsolvent);
intraBCD();
//if(tc%dataRatex==0)
//fprintf(stderr,"tor: %d, bond:: %d, bend: %d\n",torsions, bonds, bends);
//fprintf(stderr,"leaving BCDForce.c\n");
}
/*
*/
intraBCD(){
int i,j,/*k,l,*/m,n,n1,n2,n3,nw,nsolvent;
double dx,dy,dz,r2,eg,dedr;
double calcstrch(), calcbend(), calctorq(), ttraljq14(), ttraljq();  
double bond_e, tors_e, bend_e, lj14_e, ljqNB_e;
tripd sdist;
bond_e = 0.;
tors_e = 0.;
bend_e = 0.;
lj14_e = 0.;
ljqNB_e = 0.;

nsolvent = natoms - nBCD*BCDs - nsolute;
nw = nsolvent - nBrO*BrOs;
//fprintf(stderr,"enter intraBCD. nsolvent = %d\n",nsolvent);

/*fprintf(stderr,"%d ",k);*/
for(i=0;i<7;i++){   //loop over the 7 sub-units of the b-CD
//	k = 21*i+(natoms-BCDs); // after it is working, s/'ix21'/k/g -jjk
//	l = ((1+i)%7)*21+(natoms-BCDs); // also s/ ... /l/g -jjk
/* loop over all bond stretches*/
	/*carbon-carbon*/
	for (j=0;j<5;j++){
		n = i*21+(2*j+1);  // 1,3,5,7,9
		m = i*21+(2*j+3);  // 3,5,7,9,11
    		bond_e += calcstrch(n,m,0);
//fprintf(stderr,"%f ",bond_e*KCAL);
	}
	/*carbon-hydrogen 1*/
	for (j=0;j<5;j++){
		n = i*21 +(2*j+3);  // 3,5,7,9,11
		m = i*21+(2*j+4);  // 4,6,8,10,12
    		bond_e += calcstrch(n,m,1);
	}
	n = i*21+11;
	m = i*21+13;
    	bond_e += calcstrch(n,m,1);
	/*fprintf(stderr,"%f ",bond_e*KCAL);*/
	
	/*carbon-hydrogen 2*/
	n = i*21+1;
	m = i*21+2;
    	bond_e += calcstrch(n,m,2);
	/*fprintf(stderr,"%f ",bond_e*KCAL);*/
	/*carbon-oxygen 3*/
	n = i*21+3;
	m = i*21+19;
    	bond_e += calcstrch(n,m,3);
	n = i*21+5;
	m = i*21+17;
    	bond_e += calcstrch(n,m,3);
	n = i*21+11;
	m = i*21+14;
    	bond_e += calcstrch(n,m,3);
	/*fprintf(stderr,"%f ",bond_e*KCAL);*/
	/*carbon-oxygen 4*/
	n = i*21+1;
	m = i*21+16;
    	bond_e += calcstrch(n,m,4);
	n = i*21+7;
	m = i*21+0;
    	bond_e += calcstrch(n,m,4);
	n = i*21+9;
	m = i*21+16;
    	bond_e += calcstrch(n,m,4);
	/*inter-unit*/
	n = i*21+0;
	m = ((i+1)%7)*21+1;
    	bond_e += calcstrch(n,m,4);
	/*fprintf(stderr,"%f ",bond_e*KCAL);*/
	/*oxygen-hydrogen 5*/
	n = i*21+14;
	m = i*21+15;
    	bond_e += calcstrch(n,m,5);
	n = i*21+17;
	m = i*21+18;
    	bond_e += calcstrch(n,m,5);
	n = i*21+19;
	m = i*21+20;
    	bond_e += calcstrch(n,m,5);
//fprintf(stderr,"stretches: %f ",bond_e*KCAL);
/*fprintf(stderr,"%f ",bond_e*KCAL);*/
/*fprintf(stderr,"\n");*/
/* loop over all torsions*/
//	tors_e=0.;
	/*O-C-C-O 0*/
	n = i*21 + 14;/*O*/
	n1 = i*21 + 11;/*C*/
	n2 = i*21 + 9;/*C*/
	n3 = i*21 + 16;/*O*/
    	tors_e += calctorq(n,n1,n2,n3,0);
	n = i*21 + 16;/*O*/
	n1 = i*21 + 1;/*C*/
	n2 = i*21 + 3;/*C*/
	n3 = i*21 + 19;/*O*/
    	tors_e += calctorq(n,n1,n2,n3,0);
	n1 = i*21 + 9;/*C*/
	n2 = i*21 + 7;/*C*/
	n3 = i*21 + 0;/*O*/
    	tors_e += calctorq(n,n1,n2,n3,0);
	n = i*21 + 17;/*O*/
	n1 = i*21 + 5;/*C*/
	n2 = i*21 + 3;/*C*/
	n3 = i*21 + 19;/*O*/
    	tors_e += calctorq(n,n1,n2,n3,0);
	n = i*21 + 0;/*O*/
	n1 = i*21 + 7;/*C*/
	n2 = i*21 + 5;/*C*/
	n3 = i*21 + 17;/*O*/
    	tors_e += calctorq(n,n1,n2,n3,0);
	n = i*21 + 0;/*O 0*/
	n1 = ((i+1)%7)*21 + 1;/*C 22*/
	n2 = n1 + 2;/*C 24*/
	n3 = n1 + 18;/*O 40*/
    	tors_e += calctorq(n,n1,n2,n3,0);
	/*C-C-C-H 1*/
	for(j=0;j<4;j++){
		n = i*21 + 2*j + 1;/*C 1,3,5,7*/
		n1 = n + 2;/*C 3,5,7,9*/
		n2 = n + 4;/*C 5,7,9,11*/
		n3 = n + 5;/*O 6,8,10,12*/
		tors_e += calctorq(n,n1,n2,n3,1);
	}
	for(j=0;j<5;j++){
		n = i*21 + 2*j + 2;/*H 2,4,6,8,10*/
		n1 = n - 1;/*C 1,3,5,7,9*/
		n2 = n + 1;/*C 3,5,7,9,11*/
		n3 = n + 2;/*H 4,6,8,10,12*/
		tors_e += calctorq(n,n1,n2,n3,1);
		n3 = n + 3;/*C 5,7,9,11,13*/
		tors_e += calctorq(n,n1,n2,n3,1);
	}
	/*C-C-C-O*/
	n = i*21 + 1;/*C*/
	n1 = i*21 + 3;/*C*/
	n2 = i*21 + 5;/*C*/
	n3 = i*21 + 17;/*O*/
    	tors_e += calctorq(n,n1,n2,n3,1);
	/*H-C-C-O*/
	n = i*21 + 2;/*H*/
	n1 = i*21 + 1;/*C*/
	n2 = i*21 + 3;/*C*/
	n3 = i*21 + 19;/*O*/
    	tors_e += calctorq(n,n1,n2,n3,1);
	/*C-C-C-O*/
	n = i*21 + 3;/*C*/
	n1 = i*21 + 5;/*C*/
	n2 = i*21 + 7;/*C*/
	n3 = i*21 + 0;/*O*/
    	tors_e += calctorq(n,n1,n2,n3,1);
	/*C-C-C-O*/
	n = i*21 + 5;/*C*/
	n1 = i*21 + 7;/*C*/
	n2 = i*21 + 9;/*C*/
	n3 = i*21 + 16;/*O*/
    	tors_e += calctorq(n,n1,n2,n3,1);
	/*C-C-C-O*/
	n = i*21 + 7;/*C*/
	n1 = i*21 + 5;/*C*/
	n2 = i*21 + 3;/*C*/
	n3 = i*21 + 19;/*O*/
    	tors_e += calctorq(n,n1,n2,n3,1);
	/*C-C-C-H/O*/
	n = i*21 + 7;/*C*/
	n1 = i*21 + 9;/*C*/
	n2 = i*21 + 11;/*C*/
	n3 = i*21 + 13;/*H*/
    	tors_e += calctorq(n,n1,n2,n3,1);
	n3 = i*21 + 14;/*O*/
    	tors_e += calctorq(n,n1,n2,n3,1);
	/*C-C-C-O*/
	n = i*21 + 9;/*C*/
	n1 = i*21 + 7;/*C*/
	n2 = i*21 + 5;/*C*/
	n3 = i*21 + 17;/*O*/
    	tors_e += calctorq(n,n1,n2,n3,1);
	/*C-C-C-O*/
	n = i*21 + 11;/*C*/
	n1 = i*21 + 9;/*C*/
	n2 = i*21 + 7;/*C*/
	n3 = i*21 + 0;/*O*/
    	tors_e += calctorq(n,n1,n2,n3,1);
	/*O-C-C-O*/
	n = i*21 + 16;/*O*/
	n1 = i*21 + 1;/*C*/
	n2 = i*21 + 3;/*C*/
	n3 = i*21 + 5;/*C*/
    	tors_e += calctorq(n,n1,n2,n3,1);
	/*O-C-C-C*/
	n = i*21 + 0;/*O 0*/
	n1 = ((i+1)%7)*21 + 1;/*C 22*/
	n2 = n1 + 2;/*C 24*/ 
	n3 = n1 + 4;/*C 26*/
    	tors_e += calctorq(n,n1,n2,n3,1);
	/*H-C-O-H 2*/
	n = i*21 + 4;/*H*/
	n1 = i*21 + 3;/*C*/
	n2 = i*21 + 19;/*O*/
	n3 = i*21 + 20;/*H*/
    	tors_e += calctorq(n,n1,n2,n3,2);
	/*H-C-O-H 2*/
	n = i*21 + 6;/*H*/
	n1 = i*21 + 5;/*C*/
	n2 = i*21 + 17;/*O*/
	n3 = i*21 + 18;/*H*/
    	tors_e += calctorq(n,n1,n2,n3,2);
	/*H-C-O-H 2*/
	n = i*21 + 12;/*H*/
	n1 = i*21 + 11;/*C*/
	n2 = i*21 + 14;/*O*/
	n3 = i*21 + 15;/*H*/
    	tors_e += calctorq(n,n1,n2,n3,2);
	n = i*21 + 13;/*O*/
    	tors_e += calctorq(n,n1,n2,n3,2);
	/*H-C-C-O 3*/
	n = i*21 + 4;/*H*/
	n1 = i*21 + 3;/*C*/
	n2 = i*21 + 5;/*C*/
	n3 = i*21 + 17;/*O*/
    	tors_e += calctorq(n,n1,n2,n3,3);
	n = i*21 + 8;/*H*/
	n1 = i*21 + 7;/*C*/
    	tors_e += calctorq(n,n1,n2,n3,3);
	/*H-C-C-O 3*/
	n = i*21 + 6;/*H*/
	n1 = i*21 + 5;/*C*/
	n2 = i*21 + 7;/*C*/
	n3 = i*21 + 0;/*O*/
    	tors_e += calctorq(n,n1,n2,n3,3);
	n = i*21 + 10;/*H*/
	n1 = i*21 + 9;/*C*/
    	tors_e += calctorq(n,n1,n2,n3,3);
	/*H-C-C-O 3*/
	n = i*21 + 8;/*H*/
	n1 = i*21 + 7;/*C*/
	n2 = i*21 + 9;/*C*/
	n3 = i*21 + 16;/*O*/
    	tors_e += calctorq(n,n1,n2,n3,3);
	n = i*21 + 12;/*H*/
	n1 = i*21 + 11;/*C*/
    	tors_e += calctorq(n,n1,n2,n3,3);
	n = i*21 + 13;/*H*/
    	tors_e += calctorq(n,n1,n2,n3,3);
	/*H-C-C-O 3*/
	n = i*21 + 6;/*H*/
	n1 = i*21 + 5;/*C*/
	n2 = i*21 + 3;/*C*/
	n3 = i*21 + 19;/*O*/
    	tors_e += calctorq(n,n1,n2,n3,3);
	/*H-C-C-O 3*/
	n = i*21 + 10;/*H*/
	n1 = i*21 + 9;/*C*/
	n2 = i*21 + 11;/*C*/
	n3 = i*21 + 14;/*O*/
    	tors_e += calctorq(n,n1,n2,n3,3);
	/*H-C-C-O 3*/
	n = i*21 + 16;/*H*/
	n1 = i*21 + 1;/*C*/
	n2 = i*21 + 3;/*C*/
	n3 = i*21 + 4;/*O*/
    	tors_e += calctorq(n,n1,n2,n3,3);
	/*O-C-C-H 3*/
	n = i*21 + 0;/*O 0*/
	n1 = ((i+1)%7)*21 + 1;/*C 22*/
	n2 = n1 + 2;/*C 24*/
	n3 = n1 + 3;/*H 25*/
    	tors_e += calctorq(n,n1,n2,n3,3);
	/*C-O-C-H 4*/
	n = i*21 + 1;/*C*/
	n1 = i*21 + 16;/*O*/
	n2 = i*21 + 9;/*C*/
	n3 = i*21 + 10;/*H*/
    	tors_e += calctorq(n,n1,n2,n3,4);
	/*H-C-O-C 4*/
	n = i*21 + 2;/*H*/
	n1 = i*21 + 1;/*C*/
	n2 = i*21 + 16;/*O*/
	n3 = i*21 + 9;/*C*/
    	tors_e += calctorq(n,n1,n2,n3,4);
	/*C-O-C-H 4*/
	n = i*21 + 7;/*C*/
	n1 = i*21 + 0;/*O*/
	n2 = ((i+1)%7)*21 + 1;/*C*/
	n3 = n2 + 1;/*H*/
//	fprintf(stderr,"torque: n1=%d, n2=%d, n3=%d, n4=%d, index=%d\n",n,n1,n2,n3,4);
    	tors_e += calctorq(n,n1,n2,n3,4);
	/*H-C-O-C 4*/
	n = i*21 + 8;/*H*/
	n1 = i*21 + 7;/*C*/
	n2 = i*21 + 0;/*O*/
	n3 = ((i+1)%7)*21 + 1;/*C*/
    	tors_e += calctorq(n,n1,n2,n3,4);
	/*C-C-O-H 5*/
	n = i*21 + 1;/*C*/
	n1 = i*21 + 3;/*C*/
	n2 = i*21 + 19;/*O*/
	n3 = i*21 + 20;/*H*/
    	tors_e += calctorq(n,n1,n2,n3,5);
	n = i*21 + 5;/*C*/
    	tors_e += calctorq(n,n1,n2,n3,5);
	/*C-C-O-H 5*/
	n = i*21 + 3;/*C*/
	n1 = i*21 + 5;/*C*/
	n2 = i*21 + 17;/*O*/
	n3 = i*21 + 18;/*H*/
    	tors_e += calctorq(n,n1,n2,n3,5);
	n = i*21 + 7;/*C*/
    	tors_e += calctorq(n,n1,n2,n3,5);
	/*C-C-O-H 5*/
	n = i*21 + 9;/*C*/
	n1 = i*21 + 11;/*C*/
	n2 = i*21 + 14;/*O*/
	n3 = i*21 + 15;/*H*/
    	tors_e += calctorq(n,n1,n2,n3,5);
	/*C-C-O-C 6*/
	n = i*21 + 1;/*C*/
	n1 = i*21 + 16;/*C*/
	n2 = i*21 + 9;/*O*/
	n3 = i*21 + 7;/*C*/
    	tors_e += calctorq(n,n1,n2,n3,6);
	n3 = i*21 + 11;/*C*/
    	tors_e += calctorq(n,n1,n2,n3,6);
	/*C-C-O-C 6*/
	n = i*21 + 3;/*C*/
	n1 = i*21 + 1;/*C*/
	n2 = i*21 + 16;/*O*/
	n3 = i*21 + 9;/*C*/
    	tors_e += calctorq(n,n1,n2,n3,6);
	/*C-C-O-C 6*/
	n = i*21 + 5;/*C*/
	n1 = i*21 + 7;/*C*/
	n2 = i*21 + 0;/*O*/
	n3 = ((i+1)%7)*21 + 1;/*C*/
    	tors_e += calctorq(n,n1,n2,n3,6);
	n = i*21 + 9;/*C*/
    	tors_e += calctorq(n,n1,n2,n3,6);
	/*C-O-C-C 6*/
	n = i*21 + 7;/*C*/
	n1 = i*21 + 0;/*O*/
	n2 = ((i+1)%7)*21 + 1;/*C*/
	n3 = n2 + 2;/*C*/
    	tors_e += calctorq(n,n1,n2,n3,6);
	/*C-C-C-C 7*/
	for(j=0;j<3;j++){
		n = i*21 + (2*j) +1;/*C*/
		n1 = n + 2;/*C*/
		n2 = n + 4;/*C*/
		n3 = n + 6;/*C*/
    		tors_e += calctorq(n,n1,n2,n3,7);
	}
	/*C-O-C-O 8*/
	n = i*21 + 7;/*C 7*/
	n1 = i*21 + 0;/*O 0*/
	n2 = ((i+1)%7)*21 + 1;/*C 22*/
	n3 = n2 + 15;/*O 37*/
    	tors_e += calctorq(n,n1,n2,n3,8);
	/*C-O-C-O 8*/
	n = i*21 + 0;/*C 0*/
	n1 = ((i+1)%7)*21 + 1;/*O 22*/
	n2 = n1 + 15;/*C 37*/
	n3 = n1 + 8;/*O 30*/
    	tors_e += calctorq(n,n1,n2,n3,8);
//fprintf(stderr,"torsions: %f ",tors_e*KCAL);
	
/* loop over all angle bends*/
//	bend_e = 0.;
	/*H-C-C 0 */
	for (j=0;j<5;j++){
		n = i*21+(j*2)+1;/*C 1,3,5,7,9*/
		n1 = n+2;/*C 3,5,7,9,11*/
		n2 = n+3;/*H 4,6,8,10,12*/
    		bend_e += calcbend(n,n1,n2,0);
	}
	for (j=0;j<4;j++){
		n = i*21+(j*2)+4;/*H 4,6,8,10*/
		n1 = n-1;/*C 3,5,7,9*/
		n2 = n+1;/*C 5,7,9,11*/
    		bend_e += calcbend(n,n1,n2,0);
	}
	n = i*21+9;/*C*/
	n1 = i*21+11;/*C*/
	n2 = i*21+13;/*H*/
    	bend_e += calcbend(n,n1,n2,0);
	/*H-C-C 1 */
	n = i*21+2;/*H*/
	n1 = i*21+1;/*C*/
	n2 = i*21+3;/*C*/
    	bend_e += calcbend(n,n1,n2,1);
	/*C-C-C 2*/
	for (j=0;j<4;j++){
		n = i*21+2*j+1;/*1,3,5,7*/
		n1 = n+2;
		n2 = n+4;
    		bend_e += calcbend(n,n1,n2,2);
	}
	/*H-C-H 3*/
	n = i*21+12;/*H*/
	n1 = i*21+11;/*C*/
	n2 = i*21+13;/*H*/
    	bend_e += calcbend(n,n1,n2,3);
	/*C-C-O 4*/
	n = i*21+5;/*C*/
	n1 = i*21+7;/*C*/
	n2 = i*21+0;/*O*/
    	bend_e += calcbend(n,n1,n2,4);
	n = n + 4;/*C (9)*/
    	bend_e += calcbend(n,n1,n2,4);
	n = i*21+3;/*C*/
	n1 = i*21+1;/*C*/
	n2 = i*21+16;/*O*/
    	bend_e += calcbend(n,n1,n2,4);
	n = i*21+7;/*C*/
	n1 = i*21+9;/*C*/
    	bend_e += calcbend(n,n1,n2,4);
	n = i*21+11;/*C*/
    	bend_e += calcbend(n,n1,n2,4);
	n = i*21+0;/*O*/
	n1 = ((i+1)%7)*21+1;/*C 22*/
	n2 = n1+2;/*C 24*/
    	bend_e += calcbend(n,n1,n2,4);
	/*H-C-O 5*/
	n = i*21+2;/*H*/
	n1 = i*21+1;/*C*/
	n2 = i*21+16;/*O*/
    	bend_e += calcbend(n,n1,n2,5);
	n = i*21+10;/*H*/
	n1 = i*21+9;/*C*/
    	bend_e += calcbend(n,n1,n2,5);
	n = i*21+8;/*H*/
	n1 = i*21+7;/*C*/
	n2 = i*21+0;/*O*/
    	bend_e += calcbend(n,n1,n2,5);
	n = i*21+0;/*H*/
	n1 = ((i+1)%7)*21+1;/*C 22*/
	n2 = n1+1;/*O 23*/
    	bend_e += calcbend(n,n1,n2,5);
	/*C-C-O 6*/
	n1 = i*21+5;/*C*/
	n2 = i*21+17;/*O*/
	for(j=0;j<2;j++){
		n = i*21+4*j+3;/*C 3,7*/
    		bend_e += calcbend(n,n1,n2,6);
	}
	n1 = i*21+3;/*C*/
	n2 = i*21+19;/*O*/
	for(j=0;j<2;j++){
		n = i*21+4*j+1;/*C 1,5*/
    		bend_e += calcbend(n,n1,n2,6);
	}
	n = i*21+9;/*C*/
	n1 = i*21+11;/*C*/
	n2 = i*21+14;/*O*/
    	bend_e += calcbend(n,n1,n2,6);
	/*C-O-H 7*/
	n = i*21+3;/*C*/
	n1 = i*21+19;/*O*/
	n2 = i*21+20;/*H*/
    	bend_e += calcbend(n,n1,n2,7);
	n = i*21+5;/*C*/
	n1 = i*21+17;/*O*/
	n2 = i*21+18;/*H*/
    	bend_e += calcbend(n,n1,n2,7);
	n = i*21+11;/*C*/
	n1 = i*21+14;/*O*/
	n2 = i*21+15;/*H*/
    	bend_e += calcbend(n,n1,n2,7);
	/*H-C-O 8*/
	n = i*21+4;/*H*/
	n1 = i*21+3;/*C*/
	n2 = i*21+19;/*O*/
    	bend_e += calcbend(n,n1,n2,8);
	n = i*21+6;/*H*/
	n1 = i*21+5;/*C*/
	n2 = i*21+17;/*O*/
    	bend_e += calcbend(n,n1,n2,8);
	n = i*21+12;/*H*/
	n1 = i*21+11;/*C*/
	n2 = i*21+14;/*O*/
    	bend_e += calcbend(n,n1,n2,8);
	n = i*21+13;/*H*/
    	bend_e += calcbend(n,n1,n2,8);
	/*C-O-C 9*/
	n = i*21+1;/*C*/
	n1 = i*21+16;/*O*/
	n2 = i*21+9;/*C*/
    	bend_e += calcbend(n,n1,n2,9);
	n = i*21+7;/*C*/
	n1 = i*21+0;/*O*/
	n2 = ((i+1)%7)*21+1;/*C*/
    	bend_e += calcbend(n,n1,n2,9);
	/*O-C-O 10*/
	n = i*21+0;/*O*/
	n1 = ((i+1)%7)*21+1;/*C*/
	n2 = n1+15;/*O*/
    	bend_e += calcbend(n,n1,n2,10);
//fprintf(stderr,"bends: %f /n",bend_e*KCAL);
	/* 1-4 L-J pairs - intra-unit */
//	lj14_e=0.0;
//fprintf(stderr,"test fprintf\n");  // this causes a segfault (?)
//	for(j=0;j<0;j++){
	for(j=0;j<53;j++){
		eg = 0.0;
		m=i*21+pair[j][0]+nsolvent;
		n=i*21+pair[j][1]+nsolvent;
//if(m>147 || n>147) fprintf(stderr,"BCDForce.c m = %d, n = %d, pair = %d\n",m,n,j);
//if(m>147 || n>147) fprintf(stderr,"BCDForce.c m = %d, n = %d, pair = %d\n",pair[j][0],pair[j][1],j);
		sdist.fx = pos[m].fx - pos[n].fx;
		sdist.fy = pos[m].fy - pos[n].fy;
		sdist.fz = pos[m].fz - pos[n].fz;
		mvimage(&sdist);
		dx = sdist.fx;
		dy = sdist.fy;
		dz = sdist.fz;
		r2 = dx*dx+dy*dy+dz*dz;
//fprintf(stderr,"test fprintf\n");  // this causes a segfault (?)
		dedr = ttraljq14(r2,pair[j][0],pair[j][1],&eg); // pass index 0-21
//fprintf(stderr,"BCDForce 1-4LJ: passed m = %d, n = %d \n",pair[j][0],pair[j][1]);
//if(m>147 || n>147) fprintf(stderr,"BCDForce.c got pair dedr\n",m,n,j);
		force[m].fx-=dedr*dx;
		force[m].fy-=dedr*dy;
		force[m].fz-=dedr*dz;
		force[n].fx+=dedr*dx;
		force[n].fy+=dedr*dy;
		force[n].fz+=dedr*dz;
//		fprintf(stderr,"BCD ljq14 m = %d, n = %d. r = %f eg = %f\n",m,n,sqrt(r2),eg*KCAL);
		lj14_e += eg;
	}
	/* 1-4 L-J pairs - inter-unit */
//	for(j=0;j<0;j++){
	for(j=53;j<63;j++){
		eg = 0.0;
		m=i*21+pair[j][0]+nsolvent;
		n=((i+1)%7)*21-21+pair[j][1]+nsolvent;
//if(m>147 || n>147) fprintf(stderr,"BCDForce.c m = %d, n = %d, pair = %d\n",m,n,j);
		sdist.fx = pos[m].fx - pos[n].fx;
		sdist.fy = pos[m].fy - pos[n].fy;
		sdist.fz = pos[m].fz - pos[n].fz;
		mvimage(&sdist);
		dx = sdist.fx;
		dy = sdist.fy;
		dz = sdist.fz;
		r2 = dx*dx+dy*dy+dz*dz;
//		n1 = pair[j][1]-21;
//fprintf(stderr,"BCDForce 1-4LJ: passed m = %d, n = %d\n",pair[j][0],n1);
		dedr = ttraljq14(r2,pair[j][0],pair[j][1]-21,&eg); //pass index 0-21
//if(m>147 || n>147) fprintf(stderr,"BCDForce.c got pair dedr\n",m,n,j);
		force[m].fx-=dedr*dx;
		force[m].fy-=dedr*dy;
		force[m].fz-=dedr*dz;
		force[n].fx+=dedr*dx;
		force[n].fy+=dedr*dy;
		force[n].fz+=dedr*dz;
//		fprintf(stderr,"BCD ljq14 m = %d, n = %d. r = %f eg = %f\n",m,n,sqrt(r2),eg*KCAL);
		lj14_e += eg;
	}
}
// non-bonding L-J q interactions
ljqNB_e = 0.0;
//fprintf(stderr,"nsolvent = %d\n",nsolvent);
for(i=0;i<9849;i++){
	eg = 0.0;
	m=pairNB[i][0]+nsolvent;
	n=pairNB[i][1]+nsolvent;

	sdist.fx = pos[m].fx - pos[n].fx;
	sdist.fy = pos[m].fy - pos[n].fy;
	sdist.fz = pos[m].fz - pos[n].fz;
	mvimage(&sdist);
	dx = sdist.fx;
	dy = sdist.fy;
	dz = sdist.fz;

	r2 = dx*dx+dy*dy+dz*dz;
	dedr = ttraljq(r2,((pairNB[i][0])%21),((pairNB[i][1])%21),&eg); //pass index 0-21
	force[m].fx-=dedr*dx;
	force[m].fy-=dedr*dy;
	force[m].fz-=dedr*dz;
	force[n].fx+=dedr*dx;
	force[n].fy+=dedr*dy;
	force[n].fz+=dedr*dz;
//if(tc==0)
//   fprintf(stderr,"1-5 energy: %f\n",eg*KCAL);
	ljqNB_e += eg;
}

/* total sum for all molecules*/
	bondE += bond_e;
	bendE += bend_e;
	torsE += tors_e;
	nb14E += lj14_e;
	nb15E += ljqNB_e;
//if(tc==0)
//   fprintf(stderr,"1-5 total energy: %f 1-4: %f\n",ljqNB_e*KCAL,lj14_e*KCAL);
/*
fprintf(stderr,"bonds: %f \n",bond_e*KCAL);
fprintf(stderr,"bends: %f \n",bend_e*KCAL);
fprintf(stderr,"tors:  %f \n",tors_e*KCAL);
fprintf(stderr,"LJ1-4: %f \n",lj14_e*KCAL);
fprintf(stderr,"LJqNB: %f \n",ljqNB_e*KCAL);
*/
/* total intramolecular potential energy*/
    	INTRABCD += bond_e + tors_e + bend_e + lj14_e + ljqNB_e;
}

double calcstrch(i,j,k)
int i,j,k;
{
int nsolvent;
nsolvent = natoms - nBCD*BCDs - nsolute;
i+=nsolvent;
j+=nsolvent;
tripd 	d, sdist; 
double	bond, dedr, pe;
//bonds++;
sdist.fx = pos[i].fx - pos[j].fx;
sdist.fy = pos[i].fy - pos[j].fy;
sdist.fz = pos[i].fz - pos[j].fz;
mvimage(&sdist);
d.fx = sdist.fx;
d.fy = sdist.fy;
d.fz = sdist.fz;
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
double calctorq (i1,i2,i3,i4,m)
int i1,i2,i3,i4,m;
{
int nsolvent;
nsolvent = natoms - nBCD*BCDs - nsolute;
i1+=nsolvent;
i2+=nsolvent;
i3+=nsolvent;
i4+=nsolvent;
//torsions++;
tripd sdist;
int i,j,k,index,ats[4];
double dmat[3][3], 
      cmat[3][3],
      grad[4][3], 
      d[3][3], 
      da,cosa, 
      K1,K10,K11,K12,K13,K21,K22,K23, 
      D1201, C1100, 
      dudcos, 
      pe,V1,V2;
ats[0] = i1;
ats[1] = i2;
ats[2] = i3;
ats[3] = i4;
    for (i=0; i<3; i++){/* get vectors for atoms */
        sdist.fx = pos[ats[i+1]].fx-pos[ats[i]].fx;
        sdist.fy = pos[ats[i+1]].fy-pos[ats[i]].fy;
        sdist.fz = pos[ats[i+1]].fz-pos[ats[i]].fz;   
	mvimage(&sdist);
	d[i][0] = sdist.fx;
	d[i][1] = sdist.fy;
	d[i][2] = sdist.fz;	
    }
    for (i=0; i<3; i++){
        for (j=0; j<3; j++){
	    cmat[i][j]=0.0;
            for (k=0;  k<3;  k++){
                cmat[i][j] += d[i][k]*d[j][k];
            }
        }
    }
    for (i=0; i<3; i++){
        for (j=0; j<3; j++){
            dmat[i][j] = cmat[i][i]*cmat[j][j]-cmat[i][j]*cmat[i][j];
        }
    }

    
    K1 = cmat[0][1]*cmat[1][2] - cmat[0][2]*cmat[1][1];
    D1201 = sqrt (dmat[1][2]*dmat[0][1]);

    if (D1201 == 0){
	fprintf(stderr,"No torsion defined\n");
	exit(1);
    }


    /*
    **	    minus gradient at atom1 with respect to the cosine of the angle
    */

    for (i=0;  i<3;  i++){
        grad[0][i]=( ( (K1* (cmat[0][1] * d[1][i]-cmat[1][1] * d[0][i]) )
	/dmat[0][1]) + cmat[1][2]*d[1][i] - cmat[1][1]*d[2][i])/D1201;
    }

    /*
    **	    minus gradient at atom2 with respect to the cosine of the angle
    */

    K11 = cmat[0][1] + cmat[1][1];
    K12 = cmat[1][2] + 2*cmat[0][2];
    K13 = cmat[0][0] + cmat[0][1];
   
    for (i=0;  i<3;  i++){
	V1 = cmat[1][2]*d[2][i] - cmat[2][2]*d[1][i];
	V2 = K11*d[0][i] - K13*d[1][i];
        grad[1][i] = (K11*d[2][i] - K12*d[1][i] + cmat[1][2]*d[0][i] + 
			K1*(V1/dmat[1][2] + V2/dmat[0][1]))/D1201;
    }

    /*
    **	   minus gradient at atom3 with respect to the cosine of the angle
    */
    K21= cmat[0][1] + 2*cmat[0][2];
    K22= cmat[1][2] + cmat[1][1];
    K23= cmat[2][2] + cmat[1][2];
    for (i=0;  i<3;  i++){
        V1 = cmat[0][0]*d[1][i] - cmat[0][1]*d[0][i];
	V2 = K23*d[1][i] - K22*d[2][i];
	grad[2][i] = (K21*d[1][i] - K22*d[0][i] - cmat[0][1]*d[2][i] +
			K1*(V1/dmat[0][1] + V2/dmat[1][2]))/D1201;
    }

    /*
    **	    minus gradient at atom4 with respect to the cosine of the angle
    */
    for (i=0;  i<3;  i++){
        grad[3][i] = (-cmat[0][1]*d[1][i] + cmat[1][1]*d[0][i] + 
	    K1*(cmat[1][1]*d[2][i] - cmat[1][2]*d[1][i])/dmat[1][2])/D1201;
    }
    /*
    **	    calculate cos angle
    */
    cosa = K1/D1201;
    if(m<6){
       /*
       **  evaluate dU/dcos(a)
       **  U = pTors[l][0]*(1+cos(a))+pTors[l][1]*(1-cos(2*a))+pTors[l][2]*(1+cos(3*a))
       */
       pe = (1+cosa)*(pTors[m][0]+2*pTors[m][1]*(1-cosa)+pTors[m][2]*sq(2*cosa-1));
       dudcos = pTors[m][0] - 4*pTors[m][1]*cosa + 3*pTors[m][2]*(4*cosa*cosa-1);
    }
    else{ 
       /*
       **  evaluate dU/dcos(a) -- when m= 6,7,8... phase n1 and n2 = PI (n3 still = 0)
       **  U = pTors[l][0]*(1+cos(a-PI))+pTors[l][1]*(1-cos(2*a-PI))+pTors[l][2]*(1+cos(3*a))
       */
       pe = pTors[m][0]*(1-cosa)+pTors[m][1]*(2*cosa*cosa)+pTors[m][2]*(4.0*sq(cosa)*sq(cosa)-3.0*cosa+1.0);
       dudcos = -pTors[m][0] + 4*pTors[m][1]*cosa + pTors[m][2]*(12.0*sq(cosa)*cosa-3.0);
    }
    /*
    **	    add resulting force to atoms for angle
    */
    for (i=0;  i<4;  i++){
        force[ats[i]].fx += dudcos*grad[i][0];
	force[ats[i]].fy += dudcos*grad[i][1];
	force[ats[i]].fz += dudcos*grad[i][2];
    }
    return(pe);
}

double calcbend(i,j,k,l)
int i,j,k,l;
{
int nsolvent;
nsolvent = natoms - nBCD*BCDs - nsolute;
tripd sdist;
i+=nsolvent;
j+=nsolvent;
k+=nsolvent;
//bends++;
int i0,i1,i2;
tripd r[2], grad[3];
double pe,r1,r2,r1r2,n1n2,da,dedd;

i0 = j;/*center atom*/
i1 = i;
i2 = k;

sdist.fx = pos[i0].fx - pos[i1].fx;
sdist.fy = pos[i0].fy - pos[i1].fy;
sdist.fz = pos[i0].fz - pos[i1].fz;
mvimage(&sdist);
r[0].fx = sdist.fx;
r[0].fy = sdist.fy;
r[0].fz = sdist.fz;
sdist.fx = pos[i0].fx - pos[i2].fx;
sdist.fy = pos[i0].fy - pos[i2].fy;
sdist.fz = pos[i0].fz - pos[i2].fz;
mvimage(&sdist);
r[1].fx = sdist.fx;
r[1].fy = sdist.fy;
r[1].fz = sdist.fz;
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
/* intra 1-4 L-J */
double
ttraljq14(r2,m,n,eg)
double r2, *eg;
int m,n;
{
double der,r,ir,ir6,ec,a,b,q,aa,bb;
	ir = 1. / r2;
	ir6 = ir * ir * ir;
	a = BCDlj[m][n].a;
	b = BCDlj[m][n].b;
	q = BCDlj[m][n].q;
	aa =12*a;
	bb = 6*b;
	r = sqrt(r2);
	ec = q/r;
	*eg +=  BCDfactorlj*( a * ir6 - b ) * ir6 + BCDfactorq*ec;
	der = BCDfactorlj*((bb - aa * ir6) * ir6 * ir) - BCDfactorq*ec/r2;
	return(der);
}
/* intra L-J q NB */
double
ttraljq(r2,m,n,eg)
double r2, *eg;
int m,n;
{
double der,r,ir,ir6,ec,a,b,q,aa,bb;

	ir = 1. / r2;
	ir6 = ir * ir * ir;
	a = BCDlj[m][n].a;
	b = BCDlj[m][n].b;
	q = BCDlj[m][n].q;
	aa =12*a;
	bb = 6*b;
	r = sqrt(r2);
	ec = q/r;
	*eg +=  ( a * ir6 - b ) * ir6 + ec;
	der = (bb - aa * ir6) * ir6 * ir - ec/r2;
	return(der);
}

/* windoing potential */
void window(){
   double zeta;
   double V_w,deriv,solMass,BCDMass;
   int i,k;
   V_w = 0.0;
   zeta = fabs(BCDcom.fz - center_w) - width_w/2.0;
//fprintf(stderr,"windowing potential: zeta = %f, V_w = %f\n",zeta,V_w*KCAL);
   if(zeta > 0.0){
	V_w = pot_w*pow(zeta,power_w-1);
	deriv = V_w*power_w*sgn(BCDcom.fz - center_w);
	V_w *= zeta;
	BCDMass = 0.0;
	for(i=0;i<BCDs;i++){
	   k = (natoms - nBCD*BCDs) + i;
	   BCDMass += mass[k];
	}
	for(i=0;i<BCDs;i++){
	   k = (natoms - nBCD*BCDs) + i;
	   force[k].fz -= (deriv*mass[k]/BCDMass);
	}
	if(KillFinger==0){
	   solMass = 0.0;
           for(i=0;i<natoms-nBCD*BCDs;i++){
	      solMass += mass[i];
	   }
	   for(i=0;i<natoms-nBCD*BCDs;i++){
	      force[i].fz += (deriv * mass[i]/solMass);
	   }
	}
//fprintf(stderr,"windowing potential done\n");
   }
   VINT += V_w;
   V_window = V_w;
//fprintf(stderr,"windowing potential: zeta = %f, V_w = %f\n",zeta,V_w*KCAL);
}
