UDIR = ../../util/
DDIR = ../../data/
MIDIR = ../../min/
WDIR = ../../water/
SDIR = ../
IDIR = ../../include/
ALLD  = $(UDIR) $(DDIR) $(FDIR) $(SDIR) $(MIDIR) $(ADIR) $(WDIR)

SYS = 	$(SDIR)sysInit.o $(SDIR)liqForce.o $(SDIR)intForce.o $(SDIR)wxtr.o\
	$(SDIR)sysForce.o $(SDIR)doSpec.o $(SDIR)wstat.o $(SDIR)ElecPot.o\
	$(SDIR)moveZW.o $(SDIR)NITROForce.o  $(SDIR)waterNITRO.o $(SDIR)waterForceWNIT.o\
	$(SDIR)wdatNITRO.o

UTIL =  $(UDIR)getmass.o $(UDIR)life.o $(UDIR)ranvel.o\
	$(UDIR)rinp.o $(UDIR)settime.o $(UDIR)swtch.o $(UDIR)offset.o\
	$(UDIR)doTraj.o $(UDIR)getforce.o $(UDIR)crossp.o $(UDIR)save.o\
	$(UDIR)kill.o $(UDIR)main.o $(UDIR)doEqu.o $(UDIR)winp.o\
	$(UDIR)init.o $(UDIR)mvimage.o $(UDIR)verlet.o $(UDIR)constrain.o\
	$(UDIR)reStart.o $(UDIR)rvel.o $(UDIR)wrsv.o $(UDIR)equil.o\
	$(UDIR)dist.o $(UDIR)pimage.o

MIN =	$(MIDIR)conjug.o $(MIDIR)steepest.o $(UDIR)doMin.o $(SDIR)wstat_min.o

WATER = $(WDIR)waterInit.o

DATA =	$(DDIR)collect.o $(DDIR)wfrc.o $(DDIR)wfsolute.o\
	$(DDIR)wpos.o $(DDIR)wvel.o $(DDIR)wshel.o $(DDIR)wpsolute.o $(DDIR)wvsolute.o

FILES_O = $(SYS) $(UTIL) $(DATA) $(MIN) $(WATER)

CLIBES = -lm 
CFLAGS = -I../../include -I../ -I/usr/include/bsd
DFLAGS = -O2 -DNEWSPC  

INCLUDE = 	$(IDIR)atom.h $(IDIR)md.h $(IDIR)potentials.h\
		$(IDIR)atomtypes.h $(IDIR)typedefs.h $(IDIR)units.h\
		$(SDIR)system.h $(IDIR)water.h

.c.o:	;cc -c -DFLAGS='"@$@: $(DFLAGS)"' $(CFLAGS)\
		$(DFLAGS) $< -lm
		mv $(*F).o $@
		touch $@

run:	$(FILES_O)
	cc $(DFLAGS) $(FILES_O) $(CLIBES) -o run
	size run
	mv run ../run/run

#
#	implicit dependencies
#
$(SYS):		$(INCLUDE)
$(UTIL):	$(INCLUDE)
$(DATA):	$(INCLUDE)
$(MIN):		$(INCLUDE)

cleanup_run:	;-rm $(FILES_O)
