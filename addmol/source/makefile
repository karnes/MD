Cl2TS = Cl2toTS.o getmass.o globals.o
B2Cl = BrtoCl2.o getmass.o globals.o
L2C = ll_to_cub.o globals.o
SCZ = setCOMz.o getmass.o globals.o 
H2CL = HtoCl.o globals.o
BBIG = build_big_LL.o globals.o
CWB = chopWatBox.o globals.o
CWBR = chopWatBrO.o globals.o
CUT = cutLL.o globals.o
C2O = cub_to_oct.o globals.o
READF = readbin.o globals.o
WNFC = wnFC.o globals.o
TSANG = tsAngs.o globals.o
TSANG2 = tsAngs2.o globals.o
WRITF = writebin.o globals.o
ADDMF = addmol.o getmass.o rotcoord.o wfile.o globals.o
INI = insertIon.o getmass.o globals.o
INW = insertW.o getmass.o globals.o
IN3S = insert3sol.o getmass.o globals.o
INSN2 = insertsn2.o getmass.o globals.o
INSC = insertCenter.o getmass.o globals.o
INBW = insertBCDw.o getmass.o globals.o
INBWB = insertBCDwb.o getmass.o globals.o
INBLL = insertBCDll.o getmass.o globals.o
INSW = insertH2O.o globals.o
ADDW = addW.o getmass.o globals.o
XYZ = psxyz.o globals.o
XYZ0 = psxyz0.o globals.o
SILICA = silica.o 
ADDM = addMethanol.o globals.o
ZSHFT = zshift.o globals.o
MIXM = mixMeOH.o getmass.o globals.o
MIX2M = mixMeOH2.o getmass.o globals.o
MIX3M = mixMeOH3.o getmass.o globals.o
MIXA = mixMCN.o getmass.o globals.o
MIXA3 = mixMCN3.o getmass.o globals.o
FHBM = findHBM.o globals.o
FHBA = findHBA.o globals.o
H2B8 = hex2br8.o globals.o
RMM = rmMeOH.o globals.o
RMA = rmACN.o globals.o
SHR = shrinkOCT.o globals.o
MEETOH = MeOHEtOH.o globals.o
CLIBES = -lm
CFLAGS = -I../include -I/usr/include/bsd
DFLAGS =  -g

INCLUDE = ../include/globals.h ../include/atomtypes.h ../include/typedefs.h

.c.o:	;cc -c -DFLAGS='"@$@: $(DFLAGS)"' $(CFLAGS)\
		$(DFLAGS) $< -lm
		touch $@

Cl2toTS:$(Cl2TS)
	cc $(DFLAGS) $(Cl2TS) $(CLIBES) -o Cl2toTS

BrtoCl2:$(B2Cl)
	cc $(DFLAGS) $(B2Cl) $(CLIBES) -o BrtoCl2

ll_to_cub:$(L2C)
	cc $(DFLAGS) $(L2C) $(CLIBES) -o ll_to_cub

setCOMz:$(SCZ)
	cc $(DFLAGS) $(SCZ) $(CLIBES) -o setCOMz

HtoCl:$(H2CL)
	cc $(DFLAGS) $(H2CL) $(CLIBES) -o HtoCl

build_big_LL:$(BBIG)
	cc $(DFLAGS) $(BBIG) $(CLIBES) -o build_big_LL

chopWatBox:$(CWB)
	cc $(DFLAGS) $(CWB) $(CLIBES) -o chopWatBox

chopWatBrO:$(CWBR)
	cc $(DFLAGS) $(CWBR) $(CLIBES) -o chopWatBrO

cutLL:$(CUT)
	cc $(DFLAGS) $(CUT) $(CLIBES) -o cutLL

cub_to_oct:$(C2O)
	cc $(DFLAGS) $(C2O) $(CLIBES) -o cub_to_oct

wnFC:$(WNFC)
	cc $(DFLAGS) $(WNFC) $(CLIBES) -o wnFC

hex2br8:$(H2B8)
	cc $(DFLAGS) $(H2B8) $(CLIBES) -o hex2br8

rmMeOH:$(RMM)
	cc $(DFLAGS) $(RMM) $(CLIBES) -o rmMeOH

rmACN:$(RMA)
	cc $(DFLAGS) $(RMA) $(CLIBES) -o rmACN

shrinkOCT:$(SHR)
	cc $(DFLAGS) $(SHR) $(CLIBES) -o shrinkOCT

findHBM:$(FHBM)
	cc $(DFLAGS) $(FHBM) $(CLIBES) -o findHBM

findHBA:$(FHBA)
	cc $(DFLAGS) $(FHBA) $(CLIBES) -o findHBA

mixMCN:$(MIXA)
	cc $(DFLAGS) $(MIXA) $(CLIBES) -o mixMCN

mixMCN3:$(MIXA3)
	cc $(DFLAGS) $(MIXA3) $(CLIBES) -o mixMCN3

mixMeOH:$(MIXM)
	cc $(DFLAGS) $(MIXM) $(CLIBES) -o mixMeOH

mixMeOH2:$(MIX2M)
	cc $(DFLAGS) $(MIX2M) $(CLIBES) -o mixMeOH2

mixMeOH3:$(MIX3M)
	cc $(DFLAGS) $(MIX3M) $(CLIBES) -o mixMeOH3

insertH2O:$(INSW)
	cc $(DFLAGS) $(INSW) $(CLIBES) -o insertH2O

addW:$(ADDW)
	cc $(DFLAGS) $(ADDW) $(CLIBES) -o addW

insertW:$(INW)
	cc $(DFLAGS) $(INW) $(CLIBES) -o insertW

insert3sol:$(IN3S)
	cc $(DFLAGS) $(IN3S) $(CLIBES) -o insert3sol

insertsn2:$(INSN2)
	cc $(DFLAGS) $(INSN2) $(CLIBES) -o insertsn2

insertCenter:$(INSC)
	cc $(DFLAGS) $(INSC) $(CLIBES) -o inscent

insertBCDw:$(INBW)
	cc $(DFLAGS) $(INBW) $(CLIBES) -o insBCDw

insertBCDwb:$(INBWB)
	cc $(DFLAGS) $(INBWB) $(CLIBES) -o insBCDwb

insertBCDll:$(INBLL)
	cc $(DFLAGS) $(INBLL) $(CLIBES) -o insBCDll

insIon:$(INI)
	cc $(DFLAGS) $(INI) $(CLIBES) -o insertIon

readbin:$(READF)
	cc $(DFLAGS) $(READF) $(CLIBES) -o readbin

tsAngs:$(TSANG)
	cc $(DFLAGS) $(TSANG) $(CLIBES) -o tsAngs

tsAngs2:$(TSANG2)
	cc $(DFLAGS) $(TSANG2) $(CLIBES) -o tsAngs2

writebin:$(WRITF)
	cc $(DFLAGS) $(WRITF) $(CLIBES) -o writebin

addmol:$(ADDMF)
	cc $(DFLAGS) $(ADDMF) $(CLIBES) -o addmol

psxyz:$(XYZ)
	cc $(DFLAGS) $(XYZ) $(CLIBES) -o psxyz

psxyz0:$(XYZ0)
	cc $(DFLAGS) $(XYZ0) $(CLIBES) -o psxyz0

silica:$(SILICA)
	cc $((DFLAGS) $(SILICA) $(CLIBES) -o silica

addM:$(ADDM)
	cc $(DFLAGS) $(ADDM) $(CLIBES) -o addM

zshift:$(ZSHFT)
	cc $(DFLAGS) $(ZSHFT) $(CLIBES) -o zshift

MeOHEtOH:$(MEETOH)
	cc $(DFLAGS) $(MEETOH) $(CLIBES) -o MeOHEtOH

#
#	implicit dependencies
#
$(Cl2TS):	$(INCLUDE)
$(B2Cl):	$(INCLUDE)
$(L2C):	$(INCLUDE)
$(SCZ):	$(INCLUDE)
$(H2CL):	$(INCLUDE)
$(BBIG):	$(INCLUDE)
$(CWB):	$(INCLUDE)
$(CWBR):	$(INCLUDE)
$(C2O):	$(INCLUDE)
$(CUT):	$(INCLUDE)
$(WNFC):	$(INCLUDE)
$(H2B8):	$(INCLUDE)
$(READF):	$(INCLUDE)
$(TSANG):	$(INCLUDE)
$(TSANG2):	$(INCLUDE)
$(WRITF):	$(INCLUDE)
$(ADDMF):	$(INCLUDE)
$(INI):	$(INCLUDE)
$(INW):	$(INCLUDE)
$(IN3S):	$(INCLUDE)
$(INSC):	$(INCLUDE)
$(INSN2):	$(INCLUDE)
$(INBW):	$(INCLUDE)
$(INBWB):	$(INCLUDE)
$(INBLL):	$(INCLUDE)
$(ADDW):	$(INCLUDE)
$(INSW):	$(INCLUDE)
$(XYZ):	$(INCLUDE)
$(XYZ0):	$(INCLUDE)
$(ADDM):	$(INCLUDE)
$(ZSHFT):	$(INCLUDE)
$(MIXM):	$(INCLUDE)
$(MIX2M):	$(INCLUDE)
$(MIX3M):	$(INCLUDE)
$(MIXA):	$(INCLUDE)
$(MIXA3):	$(INCLUDE)
$(FHBM):	$(INCLUDE)
$(FHBA):	$(INCLUDE)
$(RMM):	$(INCLUDE)
$(RMA):	$(INCLUDE)
$(SHR):	$(INCLUDE)
$(MEETOH):	$(INCLUDE)

clean:	;-rm $(READF) $(WRITF) $(ADDMF) $(XYZ) readbin writebin addmol
