# UNIX HTK TOOLS MAKEFILE
#
# 
# This makefile will build executables for each HTK Basic
# tools.  To build a specific tool type
#       make HTool
# To build all tools type
#       make all
# When compiled on a machine XXX, the executables are stored in
# a directory called bin.XXX which must already exist.
#
# When invoked it expects to have the following environment variables
# set
#       CPU             - set to the machine name
#       HTKCC           - name of C compiler (either cc or gcc)
#       HTKCF           - compiler flags
#       HTKLF           - linker flags
#       HBIN            - the directory containing "bin.$CPU"
#
# Suggested settings are
#       CPU           HTKCC  HTKCF                                                            HTKLF
# SGI:  sgi           cc     -O -DSGI_AUDIO                                                   -lmalloc -laudio
# HP:   hp700         cc     +e -O -I/usr/include/X11R5 -I/usr/include/audio -I/opt/audio/include -D_HPUX_SOURCE -Aa -DHPUX_AUDIO -DUSE_TIMER   -lAlib
# Sun:  sun4_solaris  gcc    -ansi -I/usr/openwin/include -O2 -DSOLARIS16_AUDIO -DBSD_COMP    -ldl -lsocket -lnsl -lintl
# Linux:linux	      gcc    -ansi -O2 -DOSS_AUDIO					      -L/usr/X11R6/lib
#
# When running in a multiple machine environment, type
#        make clean
# before recompiling for a new machine and
#        make cleanup
# afterwards.
#

hlib = ../HTKLib

HLIBS = $(hlib)/HTKLib.$(CPU).a 

CC = $(HTKCC)
CFLAGS = -I$(hlib) $(HTKCF) 

all:    HEAdapt	HList	HCompV HCopy   HDMan HERest HHEd    HInit  \
        HParse HBuild HQuant  HRest  HResults HSmooth HVite HLEd   \
	HLRescore HLStats HSGen  HSLab

HEAdapt :  $(hlib)/HTKLib.$(CPU).a HEAdapt.o
	$(CC) HEAdapt.o $(HLIBS) -lm $(HTKLF)
	mv a.out $(HBIN)/bin.$(CPU)/HEAdapt

HList:  $(hlib)/HTKLib.$(CPU).a HList.o
	$(CC) HList.o $(HLIBS) -lm $(HTKLF)
	mv a.out $(HBIN)/bin.$(CPU)/HList

HCompV:  $(hlib)/HTKLib.$(CPU).a HCompV.o
	$(CC) HCompV.o $(HLIBS) -lm $(HTKLF)
	mv a.out $(HBIN)/bin.$(CPU)/HCompV

HCopy:  $(hlib)/HTKLib.$(CPU).a HCopy.o
	$(CC) HCopy.o $(HLIBS) -lm $(HTKLF)
	mv a.out $(HBIN)/bin.$(CPU)/HCopy

HDMan:  $(hlib)/HTKLib.$(CPU).a HDMan.o
	$(CC) HDMan.o $(HLIBS) -lm $(HTKLF)
	mv a.out $(HBIN)/bin.$(CPU)/HDMan

HERest: $(hlib)/HTKLib.$(CPU).a HERest.o
	$(CC) HERest.o $(HLIBS) -lm $(HTKLF)
	mv a.out $(HBIN)/bin.$(CPU)/HERest

HLEd:   $(hlib)/HTKLib.$(CPU).a HLEd.o
	$(CC) HLEd.o $(HLIBS) -lm $(HTKLF)
	mv a.out $(HBIN)/bin.$(CPU)/HLEd

HLRescore:   $(hlib)/HTKLib.$(CPU).a HLRescore.o
	$(CC) HLRescore.o $(HLIBS) -lm $(HTKLF)
	mv a.out $(HBIN)/bin.$(CPU)/HLRescore

HHEd:   $(hlib)/HTKLib.$(CPU).a HHEd.o
	$(CC) HHEd.o $(HLIBS) -lm $(HTKLF)
	mv a.out $(HBIN)/bin.$(CPU)/HHEd

HInit:  $(hlib)/HTKLib.$(CPU).a HInit.o
	$(CC) HInit.o $(HLIBS) -lm $(HTKLF)
	mv a.out $(HBIN)/bin.$(CPU)/HInit

HParse: $(hlib)/HTKLib.$(CPU).a HParse.o
	$(CC) HParse.o $(HLIBS) -lm $(HTKLF)
	mv a.out $(HBIN)/bin.$(CPU)/HParse

HBuild: $(hlib)/HTKLib.$(CPU).a HBuild.o
	$(CC) HBuild.o $(HLIBS) -lm $(HTKLF)
	mv a.out $(HBIN)/bin.$(CPU)/HBuild

HQuant: $(hlib)/HTKLib.$(CPU).a HQuant.o
	$(CC) HQuant.o $(HLIBS) -lm $(HTKLF)
	mv a.out $(HBIN)/bin.$(CPU)/HQuant

HRest:  $(hlib)/HTKLib.$(CPU).a HRest.o
	$(CC) HRest.o $(HLIBS) -lm $(HTKLF)
	mv a.out $(HBIN)/bin.$(CPU)/HRest

HResults:       $(hlib)/HTKLib.$(CPU).a HResults.o
	$(CC) HResults.o $(HLIBS) -lm $(HTKLF)
	mv a.out $(HBIN)/bin.$(CPU)/HResults

HSmooth:	$(hlib)/HTKLib.$(CPU).a HSmooth.o
	$(CC) HSmooth.o $(HLIBS) -lm $(HTKLF)
	mv a.out $(HBIN)/bin.$(CPU)/HSmooth

HVite:  $(hlib)/HTKLib.$(CPU).a HVite.o
	$(CC) HVite.o $(HLIBS) -lm $(HTKLF)
	mv a.out $(HBIN)/bin.$(CPU)/HVite

HLStats:	$(hlib)/HTKLib.$(CPU).a HLStats.o
	$(CC) HLStats.o $(HLIBS) -lm $(HTKLF)
	mv a.out $(HBIN)/bin.$(CPU)/HLStats

HSGen:  $(hlib)/HTKLib.$(CPU).a HSGen.o
	$(CC) HSGen.o $(HLIBS) -lm $(HTKLF)
	mv a.out $(HBIN)/bin.$(CPU)/HSGen

HSLab:  $(hlib)/HTKLib.$(CPU).a HSLab.o
	$(CC) HSLab.o $(HLIBS) -lm $(HTKLF) -lX11
	mv a.out $(HBIN)/bin.$(CPU)/HSLab

clean:
	-rm -f *.o *.cpu
	touch $(CPU).cpu

cleanup:
	-rm -f *.o *.cpu

removetools:
	-rm -f $(HBIN)/bin.$(CPU)/HEAdapt
	-rm -f $(HBIN)/bin.$(CPU)/HList
	-rm -f $(HBIN)/bin.$(CPU)/HCompV
	-rm -f $(HBIN)/bin.$(CPU)/HCopy
	-rm -f $(HBIN)/bin.$(CPU)/HDMan
	-rm -f $(HBIN)/bin.$(CPU)/HERest
	-rm -f $(HBIN)/bin.$(CPU)/HLEd
	-rm -f $(HBIN)/bin.$(CPU)/HLRescore
	-rm -f $(HBIN)/bin.$(CPU)/HHEd
	-rm -f $(HBIN)/bin.$(CPU)/HInit
	-rm -f $(HBIN)/bin.$(CPU)/HParse
	-rm -f $(HBIN)/bin.$(CPU)/HBuild
	-rm -f $(HBIN)/bin.$(CPU)/HQuant
	-rm -f $(HBIN)/bin.$(CPU)/HRest
	-rm -f $(HBIN)/bin.$(CPU)/HResults
	-rm -f $(HBIN)/bin.$(CPU)/HSmooth
	-rm -f $(HBIN)/bin.$(CPU)/HVite
	-rm -f $(HBIN)/bin.$(CPU)/HLStats
	-rm -f $(HBIN)/bin.$(CPU)/HSGen
	-rm -f $(HBIN)/bin.$(CPU)/HSLab
