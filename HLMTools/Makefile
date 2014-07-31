# UNIX HTK LANGUAGE MODELLING TOOLS MAKEFILE
#
# 
# This makefile will build executables for each HTK LM
# tool.  To build a specific tool type
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
# afterwards.  Or type
#        make removetools
# to delete just the compiled tools.
#

hlib = ../HTKLib
llib = ../HLMLib

HLIBS = $(hlib)/HTKLib.$(CPU).a $(llib)/HLMLib.$(CPU).a

CC = $(HTKCC)
CFLAGS = -I$(hlib) -I$(llib) $(HTKCF) 

all:    Cluster LAdapt LFoF LGList LMerge LNorm HLMCopy \
	LBuild LGCopy LGPrep LPlex LSubset LLink LNewMap

Cluster:	$(hlib)/HTKLib.$(CPU).a $(llib)/HLMLib.$(CPU).a Cluster.o
	$(CC) Cluster.o $(HLIBS) -lm $(HTKLF)
	mv a.out $(HBIN)/bin.$(CPU)/Cluster

LAdapt:	$(hlib)/HTKLib.$(CPU).a $(llib)/HLMLib.$(CPU).a LAdapt.o
	$(CC) LAdapt.o $(HLIBS) -lm $(HTKLF)
	mv a.out $(HBIN)/bin.$(CPU)/LAdapt

LFoF:	$(hlib)/HTKLib.$(CPU).a $(llib)/HLMLib.$(CPU).a LFoF.o
	$(CC) LFoF.o $(HLIBS) -lm $(HTKLF)
	mv a.out $(HBIN)/bin.$(CPU)/LFoF

LGList:	$(hlib)/HTKLib.$(CPU).a $(llib)/HLMLib.$(CPU).a LGList.o
	$(CC) LGList.o $(HLIBS) -lm $(HTKLF)
	mv a.out $(HBIN)/bin.$(CPU)/LGList

LMerge:	$(hlib)/HTKLib.$(CPU).a $(llib)/HLMLib.$(CPU).a LMerge.o
	$(CC) LMerge.o $(HLIBS) -lm $(HTKLF)
	mv a.out $(HBIN)/bin.$(CPU)/LMerge

LNorm:	$(hlib)/HTKLib.$(CPU).a $(llib)/HLMLib.$(CPU).a LNorm.o
	$(CC) LNorm.o $(HLIBS) -lm $(HTKLF)
	mv a.out $(HBIN)/bin.$(CPU)/LNorm

HLMCopy:	$(hlib)/HTKLib.$(CPU).a $(llib)/HLMLib.$(CPU).a HLMCopy.o
	$(CC) HLMCopy.o $(HLIBS) -lm $(HTKLF)
	mv a.out $(HBIN)/bin.$(CPU)/HLMCopy

LBuild:	$(hlib)/HTKLib.$(CPU).a $(llib)/HLMLib.$(CPU).a LBuild.o
	$(CC) LBuild.o $(HLIBS) -lm $(HTKLF)
	mv a.out $(HBIN)/bin.$(CPU)/LBuild

LGCopy:	$(hlib)/HTKLib.$(CPU).a $(llib)/HLMLib.$(CPU).a LGCopy.o
	$(CC) LGCopy.o $(HLIBS) -lm $(HTKLF)
	mv a.out $(HBIN)/bin.$(CPU)/LGCopy

LGPrep:	$(hlib)/HTKLib.$(CPU).a $(llib)/HLMLib.$(CPU).a LGPrep.o
	$(CC) LGPrep.o $(HLIBS) -lm $(HTKLF)
	mv a.out $(HBIN)/bin.$(CPU)/LGPrep

LPlex:	$(hlib)/HTKLib.$(CPU).a $(llib)/HLMLib.$(CPU).a LPlex.o
	$(CC) LPlex.o $(HLIBS) -lm $(HTKLF)
	mv a.out $(HBIN)/bin.$(CPU)/LPlex

LSubset:	$(hlib)/HTKLib.$(CPU).a $(llib)/HLMLib.$(CPU).a LSubset.o
	$(CC) LSubset.o $(HLIBS) -lm $(HTKLF)
	mv a.out $(HBIN)/bin.$(CPU)/LSubset

LLink:	$(hlib)/HTKLib.$(CPU).a $(llib)/HLMLib.$(CPU).a LLink.o
	$(CC) LLink.o $(HLIBS) -lm $(HTKLF)
	mv a.out $(HBIN)/bin.$(CPU)/LLink

LNewMap:	$(hlib)/HTKLib.$(CPU).a $(llib)/HLMLib.$(CPU).a LNewMap.o
	$(CC) LNewMap.o $(HLIBS) -lm $(HTKLF)
	mv a.out $(HBIN)/bin.$(CPU)/LNewMap

.PHONY:	clean cleanup removetools

clean:
	-rm -f *.o *.cpu
	touch $(CPU).cpu

cleanup:
	-rm -f *.o *.cpu

removetools:
	-rm -f $(HBIN)/bin.$(CPU)/Cluster
	-rm -f $(HBIN)/bin.$(CPU)/LAdapt
	-rm -f $(HBIN)/bin.$(CPU)/LFoF
	-rm -f $(HBIN)/bin.$(CPU)/LGList
	-rm -f $(HBIN)/bin.$(CPU)/LMerge
	-rm -f $(HBIN)/bin.$(CPU)/LNorm
	-rm -f $(HBIN)/bin.$(CPU)/HLMCopy
	-rm -f $(HBIN)/bin.$(CPU)/LBuild
	-rm -f $(HBIN)/bin.$(CPU)/LGCopy
	-rm -f $(HBIN)/bin.$(CPU)/LGPrep
	-rm -f $(HBIN)/bin.$(CPU)/LPlex
	-rm -f $(HBIN)/bin.$(CPU)/LSubset
	-rm -f $(HBIN)/bin.$(CPU)/LLink
	-rm -f $(HBIN)/bin.$(CPU)/LNewMap
