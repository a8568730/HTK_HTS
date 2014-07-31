# UNIX HTK LIBRARY MAKEFILE

# This makefile will build an archive file called HTKLib.XXX.a
# and/or a specific library test program
# containing the HTK library object modules where XXX is the
# type of machine.  
# 
# When invoked it expects to have the following environment variables
# set
#       CPU             - set to the machine name
#       HTKCC           - name of C compiler (either cc or gcc)
#       HTKCF           - compiler flags
#
# Suggested settings are
#		CPU	      HTKCC  HTKCF
# SGI:		sgi           cc     -O -DSGI_AUDIO
# HP:		hp700         cc     +e -O -I/usr/include/X11R5 -I/usr/include/audio -I/opt/audio/include -D_HPUX_SOURCE -Aa -DHPUX_AUDIO -DUSE_TIMER 
# Sun:		sun4_solaris  gcc    -ansi -I/usr/openwin/include -O2 -DSOLARIS16_AUDIO -DBSD_COMP
# Linux:	linux	      gcc    -ansi -O2 -DOSS_AUDIO
#
# When running in a multiple machine environment, type
#        make clean
# before recompiling for a new machine.

hlib = .
HLIBS = $(hlib)/HTKLib.$(CPU).a

CC = $(HTKCC)
CFLAGS = $(HTKCF) -D'ARCH="$(CPU)"' -I.

modules = HShell.o HMath.o  HSigP.o  HWave.o HAudio.o HParm.o HVQ.o  HGraf.o\
          HLabel.o HModel.o HUtil.o HTrain.o HDict.o  HLM.o   HRec.o HNet.o \
          HAdapt.o HFB.o HMem.o HLat.o \
          esignal.o esig_asc.o esig_edr.o esig_nat.o strarr.o

all:    HTKLib.$(CPU).a

HTKLib.$(CPU).a: $(modules)
	-rm HTKLib.$(CPU).a
	ar rv HTKLib.$(CPU).a $(modules)
	-ranlib HTKLib.$(CPU).a

.PHONY: clean cleanup depend
clean:
	-rm -f *.o HTKLib.$(CPU).a *.cpu
	touch $(CPU).cpu

cleanup:
	-rm -f *.o *.cpu

depend:
	-makedepend -Y *.c

# DO NOT DELETE THIS LINE -- make depend depends on it.

HAdapt.o: HShell.h HMem.h HMath.h HSigP.h HWave.h HAudio.h HParm.h HLabel.h
HAdapt.o: HModel.h HTrain.h HUtil.h HAdapt.h
HAudio.o: HShell.h HMem.h HMath.h HWave.h HAudio.h
HDict.o: HShell.h HMem.h HMath.h HWave.h HLabel.h HDict.h
HFB.o: HShell.h HMem.h HMath.h HSigP.h HAudio.h HWave.h HVQ.h HParm.h
HFB.o: HLabel.h HModel.h HTrain.h HUtil.h HAdapt.h HFB.h
HGraf.o: HShell.h HMem.h HMath.h HGraf.h
HGraf.null.o: HShell.h HMem.h HMath.h HGraf.h
HGraf_WIN32.o: HShell.h HMem.h HMath.h HGraf.h
HLM.o: HShell.h HMem.h HMath.h HWave.h HLabel.h HLM.h
HLabel.o: HShell.h HMem.h HMath.h HWave.h HLabel.h
HLat.o: HShell.h HMem.h HMath.h HWave.h HAudio.h HParm.h HLabel.h HModel.h
HLat.o: HUtil.h HDict.h HNet.h HLM.h HLat.h
HMath.o: HShell.h HMem.h HMath.h
HMem.o: HShell.h HMem.h
HModel.o: HShell.h HMem.h HMath.h HWave.h HAudio.h HParm.h HLabel.h HModel.h
HModel.o: HUtil.h HTrain.h HAdapt.h
HNet.o: HShell.h HMem.h HMath.h HWave.h HAudio.h HParm.h HLabel.h HModel.h
HNet.o: HUtil.h HDict.h HNet.h
HParm.o: HShell.h HMem.h HMath.h HSigP.h HAudio.h HWave.h HVQ.h HParm.h
HParm.o: HLabel.h HModel.h esignal.h
HRec.o: HShell.h HMem.h HMath.h HWave.h HAudio.h HParm.h HLabel.h HModel.h
HRec.o: HDict.h HNet.h HRec.h HUtil.h
HShell.o: HShell.h
HSigP.o: HShell.h HMem.h HMath.h HSigP.h
HTrain.o: HShell.h HMem.h HMath.h HSigP.h HAudio.h HWave.h HVQ.h HParm.h
HTrain.o: HLabel.h HModel.h HUtil.h HTrain.h
HUtil.o: HShell.h HMem.h HMath.h HSigP.h HAudio.h HWave.h HVQ.h HParm.h
HUtil.o: HLabel.h HModel.h HUtil.h
HVQ.o: HShell.h HMem.h HMath.h HVQ.h
HWave.o: HShell.h HMem.h HMath.h HWave.h esignal.h HAudio.h HParm.h
esig_asc.o: esignal.h
esig_edr.o: esignal.h
esig_nat.o: esignal.h
esignal.o: esignal.h
strarr.o: esignal.h strarr.h
