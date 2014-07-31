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
          HAdapt.o HFB.o HMem.o \
          esignal.o esig_asc.o esig_edr.o esig_nat.o strarr.o

all:    HTKLib.$(CPU).a

HTKLib.$(CPU).a: $(modules)
	ar rv HTKLib.$(CPU).a $?
	-ranlib HTKLib.$(CPU).a

clean:
	-rm -f *.o HTKLib.$(CPU).a *.cpu
	touch $(CPU).cpu

cleanup:
	-rm -f *.o *.cpu

