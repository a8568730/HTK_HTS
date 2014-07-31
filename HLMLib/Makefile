# UNIX HLM LIBRARY MAKEFILE


HLMLIB = HLMLib.$(CPU).a
HTKLIB = ../HTKLib/HTKLib.$(CPU).a

CC = $(HTKCC)
CFLAGS = $(HTKCF) -I. -I../HTKLib/ \
	-DSANITY $(HAPI_DEFS) -g

hlmmod = LModel.${CPU}.o LPMerge.${CPU}.o LPCalc.${CPU}.o LUtil.${CPU}.o \
	LWMap.${CPU}.o LCMap.${CPU}.o LGBase.${CPU}.o

HEADERS = LCMap.h LGBase.h LModel.h LPCalc.h LPMerge.h LUtil.h LWMap.h

# the following builds the HLM library
$(HLMLIB):	$(hlmmod)
		ar rv $(HLMLIB) $?
		ranlib $(HLMLIB)

LUtil.${CPU}.o:		LUtil.c $(HEADERS)
			$(CC) $(CFLAGS) -c LUtil.c -o LUtil.$(CPU).o
LWMap.${CPU}.o:		LWMap.c $(HEADERS)
			$(CC) $(CFLAGS) -c LWMap.c -o LWMap.$(CPU).o
LCMap.${CPU}.o:		LCMap.c $(HEADERS)
			$(CC) $(CFLAGS) -c LCMap.c -o LCMap.$(CPU).o
LGBase.${CPU}.o:	LGBase.c $(HEADERS)
			$(CC) $(CFLAGS) -c LGBase.c -o LGBase.$(CPU).o
LPCalc.${CPU}.o:	LPCalc.c $(HEADERS)
			$(CC) $(CFLAGS) -c LPCalc.c -o LPCalc.$(CPU).o
LPMerge.${CPU}.o:	LPMerge.c $(HEADERS)
			$(CC) $(CFLAGS) -c LPMerge.c -o LPMerge.$(CPU).o
LModel.${CPU}.o:	LModel.c $(HEADERS)
			$(CC) $(CFLAGS) -c LModel.c -o LModel.$(CPU).o

all: 			$(HLMLIB)

clean:
			rm -f *.${CPU}.o $(HLMLIB)
cleanup:
			rm -f *.${CPU}.o
