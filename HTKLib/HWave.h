/* ----------------------------------------------------------- */
/*                                                             */
/*                          ___                                */
/*                       |_| | |_/   SPEECH                    */
/*                       | | | | \   RECOGNITION               */
/*                       =========   SOFTWARE                  */ 
/*                                                             */
/*                                                             */
/* ----------------------------------------------------------- */
/*         Copyright: Microsoft Corporation                    */
/*          1995-2000 Redmond, Washington USA                  */
/*                    http://www.microsoft.com                 */
/*                                                             */
/*   Use of this software is governed by a License Agreement   */
/*    ** See the file License for the Conditions of Use  **    */
/*    **     This banner notice must not be removed      **    */
/*                                                             */
/* ----------------------------------------------------------- */
/*         File: HWave.h: Speech Waveform File Input           */
/* ----------------------------------------------------------- */

/*  *** THIS IS A MODIFIED VERSION OF HTK ***                        */
/*  ---------------------------------------------------------------  */
/*           The HMM-Based Speech Synthesis System (HTS)             */
/*                       HTS Working Group                           */
/*                                                                   */
/*                  Department of Computer Science                   */
/*                  Nagoya Institute of Technology                   */
/*                               and                                 */
/*   Interdisciplinary Graduate School of Science and Engineering    */
/*                  Tokyo Institute of Technology                    */
/*                                                                   */
/*                     Copyright (c) 2001-2006                       */
/*                       All Rights Reserved.                        */
/*                                                                   */
/*  Permission is hereby granted, free of charge, to use and         */
/*  distribute this software in the form of patch code to HTK and    */
/*  its documentation without restriction, including without         */
/*  limitation the rights to use, copy, modify, merge, publish,      */
/*  distribute, sublicense, and/or sell copies of this work, and to  */
/*  permit persons to whom this work is furnished to do so, subject  */
/*  to the following conditions:                                     */
/*                                                                   */
/*    1. Once you apply the HTS patch to HTK, you must obey the      */
/*       license of HTK.                                             */
/*                                                                   */
/*    2. The source code must retain the above copyright notice,     */
/*       this list of conditions and the following disclaimer.       */
/*                                                                   */
/*    3. Any modifications to the source code must be clearly        */
/*       marked as such.                                             */
/*                                                                   */
/*  NAGOYA INSTITUTE OF TECHNOLOGY, TOKYO INSTITUTE OF TECHNOLOGY,   */
/*  HTS WORKING GROUP, AND THE CONTRIBUTORS TO THIS WORK DISCLAIM    */
/*  ALL WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL       */
/*  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT   */
/*  SHALL NAGOYA INSTITUTE OF TECHNOLOGY, TOKYO INSTITUTE OF         */
/*  TECHNOLOGY, HTS WORKING GROUP, NOR THE CONTRIBUTORS BE LIABLE    */
/*  FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY        */
/*  DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,  */
/*  WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTUOUS   */
/*  ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR          */
/*  PERFORMANCE OF THIS SOFTWARE.                                    */
/*                                                                   */
/*  ---------------------------------------------------------------  */

/* !HVER!HWave:   3.3 [CUED 28/04/05] */

/*  Configuration Parameters:
   NSAMPLES       - num samples in alien file input via a pipe
   HEADERSIZE     - size of header in alien file
   BYTEORDER      - define byte order VAX or other
   STEREOMODE     - select LEFT or RIGHT chan (default both)
   SOURCERATE     - sample period of source
   TARGETRATE     - sample period of target
   SOURCEFORMAT   - input file format
   TARGETFORMAT   - output file format
   TRACE          - trace level
*/

#ifndef _HWAVE_H_
#define _HWAVE_H_

#ifdef __cplusplus
extern "C" {
#endif

/* Added for esignal supprt */
typedef struct FieldSpec **HFieldList;

typedef enum {
        NOHEAD,            /* Headerless File */
        HAUDIO,            /* Direct Audio Input */
        HTK,               /* used for both wave and parm files */
        TIMIT,             /* Prototype TIMIT database */
        NIST,              /* NIST databases eg RM1,TIMIT */
        SCRIBE,            /* UK Scribe databases */
        AIFF,              /* Apple Audio Interchange format */
        SDES1,             /* Sound Designer I format */
        SUNAU8,            /* Sun 8 bit MuLaw .au format */
        OGI,               /* Oregon Institute format (similar to TIMIT) */
        ESPS,              /* used for both wave and parm files */
	ESIG,              /* used for both wave and parm files */
	WAV,               /* Microsoft WAVE format */
        UNUSED,
        ALIEN,             /* Unknown */
        UNDEFF
} FileFormat;

typedef struct _Wave *Wave;  /* Abstract type representing waveform file */

void InitWave(void);
/*
   Initialise module
*/

void ResetWave(void);
/*
   Reset module
*/

Wave OpenWaveInput(MemHeap *x, char *fname, FileFormat fmt, HTime winDur, 
                   HTime frPeriod, HTime *sampPeriod);
/*
   Open the named input file with the given format and return a
   Wave object. If fmt==UNDEFF then the value of the configuration
   parameter SOURCEFORMAT is used.  If this is not set, then the format
   HTK is assumed. Samples are returned in frames of duration winDur.  The 
   period between successive frames is given by frPeriod.  If the value of 
   sampPeriod is not 0.0, then it overrides the sample period specified in
   the file, otherwise the actual value is returned. Returns NULL on error.
*/

void CloseWaveInput(Wave w);
/* 
   Terminate Wave input and free any resources associated with w
*/

void ZeroMeanWave(Wave w);
/*
   Ensure that mean of wave w is zero
*/

int FramesInWave(Wave w);
/*
   Return number of whole frames which are currently
   available in the given Wave
*/

int SampsInWaveFrame(Wave w);
/*
   Return number of samples in each frame of the given Wave
*/

void GetWave(Wave w, int nFrames, float *buf);
/* 
   Get next nFrames from Wave input buffer and store sequentially in
   buf as floats.  If a frame overlap has been set then samples will be
   duplicated in buf.  It is a fatal error to request more frames
   than exist in the Wave (as determined by FramesInWave.
*/

short *GetWaveDirect(Wave w, long *nSamples);
/* 
   Returns a pointer to the waveform stored in w.
*/

Wave OpenWaveOutput(MemHeap *x, HTime *sampPeriod, long bufSize);
/*
   Initialise a Wave object to store waveform data at the given 
   sample period, using buffer of bufSize shorts.  
*/

void PutWaveSample(Wave w, long nSamples, short *buf);
/*
   Append given nSamples in buf to wave w.
*/

ReturnStatus CloseWaveOutput(Wave w, FileFormat fmt, char *fname);
/* 
   Output wave w to file fname in given fmt and free any 
   associated resources.  If fmt==UNDEFF then value of
   configuration variable TARGETFORMAT is used, if any,
   otherwise the HTK format is used. If an error then 
   returns FAIL and does not free any memory. 
*/

FileFormat WaveFormat(Wave w);
/* 
   Return format of given wave
*/

char *Format2Str(FileFormat format);
FileFormat Str2Format(char *fmt);
/*
   Convert between FileFormat enum type & string.
*/

/* --------------------- HTK Header Routines --------------------- */

Boolean ReadHTKHeader(FILE *f,long *nSamp,long *sampP,short *sampS,
                      short *kind, Boolean *bSwap);
/* 
   Get header info from HTK file f, return false if apparently not
   a HTK file.  If byte-swapped bswap returns true.  NB only
   the user can specify required byte order via NATREADORDER config var 
   since it is not defined for HTK files)
*/

void WriteHTKHeader(FILE *f, long nSamp, long sampP, short sampS, 
		    short kind, Boolean *bSwap);
/* 
   Write header info to HTK file f.  
   Sets bSwap to indicate whether header was byte swapped before writing.
*/

void StoreESIGFieldList(HFieldList fList);
/*
   Store the field list of an ESIG input file 
*/
void RetrieveESIGFieldList(HFieldList *fList);
/*
   Retrieve the field list of an ESIG input file 
*/

Boolean ReadEsignalHeader(FILE *f, long *nSamp, long *sampP, short *sampS,
 			  short *kind, Boolean *bSwap, long *hdrS,
 			  Boolean isPipe);
/*
    Get header from Esignal file f; return FALSE in case of failure.
*/

#ifdef __cplusplus
}
#endif

#endif  /* _HWAVE_H_ */

/* ------------------------ End of HWave.h ----------------------- */

