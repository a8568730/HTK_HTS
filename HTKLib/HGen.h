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
/*  distribute this software and its documentation without           */
/*  restriction, including without limitation the rights to use,     */
/*  copy, modify, merge, publish, distribute, sublicense, and/or     */
/*  sell copies of this work, and to permit persons to whom this     */
/*  work is furnished to do so, subject to the following conditions: */
/*                                                                   */
/*    1. The source code must retain the above copyright notice,     */
/*       this list of conditions and the following disclaimer.       */
/*                                                                   */
/*    2. Any modifications to the source code must be clearly        */
/*       marked as such.                                             */
/*                                                                   */
/*    3. Redistributions in binary form must reproduce the above     */
/*       copyright notice, this list of conditions and the           */
/*       following disclaimer in the documentation and/or other      */
/*       materials provided with the distribution.  Otherwise, one   */
/*       must contact the HTS working group.                         */
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
/*         File: HGen.c: Generate parameter sequence from HMM        */
/*  ---------------------------------------------------------------  */


/* !HVER!HGen:   2.0 [NIT 06/11/06] */

#ifndef _HGEN_H_
#define _HGEN_H_

#ifdef __cplusplus
extern "C" {
#endif

#define MAXWINNUM 10   /* maximum number of static + deltas */

typedef enum {WLEFT=0, WRIGHT=1} WINWIDTH;
typedef enum {CHOLESKY=0, MIX=1, FB=2} ParmGenType;

typedef struct {
   char fn[MAXWINNUM][MAXFNAMELEN];  /* window coefficient file(s) */
   int num;                          /* number of windows */
   int maxw[2];                      /* max width [0(left) 1(right)] */
   int  **width;                     /* width [0..num-1][0(left) 1(right)] */
   float **coef;                     /* window coefficient [0..num][length[1]..length[2]] */
} Window;

typedef struct {
   char ext[MAXSTRLEN];  /* filename extension for this PdfStream */
   Boolean *ContSpace;   /* space indexes */
   Boolean fullCov;      /* full covariance flag */
   int vSize;            /* vector size of observation vector */
   int order;            /* vector size of static feature vector */
   int t;                /* time counter */
   int T;                /* number of frames */
   int width;            /* band width */
   Window win;           /* window coefficients */
   Matrix mseq;          /* sequence of mean vector */
   Covariance *vseq;     /* sequence of covariance matrices */
   Matrix C;             /* generated parameter c */
   DVector g;   
   DMatrix WUW;
   DVector WUM;
} PdfStream;

typedef struct {
   MemHeap *genMem;       /* generation memory heap */
   
   float speakRate;       /* speaking rate */
   float MSDthresh;       /* MSD threshold */
   Boolean modelAlign;    /* use model-level alignment from label */
   Boolean stateAlign;    /* use state-level alignment from label */
   HTime frameRate;       /* frame rate in 100ns */
   
   HMMSet *hset;          /* set of HMMs */
   HMMSet *dset;          /* set of duration models */
   int maxStates;         /* max # of states in hset */
   
   PdfStream pst[SMAX];   /* PdfStream for generation */
   int nPdfStream[SMAX];  /* # of PdfStreams and its size */ 
   
   Transcription *labseq; /* input label sequence */
   int labseqlen;         /* # of labels */
   Label **label;         /* labels sequence */
   
   HLink *hmm;            /* a sentence HMM for given label */
   HLink *dm;             /* a sentence duration models for given label */
   
   IMatrix sindex;        /* state sequence indexes */
   IMatrix durations;     /* state durations */
   int tframe;            /* total # of frames */
} GenInfo;

/* EXPORTED functions ------------------ */ 

void InitGen(void);
/*
   Initialise module
*/

void ResetGen(void);
/*
   Reset module
*/

void SetTraceGen(void);
/* 
   Set trace level 
*/

void InitialiseGenInfo (GenInfo *, Transcription *);
/*
 * Initialise GenInfo 
 */

void ResetGenInfo (GenInfo *);
/*
 * Reset GenInfo 
 */
 
void ParamGen (GenInfo *, UttInfo *, FBInfo *, ParmGenType);
/*
   Generate parameter sequence 
 */
   
#ifdef __cplusplus
}
#endif

#endif  /* _HGEN_H_ */

/* ------------------------- End of HGen.h --------------------------- */
