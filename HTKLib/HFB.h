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
/*         File: HFB.h: Forward Backward routines module       */
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
/*                     Copyright (c) 2001-2007                       */
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

/* !HVER!HFB:   3.3 [CUED 28/04/05] */

#ifndef _HFB_H_
#define _HFB_H_

#ifdef __cplusplus
extern "C" {
#endif

#define NOPRUNE 1.0E20

/* structure for the utterance information */
typedef struct {

  MemHeap transStack; /* utterance transcript information heap */
  MemHeap dataStack;  /* utterance data information heap */
  MemHeap dataStack2; /* utterance data2 information heap */

  int Q;              /* number of models in transcription */
  Transcription *tr;  /* current transcription */

  Boolean twoDataFiles; /* Using two data files */
   
  int S;              /* number of data streams */
  int T;              /* number of frames in utterance */
  ParmBuf pbuf;       /* parameter buffer */
  ParmBuf pbuf2;      /* a second parameter buffer (if required) */

   HTime tgtSampRate;  /* frame rate */
   HTime tgtSampRate2; /* second frame rate */
    
   Observation *o;     /* Observations */
   Observation *o2;    /* Cepstral Mean Normalised obervation, used in
                               single pass re-training */

  LogDouble pr;        /* log prob of current utterance */

} UttInfo;


/* structure for the forward-backward  pruning information */
typedef struct {
 
  short *qHi;               /* array[1..T] of top of pruning beam */
  short *qLo;               /* array[1..T] of bottom of pruning beam */
  int maxBeamWidth;         /* max width of beam in model units */
  LogDouble maxAlphaBeta;   /* max alpha/beta product along beam ridge */
  LogDouble minAlphaBeta;   /* min alpha/beta product along beam ridge */
  LogDouble pruneThresh;    /* pruning threshold currently */

} PruneInfo;

/* structure for the forward-backward alpha-beta structures */
typedef struct {
  
  MemHeap abMem;      /* alpha beta memory heap */
  PruneInfo *pInfo;   /* pruning information */
  HLink *up_qList;    /* array[1..Q] of active HMM defs */
  HLink *al_qList;    /* array[1..Q] of active align HMM defs */
   HLink *up_dList;    /* array[1..Q] of active dur model defs */
   HLink *al_dList;    /* array[1..Q] of active align dur model defs */
   MLink *qLink;       /* array[1..Q] of link to active HMM defs */
  LabId  *qIds;       /* array[1..Q] of logical HMM names (in qList) */
  short *qDms;        /* array[1..Q] of minimum model duration */
  DVector *alphat;    /* array[1..Q][1..Nq] of prob */
  DVector *alphat1;   /* alpha[t-1] */
   DVector **alpha;    /* array[1..T][1..Q][1..Nq] of prob */
  DVector **beta;     /* array[1..T][1..Q][1..Nq] of prob */
   float *****otprob;  /* array[1..T][1..Q][2..Nq-1][0..S] of prob */
  LogDouble pr;       /* log prob of current utterance */
  Vector occt;        /* occ probs for current time t */
   Vector **occa;      /* array[1..T][1..Q][1..Nq] of occ probs (trace only) */
   Vector ****occm;    /* array[1..T][1..Q][1..Nq][1..S][1..M] of occ probs (param gen only) */

} AlphaBeta;

/* structure storing the model set and a pointer to it's alpha-beta pass structure */
typedef struct {
  Boolean twoModels;  /* Enable two model reestimation */
   Boolean useAlign;   /* Using model alignment */
  HMMSet *up_hset;    /* set of HMMs to be re-estimated */
  HMMSet *al_hset;    /* HMMs to use for alignment */
                      /* these are equal unless 2 model reest */
   HMMSet *up_dset;    /* set of duration models to be estimated */
   HMMSet *al_dset;    /* duration models to use for alignment */
                       /* these are equal unless 2 model reest */
  HSetKind hsKind;    /* kind of the alignment HMM system */
   UPDSet uFlags;      /* parameter update flags for HMMs */
   UPDSet dur_uFlags;  /* parameter update flags for duration models */
  int skipstart;      /* Skipover region - debugging only */
  int skipend;
  int maxM;           /* maximum number of mixtures in hmmset */
  int maxMixInS[SMAX];/* array[1..swidth[0]] of max mixes */
  AlphaBeta *ab;      /* Alpha-beta structure for this model */
   
   AdaptXForm *inXForm;        /* current input transform for HMMs (if any) */
  AdaptXForm *al_inXForm;/* current input transform for al_hset (if any) */
   AdaptXForm *paXForm;        /* current parent transform for HMMs (if any) */
   AdaptXForm *dur_inXForm;    /* current input transform for duration models (if any) */
   AdaptXForm *dur_al_inXForm; /* current input transform for al_dset (if any) */
   AdaptXForm *dur_paXForm;    /* current parent transform for duration models (if any) */
} FBInfo;

/* EXPORTED FUNCTIONS-------------------------------------------------*/

/* Initialise HFB module */
void InitFB(void) ;

/* Reset HFB module */
void ResetFB(void);

/* Allow tools to enable top-level tracing in HFB. Only here for historical reasons */
void SetTraceFB(void);

/* Initialise the forward backward memory stacks etc */
void InitialiseForBack(FBInfo *fbInfo, MemHeap *x, HMMSet *hset, UPDSet uFlags, HMMSet *dset, UPDSet dur_uFlags,
                       LogDouble pruneInit, LogDouble pruneInc, LogDouble pruneLim, 
                       float minFrwdP, Boolean useAlign);

/* Use a different model set for alignment */
void UseAlignHMMSet(FBInfo* fbInfo, MemHeap* x, HMMSet *al_hset, HMMSet *al_dset);

/* Initialise the utterance Information */
void InitUttInfo(UttInfo *utt, Boolean twoFiles );

/* GetInputObs: Get input Observations for t */
void GetInputObs( UttInfo *utt, int t, HSetKind hsKind );

/* load the labels into the UttInfo structure from file */
void LoadLabs(UttInfo *utt, FileFormat lff, char * datafn,
	      char *labDir, char *labExt);

/* load the data file(s) into the UttInfo structure */
void LoadData(HMMSet *hset, UttInfo *utt, FileFormat dff, 
	      char * datafn, char * datafn2);

/* Initialise the observation structures within UttInfo */
void InitUttObservations(UttInfo *utt, HMMSet *hset,
			 char * datafn, int * maxmixInS);

/* reset the observation structures within UttInfo */ 
void ResetUttObservations (UttInfo *utt, HMMSet *hset);

/* FBUtt: apply forward-backward to given utterance */
Boolean FBUtt(FBInfo *fbInfo, UttInfo *utt);

/* PrLog: print a log value */
void PrLog(LogDouble x);

#ifdef __cplusplus
}
#endif

#endif  /* _HFB_H_ */
