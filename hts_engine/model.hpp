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
/*    model.h: model definition                                      */
/*  ---------------------------------------------------------------  */

/* Model: structure for individual HMM */
typedef struct _Model {    
   char *name;             /* the name of this HMM */
   int durpdf;             /* duration pdf index for this HMM */
   int *lf0pdf;            /* mel-cepstrum pdf indexes for each state of this HMM */  
   int *mceppdf;           /* log f0 pdf indexes for each state of this HMM */
   int *dur;               /* duration for each state of this HMM */
   int totaldur;           /* total duration of this HMM */
   double **lf0mean;       /* mean vector of log f0 pdfs for each state of this HMM */
   double **lf0variance;   /* variance (diag) elements of log f0 for each state of this HMM */
   double **mcepmean;      /* mean vector of mel-cepstrum pdfs for each state of this HMM */
   double **mcepvariance;  /* variance (diag) elements of mel-cepstrum for each state of this HMM */
   HTS_Boolean *voiced;    /* voiced/unvoiced decision for each state of this HMM */
   struct _Model *next;    /* pointer to next HMM */
} Model; 

/* UttMode: structure for utterance HMM */
typedef struct _UttModel {
   Model *mhead;
   Model *mtail;
   int nModel;     /* # of models for current utterance     */
   int nState;     /* # of HMM states for current utterance */ 
   int totalframe; /* # of frames for current utterance     */
} UttModel;

/* ModelSet: structure for HMM set */
typedef struct _ModelSet {
   int nstate;               /* # of HMM states for individual HMM */
   int lf0stream;            /* # of stream for log f0 modeling */
   int mcepvsize;            /* vector size for mcep modeling */
   int *nlf0pdf;             /* # of pdfs for each state position (log F0) */
   int *nmceppdf;            /* # of pdfs for each state position (mcep) */
   int ndurpdf;              /* # of pdfs for duration */
   double **durpdf;          /* pdfs for duration */
   double ***mceppdf;        /* pdfs for mcep     */
   double ****lf0pdf;        /* pdfs for lf0      */
   FILE *fp[HTS_NUMMTYPE];   /* file pointer for mcep, logF0, and duration model */
} ModelSet;

/* Function prototypes for tree handling */
void LoadModelSet (ModelSet *);
void ClearModelSet (ModelSet *);

void FindDurPDF (Model *, ModelSet *, const double, double *);
void FindLF0PDF (const int, Model *, ModelSet *, const double);
void FindMcpPDF (const int, Model *, ModelSet *);

void InitModelSet (ModelSet *);

/* -------------------- End of "model.h" -------------------- */
