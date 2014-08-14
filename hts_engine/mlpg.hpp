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
/*    mlpg.h: speech parameter generation from pdf sequence          */
/*  ---------------------------------------------------------------  */

#define INFTY   ((double) 1.0e+38)
#define INFTY2  ((double) 1.0e+19)
#define INVINF  ((double) 1.0e-38)
#define INVINF2 ((double) 1.0e-19)
#define LTPI    1.83787706640935     /* log(2*PI) */

#define WLEFT  0
#define WRIGHT 1

/* DWin: structure for regression window */
typedef struct _DWin {
   int num;         /* number of static + deltas */
   char **fn;       /* delta window coefficient file */
   int **width;     /* width [0..num-1][0(left) 1(right)] */
   double **coef;   /* coefficient [0..num-1][length[0]..length[1]] */
   int maxw[2];     /* max width [0(left) 1(right)] */
   int max_L;       /* max {maxw[0], maxw[1]} */
} DWin;

/* SMatrices: structure for matrices to generate sequence of speech parameter vectors */
typedef struct _SMatrices {
   double **mseq;   /* sequence of mean vector */
   double **ivseq;  /* sequence of invarsed variance vector */
   double *g;       /* for forward substitution */
   double **WUW;    /* W' U^-1 W  */
   double *WUM;     /* W' U^-1 mu */
} SMatrices;

/* PStream: structure for parameter generation setting */
typedef struct _PStream {
   int vSize;       /* vector size of observation vector (include static and dynamic features) */
   int order;       /* vector size of static features */
   int T;           /* length */
   int width;       /* width of dynamic window */
   DWin dw;         /* dynamic window */
   double **par;     /* output parameter vector */
   SMatrices sm;    /* matrices for parameter generation */
} PStream;

/* Function prototype for HMM-based speech synthesis */
void pdf2speech (FILE *, FILE *, FILE *, 
                 PStream *, PStream *,  
                 globalP *, ModelSet *, UttModel *, VocoderSetup *);

/* -------------------- End of "mlpg.h" -------------------- */
