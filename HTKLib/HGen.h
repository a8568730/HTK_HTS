/* ---------------------------------------------------------------- */
/*                                                                  */
/*     The HMM-Based Speech Synthesis System (HTS): version 1.0     */
/* 		      HTS Working Group                                     */
/*                                                                  */
/* 		 Department of Computer Science                             */
/* 		 Nagoya Institute of Technology                             */
/* 			      and                                               */
/*   Interdisciplinary Graduate School of Science and Engineering   */
/* 		 Tokyo Institute of Technology                              */
/* 		    Copyright (c) 2001-2002                                 */
/* 		      All Rights Reserved.                                  */
/*                                                                  */
/* Permission is hereby granted, free of charge, to use and         */
/* distribute this software in the form of patch code to HTK and    */
/* its documentation without restriction, including without         */
/* limitation the rights to use, copy, modify, merge, publish,      */
/* distribute, sublicense, and/or sell copies of this work, and to  */
/* permit persons to whom this work is furnished to do so, subject  */
/* to the following conditions:                                     */
/*                                                                  */
/*   1. Once you apply the HTS patch to HTK, you must obey the      */
/*      license of HTK.                                             */
/*                                                                  */
/*   2. The code must retain the above copyright notice, this list  */
/*      of conditions and the following disclaimer.                 */
/*                                                                  */
/*   3. Any modifications must be clearly marked as such.           */
/*                                                                  */
/* NAGOYA INSTITUTE OF TECHNOLOGY, TOKYO INSTITUTE OF TECHNOLOGY,   */
/* HTS WORKING GROUP, AND THE CONTRIBUTORS TO THIS WORK DISCLAIM    */
/* ALL WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL       */
/* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT   */
/* SHALL NAGOYA INSTITUTE OF TECHNOLOGY, TOKYO INSTITUTE OF         */
/* TECHNOLOGY, SPTK WORKING GROUP, NOR THE CONTRIBUTORS BE LIABLE   */
/* FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY        */
/* DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,  */
/* WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTUOUS   */
/* ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR          */
/* PERFORMANCE OF THIS SOFTWARE.                                    */
/*                                                                  */
/* ---------------------------------------------------------------- */ 
/*  HGen.h : parameter generation from pdf sequence based on        */
/*           Maximum Likelihood criterion with dynamic feature      */
/*           window constraints                                     */
/*                                                                  */
/*                                   2002/12/25 by Heiga Zen        */
/* ---------------------------------------------------------------- */

#define	WLEFT    0
#define	WRIGHT   1
#define	MAXDNUM 20   /* maximum number of static + deltas */
#define INFTY   (1.0E10)
#define INVINF  (1.0E-10)

typedef struct {
   int num;            /* number of static + deltas */
   char	*fn[MAXDNUM];  /* delta window coefficient file */
   int   **width;      /* width [0..num-1][0(left) 1(right)] */
   double **coef;      /* coefficient [1..num][length[1]..length[2]] */
   int maxw[2];	       /* max width [0(left) 1(right)] */
   int max_L;
} DWin;

typedef struct {
   SVector *mseq;     /* sequence of mean vector */
   SVector *ivseq;    /* sequence of invarsed variance vector */
   Matrix C;          /* generated parameter c */
   DVector g;			
   DMatrix WUW;
   DVector WUM;
} SMatrices;

typedef struct {
   int vSize;
   int order;
   int T;
   int width;
   DWin dw;
   SMatrices sm;
} PStream;

void InitPStream(MemHeap *x, PStream *);
void pdf2par(PStream *);
void FreePStream(MemHeap *x, PStream *);

