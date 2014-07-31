/*  ---------------------------------------------------------------  */
/*      The HMM-Based Speech Synthesis System (HTS): version 1.1,1   */
/*                        HTS Working Group                          */
/*                                                                   */
/*                   Department of Computer Science                  */
/*                   Nagoya Institute of Technology                  */
/*                                and                                */
/*    Interdisciplinary Graduate School of Science and Engineering   */
/*                   Tokyo Institute of Technology                   */
/*                      Copyright (c) 2001-2003                      */
/*                        All Rights Reserved.                       */
/*                                                                   */
/*  Permission is hereby granted, free of charge, to use and         */
/*  distribute this software and its documentation without           */
/*  restriction, including without limitation the rights to use,     */
/*  copy, modify, merge, publish, distribute, sublicense, and/or     */
/*  sell copies of this work, and to permit persons to whom this     */
/*  work is furnished to do so, subject to the following conditions: */
/*                                                                   */
/*    1. The code must retain the above copyright notice, this list  */
/*       of conditions and the following disclaimer.                 */
/*                                                                   */
/*    2. Any modifications must be clearly marked as such.           */
/*                                                                   */
/*  NAGOYA INSTITUTE OF TECHNOLOGY, TOKYO INSITITUTE OF TECHNOLOGY,  */
/*  HTS WORKING GROUP, AND THE CONTRIBUTORS TO THIS WORK DISCLAIM    */
/*  ALL WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL       */
/*  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT   */
/*  SHALL NAGOYA INSTITUTE OF TECHNOLOGY, TOKYO INSITITUTE OF        */
/*  TECHNOLOGY, HTS WORKING GROUP, NOR THE CONTRIBUTORS BE LIABLE    */
/*  FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY        */
/*  DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,  */
/*  WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTUOUS   */
/*  ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR          */
/*  PERFORMANCE OF THIS SOFTWARE.                                    */
/*                                                                   */
/*  ---------------------------------------------------------------  */
/*     HGen.c : parameter generation from pdf sequence based on      */
/*              Maximum Likelihood criterion with dynamic feature    */
/*              window constraints                                   */
/*                                                                   */
/*                                   2003/12/26 by Heiga Zen         */
/* ----------------------------------------------------------------  */

char *hgen_version = "!HVER!HGen:   1.1.1 [NIT 26/12/03]";
char *hgen_vc_id = "$Id: HGen.c,v 1.1.1 2003/12/26 15:40:13 zen Exp $";

#include "HShell.h"
#include "HMem.h"
#include "HMath.h"
#include "HGen.h"
#include "SPTK.h"

/* InitDWin: Initialise window coefficients used for calcuration of delta parameter */
static void InitDWin(MemHeap *x, PStream *pst)
{
   int i;
   int fsize, leng;
   FILE *fp;

   /* memory allocation */
   if ((pst->dw.width = (int **) New(x,pst->dw.num*sizeof(int *))) == NULL) {
      HError(9905, "InitDWin: Cannot Allocate Memory\n");
      exit(1);
   }
   for (i=0; i<pst->dw.num; i++)
      if ( (pst->dw.width[i] = (int *) New(x,2*sizeof(int))) == NULL) {
         HError(9905, "InitDWin: Cannot Allocate Memory\n");
         exit(1);
      }
   if ((pst->dw.coef=(double **)New(x,pst->dw.num*sizeof(double *))) == NULL) {
      HError(9905, "InitDWin: Cannot Allocate Memory\n");
      exit(1);
   }
  
   /* window for static parameter */
   pst->dw.width[0][WLEFT] = pst->dw.width[0][WRIGHT] = 0;
   pst->dw.coef[0] = (double *) New(x, sizeof(double));
   pst->dw.coef[0][0] = 1;

   /* read delta coefficients */
   for (i=1;i<pst->dw.num;i++) {
      if ((fp = fopen(pst->dw.fn[i], "r")) == NULL) {
         HError(9901, "InitDWin: Cannot open file  <%s>\n", pst->dw.fn[i]);
         exit(1);
      }

      /* check the number of coefficients */
      fseek(fp, 0L, 2);
      fsize = ftell(fp) / sizeof(float);
      fseek(fp, 0L, 0);

      /* read coefficients */
      pst->dw.coef[i] = (double *) New(x,fsize*sizeof(double));
      freadf(pst->dw.coef[i], sizeof(**(pst->dw.coef)), fsize, fp);

      /* set pointer */
      leng = fsize / 2;
      pst->dw.coef[i] += leng;
      pst->dw.width[i][WLEFT] = -leng;
      pst->dw.width[i][WRIGHT] = leng;
      if (fsize % 2 == 0)
         pst->dw.width[i][WRIGHT]--;
      
      fclose(fp);
   }

   pst->dw.maxw[WLEFT] = pst->dw.maxw[WRIGHT] = 0;
   for (i = 0; i < pst->dw.num; i++) {
      if (pst->dw.maxw[WLEFT] > pst->dw.width[i][WLEFT])
         pst->dw.maxw[WLEFT] = pst->dw.width[i][WLEFT];
      if (pst->dw.maxw[WRIGHT] < pst->dw.width[i][WRIGHT])
         pst->dw.maxw[WRIGHT] = pst->dw.width[i][WRIGHT];
   }

   /* calcurate max_L to determine size of band matrix */
   if( pst->dw.maxw[WLEFT] >= pst->dw.maxw[WRIGHT] )
      pst->dw.max_L = pst->dw.maxw[WLEFT];
   else
      pst->dw.max_L = pst->dw.maxw[WRIGHT];
  
   return;
}

/* EXPORT->InitPStream: Initialize PStream which contain all parameter for ML parameter generation */
void InitPStream(MemHeap *x, PStream *pst)
{
   InitDWin(x,pst);
   
   pst->vSize = pst->order*pst->dw.num;        /* size of vector */
   pst->width = pst->dw.max_L*2+1;             /* band width of WUW */

   pst->sm.mseq  = (SVector *)New(x,pst->T*sizeof(SVector));  /* mean vector sequence */ 
   pst->sm.ivseq = (SVector *)New(x,pst->T*sizeof(SVector));  /* inversed covariance(diagonal) sequeance */
   pst->sm.C     = CreateMatrix(x,pst->T,pst->order);   /* generated parameter */
   pst->sm.g   = CreateDVector(x,pst->T);               /* temporary */
   pst->sm.WUM = CreateDVector(x,pst->T);               /* W^{T}U^{-1}W */
   pst->sm.WUW = CreateDMatrix(x,pst->T,pst->width);    /* W^{T}U^{-1}M */

   pst->sm.mseq--; pst->sm.ivseq--;
   
   return;
}

/* calc_WUM_and_WUW: calcurate W'* U^-1 * M and W'* U^-1 * W */
static void calc_WUM_and_WUW(PStream *pst, const int m)
{
   register int t, i, j, k;
   double WU;
   
   for (t=1;t<=pst->T;t++) {
      pst->sm.WUM[t] = pst->sm.ivseq[t][m] * pst->sm.mseq[t][m];
      pst->sm.WUW[t][1] = pst->sm.ivseq[t][m]; 
      
      for (i=2;i<=pst->width;i++)  
         pst->sm.WUW[t][i]=0.0;      
      
      for (i=1;i<pst->dw.num;i++)
         for (j=pst->dw.width[i][0];j<=pst->dw.width[i][1];j++)
            if ((t+j>0) && (t+j<=pst->T) && (pst->dw.coef[i][-j]!=0.0)) {
               WU = pst->dw.coef[i][-j]*pst->sm.ivseq[t+j][i*pst->order+m];
               pst->sm.WUM[t] += WU*pst->sm.mseq[t+j][i*pst->order+m];
               
               for (k=0;k<pst->width;k++)
                  if ((k-j<=pst->dw.width[i][1]) && (t+k<=pst->T) && (pst->dw.coef[i][k-j]!=0.0))
                     pst->sm.WUW[t][k+1] += WU*pst->dw.coef[i][k-j];
            }
   }
   
   return;
}

/* Cholesky: Cholesky factorization for band matrix W'* U^-1 * W */
static void Cholesky(PStream *pst)
{
   int t,i,j;
   
   pst->sm.WUW[1][1] = sqrt(pst->sm.WUW[1][1]);

   for (i=2;i<=pst->width;i++)
      pst->sm.WUW[1][i] /= pst->sm.WUW[1][1];

   for (t=2;t<=pst->T;t++) {
      for (i=1;i<pst->width;i++)
         if(t-i>0)
            pst->sm.WUW[t][1] -= pst->sm.WUW[t-i][i+1]*pst->sm.WUW[t-i][i+1];

      pst->sm.WUW[t][1] = sqrt(pst->sm.WUW[t][1]);

      for (i=2;i<=pst->width;i++) {
         for (j=1;j<=pst->dw.max_L;j++)
            if (i<pst->width && t-j>0)
               pst->sm.WUW[t][i] -= pst->sm.WUW[t-j][i-j+1]*pst->sm.WUW[t-j][i+1];

         pst->sm.WUW[t][i] /= pst->sm.WUW[t][1];
      }
   }
   
   return;
}

/* Cholesky_forward : forward substitution */
static void Cholesky_forward(PStream *pst)
{
   int t,i;
   double hold;

   pst->sm.g[1] = pst->sm.WUM[1]/pst->sm.WUW[1][1];

   for (t=2;t<=pst->T;t++) {
      hold = 0.0;
      for (i=1;i<pst->width;i++)
         if (t-i>0)
            hold += pst->sm.WUW[t-i][i+1] * pst->sm.g[t-i];
      pst->sm.g[t] = (pst->sm.WUM[t]-hold)/pst->sm.WUW[t][1];
   }
   
   return;
}

/* Cholesky_backward : backward substitution */
static void Cholesky_backward(PStream *pst, const int m)
{
   int t,i;
   double hold;
   
   pst->sm.C[pst->T][m] = pst->sm.g[pst->T]/pst->sm.WUW[pst->T][1];

   for (t=pst->T-1;t>0;t--) {
      hold = 0.0;
      for (i=1;i<pst->width;i++)
         if (t+i<=pst->T)
            hold += pst->sm.WUW[t][i+1]*pst->sm.C[t+i][m];
      pst->sm.C[t][m] = (pst->sm.g[t]-hold)/pst->sm.WUW[t][1];
   }
   
   return;
}

/* EXPORT->pdf2par: ML parameter generation */
void pdf2par(PStream *pst)
{   
   int m;

   for (m=1;m<=pst->order;m++) {
      calc_WUM_and_WUW(pst,m);
      Cholesky(pst);            /* cholesky decomposition */
      Cholesky_forward(pst);    /* forward substitution */
      Cholesky_backward(pst,m); /* backward substitution */
   }
   
   return;
}

/* FreeDWin: Free delta window coefficients */
static void FreeDWin(MemHeap *x, PStream *pst)
{
   int i;
    
   /* free pst->dw.coef */
   for (i=pst->dw.num-1; i>=0; i--)
      Dispose(x,pst->dw.coef[i]);
   Dispose(x,pst->dw.coef);
   
   /* free pst->dw.width */
   for (i=pst->dw.num-1; i>=0; i--)
      Dispose(x,pst->dw.width[i]);
   Dispose(x,pst->dw.width);
   
   return;
}

/* EXPORT->FreePStream: Free PStream */
void FreePStream(MemHeap *x, PStream *pst)
{
   FreeDMatrix(x,pst->sm.WUW); 
   FreeDVector(x,pst->sm.WUM);
   FreeDVector(x,pst->sm.g);
   FreeMatrix(x,pst->sm.C);

   pst->sm.ivseq++;
   pst->sm.mseq++;

   Dispose(x,pst->sm.ivseq);
   Dispose(x,pst->sm.mseq);
   
   FreeDWin(x,pst);
   
   return;
}

/* ---------------------------- End of HGen.c --------------------------- */

