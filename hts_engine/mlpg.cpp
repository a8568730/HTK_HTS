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
/*    mlpg.c: speech parameter generation from pdf sequence          */
/*  ---------------------------------------------------------------  */

/* hts_engine libraries */
#include "misc.hpp"
#include "model.hpp"
#include "global.hpp"
#include "vocoder.hpp"
#include "mlpg.hpp"

double finv (const double x)
{
   if (x >=  INFTY2) return 0.0;
   if (x <= -INFTY2) return 0.0;
   if (x <= INVINF2 && x >= 0) return INFTY;
   if (x >= -INVINF2 && x < 0) return -INFTY;
   
   return(1.0/x);
}

/*------ parameter generation fuctions */
/* Calc_WUW_and_WUM: calcurate W'U^{-1}W and W'U^{-1}M */
void Calc_WUW_and_WUM (PStream *pst, const int m)
{
   int t, i, j, k;
   double WU;
   
   for (t=1; t<=pst->T; t++) {
      /* initialize */
      pst->sm.WUM[t] = 0.0;
      for (i=1; i<=pst->width; i++)  
         pst->sm.WUW[t][i] = 0.0;      
      
      /* calc WUW & WUM */
      for (i=0; i<pst->dw.num; i++)
         for (j=pst->dw.width[i][0]; j<=pst->dw.width[i][1]; j++)
            if ((t+j>0) && (t+j<=pst->T) && (pst->dw.coef[i][-j]!=0.0)) {
               WU = pst->dw.coef[i][-j]*pst->sm.ivseq[t+j][i*pst->order+m];
               pst->sm.WUM[t] += WU*pst->sm.mseq[t+j][i*pst->order+m];
               
               for (k=0; (k<pst->width) && (t+k<=pst->T); k++)
                  if ((k-j<=pst->dw.width[i][1]) && (pst->dw.coef[i][k-j]!=0.0))
                     pst->sm.WUW[t][k+1] += WU*pst->dw.coef[i][k-j];
            }
   }
   
   return;
}

/* LDL_Factorization: Factorize W'*U^{-1}*W to L*D*L' (L: lower triangular, D: diagonal) */
void LDL_Factorization (PStream *pst)
{
   int t,i,j;
    
   for (t=1; t<=pst->T; t++) {
      for (i=1; (i<pst->width) && (t-i>0); i++)
         pst->sm.WUW[t][1] -= pst->sm.WUW[t-i][i+1]*pst->sm.WUW[t-i][i+1]*pst->sm.WUW[t-i][1];
    
      for (i=2; i<=pst->width; i++) {
         for (j=1; (i+j<=pst->width) && (t-j>0); j++)
            pst->sm.WUW[t][i] -= pst->sm.WUW[t-j][j+1]*pst->sm.WUW[t-j][i+j]*pst->sm.WUW[t-j][1];
         pst->sm.WUW[t][i] /= pst->sm.WUW[t][1];
      }
   }
     
   return;
}

/* Forward_Substitution */
void Forward_Substitution (PStream *pst)
{
   int t,i;
   
   for (t=1; t<=pst->T; t++) {
      pst->sm.g[t] = pst->sm.WUM[t];
      for (i=1; (i<pst->width) && (t-i>0); i++)
         pst->sm.g[t] -= pst->sm.WUW[t-i][i+1]*pst->sm.g[t-i];
   }
   
   return;
}

/* Backward_Substitution */
void Backward_Substitution (PStream *pst, const int m)
{
   int t,i;

   for (t=pst->T; t>0; t--) {
      pst->par[t][m] = pst->sm.g[t]/pst->sm.WUW[t][1];
      for (i=1; (i<pst->width) && (t+i<=pst->T); i++)
         pst->par[t][m] -= pst->sm.WUW[t][i+1]*pst->par[t+i][m];
   }
   
   return;
}

/* mlpg: generate sequence of speech parameter vector maximizing its output probability for given pdf sequence */
void mlpg (PStream *pst)
{
   int m;
   const int M = pst->order;

   for (m=1; m<=M; m++) {
      Calc_WUW_and_WUM(pst,m);
      LDL_Factorization(pst);       /* LDL factorization */
      Forward_Substitution(pst);    /* forward substitution   */
      Backward_Substitution(pst,m); /* backward substitution  */
   }
}

/* ReadDouble: read one double */
double ReadDouble (FILE *fp)
{
   double d=0.0;

   if (fscanf(fp,"%lf",&d) != 1) {
      return 0.0;
   }

   return d;
}

/* ReadInt: read one integer from file */
int ReadInt (FILE *fp) 
{
   int i=0;
  
   if (fscanf(fp,"%d",&i) != 1){
      return 0;
   }

   return i;
}

/* InitDWin: Initialise window coefficients */
void InitDWin (PStream *pst)
{   
   int i,j;
   int fsize, leng;
   FILE *fp;

   /* memory allocation */
   pst->dw.width = (int **) HTS_Calloc(pst->dw.num, sizeof(int *));
   for (i=0; i<pst->dw.num; i++)
      pst->dw.width[i] = (int *) HTS_Calloc(2, sizeof(int));   
   pst->dw.coef = (double **) HTS_Calloc(pst->dw.num, sizeof(double *));
   
   /* set window coefficients */
   for (i=0; i<pst->dw.num; i++) {
      fp = HTS_Getfp(pst->dw.fn[i], "rt");      

      /* check the number of coefficients */
      fscanf(fp,"%d",&fsize);
      if (fsize<1)
         HTS_Error(1, "InitDWIn: number of coefficients in %s is invalid", pst->dw.fn[i]);

      /* read coefficients */
      pst->dw.coef[i] = (double *) HTS_Calloc(fsize, sizeof(double));

      for (j=0;j<fsize;j++) {
         fscanf(fp, "%lf", &(pst->dw.coef[i][j]));
      }

      fclose(fp);

      /* set pointer */
      leng = fsize / 2;
      pst->dw.coef[i] += leng;
      pst->dw.width[i][WLEFT] = -leng;
      pst->dw.width[i][WRIGHT] = leng;
         
      if (fsize % 2 == 0)
         pst->dw.width[i][WRIGHT]--;
   }
 
   pst->dw.maxw[WLEFT] = pst->dw.maxw[WRIGHT] = 0;
      
   for (i=0; i<pst->dw.num; i++) {
      if (pst->dw.maxw[WLEFT] > pst->dw.width[i][WLEFT])
         pst->dw.maxw[WLEFT] = pst->dw.width[i][WLEFT];
      if (pst->dw.maxw[WRIGHT] < pst->dw.width[i][WRIGHT])
         pst->dw.maxw[WRIGHT] = pst->dw.width[i][WRIGHT];
   }

   /* calcurate max_L to determine size of band matrix */
   if ( pst->dw.maxw[WLEFT] >= pst->dw.maxw[WRIGHT] )
      pst->dw.max_L = pst->dw.maxw[WLEFT];
   else
      pst->dw.max_L = pst->dw.maxw[WRIGHT];

   return;
}

/* FreeDWin: free regression window */
void FreeDWin (PStream *pst)
{   
   int i;
   
   /* free window */
   for (i=pst->dw.num-1; i>=0; i--) {
      pst->dw.coef[i] += pst->dw.width[i][WLEFT];
      HTS_Free(pst->dw.coef[i]);
   }
   HTS_Free(pst->dw.coef);
   
   for (i=pst->dw.num-1; i>=0; i--)
      HTS_Free(pst->dw.width[i]);
   HTS_Free(pst->dw.width);
      
   return;
}

/* InitPStream: Initialise PStream for parameter generation */
void InitPStream (PStream *pst)
{
   pst->width    = pst->dw.max_L*2+1;  /* band width of R */

   pst->sm.mseq  = HTS_AllocMatrix(pst->T, pst->vSize);
   pst->sm.ivseq = HTS_AllocMatrix(pst->T, pst->vSize);
   pst->sm.WUW   = HTS_AllocMatrix(pst->T, pst->width);
   pst->par      = HTS_AllocMatrix(pst->T, pst->order);
   
   pst->sm.g     = HTS_AllocVector(pst->T);
   pst->sm.WUM   = HTS_AllocVector(pst->T);
   
   return;
}

/* FreePStream: Free PStream */
void FreePStream (PStream *pst)
{
   HTS_FreeVector(pst->sm.WUM);
   HTS_FreeVector(pst->sm.g);
   
   HTS_FreeMatrix(pst->par,      pst->T);
   HTS_FreeMatrix(pst->sm.WUW,   pst->T);
   HTS_FreeMatrix(pst->sm.ivseq, pst->T);
   HTS_FreeMatrix(pst->sm.mseq,  pst->T);
   
   return;
}

/* pdf2speech: parameter generation and waveform synthesis */
void pdf2speech (FILE *rawfp, FILE *lf0fp, FILE *mcepfp, 
                 PStream *mceppst, PStream *lf0pst, 
                 globalP *gp, ModelSet *ms, UttModel *um, VocoderSetup *vs)
{
   int frame, mcepframe, lf0frame;
   int state, lw, rw, k, n;
   Model *m;
   HTS_Boolean nobound, *voiced;
   double f0;
   float temp;

   lf0pst->vSize  = ms->lf0stream;
   lf0pst->order  = 1;
   mceppst->vSize = ms->mcepvsize;
   mceppst->order = mceppst->vSize / mceppst->dw.num;

   InitDWin(mceppst);
   InitDWin(lf0pst);

   mcepframe = 0;
   lf0frame  = 0;
 
   /* voiced/unvoiced decision */
   voiced = (HTS_Boolean *) HTS_Calloc(um->totalframe, sizeof(HTS_Boolean));
   voiced--;
   
   for (m=um->mhead; m!=um->mtail ; m=m->next) {
      for (state=2; state<=ms->nstate+1; state++) {
         for (frame=1; frame<=m->dur[state]; frame++) {
            voiced[++mcepframe] = m->voiced[state];
            if (m->voiced[state]) {
               ++lf0frame;
            }
         }
      }
   }
   
   /* set the number of frames for mcep and lf0 */
   mceppst->T = mcepframe;
   lf0pst->T  = lf0frame;
  
   /* initialise parameter generation */
   InitPStream(mceppst);      
   InitPStream(lf0pst);

   mcepframe = 1;
   lf0frame  = 1;

   /* copy pdfs */
   for (m=um->mhead; m!=um->mtail; m=m->next) {
      for (state=2; state<=ms->nstate+1; state++) {
         for (frame=1; frame<=m->dur[state]; frame++) {
            /* copy pdf for mcep */
            for (k=0; k<ms->mcepvsize; k++) {
               mceppst->sm.mseq[mcepframe][k+1]  = m->mcepmean[state][k];
               mceppst->sm.ivseq[mcepframe][k+1] = finv(m->mcepvariance[state][k]);
            }
            /* copy pdfs for lf0 */ 
            for (k=0; k<ms->lf0stream; k++) {
               lw = lf0pst->dw.width[k][WLEFT];
               rw = lf0pst->dw.width[k][WRIGHT];
               nobound = (HTS_Boolean)1;
               /* check current frame is voiced/unvoiced boundary or not */
               for (n=lw; n<=rw;n++)
                  if (mcepframe+n<=0 || um->totalframe<mcepframe+n)
                     nobound = (HTS_Boolean) 0;
                  else
                     nobound = (HTS_Boolean) ((int)nobound & (int)voiced[mcepframe+n]);
                  
               /* copy pdfs */
               if (voiced[mcepframe]) {
                  lf0pst->sm.mseq[lf0frame][k+1] = m->lf0mean[state][k+1];
                  if (nobound || k==0) 
                     lf0pst->sm.ivseq[lf0frame][k+1] = finv(m->lf0variance[state][k+1]);
                  else   /* the variances for dynamic feature are set to inf on v/uv boundary */
                     lf0pst->sm.ivseq[lf0frame][k+1] = 0.0;
               }
            }
            if (voiced[mcepframe])
               lf0frame++;
            mcepframe++;
         }
      }
   }

   /* parameter generation for mcep */
   mlpg(mceppst);
   
   /* parameter generation for lf0 */
   if (lf0frame>0)
      mlpg(lf0pst);

   lf0frame = 1;
   
   /* synthesize speech waveforms */
   for (mcepframe=1; mcepframe<=mceppst->T; mcepframe++) {
      /* f0 modification */
      if (voiced[mcepframe])
         f0 = gp->F0_STD*exp(lf0pst->par[lf0frame++][1])+gp->F0_MEAN;  
      else                  
         f0 = 0.0;

      /* synthesize waveforms by MLSA filter */
      if (rawfp!=NULL)
         Vocoder(f0, &mceppst->par[mcepframe][1], mceppst->order-1, rawfp, gp, vs);
         
      /* output mcep sequence */
      if (mcepfp != NULL) {
         for (k=1;k<=mceppst->order;k++) {
            temp = (float)mceppst->par[mcepframe][k];
            fwrite(&temp, sizeof(float), 1, mcepfp);
         }
      }
      
      /* output f0 sequence */
      if (lf0fp != NULL) {
         temp = (float)f0;
         fwrite(&temp, sizeof(float), 1, lf0fp);
      }
   }
   
   /* close files */
   if (mcepfp != NULL) fclose(mcepfp);
   if (lf0fp  != NULL) fclose(lf0fp);
   if (rawfp  != NULL) fclose(rawfp);

   /* free memory */
   FreePStream(lf0pst);
   FreePStream(mceppst);
   
   HTS_Free(++voiced);
   
   FreeDWin(lf0pst);
   FreeDWin(mceppst);
   
   return;
}

/* -------------------- End of "mlpg.cc" -------------------- */
