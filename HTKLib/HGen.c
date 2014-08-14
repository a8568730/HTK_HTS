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
/*         File: HGen.c: Generate parameter sequence from HMM        */
/*  ---------------------------------------------------------------  */

char *hgen_version = "!HVER!HGen:   2.0.1 [NIT 01/10/07]";
char *hgen_vc_id = "$Id: HGen.c,v 1.23 2007/09/17 11:26:41 zen Exp $";

#include "HShell.h"     /* HMM ToolKit Modules */
#include "HMem.h"
#include "HMath.h"
#include "HSigP.h"
#include "HAudio.h"
#include "HWave.h"
#include "HVQ.h"
#include "HParm.h"
#include "HLabel.h"
#include "HModel.h"
#include "HTrain.h"
#include "HUtil.h"
#include "HAdapt.h"
#include "HFB.h"
#include "HGen.h"

/* ------------------- Trace Information ------------------------------ */
/* Trace Flags */
#define T_TOP   0001  /* Top level tracing */
#define T_MAT   0002  /* trace matrices */
#define T_STA   0004  /* trace state sequence */

static int trace = 0;

static ConfParam *cParm[MAXGLOBS];  /* config parameters */
static int nParm = 0;

static int maxIter   = 20;          /* max iterations in EM-based parameter generation */
static float epsilon = 1.0E-4;      /* convergence criterion */
static Boolean rndpg = FALSE;       /* random generation instead of ML one */
static float rndMean = 0.0;         /* mean of Gaussian noise for random generation */
static float rndVar  = 1.0;         /* variance of Gaussian noise for random generation */

/* -------------------------- Initialisation ----------------------- */

/* EXPORT->InitGen: initialise module */
void InitGen(void)
{
   int m;
   int i;
   double d;
   Boolean b;

   Register(hgen_version,hgen_vc_id);

   for (m=0; m<2; m++) {
      nParm = GetConfig("HGEN", TRUE, cParm, MAXGLOBS);
      if (nParm>0){
         if (GetConfInt (cParm,nParm,"TRACE",     &i)) trace   = i;
         if (GetConfInt (cParm,nParm,"MAXITER",   &i)) maxIter = i;
         if (GetConfFlt (cParm,nParm,"EPSILON",   &d)) epsilon = d;
         if (GetConfBool(cParm,nParm,"RNDPG",     &b)) rndpg   = b;
         if (GetConfFlt (cParm,nParm,"RNDMEAN",   &d)) rndMean = d;
         if (GetConfFlt (cParm,nParm,"RNDVAR",    &d)) rndVar  = d;
      }
   }
}

/* EXPORT->ResetGen: reset module */
void ResetGen(void)
{
   return;  /* do nothing */
}

/* EXPORT->SetConfGen: set trace flag for this module */
void SetTraceGen(void)
{
   trace |= T_TOP;
}

/* -------------------------- Window coefficients handling -------------------------- */

/* InitWindow: Initialise and load window coefficients */
static void InitWindow (MemHeap *x, PdfStream *pst)
{
   Source src;
   int i, winlen;

   /* Allocate array for window coefficients */
   pst->win.width = (int **) New(x, pst->win.num*sizeof(int *));
   for (i=0; i<pst->win.num; i++)
      pst->win.width[i] = (int *) New(x,2*sizeof(int));
   pst->win.coef = (float **)New(x, pst->win.num*sizeof(float *));

   /* Load window coefficients */
   for (i=0; i<pst->win.num; i++) {
      if (InitSource(pst->win.fn[i],&src,NoFilter)<SUCCESS)
         HError(2610,"InitWindow: Can't open file %s", pst->win.fn[i]);

      /* Read window length */
      ReadInt(&src, &winlen, 1, FALSE);

      /* Read window coefficients */
      pst->win.coef[i] = (float *) New(x,winlen*sizeof(float));
      /* for (j=0; j<winlen; j++)
         pst->win.coef[i][j] = 0.0; */
      ReadFloat(&src, pst->win.coef[i], winlen, FALSE);

      /* Set pointer */
      pst->win.coef[i] += winlen/2;
      pst->win.width[i][WLEFT]  = -winlen/2;
      pst->win.width[i][WRIGHT] =  winlen/2;
      if (winlen%2 == 0)
         pst->win.width[i][WRIGHT]--;

      CloseSource(&src);
   }

   /* Set pst->width */
   pst->win.maxw[WLEFT] = pst->win.maxw[WRIGHT] = 0;
   for (i=0; i<pst->win.num; i++) {
      if (pst->win.maxw[WLEFT] > pst->win.width[i][WLEFT])
         pst->win.maxw[WLEFT] = pst->win.width[i][WLEFT];
      if (pst->win.maxw[WRIGHT] < pst->win.width[i][WRIGHT])
         pst->win.maxw[WRIGHT] = pst->win.width[i][WRIGHT];
   }
   pst->width = (abs(pst->win.maxw[WLEFT]) > abs(pst->win.maxw[WRIGHT])) ? abs(pst->win.maxw[WLEFT]) : abs(pst->win.maxw[WRIGHT]);
   pst->width = 2*pst->width+1;

   return;
}

/* -------------------------- PdfStream handling -------------------------- */

/* CreatePdfStreams: create PdfStream in GenInfo */
static void CreatePdfStreams (GenInfo *genInfo)
{
   int p,t,M;
   Boolean fullCov=FALSE;
   PdfStream *pst;
   MemHeap *x = genInfo->genMem;

   switch(genInfo->hset->ckind) {
   case DIAGC:
   case INVDIAGC:
      fullCov = FALSE;
      break;
   case FULLC:
      fullCov = TRUE;
      break;
   default:
      HError(9999,"CreatePdfStreams: Currently only DIAGC, INVDIAGC and FULLC are supported.");
   }

   for (p=1; p<=genInfo->nPdfStream[0]; p++) {
      /* p-th PdfStream */
      pst = &(genInfo->pst[p]);
      
      /* Load window file for p-th PdfStream */
      InitWindow(x, pst);

      /* full or diag covariance */
      pst->fullCov = fullCov;

      /* prepare mean vector and inverse covariance matrix sequences */
      pst->mseq = CreateMatrix(x, pst->T, pst->vSize);               /* mean vectors */
      pst->vseq = (Covariance *) New(x, pst->T*sizeof(Covariance));  /* inverse covariance matrices */
      pst->vseq--;

      for (t=1; t<=pst->T; t++) {
         if (fullCov)
            pst->vseq[t].inv = CreateSTriMat(x, pst->vSize);  /* inverse covariance (precision) matrices */
         else
            pst->vseq[t].var = CreateVector (x, pst->vSize);  /* inverse variances */
      }

      /* prepare vector and matrices for a set of linear equations */
      M = (fullCov) ? pst->order : 1;
      pst->C    = CreateMatrix (x, pst->T, pst->order);     /* generated parameter sequence */
      pst->g    = CreateDVector(x, M*pst->T);               /* vector for forward and backward substitution */
      pst->WUM  = CreateDVector(x, M*pst->T);               /* W'*U^{-1}*W */
      pst->WUW  = CreateDMatrix(x, M*pst->T, M*pst->width); /* W'*U^{-1}*M */
   }

   return;
}

/* ChkBoundary: check k-th dimension at absolute_t-th frame is on the boundary where switching MSD space */
static Boolean ChkBoundary (PdfStream *pst, const int k, const int absolute_t, const int absolute_T)
{
   int j;

   const int t = pst->t;  /* time counter in this PdfStream */
   const int T = pst->T;  /* total frame in this PdfStream */
   const int i = (k-1) / pst->order;  /* window index for k-th dimension */

   /* check */
   for (j=pst->win.width[i][0]; j<=pst->win.width[i][1]; j++)
      if ( t+j<0 || t+j>T || (pst->win.coef[i][j]!=0.0 && absolute_t+j>=1 && absolute_t+j<=absolute_T && !pst->ContSpace[absolute_t+j]) )
         return TRUE;

   return FALSE;
}

/* SetupPdfStreams: setup PdfStreams for parameter generation */
static void SetupPdfStreams (GenInfo *genInfo)
{
   int i,j,k,l,d,t,p,s,stream,m,v,max;
   float weight;
   Boolean bound;
   PdfStream *pst;
   StateInfo *si;
   MixPDF *mpdf;
   Vector mseq_t;
   Covariance vseq_t;

   /* absolute time */
   const int T = genInfo->tframe;

   /* initialize time counter and statistics of each PdfStream */
   for (p=1; p<=genInfo->nPdfStream[0]; p++) {
      /* p-th PdfStream */
      pst = &(genInfo->pst[p]);

      /* initialize time counter, mseq, and vseq */
      pst->t = 1;
      for (t=1; t<=pst->T; t++) {
         ZeroVector(pst->mseq[t]);
         if (pst->fullCov)
            ZeroTriMat(pst->vseq[t].inv);
         else
            ZeroVector(pst->vseq[t].var);
      }
   }

   /* load mean and variance and set mseq and vseq */
   for (i=1,t=1; i<=genInfo->labseqlen; i++) {
      for (j=1; genInfo->sindex[i][j]!=0; j++) {
         si = genInfo->hmm[i]->svec[genInfo->sindex[i][j]].info;
         for (d=0; d<genInfo->durations[i][j]; d++,t++) {
            for (p=stream=1; p<=genInfo->nPdfStream[0]; stream+=genInfo->nPdfStream[p++]) {
               /* p-th PdfStream */
               pst = &(genInfo->pst[p]);
            
               /* continuous space or not */
               if (pst->ContSpace[t]) {
                  /* mean vector and covariance matrix at this frame */
                  mseq_t = pst->mseq[ pst->t ];
                  vseq_t = pst->vseq[ pst->t ];
                  
                  /* calculate mean_jp and cov_jp */
                  for (s=stream,v=1; s<stream+genInfo->nPdfStream[p]; v+=genInfo->hset->swidth[s++]) {
                     /* calculate total weight of continuous spaces */
                     for (m=1,max=0,weight=LZERO; m<=si->pdf[s].info->nMix; m++) {
                        mpdf = si->pdf[s].info->spdf.cpdf[m].mpdf;
                        if (VectorSize(mpdf->mean)==genInfo->hset->swidth[s]) {
                           /* use the best mixture only */
                           if (si->pdf[s].info->spdf.cpdf[m].weight >= weight) {
                              weight = si->pdf[s].info->spdf.cpdf[m].weight;
                              max = m;
                           }
                        }
                     }
                     if (max==0)
                        HError(999,"SetupPdfStreams: no mix found");
                  
                     /* set pdf streams */
                     mpdf = si->pdf[s].info->spdf.cpdf[max].mpdf;
                     for (k=1; k<=genInfo->hset->swidth[s]; k++) {
                        bound = ChkBoundary(pst, v+k-1, t, T);
                        switch (mpdf->ckind) {
                        case DIAGC:
                           mseq_t    [v+k-1] += (bound) ? 0.0 : 1.0 / mpdf->cov.var[k] * mpdf->mean[k];
                           vseq_t.var[v+k-1] += (bound) ? 0.0 : 1.0 / mpdf->cov.var[k];
                           break;
                        case INVDIAGC:
                           mseq_t    [v+k-1] += (bound) ? 0.0 : mpdf->cov.var[k] * mpdf->mean[k];
                           vseq_t.var[v+k-1] += (bound) ? 0.0 : mpdf->cov.var[k];
                           break;
                        case FULLC:
                           /* diagonal elements */
                           mseq_t    [v+k-1]        += (bound) ? 0.0 : mpdf->cov.inv[k][k] * mpdf->mean[k];
                           vseq_t.inv[v+k-1][v+k-1] += (bound) ? 0.0 : mpdf->cov.inv[k][k];

                           /* off-diagonal elements */
                           for (l=1; l<k; l++) {
                              mseq_t    [v+k-1]        += (bound) ? 0.0 : mpdf->cov.inv[k][l] * mpdf->mean[l];
                              mseq_t    [v+l-1]        += (bound) ? 0.0 : mpdf->cov.inv[k][l] * mpdf->mean[k];
                              vseq_t.inv[v+k-1][v+l-1] += (bound) ? 0.0 : mpdf->cov.inv[k][l];
                           }
                           break;
                        default:
                           HError(9999, "SetupPdfStreams: Only DIAGC, INVDIAGC and FULLC are supported.");
                        }
                     }
                  }
                  /* update counter */
                  pst->t++;
               }
            }
         }
      }
   }

   return;
}

/* ----------------------- Sentence model initization/reset routines ----------------------- */

/* GetStateIndex: get state index from name */
static int GetStateIndex (LabId id)
{
   return(atoi(strrchr(id->name,'[')+1));
}

/* SetStateSequence: set state sequence which maximize its state sequence prob */
static void SetStateSequence (GenInfo *genInfo)
{
   int i, j, k, n, best;
   LogFloat trans;
   HLink hmm;
   Label *label;
   IntSet acyclic;

   /* set state sequence */
   if (genInfo->stateAlign) {        /* state alignments are given */
      for (i=1,j=0,n=0; i<=CountLabs(genInfo->labseq->head); i++) {
         /* get i-th label */
         label = GetLabN(genInfo->labseq->head, i);
      
         /* prepare an array to store state durations in the j-th model */
         if (label->auxLab[1]!=NULL) {
            j++; n=0;
         }
      
         /* get state duration from the given label */
         k = GetStateIndex(label->labid);
         if (k<=1 || k>=genInfo->hmm[j]->numStates)
            HError(9999, "SetStateSequence: #state in the %d-th label is out of range", i);
         
         /* set the n-th state in this model */
         genInfo->sindex[j][++n] = k;
      }
   }
   else {      /* state alignment is not given */
      /* IntSet to detect acyclic graph */
      acyclic = CreateSet(genInfo->maxStates);

      if (trace&T_STA)
         printf(" State sequence:\n");
               
      /* get state durations from given label sequence */
      for (i=1; i<=genInfo->labseqlen; i++) {
         hmm = genInfo->hmm[i];
         ClearSet(acyclic);
         
         if (trace&T_STA)
            printf("  %d-th model: %d",i,1);
         
         /* trace most likely state sequence in this model */
         for (j=1,n=0; j!=hmm->numStates; j=best) {
            for (k=2,best=1,trans=LZERO; k<=hmm->numStates; k++) {
               if (k==j) continue;  /* exclude self-transition */
               if (hmm->transP[j][k]>trans) {
                  trans = hmm->transP[j][k];
                  best  = k;
               }
            }
            
            /* check acyclic transition */
            if (IsMember(acyclic,best))
               HError(9999, "SetStateSeuqnce: acyclic HMM transition is detected.");
            AddMember(acyclic,best);
            
            if (trace&T_STA)
               printf(" -> %3d (%1.2f)", best, L2F(trans));
            
            /* set n-th state to best */
            if (best!=hmm->numStates)
               genInfo->sindex[i][++n] = best;
         }
         
         if (trace&T_STA) {
            printf("\n");
            fflush(stdout);
         }
      }
      
      /* free acyclic */
      FreeSet(acyclic);
   }
   
   return;
}

/* CountDurStat: count duration statistics */
static void CountDurStat (Vector mean, Vector ivar, float *sum, float *sqr, IntVec sindex)
{
   int i,j;
   
   for (i=1; sindex[i]!=0; i++) {
      j = sindex[i] - 1;
      *sum += mean[j];
      *sqr += 1.0/ivar[j];
   }
   
   return;
}

/* SetStateDurations: set state durations */
static void SetStateDurations (GenInfo *genInfo)
{
   int i,j,k,s,cnt,nStates,modeldur,start=0,tframe=0;
   float sum,sqr,rho,diff=0.0;
   Vector *mean, *ivar;
   Label *label;
   HLink dm;
   
   /* state duration statistics storage */
   if ((mean = (Vector *) New(genInfo->genMem, genInfo->labseqlen*sizeof(Vector)))== NULL)
      HError(9905,"SetStateDurations: Cannot allocate memory for mean");
   if ((ivar = (Vector *) New(genInfo->genMem, genInfo->labseqlen*sizeof(Vector)))== NULL)
      HError(9905,"SetStateDurations: Cannot allocate memory for inverse variance");
   mean--; ivar--;
   
   /* prepare duration and calculate statistics to set speaking rate control parameter, rho */
   sum = sqr = 0.0;
   for (i=1; i<=genInfo->labseqlen; i++) {
      /* duration model for the i-th state */
      dm = genInfo->dm[i];
      
      /* # of states in the i-th model */
      nStates = genInfo->hmm[i]->numStates-2;
      mean[i] = CreateVector(genInfo->genMem, nStates);
      ivar[i] = CreateVector(genInfo->genMem, nStates);
      
      /* set statistics of the i-th state */
      for (s=cnt=1; s<=genInfo->dset->swidth[0]; s++) {
         for (k=1; k<=genInfo->dset->swidth[s]; k++,cnt++) {
            mean[i][cnt] = dm->svec[2].info->pdf[s].info->spdf.cpdf[1].mpdf->mean[k];
            ivar[i][cnt] = dm->svec[2].info->pdf[s].info->spdf.cpdf[1].mpdf->cov.var[k];  /* inverse variance (ConvDiagC was called when duration model set was loaded) */
         }
      }
      
      /* acc duration statistics to set rho */
      CountDurStat(mean[i], ivar[i], &sum, &sqr, genInfo->sindex[i]);
   }
   
   /* set rho, please refer to
    * T. Yoshimura, et al. "Duration Modeling in HMM-based Speech Synthesis System", 
    * Proc. of ICSLP, vol.2, pp.29-32, 1998, for detail 
    * */
   rho = (genInfo->speakRate*sum-sum)/sqr;
   
   /* set state durations of given label sequence */ 
   for (i=1; i<=genInfo->labseqlen; i++) {
      /* # of states in the i-th model */
      nStates = genInfo->hmm[i]->numStates-2;

      /* i-th label */
      label = genInfo->label[i];
      
      /* use model-level aligment */
      if (genInfo->modelAlign) {
         if (label->start<=0 && label->end<=0) {   /* model-level alignment of the i-th label is not specified */
            HError(-9999,"SetStateDurations: model duration is not specified in %d-th label",i); 
            rho = 1.0;
         }
         else {  /* model-level alignment of the i-th label is specified */
            modeldur = (int) (label->end - label->start)/genInfo->frameRate;
            sum = sqr = 0.0;
            CountDurStat(mean[i], ivar[i], &sum, &sqr, genInfo->sindex[i]);
            rho = (modeldur-sum)/sqr;
         }
      }
      
      /* calculate state durations for the i-th label */
      modeldur = 0;
      for (j=1; genInfo->sindex[i][j]!=0; j++) {
         genInfo->durations[i][j] = (int)(mean[i][genInfo->sindex[i][j]-1]+rho/ivar[i][genInfo->sindex[i][j]-1]+diff+0.5);
         
         /* set minimum duration -> 1 */
         if (genInfo->durations[i][j]<1)
            genInfo->durations[i][j] = 1;

         diff += mean[i][genInfo->sindex[i][j]-1]+rho/ivar[i][genInfo->sindex[i][j]-1]-(float)genInfo->durations[i][j];
         tframe += genInfo->durations[i][j];
         modeldur += genInfo->durations[i][j];
      }
      
      /* assign model duration */
      label->start = (HTime)start*genInfo->frameRate;
      label->end   = (HTime)(start+modeldur)*genInfo->frameRate;
      start += modeldur;
   }
   genInfo->tframe = tframe;
   
   /* free memory */
   Dispose(genInfo->genMem, ++mean);
   
   return;
}

/* GetLabStateDurations: parse state durations from label */
static void GetLabStateDurations (GenInfo *genInfo)
{
   int i,j,k,tframe=0;
   float diff=0.0;
   Label *label;
   
   /* get state durations from given label sequence */ 
   for (i=1,j=0,k=1; i<=CountLabs(genInfo->labseq->head); i++) {
      label = GetLabN(genInfo->labseq->head, i);
      if (label->auxLab[1]!=NULL) {
         j++; k=1;
      }
      
      /* get state duration from label */
      genInfo->durations[j][k] = (int)((label->end - label->start)/genInfo->frameRate+diff+0.5);
      diff += (label->end - label->start)/genInfo->frameRate - (float)genInfo->durations[j][k];
      
      /* count total frame */
      tframe += genInfo->durations[j][k++];
   }
   
   genInfo->tframe = tframe;
   
   return;
}

/* SetSpaceIndexes: set space indexes for each PdfStream according to MSD threshold */
static void SetSpaceIndexes (MemHeap *x, GenInfo *genInfo)
{
   int i,j,k,s,m,t,p,stream;
   float ContSpaceWeight;
   Boolean ContSpace;
   PdfStream *pst;
   StateInfo *si;
   
   /* initialise space indexes and total frame for each PdfStream  */
   for (p=1; p<=genInfo->nPdfStream[0]; p++) {
      /* p-th PdfStream */
      pst = &(genInfo->pst[p]);
      
      /* prepare space indexes */
      pst->ContSpace = (Boolean *) New(x, genInfo->tframe*sizeof(Boolean));
      pst->ContSpace--;
      
      /* initialize space indexes and number of continuous space */
      for (t=1; t<=genInfo->tframe; t++)
         genInfo->pst[p].ContSpace[t] = FALSE;
      genInfo->pst[p].T = 0;
   }
   
   /* select space according to MSD threshold */
   for (i=1,t=1; i<=genInfo->labseqlen; i++) {
      /* determine continuous space or not */ 
      for (j=1; genInfo->sindex[i][j]!=0; j++) {
         si = genInfo->hmm[i]->svec[genInfo->sindex[i][j]].info;
         for (p=stream=1; p<=genInfo->nPdfStream[0]; stream+=genInfo->nPdfStream[p++]) {
            ContSpace = FALSE;
            for (s=stream; s<stream+genInfo->nPdfStream[p]; s++) {
               if (genInfo->hset->msdflag[s]) {
                  ContSpaceWeight = 0.0;
                  for (m=1; m<=si->pdf[s].info->nMix; m++)
                     if (VectorSize(si->pdf[s].info->spdf.cpdf[m].mpdf->mean)==genInfo->hset->swidth[s])  /* total weight of all continuous mixtures */
                        ContSpaceWeight += MixWeight(genInfo->hset, si->pdf[s].info->spdf.cpdf[m].weight);
                  
                  /* if any streams in the p-th PdfStream is determined to continuous, this frame determined to be continous */
                  if (ContSpaceWeight > genInfo->MSDthresh) {
                     ContSpace = TRUE;
                     break;
                  }
               }
               else {
                  ContSpace = TRUE;
                  break;
               }
            }
               
            /* set frames belonging to this state continuous */
            if (ContSpace) {
               /* set ContSpace flag */
               for (k=0; k<genInfo->durations[i][j]; k++)
                  genInfo->pst[p].ContSpace[t+k] = TRUE;
                  
               /* update total number of frame in this PdfStream */
               genInfo->pst[p].T += genInfo->durations[i][j];
            }
         }
         t += genInfo->durations[i][j];
      }
   }
   
   if (trace & T_TOP) { 
      printf(" Total number of frames = %d\n", genInfo->tframe);
      for (p=1; p<=genInfo->nPdfStream[0]; p++) {
         printf("  PdfStream[%d]: %d frames\n", p, genInfo->pst[p].T);
      }
      fflush(stdout);
   }  
   
   return;
}
   
/* EXPORT->InitialiseGenInfo: initialize a genInfoence HMM according to the given label */
void InitialiseGenInfo (GenInfo *genInfo, Transcription *tr)
{
   int i, j, n=0, max=0;
   MLink hmacro, dmacro;
   Label *label;
   LabId id;
      
   /* set input label file */
   genInfo->labseq = tr;
   
   /* count # of individual models in this label */
   if (genInfo->stateAlign) {
      for (i=1,genInfo->labseqlen=0; i<=CountLabs(genInfo->labseq->head); i++) {
         label = GetLabN(genInfo->labseq->head, i);
         if (i==1 && label->auxLab==NULL)
            HError(9999, "InitGenInfo: Invalid label format for parameter generation with state alignments");      
         if (label->auxLab[1]!=NULL)  /* usually model name is written in auxLab[1] */
            genInfo->labseqlen++;
      }
   }
   else
      genInfo->labseqlen = CountLabs(genInfo->labseq->head);
   
   /* create label storage */
   if ((genInfo->label = (Label **) New(genInfo->genMem, genInfo->labseqlen*sizeof(Label *)))==NULL)
      HError(9905,"InitGenInfo: Cannot allocate memory for labels");
   genInfo->label--;
    
    
   /* create hmm storage */
   if ((genInfo->hmm = (HLink *) New(genInfo->genMem, genInfo->labseqlen*sizeof(HLink)))==NULL)
      HError(9905,"InitGenInfo: Cannot allocate memory for HMMs");
   genInfo->hmm--;

   /* create duration model storage */   
   if (!genInfo->stateAlign) {
      if ((genInfo->dm = (HLink *) New(genInfo->genMem, genInfo->labseqlen*sizeof(HLink)))== NULL)
         HError(9905,"InitGenInfo: Cannot allocate memory for duration models");    
      genInfo->dm--;
   }
   
   /* parse label and compose a genInfoence HMM with state duration models */
   for (i=1,j=1; i<=CountLabs(genInfo->labseq->head); i++) {
      /* get label */
      label = GetLabN(genInfo->labseq->head, i);
      
      /* LabId of this model */
      if (genInfo->stateAlign) {
         n++;
         if ((id=label->auxLab[1]) == NULL) {  /* auxLab[1] == NULL -> ignore */
            continue;
         }
         else {
            max = (n>max) ? n : max;
            n = 0;
         }
      }
      else {
         /* find state duration model */
         id = label->labid;
         if ((dmacro = FindMacroName(genInfo->dset, 'l', id)) == NULL)
            HError(9935,"Generator: Cannot find duration model %s in current list", id->name);
         genInfo->dm[j]  = (HLink) dmacro->structure;
      }
      
      /* find HMM */
      if ((hmacro = FindMacroName(genInfo->hset, 'l', id)) == NULL)
         HError(9935,"Generator: Cannot find hmm %s in current model list", id->name);
      genInfo->hmm[j] = (HLink) hmacro->structure;
      
      /* set label */
      genInfo->label[j] = label;
      
      j++;
   }
   
   
   /* set state sequence which maximizes its state sequence prob */
   genInfo->sindex = CreateIMatrix(genInfo->genMem, genInfo->labseqlen, ((genInfo->stateAlign)?max+1:genInfo->maxStates));
   ZeroIMatrix(genInfo->sindex);
   SetStateSequence(genInfo);
   
   /* set state durations which maximize their state duration prob */
   genInfo->durations = CreateIMatrix(genInfo->genMem, genInfo->labseqlen, ((genInfo->stateAlign)?max+1:genInfo->maxStates));
   ZeroIMatrix(genInfo->durations);
   if (genInfo->stateAlign)
      GetLabStateDurations(genInfo);
   else 
      SetStateDurations(genInfo);

      
   /* set MSD space indexes for each PdfStream */
   SetSpaceIndexes(genInfo->genMem, genInfo);

   
   /* Create PdfStreams */
   CreatePdfStreams(genInfo);
  
   return;
}

/* EXPORT->ResetGenInfo: reset a genInfoence HMM according to the given label */
void ResetGenInfo (GenInfo *genInfo)
{
   Dispose(genInfo->genMem, ++genInfo->label);
   
   return;
}

/* JointProb: joint probability of given observations and state sequence */ 
static void JointProb (GenInfo *genInfo, UttInfo *utt)
{
   int i, j, d, t;
   LogFloat prob=0.0;
   StateInfo *si;
      
   /* compute output probability of given observation sequence 
    * according to state durations */
   for (i=1,t=1; i<=genInfo->labseqlen; i++) {
      /* initial state prob */
      prob += genInfo->hmm[i]->transP[1][ genInfo->sindex[i][1] ];
      
      /* state output and self-transition prob */
      for (j=1; genInfo->sindex[i][j]!=0; j++) {
         si = genInfo->hmm[i]->svec[genInfo->sindex[i][j]].info;
         
         /* state output prob */
         for (d=1; d<=genInfo->durations[i][j]; d++,t++)
            prob += POutP(genInfo->hset, &(utt->o[t]), si);
         
         /* state transition prob */
         prob += (genInfo->durations[i][j] - 1) * genInfo->hmm[i]->transP[ genInfo->sindex[i][j] ][ genInfo->sindex[i][j] ]; 
         if (genInfo->sindex[i][j+1]!=0)
            prob += genInfo->hmm[i]->transP[ genInfo->sindex[i][j] ][ genInfo->sindex[i][j+1] ]; 
      }
      
      /* final state prob */
      prob += genInfo->hmm[i]->transP[ genInfo->sindex[i][j-1] ][ genInfo->hmm[i]->numStates ];
   }
   
   utt->pr = prob;
}

/* OutProb: output probability of an observation sequence for a given state sequence */ 
static void OutProb (GenInfo *genInfo, UttInfo *utt)
{
   int i, j, d, t;
   LogFloat prob=0.0;
   StateInfo *si;
     
   for (i=1,t=1; i<=genInfo->labseqlen; i++) {
      for (j=1; genInfo->sindex[i][j]!=0; j++) {
         si = genInfo->hmm[i]->svec[genInfo->sindex[i][j]].info;
         
         /* state output prob */
         for (d=1; d<=genInfo->durations[i][j]; d++,t++)
            prob += POutP(genInfo->hset, &(utt->o[t]), si);
      }
   }
   
   utt->pr = prob;
}

/* -------------------------- Cholesky decomposition-based parameter generation -------------------------- */

/* Calc_WUM_and_WUW: calcurate W'*U^{-1}*M and W'*U^{-1}*W */
static void Calc_WUM_and_WUW (PdfStream *pst, const int bias)
{
   int t, m, n, d, j, l, k;
   double cov, WU;

   const Boolean full = pst->fullCov;
   const int M = (full) ? pst->order : 1;

   /* initialization */
   ZeroDMatrix(pst->WUW);
   ZeroDVector(pst->WUM);

   /* computation */
   #pragma omp parallel for private(m,n,d,j,l,cov,WU,k)
   for (t=1; t<=pst->T; t++) {
      for (m=1; m<=M; m++) {
         for (n=1; n<=M; n++) {
            for (d=0; d<pst->win.num; d++) {
               for (j=pst->win.maxw[WLEFT]; j<=pst->win.maxw[WRIGHT]; j++) {
                  if ((t+j>0) && (t+j<=pst->T)) {
                     /* accumulate W'*U^{-1}*M */
                     if ((n==1) && (pst->win.width[d][WLEFT]<=j) && (j<=pst->win.width[d][WRIGHT]) && (pst->win.coef[d][-j]!=0.0))
                        pst->WUM[M*(t-1)+m] += ((double)pst->win.coef[d][-j]) * pst->mseq[t+j][d*pst->order+m+bias];

                     /* accumulate W'*U^{-1}*W */
                     /* W'U^{-1} */
                     for (l=((full)?0:d),WU=0.0; l<=((full)?pst->win.num-1:d); l++) {
                        cov = (!full) ? pst->vseq[t+j].var[pst->order*l+m+bias] :
                              ((pst->order*l+m > pst->order*d+n) ? pst->vseq[t+j].inv[pst->order*l+m][pst->order*d+n]
                                                                 : pst->vseq[t+j].inv[pst->order*d+n][pst->order*l+m]);

                        if (cov!=0.0 && pst->win.width[l][WLEFT] <=j && j<=pst->win.width[l][WRIGHT] && pst->win.coef[l][-j]!=0.0)
                           WU += cov * (double)pst->win.coef[l][-j];
                     }

                     /* W'*U^{-1}*W */
                     for (k=0; (WU!=0.0) && (k<pst->width) && (t+k<=pst->T); k++)
                        if ((pst->win.width[d][WLEFT]<=k-j) && (k-j<=pst->win.width[d][WRIGHT]) && (M*k+n-m+1>0) && (pst->win.coef[d][k-j]!=0.0))
                           pst->WUW[M*(t-1)+m][M*k+n-m+1] += WU*(double)pst->win.coef[d][k-j];
                  }
               }
            }
         }
      }
   }

   if (trace & T_MAT) {
      ShowDMatrix  ("  WUW", pst->WUW, pst->width, pst->T);
      ShowDVector("\n  WUM", pst->WUM, pst->T);
      fflush(stdout);
   }
}

/* Cholesky_Factorization: Compute Cholesky factor of matrix W'*U^{-1}*W */
static void Cholesky_Factorization (PdfStream *pst)
{
   int t,i,j;

   DMatrix WUW = pst->WUW;
   DMatrix U   = pst->WUW;

   /* sizes of matrix */
   const int M = (pst->fullCov) ? pst->order : 1;
   const int T = M*pst->T;
   const int width = M*pst->width;

   /* Cholesky decomposition */
   for (t=1; t<=T; t++) {
      for (i=1; (i<width) && (t-i>0); i++)
         U[t][1] -= U[t-i][i+1]*U[t-i][i+1];

      if (WUW[t][1]<0.0)
         HError(9999,"Cholesky_Factorization: (%d,%d)-th element of W'*U^{-1}*W is negative.\n",t,t);

      U[t][1] = sqrt(WUW[t][1]);

      for (i=2; i<=width; i++) {
         for (j=1; (i+j<=width) && (t-j>0); j++)
            U[t][i] -= U[t-j][j+1]*U[t-j][i+j];
         U[t][i] /= U[t][1];
      }
   }

   if (trace & T_MAT) {
      ShowDMatrix("\n  Cholesky factor", U, width, T);
      fflush(stdout);
   }

   return;
}

/* Forward_Substitution: forward substitution to solve set of linear equations */
static void Forward_Substitution (PdfStream *pst)
{
   int t,i;

   DMatrix U = pst->WUW;
   DVector r = pst->WUM;
   DVector g = pst->g;

   /* sizes of matrix and vector */
   const int M = (pst->fullCov) ? pst->order : 1;
   const int T = M*pst->T;
   const int width = M*pst->width;

   /* forward substitution */
   for (t=1; t<=T; t++) {
      g[t] = r[t];
      for (i=1; (i<width) && (t-i>0); i++)
         g[t] -= U[t-i][i+1] * g[t-i];
      g[t] /= U[t][1];
      
      /* random generation */
      if (rndpg)
         g[t] += GaussDeviate(rndMean, rndVar);
   }
   
   if (trace & T_MAT) {
      ShowDVector("\n  g", g, T);
      fflush(stdout);
   }

   return;
}

/* Backward_Substitution: backward substitution to solve set of linear equations */
static void Backward_Substitution (PdfStream *pst, const int bias)
{
   int t,i;
   double c;

   Matrix  C = pst->C;
   DMatrix U = pst->WUW;
   DVector g = pst->g;

   /* sizes of matrix and vector */
   const int M = (pst->fullCov) ? pst->order : 1;
   const int T = M*pst->T;
   const int width = M*pst->width;

   if (trace & T_MAT)
      printf("\n  solution\n   ");

   /* backward substitution */
   for (t=T; t>0; t--) {
      c = g[t];
      for (i=1; (i<width) && (t+i<=T); i++)
         c -= U[t][i+1]*C[(t+i+M-1)/M][(t+i+M-1)%M+1+bias];
      c /= U[t][1];
      C[(t+M-1)/M][(t+M-1)%M+1+bias] = (float) c;

      if (trace & T_MAT) {
         printf("%8.2f ", C[(t+M-1)/M][(t+M-1)%M+1+bias]);
         fflush(stdout);
      }
   }

   return;
}

/* Cholesky_ParmGen: Generate parameter sequence using Cholesky decomposition */
static void Cholesky_ParmGen (GenInfo *genInfo)
{
   int p,m;
   PdfStream *pst;

   for (p=1; p<=genInfo->nPdfStream[0]; p++) {
      /* p-th PdfStream */
      pst = &(genInfo->pst[p]);
      
      if (pst->T<1)
         continue;

      if (pst->fullCov) {
         /* full covariance */
         /* generate all feature simultaneously */
         Calc_WUM_and_WUW(pst, 0);
         Cholesky_Factorization(pst);    /* Cholesky decomposition */
         Forward_Substitution(pst);      /* forward substitution   */
         Backward_Substitution(pst, 0);  /* backward substitution  */
      }
      else {
         /* diagonal covariance */
         for (m=1; m<=pst->order; m++) {
            if (trace & T_MAT) {
               printf("  Feature %d:\n", m);
               fflush(stdout);
            }

            /* generate m-th feature */
            Calc_WUM_and_WUW(pst, m-1);
            Cholesky_Factorization(pst);     /* Cholesky decomposition */
            Forward_Substitution(pst);       /* forward substitution   */
            Backward_Substitution(pst, m-1); /* backward substitution  */
         }
      }
   }

   return;
}

/* -------------------------- Observation generation routines --------------------------------- */
/* ApplyWindow: check t-th frame is on swithing space boundary */
static float ApplyWindow (PdfStream *pst, const short msdflag, const int absolute_t, const int v)
{
   int j;
   float otsk;

   /* constants */
   const int t = pst->t;  /* time counter in this PdfStream */
   const int T = pst->T;  /* total number of frames in this PdfStream */
   const int i = (v-1) / pst->order;      /* window for this dimension */
   const int m = (v-1) % pst->order + 1;  /* static feature for this dimension */

   /* applying window */
   otsk = 0.0;
   for (j=pst->win.width[i][0]; j<=pst->win.width[i][1]; j++) {
      if ( msdflag && (t+j<1 || t+j>T || (pst->win.coef[i][j]!=0.0 && !pst->ContSpace[absolute_t+j])) ) {
         otsk = ReturnIgnoreValue();
         break;
      }
      else {
         if (t+j>=1 && t+j<=T)
            otsk += pst->win.coef[i][j] * pst->C[t+j][m];
      }
   }

   return(otsk);
}

/* UpdateUttObs: update observations in UttInfo using generated static feature vector sequence */
static void UpdateUttObs (GenInfo *genInfo, UttInfo *utt)
{
   int t,s,k,p,stream,v;
   PdfStream *pst;

   const int T = utt->T;

   /* initialize PdfStream time counter */
   for (p=1; p<=genInfo->nPdfStream[0]; p++)
      genInfo->pst[p].t=1;

   /* construct observation vector sequence */
   for (t=1; t<=T; t++) {
      for (p=stream=1; p<=genInfo->nPdfStream[0]; stream+=genInfo->nPdfStream[p++]) {
         /* p-th PdfStream */
         pst = &(genInfo->pst[p]);
         
         /* whether t-th frame is continuous or not */
         if (pst->ContSpace[t]) {
            for (s=stream,v=1; s<stream+genInfo->nPdfStream[p]; v+=genInfo->hset->swidth[s++]) {
               for (k=1; k<=genInfo->hset->swidth[s]; k++)
                  utt->o[t].fv[s][k] = ApplyWindow(pst, genInfo->hset->msdflag[s], t, v+k-1);
            }
            pst->t++;
         }
         else {
            for (s=stream; s<stream+genInfo->nPdfStream[p]; s++)
               for (k=1; k<=genInfo->hset->swidth[s]; k++)
                  utt->o[t].fv[s][k] = ReturnIgnoreValue();   /* ignoreValue is used for repregenInfoing MSD discrete symbol */
         }
      }
   }

   return;
}

/* -------------------------- EM-based parameter generation algorithm -------------------------- */

/* UpdatePdfStreams: update PdfStreams according to occ prob */
static void UpdatePdfStreams (GenInfo *genInfo, FBInfo *fbInfo, UttInfo *utt)
{
   int t, p, q, j, N, k, l, m, s, v, stream;
   float Lr;
   AlphaBeta *ab = fbInfo->ab;
   HLink hmm;
   MixPDF *mpdf;
   Vector mseq_t;
   Covariance vseq_t;
   Boolean bound;
   PdfStream *pst;

   /* absolute time */
   const int T = utt->T;

   /* initialize time counter and statistics of each PdfStream */
   for (p=1; p<=genInfo->nPdfStream[0]; p++) {
      /* p-th PdfStream */
      pst = &(genInfo->pst[p]);
      
      /* initialize time counter, mseq, and vseq */
      pst->t = 1;
      for (t=1; t<=pst->T; t++) {
         ZeroVector(pst->mseq[t]);
         if (pst->fullCov)
            ZeroTriMat(pst->vseq[t].inv);
         else
            ZeroVector(pst->vseq[t].var);
      }
   }

   /* update statistics */
   for (t=1; t<=T; t++) {   /* absolute time */
      for (p=stream=1; p<=genInfo->nPdfStream[0]; stream+=genInfo->nPdfStream[p++]) {
         /* p-th PdfStream */
         pst = &(genInfo->pst[p]);
      
         /* whether t-th frame is continuous or not */
         if (pst->ContSpace[t]) {
            /* mean vector and covariance matrix at this frame */
            mseq_t = pst->mseq[ pst->t ];
            vseq_t = pst->vseq[ pst->t ];

            /* update mseq and vseq at time t */
            for (q=ab->pInfo->qLo[t]; q<=ab->pInfo->qHi[t]; q++) {
               hmm = genInfo->hmm[q];
               N = hmm->numStates;

               for (j=2; j<N; j++) {
                  /* update statistics */
                  for (s=stream,v=1; s<stream+genInfo->nPdfStream[p]; v+=fbInfo->al_hset->swidth[s++]) {
                     for (m=1; m<=hmm->svec[j].info->pdf[s].info->nMix; m++) {
                        mpdf = hmm->svec[j].info->pdf[s].info->spdf.cpdf[m].mpdf;
                        if (VectorSize(mpdf->mean) == fbInfo->al_hset->swidth[s]) {  /* check MSD */
                           Lr = ab->occm[t][q][j][s][m];  /* absolute time */
                           if (Lr > 0.0) {
                              for (k=1; k<=fbInfo->al_hset->swidth[s]; k++) {
                                 bound = ChkBoundary(pst, v+k-1, t, T);
                                 switch(mpdf->ckind) {
                                 case DIAGC:
                                    mseq_t    [v+k-1] += (bound) ? 0.0 : Lr / mpdf->cov.var[k] * mpdf->mean[k];
                                    vseq_t.var[v+k-1] += (bound) ? 0.0 : Lr / mpdf->cov.var[k];
                                    break;
                                 case INVDIAGC:
                                    mseq_t    [v+k-1] += (bound) ? 0.0 : Lr * mpdf->cov.var[k] * mpdf->mean[k];
                                    vseq_t.var[v+k-1] += (bound) ? 0.0 : Lr * mpdf->cov.var[k];
                                    break;
                                 case FULLC:
                                    /* diagonal elements */
                                    mseq_t    [v+k-1]        += (bound) ? 0.0 : Lr * mpdf->cov.inv[k][k] * mpdf->mean[k];
                                    vseq_t.inv[v+k-1][v+k-1] += (bound) ? 0.0 : Lr * mpdf->cov.inv[k][k];

                                    /* off-diagonal elements */
                                    for (l=1; l<k; l++) {
                                       mseq_t    [v+k-1]        += (bound) ? 0.0 : Lr * mpdf->cov.inv[k][l] * mpdf->mean[l];
                                       mseq_t    [v+l-1]        += (bound) ? 0.0 : Lr * mpdf->cov.inv[k][l] * mpdf->mean[k];
                                       vseq_t.inv[v+k-1][v+l-1] += (bound) ? 0.0 : Lr * mpdf->cov.inv[k][l];
                                    }
                                    break;
                                 default:
                                    HError(9999, "UpPdfStreams: Only DIAGC, INVDIAGC and FULLC are supported.");
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
            /* update counter */
            pst->t++;
         }
      }
   }

   return;
}

/* MixUtt: compute mixture-level posterior prob */
static Boolean MixUtt (GenInfo *genInfo, FBInfo *fbInfo, UttInfo *utt)
{
   int t,d,s,i,j,k,l,m;
   LogDouble ax,bx;
   AlphaBeta *ab = fbInfo->ab;
   PruneInfo *p;
   StateInfo *si;
   MixPDF *mp;
      
   ab->pInfo = (PruneInfo *) New(&ab->abMem, sizeof(PruneInfo));
   p = ab->pInfo;
   
   p->qHi = CreateShortVec(&ab->abMem, genInfo->tframe);
   p->qLo = CreateShortVec(&ab->abMem, genInfo->tframe);
   
   ab->occm = (Vector ****) New (&ab->abMem, genInfo->tframe*sizeof(Vector ***));
   ab->occm--;
   
   for (i=1,t=1; i<=genInfo->labseqlen; i++) {
      for (j=1; genInfo->sindex[i][j]!=0; j++) {
         for (d=1; d<=genInfo->durations[i][j]; d++,t++) {
            /* set qLo & qHi */
            p->qLo[t] = p->qHi[t] = i; 

            /* occupancy prob storage */
            ab->occm[t] = (Vector ***) New(&ab->abMem, sizeof(Vector **));
            ab->occm[t] -= i;
      
            ab->occm[t][i] = (Vector **) New(&ab->abMem, genInfo->hmm[i]->numStates*sizeof(Vector *));
            ab->occm[t][i]--;
      
            for (k=2; k<genInfo->hmm[i]->numStates; k++) {
               ab->occm[t][i][k] = (Vector *) New(&ab->abMem, genInfo->hset->swidth[0]*sizeof(Vector));
               ab->occm[t][i][k]--;
            
               /* first compute output prob of each mix */
               si = genInfo->hmm[i]->svec[k].info; 
               for (s=1; s<=genInfo->hset->swidth[0]; s++) {
                  ab->occm[t][i][k][s] = CreateVector(&ab->abMem, si->pdf[s].info->nMix);
                  for (m=1; m<=si->pdf[s].info->nMix; m++) {
                     mp = si->pdf[s].info->spdf.cpdf[m].mpdf;
                     if (k==genInfo->sindex[i][j])
                        ab->occm[t][i][k][s][m] = MixLogWeight(genInfo->hset, si->pdf[s].info->spdf.cpdf[m].weight) 
                                                + MOutP(utt->o[t].fv[s], mp);   
                     else
                        ZeroVector(ab->occm[t][i][k][s]);
                  }
               }
            }
            
            /* then set mix-level posterior prob */
            k = genInfo->sindex[i][j];
            si = genInfo->hmm[i]->svec[k].info; 
            for (s=1; s<=genInfo->hset->swidth[0]; s++) {
               for (l=1,ax=0.0; l<=genInfo->hset->swidth[0]; l++) {
                  if (l!=s) {
                     for (m=1,bx=LZERO; m<=si->pdf[l].info->nMix; m++)
                        bx = LAdd(bx, ab->occm[t][i][k][l][m]);
                     ax += si->weights[l] * bx;
                  }
               }
               
               for (m=1,bx=LZERO; m<=si->pdf[s].info->nMix; m++)
                  bx = LAdd(bx, ab->occm[t][i][k][s][m]);
               bx += ax;
               
               for (m=1; m<=si->pdf[s].info->nMix; m++)
                  ab->occm[t][i][k][s][m] = L2F(ab->occm[t][i][k][s][m] + ax - bx);
            }
         }
      }
   }  
   
   OutProb(genInfo, utt);
   
   return ((utt->pr<LSMALL) ? FALSE : TRUE);
}

/* EXPORT->ParamGen: Generate parameter sequence */
void ParamGen (GenInfo *genInfo, UttInfo *utt, FBInfo *fbInfo, const ParmGenType type)
{
   int n=1;
   Boolean success=TRUE;
   LogFloat prev, curr=LZERO;

   /* # of frames */
   const int T = genInfo->tframe;
   
   /* UttInfo settings */   
   utt->tgtSampRate = genInfo->frameRate;
   utt->S = genInfo->hset->swidth[0];
   
   /* First perform Cholesky-based parameter generation */
   SetupPdfStreams(genInfo);
   Cholesky_ParmGen(genInfo);
   UpdateUttObs(genInfo, utt);

   /* Cholesky case */
   if (type==CHOLESKY) {
      OutProb(genInfo, utt);
      if (trace&T_TOP)
         printf("  Average LogP = %e\n", utt->pr/T);
      return;
   }

   /* EM-based parameter generation */
   if (trace & T_TOP) {
      printf(" EM-based parameter generation ");
      switch(type) {
      case MIX: printf("(Hidden: MIX, Given: STATE)\n"); break;
      case FB:  printf("(Hidden: MIX, STATE)\n");
      }
      fflush(stdout);
   }

   /* Optimize pst->C using EM algorithm */
   do {
      switch(type) {
      case MIX: 
         success = MixUtt(genInfo, fbInfo, utt);  /* compute mixture-level posterior */
         break;
      case FB:
         success = FBUtt(fbInfo, utt);  /* perform forward-backward */
         break;
      default:
         HError(9999,"ParamGen: not supported parameter generation type");
      } 
      if (!success)
         HError(9999,"ParamGen: failed to compute output prob");

      /* output prob */
      prev = curr;
      curr = utt->pr;

      /* generate parameters */
      UpdatePdfStreams(genInfo, fbInfo, utt);  /* update PdfStreams */
      Cholesky_ParmGen(genInfo);
      UpdateUttObs(genInfo, utt);
      
      /* reset heap */
      ResetHeap(&fbInfo->ab->abMem);
      
      if (trace & T_TOP) {
         printf("  Iteration %d: Average LogP = %e", n, curr/T);
         if (n>1)
            printf("  Change = %f", (curr-prev)/T);
         printf("\n");
         fflush(stdout);
      }
   } while(fabs((curr-prev)/T)>epsilon && n++<=maxIter);

   return;
}

/* ---------------------------- End of HGen.c --------------------------- */
