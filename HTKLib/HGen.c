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
char *hgen_vc_id = "$Id: HGen.c,v 1.38 2008/01/11 02:58:08 zen Exp $";

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
#include "lbfgs.h"

#define LBFGSMEM 7
#define MAXWINLEN 10

/* ------------------- Trace Information ------------------------------ */
/* Trace Flags */
#define T_TOP   0001  /* Top level tracing */
#define T_MAT   0002  /* trace matrices */
#define T_STA   0004  /* trace state sequence */
#define T_GV    0010  /* trace gv param gen */

static int trace =  0010;

static ConfParam *cParm[MAXGLOBS];  /* config parameters */
static int nParm = 0;

static int maxEMIter = 20;          /* max iterations in the EM-based parameter generation */
static double EMepsilon = 1.0E-4;   /* convergence factor for EM iteration */

typedef enum { RNDTRANS=1,RNDDUR=2,RNDSPACE=4,RNDMIX=8,RNDPAR=16} RNDSet;
static RNDSet rFlags  = (RNDSet)0;  /* random generation flag */
static double rndParMean = 0.0;     /* mean of Gaussian noise for random generation */
static double rndParVar  = 1.0;     /* variance of Gaussian noise for random generation */

/* GV related variables */
static Boolean useGV = FALSE;
static int maxGVIter = 50;          /* max iterations in the speech parameter generation considering GV */
static double GVepsilon = 1.0E-4;   /* convergence factor for GV iteration */
static double minEucNorm = 1.0E-2;  /* minimum Euclid norm of a gradient vector */ 
static double stepInit = 1.0;       /* initial step size */
static double stepDec = 0.5;        /* step size deceralation factor */
static double stepInc = 1.2;        /* step size acceleration factor */
static double w1 = 1.0;             /* weight for HMM output prob. */
static double w2 = 1.0;             /* weight for GV output prob. */

typedef enum {STEEPEST=0,NEWTON=1,LBFGS=2} OptKind;
static OptKind optKind = LBFGS;     /* optimization method */
  

static char gvDir[MAXFNAMELEN];  /* dir to look for GV defs */
static char gvExt[MAXSTRLEN];    /* GV def file extension */
static char gvMMF[MAXFNAMELEN];  /* GV MMF */
static char gvLst[MAXFNAMELEN];  /* GV list */
static HMMSet gvset;             /* GV set */
static MemHeap gvStack;          /* Stack holds all GV related info */

/* -------------------------- Initialisation ----------------------- */

/* Str2OptKind: parse the string into the correct OptKind */
static OptKind Str2OptKind(char *str)
{
   OptKind okind = STEEPEST;
   if (!(strcmp(str,"STEEPEST"))) okind = STEEPEST;
   else if (!(strcmp(str,"NEWTON"))) okind = NEWTON;
   else if (!(strcmp(str,"LBFGS")))  okind = LBFGS;
   else HError(9999,"Str2OptKind: Unknown OptKind kind");
   return okind;
}

/* OptKind2Str: Return string representation of enum OptKind */
static char *OptKind2Str(OptKind okind, char *buf)
{
   static char *okindmap[] = {"STEEPEST","NEWTON","LBFGS"};
   return strcpy(buf,okindmap[okind]);
}

/* EXPORT->SetrFlags: set random generation flags */
void SetrFlags(char *s)
{
   rFlags=(RNDSet) 0;        
   while (*s != '\0')
      switch (*s++) {
      case 't': rFlags = (RNDSet) (rFlags+RNDTRANS); break;
      case 'd': rFlags = (RNDSet) (rFlags+RNDDUR); break;
      case 's': rFlags = (RNDSet) (rFlags+RNDSPACE); break;
      case 'm': rFlags = (RNDSet) (rFlags+RNDMIX); break;
      case 'p': rFlags = (RNDSet) (rFlags+RNDPAR); break;
      default: HError(2320,"SetrFlags: Unknown random generation flag %c",*s);
         break;
      }
}

/* EXPORT->PrintrFlags: print random generation flags */
void PrintrFlags (void)
{
   printf("HMGenS  ");
   
   if (!rFlags) printf("ML Generating: ");
   if (!(rFlags&RNDTRANS)) printf("Transitions "); 
   if (!(rFlags&RNDDUR))   printf("Durations ");
   if (!(rFlags&RNDSPACE)) printf("Spaces ");
   if (!(rFlags&RNDMIX))   printf("Mixtures ");
   if (!(rFlags&RNDPAR))   printf("Parameters ");
   if (!rFlags) printf("   ");

   if (rFlags) printf("Random Generating: ");
   if (rFlags&RNDTRANS) printf("Transitions "); 
   if (rFlags&RNDDUR)   printf("Durations ");
   if (rFlags&RNDSPACE) printf("Spaces ");
   if (rFlags&RNDMIX)   printf("Mixtures ");
   if (rFlags&RNDPAR)   printf("Parameters ");
   printf("\n\n");
}

/* EXPORT->InitGen: initialise module */
void InitGen(void)
{
   int i;
   double d;
   Boolean b;
   char buf[MAXSTRLEN];

   Register(hgen_version,hgen_vc_id);

   nParm = GetConfig("HGEN", TRUE, cParm, MAXGLOBS);
   if (nParm>0){
      if (GetConfInt (cParm,nParm,"TRACE",     &i))  trace      = i;
      if (GetConfInt (cParm,nParm,"MAXEMITER", &i))  maxEMIter  = i;
      if (GetConfFlt (cParm,nParm,"EMEPSILON", &d))  EMepsilon  = d;
      if (GetConfFlt (cParm,nParm,"RNDPARMEAN",&d))  rndParMean = d;
      if (GetConfFlt (cParm,nParm,"RNDPARVAR", &d))  rndParVar  = d;
      if (GetConfBool(cParm,nParm,"USEGV",     &b))  useGV      = b;
      if (GetConfInt (cParm,nParm,"MAXGVITER", &i))  maxGVIter  = i;
      if (GetConfFlt (cParm,nParm,"GVEPSILON", &d))  GVepsilon  = d;
      if (GetConfFlt (cParm,nParm,"MINEUCNORM",&d))  minEucNorm = d;
      if (GetConfFlt (cParm,nParm,"STEPINIT" , &d))  stepInit   = d;
      if (GetConfFlt (cParm,nParm,"STEPDEC",   &d))  stepDec    = d;
      if (GetConfFlt (cParm,nParm,"STEPINC",   &d))  stepInc    = d;
      if (GetConfFlt (cParm,nParm,"HMMWEIGHT", &d))  w1         = d;
      if (GetConfFlt (cParm,nParm,"GVWEIGHT",  &d))  w2         = d;
      if (GetConfStr (cParm,nParm,"OPTKIND",   buf)) optKind    = Str2OptKind(buf);
      if (GetConfStr (cParm,nParm,"RNDFLAGS",  buf)) SetrFlags(buf);
      if (GetConfStr (cParm,nParm,"GVMODELMMF",buf)) strcpy(gvMMF,buf);
      if (GetConfStr (cParm,nParm,"GVHMMLIST", buf)) strcpy(gvLst,buf);
      if (GetConfStr (cParm,nParm,"GVMODELDIR",buf)) strcpy(gvDir,buf);
      if (GetConfStr (cParm,nParm,"GVMODELEXT",buf)) strcpy(gvExt,buf);
   }

   if (useGV) {
      /* Stack for GV */
      CreateHeap(&gvStack, "gvStore", MSTAK, 1, 1.0, 50000, 500000);

      /* load GV set */
      CreateHMMSet(&gvset,&gvStack,TRUE);
      AddMMF(&gvset,gvMMF);
      if (MakeHMMSet(&gvset, gvLst)<SUCCESS)
         HError(9928,"InitGen: MakeHMMSet failed");
      if (LoadHMMSet(&gvset, gvDir, gvExt)<SUCCESS)
         HError(9928,"InitGen: LoadHMMSet failed");
      ConvDiagC(&gvset,TRUE);
      ConvLogWt(&gvset);

      /* and echo status */
      if (trace&T_TOP) {
         printf("GV enabled\n");
         printf(" %d Logical/%d Physical Models Loaded, VecSize=%d\n", gvset.numLogHMM, gvset.numPhyHMM, gvset.vecSize);
      }
   }
   
   return;
}

/* EXPORT->ResetGen: reset module */
void ResetGen(void)
{
   if (useGV) 
      ResetHeap(&gvStack);

   return;
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
      if (winlen<0)
         HError(9999,"InitWindow: Length of %d-th window (%d) is invalid", i+1, winlen);
      if (winlen>MAXWINLEN)
         HError(9999,"InitWindow: Length of %d-th window (%d) exceeds maximum length (%d)", i+1, winlen, MAXWINLEN);

      /* Read window coefficients */
      pst->win.coef[i] = (float *) New(x,winlen*sizeof(float));
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
   PdfStream *pst;
   MemHeap *x = genInfo->genMem;

   for (p=1; p<=genInfo->nPdfStream[0]; p++) {
      /* p-th PdfStream */
      pst = &(genInfo->pst[p]);
      
      /* Load window file for p-th PdfStream */
      InitWindow(x, pst);

      /* full or diag covariance */
      pst->fullCov = ((genInfo->hset->ckind==FULLC) || (useGV && (gvset.ckind==FULLC))) ? TRUE : FALSE;

      /* prepare mean vector and inverse covariance matrix sequences */
      pst->mseq = CreateMatrix(x, pst->T, pst->vSize);               /* mean vectors */
      pst->vseq = (Covariance *) New(x, pst->T*sizeof(Covariance));  /* inverse covariance matrices */
      pst->vseq--;

      for (t=1; t<=pst->T; t++) {
         if (pst->fullCov)
            pst->vseq[t].inv = CreateSTriMat(x, pst->vSize);  /* inverse covariance (precision) matrices */
         else
            pst->vseq[t].var = CreateVector (x, pst->vSize);  /* inverse variances */
      }

      /* prepare vector and matrices for a set of linear equations */
      M = (pst->fullCov) ? pst->order : 1;
      pst->C    = CreateMatrix (x, pst->T, pst->order);     /* generated parameter sequence */
      pst->g    = CreateDVector(x, M*pst->T);               /* vector for forward and backward substitution and gradient */
      pst->c    = CreateDVector(x, M*pst->T);               /* vector for generated parameter sequence */
      pst->WUM  = CreateDVector(x, M*pst->T);               /* W'*U^{-1}*W */
      pst->WUW  = CreateDMatrix(x, M*pst->T, M*pst->width); /* W'*U^{-1}*M */
      
      if (useGV) {
         /* gv mean */
         pst->gvmean = CreateVector(x, pst->vSize);
         ZeroVector(pst->gvmean);

         /* gv covariance */
         if (gvset.ckind==FULLC) {
            pst->gvcov.inv = CreateSTriMat(x, pst->vSize); 
            ZeroTriMat(pst->gvcov.inv);
         }
         else {
            pst->gvcov.var = CreateVector(x, pst->vSize);  
            ZeroVector(pst->gvcov.var);
         }
      }
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
      if (t+j<1 || t+j>T || (pst->win.coef[i][j]!=0.0 && absolute_t+j>=1 && absolute_t+j<=absolute_T && !pst->ContSpace[absolute_t+j]))
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
   Vector mseq_t,mnomial;
   Covariance vseq_t;
   LabId id;
   MLink macro;

   /* multinomial distribution for random generation */
   mnomial = CreateVector(&gstack,genInfo->maxMixes);
   
   /* absolute time */
   const int T = genInfo->tframe;

   /* initialize time counter and statistics of each PdfStream */
   for (p=stream=1; p<=genInfo->nPdfStream[0]; stream+=genInfo->nPdfStream[p++]) {
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
      
      /* get gv mean and covariance */
      if (useGV) {
         id = GetLabId("gv",TRUE);
         /* id = genInfo->label[1]->labid; */
         if ((macro=FindMacroName(&gvset,'h',id))!=NULL) {
            for (s=stream,v=1; s<stream+genInfo->nPdfStream[p]; v+=genInfo->hset->swidth[s++]) {
               mpdf = ((HLink)macro->structure)->svec[2].info->pdf[p].info->spdf.cpdf[1].mpdf;  /* currently only sigle-mixture GV pdf is supported */
               for (k=1; k<=VectorSize(mpdf->mean); k++) {
                  /* gv mean */
                  pst->gvmean[v+k-1] = mpdf->mean[k];
                  /* gv covariance */
                  switch(mpdf->ckind) {  /* store inverse variance */
                  case DIAGC:
                     pst->gvcov.var[v+k-1] = 1.0/mpdf->cov.var[k]; 
                     break;
                  case INVDIAGC:
                     pst->gvcov.var[v+k-1] = mpdf->cov.var[k]; 
                     break;
                  case FULLC:
                     /* diagonal elements */
                     pst->gvcov.inv[v+k-1][v+k-1] = mpdf->cov.inv[k][k];
                     /* off-diagonal elements */
                     for (l=1; l<k; l++)
                        pst->gvcov.inv[v+k-1][v+l-1] = mpdf->cov.inv[k][l];
                     break;
                  default: HError(9999,"SetupPdfStreams: Only DIAGC, INVDIAGC, or FULLC is supported for GV Pdf.");
                  }
               }
            }
         }
         else
            HError(9935,"Generator: Cannot find GV model %s in current list", id->name);
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
                     ZeroVector(mnomial);
                     for (m=1,max=0,weight=LZERO; m<=si->pdf[s].info->nMix; m++) {
                        mpdf = si->pdf[s].info->spdf.cpdf[m].mpdf;
                        if (VectorSize(mpdf->mean)==genInfo->hset->swidth[s]) {
                           /* use the best mixture only */
                           if (si->pdf[s].info->spdf.cpdf[m].weight >= weight) {
                              weight = si->pdf[s].info->spdf.cpdf[m].weight;
                              max = m;
                           }
                           mnomial[m] = MixWeight(genInfo->hset,si->pdf[s].info->spdf.cpdf[m].weight);
                        }
                     }
                     if (max==0)
                        HError(999,"SetupPdfStreams: no mix found");
                     if (rFlags&RNDMIX)  /* randomly select one mixture component */
                        max = MultiNomial(mnomial,si->pdf[s].info->nMix);
                        
                     /* set pdf streams */
                     mpdf = si->pdf[s].info->spdf.cpdf[max].mpdf;
                     for (k=1; k<=genInfo->hset->swidth[s]; k++) {
                        bound = ChkBoundary(pst, v+k-1, t, T);
                        if (pst->fullCov) {
                           /* full case */
                           switch (mpdf->ckind) {
                           case DIAGC:
                              mseq_t    [v+k-1]        += (bound) ? 0.0 : 1.0/mpdf->cov.var[k] * mpdf->mean[k];
                              vseq_t.inv[v+k-1][v+k-1] += (bound) ? 0.0 : 1.0/mpdf->cov.var[k];
                              break;
                           case INVDIAGC:
                              mseq_t    [v+k-1]        += (bound) ? 0.0 : mpdf->cov.var[k] * mpdf->mean[k];
                              vseq_t.inv[v+k-1][v+k-1] += (bound) ? 0.0 : mpdf->cov.var[k];
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
                              HError(9999, "SetupPdfStreams: Only DIAGC, INVDIAGC or FULLC is supported.");
                           }
                        }
                        else {
                           switch (mpdf->ckind) {
                           case DIAGC:
                              mseq_t    [v+k-1] += (bound) ? 0.0 : 1.0/mpdf->cov.var[k] * mpdf->mean[k];
                              vseq_t.var[v+k-1] += (bound) ? 0.0 : 1.0/mpdf->cov.var[k];
                              break;
                           case INVDIAGC:
                              mseq_t    [v+k-1] += (bound) ? 0.0 : mpdf->cov.var[k] * mpdf->mean[k];
                              vseq_t.var[v+k-1] += (bound) ? 0.0 : mpdf->cov.var[k];
                              break;
                           default:
                              HError(9999, "SetupPdfStreams: Only DIAGC or INVDIAGC is supported.");
                           }
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
   
   FreeVector(&gstack,mnomial);

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
   Vector mnomial;

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
      
      /* transition probs for random generation */
      mnomial = CreateVector(&gstack,genInfo->maxStates);

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
            ZeroVector(mnomial);
            for (k=2,best=1,trans=LZERO; k<=hmm->numStates; k++) {
               if (k==j) continue;  /* exclude self-transition */
               if (hmm->transP[j][k]>trans) {  /* select ML one */
                  trans = hmm->transP[j][k];
                  best  = k;
               }
               mnomial[k] = L2F(hmm->transP[j][k]);
            }
            if (rFlags&RNDTRANS) {  /* randomly select next state */
               best  = MultiNomial(mnomial, hmm->numStates);
               trans = hmm->transP[j][best];
            }
            
            /* check acyclic transition */
            if (!(rFlags&RNDTRANS) && IsMember(acyclic,best))
               HError(9999, "SetStateSeuqnce: acyclic transition is detected.");
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
      
      /* free */
      FreeVector(&gstack,mnomial);
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
   float sum,sqr,dur,rho=0.0,diff=0.0;
   Vector *mean, *ivar;
   Label *label;
   HLink dm;
   
   if (genInfo->speakRate!=1.0 && rFlags&RNDDUR)
      HError(9999,"SetStateDurations: Cannot change speaking rate in random duration generation");
   if (genInfo->modelAlign && rFlags&RNDDUR)
      HError(9999,"SetStateDurations: Cannot use model-level alignments in random duration generation");
   
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
            switch(dm->svec[2].info->pdf[s].info->spdf.cpdf[1].mpdf->ckind) {
            case DIAGC:    ivar[i][cnt] = 1.0 / dm->svec[2].info->pdf[s].info->spdf.cpdf[1].mpdf->cov.var[k]; break;
            case INVDIAGC: ivar[i][cnt] = dm->svec[2].info->pdf[s].info->spdf.cpdf[1].mpdf->cov.var[k]; break;
            case FULLC:    ivar[i][cnt] = dm->svec[2].info->pdf[s].info->spdf.cpdf[1].mpdf->cov.inv[k][k]; break;
            }
         }
      }
      
      /* acc duration statistics to set rho */
      if (genInfo->speakRate!=1.0)
         CountDurStat(mean[i], ivar[i], &sum, &sqr, genInfo->sindex[i]);
   }
   
   /* set rho, please refer to
    * T. Yoshimura, et al. "Duration Modeling in HMM-based Speech Synthesis System", 
    * Proc. of ICSLP, vol.2, pp.29-32, 1998, for detail 
    * */
   if (genInfo->speakRate!=1.0)
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
            rho = 0.0;
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
         k = genInfo->sindex[i][j]-1;
         dur = (rFlags&RNDDUR) ? GaussDeviate(mean[i][k],sqrt(1.0/ivar[i][k])) /* random duration sampling */
                               : mean[i][k]+rho/ivar[i][k];
         genInfo->durations[i][j] = (int)(dur+diff+0.5);
         
         /* set minimum duration -> 1 */
         if (genInfo->durations[i][j]<1)
            genInfo->durations[i][j] = 1;

         diff += dur-(float)genInfo->durations[i][j];
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
   int i,j,d,s,m,t,p,stream;
   float ContSpaceWeight,binomial[3];
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
         for (d=0; d<genInfo->durations[i][j]; d++) {
            for (p=stream=1; p<=genInfo->nPdfStream[0]; stream+=genInfo->nPdfStream[p++]) {
               ContSpace = FALSE;
               for (s=stream; s<stream+genInfo->nPdfStream[p]; s++) {
                  if (genInfo->hset->msdflag[s]) {
                     ContSpaceWeight = 0.0;
                     for (m=1; m<=si->pdf[s].info->nMix; m++)
                        if (VectorSize(si->pdf[s].info->spdf.cpdf[m].mpdf->mean)==genInfo->hset->swidth[s])  /* total weight of all continuous mixtures */
                           ContSpaceWeight += MixWeight(genInfo->hset, si->pdf[s].info->spdf.cpdf[m].weight);
                     
                     /* binomial distribution to randomly select continuous/discrete space */ 
                     binomial[1] = ContSpaceWeight;  
                     binomial[2] = 1.0 - ContSpaceWeight;
                     
                     /* if any streams in the p-th PdfStream is determined to continuous, this frame determined to be continous */
                     if (  ((rFlags&RNDSPACE) && MultiNomial(binomial,2)==1) ||
                          (!(rFlags&RNDSPACE) && ContSpaceWeight>genInfo->MSDthresh) ) {
                        ContSpace = TRUE;
                        break;
                     }
                  }
                  else {
                     ContSpace = TRUE;
                     break;
                  }
               }
               
               if (ContSpace) {
                  genInfo->pst[p].ContSpace[t+d] = TRUE;  /* set ContSpace flag */
                  genInfo->pst[p].T++;  /* increment total number of frame in this PdfStream */
               }
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

/* EXPORT->JointProb: joint probability of given observations and state sequence */ 
void JointProb (GenInfo *genInfo, UttInfo *utt)
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
      if (rFlags&RNDPAR)
         g[t] += GaussDeviate(rndParMean, rndParVar);
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

   Matrix  C = pst->C;
   DMatrix U = pst->WUW;
   DVector g = pst->g;
   DVector c = pst->c;

   /* sizes of matrix and vector */
   const int M = (pst->fullCov) ? pst->order : 1;
   const int T = M*pst->T;
   const int width = M*pst->width;

   if (trace & T_MAT)
      printf("\n  solution\n   ");

   /* backward substitution */
   for (t=T; t>0; t--) {
      c[t] = g[t];
      for (i=1; (i<width) && (t+i<=T); i++)
         c[t] -= U[t][i+1]*c[t+i];
      c[t] /= U[t][1];

      if (trace & T_MAT) {
         printf("%8.2lf ", c[t]);
         fflush(stdout);
      }
   }
   
   /* store generated parameters */
   for (t=1; t<=T; t++)
      C[(t+M-1)/M][(t+M-1)%M+1+bias] = (float) c[t];
   
   return;
}

/* -------------------------------- GV related functions ------------------------------------- */
/* Calc_GV: calculate GV of the current c */
static void Calc_GV (PdfStream *pst, const int bias, DVector mean, DVector var)
{
   int m,t;
   DVector c = pst->c;
   
   /* constant values */
   const int M = (pst->fullCov) ? pst->order : 1;
   const int T = M*pst->T;
   
   /* mean */
   ZeroDVector(mean);
   for (t=1; t<=T; t++) {
      m = (t+M-1)%M+1+bias;
      mean[m] += c[t];
   }
   for (m=1; m<=M; m++)
      mean[m+bias] /= pst->T;

   /* variance */
   ZeroDVector(var);
   for (t=1; t<=T; t++) {
      m = (t+M-1)%M+1+bias;
      var[m] += (c[t]-mean[m])*(c[t]-mean[m]);
   }
   for (m=1; m<=M; m++)
      var[m+bias] /= pst->T;

   return;
}

/* Conv_GV: expand c according to mean value of a given GV pdf */
static void Conv_GV (PdfStream *pst, const int bias, DVector mean, DVector var)
{
   int m,t;
   
   /* constant values */
   const int M = (pst->fullCov) ? pst->order : 1;
   const int T = M*pst->T;

   /* Vector/Matrices */
   DVector c = pst->c;
   DVector ratio = CreateDVector(&gvStack, pst->order);
      
   /* calculate GV of c */
   Calc_GV(pst, bias, mean, var);
   
   /* ratio between GV mean and variance of c */
   for (m=1; m<=M; m++)
      ratio [m+bias] = sqrt(pst->gvmean[m+bias]/var[m+bias]);
   
   /* expand c */
   for (t=1; t<=T; t++) {
      m = (t+M-1)%M+1+bias;
      c[t] = ratio[m] * (c[t]-mean[m]) + mean[m];
   }
   
   FreeDVector(&gvStack, ratio);
   
   return;
}

/* Calc_Gradient: calculate a gradient vector of the GV objective function with respect to c */
static LogDouble Calc_Gradient (PdfStream *pst, const int bias, DVector mean, DVector var, 
                                LogDouble *GVobj, LogDouble *HMMobj, double *norm)
{
   int m,l,t,i;
   double inv,h;
   
   /* constant values */
   const int M = (pst->fullCov) ? pst->order : 1;
   const int T = M*pst->T;
   const int width = M*pst->width;
   const double w = 1.0/(pst->win.num*pst->T);

   /* Vector/Matrices */
   DMatrix R = pst->WUW;
   DVector r = pst->WUM;
   DVector c = pst->c;
   DVector g = pst->g;

   /* GV pdf statistics */
   DVector vd = CreateDVector(&gvStack, pst->order);
   Vector gvmean = pst->gvmean;
   Covariance gvcov = pst->gvcov;
   Boolean fullGV = (gvset.ckind==FULLC) ? TRUE : FALSE;

   /* recalculate GV of the current c */
   Calc_GV(pst, bias, mean, var);
   
   /* GV objective function and its derivative with respect to c */
   ZeroDVector(vd);
   for (m=1,*GVobj=0.0; m<=M; m++) {
      for (l=((fullGV)?1:m); l<=((fullGV)?M:m); l++) {
         inv = (fullGV) ? ((l<m) ? gvcov.inv[m+bias][l+bias] : gvcov.inv[l+bias][m+bias]) : gvcov.var[l+bias];
         *GVobj -= 0.5*w2*(var[m+bias]-gvmean[m+bias])*inv*(var[l+bias]-gvmean[l+bias]);
         vd[m+bias] += inv*(var[l+bias]-gvmean[l+bias]);
      }
   }

   /* calculate g = R*c */
   for (t=1; t<=T; t++) {
      g[t] = R[t][1]*c[t];
      for (i=2; i<=width; i++) {
         if (t+i-1<=T)
            g[t] += R[t][i]*c[t+i-1];
         if (t-i+1>0)
            g[t] += R[t-i+1][i]*c[t-i+1];
      }
   }
      
   for (t=1,*HMMobj=0.0,*norm=0.0; t<=pst->T; t++) {
      for (m=1; m<=M; m++) {
         /* index */ 
         i = M*(t-1)+m;

         /* objective function */
         *HMMobj += -0.5*w1*w*c[i]*(g[i]-2.0*r[i]);

         /* Hessian */
         switch(optKind) {
         case NEWTON:
            /* only diagonal elements of Hessian matrix are used */
            inv = (fullGV) ? gvcov.inv[m+bias][m+bias] : gvcov.var[m+bias];
            h = -w1*w*R[i][1] - w2*2.0/(pst->T*pst->T)*((pst->T-1)*vd[m+bias]+2.0*inv*(c[i]-mean[m+bias])*(c[i]-mean[m+bias]));
            h = -1.0/h;
            break;
         case LBFGS:
         case STEEPEST:
            /* do not use hessian */
            h = 1.0;
         }

         /* gradient vector */
         g[i] = h * ( w1*w*(-g[i]+r[i]) + w2*(-2.0/pst->T*(c[i]-mean[m+bias])*vd[m+bias]) );
         
         /* norm of gradient vector */
         *norm += g[i]*g[i];
      }
   }

   /* Euclidian norm of gradient vector */
   *norm = sqrt(*norm);
   
   /* free vector */
   FreeDVector(&gvStack, vd);

   return(*HMMobj+*GVobj);
}

/* GV_ParmGen: optimize C considering global variance */
static void GV_ParmGen (PdfStream *pst, const int bias)
{
   int t,iter;
   double norm, step=stepInit;
   LogDouble obj,prev,GVobj,HMMobj;
   
   /* matrix/vectors */
   Matrix  C = pst->C;
   DVector c = pst->c;
   DVector g = pst->g;
   DVector diag,w;
   
   /* constant values */
   const int M = (pst->fullCov) ? pst->order : 1;
   const int T = M*pst->T;

   /* GV pdf statistics */
   DVector mean = CreateDVector(&gvStack, pst->order); ZeroDVector(mean);
   DVector var  = CreateDVector(&gvStack, pst->order); ZeroDVector(var);
   
   /* variables for L-BFGS */
   int dim=T, mem=LBFGSMEM, diagco=0, iprint[]={-1,0}, iflag=0;
   double f,eps=1.0e-6, xtol=1.0e-15;
   diag  = CreateDVector(&gvStack, T);
   w     = CreateDVector(&gvStack, T*(2*LBFGSMEM+1)+2*LBFGSMEM);
   ZeroDVector(diag); ZeroDVector(w);

   /* first convert c according to GV pdf and use it as the initial value */
   Conv_GV(pst, bias, mean, var);

   /* recalculate R and r */
   Calc_WUM_and_WUW(pst,bias);

   /* iteratively optimize c */
   for (iter=1; iter<=maxGVIter; iter++) {
      /* calculate GV objective and its derivative with respect to c */
      obj = Calc_Gradient(pst, bias, mean, var, &GVobj, &HMMobj, &norm);

      /* accelerate/decelerate step size */
      if (iter>1 && optKind!=LBFGS) {
         /* objective function improved -> increase step size */
         if (obj>prev) step *= stepInc;
         
         /* objective function degraded -> go back c and decrese step size */
         if (obj<prev) {
            for (t=1; t<=T; t++)  /* go back c to that at the previous iteration */
               c[t] -= step * diag[t];
            step *= stepDec;
            for (t=1; t<=T; t++)  /* gradient c */
               c[t] += step * diag[t];
            iter--; continue;
         }
      }

      if (trace & T_GV) {
         printf("   Iteration %2d: GV Obj = %e (HMM:%e GV:%e)", iter, obj, HMMobj, GVobj);
         if (iter>1)
            printf("  Change = %f", obj-prev);
         printf("\n");
         fflush(stdout);
      }
            
      /* convergence check (Euclid norm, objective function, and LBFGS report) */
      if ((optKind!=LBFGS && norm<minEucNorm) || (iter>1 && fabs(obj-prev)<GVepsilon) || (iter>1 && optKind==LBFGS && iflag==0)) { 
         if (trace&T_GV) {
            if (iter>1)
               printf("   Converged (norm=%e, change=%e).\n", norm, fabs(obj-prev));
            else
               printf("   Converged (norm=%e).\n", norm);
            fflush(stdout);
         }
         break;
      }
      
      switch (optKind) {
      case STEEPEST:  /* steepest ascent */
      case NEWTON:    /* quasi Newton (only diagonal elements of Hessian matrix are used */
         for (t=1; t<=T; t++)
           c[t] += step * g[t];
         CopyDVector(g,diag);
         break;
      case LBFGS:
         for (t=1; t<=T; t++)
            g[t] = -g[t];  /* swap sign because L-BFGS minimizes objective function */
         f = -obj;
         /* call LBFGS FORTRAN routine */
         lbfgs_(&dim, &mem, &c[1], &f, &g[1], &diagco, &diag[1], iprint, &eps, &xtol, &w[1], &iflag);
         break;
      default:
         HError(9999,"GV_ParmGen: Not supported optimization kind.");
      }

      prev = obj;
   }
   
   /* convergence check (Euclid norm, objective function, and LBFGS report) */
   if (iter>maxGVIter && trace&T_GV) { 
      printf("   Optimization stopped by reaching max # of iterations (%d).\n",maxGVIter);
      fflush(stdout);
   }

   /* store generated parameters */
   for (t=1; t<=T; t++)
      C[(t+M-1)/M][(t+M-1)%M+1+bias] = (float) c[t];

   /* free vectors allocated in this function */
   FreeDVector(&gvStack, mean);
   
   return;
}
   
/* Cholesky_ParmGen: Generate parameter sequence using Cholesky decomposition */
static void Cholesky_ParmGen (GenInfo *genInfo, const Boolean GV)
{
   int p,m;
   PdfStream *pst;

  if (GV && (trace & T_GV)) {
      char buf[MAXSTRLEN];
      printf(" Parameter generation considering global variance (GV)\n");
      printf("  Optimization=%s  ", OptKind2Str(optKind,buf));
      if (optKind!=LBFGS)
         printf("step size (init=%.3f, inc=%.2f, dec=%.2f)", stepInit, stepInc, stepDec);
      printf("\n  HMM weight=%.2f, GV weight=%.2f\n", w1, w2);
      fflush(stdout);
   }

   for (p=1,pst=genInfo->pst+1; p<=genInfo->nPdfStream[0]; p++,pst++) {
      if (pst->T<1)
         continue;

      if ((trace&T_MAT) || (GV && (trace & T_GV))) {
         printf("  Stream: %d\n", p);
         fflush(stdout);
      }

      for (m=1; m<=((pst->fullCov)?1:pst->order); m++) {
         if ((trace&T_MAT) || (GV && (trace & T_GV))) {
            if (pst->fullCov) 
               printf("  Feature: all\n");
            else
               printf("  Feature: %d\n", m);
            fflush(stdout);
         }

         /* generate m-th feature */
         Calc_WUM_and_WUW(pst, m-1);
         Cholesky_Factorization(pst);     /* Cholesky decomposition */
         Forward_Substitution(pst);       /* forward substitution   */
         Backward_Substitution(pst, m-1); /* backward substitution  */
         if (GV)
            GV_ParmGen(pst, m-1);    /* iterative optimization */
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
   int iter;
   Boolean success=TRUE,converged=FALSE;
   LogFloat prev, curr;

   /* # of frames */
   const int T = genInfo->tframe;
   
   /* UttInfo settings */   
   utt->tgtSampRate = genInfo->frameRate;
   utt->S = genInfo->hset->swidth[0];
   
   /* First perform Cholesky-based parameter generation */
   SetupPdfStreams(genInfo);
   Cholesky_ParmGen(genInfo, (((type==CHOLESKY)&&useGV)?TRUE:FALSE) );
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
   for (iter=1; iter<=maxEMIter; iter++) {
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
      curr = utt->pr;
      if (trace & T_TOP) {
         printf("  Iteration %d: Average LogP = %e", iter, curr/T);
         if (iter>1)
            printf("  Change = %f", (curr-prev)/T);
         printf("\n");
         fflush(stdout);
      }
      
      /* convergence check */
      if (iter>1 && fabs(curr-prev)/T<EMepsilon) {
         converged = TRUE;
         if (trace&T_TOP) {
            printf("  Converged (change=%e).\n", fabs(curr-prev)/T);
            fflush(stdout);
         }
      }
      if (iter==maxEMIter && trace&T_TOP) {
         printf("  EM iteration stopped by reaching max # of iterations (%d).\n",maxEMIter);
      }
      
      /* generate parameters */
      UpdatePdfStreams(genInfo, fbInfo, utt);  /* update PdfStreams */
      Cholesky_ParmGen(genInfo, (((converged||iter==maxEMIter)&&useGV)?TRUE:FALSE));
      UpdateUttObs(genInfo, utt);
      
      /* reset heap */
      ResetHeap(&fbInfo->ab->abMem);
            
      if (converged)
         break;
      
      prev = curr;
   }

   return;
}

/* ---------------------------- End of HGen.c --------------------------- */
