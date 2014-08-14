/* ----------------------------------------------------------------- */
/*           The HMM-Based Speech Synthesis System (HTS)             */
/*           developed by HTS Working Group                          */
/*           http://hts.sp.nitech.ac.jp/                             */
/* ----------------------------------------------------------------- */
/*                                                                   */
/*  Copyright (c) 2008-2010  Nagoya Institute of Technology          */
/*                           Department of Computer Science          */
/*                                                                   */
/* All rights reserved.                                              */
/*                                                                   */
/* Redistribution and use in source and binary forms, with or        */
/* without modification, are permitted provided that the following   */
/* conditions are met:                                               */
/*                                                                   */
/* - Redistributions of source code must retain the above copyright  */
/*   notice, this list of conditions and the following disclaimer.   */
/* - Redistributions in binary form must reproduce the above         */
/*   copyright notice, this list of conditions and the following     */
/*   disclaimer in the documentation and/or other materials provided */
/*   with the distribution.                                          */
/* - Neither the name of the HTS working group nor the names of its  */
/*   contributors may be used to endorse or promote products derived */
/*   from this software without specific prior written permission.   */
/*                                                                   */
/* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND            */
/* CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,       */
/* INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF          */
/* MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE          */
/* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS */
/* BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,          */
/* EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED   */
/* TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,     */
/* DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON */
/* ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,   */
/* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY    */
/* OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE           */
/* POSSIBILITY OF SUCH DAMAGE.                                       */
/* ----------------------------------------------------------------- */

char *hfst_version = "!HVER!HFst:   2.1.1 beta  [CUED 25/12/09]";
char *hfst_vc_id = "$Id: HFst.c,v 1.5 2010/04/08 04:50:30 uratec Exp $";

/*
   This program is used to make WFST from HSMMs.
*/

#include "HShell.h"
#include "HMem.h"
#include "HMath.h"
#include "HAudio.h"
#include "HWave.h"
#include "HParm.h"
#include "HLabel.h"
#include "HModel.h"
#include "HTrain.h"
#include "HUtil.h"
#include "HAdapt.h"
#include "HFB.h"

/* Trace Flags */
#define T_TOP 0001              /* Top level tracing */

/* Models */
static char *hmmDir = NULL;     /* directory to look for hmm def files */
static char *durDir = NULL;     /* directory to look for duration model def files */
static char *hmmExt = NULL;     /* hmm def file extension */
static char *durExt = NULL;     /* duration model def file extension */

/* WFST */
static char *hmmFstDir = NULL;  /* directory to save Fst for HMM */
static char *durFstDir = NULL;  /* directory to save Fst for Duration */
static int prunFrame = -1;      /* prun WFST by using time information of label */
static int fstType = 0;         /* type of WFST (0: OpenFst 1: MIT FST) */
static double fstDurWeight = 1.0;       /* duration weight for WFST of duration model */
static double fstMaxWidth = 4.0;        /* max width for WFST of HSMM */

/* Label */
static char *labDir = NULL;     /* label (transcription) file directory */
static char *labExt = "lab";    /* label file extension */

/* Global settings */
static int trace = 0;           /* Trace level */
static FileFormat dff = UNDEFF; /* data file format */
static FileFormat lff = UNDEFF; /* label file format */
static ConfParam *cParm[MAXGLOBS];      /* configuration parameters */
static int nParm = 0;           /* total num params */
Boolean keepOccm = FALSE;

/* Stack */
static MemHeap hmmStack;
static MemHeap durStack;
static MemHeap uttStack;
static MemHeap tmpStack;

/* ------------------ Process Command Line ------------------------- */

/* SetConfParms: set conf parms relevant to HFst */
void SetConfParms(void)
{
   int i;

   nParm = GetConfig("HFST", TRUE, cParm, MAXGLOBS);
   if (nParm > 0) {
      if (GetConfInt(cParm, nParm, "TRACE", &i))
         trace = i;
   }
}

void ReportUsage(void)
{
   printf("\nUSAGE: HFst [options] hmmList durList dataFiles...\n\n");
   printf
       (" Option                                                    Default\n\n");
   printf(" -c      number of frame for pruning                       none\n");
   printf
       (" -d s    dir to find hmm definitions                       current\n");
   printf(" -e N    type of output WFST (0: OpenFst 1: MIT FST)       0\n");
   printf
       (" -n s    dir to find duration model definitions            current\n");
   printf
       (" -m dir  dir to write HMM WFST                             current\n");
   printf
       (" -r dir  dir to write duration WFST                        current\n");
   printf(" -v f    max width for HSMM (mean +- f * sqrt(vari))       4.0\n");
   printf(" -w f    duration weight                                   1.0\n");
   printf(" -x s    extension for hmm files                           none\n");
   printf(" -y s    extension for duration model files                none\n");
   PrintStdOpts("ACDFGHILNSTX");
   printf("\n\n");
}

int main(int argc, char *argv[])
{
   char *datafn = NULL;
   char *hmmListFn = NULL, *durListFn = NULL;
   char *s;
   UttInfo *utt;                /* utterance information storage */
   HMMSet hset;                 /* Set of HMMs to be re-estimated */
   HMMSet dset;                 /* Set of duration models to be generated */
   int *maxMixInS;

   void Initialise(HMMSet * hset, HMMSet * dset, char *hmmListFn,
                   char *durListFn, int *maxMixInS);
   void makeFst(UttInfo * utt, char *datafn, char *outHMMDir, char *outDurDir,
                HMMSet * hset, HMMSet * dset, int *maxMixInS);

   if (InitShell(argc, argv, hfst_version, hfst_vc_id) < SUCCESS)
      HError(2300, "HFst: InitShell failed");
   InitMem();
   InitMath();
   InitWave();
   InitLabel();
   InitModel();
   if (InitParm() < SUCCESS)
      HError(2300, "HFst: InitParm failed");

   if (!InfoPrinted() && NumArgs() == 0)
      ReportUsage();
   if (NumArgs() == 0)
      Exit(0);
   CreateHeap(&hmmStack, "HmmStore", MSTAK, 1, 1.0, 50000, 500000);
   CreateHeap(&durStack, "DurStore", MSTAK, 1, 1.0, 50000, 500000);
   SetConfParms();
   CreateHMMSet(&hset, &hmmStack, TRUE);
   CreateHMMSet(&dset, &durStack, TRUE);
   CreateHeap(&uttStack, "uttStore", MSTAK, 1, 0.5, 100, 1000);
   utt = (UttInfo *) New(&uttStack, sizeof(UttInfo));
   CreateHeap(&tmpStack, "tmpStore", MSTAK, 1, 0.5, 100, 1000);

   while (NextArg() == SWITCHARG) {
      s = GetSwtArg();
      if (strlen(s) != 1)
         HError(2319, "HFst: bad switch %s; must be single letter", s);
      switch (s[0]) {
      case 'c':
         prunFrame = GetIntArg();
         break;
      case 'd':
         if (NextArg() != STRINGARG)
            HError(2319, "HFst: HMM definition directory expected");
         hmmDir = GetStrArg();
         break;
      case 'e':
         fstType = GetChkedInt(0, 1, s);
         break;
      case 'm':
         if (NextArg() != STRINGARG)
            HError(2319, "HFst: output WFST file directory expected");
         hmmFstDir = GetStrArg();
         break;
      case 'n':
         if (NextArg() != STRINGARG)
            HError(2319, "HFst: duration model definition directory expected");
         durDir = GetStrArg();
         break;
      case 'r':
         if (NextArg() != STRINGARG)
            HError(2319, "HFst: output duration WFST file directory expected");
         durFstDir = GetStrArg();
         break;
      case 'v':
         fstMaxWidth = GetChkedFlt(0.0, 100000.0, s);
         break;
      case 'w':
         fstDurWeight = GetChkedFlt(0.0, 100000.0, s);
         break;
      case 'x':
         if (NextArg() != STRINGARG)
            HError(2319, "HFst: HMM file extension expected");
         hmmExt = GetStrArg();
         break;
      case 'y':
         if (NextArg() != STRINGARG)
            HError(2319, "HFst: duration model file extension expected");
         durExt = GetStrArg();
         break;
      case 'F':
         if (NextArg() != STRINGARG)
            HError(2319, "HFst: Data File format expected");
         if ((dff = Str2Format(GetStrArg())) == ALIEN)
            HError(-2389, "HFst: Warning ALIEN Data file format set");
         break;
      case 'G':
         if (NextArg() != STRINGARG)
            HError(2319, "HFst: Label File format expected");
         if ((lff = Str2Format(GetStrArg())) == ALIEN)
            HError(-2389, "HFst: Warning ALIEN Label file format set");
         break;
      case 'H':
         if (NextArg() != STRINGARG)
            HError(2319, "HFst: HMM macro file name expected");
         AddMMF(&hset, GetStrArg());
         break;
      case 'I':
         if (NextArg() != STRINGARG)
            HError(2319, "HFst: MLF file name expected");
         LoadMasterFile(GetStrArg());
         break;
      case 'L':
         if (NextArg() != STRINGARG)
            HError(2319, "HFst: Label file directory expected");
         labDir = GetStrArg();
         break;
      case 'N':
         if (NextArg() != STRINGARG)
            HError(2319, "HFst: duration model macro file name expected");
         AddMMF(&dset, GetStrArg());
         break;
      case 'T':
         trace = GetChkedInt(0, 0100000, s);
         break;
      case 'X':
         if (NextArg() != STRINGARG)
            HError(2319, "HFst: Label file extension expected");
         labExt = GetStrArg();
         break;
      default:
         HError(2319, "HFst: Unknown switch %s", s);
      }
   }
   if (NextArg() != STRINGARG)
      HError(2319, "HFst: file name of vocabulary list expected");

   hmmListFn = GetStrArg();
   durListFn = GetStrArg();
   maxMixInS = CreateIntVec(&tmpStack, hset.swidth[0]);
   Initialise(&hset, &dset, hmmListFn, durListFn, maxMixInS);
   InitUttInfo(utt, FALSE);

   do {
      if (NextArg() != STRINGARG)
         HError(2319, "HFst: data file name expected");
      datafn = GetStrArg();
      makeFst(utt, datafn, hmmFstDir, durFstDir, &hset, &dset, maxMixInS);
   } while (NumArgs() > 0);

   ResetHeap(&tmpStack);
   Dispose(&uttStack, utt);
   ResetHeap(&uttStack);
   ResetHeap(&durStack);
   ResetHeap(&hmmStack);

   ResetParm();
   ResetModel();
   ResetLabel();
   ResetWave();
   ResetMath();
   ResetMem();
   ResetShell();

   Exit(0);
   return (0);
}

/* -------------------------- Initialization ----------------------- */

void Initialise(HMMSet * hset, HMMSet * dset, char *hmmListFn, char *durListFn,
                int *maxMixInS)
{
   int L, P, S, vSize, s;

   /* HMM initialization */
   if (MakeHMMSet(hset, hmmListFn) < SUCCESS)
      HError(2321, "Initialise: MakeHMMSet failed");
   if (LoadHMMSet(hset, hmmDir, hmmExt) < SUCCESS)
      HError(2321, "Initialise: LoadHMMSet failed");
   P = hset->numPhyHMM;
   L = hset->numLogHMM;
   vSize = hset->vecSize;
   S = hset->swidth[0];
   ConvDiagC(hset, TRUE);

   if (trace & T_TOP) {
      printf("HMM:\n");
      printf("System is ");
      switch (hset->hsKind) {
      case PLAINHS:
         printf("PLAIN\n");
         break;
      case SHAREDHS:
         printf("SHARED\n");
         break;
      case TIEDHS:
         printf("TIED\n");
         break;
      case DISCRETEHS:
         printf("DISCRETE\n");
         break;
      }
      printf(" %d Logical/%d Physical Models Loaded, VecSize=%d\n", L, P,
             vSize);
      if (hset->numFiles > 0)
         printf(" %d MMF input files\n", hset->numFiles);
      fflush(stdout);
   }

   /* duration model initialization */
   if (MakeHMMSet(dset, durListFn) < SUCCESS)
      HError(2321, "Initialise: MakeHMMSet failed");
   if (LoadHMMSet(dset, durDir, durExt) < SUCCESS)
      HError(2321, "Initialise: LoadHMMSet failed");
   ConvDiagC(dset, TRUE);

   if (trace & T_TOP) {
      printf("Duration model:\n");
      printf("System is ");
      switch (dset->hsKind) {
      case PLAINHS:
         printf("PLAIN\n");
         break;
      case SHAREDHS:
         printf("SHARED\n");
         break;
      case TIEDHS:
         printf("TIED\n");
         break;
      case DISCRETEHS:
         printf("DISCRETE\n");
         break;
      }
      printf(" %d Logical/%d Physical Models Loaded, #States=%d\n",
             dset->numLogHMM, dset->numPhyHMM, dset->vecSize);
      if (dset->numFiles > 0)
         printf(" %d MMF input files\n", dset->numFiles);
      fflush(stdout);
   }

   for (s = 1; s <= S; s++)
      maxMixInS[s] = MaxMixInSetS(hset, s);
   ConvLogWt(hset);

   if (trace & T_TOP) {
      printf("\n");
      fflush(stdout);
   }
}

/* ------------------------- FST Generation ------------------------ */

void makeFst(UttInfo * utt, char *datafn, char *outHMMDir, char *outDurDir,
             HMMSet * hset, HMMSet * dset, int *maxMixInS)
{
   char labfn[MAXFNAMELEN];
   char fstfn[MAXFNAMELEN];
   char basefn[MAXFNAMELEN];
   char namefn[MAXFNAMELEN];
   int i, j, l, t, idx;
   FILE *fp;
   char *name;
   LLink llink;
   MLink mlink;
   HLink hlink;
   MixPDF *pdf;
   LogFloat p;
   double mean;
   double vari;
   int dmin;
   int dmax;
   Boolean isPipe;
   int numStates = dset->swidth[0];

   /* load utterance */
   utt->twoDataFiles = FALSE;
   utt->S = hset->swidth[0];
   strcpy(labfn, datafn);
   LoadLabs(utt, lff, labfn, labDir, labExt);
   LoadData(hset, utt, dff, datafn, NULL);
   InitUttObservations(utt, hset, datafn, maxMixInS);
   BaseOf(datafn, basefn);

   if (trace & T_TOP) {
      printf(" Processing Data: %s ;", NameOf(datafn, namefn));
      printf(" Label %s.%s\n", basefn, labExt);
      fflush(stdout);
   }

   /* save WFST for parameter */
   sprintf(fstfn, "%s%c%s.fst", outHMMDir, PATHCHAR, basefn);
   if ((fp = FOpen(fstfn, NoOFilter, &isPipe)) == NULL)
      HError(2611, "HFst: Cannot open FST file %s", fstfn);
   if (fstType == 1) {
      fprintf(fp, "#FSTBasic MinPlus\n");
      fprintf(fp, "I 0\n");
      fprintf(fp, "F %d\n", utt->T);
   }
   for (t = 1; t <= utt->T; t++) {
      /* probability calculation */
      for (llink = utt->tr->head->head, l = 0; llink != NULL;
           llink = llink->succ) {
         if (llink->labid != NULL) {
            name = llink->labid->name;
            i = (llink->start * (1.0 / utt->tgtSampRate) + 1);
            j = (llink->end * (1.0 / utt->tgtSampRate) + 1);
            if (llink->start < 0.0 || llink->end < 0.0 || prunFrame < 0 ||
                (i - prunFrame <= t && t <= j + prunFrame)) {
               mlink = FindMacroName(hset, 'l', llink->labid);
               hlink = (HLink) mlink->structure;
               for (i = 2; i < hlink->numStates; i++) {
                  p = (-1.0) * POutP(hset, &(utt->o[t]), hlink->svec[i].info);
                  if (fstType == 0)
                     fprintf(fp, "%d %d %s_f%d %s_m%d_s%d %f\n", t - 1, t, name,
                             t - 1, name, l, i, p);
                  else if (fstType == 1)
                     fprintf(fp, "T %d %d %s_f%d %s_m%d_s%d %f\n", t - 1, t,
                             name, t - 1, name, l, i, p);
               }
            }
            l++;
         }
      }
   }
   if (fstType == 0)
      fprintf(fp, "%d 0.0\n", utt->T);
   FClose(fp, isPipe);

   /* save WFST for duration */
   sprintf(fstfn, "%s%c%s.fst", outDurDir, PATHCHAR, basefn);
   if ((fp = FOpen(fstfn, NoOFilter, &isPipe)) == NULL)
      HError(2611, "HFst: Cannot open FST file %s", fstfn);
   idx = utt->Q * numStates + 1;
   if (fstType == 1) {
      fprintf(fp, "#FSTBasic MinPlus\n");
      fprintf(fp, "I 0\n");
      fprintf(fp, "F %d\n", utt->Q * numStates);
   }
   for (l = 0, llink = utt->tr->head->head; llink != NULL; llink = llink->succ) {
      if (llink->labid != NULL) {
         name = llink->labid->name;
         mlink = FindMacroName(dset, 'l', llink->labid);
         hlink = (HLink) mlink->structure;
         for (i = 1; i <= numStates; i++) {
            pdf = hlink->svec[2].info->pdf[i].info->spdf.cpdf[1].mpdf;
            mean = pdf->mean[1];
            vari = 1.0 / pdf->cov.var[1];
            dmin = mean - sqrt(vari) * fstMaxWidth;
            dmax = mean + sqrt(vari) * fstMaxWidth;
            if (dmin < 1)
               dmin = 1;
            if (dmax < dmin)
               dmax = dmin;
            for (j = 1; j <= dmax; j++) {
               p = fstDurWeight * (-1.0) * -0.5 *
                   (log(2 * PI * vari) +
                    (j - mean) * (j - mean) * (1.0 / vari));
               if (fstType == 0) {
                  if (j == 1) {
                     if (j >= dmin)
                        fprintf(fp, "%d %d %s_m%d_s%d %s_m%d_s%d %f\n",
                                l * numStates + i - 1, l * numStates + i, name,
                                l, i + 1, name, l, i + 1, p);
                  } else if (j == 2) {
                     fprintf(fp, "%d %d %s_m%d_s%d %s_m%d_s%d %f\n",
                             l * numStates + i - 1, idx, name, l, i + 1, name,
                             l, i + 1, 0.0);
                     if (j >= dmin)
                        fprintf(fp, "%d %d %s_m%d_s%d %s_m%d_s%d %f\n", idx,
                                l * numStates + i, name, l, i + 1, name, l,
                                i + 1, p);
                     idx++;
                  } else {
                     fprintf(fp, "%d %d %s_m%d_s%d %s_m%d_s%d %f\n", idx - 1,
                             idx, name, l, i + 1, name, l, i + 1, 0.0);
                     if (j >= dmin)
                        fprintf(fp, "%d %d %s_m%d_s%d %s_m%d_s%d %f\n", idx,
                                l * numStates + i, name, l, i + 1, name, l,
                                i + 1, p);
                     idx++;
                  }
               } else if (fstType == 1) {
                  if (j == 1) {
                     if (j >= dmin)
                        fprintf(fp, "T %d %d %s_m%d_s%d %s_m%d_s%d %f\n",
                                l * numStates + i - 1, l * numStates + i, name,
                                l, i + 1, name, l, i + 1, p);
                  } else if (j == 2) {
                     fprintf(fp, "T %d %d %s_m%d_s%d %s_m%d_s%d %f\n",
                             l * numStates + i - 1, idx, name, l, i + 1, name,
                             l, i + 1, 0.0);
                     if (j >= dmin)
                        fprintf(fp, "T %d %d %s_m%d_s%d %s_m%d_s%d %f\n", idx,
                                l * numStates + i, name, l, i + 1, name, l,
                                i + 1, p);
                     idx++;
                  } else {
                     fprintf(fp, "T %d %d %s_m%d_s%d %s_m%d_s%d %f\n", idx - 1,
                             idx, name, l, i + 1, name, l, i + 1, 0.0);
                     if (j >= dmin)
                        fprintf(fp, "T %d %d %s_m%d_s%d %s_m%d_s%d %f\n", idx,
                                l * numStates + i, name, l, i + 1, name, l,
                                i + 1, p);
                     idx++;
                  }
               }
            }
         }
         l++;
      }
   }
   if (fstType == 0)
      fprintf(fp, "%d 0.0\n", utt->Q * numStates);
   FClose(fp, isPipe);

   /* reset utterance */
   ResetUttObservations(utt, hset);
}

/* ----------------------------------------------------------------- */
/*                           END:  HFst.c                            */
/* ----------------------------------------------------------------- */
