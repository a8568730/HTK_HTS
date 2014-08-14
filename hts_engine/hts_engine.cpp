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
/*    hts_engine.c: a compact HMM-based speech synthesis engine      */
/*  ---------------------------------------------------------------  */

/*  hts_engine libraries */
#include "misc.hpp"
#include "tree.hpp"
#include "model.hpp"
#include "global.hpp"
#include "vocoder.hpp"
#include "mlpg.hpp"

int main(int argc, char **argv)
{
   FILE *labfp=stdin;
   FILE *lf0fp=NULL, *mcepfp=NULL, *durfp=NULL, *rawfp=NULL, *tracefp=NULL;
   double f;
   int i;
   
   void HTS_Process (FILE *, FILE *, FILE *, FILE *, FILE *, FILE *, 
                     PStream *,PStream *, 
                     globalP *,ModelSet *,TreeSet *,VocoderSetup *);
   
   ModelSet     ms;
   TreeSet      ts;
   PStream      mceppst, lf0pst;
   globalP      gp;
   VocoderSetup vs;
   
   /* default value for control parameter */
   gp.RATE     = 16000;
   gp.FPERIOD  = 80;   
   gp.RHO      = 0.0; 
   gp.ALPHA    = 0.42;
   gp.F0_STD   = 1.0;
   gp.F0_MEAN  = 0.0;
   gp.UV       = 0.5;
   gp.BETA     = 0.0;
   gp.LENGTH   = 0.0;
   gp.algnst   = 0;
   gp.algnph   = 0;
   
   /* initialise TreeSet and ModelSet */
   InitTreeSet (&ts);
   InitModelSet(&ms);
   
   /* delta window handler for log f0 */
   lf0pst.dw.fn = (char **) HTS_Calloc(argc, sizeof (char *));
   lf0pst.dw.num = 0;
   
   /* delta window handler for mel-cepstrum */
   mceppst.dw.fn = (char **) HTS_Calloc(argc, sizeof (char *));
   mceppst.dw.num = 0;

  
   /* parse command line */
   if (argc==1)
      HTS_Usage();
   
   while (--argc)
      if (**++argv == '-') {
         switch (*(*argv+1)) {
            case 'v': 
               switch (*(*argv+2)) {
                  case 's': gp.algnst = 1;  break;
                  case 'p': gp.algnph = 1;  break;
                  default:  HTS_Error(1, "hts_engine: Invalid option '-v%c'.\n", *(*argv+2));
               }
               break;
            case 't':
               switch (*(*argv+2)) {
                  case 'd': ts.fp[DUR] = HTS_Getfp(*(argv+1), "r");  break;
                  case 'f': 
                  case 'p': ts.fp[LF0] = HTS_Getfp(*(argv+1), "r");  break;
                  case 'm': ts.fp[MCP] = HTS_Getfp(*(argv+1), "r");  break;
                  default:  HTS_Error(1, "hts_engine: Invalid option '-t%c'.\n", *(*argv+2));
               }
               ++argv; --argc;
               break;
            case 'm':
               switch (*(*argv+2)) {
                  case 'd': ms.fp[DUR] = HTS_Getfp(*(argv+1), "rb");  break;
                  case 'f': 
                  case 'p': ms.fp[LF0] = HTS_Getfp(*(argv+1), "rb");  break;
                  case 'm': ms.fp[MCP] = HTS_Getfp(*(argv+1), "rb");  break;
                  default:  HTS_Error(1, "hts_engine: Invalid option '-m%c'.\n", *(*argv+2));
               }
               ++argv; --argc;
               break;
            case 'd':
               switch (*(*argv+2)) {
                  case 'm': mceppst.dw.fn[mceppst.dw.num] = *(argv+1);
                             mceppst.dw.num++;
                             break;
                  case 'f':
                  case 'p': lf0pst.dw.fn[lf0pst.dw.num] = *(argv+1);
                             lf0pst.dw.num++;
                             break;
                  default:  HTS_Error(1, "hts_engine: Invalid option '-d%c'.\n", *(*argv+2)); 
               }
               ++argv; --argc;
               break;
            case 'o':
               switch (*(*argv+2)) {
                  case 'r': rawfp   = HTS_Getfp(*(argv+1), "wb");  break;
                  case 'f': 
                  case 'p': lf0fp   = HTS_Getfp(*(argv+1), "wb");  break;
                  case 'm': mcepfp  = HTS_Getfp(*(argv+1), "wb");  break;
                  case 'd': durfp   = HTS_Getfp(*(argv+1), "wt");  break;
                  case 't': tracefp = HTS_Getfp(*(argv+1), "wt");  break;
                  default:  HTS_Error(1, "hts_engine: Invalid option '-o%c'.\n", *(*argv+2)); 
               }
               ++argv; --argc;
               break;
            case 'h': HTS_Usage(); break;
            case 's':
               i = atoi(*++argv);
               if (i>0 && i<=48000) gp.RATE = i;
               --argc;
               break;
            case 'p':
               i = atoi(*++argv);
               if (i>0 && i<=2000) gp.FPERIOD = i;
               --argc;
               break; 
            case 'a':
               f = atof(*++argv);
               if (f<=1.0 && f>=0.0) gp.ALPHA = f; 
               --argc;
               break;
            case 'b':
               f = atof(*++argv);
               if (f<=0.8 && f>=-0.8) gp.BETA = f; 
               --argc;
               break;
            case 'r':
               f = atof(*++argv);
               if (f<=1.0 && f>=-1.0) gp.RHO = f;
               --argc;
               break;
            case 'f':
               switch (*(*argv+2)) {
                  case 's': f = atof(*++argv);  
                             if (f<=5.0 && f>=0.0) gp.F0_STD=f; break;
                  case 'm': f = atof(*++argv);  
                             if (f<=100.0 && f>=0.0) gp.F0_MEAN = f;  break;
                  default:  HTS_Error(1, "hts_engine: Invalid option '-f%c'.\n", *(*argv+2)); 
               }
               --argc;
               break;
            case 'u':
               f = atof(*++argv);
               if (f<=1.0 && f>=0.0) 
                  gp.UV = f;
               --argc;
               break;
            case 'l':
               f = atof(*++argv);
               if (f<=30.0 && f>=0.0) gp.LENGTH = f;
               --argc;
               break;
            default:
               HTS_Error(1, "hts_engine: Invalid option '-%c'.\n", *(*argv+1));
         }
      }
      else 
         labfp = HTS_Getfp(*argv, "r");

   /* check command line */
   if (ts.fp[DUR] == NULL)
      HTS_Error(1, "hts_engine: file for duration trees is not specified.\n");
   if (ts.fp[LF0] == NULL)
      HTS_Error(1, "hts_engine: file for log F0 trees is not specified.\n");
   if (ts.fp[MCP] == NULL)
      HTS_Error(1, "hts_engine: file for mcep trees is not specified.\n");
   if (ms.fp[DUR] == NULL)
      HTS_Error(1, "hts_engine: file for duration pdfs is not specified.\n");
   if (ms.fp[MCP] == NULL)
      HTS_Error(1, "hts_engine: file for mcep pdfs is not specified.\n");
   if (ms.fp[LF0] == NULL)
      HTS_Error(1, "hts_engine: file for log F0 pdfs is not specified.\n");
   if (gp.algnst && gp.algnph)
      HTS_Error(1, "hts_engine: options '-vs' and '-vp' are exclusive.\n");
   if (gp.LENGTH>0.0 && gp.RHO!=0.0)
      HTS_Error(1, "hts_engine: options '-r' and '-l' are exclusive.\n");

   /* load decision trees for duration, log f0 and mel-cepstrum */
   LoadTreeSet(&ts, DUR);
   LoadTreeSet(&ts, LF0);
   LoadTreeSet(&ts, MCP);

   /* load model files for duration, log f0 and mel-cepstrum */
   LoadModelSet(&ms);

   /* check model */
   if (lf0pst.dw.num != ms.lf0stream)
      HTS_Error(1, "hts_engine: #window for log f0 is not matched to acoustic model.\n");
   if (ms.mcepvsize % mceppst.dw.num != 0 )
      HTS_Error(1, "hts_engine: #window for mcep is not matched to acoustic model.\n");


   /* initialise MLSA filter */
   InitVocoder(rawfp,ms.mcepvsize-1, &vs, gp.RATE, gp.FPERIOD);

   
   /* generate speech */
   HTS_Process(labfp, rawfp, lf0fp, mcepfp, durfp, tracefp, &mceppst, &lf0pst, &gp, &ms, &ts, &vs);


   /* free memory */
   ClearVocoder(rawfp,&vs);
   ClearModelSet(&ms);
   ClearTreeSet(&ts,MCP);
   ClearTreeSet(&ts,LF0);
   ClearTreeSet(&ts,DUR);

   HTS_Free(mceppst.dw.fn);
   HTS_Free(lf0pst.dw.fn);

   fclose(labfp);

   return 0;
}

/* OutLabel: output label with frame number or time */
void OutLabel (FILE *durfp, UttModel *um, globalP *gp) 
{
   Model *m;
   int i=0, j;
   
   for (m=um->mhead; m!=um->mtail; m=m->next) {
      j = m->totaldur;
      /* in HTK & HTS format */
      fprintf(durfp, "%d %d %s\n", (int)(i*gp->FPERIOD*1e+7/gp->RATE), 
                                   (int)((i+j)*gp->FPERIOD*1e+7/gp->RATE), m->name);
      i += j;
   }
   
   return;
}

/* OutInfo: output trace information */
void OutInfo (FILE *tracefp, UttModel *um, int nstate, globalP *gp)
{
   Model *m;
   int i, j, t;

   /* output settings */
   fprintf(tracefp, "sampring frequency                     -> %d(Hz)\n", gp->RATE);
   fprintf(tracefp, "frame period                           -> %d(point) %.2f(msec)\n", gp->FPERIOD, 1e+3*gp->FPERIOD/gp->RATE);
   fprintf(tracefp, "use state alignment for duration       -> %d\n", gp->algnst);
   fprintf(tracefp, "use phoneme alignment for duration     -> %d\n", gp->algnph);
   fprintf(tracefp, "all-pass constant                      -> %f\n", gp->ALPHA);
   fprintf(tracefp, "postfiltering coefficient              -> %f\n", gp->BETA);
   fprintf(tracefp, "control duration parameter             -> %f\n", gp->RHO);
   fprintf(tracefp, "multilply f0                           -> %f\n", gp->F0_STD);
   fprintf(tracefp, "add f0                                 -> %f\n", gp->F0_MEAN);
   fprintf(tracefp, "voiced/unvoiced threshold              -> %f\n", gp->UV);
   fprintf(tracefp, "specified utterance length             -> %f(sec.)\n\n", gp->LENGTH);
   
   /* output sentence HMM and generated utterance information */
   fprintf(tracefp, "number of HMMs        -> %d\n",um->nModel);
   fprintf(tracefp, "number of HMM states  -> %d\n",um->nState);
   fprintf(tracefp, "length of this speech -> %1.3lf sec. (%d frames)\n", (double)(um->totalframe*gp->FPERIOD)/gp->RATE,um->totalframe);
   fprintf(tracefp, "\n");
 
   /* output each state information */
   for (m=um->mhead,j=1,t=0; m!=um->mtail; m=m->next,j++) {
      fprintf(tracefp, "%d: %s \n", j, m->name);  /* model name */
      fprintf(tracefp, "            duration -> %-3d\n",m->durpdf);  /* model duration */
      for (i=2;i<nstate+2;i++) {
         /* leaf number */         
         fprintf(tracefp, " %2d-state : spectrum -> %-3d   f0 -> %-3d   ",i,m->mceppdf[i],m->lf0pdf[i]);
         /* state duration (time) */ 
         fprintf(tracefp, " %1.3lf--%1.3lf(sec)",(double)(t*gp->FPERIOD)/gp->RATE,(double)((t+m->dur[i])*gp->FPERIOD)/gp->RATE);
         /* state duration (#frame) */
         fprintf(tracefp, " %3d(frame)",m->dur[i]);
         /* voiced/unvoiced */
         fprintf(tracefp, "   %s\n", ( (m->voiced[i]) ? "voiced" : "unvoiced" ) );
         t += m->dur[i];
      }
   }
   
   return;
}

/* HTS_Process: parse label, determine state duration, generate sequence of speech parameter vector, and synthesize waveform */
void HTS_Process (FILE *labfp, FILE *rawfp, FILE *lf0fp, FILE *mcepfp, FILE *durfp, FILE *tracefp,
                   PStream *mceppst, PStream *lf0pst, 
                   globalP *gp, ModelSet *ms, TreeSet *ts, VocoderSetup *vs) 
{
   char buf[HTS_MAXBUFLEN];
   Tree *tree;
   int state;
   int start, end;
   int rate, nf;
   double f, mean, var, diffdur;
   HTS_Boolean hastime;
   Model *m,*next=NULL;
   UttModel um;
   
   rate = gp->FPERIOD*10000000/gp->RATE;
   
   mean = var = diffdur = 0.0;

   m = (Model *) HTS_Calloc(1, sizeof (Model));
   um.mtail = um.mhead = m;  
   um.totalframe = um.nState = um.nModel = 0;
   start = end = 0;
   
   /* parse label file */
   while (!feof(labfp)) {
      HTS_GetToken (labfp,buf);
      if (!isalnum(buf[0])) break;
      if (isdigit(buf[0]))
         hastime = 1;   /* label contain segmentation information */
      else 
         hastime = 0;
      
      if (hastime) {              
         if (gp->algnst) {         /* state-level segmentation */
            start = atoi(buf);     /* start time */
            HTS_GetToken(labfp, buf);  
            end = atoi(buf);       /* end time */
            HTS_GetToken(labfp, buf); 
            HTS_GetToken(labfp, buf);
         }
         else if (gp->algnph) {    /* phoneme-level segmentation */
            start = atoi(buf);     /* start time */
            HTS_GetToken(labfp, buf);
            end = atoi(buf);       /* end time */
            HTS_GetToken(labfp, buf); 
         }
         else {                    /* not use segmentation */ 
            do {
               HTS_GetToken(labfp, buf);
            } while (isdigit(buf[0]));
         }
      }
      
      /* allocate memory for current model */
      m->totaldur     = 0;
      m->name         = HTS_Strdup(buf);                                                /* model name          */
      m->dur          = (int         *) HTS_Calloc(ms->nstate+2, sizeof(int)        );  /* duration            */
      m->lf0pdf       = (int         *) HTS_Calloc(ms->nstate+2, sizeof(int)        );  /* f0 (leaf number)    */
      m->lf0mean      = (double     **) HTS_Calloc(ms->nstate+2, sizeof(double *)   );  /* f0 (mean)           */
      m->lf0variance  = (double     **) HTS_Calloc(ms->nstate+2, sizeof(double *)   );  /* f0 (variance)       */
      m->voiced       = (HTS_Boolean *) HTS_Calloc(ms->nstate+2, sizeof(HTS_Boolean));  /* voiced/unvoiced     */
      m->mceppdf      = (int         *) HTS_Calloc(ms->nstate+2, sizeof(int)        );  /* spectrum (leaf num) */
      m->mcepmean     = (double     **) HTS_Calloc(ms->nstate+2, sizeof(double *)   );  /* spectrum (mean)     */
      m->mcepvariance = (double     **) HTS_Calloc(ms->nstate+2, sizeof(double *)   );  /* spectrum (variance) */
      
      m->dur     -= 2;  
      m->lf0pdf  -= 2;  m->lf0mean  -= 2;  m->lf0variance  -= 2;  m->voiced -= 2;
      m->mceppdf -= 2;  m->mcepmean -= 2;  m->mcepvariance -= 2;  

      for (state=2; state<=ms->nstate+1; state++) {
         m->lf0mean[state]     = (double *) HTS_Calloc(ms->lf0stream, sizeof(double));
         m->lf0variance[state] = (double *) HTS_Calloc(ms->lf0stream, sizeof(double));
         m->lf0mean[state]--;  m->lf0variance[state]--;
      }
      
      /* determine state-level duration */
      if (hastime && gp->algnph) {   /* use phoneme-level segmentation */
         m->durpdf = SearchTree(m->name, ts->thead[DUR]->root);
         FindDurPDF(m, ms, gp->RHO, &diffdur);
         nf = 0;
         for (state=2; state<=ms->nstate+1; state++)
            nf += m->dur[state];
           
         fprintf(stderr, ">>>nf=%d %d\n", nf, (end-start)/rate);
         
         f = (double)(end-start)/(rate*nf);
         m->totaldur = 0;
         
         for (state=2; state<=ms->nstate+1; state++) {
            nf = (int)(f*m->dur[state]+0.5);
            if (nf<=0)  nf=1;
            fprintf(stderr, "%d: %d %f %d\n", state, m->dur[state], f, nf);
            m->dur[state] = nf;
            m->totaldur += m->dur[state];
         }
         um.totalframe += m->totaldur;
      }
      else if (hastime && gp->algnst) { /* use state-level segmentation */
         m->dur[2] = (end-start)/rate;
         m->totaldur = m->dur[2];
         um.totalframe += m->dur[2];
         
         for (state=3; state<=ms->nstate+1; state++) {
            HTS_GetToken(labfp, buf);
            start = atoi(buf);
            HTS_GetToken(labfp, buf); 
            end = atoi(buf);
            HTS_GetToken(labfp, buf);
            m->dur[state] = (end-start)/rate;
            m->totaldur += m->dur[state];
            um.totalframe  += m->dur[state];
         }
      } 
      else {  /* estimate state duration from state duration model (Gaussian) */
         m->durpdf = SearchTree(m->name, ts->thead[DUR]->root);   
         if (gp->LENGTH==0.0) {
            FindDurPDF(m, ms, gp->RHO, &diffdur);
            um.totalframe += m->totaldur;
         }
         else {   /* if total length of generated speech is specified */
            for (state=2; state<=ms->nstate+1; state++) {
               mean += ms->durpdf[m->durpdf][state];
               var  += ms->durpdf[m->durpdf][state+ms->nstate];
            }
         }
      }
     
      /* find pdf for f0 */ 
      for (tree=ts->thead[LF0],state=2; tree!=ts->ttail[LF0]; tree=tree->next,state++) {
         m->lf0pdf[state] = SearchTree(m->name, tree->root);
         FindLF0PDF(state, m, ms, gp->UV);
      }

      /* find pdf for spectrum */
      for (tree=ts->thead[MCP],state=2; tree!=ts->ttail[MCP]; tree=tree->next,state++) {
         m->mceppdf[state] = SearchTree(m->name, tree->root);
         FindMcpPDF(state, m, ms);
      }
      
      m->next = (Model *) HTS_Calloc(1, sizeof(Model));
      m = um.mtail = m->next;
      
      um.nModel++;
      um.nState+=ms->nstate;
   }
   
   /* Specified utterance length is too short */
   if (gp->LENGTH>0.0 && gp->LENGTH*gp->RATE/gp->FPERIOD<um.nState)
      HTS_Error(1, "hts_engine: specified utterance length is too short.\n");
   
   /* if total length of utterance is specified, RHO (temporal factor) have to be computed */
   if (gp->LENGTH>0.0) {
      gp->RHO = (((int)gp->LENGTH*gp->RATE/gp->FPERIOD) - mean)/var;
      /* compute state duration for each state */
      for (m=um.mhead; m!=um.mtail; m=m->next) {
         FindDurPDF(m, ms, gp->RHO, &diffdur);
         um.totalframe += m->totaldur;
      }
   }
  
   /* output trace information */
   if (tracefp!=NULL)
      OutInfo(tracefp, &um, ms->nstate, gp);
 
   /* output segment information */
   if (durfp!=NULL)
      OutLabel(durfp, &um, gp);
   
   /* generate speech parameter vector sequence and synthesize speech waveform */
   pdf2speech(rawfp, lf0fp, mcepfp, mceppst, lf0pst, gp, ms, &um, vs);
   
   /* free memory */
   for (m=um.mhead; m!=um.mtail; m=next) {
      next = m->next;
      for (state=ms->nstate+1; state>=2; state--) {
         if (m->lf0variance!=NULL) HTS_Free(++m->lf0variance[state]);
         if (m->lf0mean!=NULL)     HTS_Free(++m->lf0mean[state]);
      }
      
      m->dur     += 2;  
      m->lf0pdf  += 2;  m->lf0mean  += 2;  m->lf0variance  += 2;  m->voiced += 2;
      m->mceppdf += 2;  m->mcepmean += 2;  m->mcepvariance += 2;  
      
      HTS_Free(m->mcepvariance);
      HTS_Free(m->mcepmean);
      HTS_Free(m->mceppdf);
      HTS_Free(m->voiced);
      HTS_Free(m->lf0variance);
      HTS_Free(m->lf0mean);
      HTS_Free(m->lf0pdf);
      HTS_Free(m->dur);
      HTS_Free(m->name);
      HTS_Free(m);
   }
   HTS_Free(next);
   
   return;
}
      
/* -------------------- End of "hts_engine.cc" -------------------- */
