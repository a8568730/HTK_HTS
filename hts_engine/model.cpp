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
/*    model.c: read model and search pdf from models                 */
/*  ---------------------------------------------------------------  */

/* hts_engine libraries */
#include "misc.hpp"
#include "model.hpp"

/* LoadModelSet: load models */
void LoadModelSet (ModelSet *ms)
{
   int i, j, k, l;
   double vw, uvw;
   float temp;

   /*-------------------- load pdfs for duration --------------------*/
   /* read the number of states & the number of pdfs (leaf nodes) */

   /* read the number of HMM states */
   HTS_Fread(&ms->nstate,  sizeof(int), 1, ms->fp[DUR]);
   if (ms->nstate<0)
      HTS_Error(1, "LoadModelFiles: #HMM states must be positive value.\n");
   
   /* read the number of duration pdfs */
   HTS_Fread(&ms->ndurpdf, sizeof(int), 1, ms->fp[DUR]);
   if (ms->ndurpdf<0)
      HTS_Error(1, "LoadModelFiles: #duration pdf must be positive value.\n");

   ms->durpdf = (double **) HTS_Calloc(ms->ndurpdf, sizeof(double *));
   ms->durpdf--;
   
   /* read pdfs (mean & variance) */
   for (i=1; i<=ms->ndurpdf; i++) {
      ms->durpdf[i] = (double *) HTS_Calloc(2*ms->nstate, sizeof(double));
      ms->durpdf[i]-= 2;
      for (j=2;j<2*(ms->nstate+1);j++) {
         HTS_Fread(&temp, sizeof(float), 1, ms->fp[DUR]);
         ms->durpdf[i][j] = (double) temp;
      }
   }

   /*-------------------- load pdfs for mcep --------------------*/
   /* read vector size for spectrum */
   HTS_Fread(&ms->mcepvsize, sizeof(int), 1, ms->fp[MCP]);
   if (ms->mcepvsize<0)
      HTS_Error(1, "LoadModelFiles: vector size of mel-cepstrum part must be positive.\n");

   ms->nmceppdf = (int *) HTS_Calloc(ms->nstate, sizeof(int));
   ms->nmceppdf-= 2;

   /* read the number of pdfs for each state position */
   HTS_Fread(&ms->nmceppdf[2], sizeof(int), ms->nstate, ms->fp[MCP]);   
   for (i=2;i<=ms->nstate+1;i++) { 
      if (ms->nmceppdf[i]<0)
         HTS_Error(1, "LoadModelFiles: #mcep pdf at state %d must be positive value.\n", i);
   }
   ms->mceppdf = (double ***) HTS_Calloc(ms->nstate, sizeof(double **));
   ms->mceppdf-= 2;
   
   /* read pdfs (mean, variance) */
   for (i=2; i<=ms->nstate+1; i++) {
      ms->mceppdf[i] = (double **) HTS_Calloc(ms->nmceppdf[i], sizeof(double *));
      ms->mceppdf[i]--; 
      for (j=1; j<=ms->nmceppdf[i]; j++) {
         ms->mceppdf[i][j] = (double *) HTS_Calloc(2*ms->mcepvsize, sizeof(double));
         for (k=0;k<2*ms->mcepvsize;k++) {
            HTS_Fread(&temp, sizeof(float), 1, ms->fp[MCP]);
            ms->mceppdf[i][j][k] = (double)temp;
         }
      }
   } 

   /*-------------------- load pdfs for log F0 --------------------*/
   /* read the number of streams for f0 modeling */
   HTS_Fread(&ms->lf0stream, sizeof(int), 1, ms->fp[LF0]);
   if (ms->lf0stream<0)
      HTS_Error(1, "LoadModelFiles: #stream for log f0 part must be positive value.\n");
      
   ms->nlf0pdf = (int *) HTS_Calloc(ms->nstate, sizeof(int));
   ms->nlf0pdf-= 2;

   /* read the number of pdfs for each state position */
   HTS_Fread(&ms->nlf0pdf[2], sizeof(int), ms->nstate, ms->fp[LF0]);
   for (i=2;i<=ms->nstate+1;i++) { 
      if (ms->nlf0pdf[i]<0)
         HTS_Error(1, "LoadModelFiles: #lf0 pdf at state %d must be positive.\n", i);
   }   
   ms->lf0pdf = (double ****) HTS_Calloc(ms->nstate, sizeof(double ***));
   ms->lf0pdf-= 2;
   
   /* read pdfs (mean, variance & weight) */
   for (i=2; i<=ms->nstate+1; i++) {
      ms->lf0pdf[i] = (double ***) HTS_Calloc(ms->nlf0pdf[i], sizeof(double **));
      ms->lf0pdf[i]--; 
      for (j=1; j<=ms->nlf0pdf[i]; j++) {
         ms->lf0pdf[i][j] = (double **) HTS_Calloc(ms->lf0stream, sizeof(double *));
         ms->lf0pdf[i][j]--; 
         for (k=1; k<=ms->lf0stream; k++) {
            /* 4 -> mean, variance, voiced weight, and unvoiced weight */
            ms->lf0pdf[i][j][k] = (double *) HTS_Calloc(4, sizeof(double));
            for (l=0;l<4;l++) {
               HTS_Fread(&temp, sizeof(float), 1, ms->fp[LF0]);
               ms->lf0pdf[i][j][k][l] = (double)temp;
            }
            
            vw  = ms->lf0pdf[i][j][k][2]; /* voiced weight */
            uvw = ms->lf0pdf[i][j][k][3]; /* unvoiced weight */
            if (vw<0.0 || uvw<0.0 || vw+uvw<0.99 || vw+uvw>1.01)
               HTS_Error(1, "LoadModelFiles: voiced/unvoiced weights must be within 0.99 to 1.01.\n");
         }
      }
   }
   
   return;
}

/* ClearModelSet: clear models */
void ClearModelSet (ModelSet *ms) 
{
   int i, j, k;

   /* close */
   fclose(ms->fp[DUR]);
   fclose(ms->fp[LF0]);
   fclose(ms->fp[MCP]);
   
   /* free models for f0 */
   for (i=ms->nstate+1; i>=2; i--) {
      for (j=ms->nlf0pdf[i]; j>=1; j--) {
         for (k=ms->lf0stream; k>=1; k--) {
            HTS_Free(ms->lf0pdf[i][j][k]);
         }
         HTS_Free(++ms->lf0pdf[i][j]);
      }
      HTS_Free(++ms->lf0pdf[i]);
      fflush(stdout);
   }
   ms->lf0pdf += 2;  ms->nlf0pdf += 2;
   HTS_Free(ms->lf0pdf);
   HTS_Free(ms->nlf0pdf);

   /* free models for spectrum */
   for (i=ms->nstate+1; i>=2; i--) {
      for (j=ms->nmceppdf[i]; j>=1; j--) {
         HTS_Free(ms->mceppdf[i][j]);
      }
      HTS_Free(++ms->mceppdf[i]);
   }
   ms->mceppdf += 2;  ms->nmceppdf += 2;
   HTS_Free(ms->mceppdf);
   HTS_Free(ms->nmceppdf);

   /* free models for duration */
   for (i=ms->ndurpdf; i>=1; i--) {
      ms->durpdf[i] += 2;
      HTS_Free(ms->durpdf[i]);
   }
   HTS_Free(++ms->durpdf); 
   
   return;
}

/* FindDurPDF: find duration pdf from pdf array */
void FindDurPDF (Model *m, ModelSet *ms, const double rho, double *diffdur)
{
   double data, mean, variance;
   int s; 

   const int idx    = m->durpdf;
   const int nstate = ms->nstate;

   for (s=2; s<=nstate+1; s++) {
      mean = ms->durpdf[idx][s];
      variance = ms->durpdf[idx][nstate+s];
      data = mean + rho*variance;
      
      m->dur[s] = (int) (data+*diffdur+0.5);
      m->dur[s] = (m->dur[s]<1) ? 1 : m->dur[s]; 
      
      m->totaldur += m->dur[s];
      *diffdur += data-(double)m->dur[s];
   }
   
   return;
}

/* FindLF0PDF : find required pdf for log F0 from pdf array */
void FindLF0PDF (const int s, Model *m, ModelSet *ms, const double uvthresh)
{
   int stream;
   double *weight;

   const int idx     = m->lf0pdf[s];
   const int nstream = ms->lf0stream;

   for (stream=1; stream<=nstream; stream++) {
      m->lf0mean    [s][stream] = ms->lf0pdf[s][idx][stream][0];
      m->lf0variance[s][stream] = ms->lf0pdf[s][idx][stream][1];
      weight = ms->lf0pdf[s][idx][stream]+2;
      
      if (stream==1) {
         if (weight[0]>uvthresh)
            m->voiced[s] = 1;
         else
            m->voiced[s] = 0;
      }
   }
   
   return;
}

/* FindMcpPDF : find pdf for mel-cepstrum from pdf array */
void FindMcpPDF (const int s, Model *m, ModelSet *ms)
{
   const int idx = m->mceppdf[s];

   m->mcepmean[s] = ms->mceppdf[s][idx];
   m->mcepvariance[s] = ms->mceppdf[s][idx]+ms->mcepvsize;
   
   return;
}

/* InitModelSet: initialise model set */
void InitModelSet (ModelSet *ms)
{
   ms->fp[DUR] = NULL;
   ms->fp[LF0] = NULL;
   ms->fp[MCP] = NULL;
   
   return;
} 

/* -------------------- End of "model.cc" -------------------- */
