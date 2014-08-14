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
/*    misc.c: miscellaneous functions                                */
/*  ---------------------------------------------------------------  */

/* hts_engine libraries */
#include "misc.hpp"

#ifdef FESTIVAL
#include "EST_walloc.h"
#endif

/* ----- Routines for file input/output ----- */
/* HTS_Usage: output usage */
void HTS_Usage (void)
{
   fprintf (stderr, "\n");
   fprintf (stderr, "hts_engine - A HMM-based speech synthesis engine\n");
   fprintf (stderr, "\n");
   fprintf (stderr, "  usage:\n");
   fprintf (stderr, "       hts_engine [ options ] [ infile ] \n");
   fprintf (stderr, "  options:                                                          [def] [min--max]\n");
   fprintf (stderr, "       -td tree  : decision trees file for state duration           [N/A]\n");
   fprintf (stderr, "       -tf tree  : decision trees file for Log F0                   [N/A]\n");
   fprintf (stderr, "       -tm tree  : decision trees file for Mel-Cepstrum             [N/A]\n");
   fprintf (stderr, "       -md pdf   : model file for state duration                    [N/A]\n");
   fprintf (stderr, "       -mf pdf   : model file for Log F0                            [N/A]\n");
   fprintf (stderr, "       -mm pdf   : model file for Mel-Cepstrum                      [N/A]\n");
   fprintf (stderr, "       -df win   : window file for culcuration delta of log F0      [N/A]\n");
   fprintf (stderr, "       -dm win   : filename of delta coeffs for Mel-Cepstrum        [N/A]\n");
   fprintf (stderr, "       -od s     : filename of output duration                      [N/A]\n");
   fprintf (stderr, "       -of s     : filename of output F0                            [N/A]\n");
   fprintf (stderr, "       -om s     : filename of output mcep                          [N/A]\n");
   fprintf (stderr, "       -or s     : filename of output raw audio (generated speech)  [N/A]\n");
   fprintf (stderr, "       -ot s     : filename of output trace information             [N/A]\n");
   fprintf (stderr, "       -vs       : use state alignment for duration                 [0]\n");
   fprintf (stderr, "       -vp       : use phoneme alignment for duration               [0]\n");
   fprintf (stderr, "       -s  i     : sampring frequency                               [16000][  0--48000]\n");
   fprintf (stderr, "       -p  i     : frame period (point)                             [ 80 ][   0--2000]\n");
   fprintf (stderr, "       -a  f     : all-pass constant                                [0.42][ 0.0--1.0]\n");
   fprintf (stderr, "       -b  f     : postfiltering coefficient                        [0.00][-0.8--0.8]\n");
   fprintf (stderr, "       -r  f     : control duration parameter                       [0.00][-1.0--0.0]\n");
   fprintf (stderr, "       -fs f     : multilply f0                                     [1.00][ 0.0--5.0]\n");
   fprintf (stderr, "       -fm f     : add f0                                           [0.00][ 0.0--100.0]\n");
   fprintf (stderr, "       -u  f     : voiced/unvoiced threshold                        [0.50][ 0.0--1.0]\n");
   fprintf (stderr, "       -l  f     : length of generated speech (in second)           [N/A] [ 0.0--30.0]\n");
   fprintf (stderr, "  infile:\n");
   fprintf (stderr, "       label file\n");
   fprintf (stderr, "  note:\n");
   fprintf (stderr, "       option '-d' may be repeated to use multiple\n");
   fprintf (stderr, "       delta parameters\n");
   fprintf (stderr, "       if 'd' is an integer, delta coefficients are\n");
   fprintf (stderr, "       calculated based on a regression formula.\n\n");
   fprintf (stderr, "       generated mel-cepstrum and log F0 sequences are \n");
   fprintf (stderr, "       saved in natural endian, binary (float) format.\n");
   fprintf (stderr, "\n");
 
   exit(0);
   
   return;
}

/* HTS_Error: output error message */
void HTS_Error (const int error, char *message, ...)
{
   va_list arg;
   
   fflush(stdout);
   fflush(stderr);
   
   if (error>0)
      fprintf(stderr, "\nError: ");
   else 
      fprintf(stderr, "\nWarning: ");
         
   va_start(arg, message);
   vfprintf(stderr, message, arg);
   va_end(arg);
   
   fflush(stderr);
   
   if (error>0)
      exit(error);
   
   return;
}

/* HTS_Getfp: wrapper for fopen */
FILE *HTS_Getfp (const char *name, const char *opt)
{
   FILE *fp = fopen(name, opt);
   
   if (fp==NULL)
      HTS_Error(2, "HTS_Getfp: Cannot open %s.\n", name);
   
   return (fp);
}

/* HTS_GetToken: parser */ 
void HTS_GetToken (FILE *fp, char *buff)
{
   char c;
   int i;
   HTS_Boolean squote = 0, dquote = 0;

   c = fgetc(fp);

   while (isspace(c))
      c = fgetc(fp);
      
   if (c=='\'') {  /* single quote case */
      c = fgetc(fp);
      squote = 1;
   }
   
   if (c=='\"') {  /*double quote case */
      c = fgetc(fp);
      dquote = 1;
   }
   
   if (c==',') {   /*special character ',' */
      strcpy(buff, ",");
      return; 
   }
   
   i = 0;
   while (1) {
      buff[i++] = c;
      c = fgetc(fp);
      if (squote && c == '\'') break;
      if (dquote && c == '\"') break;
      if (!(squote || dquote || isgraph(c)) ) break;
   }
   
   buff[i]=0;
   
   return;
}

/* ----- Routines for memory allocation/free ----- */
/* HTS_Calloc: wrapper for calloc */
char *HTS_Calloc (const size_t num, const size_t size)
{
#ifdef FESTIVAL
   char *mem = (char *)safe_wcalloc(num * size);
#else
   char *mem = (char *)calloc(num, size);
#endif
   
   if (mem==NULL)
      HTS_Error(1, "HTS_calloc: Cannot allocate memory.\n");
   
   return(mem);
}

/* HTS_Free: wrapper for free */
void HTS_Free (void *ptr)
{
#ifdef FESTIVAL
   wfree(ptr);
#else
   free(ptr);
#endif
   
   return;
}

/* HTS_Strdup: wrapper for strdup */
char *HTS_Strdup (const char *in)
{
#ifdef FESTIVAL
   return (wstrdup(in));
#else
   char *tmp = (char *) HTS_Calloc(strlen(in)+1, sizeof(char));
   strcpy(tmp, in);
   return tmp;
#endif
}

/* HTS_AllocVector: allocate vector */
double *HTS_AllocVector (const int x)
{
   double *ptr = (double *) HTS_Calloc(x, sizeof(double));
   
   ptr--;
   
   return(ptr);
}

/* HTS_AllocMatrix: allocate matrix */
double **HTS_AllocMatrix (const int x, const int y)
{
   int i;
   double **ptr = (double **) HTS_Calloc(x, sizeof(double *));
 
   ptr--;
   
   for (i=1; i<=x; i++)
      ptr[i] = HTS_AllocVector(y);
   
   return(ptr);
}

/* HTS_FreeVector: free vector */
void HTS_FreeVector (double *ptr) 
{
   ptr++;
   
   HTS_Free((void *)ptr);

   return;
}

/* HTS_FreeMatrix: free matrix */
void HTS_FreeMatrix (double **ptr, const int x)
{
   int i;
   
   for (i=x;i>0;i--)
      HTS_FreeVector(ptr[i]);

   ptr++;
   
   HTS_Free((void *)ptr);
   
   return;
}

/* ----- Routines for reading from binary file ----- */ 
/* HTS_ByteSwap: byte swap */
int HTS_ByteSwap (void *p, const int size, const int blocks)
{
   char *q, tmp;
   int i, j;

   q = (char *)p;

   for (i=0; i<blocks; i++) {
      for (j=0; j<(size/2); j++) {
         tmp = *(q+j);
         *(q+j) = *(q+(size-1-j));
         *(q+(size-1-j)) = tmp;
      }
      q += size;
   }
   
   return i;
}

/* HTS_Fread: fread with byteswap */
int HTS_Fread (void *p, const int size, const int num, FILE *fp)
{
   const int block = fread(p, size, num, fp);

#ifndef WORDS_BIGENDIAN
   HTS_ByteSwap(p, size, block);
#endif /* !BIG_ENDIAN */

   return block;
}

/* -------------------- End of "misc.cc" -------------------- */
