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
/*    vocoder.h:  mel-cepstral vocoder                               */
/*               (pulse/noise excitation & MLSA filter)              */
/*  ---------------------------------------------------------------  */

#define RANDMAX 32767 

#define   IPERIOD    1
#define   SEED       1
#define   B0         0x00000001
#define   B28        0x10000000
#define   B31        0x80000000
#define   B31_       0x7fffffff
#define   Z          0x00000000

#ifdef HTS_EMBEDDED
   #define   GAUSS      0
   #define   PADEORDER  4  /* pade order (for MLSA filter) */
   #define   IRLENG     64 /* length of impulse response */
#else
   #define   GAUSS      1
   #define   PADEORDER  5
   #define   IRLENG     96
#endif

/* VocoderSetup: structure for setting of vocoder */
typedef struct _VocoderSetup {
   int fprd;
   int iprd;
   int seed;
   int pd;
   unsigned long next;
   HTS_Boolean gauss;
   double p1;
   double pc;
   double pj;
   double pade[21];
   double *ppade;
   double *c, *cc, *cinc, *d1;
   double rate;
   
   int sw;
   double r1, r2, s;
   
   int x;
   
   /* for postfiltering */
   int size;
   double *d; 
   double *g;
   double *mc;
   double *cep;
   double *ir;
   int o;
   int irleng;  
} VocoderSetup;

/* Function prototypes for vocoding */
void InitVocoder(FILE *, int , VocoderSetup *, const int, const int);
void ClearVocoder(FILE *, VocoderSetup *);
void Vocoder(double , double *, int, FILE *, globalP *, VocoderSetup *);

/* -------------------- End of "vocoder.h" -------------------- */
