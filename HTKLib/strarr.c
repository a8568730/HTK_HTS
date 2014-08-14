/* ----------------------------------------------------------- */
/*                                                             */
/*                          ___                                */
/*                       |_| | |_/   SPEECH                    */
/*                       | | | | \   RECOGNITION               */
/*                       =========   SOFTWARE                  */
/*                                                             */
/*                                                             */
/* ----------------------------------------------------------- */
/*         Copyright: Microsoft Corporation                    */
/*          1995-2000 Redmond, Washington USA                  */
/*                    http://www.microsoft.com                 */
/*                                                             */
/*   Use of this software is governed by a License Agreement   */
/*    ** See the file License for the Conditions of Use  **    */
/*    **     This banner notice must not be removed      **    */
/*                                                             */
/* ----------------------------------------------------------- */

/*  *** THIS IS A MODIFIED VERSION OF HTK ***                        */
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
/*  distribute this software in the form of patch code to HTK and    */
/*  its documentation without restriction, including without         */
/*  limitation the rights to use, copy, modify, merge, publish,      */
/*  distribute, sublicense, and/or sell copies of this work, and to  */
/*  permit persons to whom this work is furnished to do so, subject  */
/*  to the following conditions:                                     */
/*                                                                   */
/*    1. Once you apply the HTS patch to HTK, you must obey the      */
/*       license of HTK.                                             */
/*                                                                   */
/*    2. The source code must retain the above copyright notice,     */
/*       this list of conditions and the following disclaimer.       */
/*                                                                   */
/*    3. Any modifications to the source code must be clearly        */
/*       marked as such.                                             */
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

#include <malloc.h>
#include "esignal.h"
#include "strarr.h"
/*
 * "dim" should point to the dimensions of a char matrix.
 * "data" should point to the matrix data, stored in row-major order.
 * Each row is assumed to be terminated by one or more null characters.
 * This function converts the matrix to a NULL-terminated array of
 * strings---one per row of the matrix---and returns a pointer to the
 * result.
 */

char **
StrArrFromRect(long *dim, void *data)
{
   char **strarr;
   long len, wid, i;
   char *row;

   len = dim[0];

   strarr = (char **)malloc((len + 1) * sizeof(char *));

   wid = dim[1];
   row = (char *) data;
   for (i = 0; i < len; i++)
      {
         /* Assumption: each row has at least 1 terminating null */
         strarr[i] = StrDup(row);
         row += wid;
      }

   strarr[len] = NULL;

   return strarr;
}


/*
 * "strarr" should point to the beginning of a NULL-terminated string array.
 * "dimenp" and "datap" are output variables.  If not NULL, they should
 * be the addresses of variables to which results will be assigned.
 * The function converts the data in the string array to a character matrix,
 * stored in row-major order.  Each row receives the contents of one string
 * padded with null characters to bring them all up to a common length;
 * every row gets at least one terminal null.  A pointer to the dimensions
 * is returned via "dimenp".  A pointer to the matrix data is returned via
 * "datap".
 */

void
StrArrToRect(char **strarr, long **dimenp, void **datap)
{
   long len, wid, *dim, siz;
   void *data;
   long i, j;
   char *str, *row;

   len = StrArrLen(strarr);
   wid = StrArrMaxLen(strarr) + 1;

   dim = (long *) malloc(2 * sizeof(long));
   dim[0] = len;
   dim[1] = wid;

   siz = len * wid;
   data = malloc(siz * sizeof(char));
   row = (char *) data;
   for (i = 0; i < len; i++)
      {
         str = strarr[i];
         for (j = 0; str[j] != '\0'; j++)
            row[j] = str[j];
         for ( ; j < wid; j++)
            row[j] = '\0';
         row += wid;
      }

   if (dimenp)
      *dimenp = dim;
   if (datap)
      *datap = data;
}


/*
 * Number of strings in a NULL-terminated string array.
 */

int
StrArrLen(char **str_arr)
{
   int i = 0;

   if ((str_arr == NULL) || (*str_arr == NULL))
      return 0;

   while (str_arr[i] != NULL) 
      i++;

   return(i);
}


/*
 * Maximum strlen for the members of a NULL-terminated string array.
 */

int
StrArrMaxLen(char **str_arr)
{
   int maxlen = 0;
   int newmax;

   if (str_arr == NULL)
      return(0);
  
   while (*str_arr != NULL) {
      if ((newmax = strlen(*str_arr)) > maxlen)
         maxlen = newmax;
      str_arr++;
   }

   return(maxlen);
}
