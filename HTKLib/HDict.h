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
/*         File: HDict.h  Dictionary Storage                   */
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
/*                     Copyright (c) 2001-2006                       */
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

/* !HVER!HDict:   3.3 [CUED 28/04/05] */

#ifndef _HDICT_H_
#define _HDICT_H_

#ifdef __cplusplus
extern "C" {
#endif

/* size of hash table */
#define VHASHSIZE 701

/* max number of phones in a pronunciation */
#define MAXPHONES 256

/* lowest prob not floored to LZERO */
#define MINPRONPROB 1E-6

typedef struct _DictEntry *Word;
typedef struct _WordPron  *Pron;

typedef struct _WordPron{   /* storage for each pronunciation */
   short pnum;     /* Pronunciation number 1..nprons */
   short nphones;  /* Number of phones in pronuciation */
   LabId *phones;  /* Array[0..nphones-1] of phones */
   LogFloat prob;  /* Log probability of pronunciation */
   LabId outSym;   /* Output symbol generated when pronunciation recognised */
   Word word;      /* Word this is a pronuciation of */
   Pron next;      /* Next pronunciation of word */
   void *aux;      /* hook for temp info */
} WordPron;

typedef struct _DictEntry{
   LabId wordName;  /* word identifier */
   Pron pron;       /* first pronunciation */
   int nprons;      /* number of prons for this word */
   Word next;       /* next word in hash table chain */
   void *aux;       /* hook used by HTK library modules for temp info */
} DictEntry;

typedef struct {
   int nwords;          /* total number of words */
   int nprons;          /* total number of prons */
   Word nullWord;       /* dummy null word/node */
   Word subLatWord;     /* special word for HNet subLats */
   Word *wtab;          /* hash table for DictEntry's */
   MemHeap heap;        /* storage for dictionary */
   MemHeap wordHeap;    /* for DictEntry structs  */
   MemHeap pronHeap;    /* for WordPron structs   */
   MemHeap phonesHeap;  /* for arrays of phones   */
} Vocab;


void InitDict(void);
/* 
   Initialise the HDict module 
*/

void ResetDict(void);
/* 
   Rest module 
*/

void InitVocab(Vocab *voc);
/* 
   Initialise a Vocab data structure - must be called before
   and other ReadDict/GetWord routines is used
*/

void ClearVocab(Vocab *voc);
/* 
   Removes all the datstrcutures created by InitVocab
*/

ReturnStatus ReadDictWord(Source *src,LabId *labels,float *prob,int* num);
/*
   Read a single word and definitiion from source.
   Format is one pronunciation per line of form
        WORD [OUTSYM] a b c d ....
   meaning that WORD has pronunciation a b c d ... and OUTSYM is
   output when WORD is recognised.  If OUTSYM is omitted then
   OUTSYM = WORD.  (OUTSYM is used by HDMan to track source of
   pronunciations).

   With appropriate #def a modified dictionary is allowed
        WORD [OUTSYM] # a b c d ....
   with each pronunciaiton havine a probability # (assumed to be 1.0 if no
   number is present).

   *num equals no. phones read (-1 indicates EOF)
   returns ReturnStatus - FAIL if format of entry problem.
   labels[0]=WORD, labels[1]=OUTSYM, labels[2]=a ..
*/

ReturnStatus ReadDict(char *dictFn, Vocab *voc);
/* 
   Read a dictionary from dictFn and store in voc.  
   Uses ReadDictWord to read file.
*/

ReturnStatus WriteDict(char *dictFn, Vocab *voc);
/* 
   Write the given Vocab structure to the file dictFn
*/


void ShowDict(Vocab *voc);
/* 
   Print out the pronunciation dictionary voc
*/


Word GetWord(Vocab *voc, LabId wordName, Boolean insert);
/* 
   Return the Word with name wordName from Vocab voc.  If
   insert and wordName not in voc, then a new entry is
   created with a null pronunciation.
*/

void DelWord(Vocab *voc, Word word);
/* 
   Delete word and its pronunciations from voc 
*/

void NewPron(Vocab *voc, Word wid,int nphones,LabId *phones,
	     LabId outSym,float prob);
/*
   Add a pronunciation and output symbol to a given word -
   pronunciation is stored in phones, prob gives the probability
   of the pronunciation (pron->prob = log(prob)).
   Note that none of the LabId's in phones or outSym can be NULL
*/

void DelPron(Vocab *voc, Word word, Pron pron);
/* 
   DelPron: delete a specific word pronunciation 
*/

#ifdef __cplusplus
}
#endif

#endif  /* _HDICT_H_ */

/* ------------------------ End of HDict.h ----------------------- */

