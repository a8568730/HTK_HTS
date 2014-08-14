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
/*    tree.h: decision tree definition                               */
/*  ---------------------------------------------------------------  */

/* Pattern: structure for pattern */
typedef struct _Pattern{
   char *pat;               /* pattern */
   struct _Pattern *next;   /* link to next pattern */
} Pattern;

/* Question: structure for questions */
typedef struct _Question {
   char *qName;             /* name of this question */
   Pattern *phead;          /* link to head of pattern list */
   Pattern *ptail;          /* link to tail of pattern list */
   struct _Question *next;  /* link to next question */
} Question;

/* Node: structure for node of decision tree */
typedef struct _Node {
   int idx;                 /* index of this node */
   int pdf;                 /* index of pdf for this node  ( leaf node only ) */
   struct _Node *yes;       /* link to child node (yes) */
   struct _Node *no;        /* link to child node (no)  */
   struct _Node *next;        /* link to next node  */  
   Question *quest;         /* question applied at this node */
} Node;
   
/* Tree: structure for each decision tree */
typedef struct _Tree {      
   int state;               /* state position of this tree */
   Pattern *phead;          /* link to head of pattern list for this tree */
   Pattern *ptail;          /* link to tail of pattern list for this tree */
   struct _Tree *next;      /* link to next tree */
   Node *root;              /* root node of this decision tree */
   Node *leaf;              /* leaf nodes of this decision tree */
} Tree;

/* TreeSet: structure for decision tree set */
typedef struct _TreeSet {
   Question *qhead[HTS_NUMMTYPE];      /* question lists for mcep, logF0, and duration */
   Question *qtail[HTS_NUMMTYPE];     
   Tree *thead[HTS_NUMMTYPE];          /* tree lists for mcep, logF0, and duration */
   Tree *ttail[HTS_NUMMTYPE];
   int nTrees[HTS_NUMMTYPE];           /* # of trees for mcep, logF0, and duration */
   FILE *fp[HTS_NUMMTYPE];             /* file pointers for mcep, logF0, and duration */
} TreeSet;

/* Function prototypes for tree handling */
void LoadTreeSet  (TreeSet *, const HTS_Mtype);
void ClearTreeSet (TreeSet *, const HTS_Mtype);

int SearchTree (char *, Node *);

void InitTreeSet(TreeSet *);

/* -------------------- End of "tree.h" -------------------- */

