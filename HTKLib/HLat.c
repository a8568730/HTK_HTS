/* ----------------------------------------------------------- */
/*                                                             */
/*                          ___                                */
/*                       |_| | |_/   SPEECH                    */
/*                       | | | | \   RECOGNITION               */
/*                       =========   SOFTWARE                  */ 
/*                                                             */
/*                                                             */
/* ----------------------------------------------------------- */
/* developed at:                                               */
/*                                                             */
/*      Speech Vision and Robotics group                       */
/*      Cambridge University Engineering Department            */
/*      http://svr-www.eng.cam.ac.uk/                          */
/*                                                             */
/* author: Gunnar Evermann <ge204@eng.cam.ac.uk>               */
/* ----------------------------------------------------------- */
/*         Copyright:                                          */
/*         2001-2002  Cambridge University                     */
/*                    Engineering Department                   */
/*                                                             */
/*   Use of this software is governed by a License Agreement   */
/*    ** See the file License for the Conditions of Use  **    */
/*    **     This banner notice must not be removed      **    */
/*                                                             */
/* ----------------------------------------------------------- */
/*       File: HLat.c:  Lattice Manipulation                   */
/* ----------------------------------------------------------- */

/*#### todo:

     - implement lattice oracle WER calculation
     - allow batch processing?
*/


char *hlat_version = "!HVER!HLat:   3.2.1 [CUED 15/10/03]";
char *hlat_vc_id = "$Id: HLat.c,v 1.4 2003/10/15 08:10:12 ge204 Exp $";


#include "HShell.h"
#include "HMem.h"
#include "HMath.h"
#include "HWave.h"
#include "HAudio.h"
#include "HParm.h"
#include "HLabel.h"
#include "HModel.h"
#include "HUtil.h"
#include "HDict.h"
#include "HNet.h"
#include "HLM.h"
#include "HLat.h"

/* ----------------------------- Trace Flags ------------------------- */

#define T_TOP  00001
#define T_PRUN 00002
#define T_FB   00004
#define T_EXP  00010
#define T_MEM  00020
#define T_TRAN 00040

static int trace=0;
static ConfParam *cParm[MAXGLOBS];      /* config parameters */
static int nParm = 0;

/* --------------------------- Global Flags -------------------------- */

static LabId startWord;         /* word at start of Lattice (!SENT_START) */
static LabId endWord;           /* word at end of Lattice (!SENT_END) */
static LabId startLMWord;       /* word at start in LM (<s>) */
static LabId endLMWord;         /* word at end in LM (</s>) */
static LabId nullWord;          /* null word in Lattices (!NULL) */

static MemHeap slaHeap, slnHeap;/* MHEAPs for use in LatExpand() */

/* --------------------------- Prototypes ---------------------------- */


#ifndef NO_LAT_LM
typedef struct _SubLNode SubLNode;
typedef struct _SubLArc SubLArc;

struct _SubLNode {
   union {
      LMState lmstate;
      LNode *newln;
   } data;
   SubLArc *foll;
   SubLNode *next;
};

struct _SubLArc {
   LogFloat lmprob;
   SubLNode *end;
   LArc *la;
   SubLArc *next;
};
#endif

/* --------------------------- Initialisation ---------------------- */

/* EXPORT->InitLat: register module & set configuration parameters */
void InitLat(void)
{
   int i;

   Register(hlat_version,hlat_vc_id);
   nParm = GetConfig("HLAT", TRUE, cParm, MAXGLOBS);
   if (nParm>0){
      if (GetConfInt(cParm,nParm,"TRACE",&i)) trace = i;
   }

#ifndef NO_LAT_LM
   CreateHeap (&slaHeap, "LatExpand arc heap", MHEAP, sizeof (SubLArc), 1.0, 1000, 128000);
   CreateHeap (&slnHeap, "LatExpand node heap", MHEAP,sizeof (SubLNode), 1.0, 1000, 32000);
#endif
}


/* --------------------------- Lattice processing ------------------- */

/* LatCheck

     check lattice for consistency: single start & end nodes; no cycles.
*/
void LatCheck (Lattice *lat)
{
   int i, nStart, nEnd;
   LNode *ln;
   LNode **topOrder;

   /* check for unique start and end nodes */
   nStart = nEnd = 0;
   for (i = 0, ln = lat->lnodes; i < lat->nn; ++i, ++ln) {
      if (!ln->pred)
         ++nStart;
      if (!ln->foll)
         ++nEnd;
   }

   if (nStart != 1)
      HError (8622, "HLat: lattice has %d start nodes (should be 1)", nStart);
   if (nEnd != 1)
      HError (8622, "HLat: lattice has %d end nodes (should be 1)", nEnd);

   /* check wheter lat is a DAG ( <=> top order exists). */
   topOrder = (LNode **) New (&gcheap, lat->nn * sizeof(LNode *));
   if (!LatTopSort (lat, topOrder))
      HError (8622, "HLat: lattice contains cylces");

   Dispose (&gcheap, topOrder);
}


/* FixPronProbs

     replace pronunciation probabilities in lattices with values taken from
     the dictionary
*/
void FixPronProbs (Lattice *lat, Vocab *voc)
{
   int i, v;
   LNode *ln;
   Pron pron;
   LArc *la;

   for (i = 0, ln = lat->lnodes; i < lat->nn; ++i, ++ln) {
      if (ln->word != voc->nullWord) {
         if (ln->v > ln->word->nprons)
            HError (8621, "FixPronprbs: lattice refers to non-existing pron");
         pron = ln->word->pron;
         for (v = 2; v <= ln->v; ++v)
            pron = pron->next;
         if (pron->pnum != ln->v)
            HError (8621, "FixPronprbs: dict pron numbering invalid");
         
         for (la = ln->pred; la; la = la->parc)
            la->prlike = pron->prob;
      }
      else {
         for (la = ln->pred; la; la = la->parc)
            la->prlike = 0.0;
      }
   }
}


/* LatStartNode

     return lattice start node

     ###GE: maybe store this in Lattice structure to save time?
*/
LNode *LatStartNode (Lattice *lat)
{
   int i;
   LNode *ln;

   for (i = 0, ln = lat->lnodes; i < lat->nn; ++i, ++ln)
      if (!ln->pred)
         return ln;

   HError (8622, "HLat: lattice has no start node");
   return NULL;         /* make compiler happy */
}

/* LatEndNode

     return lattice end node

     ###GE: maybe store this in Lattice structure to save time?
*/
LNode *LatEndNode (Lattice *lat)
{
   int i;
   LNode *ln;

   for (i = 0, ln = lat->lnodes; i < lat->nn; ++i, ++ln)
      if (!ln->foll)
         return ln;

   HError (8622, "HLat: lattice has no end node");
   return NULL;         /* make compiler happy */
}


/* helper function for LatTopSort()
   number preceeding nodes depth first 
*/
void LatTopSortVisit (LNode *ln, int *time)
{
   LArc *la;

   ln->n = -2;          /* mark node as seen, but not numbered, yet (GRAY in CLR) */
   for (la = ln->pred; la; la = la->parc)
      if (la->start->n == -1)
         LatTopSortVisit (la->start, time);
   
   ++(*time);
   ln->n = *time;
}

/* LatTopSort

     sort lattice nodes in topological order
     returns array of LNode pointers.
     uses depth first traversal (in reverse (pred) direction) 
     based on [CLR:1990], p. 485 
*/
Boolean LatTopSort (Lattice *lat, LNode **topOrder)
{
   int time = -1;       /* new numbering will start at 0 */
   int i;
   LNode *ln;
   LArc *la;
   Boolean isDAG = TRUE;

   for (i = 0, ln = lat->lnodes; i < lat->nn; ++i, ++ln)
      ln->n = -1;        /* mark node as unseen (WHITE in CLR) */

   LatTopSortVisit (LatEndNode (lat), &time);
   
   /* we should have seen all nodes */
   assert (time+1 == lat->nn); 

   /* check topological order */
   for (i = 0, la = lat->larcs; i < lat->na; ++i, ++la)
      if (!(la->start->n < la->end->n)) {
         isDAG = FALSE;
         HError (-8622, "LatTopSort: Lattice contains cycles"); 
         break;
      }

   for (i = 0, ln = lat->lnodes; i < lat->nn; ++i, ++ln)
      topOrder[ln->n] = ln;
   
   /*#### GE: reset ln->n values? */
   for (i = 0, ln = lat->lnodes; i < lat->nn; ++i, ++ln)
      ln->n = i;

   return isDAG;
}

/* LatAttachInfo

     allocate & attach an Info structre for each node
*/
void LatAttachInfo (MemHeap *heap, size_t size, Lattice *lat)
{
   int i;
   LNode *ln;
   
   for (i = 0, ln = lat->lnodes; i < lat->nn; ++i, ++ln)
      ln->hook = (Ptr) New (heap, size);
}

/* LatDetachInfo

     free Info structre for each node
*/
void LatDetachInfo (MemHeap *heap, Lattice *lat)
{
   int i;
   LNode *ln;
   
   for (i = 0, ln = lat->lnodes; i < lat->nn; ++i, ++ln)
      Dispose (heap, ln->hook);
}


/* LatForwBackw

     perform forward-backward algorithm on lattice and store scores in
     FBInfo structre
     choice of using sum (LATFB_SUM) or max (LATFB_MAX) of scores
*/
LogDouble LatForwBackw (Lattice *lat, LatFBType type)
{
   int i;
   LNode *ln;
   LArc *la;
   LNode **topOrder;
   LogDouble score;

   /* We assume that the FBinfo structures are already allocated. */
   /* init */
   for (i = 0, ln = lat->lnodes; i < lat->nn; ++i, ++ln) {
      LNodeFw (ln) = LNodeBw (ln) = LZERO;
   }
   LNodeFw (LatStartNode (lat)) = 0.0;
   LNodeBw (LatEndNode (lat)) = 0.0;

   /* find topological order of nodes */
   topOrder = (LNode **) New (&gcheap, lat->nn * sizeof(LNode *));
   if (!LatTopSort (lat, topOrder))
      HError (8622, "LatForwBackw: cannot calculate forw/backw score on Lattice with cycles"); 
      
   /* the processing in the forward and backward directions are really
      almost the same. if the data structures had been done _slightly_
      nicer this could be done in one loop. The only readable way now 
      would be defining a couple of macros... */

   /* forward direction */
   for (i = 0; i < lat->nn; ++i) {
      ln = topOrder[i];
      for (la = ln->foll; la; la = la->farc) {
         assert (la->start == ln);
         
         score = LNodeFw (ln) + LArcTotLike (lat, la);
         switch (type) {
         case LATFB_SUM:
            LNodeFw (la->end) = LAdd (LNodeFw (la->end), score);
            break;
         case LATFB_MAX:
            if (score > LNodeFw (la->end))
               LNodeFw (la->end) = score;
            break;
         default:
            abort ();
         }
      }
   }

   /* backward direction */
   for (i = lat->nn - 1; i >= 0; --i) {
      ln = topOrder[i];
      for (la = ln->pred; la; la = la->parc) {
         assert (la->end == ln);
         
         score = LNodeBw (ln) + LArcTotLike (lat, la);
         switch (type) {
         case LATFB_SUM:
            LNodeBw (la->start) = LAdd (LNodeBw (la->start), score);
            break;
         case LATFB_MAX:
            if (score > LNodeBw (la->start))
               LNodeBw (la->start) = score;
            break;
         default:
            abort ();
         }
      }
   }

   if (trace & T_FB) {
      printf ("forward prob:  %f\n", LNodeFw (topOrder[lat->nn - 1]));
      printf ("backward prob: %f\n", LNodeBw (topOrder[0]));
   }
   score = LNodeBw (topOrder[0]);
   Dispose (&gcheap, topOrder);

   return score;
}

/* EXPORT->LatFindBest

     find the N-best paths (i.e. lowest sum of LArcTotLike()s) and generate
     Transcription.

     ####GE: currently only N=1 is supported
*/
Transcription *LatFindBest (MemHeap *heap, Lattice *lat, int N)
{
   int i;
   LNode *ln;
   LNode **topOrder;
   LArc *la;
   LogDouble score, ac, lm, pr, tot;
   Word nullWord;
   Pron pron;
   Transcription *trans;
   LabList *ll;
   LLink lab;
   
   if (N != 1)
      HError (8690, "FindBest: only 1-best supported, yet.");

   /* during the search ln->score will hold the score of the best
      path to ln (i.e. lowest score). ln->hook will point to
      the arc leading to the preceeding node in this path */

   /* init fields */
   for (i = 0, ln = lat->lnodes; i < lat->nn; ++i, ++ln) {
      ln->score = LZERO;
      ln->hook = NULL;
   }
   LatStartNode (lat)->score = 0.0;

   /* find topological order of nodes */
   topOrder = (LNode **) New (&gcheap, lat->nn * sizeof(LNode *));
   if (!LatTopSort (lat, topOrder))
      HError (8690, "LatFindBest: cannot find best path in Lattice with cycles");

   assert (topOrder[0] == LatStartNode (lat));

   /* traverse nodes in top order */
   for (i = 0; i < lat->nn; ++i) {
      ln = topOrder[i];
      
      /* for all outgoing arcs: propagate scores forward */
      for (la = ln->foll; la; la = la->farc) {
         assert (la->start == ln);
         score = la->start->score + LArcTotLike (lat, la);
         if (score > la->end->score) {
            la->end->score = score;
            la->end->hook = (Ptr) la;
         }
      }
   }

   /* create traqnscription */
   trans = CreateTranscription (heap);
   ll = CreateLabelList (heap, 0);
   
   nullWord = lat->voc->nullWord;

   ac = lm = pr = tot = 0;
   if (trace & T_TRAN)
      printf ("1best trans (wordPron t ac lm pr tot): ");

   /* backtrack from end node along best path and generate transcription */
   ln = LatEndNode (lat);
   la = (LArc *) ln->hook;
   while (la) {
      LabId outlab;

      if (ln->word != nullWord) {
         for (pron = ln->word->pron; pron; pron = pron->next)
            if (pron->pnum == ln->v) 
               break;
         if (pron)
            outlab = pron->outSym;
         else   /* if we can't find pronvar (e.g. wlist), fall back to word */
            outlab = ln->word->wordName;
         if (outlab) {
            lab = CreateLabel (heap, ll->maxAuxLab);
            lab->labid = outlab;
            lab->start = la->start->time * 1.0e7;
            lab->end = ln->time * 1.0e7;
            lab->score = LArcTotLike (lat, la);
            
            lab->succ = ll->head->succ;
            lab->pred = ll->head;
            lab->succ->pred = lab->pred->succ = lab;
         }
      }

      ac += la->aclike;
      lm += la->lmlike;
      pr += la->prlike;
      tot += LArcTotLike (lat, la);
      if (trace & T_TRAN)
         printf ("(%s%d %.2f %.3f %.3f %.3f %.3f) ", ln->word->wordName->name, ln->v,
                 ln->time, la->aclike, la->lmlike, la->prlike, LArcTotLike (lat, la));


      ln = la->start;
      la = (LArc *) ln->hook;
   }
   AddLabelList (ll, trans);

   if (trace & T_TRAN) {
      printf ("\n");
      printf ("ac lm pr tot: %.3f %.3f %.3f %.3f\n", ac, lm, pr, tot);
   }

   Dispose (&gcheap, topOrder);
   return trans;
}


/* EXPORT->LatPrune

     prune lattice by removing all paths with score more than 
     thresh lower than the best path.
*/
Lattice *LatPrune (MemHeap *heap, Lattice *lat, LogDouble thresh, float arcsPerSec)
{
   LogDouble best, limit, score;
   LNode *ln, *newln;
   LArc *la, *newla;
   int i, nn, na;
   Lattice *newlat;

   LatAttachInfo (&gcheap, sizeof (FBinfo), lat);

   best = LatForwBackw (lat, LATFB_MAX);
   limit = best - thresh;
   
   /* modify thresh according to arcPerSec limit */
   if (arcsPerSec > 0) {
#define NBIN 1000
      HTime length;
      int nArc, nArcLimit, bin;
      int hist[NBIN];
      float binWidth;

      length = LatEndNode (lat)->time;
      nArcLimit = length * arcsPerSec;

      binWidth = thresh / NBIN;

      for (i = 0; i < NBIN; ++i)
         hist[i] = 0;

      nArc = 0;
      for (i = 0, la = lat->larcs; i < lat->na; ++i, ++la) {
         score = LNodeFw (la->start) + LArcTotLike (lat, la) + LNodeBw (la->end);
         bin = (best - score) / binWidth;
         assert (bin >= 0);
         if (bin < NBIN) {     /* keep */
            ++hist[bin];
            ++nArc;
         }
      }
      
      if (nArc > nArcLimit) {
         nArc = 0;
         for (i = 0; i < NBIN; ++i) {
            nArc += hist[i];
            if (nArc > nArcLimit)
               break;
         }
         thresh = i * binWidth;
         limit = best - thresh;
         if (trace &T_PRUN)
            printf ("length %.2f nArcLimit %d beam adjusted to %f\n", length, nArcLimit, thresh);
      }
   }

   nn = na = 0;
   /* scan arcs and count survivors */
   for (i = 0, la = lat->larcs; i < lat->na; ++i, ++la) {
      score = LNodeFw (la->start) + LArcTotLike (lat, la) + LNodeBw (la->end);
      if (score >= limit) {     /* keep */
         la->score = (float) na;
         ++na;
      }
      else                   /* remove */
         la->score = -1.0;
   }

   /* scan nodes, count survivors and verify consistency */
   for (i = 0, ln = lat->lnodes; i < lat->nn; ++i, ++ln) {
      score = LNodeFw (ln) + LNodeBw (ln);
      if (score >= limit) {       /* keep */
         ln->n = nn;
         ++nn;
      }
      else {                    
         ln->n = -1;
         /* check that there are no arcs to keep attached to ln */
         for (la = ln->foll; la; la = la->farc) {
            if (la->score >= 0.0)
               HError (8691, "LatPrune: arc score (%f) better than node score (%f)\n", 
                       LNodeFw (la->start) + LArcTotLike (lat, la) +
                       LNodeBw (la->end), score);
         }
      }
   }
   
   if (trace & T_PRUN)
      printf ("lattice pruned from %d/%d to %d/%d\n", lat->nn, lat->na, nn, na);

   /* build new lattice from remaining nodes/arcs */
   newlat = NewILattice (heap, nn, na, lat);

   /* copy all the nodes we want to keep */
   newln = newlat->lnodes;
   for (i = 0, ln = lat->lnodes; i < lat->nn; ++i, ++ln) {
      if (ln->n != -1) {
         *newln = *ln;
         /* init linnked lists */
         newln->foll = newln->pred = NULL;

         ++newln;
      }
   }
   assert (newln == newlat->lnodes + nn);

   newla = newlat->larcs; 
   for (i = 0, la = lat->larcs; i < lat->na; ++i, ++la) {
      if (la->score >= 0.0) {     /* keep */
         *newla = *la;
         newla->start = newlat->lnodes + ((int) la->start->n);
         newla->end = newlat->lnodes + ((int) la->end->n);
         
         /* insert newla into foll list */
         newla->farc = newla->start->foll;
         newla->start->foll = newla;
         /* ...and into pred list */
         newla->parc = newla->end->pred;
         newla->end->pred = newla;

         ++newla;
      }
   }
   assert (newla == newlat->larcs + na);

   LatDetachInfo (&gcheap, lat);

   return newlat;
}


/* StatsInfo structure attached to all LNodes */
typedef struct _StatsInfo {
   LogDouble nPaths;     /* number of paths from start node */
} StatsInfo;

#define LNodeStats(ln)  (((StatsInfo *) (ln)->hook))


/* CalcStats

     calculate and output some global statistics for a lattice
*/
void CalcStats (Lattice *lat)
{
   LNode **topOrder;
   LNode *ln;
   LArc *la;
   int i, d, max_inDegree, max_outDegree, nWords;
   LogDouble nPaths;
   Boolean isDAG;
   LNode *lnStart, *lnEnd;
   Word word;

   /* attach StatsInfo structre to nodes */
   LatAttachInfo (&gcheap, sizeof (StatsInfo), lat);
   for (i = 0, ln = lat->lnodes; i < lat->nn; ++i, ++ln)
      LNodeStats(ln)->nPaths = 0;

   /* find topological order of nodes */
   topOrder = (LNode **) New (&gcheap, lat->nn * sizeof(LNode *));
   isDAG = LatTopSort (lat, topOrder);

   lnStart = isDAG ? topOrder[0] : LatStartNode (lat);
   lnEnd = isDAG ? topOrder[lat->nn-1] : LatEndNode (lat);

   max_inDegree = max_outDegree = 0;
   LNodeStats(lnStart)->nPaths = 1;

   /* reset word counters */
   for (i = 0; i< VHASHSIZE; i++)
      for (word = lat->voc->wtab[i]; word != NULL; word = word->next)
         word->aux = (Ptr) 0;

   /* iterate over all nodes */
   for (i = 0; i < lat->nn; ++i) {
      ln = topOrder[i];

      /* count words */
      ln->word->aux = (Ptr) (((int)ln->word->aux) + 1);

      /* count incoming and outgoing arcs */
      d = 0;
      for (la = ln->pred; la; la = la->parc)
         ++d;
      if (d > max_inDegree)
         max_inDegree = d;

      d = 0;
      for (la = ln->foll; la; la = la->farc)
         ++d;
      if (d > max_outDegree)
         max_outDegree = d;
      
      if (isDAG) {
         /* propagate nPaths forward */
         nPaths = LNodeStats(ln)->nPaths;
         for (la = ln->foll; la; la = la->farc)
            LNodeStats(la->end)->nPaths += nPaths;
      }
   }

   /* find number of words seen in lattice */
   nWords = 0;
   for (i = 0; i < VHASHSIZE; i++)
      for (word = lat->voc->wtab[i]; word != NULL; word = word->next) {
         if (word->aux)
            ++nWords;
         word->aux = (Ptr) 0;
      }


   nPaths = LNodeStats(topOrder[lat->nn-1])->nPaths;

   printf("length[s]      %.2f\n", lnEnd->time);
   printf("ArcsPerSec     %.2f\n", lat->na/lnEnd->time);
   printf("nodes          %d\n", lat->nn);
   printf("arcs:          %d\n", lat->na);
   printf("max_inDegree   %d\n", max_inDegree);
   printf("max_outDegree  %d\n", max_outDegree);
   printf("nWords         %d\n", nWords);
   if (isDAG)
      printf("nPaths         %e.0\n", nPaths);
   else
      printf("nPaths         inf\n");

   Dispose (&gcheap, topOrder);
}


/* LatSetBoundaryWords

     set start and end words for use in lattices and LM
*/
void LatSetBoundaryWords (char *start, char *end, char  *startLM, char *endLM)
{
   startWord = GetLabId (start ? start : "!SENT_START", TRUE);
   startLMWord = GetLabId (startLM ? startLM : "<s>", TRUE);
   endWord = GetLabId (end ? end : "!SENT_END", TRUE);
   endLMWord = GetLabId (endLM ? endLM : "</s>", TRUE);
   nullWord = GetLabId ("!NULL", TRUE);
}

/* FixBadLat

     if final word is !NULL replace it with endWord. This happens in lattices
     produced with some older decoders.
*/
void FixBadLat (Lattice *lat)
{
   LNode *ln;
   LArc *la;
   Boolean seenSE, seenOther;
   Word end;

   LatCheck (lat);

   end = GetWord (lat->voc, endWord, FALSE);
   if (!end)
      HError (8623, "HLRescore: SentEnd word (%d) not in vocabulary", endWord->name);

   ln = LatEndNode (lat);

   if (ln->word == end)
      return;

   if (ln->word && ln->word != lat->voc->nullWord)
      HError (8624, "HLRescore: final word in lattice (%s) is not !NULL! or sentend", 
              ln->word->wordName->name);

   /* now we know that we have !NULL at the end */

   seenSE = seenOther = FALSE;

   for (la = ln->pred; la; la = la->parc) {
      if (la->start->word == end)
         seenSE = TRUE;
      else
         seenOther = TRUE;
   }

   if (seenSE && seenOther)
      HError (8624, "HLRescore: lattice contains both SentEnd and other words in final position");

   if (!seenOther)
      return;

   ln->word = end;
   ln->v = 1;
}


#ifndef NO_LAT_LM
/* LatLMTrans

     wrapper around HLM:LMTrans() to take start and end words into account
*/
static LogFloat  LatLMTrans (LModel *lm, LMState src, LabId wordId, LMState *dest)
{
   LogFloat prob;
   
   if (wordId == nullWord) {
      dest = NULL;
      return 0.0;
   }
   else if (wordId == startWord && src == NULL) {
      /* always return P(<s>) = 1.0 at start of sentence */
      wordId = startLMWord;
      prob = LMTrans (lm, src, wordId, dest);
      return 0.0;
   }
   else if (wordId == endWord) {
      /* destination state of </s> transition is always NULL */
      wordId = endLMWord;
      prob = LMTrans (lm, src, wordId, dest);
      *dest = NULL;
      return prob;
   }
   else {
      prob = LMTrans (lm, src, wordId, dest);
   }

   return prob;
}

/* FindAddSubLNode

     Search for SubLNode in chain and add it if necessary
*/
static SubLNode *FindAddSubLNode (MemHeap *heap, LNode *ln, LMState lmstate, int *nsln)
{
   SubLNode *subln;
   
   for (subln = (SubLNode *) ln->hook; subln; subln = subln->next) {
      if (subln->data.lmstate == lmstate)
         return subln;
   }
   if (!subln) {
      ++*nsln;
      subln = New (heap, sizeof (SubLNode));
      subln->data.lmstate = lmstate;
      subln->foll = NULL;
      subln->next = (SubLNode *) ln->hook;
      ln->hook = (Ptr) subln;
   }
   return subln;
}


/* EXPORT->LatExpand

     expand lattice using new (typically higher-order) language Model
*/
Lattice *LatExpand (MemHeap *heap, Lattice *lat, LModel *lm)
{
   int i, nsln, nsla;
   LNode *ln, *newln;
   LNode **topOrder;
   LArc *la, *newla;
   SubLNode *startSLN, *endSLN, *sln;
   SubLArc *sla;
   LogFloat lmprob;
   LMState dest;
   Lattice *newlat;

   nsln = nsla = 0;

   /* The idea of this algorithm is that we will split each node and arc
      in the lattice into multiple sub-nodes and sub-arcs as required
      by the new LM. 
      N.B.We never join existing nodes, i.e. Calling LatExpand() with a
      bigram on a fourgram lattice will not reduce the size of the lattice.
   */


   /* for each node in the lattice we keep a linked list (hung of
      ln->hook) of sub-nodes (corresponding to LMStates in the new
      LM). */

   /* init sub-node linked lists */
   for (i = 0, ln = lat->lnodes; i < lat->nn; ++i, ++ln) {
      ln->hook = NULL;
   }

   /* create one sub-node for lattice start node with LMState = NULL */
   FindAddSubLNode (&slnHeap, LatStartNode (lat), NULL, &nsln);
   
   /* find topological order of nodes */
   topOrder = (LNode **) New (&gcheap, lat->nn * sizeof(LNode *));
   LatTopSort (lat, topOrder);

   
   /* create lists of sub-nodes and sub-arcs and count them as we go along */
   for (i = 0; i < lat->nn; ++i) {
      ln = topOrder[i];
      for (startSLN = (SubLNode *) ln->hook; startSLN; startSLN = startSLN->next) {
         /* for each outgoing arc from current subLNode */
         for (la = ln->foll; la; la = la->farc) {
            assert (la->start == ln);
            lmprob = LatLMTrans (lm, startSLN->data.lmstate, la->end->word->wordName, &dest);
            endSLN = FindAddSubLNode (&slnHeap, la->end, dest, &nsln);

            /* add new subLArc */
            ++nsla;
            sla = New (&slaHeap, sizeof (SubLArc));
            sla->lmprob = lmprob;
            sla->end = endSLN;
            sla->la = la;
            /* add to list of arcs leaving startSLN */
            sla->next = startSLN->foll;
            startSLN->foll = sla;
         }
      }
   }

   if (trace & T_EXP)
      printf ("expanded lattice from %d/%d  to %d/%d\n", lat->nn, lat->na, nsln, nsla);
   
   /* build new lattice from sub-node/-arc lists */
   newlat = NewILattice (heap, nsln, nsla, lat);
   newlat->net = CopyString (heap, lm->name);

   /* create one node in new lattice for each sub node in old lattice */
   newln = newlat->lnodes;
   for (i = 0; i < lat->nn; ++i) {
      ln = topOrder[i];
      for (sln = (SubLNode *) ln->hook; sln; sln = sln->next) {
         *newln = *ln;
         newln->foll = newln->pred = NULL;
         newln->n = 0;
         newln->hook = NULL;

         sln->data.newln = newln;
         ++newln;
      }
   }
   assert (newln = newlat->lnodes + newlat->nn);

   /* create arcs in new lattice */
   newla = newlat->larcs;
   for (i = 0; i < lat->nn; ++i) {
      ln = topOrder[i];
      for (sln = (SubLNode *) ln->hook; sln; sln = sln->next) {
         newln = sln->data.newln;
         for (sla = sln->foll; sla; sla = sla->next) {
            *newla = *sla->la;
            newla->start = newln;
            newla->end = sla->end->data.newln;
            newla->lmlike = sla->lmprob;
            
            /* add to start node foll list */
            newla->farc = newla->start->foll;
            newla->start->foll = newla;
            /* add to end node pred list */
            newla->parc = newla->end->pred;
            newla->end->pred = newla;
            
            ++newla;
         }
      }
   }
   assert (newla == newlat->larcs + newlat->na);

   if (trace & T_MEM) {
      printf("Memory State after expanding\n");
      PrintAllHeapStats();
   }

   Dispose (&gcheap, topOrder);
   ResetHeap (&slaHeap);
   ResetHeap (&slnHeap);

   return newlat;
}

#endif

