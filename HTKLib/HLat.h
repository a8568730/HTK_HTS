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
/*       File: HLat.h:  Lattice Manipulation                   */
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

/* !HVER!HLat:   3.3 [CUED 28/04/05] */


#ifndef _HLAT_H_
#define _HLAT_H_

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------ Initialisation --------------------------- */

void InitLat (void);
/*
   register module & set configuration parameters
*/

void ResetLat (void);
/*
   reset the module
*/

/* ------------------------ Datatype ----------------------------- */

/* Forward-Backward info structre attached to LNodes
     use doubles for imporved accuracy
*/
typedef struct FBlnodeInfo {
   LogDouble fwlike;     /* forward likelihood */
   LogDouble bwlike;     /* backward likelihood */
} FBinfo;

#define LNodeFw(ln)  (((FBinfo *) (ln)->hook)->fwlike)
#define LNodeBw(ln)  (((FBinfo *) (ln)->hook)->bwlike)

typedef enum {LATFB_SUM, LATFB_MAX} LatFBType;


/* ------------------------ Prototypes --------------------------- */

Transcription *LatFindBest (MemHeap *heap, Lattice *lat, int N);

void LatSetScores (Lattice *lat);

Lattice *LatPrune (MemHeap *heap, Lattice *lat, LogDouble thresh, float arcsPerSec);

void CalcStats (Lattice *lat);

LNode *LatStartNode (Lattice *lat);

LNode *LatEndNode (Lattice *lat);

void LatSetBoundaryWords (char *start, char *end, char  *startLM, char *endLM);

void LatCheck (Lattice *lat);

void FixBadLat (Lattice *lat);

void FixPronProbs (Lattice *lat, Vocab *voc);

Boolean LatTopSort (Lattice *lat, LNode **topOrder);

void LatAttachInfo (MemHeap *heap, size_t size, Lattice *lat);
void LatDetachInfo (MemHeap *heap, Lattice *lat);

LogDouble LatForwBackw (Lattice *lat, LatFBType type);


#ifndef NO_LAT_LM
Lattice *LatExpand (MemHeap *heap, Lattice *lat, LModel *lm);

void ApplyNGram2LabLat(Lattice *lat, LModel *lm);
#endif


Lattice *GetLattice (char *fn, char *path, char *ext,
                     /* arguments of ReadLattice() below */
                     MemHeap *heap, Vocab *voc, 
                     Boolean shortArc, Boolean add2Dict);

Lattice *MergeLatNodesArcs(Lattice *lat, MemHeap *heap, Boolean mergeFwd);

void ApplyWPNet2LabLat(Lattice *lat, Lattice *wdNet);

#ifdef __cplusplus
}
#endif

#endif  /* _HLAT_H_ */
