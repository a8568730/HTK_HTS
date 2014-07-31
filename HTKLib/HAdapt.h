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
/*         File: HAdapt.h      Adaptation Library module       */
/* ----------------------------------------------------------- */

/* !HVER!HAdapt:   3.2.1 [CUED 15/10/03] */

#ifndef _HADAPT_H_
#define _HADAPT_H_

#ifdef __cplusplus
extern "C" {
#endif

enum _RegClassType{ADPTUNDEF=0,ADPTFIXED=1,ADPTTREE=2};
typedef enum _RegClassType RegClassType;

enum _RegTransType{TRANSUNDEF=0,MEANONLY=1,MEANVAR=2};
typedef enum _RegTransType RegTransType;


#define DEF_BLOCKSIZE          1   /* define the default blocksize */
#define DEF_REGCLASS    ADPTTREE   /* define the default Reg Class Type */

/* ------------- Block Matrix Definitions -------------- */

/* A block matrix is a array of square matrices */
/* There are nBlocks matrices, 
   and each matrix is a square matrix of blockSize */

/* [1..nBlocks] of matrix[blockSize][blockSize] */

typedef Matrix  *BlockMatrix;
typedef DMatrix *BlockDMatrix;
typedef TriMat  *BlockTriMatrix;

/* ------------ End Block Matrix Definitions ------------ */

/* structure to hold the mean regression matrix & bias vector: x' = ax + b
   plus the Z accumumulation structure */
typedef struct {

  BlockMatrix a;            /* block diagonal matrix */
  Vector b;                 /* bias vector           */

} OffsetBMat;

/* specialist structure for the G matrix ! */
typedef struct {

  BlockTriMatrix a;         /* block tri-diagonal matrix */
  Vector b;                 /* bias vector               */

} OffsetTriBMat;


typedef struct {

  Boolean speechFlag;       /* flag to see if class is speech/non-speech */
  short baseformIndex;      /* regression class index number */
  int nComponents;          /* number of components in this class */
  float  occ;               /* baseform's occ count for ALL incoming data */  
  Vector H;                 /* the H accumulates for base class (covar) */
  OffsetTriBMat **G;        /* the G accumulates for base class */
  OffsetBMat *Z;            /* the Z accumulates for base class */
  MixPDF **pdfList;         /* array list of pointers to pdfs in base class */

} BaseformAcc;


typedef struct {

  short        nBases;      /* number of baseforms for this node */
  short     nodeIndex;      /* node's index number */
  short        *bases;      /* list of baseform indices for this node */
  int       nodeComps;      /* the number of components present at a node */  
  float       nodeOcc;      /* node's occupation count */    
  OffsetBMat  *WTrans;      /* mean matrix transform for node */ 
  Vector       HTrans;      /* diag cov vector transform for node */
  Matrix    backTrans;      /* inverse mean transform matrix */

} RegNode ;

/* structure to hold the regression tree :-
   a) Global transform has all stats stored at the root node 
   b) Fixed transform has only left children (all right children are NULL)
   c) Standard regression tree */
typedef struct _RegTree {
  
  struct _RegTree *left;
  struct _RegTree *right;
  RegNode *nodeInfo;
  Ptr node;

} RegTree;

/* structure to hold information pertaining to this transform set */
typedef struct {

  char *uid;                /* user identifier */ 
  char *name;               /* user's full name */
  char *hmmSetId;           /* hmm set identifier */
  char *rcId;               /* regression class type identifier */
  char *chan;               /* desc of channel used for data collection */
  char *desc;               /* general description */

} TransformId;

/* structure to hold the regression classes and accumulation information */
typedef struct {

  TransformId *transId;     /* information identifier for this transform */
  HMMSet *hset;             /* HMMSet being transformed */
  MemHeap *hmem;            /* mem heap for HAdapt module */
  MemHeap pdfStatsMem;      /* mem heap solely for the mixPDF accumulations */
  MemHeap transMem;         /* mem heap solely for the transform structures */
  short vSize;              /* vector size */
  short nBlocks;            /* number of blocks in each transform */ 
  short nBaseTransforms;    /* total number of possible transforms   */
  RegClassType classKind;   /* type of regression class */
  RegTransType transKind;   /* type of regression transformations */
  TriState adptSil;         /* Flag for adapting silence */
  float nodeOccThresh;      /* Node occupation threshold for 
			       generating a transform */
  /* Boolean transWithBias; Use of loaded bias */
  /*  TriState useBias;     Use of bias */
  RegTree *rtree;           /* regression tree */
  BaseformAcc **baseforms;  /* list of baseform accumulates */

} RegTransInfo ;

/* structure to store the regression accumulations */
typedef struct {

  /* in place of PreComp -- necessary for f-b routines */
  /* DO NOT MOVE or REMOVE !!! */
  PreComp  pre;             /* PreComputed Mixture calculation structure */

  /* the rest */
  float    occ;             /* occupation count */
  Vector   obsSum;          /* summed speech vectors, scaled by posterior */
  Vector   obsSqSum;        /* summed sqared sp. vecs, scaled by posterior */
  /* MixPDF   *nextMix;      next mixture in the regression class list */

} RegAcc ;

/* -------------------- Initialisation Functions -------------------------- */

void InitAdapt(void);
/*
   Initialise configuration parameters
*/

void InitialiseTransform( HMMSet *hset, MemHeap *x, RegTransInfo *rt,
			  Boolean adapt );
/*
   Initialise transforms storage and grouping of hmmset components
   to regression base class using linked lists 
*/

void InitialiseAdapt(HMMSet *hset,MemHeap *x, RegTransInfo *rt);
/* 
   Initialise adaptation storage
*/



/* -------------------- Save/Load Transform Functions --------------------- */

void SaveTransformSet(HMMSet *hset, RegTransInfo *rt, char *saveFile, 
		      char *transPath, char *uid, char *uname, char *chan, 
		      char *desc, Boolean saveStats, Boolean global,
		      Boolean saveBinary);
/* 
   Save the transforms to a file, and also optionally save
   transform statistics allowing continued adaptation
*/

Boolean LoadTransformSet(HMMSet *hset, char *tfile, char *uid,
			 RegTransInfo *rt, Boolean *useStats);
/* 
   Load the transforms into memory; if useStats is requested and 
   transform does not contain any statistics, then useStats is
   reset to FALSE -- returns FALSE if tmf is unsuitable
   for current HMMSet
*/


/* ---------------- Accumulation Control Functions ------------------------ */

void AccAdaptFrame(double Lr, Vector speechVec, MixPDF *mp, RegTransInfo *rt);
/* 
   Accumulate frame stats into specific mixture comp
*/

void ClearRegCompStats(HMMSet *hset, RegTransInfo *rt);
/* 
   Clear for regression level accumulated stats 
*/

void ClearBaseClassStats(RegTransInfo *rt);
/* 
   Clear for base class level accumulated stats 
*/

/* ----------------- Mixture Transformation Functions --------------------- */

/* ----- Transform application functions ---- */

void ApplyMeanTransforms(RegTransInfo *rt, RegTree *t);
/* 
   Apply the mean transform to the regression classes
*/

void ApplyCovarTransforms(RegTransInfo *rt, RegTree *t);
/*
   Apply the covariance transform to the regression classes
*/

Boolean ApplyMeanGlobalTransform(RegTransInfo *rt);
/* 
   Apply the mean Global transformation -- returns true if global is found 
*/

Boolean ApplyCovarGlobalTransform(RegTransInfo *rt);
/*
   Apply the covar Global transformation -- returns true if global is found 
*/

void ApplyTransforms(RegTransInfo *rt);
/* 
   Apply the mean and possibly the covar MLLR transformations
*/

Boolean ApplyBackwardGlobalTransform(RegTransInfo *rt);
/* 
   Apply the backward global transformation to the model set 
   to obtain the "original" models -- returns TRUE if global is found
*/

void ApplyBackwardTransforms(RegTransInfo *rt);
/*
   Apply the backward transformation to the model set 
   to obtain the "original" models
*/

/* ----- Transform calculation/update ----- */

void CalcMeanTransforms(RegTransInfo *rt);
/* 
   Calculate the mean transforms for every regression class 
*/

void CalcCovarTransforms(RegTransInfo *rt);
/* 
   Calculate the covariance transforms for every regression class
*/

void GetBackwardTransforms(RegTransInfo *rt);
/* 
   Calculates the backward transformation by inverting 
   the original transform
*/

void UpdateMAP(RegTransInfo *rt, float tau);
/* 
   Update the models using Maximum A-Posteriori Adaptation (MAP) 
*/

/* ---------------- General adaptation Call ------------------------------- */

/* EXPORT->DoAdaptation: Given the initialisation to the adaptation */

/* MLLR Adaptation 
   Before this function is called HFB must be called via "FBFile"
   to perfrom a frame/state alignment and accumulate regression statistics
   at the mixture level 

   DoAdaptation does:-
   1) Transform current models back to their original state (as loaded)
   2) Calculation mean transform and apply the transform
   3) Calculate the covariance transform (if requested) and apply
*/

void DoAdaptation(RegTransInfo *rt, Boolean global);

#ifdef __cplusplus
}
#endif

#endif  /* _HADAPT_H_ */

/* ---------------------------- END HAdapt.h ------------------------------ */

