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
/*         File: HAdapt.c      Adaptation Library module       */
/* ----------------------------------------------------------- */

char *hadapt_version = "!HVER!HAdapt:   3.2.1 [CUED 15/10/03]";
char *hadapt_vc_id =  "$Id: HAdapt.c,v 1.10 2003/10/15 08:10:12 ge204 Exp $";


#include <stdio.h>      /* Standard C Libraries */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "HShell.h"     /* HMM ToolKit Modules */
#include "HMem.h"
#include "HMath.h"
#include "HSigP.h"
#include "HWave.h"
#include "HAudio.h"
#include "HParm.h"
#include "HLabel.h"
#include "HModel.h"
#include "HTrain.h"
#include "HUtil.h"
#include "HAdapt.h"

/* trace flags */
#define T_TOP   00001    /* Top level tracing */
#define T_NAC   00002    /* trace on accumulates */
#define T_TRA   00004    /* trace on transformations */
#define T_AUX   00010    /* output the auxilliary function */
#define T_RIO   00020    /* regression classes input output */
#define T_USE   00040    /* regression classes tree usage */
#define T_DET   00200    /* detailed trace for class level accumulates */
#define T_MEM   00400    /* trace to print heap stats */

#define MINOCC        1E-08
#define OCC_THRESH    700.0

static int trace = 0;               /* trace info */
static Boolean initVar=FALSE;       /* use a variance transform */
static Boolean initAdptSil=TRUE;    /* adapt the silence */
static Boolean writeBin=FALSE;      /* write file in binary format */
static int blocks = DEF_BLOCKSIZE;  /* number of blocks initialisation */
static RegClassType regClass=DEF_REGCLASS; /* regClass type to be used */
static float occThreshInit=OCC_THRESH;   /* tree ocupation count minimum */
static char *undefined="Undefined"; /* undefined tmf header field variable */

/* matrices and vectors for SVD */
static DMatrix svdMat;             /* general double SVD matrix */
static DMatrix u, v;               /* temporary matrices for SVD */
static DVector w;                  /* temporary vector for SVD */
static DMatrix svdMat1;            /* general double SVD matrix for btrans */
static DMatrix u1, v1;             /* temp matrices for SVD for backtrans*/
static DVector w1;                 /* temp vector for SVD for backtrans */

/* matrices and vectors for accumulating regression
   statistics at the base class level. */
static Matrix *Gg;
static Matrix Zg;
static Vector Hg;
static Matrix Wg;

static DMatrix blkZ;
static DMatrix blkG;
static DMatrix blku;
static DMatrix blkv;
static DVector blkw;

ConfParam *cParm[MAXGLOBS];      /* config parameters */
int nParm = 0;

/* ----------------------------------------------------------------------*/
/*                 Block Diagonal Matrix Functions                       */
/* ----------------------------------------------------------------------*/

/* --------------- Sizing funtions ---------------------- */

/* Return number of blocks in a block matrix */
static int GetMatBlockSize(BlockMatrix B)
{
   return ((int) B[0]);
}

/* Return number of blocks in a block double matrix */
static int GetDMatBlockSize(BlockDMatrix B)
{
   return ((int) B[0]);
}

/* Return number of blocks in a block tri diagonal matrix */
static int GetTriMatBlockSize(BlockTriMatrix B)
{
   return ((int) B[0]);
}


/*------------- Create the block matrix structures --------------------*/

/* Allocate space for a block matrtix 
   B[1..nBlocks][1..blockSize][1..blockSize] */
static BlockMatrix CreateBlockMat(MemHeap *x, int vSize, int nBlocks) 
{

   int blockSize = 0;
   BlockMatrix B;
   int i;

   /* firstly do check */
   blockSize = vSize / nBlocks;
   if (blockSize * nBlocks != vSize)
      HError(7460, "CreateBlockMat: Number of blocks incompatible with vector size!");
  
   B = (BlockMatrix) New(x, (nBlocks+1)*sizeof(Matrix));
   B[0] = (Matrix) nBlocks;
   for (i = 1; i <= nBlocks; i++)
      B[i] = CreateMatrix(x, blockSize, blockSize);
  
   return B;

}

/* Allocate space for blockdouble matrix 
   B[1..nBlocks][1..blockSize][1..blockSize] */
static BlockDMatrix CreateBlockDMat(MemHeap *x, int vSize, int nBlocks) 
{

   int blockSize = 0;
   BlockDMatrix B;
   int i;

   /* firstly do check */
   blockSize = vSize / nBlocks;
   if (blockSize * nBlocks != vSize)
      HError(7460, "CreateBlockDMat: Number of blocks incompatible with vector size!");
  
   B = (BlockDMatrix) New(x, (nBlocks+1)*sizeof(DMatrix));
   B[0] = (DMatrix) nBlocks;
   for (i = 1; i <= nBlocks; i++)
      B[i] = CreateDMatrix(x, blockSize, blockSize);
  
   return B;

}

/* Allocate space for tridiagonal matrix 
   B[1..nBlocks][1..blockSize][1..blockSize] */
static BlockTriMatrix CreateBlockTriMat(MemHeap *x, int vSize, int nBlocks) 
{

   int blockSize = 0;
   BlockTriMatrix B;
   int i;

   /* firstly do check */
   blockSize = vSize / nBlocks;
   if (blockSize * nBlocks != vSize)
      HError(7460, "CreateBlockTriMat: Number of blocks incompatible with vector size!");
  
   B = (BlockTriMatrix) New(x, (nBlocks+1)*sizeof(TriMat));
   B[0] = (TriMat) nBlocks;
   for (i = 1; i <= nBlocks; i++)
      B[i] = CreateTriMat(x, blockSize);
  
   return B;

}

/* ----------- Free / delete block matrices ------------- */ 

/* Free memory allocated for a block matrix */
static void FreeBlockMat(MemHeap *x, BlockMatrix B)
{  
   int i, nBlocks;

   nBlocks = GetMatBlockSize(B);
   for (i = 1; i <= nBlocks; i++)
      FreeMatrix(x, B[i]);

   Dispose(x, B);
}

/* Free memory  allocated for block double matrix */
static void FreeBlockDMat(MemHeap *x, BlockDMatrix B)
{
   int i, nBlocks;

   nBlocks = GetDMatBlockSize(B);
   for (i = 1; i <= nBlocks; i++)
      FreeDMatrix(x, B[i]);

   Dispose(x, B);
}

/* Free memory  allocated for block tridiagonal matrix */
static void FreeBlockTriMat(MemHeap *x, BlockTriMatrix B)
{
   int i, nBlocks;

   nBlocks = GetTriMatBlockSize(B);
   for (i = 1; i <= nBlocks; i++)
      FreeTriMat(x, B[i]);

   Dispose(x, B);

}

/* ------------- Zero block matrices ------------- */

/*  Zero a block matrix */
static void ZeroBlockMat(BlockMatrix B) 
{

   int i, nBlocks;

   nBlocks = GetMatBlockSize(B);
   for (i = 1; i <= nBlocks; i++)
      ZeroMatrix(B[i]);
  
}

/* Zero a block double matrix */
static void ZeroBlockDMat(BlockDMatrix B) 
{

   int i, nBlocks;

   nBlocks = GetDMatBlockSize(B);
   for (i = 1; i <= nBlocks; i++)
      ZeroDMatrix(B[i]);
  
}

/* Zero a block tri-diagonal matrix */
static void ZeroBlockTriMat(BlockTriMatrix B) 
{

   int i, nBlocks;

   nBlocks = GetTriMatBlockSize(B);
   for (i = 1; i <= nBlocks; i++)
      ZeroTriMat(B[i]);
  
}

/* -------- Arithmetic functions with block matrices ------ */

/* Multiply  a vector by a block matrix */ 
static void MultBlockMat_Vec(BlockMatrix B, Vector in, Vector out) 
{
   int vSize, bSize, nBlocks;
   int i, j, k, fj, fk, blockStart;
   Matrix m;
   Vector v;


   /* do check */
   vSize = VectorSize(in);
   nBlocks = GetMatBlockSize(B);
   bSize = 1;
   if (nBlocks >= 1)
      bSize = NumRows(B[1]);
   else
      HError(7460, "MultBlockMat_Vec: Block matrix has less than 1 block!");

   if (bSize * nBlocks != vSize)
      HError(7460, "MultBlockMat_Vec: Incompatible multiplication %d blocks, %d block size, %d full vector size\n", nBlocks, bSize, vSize);
  
   /* allocate space for temporary vector v */
   v = CreateVector(&gstack, vSize);
   ZeroVector(v);

   blockStart = 0;
   for (i = 1; i <= nBlocks; i++) {
      m = B[i];
      for (j = 1; j <= bSize; j++) {
         fj = j + blockStart;
         for (k = 1; k <= bSize; k++) {
            fk = k + blockStart;
            v[fj] += m[j][k] * in[fk];
         }
      }
      blockStart += bSize;
   }

   CopyVector(v, out);
   FreeVector(&gstack, v);

}

/* Multiply a block matrix by a block matrix */ 
static void MultBlockMat_BlockMat(BlockMatrix B1, BlockMatrix B2, 
                                  BlockMatrix out) 
{
   int vSize, bSize, nBlocks;
   int n, i, j, k;
   BlockMatrix m;


   /* do check */
   nBlocks = GetMatBlockSize(B1);
   bSize = NumRows(B1[1]);
   vSize = nBlocks * bSize;

   /* allocate space for temporary vector v */
   m = CreateBlockMat(&gstack, vSize, nBlocks);
   ZeroBlockMat(m);

   for (n = 1; n <= nBlocks; n++) {
      for (i = 1; i <= bSize; i++)
         for (j = 1; j <= bSize; j++)
            for (k = 1; k <= bSize; k++)
               m[n][i][j] += B1[n][i][k] * B2[n][k][j];
      for (i = 1; i <= bSize; i++)
         for (j = 1; j <= bSize; j++)
            out[n][i][j] = m[n][i][j];
   }
  
   FreeBlockMat(&gstack, m);

}

/* Multiply a vector by a block double matrix */ 
static void MultBlockDMat_Vec(BlockDMatrix B, Vector in, Vector out) 
{
   int vSize, bSize, nBlocks;
   int i, j, k, fj, fk, blockStart;
   DMatrix m;
   Vector v;

   bSize = 1;
   /* do check */
   vSize = VectorSize(in);
   nBlocks = GetDMatBlockSize(B);
   if (nBlocks >= 1)
      bSize = NumDRows(B[1]);
   else
      HError(7460, "MultBlockDMat_Vec: Block matrix has less than 1 block!");

   if (bSize * nBlocks != vSize)
      HError(7460, "MultBlockDMat_Vec: Incompatible multiplication %d blocks, %d block size, %d full vector size\n", nBlocks, bSize, vSize);
  
   /* allocate space for temporary vector v */
   v = CreateVector(&gstack, vSize);
   ZeroVector(v);

   blockStart = 0;
   for (i = 1; i <= nBlocks; i++) {
      m = B[i];
      for (j = 1; j <= bSize; j++) {
         fj = j + blockStart;
         for (k = 1; k <= bSize; k++) {
            fk = k + blockStart;
            v[fj] += m[j][k] * in[fk];
         }
      }
      blockStart += bSize;
   }

   CopyVector(v, out);
   FreeVector(&gstack, v);

}


/*------------------------------------------------------------------------*/
/*            Initialisations and general structures allocation           */
/*------------------------------------------------------------------------*/

/* EXPORT->InitAdapt: initialise configuration parameters */
void InitAdapt (void) 
{
   int i;
   Boolean b;
   double d;

   Register(hadapt_version,hadapt_vc_id);
   nParm = GetConfig("HADAPT", TRUE, cParm, MAXGLOBS);
  
   if (nParm>0){
      if (GetConfInt(cParm,nParm,"TRACE",&i)) trace = i;
      if (GetConfBool(cParm,nParm,"USEVAR",&b)) initVar = b;
      if (GetConfBool(cParm,nParm,"ADPTSIL",&b)) initAdptSil = b;
      if (GetConfInt(cParm,nParm,"BLOCKS",&i)) blocks = i;
      if (GetConfBool(cParm,nParm,"SAVEBINARY",&b)) writeBin = b;
      if (GetConfFlt(cParm,nParm,"OCCTHRESH",&d)) occThreshInit = (float) d;
   }

}

/* Create space for the mean transformation matrix and offset vector */
static OffsetBMat *CreateMeanTransform(MemHeap *x, int size, int nBlocks) 
{
   OffsetBMat *m;
  
   m = (OffsetBMat *) New(x, sizeof(OffsetBMat));
  
   m->a = CreateBlockMat(x, size, nBlocks);
   m->b = CreateVector(x, size);
   ZeroBlockMat(m->a);
   ZeroVector(m->b);
  
   return m;
  
}

/* create space for a transformation */
static void CreateNodeTransform(RegNode *n, MemHeap *x, 
                                RegTransType transKind, int vSize, int nBlocks)
{
   int i;

   n->WTrans = CreateMeanTransform(x, vSize, nBlocks);

   n->HTrans = CreateVector(x, vSize);
   for (i = 1; i <= vSize; i++)
      n->HTrans[i] = 1.0; 

}

/* Allocate space for a regression node of the regression tree */
static void CreateNodeInfo(RegNode *n, RegTransInfo *rt, Boolean fixed)
{

   n->nodeOcc = 0.0;
   n->nodeComps = 0;
   n->WTrans    = NULL;
   n->backTrans = NULL;
   n->HTrans    = NULL;

   /* check if fixed or global adaptation */
   if (fixed) {
      n->nBases = 1;
      n->nodeIndex = 1;
      n->bases = CreateShortVec(rt->hmem, n->nBases);
   }

}

/* Allocate space for a baseform (regression class) structure */
static BaseformAcc *CreateBaseform(RegTransInfo *rt) {

   BaseformAcc *b;
   int vSize;
  
   vSize = rt->vSize;
   b = (BaseformAcc *) New(rt->hmem, sizeof(BaseformAcc));
   b->nComponents = 0;
   b->baseformIndex = 0;
   b->speechFlag = TRUE;
   b->occ = 0.0;
   b->H = NULL;
   b->Z = NULL;
   b->G = NULL;
   b->pdfList = NULL;
 
   return(b);
}

/* Allocate space for a baseform (regression class) accumulators G, Z and H */
static void CreateBaseformAccStorage(BaseformAcc *b, RegTransInfo *rt)
{
   int j, vSize;
   OffsetBMat *m;
   OffsetTriBMat *tm;

   vSize = rt->vSize;

   b->H = CreateVector(rt->hmem, vSize);
   ZeroVector(b->H);

   m = (OffsetBMat *) New(rt->hmem, sizeof(OffsetBMat));
   m->a = CreateBlockMat(rt->hmem, vSize, rt->nBlocks);
   m->b = CreateVector(rt->hmem, vSize);
   ZeroBlockMat(m->a);
   ZeroVector(m->b);
   b->Z = m;

   b->G = (OffsetTriBMat **) New(rt->hmem, sizeof(OffsetTriBMat *)*vSize);
   --(b->G);
   for (j = 1; j <= vSize; j++) {
      tm = (OffsetTriBMat *) New(rt->hmem, sizeof(OffsetTriBMat));
      tm->a = CreateBlockTriMat(rt->hmem, vSize, rt->nBlocks);
      tm->b = CreateVector(rt->hmem, vSize+1);
      ZeroBlockTriMat(tm->a);
      ZeroVector(tm->b);
      b->G[j] = tm;
   }
}

/* allocate memory for the transforms that are stored at the nodes */
static void AllocNodeTransforms(RegTree *t, RegTransInfo *rt)
{
 
   if (t != NULL) {
      AllocNodeTransforms(t->left, rt);
      CreateNodeInfo(t->nodeInfo, rt, FALSE);
      AllocNodeTransforms(t->right, rt);
   }

}

/*------------------------------------------------------------------------*/
/*                    General purpose functions                           */
/*------------------------------------------------------------------------*/

/* add symmetric OffsetTriBMat G to Matrix m2 (lower triangle!) */
static void AddGSymMatrix(OffsetTriBMat *G,  Matrix m2) {

   int i, j, k, fk, nBlocks, bSize, bStart;
   int m;
   Vector v1, v2;
   Vector b;
   TriMat m1;

   m = NumRows(m2);
   b = G->b;
   nBlocks = GetTriMatBlockSize(G->a);
   bSize = (m-1)/nBlocks;
   if (bSize * nBlocks != m-1)
      HError(7460, "AddGMatrix: Number of blocks incompatible with vector size!");

   /* add the offset */
   for (i = 1; i <= m; i++)
      m2[i][1] += b[i];

   /* add the matrix */
   bStart = 0;
   for (i = 1; i <= nBlocks; i++) {
      m1 = G->a[i];
      for (j = 1; j <= bSize; j++) {
         v1 = m1[j]; 
         v2 = m2[bStart+j+1];
         for (k = 1; k <= j; k++) {
            fk = bStart+k+1;
            v2[fk] += v1[k];
         }
      }
      bStart+=bSize;
   }
}

/* add OffsetBMat Z to Matrix m2 */
static void AddZMatrix(OffsetBMat *Z,  Matrix m2) {

   int i, j, k, fk, nBlocks, bSize, bStart;
   int m;
   Vector v1, v2;
   Vector b;
   Matrix m1;

   b = Z->b;
   nBlocks = GetMatBlockSize(Z->a);
   m = NumRows(m2);
   bSize = m/nBlocks;
   if (bSize * nBlocks != m)
      HError(7460, "AddZMatrix: Number of blocks incompatible with vector size!");

   /* add the offset */
   for (i = 1; i <= m; i++)
      m2[i][1] += b[i];

   /* add the matrix */
   bStart = 0;
   for (i = 1; i <= nBlocks; i++) {
      m1 = Z->a[i];
    
      for (j = 1; j <= bSize; j++) {
         v1 = m1[j]; 
         v2 = m2[bStart+j];
         for (k = 1; k <= bSize; k++) {
            fk = bStart+k+1;
            v2[fk] += v1[k];
         }
      }
      bStart+=bSize;
   }

}

/* Convert a standard matrix into a block diagonal matrix */
static void ConvertMat2BlockMat(Matrix W, OffsetBMat *M, 
                                int vSize, int nBlocks) 
{
   int i, j, fi, fj, n, bStart, bSize;
   BlockMatrix a;
   Matrix m;
  
   a = M->a;
   bSize = vSize/nBlocks;
   bStart = 0;
   for (n = 1; n <= nBlocks; n++) {
      m = a[n];
      for (i = 1; i<= bSize; i++) {
         fi = i+bStart;
         for (j = 1; j<= bSize; j++) {
            fj = j+bStart;
            m[i][j] = W[fi][fj];
         }
      }
      bStart += bSize;
   }
  
}


/* return the regression accumulator hook (on the mixture PDF) 
   used commonly by many functions in this module */
static RegAcc *GetRegAcc(MixPDF *mp)
{

   return ((RegAcc *) mp->hook);

}

/* removes and frees a mean transform structure */
static void RemoveMeanTransform(MemHeap *x, OffsetBMat *m) 
{

   FreeBlockMat(x, m->a);
   FreeVector(x, m->b);
   Dispose(x, m);

}

/* function that frees up unsued transforms */
static void RemoveUnusedTransforms(RegTree *t, RegTransInfo *rt) 
{
   RegNode *n, *n1, *n2;

   if (t->left != NULL) {

      RemoveUnusedTransforms(t->left, rt);

      n = t->nodeInfo;
      n1 = t->left->nodeInfo;
      n2 = t->right->nodeInfo;

      if (n->WTrans != NULL && t != rt->rtree) {
         if (n1->WTrans != NULL && n2->WTrans != NULL) {
            RemoveMeanTransform(&rt->transMem, n->WTrans);
            n->WTrans = NULL;
            Dispose(&rt->transMem, n->HTrans);
            n->HTrans = NULL;
         }
      }

      RemoveUnusedTransforms(t->right, rt);
    
   }
}

/* function that frees up ALL transforms */
static void RemoveAllTransforms(RegTree *t, RegTransInfo *rt) 
{
   RegNode *n;

   if (t != NULL) {

      RemoveAllTransforms(t->left, rt);

      n = t->nodeInfo;

      if (n->WTrans != NULL) {
         RemoveMeanTransform(&rt->transMem, n->WTrans);
         n->WTrans = NULL;
         Dispose(&rt->transMem, n->HTrans);
         n->HTrans = NULL;
      }

      RemoveAllTransforms(t->right, rt);
   }
}

/* Print the regression tree -- Used by tracing functions */
static void PrintTree(RegTree *t, RegTransInfo *rt)
{
   short idx;
   int i;

   if (t != NULL) {
      PrintTree(t->left, rt);

      if (t->left == NULL) {
         idx = t->nodeInfo->bases[1];
         printf("%d : Terminal (%d, %f)\n", t->nodeInfo->nodeIndex,
                rt->baseforms[idx]->nComponents, t->nodeInfo->nodeOcc);
      }
      else {
         printf("%d : { %d %d } (%f)\n", t->nodeInfo->nodeIndex,
                t->left->nodeInfo->nodeIndex, t->right->nodeInfo->nodeIndex,
                t->nodeInfo->nodeOcc); 
         printf("Bases are ");
         for (i = 1; i <= t->nodeInfo->nBases; i++) {
            idx = t->nodeInfo->bases[i];
            printf("%d ", rt->baseforms[idx]->baseformIndex);
         }
         printf("\n");
         fflush(stdout);
      }
    
      PrintTree(t->right, rt);
   }
  
}

/* Find a node in the regression tree given the index and return the node */
static RegTree *FindNode(RegTree *t, RegTree *r, int i)
{
 
   if (t != NULL) {
      r = FindNode(t->left, r, i);
      if (t->nodeInfo->nodeIndex == i)
         return t;
      r = FindNode(t->right, r, i);
   }

   return r;

}

/* find the number of terminals in the tree (returned in nTerms) */
static void FindNumTerminals(RegTree *t, int *nTerms) 
{
  
   if (t != NULL) {
      FindNumTerminals(t->left, nTerms);
      if (t->left == NULL)
         *nTerms += 1;
      FindNumTerminals(t->right, nTerms);
   }

}

/* Given a base class number (terminal node index), find
   the baseforms index number */
static int FindClassIndex(RegTransInfo *rt, int base) 
{
   int i;

   if (base == 0) 
      HError(7430, "FindClassIndex: Regression class is zero;\nNeed to run HHEd to build a regression class tree  before running adaptation!");

   for (i = 1; i <= rt->nBaseTransforms; i++)
      if (rt->baseforms[i]->baseformIndex == base)
         break;

   if (i > rt->nBaseTransforms)
      HError(7430, "FindClassIndex: Can't find base form class %d\n", base);
  
   return i;

}

/* Update the node occupation count */
static void GetNodeOcc(RegTree *t, RegTransInfo *rt) 
{
   short idx, i;
   RegNode *n;

   if (t != NULL) {
    
      GetNodeOcc(t->left, rt);

      n = t->nodeInfo;
      n->nodeOcc = 0.0;
      for (i = 1; i <= n->nBases; i++) {
         idx = n->bases[i];
         if (rt->baseforms[idx] == NULL)
            HError(7430, "GetNodeOcc:Can't find base class %d for node %d", i, 
                   n->nodeIndex);
         n->nodeOcc += rt->baseforms[idx]->occ;
      }
      if (trace & T_DET) {
         printf("Occupation count for node %d is %f\n",
                n->nodeIndex, n->nodeOcc);
         fflush(stdout);
      }
      GetNodeOcc(t->right, rt);

   }
}

/* Count of the number mixture components present at each node */
static void GetNodeComps(RegTree *t, RegTransInfo *rt) 
{
   short idx, i;
   RegNode *n;

   if (t != NULL) {
    
      GetNodeComps(t->left, rt);

      n = t->nodeInfo;
      n->nodeComps = 0;
      for (i = 1; i <= n->nBases; i++) {
         idx = n->bases[i];
         if (rt->baseforms[idx] == NULL)
            HError(7430, "GetNodeComps:Can't find base class %d for node %d", i, 
                   n->nodeIndex);
         n->nodeComps += rt->baseforms[idx]->nComponents;
      }
      if (trace & T_DET) {
         printf("Mixture component count for node %d is %d\n",
                n->nodeIndex, n->nodeComps);
         fflush(stdout);
      }
      GetNodeComps(t->right, rt);

   }
}

/* General function to check a node's occupation, to see
   if a transform is available or whether one needs constructing */
static Boolean CheckNodeOcc(RegTree *left, RegTree *right,
                            RegTransInfo *rt)
{
   if (left == NULL) {
      if (right->nodeInfo->nodeOcc < rt->nodeOccThresh || 
          right->nodeInfo->nodeComps <= rt->vSize)
         return TRUE;
      else
         return FALSE;
   }

   if (right == NULL) {
      if (left->nodeInfo->nodeOcc < rt->nodeOccThresh || 
          left->nodeInfo->nodeComps <= rt->vSize)
         return TRUE;
      else
         return FALSE;
   }
  
   if (left->nodeInfo->nodeOcc < rt->nodeOccThresh || 
       right->nodeInfo->nodeOcc < rt->nodeOccThresh || 
       left->nodeInfo->nodeComps <= rt->vSize ||
       right->nodeInfo->nodeComps <= rt->vSize)
      return TRUE;
   else
      return FALSE;

}

/* ----------------------------------------------------------------------*/
/*                 Auxilliary Calculation Functions                      */
/* ----------------------------------------------------------------------*/

/* note all auxilliary function calculation is used for tracing only */

/* Calculate the auxilliary function for a regression class */
static float CalcRegClassAuxilliary(BaseformAcc *base, float *frames) 
{
   int i, k, n, vSize;
   MixPDF *mp=NULL;
   RegAcc *ra=NULL;
   float sum = 0.0;
   Vector mean;
   Vector var;
   float gconst;
   float ivar;

   vSize = VectorSize(base->pdfList[1]->mean);
   n = base->nComponents;

   for (i = 1; i <= n; i++) {
      mp = base->pdfList[i];
      ra = GetRegAcc(mp);
      mean = mp->mean;
      gconst = mp->gConst;
      var  = mp->cov.var;
      if (ra->occ > MINOCC) {
         for(k = 1; k <= vSize; k++) {
            if (mp->ckind == INVDIAGC)
               ivar = var[k];
            else
               ivar = 1.0/var[k];
            sum += ra->obsSqSum[k] * ivar;
            sum -= 2.0 * ra->obsSum[k] * mean[k] * ivar;
            sum += ra->occ * mean[k] * mean[k] * ivar;
         }
         sum += ra->occ * gconst;
         *frames += ra->occ;
      }
   }
   return sum;
}
      

/* CalcAuxilliary: Calculate the auxilliary function */
static float CalcAuxilliary(RegTransInfo *rt) 
{
   float prob = 0.0;
   float frames = 0.0;
   int i;

   if (rt->transKind & MEANVAR) {
      for (i = 1; i <= rt->nBaseTransforms; i++)
         prob += CalcRegClassAuxilliary(rt->baseforms[i], &frames);
      printf("\n\tFinal auxilliary = %e, for %d frames\n",
             -0.5*prob/frames, (int) floor(frames + 0.5));
      fflush(stdout); 
  
   }
   else
      HError(-7420, "CalcAuxilliary: Only possible with mean & variance transforms\n");

   return prob;
}

/*------------------------------------------------------------------------*/
/*                    Input / Output functions                            */
/*------------------------------------------------------------------------*/

/* -------------------------- Lexical Scanner --------------------------- */
/* --------------------- Mirrors HModel Format -------------------------- */
#define MAXSYMLEN 40

/* Verbose keywords for readable external MLLR transform format */
static char *symMap[] = {
   "TRANSFORM", "MEAN_TR", "VARIANCE_TR", "BLOCK", "BIASOFFSET",
   "NODEOCC", "NODETHRESH", "NBLOCKS", "BASESTATS", 
   "LHSACC", "RHSACC", "VARACC", "EOF", ""
};

/* Internal and binary keyword representation */
typedef enum { 
   TRANSFORM, MEAN_TR, VARIANCE_TR, BLOCK, BIASOFFSET,
   NODEOCC, NODETHRESH, NBLOCKS, BASESTATS, 
   LHSACC, RHSACC, VARACC, EOFSYM, NULLSYM
} Symbol;

/* ---------------------------------------------------------------------- */

/* PutSymbol: output symbol to f in ascii or binary form */
static void PutSymbol(FILE *f, Symbol sym, Boolean binary)
{
   if (binary){
      fputc(':',f);fputc(sym,f);
   } else
      fprintf(f,"<%s>",symMap[sym]);
}

/* GetSymbol: read in a symbol in ascii or binary form */
static Boolean GetSymbol(Source *src, Symbol *sym)
{
   int i, c, imax;
   char buf[MAXSYMLEN];
   Symbol s;
  
   while (isspace(c=GetCh(src)));     /* Look for symbol */
   if (c != '<' && c != ':') {
      if (c == EOF) {
         *sym=EOFSYM; 
         return FALSE;
      }
      HError(7440, "GetSymbol: Symbol expected");
   }

   i=0; imax = MAXSYMLEN-1;
   if (c=='<') {                 /* Read verbose symbol string into buf */
      while ((c=GetCh(src)) != '>' && i<imax)
         buf[i++] = islower(c)?toupper(c):c;
      buf[i] = '\0';
      if (c != '>')
         HError(7440,"GetSymbol: > missing in symbol");
      for (s=TRANSFORM; s<=VARACC; s=(Symbol) (s+1)) /* Look symbol up in symMap */
         if (strcmp(symMap[s],buf) == 0) {
            *sym = s;
            return FALSE;                                /* and return */  
         }
   }
   else {
      /* Read binary symbol into buf */
      s = (Symbol) GetCh(src);
      if (s>=TRANSFORM && s<=VARACC) {
         *sym = s;
         return TRUE;                                /* and return */  
      }     
   }

   return FALSE;
}



static void SaveNodeTransform(FILE *fp, RegTransInfo *rt,
                              OffsetBMat *MTrans, Vector HTrans, short nodeId)
{
   short n, i;
   short bSize;
   short size;
   Matrix m;

   size = (short) rt->vSize;
   bSize = size / rt->nBlocks;

   /* print the node id */
   PutSymbol(fp, TRANSFORM, writeBin);
   WriteShort(fp, &nodeId, 1, writeBin);
   if (!writeBin) fputc('\n', fp);
   /* save the mean transformation */
   PutSymbol(fp, MEAN_TR, writeBin);
   WriteShort(fp, &bSize, 1, writeBin);
   if (!writeBin) fputc('\n', fp);

   if (bSize * rt->nBlocks != size)
      HError(7440, "Illegal number of blocks! Vector size %d, blocksize %d!",
             size, rt->nBlocks);

   /* save the mean matrix transform */
   for (n = 1; n <= rt->nBlocks; n++) {
      m = MTrans->a[n];
      PutSymbol(fp, BLOCK, writeBin); 
      WriteShort(fp, &n, 1, writeBin);
      if (!writeBin) fputc('\n', fp);
      for (i = 1; i <= bSize; i++)
         WriteVector(fp, m[i], writeBin);
   }
  
   /* save the mean bias vector */
   PutSymbol(fp, BIASOFFSET, writeBin);
   WriteShort(fp, &size, 1, writeBin);
   if (!writeBin) fputc('\n', fp);
   WriteVector(fp, MTrans->b, writeBin);

   /* save the covariance transform vector */
   if (rt->transKind & MEANVAR) {
      PutSymbol(fp, VARIANCE_TR, writeBin);
      WriteShort(fp, &size, 1, writeBin);
      if (!writeBin) fputc('\n', fp);
      WriteVector(fp, HTrans, writeBin);
   }

}

static void SaveTransforms(FILE *fp, RegTree *t, RegTransInfo *rt)
{
   RegNode *n;

   if (t != NULL) {

      SaveTransforms(fp, t->left,rt);

      if ( t->nodeInfo->nodeOcc > rt->nodeOccThresh && 
           t->nodeInfo->nodeComps > rt->vSize ) {
         n =  t->nodeInfo;
         if (t->left != NULL) {
            if (CheckNodeOcc(t->left, t->right, rt))
               SaveNodeTransform(fp, rt, n->WTrans, n->HTrans, n->nodeIndex);
         }
         else 
            SaveNodeTransform(fp, rt, n->WTrans, n->HTrans, n->nodeIndex);
      }

      SaveTransforms(fp, t->right, rt);

   }

}

static void SaveNodeOcc(FILE *fp, RegTree *t)
{
   RegNode *n;
  
   if (t != NULL) {

      SaveNodeOcc(fp, t->left);
      n =  t->nodeInfo;
      PutSymbol(fp, NODEOCC, writeBin);
      WriteShort(fp, &(n->nodeIndex), 1, writeBin);
      WriteFloat(fp, &(n->nodeOcc), 1, writeBin);
      if (!writeBin) fputc('\n', fp);
      SaveNodeOcc(fp, t->right);

   }
}


/* SaveTransformSetHeader: Save the header information */
void SaveTransformSetHeader(FILE *file, HMMSet *hset, RegTransInfo *rt,
                            char *uid, char *uname, char *chan, char *desc,
                            Boolean global, Boolean writeBin)
{
   MLink m;
   int n, h;
  
   /* get hold of the regression tree macro identifier (if it exists) */
   for (n=0,h=0; h<MACHASHSIZE; h++) {
      for (m=hset->mtab[h]; m!=NULL; m=m->next)
         if (m->type == 'r')
            break;
      if (m != NULL)
         break;
   }
  
   /* save the transform headers */
   /* user id */
   if (uid == NULL || strlen(uid) == 0)
      HError(7440, "SaveTransformSetHeader:Must specify a user id for the tmf header field <UID>!");
   fprintf(file, "<UID> %s\n", uid);

   /* save full name if available */
   if (uname == NULL || strlen(uname) == 0)
      fprintf(file, "<NAME> %s\n", undefined);
   else
      fprintf(file, "<NAME> %s\n", uname);

   /* mmf id */
   fprintf(file, "<MMFID> %s\n", hset->hmmSetId==NULL?"Unset":hset->hmmSetId);

   /* regression class type identifier */
   if (m == NULL || global)
      fprintf(file, "<RCID> global\n");
   else
      fprintf(file, "<RCID> %s\n", m->id->name);
  
   /* channel/microphone type */
   if (chan == NULL || strlen(chan) == 0)
      fprintf(file, "<CHAN> %s\n", undefined);
   else
      fprintf(file, "<CHAN> %s\n", chan);
 
   /* Store a general description field */
   if (desc == NULL || strlen(desc) == 0)
      fprintf(file, "<DESC> %s\n", undefined);
   else
      fprintf(file, "<DESC> %s\n", desc);

   PutSymbol(file, NBLOCKS, writeBin);
   WriteShort(file, &(rt->nBlocks), 1, writeBin);
   if (!writeBin) fputc('\n', file);
   /* fprintf(file, "<NBLOCKS> %d\n", rt->nBlocks); */
   PutSymbol(file, NODETHRESH, writeBin);
   WriteFloat(file, &rt->nodeOccThresh, 1, writeBin);
   if (!writeBin) fputc('\n', file);

}

/* EXPORT->SaveTransformSet: Save the transforms to a file */
void SaveTransformSet(HMMSet *hset, RegTransInfo *rt, char *saveFile, 
                      char *transPath, char *uid, char *uname, char *chan, 
                      char *desc, Boolean saveStats, Boolean global,
                      Boolean saveBinary)
{
   FILE *file;
   RegTree *rtree;
   char tfile[256];

   rtree = rt->rtree;
   writeBin |= saveBinary;

   /* create a transform file, made up of MMF name and a user's id */
   MakeFN(saveFile, transPath, NULL, tfile);
  
   if ((file=fopen(tfile,"w"))==NULL)
      HError(7440, "SaveTransformsSet: Cannot save transform file %s", tfile);

   if (trace & T_TOP) {
      printf("Saving transforms to %s\n", tfile);
   }

   SaveTransformSetHeader(file, hset, rt, uid, uname, chan, desc, global, 
                          writeBin);
   SaveNodeOcc(file, rtree);
   if (global)
      SaveNodeTransform(file, rt, rt->rtree->nodeInfo->WTrans, 
                        rt->rtree->nodeInfo->HTrans, 
                        rt->rtree->nodeInfo->nodeIndex);
   else {
      /* always save a root node transform, so see if it will be produced 
         by the  SaveTransforms function -- if not then call SaveNodeTransform */
      if (rt->rtree->left != NULL) {  
         if (rt->rtree->left->nodeInfo->nodeOcc > rt->nodeOccThresh && 
             rt->rtree->right->nodeInfo->nodeOcc > rt->nodeOccThresh && 
             rt->rtree->left->nodeInfo->nodeComps >= rt->vSize &&
             rt->rtree->right->nodeInfo->nodeComps >= rt->vSize)
            SaveNodeTransform(file, rt, rt->rtree->nodeInfo->WTrans, 
                              rt->rtree->nodeInfo->HTrans, 
                              rt->rtree->nodeInfo->nodeIndex);
      }
      SaveTransforms(file, rtree, rt);
   }

   if (saveStats)
      HError(-7440, "SaveTransformsSet: Saving adaptation statistics is not available!");
    
   fclose(file);

}


static void GetMeanTransformSet(Source *src,  RegNode *node, MemHeap *x, 
                                short vSize, short nBlocks,
                                RegTransType transKind, Boolean binary)

{
   int i, n;
   short bSize, bNum;
   BlockMatrix a;
   Matrix m = NULL;
   Symbol sym;
  
   if (!ReadShort(src, &bSize, 1, binary))
      HError(7440, "GetMeanTransformSet: expected mean vector block size");
  
   /* check that the transform has been created -- if not create! */
   if (node->WTrans == NULL)
      CreateNodeTransform(node, x, transKind, (int) vSize, (int) nBlocks);

   a = node->WTrans->a;

   for (n = 1; n <= nBlocks; n++) {
      binary = GetSymbol(src, &sym);
      if (sym != BLOCK)
         HError(7440, "GetMeanTransformSet: expected BLOCK symbol\n");
      if (!ReadShort(src, &bNum, 1, binary))
         HError(7440, "GetMeanTransformSet: expected block number");
      if (n == bNum)
         m = a[n];
      else
         HError(7440, "GetMeanTransformSet: Expected the block %d", n);
      for (i = 1; i <= bSize; i++) {
         if (!ReadVector(src, m[i], binary))
            HError(7440, "GetMeanTransformSet: Mean transform row %d expected", i+n*bSize);
      }
   }
}

static void GetBiasOffset(Source *src, RegNode *node, Boolean binary)
{
   short size;
   Vector v;

   if (!ReadShort(src, &size, 1, binary))
      HError(7440, "GetBiasOffset: expected bias offset size");

   v = node->WTrans->b;
   if (!ReadVector(src, v, binary))
      HError(7440, "GetBiasOffset: expected bias offset vector");

}

static void GetCovarTransformSet(Source *src, RegNode *node, Boolean useVar, 
                                 Boolean binary)
{
   int i;
   short size;
   Vector v, covar;
  
   if (!ReadShort(src, &size, 1, binary))
      HError(7440, "GetCovarTransformSet: expected covariance transform vector size");

   v = CreateVector(&gstack, size);

   covar = node->HTrans;
   if (!ReadVector(src, v, binary))
      HError(7440, "GetCovarTransformSet: expected covariance transform vector");
  
   if (useVar)
      for (i = 1; i <= size; i++)
         covar[i] = v[i];

   FreeVector(&gstack, v);
}


static Boolean LoadTransforms(Source *src, RegTransInfo *rt, Symbol *s)
{
   RegTree *t;
   RegNode *n;
   float tmp;
   Symbol sym;
   Boolean useVar=FALSE, varThere=FALSE;
   Boolean binary=FALSE;
   short nodeId, nblocks, bsize;

   /* load the number of blocks */
   binary = GetSymbol(src, &sym);
   if (sym != NBLOCKS) {
      HError(7440, "LoadTransforms: expected NBLOCKS symbol, not %s", symMap[sym]);
   }
   if (!ReadShort(src, &nblocks, 1, binary))
      HError(7440, "LoadTransforms: expected the number of blocks in block diagonal matrices");

   /* check blockSize */
   bsize = rt->vSize/nblocks;
   if (rt->vSize != bsize*nblocks)
      HError(7440, "LoadTransforms: Blocksize in tmf of %d incompatible with vector size %d\n", bsize, rt->vSize);


   if (nblocks != rt->nBlocks)
      HError(-7440, "LoadTransforms: Using no. of blocks specified in tmf (%d) which differs from the default, command-line or config setting!", nblocks);

   /* override any setting of number of blocks with number of blocks in tmf */
   rt->nBlocks = nblocks ;

   /* Load the occupation threshold and the node occupation counts */
   binary = GetSymbol(src, &sym);
   if (sym != NODETHRESH)
      HError(7440, "LoadTransforms: expected NODETHRESH symbol, not %s", symMap[sym]);

   /* if the occupation threshold is not set in the config file use saved one */
   if (strcmp(rt->transId->rcId, "global") != 0) {
      if (!ReadFloat(src, &rt->nodeOccThresh, 1, binary))
         HError(7440, "LoadTransforms: expected a node occ threshold value");
   }
   else
      if (!ReadFloat(src, &tmp, 1, binary))
         HError(7440, "LoadTransforms: expected a node occ threshold value");

   binary = GetSymbol(src, &sym);
   while (sym == NODEOCC) {
      if (!ReadShort(src, &nodeId, 1, binary))
         HError(7440, "LoadTransforms: expected a node index");
      t = FindNode(rt->rtree, NULL, nodeId);
      if (t==NULL)
         HError(7430, "LoadTransforms: Can't find node %d in tree\n", nodeId);
      if (!ReadFloat(src, &(t->nodeInfo->nodeOcc), 1, binary))
         HError(7440, "LoadTransforms: expected a node occ count for node %d",
                nodeId);
      binary = GetSymbol(src, &sym);
   }
  
   /* set some variables */
   if (rt->transKind & MEANVAR)
      useVar = TRUE;

   /* clear all transforms previously stored in memory */
   RemoveAllTransforms(rt->rtree, rt);

   /* Load the MLLR transforms */
   while(sym == TRANSFORM) {
      if (!ReadShort(src, &nodeId, 1, binary))
         HError(7440, "LoadTransforms: expected a node index");
      /* find node with index nodeId */
      t = FindNode(rt->rtree, NULL, nodeId);
      n = t->nodeInfo;
      sym = NULLSYM;
      while(sym != TRANSFORM && sym != BASESTATS && sym != EOFSYM) {
         binary = GetSymbol(src, &sym);
         switch(sym) {
         case MEAN_TR:
            GetMeanTransformSet(src, n, &rt->transMem, rt->vSize, 
                                rt->nBlocks, rt->transKind, binary);
            break;
         case BIASOFFSET:
            GetBiasOffset(src, n, binary);
            break;
         case VARIANCE_TR:
            GetCovarTransformSet(src, n, useVar, binary);
            varThere=TRUE;
            break;
         case TRANSFORM: case BASESTATS: case EOFSYM:
            break;
         default:
            HError(7440, "Unexpected symbol %s", symMap[sym]);
            break;
         }
      }
   }

   if (varThere == FALSE)
      rt->transKind = MEANONLY;
  
   *s = sym;
   return binary;

}

/* Load the transform set header */
static void LoadTransformHeader(Source *src, MemHeap *x, TransformId *tId,
                                char *hmmSetId, char *uid, MLink m)  
{
   char buf[256];

   /*
     Read the transform file header
   */

   /* Read the user id from the transform file and 
      check that it matches the current user */
   if (!ReadString(src, buf))
      HError(7440, "LoadTransformHeader: Expected marker <UID> in the TMF");
   if (strcmp(buf, "<UID>"))
      HError(7440, "LoadTransformHeader: Expected marker <UID> in the TMF");
   if (!ReadString(src, buf))
      HError(7440, "LoadTransformHeader: Expected a user identifier in the TMF");

   strcpy(tId->uid, buf);

   if (uid != NULL)
      if (strcmp(uid, tId->uid))
         HError(-7440, "LoadTransformHeader: Transform id does not match the user id\n\tTrans id %s, User id %s", tId->uid, uid);

   /* Read in the full user name if available (ie optional) */
   if (!ReadString(src, buf))
      HError(7440, "LoadTransformHeader: Expected marker <NAME> or <MMFID> in the TMF");
   if (strcmp(buf, "<NAME>") == 0) {
      if (!ReadString(src, buf))
         HError(7440, "LoadTransformHeader: Full user name is expected");
      if (strcmp(buf, "<MMFID>")) {
         strcpy(tId->name, buf);
         ReadLine(src, buf);
         strcat(tId->name, buf);
         SkipWhiteSpace(src);
         if (!ReadString(src, buf))
            HError(7440, "LoadTransformHeader: Expected marker <MMFID> in the TMF");
      }
      else {
         tId->name = NULL;
         HError(-7440, "LoadTransformHeader: User's full name unknown");
      } 
   }
   else {
      tId->name = NULL;
      HError(-7440, "LoadTransformHeader: Users' full name unknown");
   }

   /* read in the MMF identifier and check it */
   if (strcmp(buf, "<MMFID>"))
      HError(7440, "LoadTransformHeader: Expected marker <MMFID> in the TMF");
   if (!ReadString(src, buf))
      HError(7440, "LoadTransformHeader: Expected an MMF identifier in the TMF");
  
   /* if the hmmSet identifiers match then set the transform id hmm id */ 
   if (hmmSetId==NULL || strcmp(hmmSetId, buf))
      tId->hmmSetId = NULL;
   else
      tId->hmmSetId = hmmSetId;

   /* Read the regression class set identifier from the transform file */
   if (!ReadString(src, buf))
      HError(7440, "LoadTransformHeader: Expected marker <RCID> in the TMF");
   if (strcmp(buf, "<RCID>"))
      HError(7440, "LoadTransformHeader: Expected marker <RCID> in the TMF");  
   if (!ReadString(src, buf))
      HError(7440, "LoadTransformHeader: Expected a transform identifier in the TMF");
   strcpy(tId->rcId, buf);


   /* read in the microphone/channel description */
   if (!ReadString(src, buf))
      HError(7440, "LoadTransformHeader: Expected marker <CHAN> in the TMF");
   if (strcmp(buf, "<CHAN>"))
      HError(7440, "LoadTransformHeader: Expected marker <CHAN> in the TMF");
   /* skip the single white space */
   SkipWhiteSpace(src);
   if (!ReadLine(src, buf))
      HError(7440, "LoadTransformHeader: Expected a line describing channel type");
  
   strcpy(tId->chan, buf);
  
   /* Read in the general description */
   if (!ReadString(src, buf))
      HError(7440, "LoadTransformHeader: Expected marker <DESC> in the TMF");
   if (strcmp(buf, "<DESC>"))
      HError(7440, "LoadTransformHeader: Expected marker <DESC> in the TMF");
   /* skip the single white space */
   SkipWhiteSpace(src);
   if (!ReadLine(src, buf))
      HError(7440, "LoadTransformHeader: Expected a general description line");

   SkipWhiteSpace(src);
   strcpy(tId->desc, buf);
}


/* EXPORT->LoadTransformSet: Load the transforms into memory */
Boolean LoadTransformSet(HMMSet *hset, char *tfile, char *uid, 
                         RegTransInfo *rt, Boolean *useStats)  
{
   Source src;
   Symbol sym;
   Boolean binary;
   MLink m;
   int n, h;
  
   if(InitSource(tfile, &src, NoFilter)<SUCCESS)
      HError(7410,"LoadTransformSet: Can't open file %s", tfile);

   /* get hold of the regression tree macro identifier (if it exists) */
   for (n=0,h=0; h<MACHASHSIZE; h++) {
      for (m=hset->mtab[h]; m!=NULL; m=m->next)
         if (m->type == 'r')
            break;
      if (m != NULL)
         break;
   }
  
   /* Load the transform file header */
   LoadTransformHeader(&src, rt->hmem, rt->transId, hset->hmmSetId, uid, m);
   if (rt->transId->hmmSetId == NULL) {
      CloseSource(&src);
      return FALSE;
   }

   /* load the transforms and store in the tree */
   binary = LoadTransforms(&src,rt,&sym);

   if (*useStats)
      HError(-7410,"LoadTransformSet: Loading adaptation statistics is not available");

   if (trace & T_RIO) {
      PrintTree(rt->rtree, rt);
      fflush(stdout);
   }
 
   CloseSource(&src);
   return TRUE;
}

/*------------------------------------------------------------------------*/
/*         Initialisation of accumulate structures and hooks              */
/*------------------------------------------------------------------------*/

/* create a regression accumulation instance (used to link components
   in the same regression class together, and if adaptation to store
   regression statistics */
static RegAcc *CreateRegAccInstance(MemHeap *x, int transKind, int vSize)
{
  
   RegAcc *ra;
  
   ra = (RegAcc *) New(x, sizeof(RegAcc));

   ra->pre.time = -1;
   ra->pre.prob = LZERO;
   ra->obsSum   = NULL;
   ra->obsSqSum = NULL;
   ra->occ = 0.0;
  
   return ra;
  
}

/* allocate space in regression accumulator to store statistics */
static void CreateRegAccStorage(RegAcc *ra, RegTransInfo *rt) 
{

   ra->obsSum   =  CreateVector(&rt->pdfStatsMem,rt->vSize);
   ZeroVector(ra->obsSum);

   /* if using variances also */
   if (rt->transKind & MEANVAR) {
      ra->obsSqSum =  CreateVector(&rt->pdfStatsMem, rt->vSize);
      ZeroVector(ra->obsSqSum);
   }
}

/* AttachRegAccs: attach accumulators to hset */
static void AttachRegAccs(HMMSet *hset, MemHeap *x, RegTransType transKind, 
                          int vSize)
{
   HMMScanState hss;
   RegAcc *ra;
   int nregacc=0;
  
   NewHMMScan(hset,&hss);
   do {
      while (GoNextState(&hss,TRUE)) {
         while (GoNextStream(&hss,TRUE)) {            
            if (hss.isCont)                     /* PLAINHS or SHAREDHS */
               while (GoNextMix(&hss,TRUE)) {
                  nregacc += 1;
                  hss.mp->hook = (void *)
                     CreateRegAccInstance(x, transKind, vSize);
                  ra = GetRegAcc(hss.mp);
               }
            else
               HError(7450, "AttachRegAccs: Adaptation only available for PLAIN or SHARED systems!");
         }
      }
   } while (GoNextHMM(&hss));
   EndHMMScan(&hss);
  
   if (trace&T_NAC)
      printf("AttachRegAccs: Attached %d regression structures\n",
             nregacc);
}

/*------------------------------------------------------------------------*/
/*     Initialisation of the component and regression class structures    */
/*     and pointers for the adaptation and/or transformation process      */
/*------------------------------------------------------------------------*/

/* GroupRegMixes: For each class build a linked list of mixtures */
static void GroupRegMixes(HMMSet *hset, RegTransInfo *rt) {
  
   HMMScanState hss;
   IntVec baseMixes;
   short classNum;
   int i, j, index, n;
   char *s;

   /* sanity checking structure */
   baseMixes = CreateIntVec(&gstack, rt->nBaseTransforms);
   ZeroIntVec(baseMixes);

   /* allocate memory for the pdfList arrays and initialise */
   for (i = 1; i <= rt->nBaseTransforms; i++) {
      n = rt->baseforms[i]->nComponents;
      rt->baseforms[i]->pdfList = (MixPDF **) New(rt->hmem, 
                                                  sizeof(MixPDF *) * n);
      --(rt->baseforms[i]->pdfList);

      for (j = 1; j <= n; j++)
         rt->baseforms[i]->pdfList[j] = NULL;
      /* printf("%3d -- %3d has %4d\n", i, rt->baseforms[i]->baseformIndex, 
         rt->baseforms[i]->nComponents);
         fflush(stdout); */
   }

   /* assign mixpdfs to base class pdfLists */ 
   NewHMMScan(hset,&hss);
   do {
      while (GoNextState(&hss,TRUE)) {
         while (GoNextStream(&hss,TRUE)) {
            while (GoNextMix(&hss,TRUE)) {  
               switch(rt->classKind) {
               case ADPTTREE:
                  /* Group the components for the regression tree adaptation case */
                  classNum = hss.mp->rClass;
                  if (classNum > 0) {
                     index = FindClassIndex(rt, classNum);
                     baseMixes[index] += 1;
                     if (baseMixes[index] <= rt->baseforms[index]->nComponents)
                        rt->baseforms[index]->pdfList[baseMixes[index]] = hss.mp;
                     else
                        HError(7425, "GroupRegMixes: Mismatch num of comps for class %d -- %d, expected %d",
                               classNum, baseMixes[index], rt->baseforms[index]->nComponents);
                     s = HMMPhysName(hset, hss.hmm);
                     /* check for speechflag */
                     if (!strcmp(s, "sil") || !strcmp(s, "sp") 
                         || !strcmp(s, "h#") || !strcmp(s, "silist"))
                        rt->baseforms[index]->speechFlag = FALSE;
                     else
                        rt->baseforms[index]->speechFlag = TRUE;
                  }
                  break;
               case ADPTFIXED:
                  /* Group the components for the fixed classes adaptation case */
                  HError(7401, "GroupRegMixes: Fixed Group not yet available");
               case ADPTUNDEF:
                  HError(7419, "GroupRegMixes: Regression class type undefined!\n\tCheck config settings!");
               }
          
            }
         }
      }
   } while (GoNextHMM(&hss));
   EndHMMScan(&hss); 
  
   /* Do sanity check */
   if (rt->classKind == ADPTTREE || rt->classKind == ADPTFIXED) {
      for (i = 1; i <= rt->nBaseTransforms; i++)
         if (rt->baseforms[i]->nComponents != baseMixes[i])
            HError(7425, "GroupRegMixes: Mismatch between number of components %d and the expected number %d", 
                   baseMixes[i], rt->baseforms[i]->nComponents);
   }
   FreeIntVec(&gstack, baseMixes);
  
}

/* For this node, attach (via the bases member) the base classes */
static void AttachBaseformsToNode(RegNode *n, RegTree *t, RegTransInfo *rt) {
  
   short index, i;

   if (t != NULL) {
      AttachBaseformsToNode(n, t->left, rt);

      if (t->left == NULL) {
         for (i = 1; i <= n->nBases; i++)
            if (n->bases[i] == 0)
               break;
         /* find the index for this terminal node */
         index = FindClassIndex(rt, t->nodeInfo->nodeIndex);
         /* store a pointer to this base class for node n */
         n->bases[i] = index;
      }

      AttachBaseformsToNode(n, t->right, rt);
   }

}

/* for each node in the regression tree, attach the node's base classes */
static void AttachBaseformsToTree(RegTree *t, RegTransInfo *rt) {
  
   short i;
   RegNode *n;
   int nTerms=0;

   if (t != NULL) {
      AttachBaseformsToTree(t->left, rt);

      FindNumTerminals(t, &nTerms);
      n = t->nodeInfo;
      n->nBases = nTerms;
      n->bases = CreateShortVec(rt->hmem, n->nBases);

      for (i = 1; i <= n->nBases; i++)
         n->bases[i] = 0;
      AttachBaseformsToNode(n, t, rt);

      AttachBaseformsToTree(t->right, rt);
   }

}


static void GetBaseformComps(RegTree *t, RegTransInfo *rt) {
  
   int i;

   if (t != NULL) {
      GetBaseformComps(t->left, rt);

      if (t->left == NULL) {
         for (i = 1; i <= rt->nBaseTransforms; i++)
            if (rt->baseforms[i]->nComponents == 0) {
               rt->baseforms[i]->nComponents = t->nodeInfo->nodeComps;
               rt->baseforms[i]->baseformIndex = t->nodeInfo->nodeIndex;
               t->nodeInfo->nodeComps = 0;
               break;
            }
      }
      GetBaseformComps(t->right, rt);
   }

}

/* allocate further memory requirements to the tree */
static void AllocBaseMembers(RegTransInfo *rt)
{
   int i, nTerms = 0;
  
   FindNumTerminals(rt->rtree, &nTerms);
   rt->nBaseTransforms = nTerms;
   /* create storage for the transform matrices */
   rt->baseforms = (BaseformAcc **) New(rt->hmem, rt->nBaseTransforms * 
                                        sizeof(BaseformAcc *));
   --(rt->baseforms);

   for (i = 1; i <= rt->nBaseTransforms; i++)
      rt->baseforms[i] = CreateBaseform(rt);

   GetBaseformComps(rt->rtree, rt);

   AllocNodeTransforms(rt->rtree, rt);
  
}
/*------------------------------------------------------------------------*/
/*                    Main HAdapt initialisations                         */
/*------------------------------------------------------------------------*/
/* Initialise the regression classes (create memory storage etc) */
static void InitialiseRegClasses(HMMSet *hset, RegTransInfo *rt) 
{
   int bSize;
   RegTree *rtree=NULL;
   int h, n;
   MLink m;

   /* check blockSize */

   bSize = rt->vSize/rt->nBlocks;
   if (rt->vSize != bSize*rt->nBlocks)
      HError(7425, "InitialiseRegTrans: Blocksize of %d incompatible with vector size %d\n", bSize, rt->vSize);

   /* get hold of the regression tree from the macro structure holding it! */
   for (n=0,h=0; h<MACHASHSIZE; h++) {
      for (m=hset->mtab[h]; m!=NULL; m=m->next)
         if (m->type == 'r')
            break;
      if (m != NULL)
         break;
   }
   if (m != NULL)
      rtree = (RegTree *) m->structure;

   rt->hset = hset;

   if (hset->hmmSetId == NULL)
      HError(7421, "InitialiseRegClasses: MMF does not contain any identifier;\nUse HHEd to generate one\n");

   if (rtree == NULL)
      HError(7422, "InitialiseRegClasses: MMF does not contain a regression class tree;\nUse HHEd to generate one\n");

   /* only class kind available to set */
   rt->classKind = ADPTTREE;
   rt->rtree = rtree;

   /* allocate further memory requirements to the tree */
   AllocBaseMembers(rt);


   /* Now the regression classes are known, group the mixtures into linked
      lists for each regression class */
   GroupRegMixes(hset, rt);

   AttachBaseformsToTree(rtree, rt);
   if (trace & T_RIO) {
      PrintTree(rtree, rt);
   }

}

/* EXPORT->InitialiseTransform: Initialise transforms storage and grouping
   of hmmset components to regression base class using linked lists */
void InitialiseTransform( HMMSet *hset, MemHeap *x, RegTransInfo *rt,
                          Boolean adapt )
{  
   int vSize, bSize;
   TransformId *tId;
   HSetKind hsKind;
   int L;               /* number of logical HMM's */
   int P;               /* number of physical HMM's */

   /* set the block size and the regression class type */

   /* if block size is not set; use config variable/default */
   /* will be overriden by number of blocks in tmf if a tmf is read in */
   if (rt->nBlocks == 0)
      rt->nBlocks = blocks;

   /* if regression class type is undefined; use config variable/default */
   if (rt->classKind == ADPTUNDEF)
      rt->classKind = regClass;

   /* if regression transform type is undefined; use config variable/default */
   if (rt->transKind == TRANSUNDEF) {
      if (initVar)
         rt->transKind = MEANVAR;
      else
         rt->transKind = MEANONLY;
   }

   /* if the occupation threshold is not set; use config/default */
   if (rt->nodeOccThresh == 0.0)
      rt->nodeOccThresh = occThreshInit;

   rt->hmem = x;
   rt->transId = NULL;
   vSize = rt->vSize =(int)(hset->vecSize);
  
   if (rt->adptSil == TRI_UNDEF)
      rt->adptSil = (TriState) initAdptSil;

   /* set this for the time being */
   /* rt->adptSil = TRUE; */

   hsKind = hset->hsKind;
   P = hset->numPhyHMM;
   L = hset->numLogHMM;

   if (trace&T_TOP) {
      printf("System is ");
      switch (hsKind){
      case PLAINHS:  printf("PLAIN\n");  break;
      case SHAREDHS: printf("SHARED\n"); break;
      case TIEDHS: case DISCRETEHS:
         HError(7450, "HAdapt: Can only adapt PLAIN and SHARED systems");
      
         printf("%d Logical/%d Physical Models Loaded, VecSize=%d\n",
                L,P,vSize);
         if (hset->numFiles>0)
            printf("%d MMF input files\n",hset->numFiles);
      }

   }

   /* create an internal heap used for the transform structures */
   CreateHeap(&rt->transMem, "transMem", CHEAP, 1, 0.25, 0, 0);

   /* allocate space for SVD matrices and vector for bwd transforms */
   svdMat1 = CreateDMatrix(rt->hmem, vSize, vSize);
   u1      = CreateDMatrix(rt->hmem, vSize, vSize);
   v1      = CreateDMatrix(rt->hmem, vSize, vSize);
   w1      = CreateDVector(rt->hmem, vSize);

   /* allocate space for  blkZ & blkG */
   bSize = vSize/(rt->nBlocks);
   blkZ = CreateDMatrix(rt->hmem, bSize, bSize+1);
   blkG = CreateDMatrix(rt->hmem, bSize+1, bSize+1);
   blku = CreateDMatrix(rt->hmem, bSize+1, bSize+1);
   blkv = CreateDMatrix(rt->hmem, bSize+1, bSize+1);
   blkw = CreateDVector(rt->hmem, bSize+1);

   /* allocate space for the transform information identifier */
   tId = (TransformId *) New(x, sizeof(TransformId));

   tId->uid   = (char *) New(x, MAXSTRLEN);
   tId->name  = (char *) New(x, MAXSTRLEN);
   tId->rcId  = (char *) New(x, MAXSTRLEN);
   tId->chan   = (char *) New(x, MAXSTRLEN);
   tId->desc  = (char *) New(x, MAXSTRLEN);

   tId->uid[0] = 0;
   tId->name[0] = 0;
   tId->rcId[0] = 0;
   tId->chan[0] = 0;
   tId->desc[0] = 0;
   rt->transId = tId;

   /* attach hooks onto mix means for regression statistics accumulation */
   AttachRegAccs(hset, rt->hmem, rt->transKind, rt->vSize);

   /* Load in the regression classes group the mixtures based on classes */
   InitialiseRegClasses(hset, rt);

   if (trace & T_MEM) {
      printf("Basic transform initialisation done\n");
      PrintHeapStats(rt->hmem);
      fflush(stdout);
   }

   /* allocate space for the regression statistics at the base class level,
      if adaptation is to take place */
   /*if (adapt) 
     for (i = 1; i <= rt->nBaseTransforms; i++)
     CreateBaseformAccStorage(rt->baseforms[i], rt);*/

   /* Count of the number mixture components present at each node in tree */
   GetNodeComps(rt->rtree, rt);

}


/* EXPORT->InitialiseAdapt: Initialise adaptation storage */
void InitialiseAdapt(HMMSet *hset, MemHeap *x, RegTransInfo *rt )
{  
   int i, vSize;

   vSize = rt->vSize;

   /* create an internal heap used for accumulating statistics at the MixPDF level */
   CreateHeap(&rt->pdfStatsMem, "pdfMem", MHEAP, (vSize+1)*sizeof(float), 
              0.25, 500, 1000);

   /* allocate space for SVD matrices and vector */
   svdMat = CreateDMatrix(rt->hmem, vSize+1, vSize+1);
   u      = CreateDMatrix(rt->hmem, vSize+1, vSize+1);
   v      = CreateDMatrix(rt->hmem, vSize+1, vSize+1);
   w      = CreateDVector(rt->hmem, vSize+1);

   /* create temporary matrices to sum the regression level accumulates */
   Wg = CreateMatrix(rt->hmem, vSize, vSize);
   Gg = (Matrix *) New(rt->hmem, vSize * sizeof(Matrix));
   --Gg;
   for (i = 1; i <= vSize; i++) {
      Gg[i] = CreateMatrix(rt->hmem, vSize+1, vSize+1);
      ZeroMatrix(Gg[i]);
   }
   Zg = CreateMatrix(rt->hmem, vSize, vSize+1);
   Hg = CreateVector(rt->hmem, vSize);
   ZeroMatrix(Zg);
   ZeroVector(Hg);
   ZeroMatrix(Wg);

   /* allocate space for the regression statistics at the base class level
      if not already allocated */
   for (i = 1; i <= rt->nBaseTransforms; i++)
      if (rt->baseforms[i]->Z == NULL)
         CreateBaseformAccStorage(rt->baseforms[i], rt);

   if (trace & T_MEM) {
      printf("Transform + adaptation initialisation done\n");
      PrintHeapStats(rt->hmem);
      fflush(stdout);
   }
}


/* ----------------------------------------------------------------------*/
/*            Accumulate and Clear Adaptation Statistics                 */
/* ----------------------------------------------------------------------*/

/* EXPORT->AccAdaptFrame: Accumulate frame stats into specific mixture comp */
void AccAdaptFrame(double Lr, Vector speechVec, MixPDF *mp, RegTransInfo *rt) {

   RegAcc *ra;
   int k;
   int vSize;

   vSize = VectorSize(speechVec);

   /* get hold of the regression stats hook */
   ra = GetRegAcc(mp);
   if (ra != NULL) {
    
      if (ra->obsSum == NULL)
         CreateRegAccStorage(ra, rt); 

      ra->occ += Lr;

      for (k = 1; k <= vSize; k++) {
         ra->obsSum[k] += Lr * speechVec[k];
         if (rt->transKind & MEANVAR)
            ra->obsSqSum[k] += Lr * speechVec[k] * speechVec[k];
      }
   }
   else
      HError(7431, "AccAdaptFrame: Problem adding a frame of reg accumulation at the component level");
  
}

/* EXPORT->ClearRegCompStats: Clear for regression level accumulated stats */
void ClearRegCompStats(HMMSet *hset, RegTransInfo *rt) {

   HMMScanState hss;
   RegAcc *ra;
  

   /* zero the occupation counts */
   NewHMMScan(hset,&hss);
   do {
      while (GoNextState(&hss,TRUE)) {
         while (GoNextStream(&hss,TRUE)) {            
            if (hss.isCont)                     /* PLAINHS or SHAREDHS */
               while (GoNextMix(&hss,TRUE)) {
                  ra = GetRegAcc(hss.mp);
                  if (ra->obsSum != NULL) {
                     ra->occ = 0.0;
                     ra->obsSum = NULL;
                     ra->obsSqSum = NULL;
                  }
               }
            else
               HError(7450, "ClearRegCompStats: Only available for PLAIN or SHARED systems!");
         }
      }
   } while (GoNextHMM(&hss));
   EndHMMScan(&hss);

   /* reset the heap holding all the summed scaled speech vectors */
   ResetHeap(&rt->pdfStatsMem);

}

/* EXPORT->ClearBaseClassStats: Clear for base class level accumulated stats */
void ClearBaseClassStats(RegTransInfo *rt) {

   int i, j, vSize;
   BaseformAcc *b;

   vSize = rt->vSize;

   /* clear the base class stats */
   for (i = 1; i <= rt->nBaseTransforms; i++) {

      b = rt->baseforms[i];
      for (j = 1; j <= vSize; j++) {
         ZeroBlockTriMat(b->G[j]->a);
         ZeroVector(b->G[j]->b);
      }
      ZeroBlockMat(b->Z->a);
      ZeroVector(b->Z->b);
      ZeroVector(b->H);

      /* also need to clear the base node occupations for the dynamic tree */
      b->occ = 0.0;
   }

}

/* ----------------------------------------------------------------------*/
/*            Mean and Covariance Transformation Functions               */
/* ----------------------------------------------------------------------*/


/* ---------------------- TRANSFORM MEANS -------------------------------*/
 
/* Apply the mean tranformation to all components in this regression class */
static void ApplyMeanClassTransform(HMMSet *hset, BaseformAcc *base, 
                                    OffsetBMat *m, int nodeId) {

   int i, j, n, vSize;
   MixPDF *mp=NULL;
   int nTransformed = 0;
   SVector mean;
   BlockMatrix a;
   Vector b;

   vSize = VectorSize(base->pdfList[1]->mean);
   n = base->nComponents;

   /* apply mean transform to all means in base class that have not been seen */
   for (j = 1; j <= n; j++) {
      mp = base->pdfList[j];
      mean = mp->mean;
      if (!IsSeenV(mean)) {
         a = m->a;
         b = m->b;
         MultBlockMat_Vec(a, mean, mean);
         for (i = 1; i <= vSize; i++)
            mean[i] += b[i];
         /* mark as seen */
         TouchV(mean);
         nTransformed += 1;
      }
   }

   /* unmark all the means in this base class */
   for (i = j; j <= n; j++) {
      mp = base->pdfList[j];
      mean = mp->mean;
      if (IsSeenV(mean))
         UntouchV(mean);
   }

   if (trace&T_TRA) {
      printf("There were %d transformations for regression node class %d\n",
             nTransformed, nodeId);
      fflush(stdout);
   }

}

/* apply mean tranform at node n to the regression sub-tree t */
static void ApplyMeanTrans(RegTree *t, RegNode *n, RegTransInfo *rt) {

   short idx;

   if (t != NULL) {
      ApplyMeanTrans(t->left, n, rt);

      if (t->left == NULL) {
         /* terminal node -- corresponds to base class */
         idx = t->nodeInfo->bases[1];
         if (rt->adptSil || rt->baseforms[idx]->speechFlag) {
            if (trace&T_USE) {
               printf("Applying transform at nd %d to components at nd %d (%f)\n",
                      n->nodeIndex, t->nodeInfo->nodeIndex,  t->nodeInfo->nodeOcc);
               fflush(stdout);
            }
            ApplyMeanClassTransform(rt->hset, rt->baseforms[idx], n->WTrans,
                                    t->nodeInfo->nodeIndex);
         }
      }
      ApplyMeanTrans(t->right, n, rt);
   }
}

/* EXPORT->ApplyMeanTransforms: 
   Apply the mean transform to regression classes */
void ApplyMeanTransforms(RegTransInfo *rt, RegTree *t)
{
   short idx;
   Boolean adaptSil;
   float nodeOccThresh;
   int vSize; 
   RegNode *n;

   vSize = rt->vSize;
   nodeOccThresh = rt->nodeOccThresh;
   adaptSil = (Boolean) rt->adptSil;

   /* apply the transform, given the occupation count at each node
      if the child occupation count is above the the threshold then
      move down the tree, otherwise apply the current node's */
   if (t != NULL) {

      ApplyMeanTransforms(rt, t->left);

      if (t->nodeInfo->nodeOcc > nodeOccThresh &&
          t->nodeInfo->nodeComps > vSize) {

         n =  t->nodeInfo;
         if (t->left != NULL) {
            if (CheckNodeOcc(t->left, NULL, rt))
               ApplyMeanTrans(t->left, n, rt);
            if (CheckNodeOcc(NULL, t->right, rt))
               ApplyMeanTrans(t->right, n, rt);
        
         }
         else { 
            /* terminal node -- corresponds to base class */
            idx = t->nodeInfo->bases[1];
            if (adaptSil || rt->baseforms[idx]->speechFlag) {
          
               ApplyMeanClassTransform(rt->hset, rt->baseforms[idx], n->WTrans,
                                       n->nodeIndex);
               if (trace&T_USE) {
                  printf("Applying transform at nd %d to comps at nd %d (%f)\n",
                         n->nodeIndex, t->nodeInfo->nodeIndex, 
                         t->nodeInfo->nodeOcc);
                  fflush(stdout);
               }
            }
         }
      }
     
      ApplyMeanTransforms(rt, t->right);
  
   }
}


/* ---------------------- TRANSFORM COVARIANCES -------------------------*/

/* apply the covariance tranformation to a regression class */
static void ApplyCovarClassTransform(HMMSet *hset, BaseformAcc *base, 
                                     Vector H, int nodeId) 
{

   int i, j, n, vSize;
   MixPDF *mp=NULL;
   int nTransformed = 0;
   float hi;
   SVector var;

   vSize = VectorSize(base->pdfList[1]->cov.var);
   n = base->nComponents;

   /* apply covar transform to all covarss in base class that have not been seen */
   for (j = 1; j <= n; j++) {
      mp = base->pdfList[j];
      var = mp->cov.var;
      if (!IsSeenV(var)) {
         for (i = 1; i <= vSize; i++) {
            if (mp->ckind == INVDIAGC) {
               hi = 1.0 / H[i];
            }
            else {
               hi = H[i];
            }
            var[i] *= hi;
         }      
         /* mark as seen */
         TouchV(var);
         nTransformed += 1;
    
      }
   }

   /* unmark all the covars in this base class */
   for (j = 1; j <= n; j++) {
      mp = base->pdfList[j];
      var = mp->cov.var;
      if (IsSeenV(var))
         UntouchV(var);
      /* fix the gConst too */
      if (mp->ckind == INVDIAGC)
         FixInvDiagGConst(mp);
      else
         FixDiagGConst(mp);
   }

   if (trace&T_TRA) {
      printf("There were %d (%d) covar transformations for reg base class %d\n",
             nTransformed, n, nodeId);
      fflush(stdout);
   }

}

/* apply covariance tranform at node n to the regression sub-tree t */
static void ApplyCovarTrans(RegTree *t, RegNode *n,  RegTransInfo *rt) {
  
   short idx;

   if (t != NULL) {
      ApplyCovarTrans(t->left, n, rt);

      if (t->left == NULL) {
         /* terminal node -- corresponds to base class */
         idx = t->nodeInfo->bases[1];
         if (rt->adptSil || rt->baseforms[idx]->speechFlag) {
            if (trace&T_USE) {
               printf("Applying covartransform at node %d to components at node %d (%f)\n",
                      n->nodeIndex, t->nodeInfo->nodeIndex,  t->nodeInfo->nodeOcc);
               fflush(stdout);
            }
            ApplyCovarClassTransform(rt->hset, rt->baseforms[idx], n->HTrans,
                                     t->nodeInfo->nodeIndex);
         }
      }
    
      ApplyCovarTrans(t->right, n, rt);
   }
}

/* EXPORT->ApplyCovarTransforms: 
   Apply the covariance transform to the regression classes */
void ApplyCovarTransforms(RegTransInfo *rt, RegTree *t)
{
   short idx;
   Boolean adaptSil;
   int vSize; 
   float nodeOccThresh;
   RegNode *n;

   vSize = rt->vSize;
   nodeOccThresh = rt->nodeOccThresh;
   adaptSil = (Boolean) rt->adptSil;

   /* apply the transform, given the occupation count at each node
      if the child occupation count is above the the threshold then
      move down the tree, otherwise apply the current node's */
   if (t != NULL) {

      ApplyCovarTransforms(rt, t->left);

      if ( t->nodeInfo->nodeOcc > nodeOccThresh &&
           t->nodeInfo->nodeComps > vSize ) {
         n =  t->nodeInfo;
         if (t->left != NULL) {
            if (CheckNodeOcc(t->left, NULL, rt))
               ApplyCovarTrans(t->left, n, rt);
            if (CheckNodeOcc(NULL, t->right, rt))
               ApplyCovarTrans(t->right, n, rt);
         }
         else { 
            /* terminal node -- corresponds to base class */
            idx = t->nodeInfo->bases[1];
            if (adaptSil || rt->baseforms[idx]->speechFlag)
               ApplyCovarClassTransform(rt->hset, rt->baseforms[idx], n->HTrans, 
                                        n->nodeIndex);
         }
      }

 
      ApplyCovarTransforms(rt, t->right);
   }
}

/* EXPORT->ApplyMeanGlobalTransform: Apply the mean Global transformation 
   -- returns true if global is found */
Boolean ApplyMeanGlobalTransform(RegTransInfo *rt) {

   short i, idx;

   /* check to see if root transform is available */
   if (rt->rtree->nodeInfo->WTrans == NULL)
      return FALSE;

   /* apply mean transform at root node only */
   if (rt->rtree->nodeInfo->nodeOcc > rt->nodeOccThresh)
      for (i = 1; i <= rt->rtree->nodeInfo->nBases; i++) {
         idx = rt->rtree->nodeInfo->bases[i];
         ApplyMeanClassTransform(rt->hset, 
                                 rt->baseforms[idx],
                                 rt->rtree->nodeInfo->WTrans, 
                                 rt->rtree->nodeInfo->nodeIndex);
      }

   return TRUE;
}

/* EXPORT->ApplyCovarGlobalTransform: Apply the covar Global transformation 
   -- returns true if global is found */
Boolean ApplyCovarGlobalTransform(RegTransInfo *rt) {

   short i, idx;

   /* check to see if root transform is available */
   if (rt->rtree->nodeInfo->WTrans == NULL)
      return FALSE;

   /* apply covar transform at root node only */
   if (rt->transKind & MEANVAR)
      if (rt->rtree->nodeInfo->nodeOcc > rt->nodeOccThresh)
         for (i = 1; i <= rt->rtree->nodeInfo->nBases; i++) {
            idx = rt->rtree->nodeInfo->bases[i];
            ApplyCovarClassTransform(rt->hset, 
                                     rt->baseforms[idx],
                                     rt->rtree->nodeInfo->HTrans, 
                                     rt->rtree->nodeInfo->nodeIndex);
         }

   return TRUE;
}

/* EXPORT->ApplyTransforms: Apply the mean and possibly 
   the covar MLLR transformations */
void ApplyTransforms(RegTransInfo *rt) {

   Boolean global=FALSE;

   /* check to see if transform is a global one */
   if (!strcmp(rt->transId->rcId, "global"))
      global = TRUE;

   if (global)
      ApplyMeanGlobalTransform(rt);
   else
      ApplyMeanTransforms(rt, rt->rtree);
  
   if (rt->transKind & MEANVAR) {
      if (global)
         ApplyCovarGlobalTransform(rt);
      else
         ApplyCovarTransforms(rt, rt->rtree);
   }
   
   if (trace & T_TOP) {
      printf("Applied mean ");
      if (rt->transKind & MEANVAR)
         printf("and variance ");
      printf("transforms\n");
      fflush(stdout);
   }
}

/* ----------------------------------------------------------------------*/
/*     Backward Transformation Functions (so only 1 HMMset needed)       */
/* ----------------------------------------------------------------------*/

/* Get the backward transform for a node in the regression tree */
static void GetBackwardTransform(RegNode *node, int vSize, int nBlocks) {

   int i, j, n, fi, fj;
   int bStart, bSize;
   Matrix A;
   Matrix m;

   A = node->backTrans;
   ZeroMatrix(A);
   ZeroDMatrix(svdMat1);

   bStart = 0;
   bSize = vSize / nBlocks;

   /* convert block matrix to a double matrix so it can be inverted */
   for(n = 1; n <= nBlocks; n++) {
      m = node->WTrans->a[n];
      for(i = 1; i<= bSize; i++) {
         fi = i + bStart;
         for(j = 1; j <= bSize; j++) {
            fj = j + bStart;
            svdMat1[fi][fj] = (double) m[i][j];
         }
      }
      bStart += bSize;
   }
  
   /* invert matrix */
   InvSVD(svdMat1, u1, w1, v1, svdMat1);

   /* convert double to float */
   DMat2Mat(svdMat1, A);

}

/* Storing only one set of gaussians -- transform means and covar back
   to their original SI model form */
static void TransformMixesBack(BaseformAcc *base, HMMSet *hset, RegNode *n, 
                               int transKind, int vSize, int i) {
  
   MixPDF *mp=NULL;
   SVector mean, var;
   int j, k, m, nComps;
  
   nComps = base->nComponents;

   for (m = 1; m <= nComps; m++) {
      mp = base->pdfList[m];
      ZeroDVector(w1);
      mean     = mp->mean;
      var      = mp->cov.var;

      /* do mean only */
      if (!IsSeenV(mean)) {
         for (j = 1; j <= vSize; j++)
            mean[j] -= n->WTrans->b[j];

         for (j = 1; j <= vSize; j++)
            for (k = 1; k <= vSize; k++)
               w1[j] += n->backTrans[j][k] * mean[k];
         for (j = 1; j <= vSize; j++) {
            mean[j] = (float) w1[j];
         }
         TouchV(mean);
      }

      /* do var only */
      if (transKind & MEANVAR) {
         if (!IsSeenV(var)) {
            for (j = 1; j <= vSize; j++) {
               if (mp->ckind == INVDIAGC)
                  var[j] *= n->HTrans[j];
               else
                  var[j] /= n->HTrans[j];
            }

            TouchV(var);
         }
      }
   }

   /* unmark means and vars for this base class */
   for (m = 1; m <= nComps; m++) {
      mp = base->pdfList[m];
      mean     = mp->mean;
      var      = mp->cov.var;  
      if (IsSeenV(mean))
         UntouchV(mean);
      if (IsSeenV(var))
         UntouchV(var);
      /* fix the gConst too */
      if (mp->ckind == INVDIAGC)
         FixInvDiagGConst(mp);
      else
         FixDiagGConst(mp);
   }

}

/* Apply all the backward transformation at node n to the rest of tree t */
static void ApplyBwdTransform(RegTree *t, RegNode *n, RegTransInfo *rt) {

   short idx;

   if (t != NULL) {

      ApplyBwdTransform(t->left, n, rt);

      if (t->left == NULL) {
         idx = t->nodeInfo->bases[1];
         /* this is a terminal node! */
         if (rt->adptSil || rt->baseforms[idx]->speechFlag) {
            TransformMixesBack(rt->baseforms[idx], rt->hset, n, 
                               rt->transKind, rt->vSize, 
                               t->nodeInfo->nodeIndex);
         }
      }

      ApplyBwdTransform(t->right, n, rt);
   }

}

/* Apply all the backward transformations for the regression tree */
static void ApplyBwdTransforms(RegTree *t, RegTransInfo *rt) {

   short idx;
   int vSize;
   Boolean backTrans = FALSE;
  
   vSize = rt->vSize;

   if (t != NULL) {

      ApplyBwdTransforms(t->left, rt);

      /* apply back transform at node */
      if (t->nodeInfo->nodeOcc > rt->nodeOccThresh &&
          t->nodeInfo->nodeComps > vSize) {
         if (t->left != NULL) {
            if (CheckNodeOcc(t->left, NULL, rt)) {
               t->nodeInfo->backTrans = CreateMatrix(rt->hmem, vSize, vSize);
               GetBackwardTransform(t->nodeInfo, vSize, rt->nBlocks);
               ApplyBwdTransform(t->left, t->nodeInfo, rt);
               backTrans = TRUE;
            }
            if (CheckNodeOcc(NULL, t->right, rt)) {
               if (!backTrans) {
                  t->nodeInfo->backTrans = CreateMatrix(rt->hmem, vSize, vSize);
                  GetBackwardTransform(t->nodeInfo, vSize, rt->nBlocks);
                  backTrans=TRUE;
               }
               ApplyBwdTransform(t->right, t->nodeInfo, rt);
            }
            if (backTrans) {
               Dispose(rt->hmem, t->nodeInfo->backTrans);
               backTrans=FALSE;
            }
         }
         else {
            /* this is a terminal node! */
            idx = t->nodeInfo->bases[1];
            t->nodeInfo->backTrans = CreateMatrix(rt->hmem, vSize, vSize);
            ZeroMatrix(t->nodeInfo->backTrans);
            GetBackwardTransform(t->nodeInfo, vSize, rt->nBlocks);
            TransformMixesBack(rt->baseforms[idx], rt->hset, 
                               t->nodeInfo, rt->transKind,
                               vSize, t->nodeInfo->nodeIndex);
            Dispose(rt->hmem, t->nodeInfo->backTrans);
         }
      }

      ApplyBwdTransforms(t->right, rt);
   }
  
}


/* EXPORT->ApplyBackwardGlobalTransform:
   Apply the backward globale transformation to the model set 
   to obtain the "original" models -- returns true if global is found */
Boolean ApplyBackwardGlobalTransform(RegTransInfo *rt) {

   short i, idx;

   if (rt->rtree->nodeInfo->WTrans == NULL)
      return FALSE;

   if (rt->rtree->nodeInfo->nodeOcc > rt->nodeOccThresh && 
       rt->rtree->nodeInfo->WTrans != NULL) {
      rt->rtree->nodeInfo->backTrans = CreateMatrix(rt->hmem, 
                                                    rt->vSize,
                                                    rt->vSize);
      ZeroMatrix(rt->rtree->nodeInfo->backTrans);
      GetBackwardTransform(rt->rtree->nodeInfo, rt->vSize, rt->nBlocks);
      for (i = 1; i <= rt->rtree->nodeInfo->nBases; i++) {
         idx = rt->rtree->nodeInfo->bases[i];
         TransformMixesBack(rt->baseforms[idx], rt->hset, 
                            rt->rtree->nodeInfo, rt->transKind,
                            rt->vSize, rt->rtree->nodeInfo->nodeIndex);
      }
      Dispose(rt->hmem, rt->rtree->nodeInfo->backTrans);
   }
  
   return TRUE;

}


/* EXPORT->AppplyBackwardTransforms: Interface to ApplyBwdTransforms
   Apply the backward transformation to the model set to obtain the "original" models */
void ApplyBackwardTransforms(RegTransInfo *rt) {

   ApplyBwdTransforms(rt->rtree, rt);

}

/* ----------------------------------------------------------------------*/
/*            Regression Class Accumulation Functions                    */
/* ----------------------------------------------------------------------*/

/* accumulate the LHS of the regression equation into Z */
static void AccTransLHS(OffsetBMat *Z, Vector mean, Vector covar, 
                        CovKind ckind, RegTransInfo *rt, Vector obs)
{
   int i, n, k, fi, fk, bStart, bSize, vSize, nBlocks;
   float scaledObs, ivar;
   Vector b;
   Matrix a;

   vSize   = rt->vSize;
   nBlocks = rt->nBlocks;
   bSize   = vSize/nBlocks;
   b       = Z->b;

   bStart = 0;
   for (n = 1; n <= nBlocks; n++) {
      a = Z->a[n];
      for (i = 1; i <= bSize; i++) {
         fi = bStart+i;
         if (ckind == INVDIAGC)
            ivar=covar[fi];
         else
            ivar=1.0/covar[fi];
         scaledObs = obs[fi] * ivar;
         /* do bias */
         b[fi] +=  scaledObs;
         for (k = 1; k <= bSize; k++) {
            fk = bStart+k;
            a[i][k] += scaledObs * mean[fk];
         }
      }
      bStart += bSize;
   }

}

/* accumulate the RHS of the regression equation into G */
static void AccTransRHS(OffsetTriBMat **G, Vector mean, Vector covar,
                        CovKind ckind, RegTransInfo *rt, float occ) 
{
   int i, j, fj, q, n, bStart, bSize, vSize, nBlocks;
   float scaledMean, scaledVar;
   TriMat A;
   Vector b;
   Vector Aj;

   /* do initial assignments */
   vSize   = rt->vSize;
   nBlocks = rt->nBlocks;
   bSize   = vSize/nBlocks;
  
   for (i = 1; i <= vSize; i++) {
      if (ckind == INVDIAGC)
         scaledVar=covar[i]*occ;
      else
         scaledVar=occ/covar[i];
      bStart = 0;
      b = G[i]->b;
      /* do first bias term */
      b[1] += scaledVar;
      for (n = 1; n <= nBlocks; n++) {
         A = G[i]->a[n];
         for (j = 1; j <= bSize; j++) {
            fj = bStart+j;
            scaledMean = scaledVar * mean[fj];
            Aj = A[j];
            /* do bias */
            b[fj+1] += scaledMean;
            for (q = 1; q <= j; q++)
               Aj[q] += scaledMean * mean[bStart+q];
         }
         bStart += bSize;
      }
   }
}

/* Accumulate the covariance (numerator) into vector Z */ 
static void AccCovarTrans(MixPDF *mp, Vector Z, Vector mean, 
                          Vector var, Vector obsSum, Vector obsSqSum)
{
   RegAcc *ra;
   float ivar, meank, occ;
   double sum;
   int k, vSize;

  
   ra = GetRegAcc(mp);
   occ      = ra->occ;

   vSize = VectorSize(mean);
   for (k = 1; k <= vSize; k++) {
      sum = 0.0;
      meank = mean[k];
      if (mp->ckind == INVDIAGC)
         ivar = var[k];
      else
         ivar = 1.0 / var[k];
    
      sum  = obsSqSum[k] * ivar;
      sum -= 2.0 * meank * obsSum[k] * ivar ;
      sum += meank * meank * occ * ivar ;
    
      Z[k] += sum ;
   }
}

/* ----------------------------------------------------------------------*/
/*              Mean Transformation Calculation Functions                */
/* ----------------------------------------------------------------------*/

static void CalcAdaptXFromAllZAllG(RegTransInfo *rt, Matrix* AllG, 
                                   Matrix AllZ, OffsetBMat* adaptX)
{
   int vSize, nBlocks, bSize, q, fullq, cntBlk, bStart, fullk, i, k, row, col, fullrow, fullcol;
   Matrix allGq, xan;
   BlockMatrix xa;
   Vector xb;
   double dbltmp;

   vSize = rt->vSize;
   nBlocks = rt->nBlocks;
   bSize = vSize / nBlocks;

   xa = adaptX->a;
   xb = adaptX->b;

   bStart=0;
   for(cntBlk = 1; cntBlk <= nBlocks; cntBlk++) {
      xan = xa[cntBlk];
      for(q=1; q<=bSize; q++) {
         /* fetch the sub Z */
         fullq = bStart + q;
         blkZ[q][1] = AllZ[fullq][1];
         for(k=1; k<=bSize; k++) {
            fullk = bStart + k;
            blkZ[q][k+1] = AllZ[fullq][fullk+1];
         }
         /* fetch the sub G */
         allGq = AllG[fullq];
         blkG[1][1] = allGq[1][1];
         for(row=1; row<=bSize; row++) {
            fullrow = bStart+row;
            blkG[row+1][1]= blkG[1][row+1] = allGq[fullrow+1][1];
            for(col=1; col<=bSize; col++) {
               fullcol = bStart + col;
               blkG[row+1][col+1]= allGq[fullrow+1][fullcol+1];
            }
         }
         /* inverse the sub G */
         InvSVD(blkG, blku, blkw, blkv, blkG);
         /* calc the current row of X */
         fullq = bStart + q;
         dbltmp = 0.0;
         for(k=1; k<=bSize+1; k++) {
            dbltmp += blkZ[q][k]*blkG[k][1];
         }
         xb[fullq] = dbltmp;
         for(i=1; i<=bSize; i++) {
            dbltmp = 0.0;
            for(k=1; k<=bSize+1; k++) {
               dbltmp += blkZ[q][k]*blkG[k][i+1];
            }
            xan[q][i] = dbltmp;
         }
         
      }
      bStart += bSize;
      
   }
}

/* Calculate the G and Z regression accumulates for each regression class */
static float GetMeanClassAccumulates(BaseformAcc *base, RegTransInfo *rt)
{
   int i, n;
   MixPDF *mp=NULL;
   RegAcc *ra;
   float sum = 0.0;
   HMMSet *hset;

   hset = rt->hset;
   n = base->nComponents;

   for (i = 1; i <= n; i++) {
      mp = base->pdfList[i];
      ra = GetRegAcc(mp);
      if (ra->occ > MINOCC) {
         AccTransLHS(base->Z, mp->mean, mp->cov.var, mp->ckind, rt, ra->obsSum);
         AccTransRHS(base->G, mp->mean, mp->cov.var, mp->ckind, rt, ra->occ);
      }
      sum += ra->occ;
   }

   return sum;
}

/* Calculate the mean transform for a regression tree node */
static void GetMeanTransform(RegNode *n, MemHeap *x, RegTransInfo *rt)
{
   short idx;
   int i, j, k, nBlocks, vSize;
   OffsetBMat *m;
   BlockMatrix a;
   Vector b;
   Matrix Gk, Gi;
   OffsetTriBMat *Gti;
   Boolean newTransform=FALSE;

   vSize = rt->vSize;
   nBlocks = rt->nBlocks;

   /* check that the transform has been created -- if not create! */
   if (n->WTrans == NULL) {
      CreateNodeTransform(n, x, rt->transKind, vSize, nBlocks);
      newTransform = TRUE;
   }
   m = n->WTrans;
   a = m->a;
   b = m->b;
  
   /* clear the transform matrix and the bias offset vector */
   for (k = 1; k<=vSize; k++)
      ZeroMatrix(Gg[k]);
   ZeroDMatrix(svdMat);
   ZeroMatrix(Wg);
   ZeroMatrix(Zg);
   ZeroVector(b);

   if (n->nBases < 1)
      HError(7430, "GetMeanTransform: No base class for node %d",
             n->nodeIndex);

   /* Sum the regression accumulates for the base classes 
      belonging to this node */
   for (i = 1; i <= n->nBases; i++) {
      idx = n->bases[i];
      if (rt->baseforms[idx] == NULL)
         HError(7430, "GetMeanTransform:Can't find base class %d for node %d", i, 
                n->nodeIndex);
      for (k = 1; k<=vSize; k++) {
         Gk = Gg[k];
         Gti = rt->baseforms[idx]->G[k];
         AddGSymMatrix(Gti, Gk);
      }
      AddZMatrix(rt->baseforms[idx]->Z, Zg);
   }
    
   /* Now do the transform calculation */

   /* G's accumulated stats only in lower triangle, 
      so fill in upper triangles */

   for (i = 1; i <= vSize; i++) {
      Gi = Gg[i];
      for (j = 1; j<=vSize+1;j++)
         for (k = 1; k<=j; k++)
            Gi[k][j] = Gi[j][k];
   }
   CalcAdaptXFromAllZAllG(rt, Gg, Zg, m);
}

/* Calculate the mean transform for each regression class */
static void GetMeanClassTransforms(RegTree *t, RegTransInfo *rt) 
{
   int vSize;

   vSize = rt->vSize;

   if (t != NULL) {
      GetMeanClassTransforms(t->left, rt);
      if (t->nodeInfo->nodeOcc > rt->nodeOccThresh &&
          t->nodeInfo->nodeComps > vSize ) {
         if (t->left != NULL) {
            if (CheckNodeOcc(t->left, t->right, rt))
               GetMeanTransform(t->nodeInfo, &rt->transMem, rt);
         }
         else { 
            GetMeanTransform(t->nodeInfo, &rt->transMem, rt);
         }
      }
      GetMeanClassTransforms(t->right, rt);
   }

}



/* EXPORT CalcMeanTransforms: Calculate the mean transforms for every 
   regression class */
void CalcMeanTransforms(RegTransInfo *rt) {

   int i;

   /* reset all the transforms */
   RemoveAllTransforms(rt->rtree, rt);

   for (i = 1; i <= rt->nBaseTransforms; i++) {
      if (trace & T_DET) {
         printf("(Means) Accumulating at the regression level for base %d\n",
                rt->baseforms[i]->baseformIndex);
         fflush(stdout);
      }

      rt->baseforms[i]->occ +=  
         GetMeanClassAccumulates(rt->baseforms[i], rt);
   }
  
   GetNodeOcc(rt->rtree, rt);
   GetMeanClassTransforms(rt->rtree, rt);

   /* always make sure to generate the root node (global) transform
      if it's not already generated */
   if (rt->rtree->left != NULL) {
      if (rt->rtree->nodeInfo->nodeOcc >= rt->nodeOccThresh &&
          rt->rtree->left->nodeInfo->nodeOcc >= rt->nodeOccThresh && 
          rt->rtree->right->nodeInfo->nodeOcc >= rt->nodeOccThresh && 
          rt->rtree->left->nodeInfo->nodeComps > rt->vSize &&
          rt->rtree->right->nodeInfo->nodeComps > rt->vSize)
         GetMeanTransform(rt->rtree->nodeInfo, &rt->transMem, rt);
   }

   if (trace & T_MEM) {
      PrintHeapStats(rt->hmem);
      PrintHeapStats(&rt->pdfStatsMem);
      PrintHeapStats(&rt->transMem);
      fflush(stdout);
   }
    

}

/* ----------------------------------------------------------------------*/
/*           Covariance Transform Calculation Functions                  */
/* ----------------------------------------------------------------------*/

/* Calculate the covariance transformation for a regression class */
static void GetCovarTransform(RegNode *n, RegTransInfo *rt) 
{
   short idx, i;
   int k, vSize;
   float sumOcc = 0.0;

   ZeroVector(Hg);
   vSize = rt->vSize;

   /* Sum the regression accumulates for the base classes 
      belonging to this node */
   for (i = 1; i <= n->nBases; i++) {
      idx = n->bases[i];
      if (rt->baseforms[idx] == NULL)
         HError(7430, "GetCovarTransform:Can't find base class %d for node %d", i, 
                n->nodeIndex);
      for (k = 1; k<=vSize; k++)
         Hg[k] += rt->baseforms[idx]->H[k];
   }

   sumOcc = 1.0 / n->nodeOcc;

   for (k = 1; k <= vSize; k++)
      n->HTrans[k] = Hg[k] * sumOcc;

}


/* Calculate the covariance transformations for a regression tree */
static void GetClassCovarTransforms(RegTree *t, RegTransInfo *rt) {
  
   int vSize;

   vSize = rt->vSize;

   if (t != NULL) {
      GetClassCovarTransforms(t->left, rt);

      if (t->nodeInfo->nodeOcc > rt->nodeOccThresh &&
          t->nodeInfo->nodeComps > vSize) {
         if (t->left != NULL) {
            if (CheckNodeOcc(t->left, t->right, rt)) {
               if (trace & T_DET) {
                  printf("(Covar) Building transform for node %d\n",
                         t->nodeInfo->nodeIndex);
                  fflush(stdout);
               }
               GetCovarTransform(t->nodeInfo, rt);
            }
         }
         else {
            if (trace & T_DET) {
               printf("(Covar) Building transform for node %d\n",
                      t->nodeInfo->nodeIndex);
               fflush(stdout);
            }
            GetCovarTransform(t->nodeInfo, rt);
         }
      }
      GetClassCovarTransforms(t->right, rt);
   }

}

/* Calculate the covariance accumulation H for a base regression class */
static int GetClassCovarAccumulates(RegTransInfo *rt, BaseformAcc *base) 
{
   int nCompsSeen=0, i, n;
   MixPDF *mp=NULL;
   RegAcc *ra=NULL;
   HMMSet *hset;

   hset = rt->hset;
   n = base->nComponents;

   for (i = 1; i <= n; i++) {
      mp = base->pdfList[i];
      ra = GetRegAcc(mp);
      if (ra->occ > MINOCC) {
         AccCovarTrans(mp, base->H, mp->mean, mp->cov.var, 
                       ra->obsSum, ra->obsSqSum);
         nCompsSeen += 1;
      }
   }


   return nCompsSeen;
}

/* EXPORT->CalcCovarTransforms: Calculate the covariance transforms for
   every regression class */
void CalcCovarTransforms(RegTransInfo *rt) {

   int i;
   int compsSeen=0;

   if (rt->transKind & MEANVAR) {

      for (i = 1; i <= rt->nBaseTransforms; i++) {
         if (trace & T_DET) {
            printf("(Covars) Accumulating at the regression level for base %d\n",
                   rt->baseforms[i]->baseformIndex);
            fflush(stdout);
         }
         compsSeen += GetClassCovarAccumulates(rt, rt->baseforms[i]);
      }
    
      GetClassCovarTransforms(rt->rtree, rt);
    
      /* always make sure to generate the root node (global) transform */
      if (rt->rtree->left != NULL) {
         if (rt->rtree->nodeInfo->nodeOcc >= rt->nodeOccThresh &&
             rt->rtree->left->nodeInfo->nodeOcc >= rt->nodeOccThresh && 
             rt->rtree->right->nodeInfo->nodeOcc >= rt->nodeOccThresh && 
             rt->rtree->left->nodeInfo->nodeComps > rt->vSize &&
             rt->rtree->right->nodeInfo->nodeComps > rt->vSize)
            GetCovarTransform(rt->rtree->nodeInfo, rt);
      }
   }
}

/* ----------------------------------------------------------------------*/
/*                        MAP Adaptation                                 */
/* ----------------------------------------------------------------------*/
void UpdateMAP(RegTransInfo *rt, float tau)
{
   HMMScanState hss;
   MixtureElem *me;
   WtAcc *wa;
   HLink hmm;
   RegAcc *ra;
   Vector mean, var;
   float scaleOcc, occ;
   float x, occi, wght;
   int i, vSize, m, M;
 
   vSize = rt->vSize;
   NewHMMScan(rt->hset,&hss);
   do {
      hmm = hss.hmm;
      while (GoNextState(&hss,TRUE)) {
         while (GoNextStream(&hss,TRUE)) {            
            if (hss.isCont) {                    /* PLAINHS or SHAREDHS */
               while (GoNextMix(&hss,TRUE)) {
                  ra = GetRegAcc(hss.mp);
                  mean = hss.mp->mean;
                  var = hss.mp->cov.var;
                  occ = ra->occ;
                  if (MixWeight(rt->hset, hss.me->weight) > MINMIX) {
                     if (occ > 0.0) {
                        scaleOcc = 1.0/(occ + tau);
                        for (i = 1; i <= vSize; i++) {
                           mean[i] *= tau;
                           mean[i] += ra->obsSum[i];
                           mean[i] *= scaleOcc;
                        } 
                     }
                  }
               }
               wa = (WtAcc *)hss.ste->hook;
               M=hss.ste->nMix;
               if (wa == NULL)
                  HError(7470, "UpdateMAP: weight hook is NULL!");
               occi = wa->occ;
               if (occi>0) {
                  scaleOcc = 1.0/(occi + tau);
                  for (m=1; m<=M; m++){
                     x = wa->c[m]/occi;
                     if (x>1.0)
                        x = 1.0;
                     me = hss.ste->spdf.cpdf+m;
                     wght = tau * MixWeight(rt->hset, hss.me->weight);
                     wght += occi * x;
                     wght *= scaleOcc;
                     me->weight = (wght>MINMIX) ? wght  : 0.0;
                  }
               }
            }
            else
               HError(7450, "UpdateMAP: Adaptation only available for PLAIN or SHARED systems!");
         }
      }
      FixGConsts(hmm);
   } while (GoNextHMM(&hss));
   EndHMMScan(&hss);

   if (trace & T_MEM) {
      PrintHeapStats(rt->hmem);
      PrintHeapStats(&rt->pdfStatsMem);
      fflush(stdout);
   }

}


/* ----------------------------------------------------------------------*/
/*                        MLLR Adaptation                                */
/* ----------------------------------------------------------------------*/


/* EXPORT->DoAdaptation: Given the initialisation to the adaptation */
void DoAdaptation(RegTransInfo *rt, Boolean global) {

   int i;
   RegTree *rtree;
   float nFrames=0.0;
   Boolean usedGlobal=FALSE;

   rtree = rt->rtree;

   if (trace&T_AUX) {
      printf("Calculating auxilliary before transformation...");
      CalcAuxilliary(rt);
      fflush(stdout);
   }

   /* check to see if transform is a global one */
   if (global || !strcmp(rt->transId->rcId, "global"))
      usedGlobal = TRUE;
  
   if (usedGlobal) {
      /* reset the usedGlobal marker for iterative approach by changing
         the rcId in transId */
      if (!global)
         rt->transId->rcId[0] = '\0';
      ApplyBackwardGlobalTransform(rt);
   }
   else
      ApplyBwdTransforms(rtree, rt);
 
   CalcMeanTransforms(rt);
  
   if (global) {
      ApplyMeanGlobalTransform(rt);
   }
   else
      ApplyMeanTransforms(rt, rt->rtree);

   if ((trace&T_AUX) && (rt->transKind == MEANONLY)) {
      printf("Calculating auxilliary after mean transformation...");
      CalcAuxilliary(rt);
      for (i = 1; i <= rt->nBaseTransforms; i++)
         nFrames += rt->baseforms[i]->occ;
      /* printf(" Total auxilliary = %e, for %d frames\n",
         -0.5*currAux/nFrames, (int) floor(nFrames + 0.5)); */
      fflush(stdout);
   }

   if (rt->transKind & MEANVAR) {
      CalcCovarTransforms(rt);
      if (global)
         ApplyCovarGlobalTransform(rt);
      else
         ApplyCovarTransforms(rt, rt->rtree);

      if (trace&T_AUX) {
         printf("Calculating auxilliary after mean & variance transformation...");
         CalcAuxilliary(rt);
         for (i = 1; i <= rt->nBaseTransforms; i++)
            nFrames += rt->baseforms[i]->occ;
         /* printf(" Total auxilliary = %e, for %d frames\n",
            -0.5*currAux/nFrames, (int) floor(nFrames + 0.5)); */
         fflush(stdout);
      }
   }

}

/* ----------------------------------------------------------------------*/
/*                       END of HAdapt.c Module                          */
/* ----------------------------------------------------------------------*/

