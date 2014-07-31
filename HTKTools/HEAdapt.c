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
/*      Entropic Cambridge Research Laboratory                 */
/*      (now part of Microsoft)                                */
/*                                                             */
/* ----------------------------------------------------------- */
/*         Copyright: Microsoft Corporation                    */
/*          1995-2000 Redmond, Washington USA                  */
/*                    http://www.microsoft.com                 */
/*                                                             */
/*              2002  Cambridge University                     */
/*                    Engineering Department                   */
/*                                                             */
/*   Use of this software is governed by a License Agreement   */
/*    ** See the file License for the Conditions of Use  **    */
/*    **     This banner notice must not be removed      **    */
/*                                                             */
/* ----------------------------------------------------------- */
/*             File: HEAdapt.c: Adaptation Tool                */
/* ----------------------------------------------------------- */

char *headapt_version = "!HVER!HEAdapt:   3.2.1 [CUED 15/10/03]";
char *headapt_vc_id = "$Id: HEAdapt.c,v 1.11 2003/10/15 08:10:13 ge204 Exp $";

#include "HShell.h"     /* HMM ToolKit Modules */
#include "HMem.h"
#include "HMath.h"
#include "HSigP.h"
#include "HAudio.h"
#include "HWave.h"
#include "HVQ.h"
#include "HParm.h"
#include "HLabel.h"
#include "HModel.h"
#include "HTrain.h"
#include "HUtil.h"
#include "HAdapt.h"
#include "HFB.h"


/* Trace Flags */
#define T_TOP   0001    /* Top level tracing */
#define T_MAP   0002    /* logical/physical hmm map */

Boolean traceHFB = FALSE;        /* pass to HFB to retain top-level tracing */

/* Global Settings */

static char * labDir = NULL;      /* label (transcription) file directory */
static char * labExt = "lab";     /* label file extension */
static char * hmmDir = NULL;      /* directory to look for hmm def files */
static char * hmmExt = NULL;      /* hmm def file extension */
static char * newDir = NULL;      /* directory to store new hmm def files */
static char * newExt = NULL;      /* extension of new reestimated hmm files */
static char * saveTransFN=NULL;   /* stats file, if any */
static char * uid ="unknown";     /* temporary uid */
static char * uname ="Unknown";   /* temporary uname */
static char * chan ="Standard";   /* temporary mic channel name */
static char * desc ="None";       /* temporary description to save in tmf */
static char * transFile = NULL;   /* file holding transformations */
static int blocks=0;              /* number of block diagonal matrices */
static RegClassType regClass=ADPTUNDEF;/* regClass type to be used */
static RegTransType regTrans=TRANSUNDEF;/* reg transform type to be used */

static float minVar  = 0.0;       /* minimum variance (diagonal only) */
static float occThresh=0.0;       /* Regression tree ocupation count minimum */
static int update    = 0;         /* update trans every "update" utterances */

static int trace     = 0;         /* Trace level */
static Boolean saveBinary=FALSE;  /* save output in binary  */
static FileFormat dff=UNDEFF;     /* data file format */
static FileFormat lff=UNDEFF;     /* label file format */

/* Global Data Structures - valid for all training utterances */
static LogDouble pruneInit = NOPRUNE;    /* pruning threshold initially */
static LogDouble pruneInc = 0.0;         /* pruning threshold increment */
static LogDouble pruneLim = NOPRUNE;     /* pruning threshold limit */
static float minFrwdP = 10.0;    /* mix prune threshold */

static float tau         = 15.0;    /* scaling for MAP adaptation */
static Boolean map       = FALSE;   /* do MAP adaptation */
static Boolean mapUseMLLR= FALSE;   /* do MAP + MLLR adaptation */
static Boolean global    = FALSE;   /* do global adaptation only */

static Boolean firstTime = TRUE;    /* Flag used to enable creation of ot */
static Boolean twoDataFiles = FALSE; /* Enables creation of ot2 for FB
                                        training using two data files */
static Vector vFloor[SMAX]; /* variance floor - default is all zero */

static LogDouble totalPr=0.0;   /* total log prob upto current utterance */
static int totalT=0;            /* total number of frames in training data */

static MemHeap hmmStack;   /*For Storage of all dynamic structures created...*/
static MemHeap uttStack;
static MemHeap regStack;
static MemHeap fbInfoStack;

/* ------------------ Process Command Line -------------------------- */
   
/* SetConfParms: set conf parms relevant to HCompV  */
void SetConfParms(void)
{
   int i;
   Boolean b;
   ConfParam *cParm[MAXGLOBS];/* configuration parameters */
   int nParm = 0;             /* total num params */
   
   nParm = GetConfig("HEADAPT", TRUE, cParm, MAXGLOBS);
   if (nParm>0) {
      if (GetConfInt(cParm,nParm,"TRACE",&i)) trace = i;
      if (GetConfBool(cParm,nParm,"SAVEBINARY",&b)) saveBinary = b;
   }
}

void ReportUsage(void)
{
   printf("\nUSAGE: HEAdapt [options] hmmList dataFiles...\n\n");
   printf(" Option                                       Default\n\n");
   printf(" -b N    set no. of blocks for transform      1\n");
   printf(" -c f    Mixture pruning threshold            10.0\n");
   printf(" -d s    dir to find hmm definitions          current\n");
   printf(" -f s1 s2  Save s2 for field s1 in tmf        defaults saved\n");
   printf(" -g      do global adaptation only            reg classes\n");
   printf(" -i N    update the transforms every N utts   static\n");
   printf(" -j f    Map adaptation with scaling factor   off\n");
   printf(" -k      Use MLLR before performing MAP       off\n");
   printf(" -m f    occ count min. to calc transform     700.0\n");
   printf(" -o s    extension for new hmm files          as src\n");
   printf(" -t f [i l] set pruning to f [inc limit]      inf\n");
   printf(" -u [mv] update means(m) vars(v)              means\n");
   printf(" -v f    set minimum variance to f            0.0\n");
   printf(" -x s    extension for hmm files              none\n");
   PrintStdOpts("BFGHIJKLMSTX");
   printf("\n\n");
}

void SetuFlags(void)
{
   char *s;
   Boolean meanFlag=FALSE;
   
   s=GetStrArg();
   while (*s != '\0')
      switch (*s++) {
      case 'm':
         meanFlag = TRUE;
         if (regTrans == TRANSUNDEF)
            regTrans=MEANONLY;
         break;
      case 'v':regTrans=MEANVAR;break;
      default: HError(2720,"SetuFlags: Unknown update flag %c",*s);
         break;
      }
   if (!meanFlag)
      HError(-2721,"SetuFlags: Mean transformation will be updated, even though unrequestred");
}

int main(int argc, char *argv[])
{
   int numUtt;
   char fmt[MAXSTRLEN];
   char *datafn, *s, *c;
   UttInfo *utt;            /* utterance information storage */
   RegTransInfo *rt;        /* regression transform information storage */
   HMMSet hset;             /* Set of HMMs to be re-estimated */
   FBInfo *fbInfo;          /* forward-backward information storage */

   void Initialise(FBInfo *fbInfo, MemHeap *x, HMMSet *hset,
                   char *hmmListFn, RegTransInfo *rt, char *transFile);
   void DoForwardBackward(FBInfo *fbInfo, UttInfo *utt, char *datafn, char *datafn2);
   
   if(InitShell(argc,argv,headapt_version,headapt_vc_id)<SUCCESS)
      HError(2700,"HEAdapt: InitShell failed");

   InitMem();
   InitMath();
   InitSigP();
   InitAudio();
   InitWave();
   InitVQ();
   InitModel();
   if(InitParm()<SUCCESS)  
      HError(2700,"HEAdapt: InitParm failed"); 
   InitLabel();
   InitTrain();
   InitUtil();
   InitFB();
   InitAdapt();  

   if (!InfoPrinted() && NumArgs() == 0)
      ReportUsage();
   if (NumArgs() == 0) Exit(0);

   SetConfParms();
   CreateHeap(&hmmStack,"HmmStore", MSTAK, 1, 1.0, 50000, 500000);
   CreateHMMSet(&hset,&hmmStack,TRUE);
   CreateHeap(&uttStack,   "uttStore",    MSTAK, 1, 0.5, 100,   1000);
   utt = (UttInfo *) New(&uttStack, sizeof(UttInfo));
   CreateHeap(&fbInfoStack,   "FBInfoStore",  MSTAK, 1, 0.5, 1000,  10000);
   fbInfo = (FBInfo *) New(&fbInfoStack, sizeof(FBInfo));
   CreateHeap(&regStack,   "regClassStore",    MSTAK, 1, 0.5, 500,   5000);
   rt = (RegTransInfo *) New(&regStack, sizeof(RegTransInfo));

   while (NextArg() == SWITCHARG) {    
      s = GetSwtArg();
      if (strlen(s)!=1) 
         HError(2719,"HEAdapt: Bad switch %s; must be single letter",s);
      switch(s[0]){
      case 'b':
         if (NextArg()!= INTARG)
            HError(2719,"HEAdapt: Number of blocks in transform matrix expected");
         blocks = GetIntArg(); break; 
      case 'c':
         minFrwdP = GetChkedFlt(0.0,1000.0,s);
         break;
      case 'd':
         if (NextArg()!=STRINGARG)
            HError(2719,"HEAdapt: HMM definition directory expected");
         hmmDir = GetStrArg(); break;   
      case 'f':
         if (NextArg()!=STRINGARG)
            HError(2719,"HEAdapt: TMF description type expected");
         strcpy(fmt,GetStrArg());
         for (c=fmt; *c!=0; *c++=toupper(*c));
         if (strcmp(fmt,"UID")==0) {
            if (NextArg()!=STRINGARG)
               HError(2719,"HEAdapt: user id string expected");
            uid = GetStrArg();
         }
         else if (strcmp(fmt,"UNAME")==0) {
            if (NextArg()!=STRINGARG)
               HError(2719,"HEAdapt: user name string expected");
            uname = GetStrArg();
         }
         else if (strcmp(fmt,"CHAN")==0) {
            if (NextArg()!=STRINGARG)
               HError(2719,"HEAdapt: channel description string expected");
            chan = GetStrArg();
         }
         else if (strcmp(fmt,"DESC")==0){
            if (NextArg()!=STRINGARG)
               HError(2719,"HEAdapt: a general description string expected");
            desc = GetStrArg();
         }
         else
            HError(2719,"HEAdapt: Unrecognised format, should be one of [uid,uname,chan,desc]");
         break;
      case 'g':
         /* global adaptation! */
         global = TRUE;
         break;
      case 'i':
         if (NextArg()!=INTARG)
            HError(2719,"HEAdapt: Number of utterances before transformation expected");
         update = GetIntArg(); break;
      case 'j':
         map = TRUE;
         if (NextArg()==FLOATARG || NextArg()==INTARG)
            tau = GetChkedFlt(0.0,1.0E20,s);
         break;
      case 'k':
         mapUseMLLR = TRUE;
         break;
      case 'm':
         occThresh = GetChkedFlt(0.0,1.0E20,s); break;
      case 'n':
         if (NextArg()!=STRINGARG)
            HError(2719,"HEAdapt: Full name string expected");
         uname = GetStrArg(); break;
      case 'o':
         if (NextArg()!=STRINGARG)
            HError(2719,"HEAdapt: HMM file extension expected");
         newExt = GetStrArg(); break;
      case 't':
         pruneInit = GetChkedFlt(0.0,1.0E20,s);
         if (NextArg()==FLOATARG || NextArg()==INTARG)
            {
               pruneInc = GetChkedFlt(0.0,1.0E20,s);
               pruneLim = GetChkedFlt(0.0,1.0E20,s);
            }
         else
            {
               pruneInc = 0.0;
               pruneLim = pruneInit;
            }
         break;
      case 'u':
         SetuFlags(); break;
      case 'v':
         minVar = GetChkedFlt(0.0,10.0,s); break;
      case 'x':
         if (NextArg()!=STRINGARG)
            HError(2319,"HERest: HMM file extension expected");
         hmmExt = GetStrArg(); break;
      case 'B':
         saveBinary=TRUE;
         break;
      case 'F':
         if (NextArg() != STRINGARG)
            HError(2719,"HEAdapt: Data File format expected");
         if((dff = Str2Format(GetStrArg())) == ALIEN)
            HError(-2789,"HEAdapt: Warning ALIEN Data file format set");
         break;
      case 'G':
         if (NextArg() != STRINGARG)
            HError(2719,"HEAdapt: Label File format expected");
         if((lff = Str2Format(GetStrArg())) == ALIEN)
            HError(-2789,"HEAdapt: Warning ALIEN Label file format set");
         break;
      case 'H':
         if (NextArg() != STRINGARG)
            HError(2719,"HEAdapt: HMM macro file name expected");
         AddMMF(&hset,GetStrArg());
         break;     
      case 'I':
         if (NextArg() != STRINGARG)
            HError(2719,"HEAdapt: MLF file name expected");
         LoadMasterFile(GetStrArg());
         break;
      case 'J':
         if (NextArg()!=STRINGARG)
            HError(2719,"HEAdapt: Transforms file name expected");
         transFile = GetStrArg(); break;
      case 'K':
         if (NextArg()!=STRINGARG)
            HError(2719,"HEAdapt: Stats file name expected");
         saveTransFN = GetStrArg(); break;
      case 'L':
         if (NextArg()!=STRINGARG)
            HError(2719,"HEAdapt: Label file directory expected");
         labDir = GetStrArg(); break;
      case 'M':
         if (NextArg()!=STRINGARG)
            HError(2719,"HEAdapt: Output model file directory expected");
         newDir = GetStrArg();
         break;     
      case 'T':
         trace = GetChkedInt(0,0100000,s);
         break;
      case 'X':
         if (NextArg()!=STRINGARG)
            HError(2719,"HEAdapt: Label file extension expected");
         labExt = GetStrArg(); break;
      default:
         HError(2719,"HEAdapt: Unknown switch %s",s);
      }
   } 
   if (NextArg() != STRINGARG)
      HError(2719,"HEAdapt: file name of vocabulary list expected");

   /* set number of blocks and regression class kind here */
   rt->nBlocks = blocks;
   rt->classKind = regClass;
   rt->transKind = regTrans;
   rt->nodeOccThresh = occThresh;
   rt->adptSil = TRI_UNDEF;

   if (map && update > 0) {
      if (mapUseMLLR) {
         HError(-2741, "HEAdapt: MLLR + MAP adaptation only available is static mode -- switching to static mode");
         update = 0; 
      }
      else {
         HError(-2740, "HEAdapt: MAP adaptation only available is static mode -- switching to static mode");
         update = 0;
      }
   } 
     
   Initialise(fbInfo, &fbInfoStack, &hset, GetStrArg(), rt, transFile);

   InitUttInfo(utt, twoDataFiles);
   numUtt = 0;

   do {
      if (NextArg()!=STRINGARG)
         HError(2719,"HEAdapt: data file name expected");
      datafn = GetStrArg();
      DoForwardBackward(fbInfo, utt, datafn, NULL) ;
      if (update > 0) {
         if (((numUtt+1) % update) == 0) {
            DoAdaptation(rt, global);
            ClearRegCompStats(&hset, rt);
         }
      }
      numUtt += 1;
   } while (NumArgs()>0);
   
   if (map) {
      if (mapUseMLLR)
         DoAdaptation(rt, global);
      UpdateMAP(rt, tau);
   }
   else
      if (update == 0 || (update > 0 && numUtt%update > 0))
         DoAdaptation(rt, global);

   if (saveTransFN == NULL || map) {    
      if (trace&T_TOP && map)
         printf("MAP Reest complete -\n\t\taverage log prob per frame = %e (%d)\n",
                totalPr/totalT, totalT); 
      ConvDiagC(&hset,TRUE);
      if(SaveHMMSet(&hset,newDir,newExt,saveBinary)<SUCCESS)
         HError(2710,"HEAdapt: Saving HMMSet failed");
   }
   
   if (saveTransFN != NULL) {
      SaveTransformSet(&hset, rt, saveTransFN, newDir, 
                       uid, uname, chan, desc, FALSE, global, saveBinary);
   }

   Exit(0);
   return (0);          /* never reached -- make compiler happy */
}

/* -------------------------- Initialisation ----------------------- */

void Initialise(FBInfo *fbInfo, MemHeap *x, HMMSet *hset,
                char *hmmListFn, RegTransInfo *rt, char *transFile)
{   

   Boolean loadTransStats=FALSE;

   /* Load HMMs and init HMMSet related global variables */
   if(MakeHMMSet( hset, hmmListFn )<SUCCESS)
      HError(2728,"Initialise: MakeHMMSet failed");
   if(LoadHMMSet( hset,hmmDir,hmmExt)<SUCCESS)
      HError(2728,"Initialise: LoadHMMSet failed");
   SetParmHMMSet(hset);
  
   /* needs to be plain or shard system for adaptation purposes */
   if (hset->hsKind == TIEDHS || hset->hsKind == DISCRETEHS)
      HError(2730, "Initialise: MLLR adaptation only available for PLAIN or SHARED systems!"); 
  
   /* needs to currently be single stream data */
   if (hset->swidth[0] != 1)
      HError(2731, "Initialise: MLLR adaptation only currently available for single stream data");

   ConvDiagC(hset,TRUE);

   InitialiseTransform(hset, &regStack, rt, TRUE);
   if (trace & T_TOP) {
      printf("Updating Mean ");
      if (rt->transKind == MEANVAR)
         printf("and Variance ");
      printf("Transforms\n");
   }
   if (transFile != NULL) {
      LoadTransformSet(hset, transFile, uid, rt, &loadTransStats);
      ApplyTransforms(rt);
   }
   InitialiseAdapt(hset, &regStack, rt);

   SetVFloor(hset, vFloor, minVar);
  
   /* initialise and pass information to the forward backward library */
   InitialiseForBack(fbInfo, x, hset, rt, (UPDSet) (UPADAPT|UPMIXES), pruneInit, 
                     pruneInc, pruneLim, minFrwdP);

}

/* Load data and call FBFile: apply forward-backward to given utterance */
void DoForwardBackward(FBInfo *fbInfo, UttInfo *utt, char * datafn, char * datafn2)
{

   utt->twoDataFiles = twoDataFiles ;
   utt->S = fbInfo->up_hset->swidth[0];

   /* Load the labels */
   LoadLabs(utt, lff, datafn, labDir, labExt);
   /* Load the data */
   LoadData(fbInfo->up_hset, utt, dff, datafn, datafn2);
   if (firstTime) {
      InitUttObservations(utt, fbInfo->up_hset, datafn, fbInfo->maxMixInS);
      firstTime = FALSE;
   }

   /* fill the alpha beta and otprobs (held in fbInfo) */
   if (FBFile(fbInfo, utt, datafn)) {
      /* update totals */
      totalT += utt->T ;
      totalPr += utt->pr ;
   }
}

