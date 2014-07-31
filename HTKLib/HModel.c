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
/*               2002 Cambridge University                     */
/*                    Engineering Department                   */
/*                                                             */
/*   Use of this software is governed by a License Agreement   */
/*    ** See the file License for the Conditions of Use  **    */
/*    **     This banner notice must not be removed      **    */
/*                                                             */
/* ----------------------------------------------------------- */
/*         File: HModel.c  HMM Model Definition Data Type      */
/* ----------------------------------------------------------- */

char *hmodel_version = "!HVER!HModel:   3.2.1 [CUED 15/10/03]";
char *hmodel_vc_id = "$Id: HModel.c,v 1.13 2003/10/15 08:10:12 ge204 Exp $";

#include "HShell.h"
#include "HMem.h"
#include "HMath.h"
#include "HWave.h"
#include "HAudio.h"
#include "HParm.h"
#include "HLabel.h"
#include "HModel.h"
#include "HUtil.h"
#include "HTrain.h"
#include "HAdapt.h"

/* --------------------------- Trace Flags ------------------------- */

static int trace = 0;

#define T_TOP  00001       /* Top Level tracing */
#define T_CHK  00002       /* Show HMM Checking */
#define T_TOK  00004       /* Scanner Token Processing */
#define T_PAR  00010       /* Trace Parsing */
#define T_PMP  00020       /* Pointer Map Tracing */
#define T_MAC  00040       /* Macro Load/Save Tracing */
#define T_ORP  00100       /* Orphan macro handling */
#define T_BTR  00200       /* Decision Tree handling */
#define T_GMX  00400       /* GMP optimisation */
#define T_XFM  01000       /* Loading of xform macros */
#define T_XFD  02000       /* Additional detail of loading of xform macros */

#define CREATEFIDX -1
#define LOADFIDX   -2

/* ------------------ Input XForm directory info ------------------- */

typedef struct _XFDirInfo *XFDirLink;

typedef struct _XFDirInfo {
   char *dirName;           /* input XForm directory name */
   XFDirLink next;          /* next directory name in list */
} XFDirInfo;

/* --------------------------- Initialisation ---------------------- */

static ConfParam *cParm[MAXGLOBS];      /* config parameters */
static int nParm = 0;
static Boolean checking   = TRUE;       /* check HMM defs */
static Boolean saveBinary = FALSE;      /* save HMM defs in binary */
static Boolean saveGlobOpts = TRUE;     /* save ~o with HMM defs */
static Boolean saveRegClass = TRUE;     /* save regression classes and tree */ 
static Boolean forceHSKind= FALSE;      /* force HMM Set Kind */
static Boolean keepDistinct=FALSE;      /* keep orphan HMMs distinct */
static Boolean discreteLZero=FALSE;     /* map DLOGZERO to LZERO */

static Boolean allowOthers=TRUE;        /* allow unseen models in files */
static HSetKind cfHSKind;
static char orphanMacFile[100]; /* last resort file for new macros */

static XFDirLink xformDirNames = NULL;  /* linked list of input transform directories */
static Boolean saveInputXForm = TRUE;   /* save input xforms with models set */
static MemHeap xformStack;              /* For Storage of xforms with no model sets ... */

void InitSymNames(void);

/* EXPORT->InitModel: initialise memory and configuration parameters */
void InitModel(void)
{
   int i;
   Boolean b;
   char buf[MAXSTRLEN];
   
   Register(hmodel_version,hmodel_vc_id);
   CreateHeap(&xformStack,"XFormStore",MSTAK, 1, 0.5, 100 ,  1000 );
   strcpy(orphanMacFile,"newMacros");
   InitSymNames();
   nParm = GetConfig("HMODEL", TRUE, cParm, MAXGLOBS);
   if (nParm>0){
      if (GetConfInt(cParm,nParm,"TRACE",&i)) trace = i;
      if (GetConfBool(cParm,nParm,"CHKHMMDEFS",&b)) checking = b;
      if (GetConfBool(cParm,nParm,"SAVEBINARY",&b)) saveBinary = b;
      if (GetConfBool(cParm,nParm,"KEEPDISTINCT",&b)) keepDistinct = b;
      if (GetConfBool(cParm,nParm,"SAVEGLOBOPTS",&b)) saveGlobOpts = b;
      if (GetConfBool(cParm,nParm,"SAVEREGCLASS",&b)) saveRegClass = b;
      if (GetConfBool(cParm,nParm,"SAVEINPUTXFORM",&b)) saveInputXForm = b;
      if (GetConfBool(cParm,nParm,"ALLOWOTHERHMMS",&b)) allowOthers = b;
      if (GetConfBool(cParm,nParm,"DISCRETELZERO",&b)) discreteLZero = b;
      if (GetConfStr (cParm,nParm,"ORPHANMACFILE",buf))
         strcpy(orphanMacFile,buf);
      if (GetConfStr (cParm,nParm,"HMMSETKIND",buf)) {
         if (strcmp(buf,"PLAIN")==0) cfHSKind = PLAINHS;
         else if (strcmp(buf,"SHARED")==0) cfHSKind = SHAREDHS;
         else if (strcmp(buf,"TIED")==0) cfHSKind = TIEDHS;
         else if (strcmp(buf,"DISCRETE")==0) cfHSKind = DISCRETEHS;
         else
            HError(7070,"InitModel: Unknown HMM kind %s",buf);
         forceHSKind = TRUE;
      }
   }
}

/* -------------------- Check Model Consistency -------------------- */

/* CheckMix: check given mixture pdf of state n, stream s, mixture m */
static ReturnStatus CheckMix(char *defName, MixPDF *mp, int n, int s, int m, int sw)
{  
   if (trace&T_CHK) {
      printf("HModel:       checking mix %d\n",m);  fflush(stdout);
   }
   if (mp->mean == NULL){
      HRError(7030,"CheckMix: %s: pdf for s=%d,j=%d,m=%d has NULL mean",
              defName,n,s,m);
      return(FAIL);
   }
   if (VectorSize(mp->mean) != sw){
      HRError(7030,"CheckMix: %s: mean for s=%d,j=%d,m=%d has bad vecSize",
              defName,n,s,m);
      return(FAIL);
   }
   if (mp->cov.var == NULL){
      HRError(7030,"CheckMix: %s: pdf for s=%d,j=%d,m=%d has NULL covariance",
              defName,n,s,m);
      return(FAIL);
   }
   
   switch(mp->ckind){
   case DIAGC:
      if (VectorSize(mp->cov.var) != sw){
         HRError(7030,"CheckMix: %s: var for s=%d,j=%d,m=%d has bad vecSize",
                 defName,n,s,m);
         return(FAIL);
      }
      break;
   case FULLC:
   case LLTC:
      if (TriMatSize(mp->cov.inv) != sw){
         HRError(7030,"CheckMix: %s: inv for s=%d,j=%d,m=%d has bad dimens",
                 defName,n,s,m);
         return(FAIL);
      }
      break;
   case XFORMC:
      if (NumCols(mp->cov.xform) != sw){
         HRError(7030,"CheckMix: %s: xform for s=%d,j=%d,m=%d has bad dimens",
                 defName,n,s,m);
         return(FAIL);
      }
      break;
   }
   if (mp->gConst == LZERO){
      if (mp->ckind==DIAGC)
         FixDiagGConst(mp);
      else if (mp->ckind==FULLC)
         FixFullGConst(mp,-CovDet(mp->cov.inv));
      else if (mp->ckind==LLTC)
         FixLLTGConst(mp);
      else if (mp->ckind==XFORMC)
         HRError(7030,"CheckMix: No GConst set for xformc covariance matrix");
   }
   Touch(&mp->nUse);
   return(SUCCESS);
}

/* CheckStream: check the stream s for state n */
static ReturnStatus CheckStream(char *defName, HLink hmm, StreamElem *se, int s, int n)
{
   double sum=0.0;
   int m,sw;
   LogFloat wt;
   MixtureElem *me;
   MixPDF *mp;
   HSetKind  hk;
   
   if (trace&T_CHK) {
      printf("HModel:    checking stream %d, se=%p\n",s,se); 
      fflush(stdout);
   }
   hk = hmm->owner->hsKind;
   sw = hmm->owner->swidth[s];
   switch(hk){
   case PLAINHS:
   case SHAREDHS:
      me = se->spdf.cpdf+1;
      for (m=1; m<=se->nMix; m++,me++){
         wt=me->weight; wt = MixWeight(hmm->owner,wt);
         sum += wt;
         if (wt > MINMIX ){
            if (me->mpdf == NULL){
               HRError(7030,"CheckStream: %s: MixDef %d missing for s=%d, j=%d",
                       defName, m, s, n);
               return(FAIL);
            }
            mp = me->mpdf;
            if((!IsSeen(mp->nUse))&&(CheckMix(defName,me->mpdf,n,s,m,sw)<SUCCESS))
               return(FAIL);
         }
      }
      break;
   case TIEDHS:
      for (m=1; m<=se->nMix; m++)
         sum += se->spdf.tpdf[m];
      break;
   case DISCRETEHS:
      for (m=1; m<=se->nMix; m++)
         sum += exp((float)se->spdf.dpdf[m]/DLOGSCALE);
      break;
   }
   if (sum<0.99 || sum>1.01){
      HRError(7031,"CheckStream: %s: Mix weights sum %e for s=%d, j=%d",
              defName,sum,s,n);
      return(FAIL);
   }
   return(SUCCESS);
}

/* CheckState: check state n */
static ReturnStatus CheckState(char *defName, HLink hmm, StateInfo *si, int n)
{
   int s,S;
   StreamElem *ste;
   
   if (trace&T_CHK) {
      printf("HModel:  checking state %d, si=%p\n",n,si); fflush(stdout);
   }
   S = hmm->owner->swidth[0];
   if (si->pdf == NULL){
      HRError(7030,"CheckState: %s: state %d has no pdf",defName,n);
      return(FAIL);
   }
   ste = si->pdf+1;
   for (s=1; s<=S; s++,ste++){
      if (ste->spdf.cpdf == NULL){
         HRError(7030,"CheckState: %s: Stream %d missing in state %d",defName,s,n);
         return(FAIL);
      }
      if(CheckStream(defName,hmm,ste,s,n)<SUCCESS)
         return(FAIL);
   }
   if (S>1)
      if (si->weights==NULL){
         HRError(7030,"CheckState: %s: Stream weights missing in state %d",defName,n);
         return(FAIL);
      }
   Touch(&si->nUse);
   return(SUCCESS);
}

/* CheckHMM: check the consistency of the given model */
static ReturnStatus CheckHMM(char *defName, HLink hmm)
{
   int i;
   StateElem *se;
   StateInfo *si;
   
   if (trace&T_CHK) {
      printf("HModel: checking HMM %s\n",defName); fflush(stdout);
   }
   for (i=2,se=hmm->svec+2; i<hmm->numStates; i++,se++) {
      si = se->info;
      if (si == NULL){
         HRError(7030,"CheckHMM: %s: state %d has NULL info",defName,i);
         return(FAIL);
      }
      else if((!IsSeen(si->nUse))&&(CheckState(defName,hmm,si,i)<SUCCESS))
         return(FAIL);
   }
   return(SUCCESS);
}

/* CheckTMRecs: check the tied mix codebook attached to a HMM Set */
static ReturnStatus CheckTMRecs(HMMSet *hset)
{
   int s,m,sw;
   MixPDF *mp;
   
   if (trace&T_CHK) {
      printf("HModel: checking tied mixture codebook\n"); fflush(stdout);
   }
   for (s=1;s<=hset->swidth[0]; s++){
      sw = hset->swidth[s];
      if (hset->tmRecs[s].mixId == NULL){
         HRError(7030,"CheckTMRecs: no mix id set in stream %d",s);
         return(FAIL);
      }
      if (hset->tmRecs[s].mixes == NULL){
         HRError(7030,"CheckTMRecs: no mixes array allocated in stream %d",s);
         return(FAIL);
      }
      for (m=1; m<=hset->tmRecs[s].nMix; m++){
         mp = hset->tmRecs[s].mixes[m];
         if(CheckMix("TMRec",mp,0,s,m,sw)<SUCCESS)
            return(FAIL);
      }
   }
   return(SUCCESS);
}

/* CheckDiscrete: check discrete HMM Set */
static ReturnStatus CheckDiscrete(HMMSet *hset)
{
   int s;

   for (s=1;s<=hset->swidth[0]; s++){
      if (hset->swidth[s] != 1) {
         HRError(7030,"CheckDiscrete: stream width not equal to 1 in discrete stream %d ",s);
         return(FAIL);
      }
   }
   return (SUCCESS);
}

/* CheckHSet: check the consistency of a complete HMM Set */
static ReturnStatus CheckHSet(HMMSet *hset)
{
   int h;
   MLink m;
   
   for (h=0; h<MACHASHSIZE; h++)
      for (m=hset->mtab[h]; m!=NULL; m=m->next)
         if (m->type == 'h')
            if(CheckHMM(m->id->name,(HLink)m->structure)<SUCCESS)
               return(FAIL);
   if ((hset->hsKind == TIEDHS) && (CheckTMRecs(hset)<SUCCESS))
      return(FAIL);
   if ((hset->hsKind == DISCRETEHS) && (CheckDiscrete(hset)<SUCCESS))
      return(FAIL);

   ClearSeenFlags(hset,CLR_ALL);
   return(SUCCESS);

}

/* ------------------------ Lexical Scanner ------------------------ */

#define MAXSYMLEN 40

/* Internal and binary keyword representation */
typedef enum {   /* Only a character big !! */
   BEGINHMM, USEMAC, ENDHMM, NUMMIXES, 
   NUMSTATES, STREAMINFO, VECSIZE, 
   NDUR, PDUR, GDUR, RELDUR, GENDUR,
   DIAGCOV,  FULLCOV, XFORMCOV,
   STATE, TMIX, MIXTURE, STREAM, SWEIGHTS,
   MEAN, VARIANCE, INVCOVAR, XFORM, GCONST,
   DURATION, INVDIAGCOV, TRANSP, DPROB, LLTCOV, LLTCOVAR,
   XFORMKIND=90, PARENTXFORM, NUMXFORMS, XFORMSET,
   LINXFORM, OFFSET, BIAS, BLOCKINFO, BLOCK, BASECLASS, 
   CLASS, XFORMWGTSET, CLASSXFORM, MMFIDMASK, PARAMETERS,
   NUMCLASSES, ADAPTKIND, PREQUAL, INPUTXFORM,
   RCLASS=110, REGTREE, NODE, TNODE,
   HMMSETID=119,
   PARMKIND=120, 
   MACRO, EOFSYM, NULLSYM   /* Special Syms - not literals */
} Symbol;

/* Mapping between verbose keywords for readable external HMM format */
/*  and internal/external binary token numbers */
static struct {
   char *name;
   Symbol sym; } symMap[] = 
   { { "BEGINHMM", BEGINHMM }, { "USE", USEMAC }, 
     { "ENDHMM", ENDHMM }, { "NUMMIXES", NUMMIXES }, 
     { "NUMSTATES", NUMSTATES }, { "STREAMINFO", STREAMINFO }, 
     { "VECSIZE", VECSIZE }, { "NULLD", NDUR }, 
     { "POISSOND", PDUR }, { "GAMMAD", GDUR }, 
     { "RELD", RELDUR }, { "GEND", GENDUR }, 
     { "DIAGC", DIAGCOV }, {  "FULLC", FULLCOV }, 
     { "XFORMC", XFORMCOV }, { "STATE", STATE }, 
     { "TMIX", TMIX }, { "MIXTURE", MIXTURE }, 
     { "STREAM", STREAM }, { "SWEIGHTS", SWEIGHTS }, 
     { "MEAN", MEAN }, { "VARIANCE", VARIANCE }, 
     { "INVCOVAR", INVCOVAR }, { "XFORM", XFORM }, 
     { "GCONST" , GCONST }, { "DURATION", DURATION }, 
     { "INVDIAGC", INVDIAGCOV }, { "TRANSP", TRANSP }, 
     { "DPROB", DPROB }, { "LLTC", LLTCOV }, 
     { "LLTCOVAR", LLTCOVAR }, 
     { "RCLASS", RCLASS }, { "REGTREE", REGTREE }, 
     { "NODE", NODE }, { "TNODE", TNODE }, 
     { "HMMSETID", HMMSETID }, 
     { "PARMKIND", PARMKIND }, { "MACRO", MACRO }, 
     { "EOF", EOFSYM }, 
     /* Transformation symbols */
     {"XFORMKIND", XFORMKIND } , {"PARENTXFORM", PARENTXFORM },
     {"NUMXFORMS", NUMXFORMS }, {"XFORMSET", XFORMSET },
     {"LINXFORM", LINXFORM }, {"OFFSET", OFFSET }, 
     {"BIAS", BIAS }, {"BLOCKINFO", BLOCKINFO }, {"BLOCK", BLOCK }, 
     {"BASECLASS", BASECLASS }, {"CLASS", CLASS }, 
     {"XFORMWGTSET", XFORMWGTSET }, {"CLASSXFORM", CLASSXFORM }, 
     {"MMFIDMASK", MMFIDMASK }, {"PARAMETERS", PARAMETERS },
     {"NUMCLASSES", NUMCLASSES }, {"ADAPTKIND", ADAPTKIND },
     {"PREQUAL", PREQUAL }, {"INPUTXFORM", INPUTXFORM },
     { "", NULLSYM } };

/* Reverse lookup table for above */
static char *symNames[NULLSYM+1];

#define NUMSYM (sizeof(symMap)/sizeof(symMap[0]))

typedef struct {     /* Token returned by the scanner */
   Symbol sym;          /* the current input symbol */
   Boolean binForm;     /* binary form of keyword symbol */
   ParmKind pkind;      /* samp kind when sym==PARMKIND */
   char macroType;      /* current macro type if sym==MACRO */
} Token;
   
void InitSymNames(void)
{
   int i;
   for (i=0;i<=NULLSYM;i++) symNames[i]="";
   for (i=0;i<NUMSYM;i++) symNames[symMap[i].sym]=symMap[i].name;
}

/* InitScanner: initialise scanner for new source */
ReturnStatus InitScanner(char *fname, Source *src, Token *tok, HMMSet *hset)
{
   if(InitSource(fname, src, HMMDefFilter)<SUCCESS){
      return(FAIL);
   }
   tok->sym = NULLSYM; tok->macroType = ' '; 
   tok->binForm = FALSE;
   return(SUCCESS);
}

/* TermScanner: terminate scanner for given source */
static void TermScanner(Source *src)
{
   CloseSource(src);
}

/* HMError: report a HMM definition error */
static void HMError(Source *src, char *message)
{
   char buf[MAXSTRLEN];
   
   fflush(stdout);
   fprintf(stderr,"HMM Def Error: %s at %s\n",
           message,SrcPosition(*src,buf));
   HRError(7050,"HMError:");
}

/* GetToken: put next symbol from given source into token */
static ReturnStatus GetToken(Source *src, Token *tok)
{
   char buf[MAXSYMLEN],tmp[MAXSTRLEN];
   int i,c,imax,sym;
   
   tok->binForm = FALSE;
   while (isspace(c=GetCh(src)));     /* Look for symbol or Macro */
   if (c != '<' && c != ':' && c != '~'  && c != '.' && c != '#') {
      if (c == EOF) {
         if (trace&T_TOK) printf("HModel:   tok=<EOF>\n");
         tok->sym=EOFSYM; return(SUCCESS);
      }
      HMError(src,"GetToken: Symbol expected");
      return(FAIL);
   }
   if (c == '~'){                    /* If macro sym return immediately */
      c = tolower(GetCh(src));
      if (c!='s' && c!='m' && c!='u' && c!='x' && c!='d' && c!='c' &&
          c!='r' && c!='a' && c!='b' && c!='g' && c!='f' && c!='y' && c!='j' &&
          c!='v' && c!='i' && c!='t' && c!='w' && c!='h' && c!='o')
         {
            HMError(src,"GetToken: Illegal macro type");
            return(FAIL);
         }
      tok->macroType = c; tok->sym = MACRO;
      if (trace&T_TOK) printf("HModel:   MACRO ~%c\n",c);
      return(SUCCESS);
   }
   i=0; imax = MAXSYMLEN-1;
   if (c=='#') {           /* if V1 mmf header convert to ~h */
      while ((c=GetCh(src)) != '#' && i<imax)
         buf[i++] = c;
      buf[i] = '\0';
      if (strcmp(buf,"!MMF!") != 0){
         HMError(src,"GetToken: expecting V1 style MMF header #!MMF!#");
         return(FAIL);
      }
      tok->sym = MACRO; tok->macroType = 'h';
      if (trace&T_TOK) printf("HModel:   MACRO ~h (#!MMF!#)\n");
      return(SUCCESS);
   }
   if (c=='.'){            /* if . and not EOF convert to ~h */
      while (isspace(c=GetCh(src)));
      if (c == EOF) {
         if (trace&T_TOK) printf("HModel:   tok=.<EOF>\n");
         tok->sym=EOFSYM;
         return(SUCCESS);
      }
      UnGetCh(c,src);
      tok->sym = MACRO; tok->macroType = 'h';
      if (trace&T_TOK) printf("HModel:   MACRO ~h (.)\n");
      return(SUCCESS);     
   }  
   if (c=='<') {                 /* Read verbose symbol string into buf */
      while ((c=GetCh(src)) != '>' && i<imax)
         buf[i++] = islower(c)?toupper(c):c;
      buf[i] = '\0';
      if (c != '>'){
         HMError(src,"GetToken: > missing in symbol");
         return(FAIL);
      }
      for (sym=0; sym<NUMSYM; sym++) /* Look symbol up in symMap */
         if (strcmp(symMap[sym].name,buf) == 0) {
            tok->sym = symMap[sym].sym;
            if (trace&T_TOK) printf("HModel:   tok=<%s>\n",buf);
            return(SUCCESS);                                /* and return */  
         }
   } else {
      /* Read binary symbol into buf */
      tok->binForm = TRUE;
      sym = GetCh(src);
      if (sym>=BEGINHMM && sym<PARMKIND) {
         if (trace&T_TOK) printf("HModel:   tok=:%s\n",symNames[sym]);
         tok->sym = (Symbol) sym;
         return(SUCCESS);                                /* and return */  
      }     
   }
         
   /* if symbol not in symMap then it may be a sampkind */
   if ((tok->pkind = Str2ParmKind(buf)) != ANON){
      tok->sym = PARMKIND;
      if (trace&T_TOK) printf("HModel:   tok=SK[%s]\n",buf);
      return(SUCCESS);
   }
   strcpy(tmp,"GetToken: Unknown symbol ");
   HMError(src,strcat(tmp,buf));
   return(FAIL);
}

/*
  Looks for the macroname in the current directory. If it is 
  notfound search through the input transform directories
  in order. FOpen is used to allow optional generation of
  Error messages using T_XFD
*/
static char *InitXFormScanner(HMMSet *hset, char *macroname, char *fname,
			      Source *src, Token *tok)
{
   static char buf[MAXSTRLEN];
   XFDirLink p;
   Boolean isPipe;
   FILE *f;

   if ((fname==NULL) || ((f=FOpen(fname,NoFilter,&isPipe)) == NULL)) {
      if ((trace&T_XFD) && (fname!=NULL))
         HRError(7010,"InitXFormScanner: Cannot open source file %s",fname);
      p = xformDirNames;
      while ((p!=NULL) && 
             ((f=FOpen(MakeFN(macroname,p->dirName,NULL,buf),NoFilter,&isPipe)) == NULL)) {
         if (trace&T_XFD) 
            HRError(7010,"InitXFormScanner: Cannot open source file %s",buf);
         p = p->next;
      }
      if (p==NULL) { /* Failed to find macroname */
         HError(7035,"Failed to find macroname %s",macroname);
      } else { /* Close file and initialise scanner */
         FClose(f,isPipe);
         InitScanner(buf,src,tok,hset);
      }
      if (trace&T_TOP) printf("Loading macro file %s\n",buf);
      return buf;
   } else {
      FClose(f,isPipe);
      if (trace&T_TOP) printf("Loading macro file %s\n",fname);
      InitScanner(fname,src,tok,hset);
      return fname;
   }
}

/* ------------------- HMM 'option' handling ----------------------- */


static void OWarn(HMMSet *hset,Boolean equal,char *opt)
{
   if (!equal && hset->optSet)
      HRError(-7032,"OWarn: change HMM Set %s",opt);
}

/* GetOption: read a HMM option specifier - value set in nState pointer */
static ReturnStatus GetOption(HMMSet *hset, Source *src, Token *tok, int *nState)
{
   DurKind dk;
   char buf[MAXSTRLEN];
   short vs,sw[SMAX],nSt=0;
   int i;
   Boolean ntok=TRUE;

   static InputXForm* GetInputXForm(HMMSet *hset, Source *src, Token *tok);

   switch (tok->sym) {
   case NUMSTATES:
      if (!ReadShort(src,&nSt,1,tok->binForm)){
         HMError(src,"NumStates Expected");
         return(FAIL);
      }
      *nState=nSt;
      break;
   case PARMKIND:    
      OWarn(hset,hset->pkind==tok->pkind,"parmKind");
      hset->pkind = tok->pkind;
      break;
   case NDUR:
   case PDUR:
   case GDUR:
   case RELDUR:
   case GENDUR:
      dk = (DurKind) (NULLD + (tok->sym-NDUR));
      OWarn(hset,hset->dkind==dk,"durKind");
      hset->dkind = dk; 
      break;
   case HMMSETID:
      if (!ReadString(src,buf)){
         HMError(src,"HMM identifier expected");
         return(FAIL);
      }
      hset->hmmSetId=CopyString(hset->hmem,buf);
      SkipWhiteSpace(src);
      break;
   case INPUTXFORM:
      if(GetToken(src,tok)<SUCCESS){
         HMError(src,"GetOption: GetToken failed");
         return(FAIL);     
      }
      if (tok->sym==MACRO && tok->macroType=='j') {
         if (!ReadString(src,buf))
            HError(7013,"GetOption: cannot read input xform macro name");
         hset->xf = LoadInputXForm(hset,buf,NULL);	
      } else {
         hset->xf = GetInputXForm(hset,src,tok);
         hset->xf->xformName = CopyString(hset->hmem,src->name);
         ntok = FALSE;
      }     
      break;
   case VECSIZE:
      if (!ReadShort(src,&vs,1,tok->binForm)){
         HMError(src,"Vector Size Expected");
         return(FAIL);
      }
      OWarn(hset,hset->vecSize==vs,"vecSize");
      hset->vecSize = vs;
      break;
   case STREAMINFO:
      if (!ReadShort(src,sw,1,tok->binForm)){
         HMError(src,"Num Streams Expected");
         return(FAIL);
      }
      if (sw[0] >= SMAX){
         HMError(src,"Stream limit exceeded");
         return(FAIL);
      }
      if (!ReadShort(src,sw+1,sw[0],tok->binForm)){
         HMError(src,"Stream Widths Expected");
         return(FAIL);
      }
      OWarn(hset,hset->swidth[0]==sw[0],"swidth[0]");       
      for (i=0; i<=sw[0]; i++) hset->swidth[i] = sw[i];
      break;
   case DIAGCOV:
      OWarn(hset,hset->ckind==DIAGC,"covKind");
      hset->ckind = DIAGC; 
      break;
   case FULLCOV:
      OWarn(hset,hset->ckind==FULLC,"covKind");
      hset->ckind = FULLC; 
      break;
   case XFORMCOV:
      OWarn(hset,hset->ckind==XFORMC,"covKind");
      hset->ckind = XFORMC; 
      break;
   case INVDIAGCOV:
      OWarn(hset,hset->ckind==INVDIAGC,"covKind");
      hset->ckind = INVDIAGC; 
      break;
   case LLTCOV:
      OWarn(hset,hset->ckind==LLTC,"covKind");
      hset->ckind = LLTC; 
      break;

   default: 
      HMError(src,"GetOption: Ilegal Option Symbol");
      return(FAIL);
   }
   if (ntok && (GetToken(src,tok)<SUCCESS)){
      HMError(src,"GetOption: GetToken failed");
      return(FAIL);     
   }
   
   return(SUCCESS);
}

/* FreezeOptions: freeze the global options in HMM set */
static ReturnStatus FreezeOptions(HMMSet *hset)
{
   int i;
   
   if (hset->optSet) return(SUCCESS);
   if (hset->vecSize == 0) {
      if (hset->swidth[0] > 0 && hset->swidth[1] > 0)
         for (i=1; i<=hset->swidth[0]; i++)
            hset->vecSize += hset->swidth[i];
      else{
         HRError(7032,"FreezeOptions: vecSize not set");
         return(FAIL);
      }
   }
   if (hset->swidth[0] == 0) {
      hset->swidth[0] = 1; hset->swidth[1] = hset->vecSize;
   }
   if (hset->pkind == 0){
      HRError(7032,"FreezeOptions: parmKind not set");
      return(FAIL);
   }
   hset->optSet = TRUE;
   return(SUCCESS);
}

/* CheckOptions: check that options are set in given HMM set */
static ReturnStatus CheckOptions(HMMSet *hset)
{
   if (!hset->optSet){
      HRError(7032,"CheckOptions: options not set in HMM Set");
      return(FAIL);
   }
   return(SUCCESS);
}

/* ---------------------- Input XForm Directory Handling ---------------------- */

/* EXPORT->AddInXFormDir: Add given file name to set */
/* Doesn't check for repeated specification! */
void AddInXFormDir(HMMSet *hset, char *dirname)
{
  XFDirLink p,q;
   
  p = (XFDirLink)New(hset->hmem,sizeof(XFDirInfo));
  p->next = NULL;
  p->dirName = CopyString(hset->hmem,dirname);
  if (xformDirNames == NULL)
    xformDirNames = p;
  else {  /* store in order of arrival */
    for (q=xformDirNames; q->next != NULL; q=q->next);
    q->next = p;
  }
}

/* ---------------------- MMF File Name Handling ---------------------- */

/* FindMMF: find given file name in HMM set */
static MILink FindMMF(HMMSet *hset, char *fname, Boolean ignorePath)
{
   MILink p;
   char buf1[MAXSTRLEN],buf2[MAXSTRLEN];
   
   for (p=hset->mmfNames; p!=NULL; p=p->next){
      if (ignorePath){
         if (strcmp(NameOf(fname,buf1),NameOf(p->fName,buf2)) == 0 ) 
            return p;
      }else{
         if (strcmp(fname,p->fName) == 0 )
            return p;
      }
   }
   return NULL;
}

/* EXPORT->AddMMF: Add given file name to set */
MILink AddMMF(HMMSet *hset, char *fname)
{
   MILink p,q;
   
   if ((p = FindMMF(hset,fname,FALSE)) != NULL) 
      return(p);
   p = (MILink)New(hset->hmem,sizeof(MMFInfo));
   ++hset->numFiles;
   p->isLoaded = FALSE; p->next = NULL;
   p->fName = CopyString(hset->hmem,fname);
   p->fidx = hset->numFiles;
   if (hset->mmfNames == NULL)
      hset->mmfNames = p;
   else {  /* store in order of arrival */
      for (q=hset->mmfNames; q->next != NULL; q=q->next);
      q->next = p;
   }
   return p;
}


/* -------------- Special Discrete/Tied Mixture Input Routines ------------ */

/* 
   Run-length encoding for verbose external form uses the notation w*n
   to mean that w is repeated n times.  For binary form, all weights
   are stored as unsigned shorts.  The msb is set to indicate that a
   repeat count follows.  Repeat counts are stored as unsigned chars.
*/

/* GetTiedWeights: parse src and get compact mixture weight def */
ReturnStatus GetTiedWeights(Source *src, Token *tok, int M, Vector tpdf)
{
   float weight=0.0;
   short repCount=0;      /* repeat counter for weights */
   int c,m;
   
   if (trace&T_PAR) printf("HModel: GetTiedWeights: M=%d\n",M);
   for (m=1; m<=M; m++) {
      if (repCount>0)         /* get mixture weight */
         --repCount;
      else {
         if (tok->binForm) {
            if (!ReadFloat(src,&weight,1,TRUE)){
               HMError(src,"Tied Weight expected");
               return(FAIL);
            }
            if (weight<0.0) {
               repCount=GetCh(src);
               --repCount; weight = weight+2.0;
            }
         } else {
            if (!ReadFloat(src,&weight,1,FALSE)){
               HMError(src,"Discrete Weight expected");
               return(FAIL);
            }
            c=GetCh(src);
            if (c == '*') {
               if (!ReadShort(src,&repCount,1,FALSE)){
                  HMError(src,"Discrete Repeat Count expected");
                  return(FAIL);
               }
               --repCount;
            } else
               UnGetCh(c,src);
         }
      }
      tpdf[m] = weight;      /* set it */
   }
   if(GetToken(src,tok)<SUCCESS){
      HMError(src,"GetToken failed");
      return(FAIL);
   }

   return(SUCCESS);
}

/* GetDiscreteWeights: parse src and get compact mixture weight def */
ReturnStatus GetDiscreteWeights(Source *src, Token *tok, int M, ShortVec dpdf)
{
   short weight=0;
   short repCount=0;      /* repeat counter for weights */
   int c,m;
   
   if (trace&T_PAR) printf("HModel: GetDiscreteWeights: M=%d\n",M);
   for (m=1; m<=M; m++) {
      if (repCount>0)         /* get mixture weight */
         --repCount;
      else {
         if (tok->binForm) {
            if (!ReadShort(src,&weight,1,TRUE)){
               HMError(src,"Discrete Weight expected");
               return(FAIL);
            }
            if (weight<0) {
               repCount=GetCh(src);
               --repCount; weight &= 077777;
            }
         } else {
            if (!ReadShort(src,&weight,1,FALSE)){
               HMError(src,"Discrete Weight expected");
               return(FAIL);
            }
            c=GetCh(src);
            if (c == '*') {
               if (!ReadShort(src,&repCount,1,FALSE)){
                  HMError(src,"Discrete Repeat Count expected");
                  return(SUCCESS);
               }
               --repCount;
            } else
               UnGetCh(c,src);
         }
      }
      dpdf[m] = weight;      /* set it */
   }
   if(GetToken(src,tok)<SUCCESS){
      HMError(src,"GetToken failed");
      return(FAIL);
   }
   
   return(SUCCESS);
}

/* InitTMixRecs: create all TMixRecs in given hset */
void InitTMixRecs(HMMSet *hset, int s, int M)
{
   TMixRec *p;
   MixPDF **mpp;
   TMProb *tm;
   
   p = hset->tmRecs+s;
   p->mixId = NULL; p->nMix = M;
   mpp = (MixPDF **)New(hset->hmem,sizeof(MixPDF *) * M);
   p->mixes = mpp-1;
   tm = (TMProb *)New(hset->hmem,sizeof(TMProb)*M);
   p->probs = tm-1;
}

/* GetTiedMixtures: parse src and get compact tied mixture def */
ReturnStatus GetTiedMixtures(HMMSet *hset, Source *src, Token *tok, 
                             int M, int s, Vector tpdf)
{
   char tmName[MAXSTRLEN],macName[2*MAXSTRLEN],intstr[20];
   LabId id,mid;
   Boolean isNew;
   int m;
   MLink q;
   
   if (trace&T_PAR) printf("HModel: GetTiedMixtures\n");
   if (!ReadString(src,tmName)){      /* read generic ~m macro name */
      HMError(src,"Tied Mix macro name expected");
      return(FAIL);
   }
   id = GetLabId(tmName,TRUE);
   isNew = hset->tmRecs[s].mixId == NULL;
   if (isNew) {
      InitTMixRecs(hset,s,M);
      hset->tmRecs[s].mixId = id;
      for (m=1; m<=M; m++){
         sprintf(intstr,"%d",m);
         strcpy(macName,tmName);
         strcat(macName,intstr);
         if((mid = GetLabId(macName,FALSE)) == NULL){
            HRError(7035,"GetTiedMixtures: Unknown tied mix macro name %s",macName);
            return(FAIL);
         }
         if ((q = FindMacroName(hset,'m',mid))==NULL){
            HRError(7035,"GetTiedMixtures: no macro %s in this set",macName);
            return(FAIL);
         }
         hset->tmRecs[s].mixes[m] = (MixPDF *)q->structure;
      }
   }else {
      if (hset->tmRecs[s].mixId != id){
         HMError(src,"Bad Generic ~m Macro Name in TMix");
         return(FAIL);
      }
      if (hset->tmRecs[s].nMix != M){
         HMError(src,"Inconsistent Num Mixtures in TMix");
         return(FAIL);
      }
   }
   if(GetTiedWeights(src,tok,M,tpdf)<SUCCESS){
      HMError(src, "GetTiedWeights failed");
      return(FAIL);
   }
   return(SUCCESS);
}

/* ------------ Standard Macro Definition Input Routines ----------------- */

/* GetOptions: read a global options macro, return numStates if set */
static ReturnStatus GetOptions(HMMSet *hset, Source *src, Token *tok, int *nState)
{
   int p=0;
   
   *nState=0; 

   if (trace&T_PAR) printf("HModel: GetOptions\n");
   if(GetToken(src,tok)<SUCCESS){
      HMError(src,"GetOptions: GetToken failed");
      return(FAIL);
   }
   while (tok->sym == PARMKIND || tok->sym == INVDIAGCOV || 
          tok->sym == HMMSETID  || tok->sym == INPUTXFORM ||
          (tok->sym >= NUMSTATES && tok->sym <= XFORMCOV)){
      if(GetOption(hset,src,tok,&p)<SUCCESS){
         HMError(src,"GetOptions: GetOption failed");
         return(FAIL);
      }
      if (p>*nState) *nState = p;
   }
   FreezeOptions(hset);
   return(SUCCESS);
}

/* GetStructure: next input token is a string containing macro name,
                 a pointer to corresponding structure is returned */
static Ptr GetStructure(HMMSet *hset, Source *src, char type)
{
   char buf[MAXSTRLEN];
   LabId id;
   MLink m;

   if (!ReadString(src,buf)){
      HRError(7013,"GetStructure: cannot read macro name");
      return(NULL);
   }
   id = GetLabId(buf,FALSE);
   if (id==NULL){
      HRError(7035,"GetStructure: undef macro name %s, type %c",buf,type);
      return(NULL);
   }
   m = FindMacroName(hset,type,id);
   if (m==NULL){
      HRError(7035,"GetStructure: no macro %s, type %c exists",buf,type);  
      return(NULL);
   }
   if (trace&T_MAC)
      printf("HModel: getting structure ~%c %s -> %p\n",
             type,buf,m->structure);
   return m->structure;
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

static RegNode *CreateNodeInfo(MemHeap *m, short nodeId)
{
   RegNode *n;
  
   n  = (RegNode *) New(m, sizeof(RegNode));
   n->nodeIndex = nodeId;
   n->nodeComps = 0;
   n->nodeOcc = 0.0;
   n->nBases = 0;
   n->WTrans = NULL;
   n->HTrans = NULL;
   n->bases  = NULL;
   n->backTrans = NULL;

   return n;
}

/* GetRegTree: parse src and return the regression tree structure */
static RegTree *GetRegTree(HMMSet *hset, Source *src, Token *tok)
{
   RegTree *rtree = NULL, *t;
   short size, i;
   int comps;
   short index;

   /* allocate space for the regression tree root node */
   rtree = (RegTree *) New(hset->hmem, sizeof(RegTree));
   rtree->left = rtree->right = NULL;
   rtree->nodeInfo = CreateNodeInfo(hset->hmem, 1);

   if (trace&T_PAR) printf("HModel: GetRegTree\n");
   if (tok->sym==REGTREE) {      
      if (!ReadShort(src,&size,1,tok->binForm)){
         HMError(src,"Size of the regression tree is expected");
         return(NULL);
      }
      size *= 2; size -= 1;
      for (i = 1; i <= size; i++) {
      
         if(GetToken(src,tok)<SUCCESS){
            HMError(src,"GetToken failed");
            return(NULL);
         }
         if (!ReadShort(src,&index,1,tok->binForm)){
            HMError(src,"Parent node index for regression tree expected");
            return(NULL);
         }
         t = FindNode(rtree, NULL, index);
         if (t == NULL){
            HRError(7085, "GetRegTree: Can't find node %d in tree", index);
            return(NULL);
         }
         switch(tok->sym) {
         case NODE:
            /* create children */
            t->left  = (RegTree *) New(hset->hmem, sizeof(RegTree));
            t->right = (RegTree *) New(hset->hmem, sizeof(RegTree));
            t->left->left = t->left->right = NULL;
            t->right->left = t->right->right = NULL;
            if (!ReadShort(src,&index,1,tok->binForm)){
               HMError(src,"Left node index for regression tree expected");
               return(NULL);
            }
            t->left->nodeInfo  = CreateNodeInfo(hset->hmem, index);
            if (!ReadShort(src,&index,1,tok->binForm)){
               HMError(src,"Right node index for regression tree expected");
               return(NULL);
            }
            t->right->nodeInfo = CreateNodeInfo(hset->hmem, index);
            break;
         case TNODE:
            /* load the number of components */
            if (!ReadInt(src,&comps,1,tok->binForm)){
               HMError(src,"Number of components for regression base class expected");
               return(NULL);
            }
            t->nodeInfo->nodeComps = comps;
            /* if (!ReadFloat(src,&tmp,1,tok->binForm))
               HMError(src,"Occ count (training) for regression base class expected");
               if (!ReadFloat(src,&tmp,1,tok->binForm))
               HMError(src,"Node score for regression base class expected"); */
            break;
         default:
            HRError(7085,"GetRegTree:Unexpected token symbol");
            return(NULL);
         }
      }
   }

   if(GetToken(src,tok)<SUCCESS){
      HMError(src,"GetToken failed");
      return(NULL);
   }
  
   return rtree; 
}


/* GetMean: parse src and return Mean structure */
static SVector GetMean(HMMSet *hset, Source *src, Token *tok)
{
   SVector m = NULL;
   short size;
   
   if (trace&T_PAR) printf("HModel: GetMean\n");
   if (tok->sym==MEAN) {      
      if (!ReadShort(src,&size,1,tok->binForm)){
         HMError(src,"Size of Mean Vector expected");
         return(NULL);
      }
      m = CreateSVector(hset->hmem,size);
      if (!ReadVector(src,m,tok->binForm)){
         HMError(src,"Mean Vector expected");
         return(NULL);
      }
   }
   
   else if (tok->sym==MACRO && tok->macroType=='u'){
      if((m=(SVector)GetStructure(hset,src,'u'))==NULL){
         HMError(src,"GetStructure Failed");
         return(NULL);
      }
      IncUse(m);
   } else{
      HMError(src,"<Mean> symbol expected in GetMean");
      return(NULL);
   }
   if(GetToken(src,tok)<SUCCESS){
      HMError(src,"GetToken failed");
      return(NULL);
   }
   return m;
}

/* GetVariance: parse src and return Variance structure */
static SVector GetVariance(HMMSet *hset, Source *src, Token *tok)
{
   SVector v = NULL;
   short size; 
   
   if (trace&T_PAR) printf("HModel: GetVariance\n");
   if (tok->sym==VARIANCE) {
      if (!ReadShort(src,&size,1,tok->binForm)){
         HMError(src,"Size of Variance Vector expected");
         return(NULL);
      }
      v = CreateSVector(hset->hmem,size);
      if (!ReadVector(src,v,tok->binForm)){
         HMError(src,"Variance Vector expected");
         return(NULL);
      }
   }
   
   else if (tok->sym==MACRO && tok->macroType=='v'){
      if((v=(SVector)GetStructure(hset,src,'v'))==NULL){
         HMError(src,"GetStructure Failed");
         return(NULL);
      }
      IncUse(v);
   } else{
      HMError(src,"<Variance> symbol expected in GetVariance");
      return(NULL);
   }
   if(GetToken(src,tok)<SUCCESS){
      HMError(src,"GetToken failed");
      return(NULL);
   }
   return v;
}

/* GetCovar: parse src and return Covariance structure */
static STriMat GetCovar(HMMSet *hset, Source *src, Token *tok)
{
   STriMat m = NULL;
   short swidth;
   
   if (trace&T_PAR) printf("HModel: GetCovar\n");
   if (tok->sym==INVCOVAR || tok->sym==LLTCOVAR) {
      if (!ReadShort(src,&swidth,1,tok->binForm)){
         HMError(src,"Size of Inv Covariance expected");
         return(NULL);
      }
      m = CreateSTriMat(hset->hmem,swidth);
      if (!ReadTriMat(src,m,tok->binForm)){
         HMError(src,"Inverse/LLT Covariance Matrix expected");
         return(NULL);
      }
   } else if (tok->sym==MACRO && 
              (tok->macroType=='i' || tok->macroType=='c') ){
      if((m=(STriMat)GetStructure(hset,src,tok->macroType))==NULL){
         HMError(src,"GetStructure Failed");
         return(NULL);
      }
      IncUse(m);
   } else{
      HMError(src,"<InvCovar>/<LLTCovar> symbol expected in GetCovar");
      return(NULL);
   }
   if(GetToken(src,tok)<SUCCESS){
      HMError(src,"GetToken failed");
      return(NULL);
   }
   return m;
}

/* GetTransform: parse src and return Transform structure */
static SMatrix GetTransform(HMMSet *hset, Source *src, Token *tok)
{
   SMatrix m = NULL;
   MemHeap *hmem;
   short xformRows,xformCols;
   
   if (trace&T_PAR) printf("HModel: GetTransform\n");
   if (tok->sym==XFORM) {
      if (hset==NULL) hmem = &xformStack;
      else hmem = hset->hmem;
      if (!ReadShort(src,&xformRows,1,tok->binForm)){
         HMError(src,"Num Rows in Xform matrix expected");
         return(NULL);
      }
      if (!ReadShort(src,&xformCols,1,tok->binForm)){
         HMError(src,"Num Cols in Xform matrix expected");
         return(NULL);
      }
      m = CreateSMatrix(hmem,xformRows,xformCols);
      if (!ReadMatrix(src,m,tok->binForm)){
         HMError(src,"Transform Matrix expected");
         return(NULL);
      }
   }  else  if ((tok->sym==MACRO && tok->macroType=='x') && (hset != NULL)) {
      if((m=(SMatrix)GetStructure(hset,src,'x'))==NULL){
         HMError(src,"GetStructure Failed");
         return(NULL);
      }
      IncUse(m);
   } else{
      HMError(src,"<Xform> symbol expected in GetTransform");
      return(NULL);
   }
   if(GetToken(src,tok)<SUCCESS){
      HMError(src,"GetToken failed");
      return(NULL);
   }
   return m;
}

/* GetDuration: parse src and return Duration structure */
static SVector GetDuration(HMMSet *hset, Source *src, Token *tok)
{
   SVector v = NULL;
   short size;
   
   if (trace&T_PAR) printf("HModel: GetDuration\n");
   if (tok->sym==DURATION) {
      if (!ReadShort(src,&size,1,tok->binForm)){
         HMError(src,"Size of Duration Vector expected");
         return(NULL);
      }
      v = CreateSVector(hset->hmem,size);
      if (!ReadVector(src,v,tok->binForm)){
         HMError(src,"Duration Vector expected");
         return(NULL);
      }
   } else  if (tok->sym==MACRO && tok->macroType=='d'){
      if((v=(SVector)GetStructure(hset,src,'d'))==NULL){
         HMError(src,"GetStructure Failed");
         return(NULL);
      }
      IncUse(v);
   } else{
      HMError(src,"<Duration> symbol expected in GetDuration");
      return(NULL);
   }
   if(GetToken(src,tok)<SUCCESS){
      HMError(src,"GetToken failed");
      return(NULL);
   }
   return v;
}

/* GetSWeights: parse src and return vector of stream weights */
static SVector GetSWeights(HMMSet *hset, Source *src, Token *tok)
{
   SVector v = NULL;
   short size;
   
   if (trace&T_PAR) printf("HModel: GetSWeights\n");
   if (tok->sym==SWEIGHTS) {
      if (!ReadShort(src,&size,1,tok->binForm)){
         HMError(src,"Num stream weights expected");
         return(NULL);
      }
      v = CreateSVector(hset->hmem,size);
      if (!ReadVector(src,v,tok->binForm)){
         HMError(src,"Stream Weights expected");
         return(NULL);
      }
   } else  if (tok->sym==MACRO && tok->macroType=='w'){
      if((v=(SVector)GetStructure(hset,src,'w'))==NULL){
         HMError(src,"GetStructure Failed");
         return(NULL);
      }
      IncUse(v);
   } else{
      HMError(src,"<SWeights> symbol expected in GetSWeights");
      return(NULL);
   }
   if(GetToken(src,tok)<SUCCESS){
      HMError(src,"GetToken failed");
      return(NULL);
   }
   return v;
}

/* GetMixPDF: parse src and return MixPDF structure */
static MixPDF *GetMixPDF(HMMSet *hset, Source *src, Token *tok)
{
   MixPDF *mp;
  
   if (trace&T_PAR) printf("HModel: GetMixPDF\n");
   if (tok->sym==MACRO && tok->macroType=='m') {
      if ((mp = (MixPDF *)GetStructure(hset,src,'m'))==NULL){
         HMError(src,"GetStructure Failed");
         return(NULL);
      }
      ++mp->nUse;
      if(GetToken(src,tok)<SUCCESS){
         HMError(src,"GetToken failed");
         return(NULL);
      }
   } else {
      mp = (MixPDF *)New(hset->hmem,sizeof(MixPDF));
      mp->nUse = 0; mp->hook = NULL; mp->gConst = LZERO;
      mp->rClass = 0;
      if (tok->sym == RCLASS) {
         short r;
         if (!ReadShort(src,&r,1,tok->binForm)){
            HMError(src,"Regression Class Number expected");
            return(NULL);
         }
         mp->rClass = r;
         if(GetToken(src,tok)<SUCCESS){
            HMError(src,"GetToken failed");
            return(NULL);
         }
      }
      if((mp->mean = GetMean(hset,src,tok))==NULL){      
         HMError(src,"GetMean Failed");
         return(NULL);
      }
      if (tok->sym==VARIANCE || (tok->sym==MACRO && tok->macroType=='v')) {
         if((mp->cov.var = GetVariance(hset,src,tok))==NULL){
            HMError(src,"GetVariance Failed");
            return(NULL);
         }
         if (hset->ckind == DIAGC || hset->ckind == NULLC)
            mp->ckind = DIAGC;
         else{
            HRError(7032,"GetMixPDF: trying to change global cov type to DiagC");
            return(NULL);
         }
      } else if (tok->sym==INVCOVAR || (tok->sym==MACRO && tok->macroType=='i')){
         if((mp->cov.inv = GetCovar(hset,src,tok))==NULL){
            HMError(src,"GetCovar Failed");
            return(NULL);
         }
         if (hset->ckind == FULLC || hset->ckind == NULLC)
            mp->ckind = FULLC;
         else{
            HRError(7032,"GetMixPDF: trying to change global cov type to FullC");
            return(NULL);
         }
      } else if (tok->sym==LLTCOVAR || (tok->sym==MACRO && tok->macroType=='c')){
         if((mp->cov.inv = GetCovar(hset,src,tok))==NULL){
            HMError(src,"GetCovar Failed");
            return(NULL);
         }
         if (hset->ckind == LLTC || hset->ckind == NULLC)
            mp->ckind = LLTC;
         else{
            HRError(7032,"GetMixPDF: trying to change global cov type to LLTC");
            return(NULL);
         }
      } else if (tok->sym==XFORM || (tok->sym==MACRO && tok->macroType=='x')){
         if((mp->cov.xform = GetTransform(hset,src,tok))==NULL){
            HMError(src,"GetTransform Failed");
            return(NULL);
         }
         if (hset->ckind == XFORMC || hset->ckind == NULLC)
            mp->ckind = XFORMC;
         else{
            HRError(7032,"GetMixPDF: trying to change global cov type to XFormC");
            return(NULL);
         }
      } else{
         HMError(src,"Variance or Xform expected in GetMixPDF");
         return(NULL);
      }
      if (tok->sym==GCONST) {
         ReadFloat(src,&mp->gConst,1,tok->binForm);
         
         if(GetToken(src,tok)<SUCCESS){
            HMError(src,"GetToken failed");
            return(NULL);
         }
      }
   }
   return mp;  
}

/* CreateCME: create an array of M Continuous MixtureElems */
static MixtureElem *CreateCME(HMMSet *hset, int M)
{
   int m;
   MixtureElem *me,*p;
   
   me = (MixtureElem *)New(hset->hmem,M*sizeof(MixtureElem));
   p = me-1;
   for (m=1;m<=M;m++,me++){
      me->weight = 0.0; me->mpdf = NULL;
   }
   return p;
}

/* CreateTME: create an array of M Tied Mix Weights (ie floats) */
static Vector CreateTME(HMMSet *hset, int M)
{
   int m;
   Vector v;
   
   v = CreateVector(hset->hmem,M);
   for (m=1;m<=M;m++)
      v[m] = 0;
   return v;
}

/* CreateDME: create an array of M Discrete Mix Weights (ie shorts) */
static ShortVec CreateDME(HMMSet *hset, int M)
{
   int m;
   ShortVec v;
   
   v = CreateShortVec(hset->hmem,M);
   for (m=1;m<=M;m++)
      v[m] = 0;
   return v;
}

/* GetMixture: parse src and store a MixtureElem in spdf array */
static ReturnStatus GetMixture(HMMSet *hset,Source *src,Token *tok,int M,MixtureElem *spdf)
{
   float w = 1.0;
   short m = 1;

   if (trace&T_PAR) printf("HModel: GetMixture\n");
   if (tok->sym == MIXTURE) {
      if (!ReadShort(src,&m,1,tok->binForm)){
         HMError(src,"Mixture Index expected");
         return(FAIL);
      }
      if (m<1 || m>M){
         HMError(src,"Mixture index out of range");
         return(FAIL);
      }
      if (!ReadFloat(src,&w,1,tok->binForm)){
         HMError(src,"Mixture Weight expected");
         return(FAIL);
      }
      if(GetToken(src,tok)<SUCCESS){
         HMError(src,"GetToken failed");
         return(FAIL);
      }
   }
   spdf[m].weight = w;
   if((spdf[m].mpdf = GetMixPDF(hset,src,tok))==NULL){
      HMError(src,"Regression Class Number expected");
      return(FAIL);
   }
   return(SUCCESS);
}
   
/* CreateSE: create an array of S StreamElems */
static StreamElem *CreateSE(HMMSet *hset, int S)
{
   int s;
   StreamElem *se,*p;
   
   se = (StreamElem *)New(hset->hmem,S*sizeof(StreamElem));
   p = se-1;
   for (s=1;s<=S;s++,se++){
      se->hook = NULL;
      se->spdf.cpdf = NULL;
   }
   return p;
}

/* EmptyMixPDF: return an empty Diag Covariance MixPDF */
static MixPDF *EmptyMixPDF(HMMSet *hset, int vSize, int s)
{
   int i;
   static Boolean isInitialised = FALSE;
   static MixPDF *t[SMAX];
   static int size[SMAX];

   if (!isInitialised){
      for (i=0; i<SMAX; i++) t[i]=NULL;
      isInitialised = TRUE;
   }
   if (t[s] != NULL) {
      if (size[s] != vSize){
         HRError(7090,"EmptyMixPDF: Size mismatch %d vs %d in EmptyDiagMixPDF",
                 vSize,size[s]);
         return(NULL);
      }
   }
   size[s] = vSize;
   t[s] = (MixPDF *)New(hset->hmem,sizeof(MixPDF));
   t[s]->ckind = DIAGC;
   t[s]->nUse = 0; t[s]->hook = NULL; t[s]->gConst = LZERO;
   t[s]->mean = CreateSVector(hset->hmem,vSize);
   ZeroVector(t[s]->mean);
   t[s]->cov.var = CreateSVector(hset->hmem,vSize);
   for (i=1; i<=vSize; i++) t[s]->cov.var[i] = 1.0;
   return t[s];
}

/* GetStream: parse src and store a StreamElem in pdf array */
static ReturnStatus GetStream(HMMSet *hset, Source *src, Token *tok,
                              StreamElem *pdf, short *nMix)
{
   int m,S,M;
   short s;
   MixtureElem *cpdf;
  
   S=hset->swidth[0];
   if (trace&T_PAR) {
      printf("HModel: GetStream - nMix =");
      for (s=1; s<=S; s++) printf(" %d",nMix[s]);
      printf("\n");
   }
   s = 1;
   if (tok->sym == STREAM) {
      if (!ReadShort(src,&s,1,tok->binForm)){
         HMError(src,"Stream Index expected");
         return(FAIL);
      }
      if (s<1 || s>S){
         HMError(src,"Stream Index out of range");
         return(FAIL);
      }
      if(GetToken(src,tok)<SUCCESS){
         HMError(src,"GetToken failed");
         return(FAIL);
      }
   }
   M = nMix[s];
   pdf[s].nMix = M;
   if (tok->sym == TMIX ) {
      if (hset->hsKind == PLAINHS)
         hset->hsKind = TIEDHS;
      else 
         if (hset->hsKind != TIEDHS){
            HRError(7032,"GetStream: change to TIEDHS from other than PLAINHS");
            return(FAIL);
         }
      pdf[s].spdf.tpdf = CreateTME(hset,M);
      if((GetTiedMixtures(hset,src,tok,M,s,pdf[s].spdf.tpdf))<SUCCESS){
         HMError(src,"GetTiedMixtures failed");
         return(FAIL);
      }
   } 
   else 
      if (tok->sym == DPROB) {
         if (hset->hsKind == PLAINHS)
            hset->hsKind = DISCRETEHS;
         else if (hset->hsKind != DISCRETEHS){
            HRError(7032,"GetStream: change to DISCRETEHS from other than PLAINHS");
            return(FAIL);
         }
         pdf[s].spdf.dpdf = CreateDME(hset,M);
         if((GetDiscreteWeights(src,tok,M,pdf[s].spdf.dpdf))<SUCCESS){
            HMError(src,"GetDiscreteWeights failed");
            return(FAIL);
         }
      } else {  /* PLAIN/SHARED Mixtures */
         cpdf = pdf[s].spdf.cpdf = CreateCME(hset,M);
         if((GetMixture(hset,src,tok,M,cpdf))<SUCCESS){
            HMError(src,"GetMixtures failed");
            return(FAIL);
         }
         while (tok->sym==MIXTURE)
            if((GetMixture(hset,src,tok,M,cpdf))<SUCCESS){
               HMError(src,"GetMixtures failed");
               return(FAIL);
            }
         for (m=1; m<=M; m++)
            if (cpdf[m].mpdf == NULL){
               if((cpdf[m].mpdf = EmptyMixPDF(hset,hset->swidth[s],s))==NULL){
                  HMError(src,"EmptyMixPDF failed");
                  return(FAIL);
               }
               cpdf[m].weight = 0.0;
            }
      }
   return(SUCCESS);
}
  
/* GetStateInfo: parse src and return StateInfo structure  */
static StateInfo *GetStateInfo(HMMSet *hset, Source *src, Token *tok)
{
   StateInfo *si;
   int i,S;
   short nMix[SMAX];

   if (trace&T_PAR) printf("HModel: GetStateInfo\n");
   S = hset->swidth[0];
   if (tok->sym==MACRO && tok->macroType=='s') {
      if((si = (StateInfo *)GetStructure(hset,src,'s'))==NULL){
         HMError(src,"GetStructure failed");
         return(NULL);
      }
      ++si->nUse;
      if(GetToken(src,tok)<SUCCESS){
         HMError(src,"GetToken failed");
         return(NULL);
      }
   } else {
      if (tok->sym == NUMMIXES){
         if (!ReadShort(src,nMix+1,S,tok->binForm)){
            HMError(src,"Num Mix in Each Stream expected");
            return(NULL);
         }
         if(GetToken(src,tok)<SUCCESS){
            HMError(src,"GetToken failed");
            return(NULL);
         }
      } else {
         for (i=1;i<=S;i++)
            nMix[i] = 1;
      }
      si = (StateInfo *)New(hset->hmem,sizeof(StateInfo));
      si->nUse = 0; si->hook = NULL; si->weights = NULL;
      si->pdf = CreateSE(hset,S);
      if (tok->sym==SWEIGHTS || (tok->sym==MACRO && tok->macroType=='w')){
         if((si->weights = GetSWeights(hset,src,tok))==NULL){
            HMError(src,"GetSWeights failed");
            return(NULL);
         }
         if (VectorSize(si->weights) != S){
            HMError(src,"Incorrect number of stream weights");
            return(NULL);
         }
      }
      if((GetStream(hset,src,tok,si->pdf,nMix))<SUCCESS){
         HMError(src,"GetStream failed");
         return(NULL);
      }
      while(tok->sym==STREAM)
         if((GetStream(hset,src,tok,si->pdf,nMix))<SUCCESS){
            HMError(src,"GetStream failed");
            return(NULL);
         }
      if (tok->sym==DURATION || (tok->sym==MACRO && tok->macroType=='d')){
         if((si->dur = GetDuration(hset,src,tok))==NULL){
            HMError(src,"GetDuration failed");
            return(NULL);
         }
      }
      else
         si->dur = NULL;
   }
   if (S>1 && si->weights == NULL) {
      si->weights = CreateSVector(hset->hmem,S);
      for (i=1;i<=S;i++)
         si->weights[i] = 1.0;
   }
   return si;  
}

/* GetTransMat: parse src and return Transition Matrix structure */
static SMatrix GetTransMat(HMMSet *hset, Source *src, Token *tok)
{
   SMatrix m;
   int i,j;
   short size;
   Vector v;
   float rSum;
   
   if (trace&T_PAR) printf("HModel: GetTransMat\n");
   if (tok->sym == TRANSP) {
      if (!ReadShort(src,&size,1,tok->binForm)){
         HMError(src,"Size of Transition matrix expected");
         return(NULL);
      }
      if (size < 1){
         HRError(7031,"GetTransMat: Bad size of transition matrix: %d\n", size);
         return(NULL);
      }
      m = CreateSMatrix(hset->hmem,size,size);
      if (!ReadMatrix(src,m,tok->binForm)){
         HMError(src,"Transition Matrix expected");
         return(NULL);
      }
      if (!hset->allowTMods && (m[1][size] > 0.0)) { /* kill teeModel */
         rSum = 0.0; v = m[1]; v[size] = 0.0;
         for (j=1;j<size;j++) rSum += v[j];
         for (j=1;j<size;j++) v[j] /= rSum;
      }     
      for (i=1;i<size;i++){                       /* convert to logs */
         v=m[i]; rSum = 0.0;
         for (j=1;j<=size;j++){
            rSum += v[j];
            v[j] = (v[j]<=MINLARG) ? LZERO : log(v[j]);
         }
         if (rSum<0.99 || rSum>1.01){
            HRError(7031,"GetTransMat: Bad Trans Mat Sum in Row %d\n",i);
            return(NULL);
         }
      }
      v = m[size];
      for (j=1;j<=size;j++)
         v[j] =  LZERO;
   } else {
      if((m = (SMatrix)GetStructure(hset,src,'t'))==NULL){
         HMError(src,"GetStructure failed");
         return(NULL);
      }
      IncUse(m);
   }
   if(GetToken(src,tok)<SUCCESS){
      HMError(src,"GetToken failed");
      return(NULL);
   }
   return m;
}

/* GetHMMDef: get a hmm def from given source */
static ReturnStatus GetHMMDef(HMMSet *hset, Source *src, Token *tok,
                              HLink hmm, int nState)
{
   short state;
   int N=0;
   StateElem *se;
   char buf[MAXSTRLEN];
   
   if (trace&T_PAR) printf("HModel: GetHMMDef\n");
   if (tok->sym != BEGINHMM){
      HMError(src,"<BeginHMM> symbol expected");
      return(FAIL);
   }
   if(GetToken(src,tok)<SUCCESS){
      HMError(src,"GetHMMDef: GetToken failed");
      return(FAIL);
   }
   if (tok->sym == USEMAC){      /* V1 style USE clause */
      if (!ReadString(src,buf)){
         HRError(7013,"GetHMMDef: cannot read USE macro name");
         return(FAIL);
      }
      if(GetToken(src,tok)<SUCCESS){
         HMError(src,"GetHMMDEf: GetToken failed");
         return(FAIL);
      }
   }
   while (tok->sym != STATE) {
      if(GetOption(hset,src,tok, &N)<SUCCESS){
         HMError(src,"GetHMMDef: GetOption failed");
         return(FAIL);
      }
      if (N>nState) nState = N;
   }
   FreezeOptions(hset);
   if (nState==0){
      HMError(src,"NumStates not set");
      return(FAIL);
   }
   hmm->numStates = N = nState;
   se = (StateElem *)New(hset->hmem,(N-2)*sizeof(StateElem));
   hmm->svec = se - 2;
   while (tok->sym == STATE) {
      if (!ReadShort(src,&state,1,tok->binForm)){
         HMError(src,"States index expected");
         return(FAIL);
      }
      if (state<2 || state >= N){
         HMError(src,"State index out of range");
         return(FAIL);
      }
      if(GetToken(src,tok)<SUCCESS){
         HMError(src,"GetHMMDef: GetToken failed");
         return(FAIL);
      }
      se = hmm->svec+state;
      if((se->info = GetStateInfo(hset,src,tok))==NULL){
         HMError(src,"GetStateInfo failed");
         return(FAIL);
      }     
   }
   if (tok->sym==TRANSP || (tok->sym==MACRO && tok->macroType=='t')){
      if((hmm->transP = GetTransMat(hset,src,tok))==NULL){
         HMError(src,"GetTransMat failed");
         return(FAIL);
      }     
      if (NumRows(hmm->transP) != N ||
          NumCols(hmm->transP) != N){
         HRError(7030,"GetHMMDef: Trans Mat Dimensions not %d x %d",N,N);
         return(FAIL);
      }
   }
   else{
      HMError(src,"Transition Matrix Missing");
      return(FAIL);
   }
   if (tok->sym==DURATION || (tok->sym==MACRO && tok->macroType=='d')){
      if((hmm->dur = GetDuration(hset,src,tok))==NULL){
         HMError(src,"GetDuration Failed");
         return(FAIL);
      }
   }
   else
      hmm->dur = NULL;
   if (tok->sym!=ENDHMM){
      HMError(src,"<EndHMM> symbol expected");
      return(FAIL);
   }
   if(GetToken(src,tok)<SUCCESS){
      HMError(src,"GetHMMDef: GetToken failed");
      return(FAIL);
   }
   return(SUCCESS);
}

/* GetBias: parse src and return Bias structure */
static SVector GetBias(HMMSet *hset, Source *src, Token *tok)
{
   SVector m = NULL;
   MemHeap *hmem;
   short size;
   
   if (trace&T_PAR) printf("HModel: GetBias\n");
   if (tok->sym==BIAS) {      
      if (hset==NULL) hmem = &xformStack;
      else hmem = hset->hmem;
      if (!ReadShort(src,&size,1,tok->binForm)){
         HMError(src,"Size of Bias Vector expected");
         return(NULL);
      }
      m = CreateSVector(hmem,size);
      if (!ReadVector(src,m,tok->binForm)){
         HMError(src,"Bias Vector expected");
         return(NULL);
      }
   }
   
   else if ((tok->sym==MACRO && tok->macroType=='y') && (hset != NULL)) {
      if((m=(SVector)GetStructure(hset,src,'y'))==NULL){
         HMError(src,"GetStructure Failed");
         return(NULL);
      }
      IncUse(m);
   } else{
      HMError(src,"<BIAS> symbol expected in GetBias");
      return(NULL);
   }
   if(GetToken(src,tok)<SUCCESS){
      HMError(src,"GetToken failed");
      return(NULL);
   }
   return m;
}

/* GetLinXForm: get a linear transformations */
static LinXForm* GetLinXForm(HMMSet *hset, Source *src, Token *tok)
{
   LinXForm *xf;
   MemHeap *hmem;
   int i,b;
   int numBlocks;

   if (trace&T_PAR) printf("HModel: GetLinXForm\n");
   if (tok->sym == VECSIZE) {
      if (hset==NULL) hmem = &xformStack;
      else hmem = hset->hmem;
      xf = (LinXForm *)New(hmem,sizeof(LinXForm));
      if (!ReadInt(src,&(xf->vecSize),1,tok->binForm)){
         HRError(7013,"GetLinXForm: cannot read vector size");
         return(NULL);
      }
      if(GetToken(src,tok)<SUCCESS){
         HMError(src,"GetToken failed");
         return(NULL);
      }
      if (tok->sym == OFFSET){
         if(GetToken(src,tok)<SUCCESS){
            HMError(src,"GetToken failed");
            return(NULL);
         }
         xf->bias = GetBias(hset,src,tok);
      } else {
         xf->bias = NULL;
      }
      if (tok->sym!=BLOCKINFO){
         HMError(src,"<BLOCKINFO> symbol expected");
         return(NULL);
      }
      if (!ReadInt(src,&numBlocks,1,tok->binForm)){
         HMError(src,"Number of transform blocks expected");
         return(NULL);
      }
      xf->blockSize = CreateIntVec(hmem,numBlocks);
      if (!ReadInt(src,xf->blockSize+1,numBlocks,tok->binForm)){
         HMError(src,"Size of blocks expected");
         return(NULL);
      }    
      xf->xform = (Matrix *)New(hmem,(numBlocks+1)*sizeof(Matrix));
      if(GetToken(src,tok)<SUCCESS){
         HMError(src,"GetToken failed");
         return(NULL);
      }
      for (i=1;i<=numBlocks;i++) {
         if (tok->sym != BLOCK){
            HMError(src,"<NUMXFORMS> symbol expected");
            return(NULL);
         }
         if (!ReadInt(src,&b,1,tok->binForm)){
            HMError(src,"Weight class expected");
            return(NULL);
         }      
         if (b != i) {
            HMError(src,"Inconsistency in transform definition");
            return(NULL);
         }
         if(GetToken(src,tok)<SUCCESS){
            HMError(src,"GetToken failed");
            return(NULL);
         }
         xf->xform[i] = GetTransform(hset,src,tok);
      }
      xf->nUse = 0;
   }
   else if ((tok->sym==MACRO && tok->macroType=='f') && (hset != NULL)){
      if((xf=(LinXForm *)GetStructure(hset,src,'f'))==NULL){
         HMError(src,"GetStructure Failed");
         return(NULL);
      }
      xf->nUse++;
      if(GetToken(src,tok)<SUCCESS){
         HMError(src,"GetToken failed");
         return(NULL);
      }
   } else{
      HMError(src,"<VECSIZE> symbol expected in GetXFormSet");
      return(NULL);
   }
   return xf;
}

/* GetInputXForm: get the input linear transformations */
static InputXForm* GetInputXForm(HMMSet *hset, Source *src, Token *tok)
{
   char buf[MAXSTRLEN];
   MemHeap *hmem;
   InputXForm *xf;

   if (trace&T_PAR) printf("HModel: GetXForm\n");
   if (tok->sym == MMFIDMASK) { 
      if (hset==NULL) hmem = &xformStack;
      else hmem = hset->hmem;
      xf = (InputXForm *)New(hmem,sizeof(InputXForm));
      xf->fname = CopyString(hmem,src->name);
      xf->xformName = NULL;
      if (!ReadString(src,buf)){
         HRError(7013,"GetInputXForm: cannot read MMFIDMASK");
         return(NULL);
      }
      xf->mmfIdMask = CopyString(hmem,buf);
      if(GetToken(src,tok)<SUCCESS){
         HMError(src,"GetToken failed");
         return(NULL);
      }
      if (tok->sym != PARMKIND) {
         HMError(src,"parameter kind symbol expected");
         return(NULL);
      }
      xf->pkind = tok->pkind;
      if(GetToken(src,tok)<SUCCESS){
         HMError(src,"GetToken failed");
         return(NULL);
      }
      if (tok->sym == PREQUAL) {
         xf->preQual = TRUE;
         if(GetToken(src,tok)<SUCCESS){
            HMError(src,"GetToken failed");
            return(NULL);
         }
      } else 
         xf->preQual = FALSE;
      if (tok->sym != LINXFORM){
         HMError(src,"<LINXFORM> symbol expected");
         return(NULL);
      }
      if(GetToken(src,tok)<SUCCESS){
         HMError(src,"GetToken failed");
         return(NULL);
      }
      xf->xform = GetLinXForm(hset,src,tok);
      xf->nUse=0;
   }
   else if ((tok->sym==MACRO && tok->macroType=='j') && (hset != NULL)){
      if((xf=(InputXForm *)GetStructure(hset,src,'j'))==NULL){
         HMError(src,"GetStructure Failed");
         return(NULL);
      }
      xf->nUse++;
      if(GetToken(src,tok)<SUCCESS){
         HMError(src,"GetToken failed");
         return(NULL);
      }
   } else{
      HMError(src,"<MMFIDMASK> symbol expected in GetInputXForm");
      return(NULL);
   }
   return xf;
}

/* ---------------------- Symbol Output Routine ------------------------ */

/* PutSymbol: output symbol to f in ascii or binary form */
static void PutSymbol(FILE *f, Symbol sym, Boolean binary)
{
   if (binary){
      fputc(':',f);fputc(sym,f);
   } else
      fprintf(f,"<%s>",symNames[sym]);
}

/* ------------- Special Tied Mixture Output Routines ------------------ */

static Boolean putWtActive;   /* false until a weight is output */

/* PutTiedWeight: output mix weight via 1-item buffer */
void PutTiedWeight(FILE *f, short repeatLast, float w, Boolean binary)
{
   static float putWtCurrent;  /* the current output weight */
   
   if (putWtActive)  {        /* not 1st time so output last w */
      if (binary){
         if (repeatLast>0) {
            putWtCurrent = putWtCurrent-2.0;
            WriteFloat(f,&putWtCurrent,1,binary);
            fputc(repeatLast,f);
         } else
            WriteFloat(f,&putWtCurrent,1,binary);
      } else {
         fprintf(f," %e",putWtCurrent);
         if (repeatLast>0) 
            fprintf(f,"*%d",repeatLast);
      }
   } else {
      putWtActive = TRUE;
   }
   putWtCurrent = w;
   if (!binary && w<0) /* all done */
      fprintf(f,"\n");
}

/* PutTiedWeights: output mixture weights in compact format */
void PutTiedWeights(FILE *f, StreamElem *se, Boolean binary)
{
   int repCount=0;
   float weight= -1;
   int m,M;
   Vector v;
   
   putWtActive = FALSE;
   M = se->nMix; v = se->spdf.tpdf;
   for (m=1; m<=M; m++){
      if (v[m] == weight && repCount < 255)
         ++repCount;
      else{
         if (repCount>0){
            PutTiedWeight(f,repCount+1,v[m],binary);
            repCount=0;
         } else
            PutTiedWeight(f,0,v[m],binary);
         weight = v[m];
      }
   }
   if (repCount>0){
      PutTiedWeight(f,repCount+1,-1,binary);
   }else
      PutTiedWeight(f,0,-1,binary);
}

/* PutMixWeight: output mix weight via 1-item buffer */
void PutMixWeight(FILE *f, short repeatLast, short w, Boolean binary)
{
   static short   putWtCurrent;  /* the current output weight */
   static short   putWtLCount;   /* num wts output on this line */
   
   if (putWtActive)  {        /* not 1st time so output last w */
      if (binary){
         if (repeatLast>0) {
            putWtCurrent |= 0100000;
            WriteShort(f,&putWtCurrent,1,binary);
            fputc(repeatLast,f);
         } else
            WriteShort(f,&putWtCurrent,1,binary);
      } else {
         fprintf(f," %d",putWtCurrent);
         if (repeatLast>0) 
            fprintf(f,"*%d",repeatLast);
         if (++putWtLCount > 8) {
            fprintf(f,"\n"); putWtLCount = 0;
         }
      }
   } else {
      putWtActive = TRUE; putWtLCount = 0;
   }
   putWtCurrent = w;
   if (!binary && putWtLCount > 0 && w<0) /* all done */
      fprintf(f,"\n");
}

/* PutDiscreteWeights: output mixture weights in compact format */
void PutDiscreteWeights(FILE *f, StreamElem *se, Boolean binary)
{
   int repCount=0;
   short weight= -1;
   int m,M;
   ShortVec v;
   
   putWtActive = FALSE;
   M = se->nMix; v = se->spdf.dpdf;
   for (m=1; m<=M; m++){
      if (v[m] == weight && repCount < 255)
         ++repCount;
      else{
         if (repCount>0){
            PutMixWeight(f,repCount+1,v[m],binary);
            repCount=0;
         } else
            PutMixWeight(f,0,v[m],binary);
         weight = v[m];
      }
   }
   if (repCount>0){
      PutMixWeight(f,repCount+1,-1,binary);
   }else
      PutMixWeight(f,0,-1,binary);
}

/* PutTiedMixtures: output mixture weights in compact tied mix format */
void PutTiedMixtures(HMMSet *hset,FILE *f,int s,StreamElem *se,Boolean binary)
{
   PutSymbol(f,TMIX,binary);
   fprintf(f," %s",ReWriteString(hset->tmRecs[s].mixId->name,NULL,DBL_QUOTE));
   if (!binary) fprintf(f,"\n");
   PutTiedWeights(f,se,binary);
}

/* PutDiscrete: output discrete weights in compact format */
void PutDiscrete(FILE *f, StreamElem *se, Boolean binary)
{
   PutSymbol(f,DPROB,binary);
   if (!binary) fprintf(f,"\n");
   PutDiscreteWeights(f,se,binary);
}

/* ---------- Standard Macro Definition Output Routines ------------------ */

/* PutMacroHdr: write macro name to f.  If m is null find structure
   corresponding to ptr */
static void PutMacroHdr(HMMSet *hset, FILE *f, MLink m, char mType, 
                        Ptr ptr, Boolean binary)
{
   if (m==NULL)   
      m = FindMacroStruct(hset,mType,ptr);
   if (m != NULL) {
      if (m->type != mType)
         HError(7091,"PutMacroHdr: Macro type error ~%c vs ~%c",m->type, mType);
      fprintf(f,"~%c %s",mType,ReWriteString(m->id->name,NULL,DBL_QUOTE));
      if (!binary) fprintf(f,"\n");
      return;
   }
   HError(7035,"PutMacroHdr: Unable to find ~%c Macro to write",mType);
}


/* PutMean: output mean vector to stream f */
static void PutMean(HMMSet *hset, FILE *f, MLink q, SVector m, 
                    Boolean inMacro, Boolean binary)
{
   int nUse;
   short size;
   
   nUse = GetUse(m);
   if (nUse > 0 || inMacro) 
      PutMacroHdr(hset,f,q,'u',m,binary);
   if (nUse == 0 || inMacro) {
      PutSymbol(f,MEAN,binary);
      size = VectorSize(m);
      WriteShort(f,&size,1,binary);
      if (!binary) fprintf(f,"\n");
      WriteVector(f,m,binary);
   }
}

/* PutVariance: output variance vector to stream f */
static void PutVariance(HMMSet *hset, FILE *f, MLink q, SVector v,
                        Boolean inMacro, Boolean binary)
{
   int nUse;
   short size;
   
   nUse = GetUse(v);
   if (nUse > 0 || inMacro) 
      PutMacroHdr(hset,f,q,'v',v,binary);
   if (nUse == 0 || inMacro){
      PutSymbol(f,VARIANCE,binary);
      size = VectorSize(v);
      WriteShort(f,&size,1,binary);
      if (!binary) fprintf(f,"\n");
      WriteVector(f,v,binary);
   }
}

/* PutCovar: output inverse/choleski cov matrix to stream f */
static void PutCovar(HMMSet *hset, FILE *f, MLink q, STriMat m,
                     Symbol sym, Boolean inMacro, Boolean binary)
{
   int nUse;
   short size;
   
   nUse = GetUse(m);
   if (nUse > 0 || inMacro) 
      PutMacroHdr(hset,f,q,sym==INVCOVAR?'i':'c',m,binary);
   if (nUse == 0 || inMacro){
      PutSymbol(f,sym,binary);
      size = NumRows(m);
      WriteShort(f,&size,1,binary);
      if (!binary) fprintf(f,"\n");
      WriteTriMat(f,m,binary);
   }
}

/* PutTransform: output xform matrix to stream f */
static void PutTransform(HMMSet *hset, FILE *f, MLink q, SMatrix m,
                         Boolean inMacro, Boolean binary)
{
   int nUse;
   short nr,nc;
   
   nUse = GetUse(m);
   if (nUse > 0 || inMacro) 
      PutMacroHdr(hset,f,q,'x',m,binary);
   if (nUse == 0 || inMacro){
      PutSymbol(f,XFORM,binary);
      nr = NumRows(m); nc = NumCols(m);
      WriteShort(f,&nr,1,binary);
      WriteShort(f,&nc,1,binary);
      if (!binary) fprintf(f,"\n");
      WriteMatrix(f,m,binary);
   }
}

/* PutDuration: output duration vector to stream f */
static void PutDuration(HMMSet *hset, FILE *f, MLink q, SVector v,
                        Boolean inMacro, Boolean binary)
{
   int nUse;
   short size;
   
   nUse = GetUse(v);
   if (nUse > 0 || inMacro) 
      PutMacroHdr(hset,f,q,'d',v,binary);
   if (nUse == 0 || inMacro){
      PutSymbol(f,DURATION,binary);
      size = VectorSize(v);
      WriteShort(f,&size,1,binary);
      if (!binary) fprintf(f,"\n");
      WriteVector(f,v,binary);
   }
}

/* PutSWeights: output stream weight vector to stream f */
static void PutSWeights(HMMSet *hset, FILE *f, MLink q, SVector v,
                        Boolean inMacro, Boolean binary)
{
   int nUse;
   short size;

   nUse = GetUse(v);
   if (nUse > 0 || inMacro) 
      PutMacroHdr(hset,f,q,'w',v,binary);
   if (nUse == 0 || inMacro){
      PutSymbol(f,SWEIGHTS,binary);
      size = VectorSize(v);
      WriteShort(f,&size,1,binary);
      if (!binary) fprintf(f,"\n");
      WriteVector(f,v,binary);
   }
}

/* PutTransMat: output transition matrix to stream f */
static void PutTransMat(HMMSet *hset, FILE *f, MLink q, SMatrix m,
                        Boolean inMacro, Boolean binary)
{
   Matrix mm;
   Vector v;
   int i,j,nUse;
   short nstates;
   float rSum;
   
   nUse = GetUse(m);
   if (nUse>0 || inMacro) 
      PutMacroHdr(hset,f,q,'t',m,binary);
   if (nUse == 0 || inMacro){
      nstates = NumRows(m);
      PutSymbol(f,TRANSP,binary);
      WriteShort(f,&nstates,1,binary);
      if (!binary) fprintf(f,"\n");
      mm=CreateMatrix(&gstack,nstates,nstates);
      CopyMatrix(m,mm);
      for (i=1;i<nstates;i++){    /* convert to real */
         v = mm[i]; rSum = 0.0;     /* then normalise */
         for (j=1; j<=nstates;j++) {
            v[j] = L2F(v[j]);
            rSum += v[j];
         }
         if (rSum==0.0)
            HError(7031,"PutTransMat: Row %d of transition mat is all zero",i);
         if (rSum < 0.99 || rSum > 1.01)
            HError(7031,"PutTransMat: Row %d of transition mat sum = %f\n",i,rSum);
         for (j=1; j<=nstates;j++)
            v[j] /= rSum;
      }
      v = mm[nstates];  /* last row is always 0.0 */
      for (j=1; j<=nstates;j++) v[j] = 0.0;
      WriteMatrix(f,mm,binary);
      FreeMatrix(&gstack,mm);
   }
}

/* ----------------- Regression Tree Output Functions -------------------- */
/* linked list of components that are grouped into a cluster */
typedef struct _CoList {
  
   struct _CoList *next;     /* next component in the linked list */
   char   *hmmName;          /* Physical hmm name for this component */
   int    state;             /* state containing this component */
   int    stream;            /* stream for this component */
   int    mix;               /* mixture for this component */
   MixPDF *mp;               /* actual component */
  
} CoList;

/* node information in the regression class tree */
typedef struct {
  
   Vector aveMean;           /* node cluster mean */
   Vector aveCovar;          /* node cluster variance */
   float  clusterScore;      /* node cluster score */
   float  clustAcc;          /* accumulates in this cluster */
   int    nComponents;       /* number of components in this cluster */
   short  nodeIndex;         /* node index number */
   CoList *list;             /* linked list of the mixture components */
  
} RNode;

/* return the node of a tree */
static RNode *GetRNode(RegTree *t) {

   return ((RNode *) t->node);

}

static void GetTreeVector(RegTree **tVec, RegTree *t) {
  
   RNode *n;
   short index;

   if (t != NULL) {
      GetTreeVector(tVec, t->left);
      if (t->nodeInfo == NULL) {
         n = GetRNode(t);
         index = n->nodeIndex;
      }
      else
         index = t->nodeInfo->nodeIndex;
      tVec[index] = t;
      GetTreeVector(tVec, t->right);
   } 
  
}

/* find the number of terminals in the tree (returned in nTerms) */
static void FindNumTerminals(RegTree *t, short *nTerms) 
{
   if (t != NULL) {
      FindNumTerminals(t->left, nTerms);
      if (t->left == NULL)
         *nTerms += 1;
      FindNumTerminals(t->right, nTerms);
   }
}


/* PutRegTree: output the regression class tree */
static void PutRegTree(HMMSet *hset, FILE *f, MLink q, RegTree *t, 
                       Boolean inMacro, Boolean binary) {
  

   RegTree **tVec=NULL;
   RNode *n=NULL;
   int i, nNodes=0;
   short nTerms=0;
   short index;
   int nComps;

   PutMacroHdr(hset,f,q,'r',t,binary);
   FindNumTerminals(t, &nTerms);
   nNodes = nTerms*2 - 1;
   tVec = (RegTree **) New(&gstack, nNodes*sizeof(RegTree *));
   --tVec;
   GetTreeVector(tVec, t);
   /* Write the regression tree header */
   PutSymbol(f, REGTREE, binary);
   WriteShort(f,&nTerms,1,binary);
   if (!binary) fprintf(f,"\n");
   /* Write the rest of the tree */
   for (i = 1; i <= nNodes; i++) {
      if (tVec[i]->nodeInfo == NULL) {
         n = GetRNode(tVec[i]);
         index = n->nodeIndex;
         nComps = n->nComponents;
      }
      else {
         index = tVec[i]->nodeInfo->nodeIndex;
         nComps = tVec[i]->nodeInfo->nodeComps;
      }
      if (tVec[i]->left == NULL) {
         PutSymbol(f, TNODE, binary);
         WriteShort(f,&index,1,binary);
         WriteInt(f,&nComps,1,binary);
         /* WriteFloat(f,&(n->clustAcc),1,binary);
            WriteFloat(f,&(n->clusterScore),1,binary); */
         if (!binary) fprintf(f,"\n");
      }
      else {
         PutSymbol(f, NODE, binary);
         WriteShort(f,&index,1,binary);
         if (tVec[i]->left->nodeInfo == NULL) {
            n = GetRNode(tVec[i]->left);
            index = n->nodeIndex;
         }
         else
            index = tVec[i]->left->nodeInfo->nodeIndex;
         WriteShort(f,&index,1,binary);      
         if (tVec[i]->right->nodeInfo == NULL) {
            n = GetRNode(tVec[i]->right);
            index = n->nodeIndex;
         }
         else
            index = tVec[i]->right->nodeInfo->nodeIndex;
         WriteShort(f,&index,1,binary);
         if (!binary) fprintf(f,"\n");
      } 
   }
   /* Dispose(&gstack, tVec); */
}


/* PutMixPDF: output mixture pdf to stream f */
static void PutMixPDF(HMMSet *hset, FILE *f, MLink q, MixPDF *mp, 
                      Boolean inMacro, Boolean binary)
{
   if (mp->nUse >0 || inMacro)  
      PutMacroHdr(hset,f,q,'m',mp,binary);
   if (mp->nUse == 0 || inMacro){
      if (saveRegClass && mp->rClass!=0) {
         PutSymbol(f,RCLASS,binary);
         WriteShort(f,&mp->rClass,1,binary);
         if (!binary) fprintf(f,"\n");
      }
      PutMean(hset,f,NULL,mp->mean,FALSE,binary);
      switch(mp->ckind){
      case DIAGC:  PutVariance(hset,f,NULL,mp->cov.var,FALSE,binary); break;
      case INVDIAGC: HError(7021,"PutMixPDF: Cannot output INVDIAGC covariance"); break;
      case FULLC:  PutCovar(hset,f,NULL,mp->cov.inv,INVCOVAR,FALSE,binary); break;
      case LLTC:   PutCovar(hset,f,NULL,mp->cov.inv,LLTCOVAR,FALSE,binary); break;
      case XFORMC: PutTransform(hset,f,NULL,mp->cov.xform,FALSE,binary); break;
      default:     HError(7033,"PutMixPDF: bad ckind %d",mp->ckind); break;
      }
      if (mp->gConst != LZERO) {
         PutSymbol(f,GCONST,binary);
         WriteFloat(f,&mp->gConst,1,binary);
         if (!binary) fprintf(f,"\n");
      }
   }  
}

/* PutStateInfo: output state info to stream f */
static void PutStateInfo(HMMSet *hset, FILE *f, MLink q, StateInfo *si, 
                         Boolean inMacro, Boolean binary)
{
   int S;
   float wt;
   StreamElem *se;
   MixtureElem *me;
   short m,s,nMix[SMAX];
   Boolean needNM = FALSE;
   
   S = hset->swidth[0];
   if (si->nUse > 0 || inMacro) 
      PutMacroHdr(hset,f,q,'s',si,binary);
   if (si->nUse == 0 || inMacro){
      se = si->pdf+1;
      for (s=1; s<=S; s++,se++) {
         if (se->nMix>1) needNM = TRUE;
         nMix[s] = se->nMix;
      }
      if (needNM){
         PutSymbol(f,NUMMIXES,binary);
         WriteShort(f,nMix+1,S,binary);
         if (!binary) fprintf(f,"\n");
      }
      if (si->weights != NULL)
         PutSWeights(hset,f,NULL,si->weights,FALSE,binary);
      se = si->pdf+1;
      for (s=1; s<=S; s++,se++){
         if (S>1){
            PutSymbol(f,STREAM,binary);
            WriteShort(f,&s,1,binary);
            if (!binary) fprintf(f,"\n");
         }
         switch (hset->hsKind) {
         case TIEDHS:
            PutTiedMixtures(hset,f,s,se,binary);
            break;
         case DISCRETEHS:
            PutDiscrete(f,se,binary);
            break;
         case PLAINHS:
         case SHAREDHS:
            me = se->spdf.cpdf+1;
            for (m=1; m<=se->nMix; m++,me++){
               wt=me->weight; wt = MixWeight(hset,wt);
               if (wt > MINMIX) {
                  if (se->nMix > 1){
                     PutSymbol(f,MIXTURE,binary);
                     WriteShort(f,&m,1,binary);
                     WriteFloat(f,&wt,1,binary);
                     if (!binary) fprintf(f,"\n");
                  }
                  PutMixPDF(hset,f,NULL,me->mpdf,FALSE,binary);
               }
            }
            break;
         }
      }
      if (hset->dkind!=NULLD && si->dur != NULL)
         PutDuration(hset,f,NULL,si->dur,FALSE,binary);
   }
}

/* PutBias: output bias vector to stream f */
static void PutBias(HMMSet *hset, FILE *f, MLink q, SVector m, 
                    Boolean inMacro, Boolean binary)
{
   int nUse;
   short size;
   
   nUse = GetUse(m);
   if (nUse > 0 || inMacro) 
      PutMacroHdr(hset,f,q,'y',m,binary);
   if (nUse == 0 || inMacro) {
      PutSymbol(f,BIAS,binary);
      size = VectorSize(m);
      WriteShort(f,&size,1,binary);
      if (!binary) fprintf(f,"\n");
      WriteVector(f,m,binary);
   }
}

static void PutLinXForm(HMMSet *hset, FILE *f, MLink q, LinXForm *xf, 
                        Boolean inMacro, Boolean binary)
{
   int i;

   if (xf->nUse > 0 || inMacro) 
      PutMacroHdr(hset,f,q,'f',xf,binary);
   if (xf->nUse == 0 || inMacro){
      PutSymbol(f,VECSIZE,binary);
      WriteInt(f,&(xf->vecSize),1,binary);
      if (!binary) fprintf(f,"\n");
      if (xf->bias != NULL) {
         PutSymbol(f,OFFSET,binary);
         if (!binary) fprintf(f,"\n");
         PutBias(hset,f,NULL,xf->bias,inMacro,binary);      
      }
      PutSymbol(f,BLOCKINFO,binary);
      WriteInt(f,xf->blockSize,1,binary);
      for (i=1;i<=IntVecSize(xf->blockSize);i++) {
         WriteInt(f,xf->blockSize+i,1,binary);
      }
      if (!binary) fprintf(f,"\n");
      for (i=1;i<=IntVecSize(xf->blockSize);i++) {
         PutSymbol(f,BLOCK,binary);
         WriteInt(f,&i,1,binary);
         if (!binary) fprintf(f,"\n");
         PutTransform(hset,f,NULL,xf->xform[i],inMacro,binary);
      }
   }
}

static void PutInputXForm(HMMSet *hset, FILE *f, MLink q, InputXForm *xf, 
			  Boolean inMacro, Boolean binary)
{
   char buf[MAXSTRLEN];
  
   if (xf->nUse > 0 || inMacro) 
      PutMacroHdr(hset,f,q,'j',xf,binary);
   if (xf->nUse == 0 || inMacro){
      PutSymbol(f,MMFIDMASK,binary);
      fprintf(f," %s ",xf->mmfIdMask);
      fprintf(f,"<%s>", ParmKind2Str(xf->pkind,buf));
      if (xf->preQual)
         PutSymbol(f,PREQUAL,binary);
      if (!binary) fprintf(f,"\n");
      PutSymbol(f,LINXFORM,binary);
      PutLinXForm(hset,f,NULL,xf->xform,FALSE,binary);
   }
}

/* PutOptions: write the current global options to f */
static void PutOptions(HMMSet *hset, FILE *f, Boolean binary)
{
   short i,S;
   char buf[64];

   if (hset->hmmSetId!=NULL) {
      PutSymbol(f,HMMSETID,binary);
      fprintf(f," %s\n", hset->hmmSetId);
   }
   S = hset->swidth[0];
   PutSymbol(f,STREAMINFO,binary);
   WriteShort(f,hset->swidth,1,binary);
   for (i=1;i<=S;i++)
      WriteShort(f,hset->swidth+i,1,binary);
   if (!binary) fprintf(f,"\n");
   PutSymbol(f,VECSIZE,binary);
   WriteShort(f,&hset->vecSize,1,binary);
   PutSymbol(f, (Symbol) (NDUR+hset->dkind), binary);
   fprintf(f,"<%s>", ParmKind2Str(hset->pkind,buf));
   if (hset->ckind != NULLC)
      fprintf(f,"<%s>", CovKind2Str(hset->ckind,buf));
   if (!binary) fprintf(f,"\n");
   if ((saveInputXForm) && (hset->xf != NULL)) {
      PutSymbol(f,INPUTXFORM, binary);
      PutInputXForm(hset,f,NULL,hset->xf,FALSE,binary);          
   }
}

/* PutHMMDef: Save the model hmm to given stream in either text or binary */
static void PutHMMDef(HMMSet *hset, FILE *f, MLink m, Boolean withHdr,
                      Boolean binary)
{
   short i;
   StateElem *se;
   HLink hmm;

   if (withHdr)
      PutMacroHdr(hset,f,m,'h',NULL,binary);
   hmm = (HLink)m->structure;
   PutSymbol(f,BEGINHMM,binary); 
   if (!binary) fprintf(f,"\n");
   PutSymbol(f,NUMSTATES,binary);
   WriteShort(f,&hmm->numStates,1,binary);
   if (!binary) fprintf(f,"\n");
   se = hmm->svec+2;
   for (i=2; i<hmm->numStates; i++,se++){
      PutSymbol(f,STATE,binary);
      WriteShort(f,&i,1,binary);
      if (!binary) fprintf(f,"\n");
      PutStateInfo(hset,f,NULL,se->info,FALSE,binary);
   }
   PutTransMat(hset,f,NULL,hmm->transP,FALSE,binary);
   if (hset->dkind != NULLD && hmm->dur != NULL)
      PutDuration(hset,f,NULL,hmm->dur,FALSE,binary);
   PutSymbol(f,ENDHMM,binary);
   if (!binary) fprintf(f,"\n");
}

/* ---------------- Macro Related Manipulations -------------------- */

/* MakeHashTab: make a macro hash table and initialise it */
void ** MakeHashTab(HMMSet *hset, int size)
{
   void **p;
   int i;
   
   p = (void **) New(hset->hmem,sizeof(void *)*size);
   for (i=0; i<size; i++)
      p[i] = NULL;
   return p;
}
   
/* Hash: return a hash value for given name */
static unsigned Hash(char *name)
{
   unsigned hashval;

   for (hashval=0; *name != '\0'; name++)
      hashval = *name + 31*hashval;
   return hashval%MACHASHSIZE;
}

/* EXPORT-> NewMacro: append a macro with given values to list */
MLink NewMacro(HMMSet *hset, short fidx, char type, LabId id, Ptr structure)
{
   unsigned int hashval;
   MLink m;
   PtrMap *p;

   if (FindMacroName(hset,type,id) != NULL){
      HError(7036,"NewMacro: macro or model name %s already exists",id->name);
   }

   m = (MLink)New(hset->hmem,sizeof(MacroDef));
   if (type == 'h')
      ++hset->numPhyHMM;
   else if (type == 'l'){
      ++hset->numLogHMM;
   }
   else
      ++hset->numMacros;
   hashval = Hash(id->name);
   m->type = type; m->fidx = fidx;
   m->id = id;
   m->structure = structure;
   m->next = hset->mtab[hashval]; hset->mtab[hashval] = m;
   if (hset->pmap != NULL) {
      hashval = (unsigned int)m->structure % PTRHASHSIZE;
      p = (PtrMap *)New(hset->hmem,sizeof(PtrMap));
      p->ptr = m->structure; p->m = m; 
      p->next = hset->pmap[hashval]; hset->pmap[hashval] = p;
   }
   return m;
}

/* EXPORT-> DeleteMacro: delete the given macro */
void DeleteMacro(HMMSet *hset, MLink p)
{
   if (p->type == 'h')
      --hset->numPhyHMM;
   else if (p->type == 'l'){
      --hset->numLogHMM;
   }
   else
      --hset->numMacros;
   p->type = '*';
}

/* EXPORT-> DeleteMacroStruct: delete the given macro */
void DeleteMacroStruct(HMMSet *hset, char type, Ptr structure)
{
   MLink p;

   
   p = FindMacroStruct(hset,type,structure);
   if (p==NULL)
      HError(7035,"DeleteMacroStruct: attempt to delete unknown macro");
   if (p->type != type)
      HError(7091,"DeleteMacroStruct: type mismatch %c/%c",p->type,type);
   DeleteMacro(hset,p);
}

/* EXPORT->FindMacroName: find macro def based on given id */
MLink FindMacroName(HMMSet *hset, char type, LabId id)
{
   unsigned int hashval;
   MLink m;
   
   
   if (id == NULL) return NULL;
   hashval = Hash(id->name);
   m = hset->mtab[hashval];
   while (m != NULL && !(m->id == id && m->type == type))
      m = m->next;
   return m;
}

/* EXPORT->FindMacroStruct: return macro for given structure */
MLink FindMacroStruct(HMMSet *hset, char type, Ptr structure)
{
   MLink m;
   int h,n;
   unsigned int i;
   PtrMap *p;
   

   if (hset->pmap == NULL){   /* 1st use so create pmap hash table */
      hset->pmap = (PtrMap **)MakeHashTab(hset,PTRHASHSIZE);
      if (trace&T_PMP) printf("HModel: creating pointer map hash table\n");
      for (n=0,h=0; h<MACHASHSIZE; h++)
         for (m=hset->mtab[h]; m!=NULL; m=m->next){
            i = (unsigned int)m->structure % PTRHASHSIZE;
            p = (PtrMap *)New(hset->hmem,sizeof(PtrMap));
            p->ptr = m->structure; p->m = m; 
            p->next = hset->pmap[i]; hset->pmap[i] = p;
            ++n;
         }
      if (trace&T_PMP) printf("HModel: %d pointers hashed\n",n);
   }
   i = (unsigned int)structure % PTRHASHSIZE;
   for (p = hset->pmap[i]; p != NULL; p = p->next) {
      m = p->m;
      if (p->ptr == structure && m->type == type)
         return m;
   }
   return NULL;
}

/* EXPORT->HasMacros: return true if numMacros>0, if types != NULL
                      create a list of all macro types used */
Boolean HasMacros(HMMSet *hset, char * types)
{
   MLink p;
   char buf[2];
   int h;
   
   if (hset->numMacros == 0) return FALSE;
   if (types!=NULL){
      types[0] = '\0'; buf[1] = '\0';
      for (h=0; h<MACHASHSIZE; h++)
         for (p=hset->mtab[h]; p!=NULL; p=p->next)
            if (strchr(types,p->type) == NULL){
               buf[0] = p->type; 
               strcat(types,buf);
            }
   }
   return TRUE;
}

/* EXPORT->SetVFloor: set vFloor[1..S] from macros "varFloorN" */
void SetVFloor(HMMSet *hset, Vector *vFloor, float minVar)
{
   int j,s,S,size;
   char mac[MAXSTRLEN], num[10];
   LabId id;
   MLink m;
   SVector v;
   
   S = hset->swidth[0];
   for (s=1; s<=S; s++){
      strcpy(mac,"varFloor");
      sprintf(num,"%d",s); strcat(mac,num);
      id = GetLabId(mac,FALSE);
      if (id != NULL  && (m=FindMacroName(hset,'v',id)) != NULL){
         v = (SVector)m->structure;
         SetUse(v,1);
         size = VectorSize(v);
         if (size != hset->swidth[s])
            HError(7023,"SetVFloor: Macro %s has vector size %d, should be %d",
                   mac,size,hset->swidth[s]);
         if (vFloor!=NULL)
            vFloor[s] = v;
      } else if (vFloor!=NULL) {
         if (id != NULL)
            HError(-7023,"SetVFLoor: %s used but not a variance floor macro",mac);
         vFloor[s] = CreateVector(hset->hmem,hset->swidth[s]);
         for (j=1; j<=hset->swidth[s]; j++)
            vFloor[s][j] = minVar;
      }
   }
}

/* EXPORT->ApplyVFloor: apply the variance floors in hset to all mix comps */
void ApplyVFloor(HMMSet *hset)
{
   int s,S,k,vSize;
   int nFloorVar,nFloorVarMix,nMix;
   Boolean mixFloored;
   Vector vFloor[SMAX],minv;
   Covariance cov;
   MLink m;
   LabId id;
   char mac[32];
   HMMScanState hss;
   LogFloat ldet;

   /* get varFloor vectors */
   S=hset->swidth[0];
   for (s=1;s<=S;s++) {
      sprintf(mac,"varFloor%d",s);
      id = GetLabId(mac,FALSE);
      if (id != NULL  && (m=FindMacroName(hset,'v',id)) != NULL) {
         vFloor[s] = (Vector)m->structure;
         vSize = VectorSize(vFloor[s]);
         if (vSize != hset->swidth[s])
            HError(7023,"SetVFloor: Macro %s has vector size %d, should be %d",
                   mac,vSize,hset->swidth[s]);
      }
      else
         HError(7023,"ApplyVFLoor: variance floor macro %s not found",mac);
   }

   /* apply the variance floors to the mixcomps */
   NewHMMScan(hset,&hss);
   nMix = nFloorVar = nFloorVarMix = 0 ;
   do {
      nMix++;
      mixFloored = FALSE;
      cov=hss.mp->cov;
      minv = vFloor[hss.s];
      vSize = VectorSize(minv);
      switch (hss.mp->ckind) {
      case DIAGC: /* diagonal covariance matrix */ 
         for (k=1; k<=vSize; k++){
            if (cov.var[k]<minv[k]) {
               cov.var[k] = minv[k];
               nFloorVar++;
               mixFloored = TRUE;
            }
         }
         FixDiagGConst(hss.mp); 
         break;
      case INVDIAGC: /* inverse diagonal covariance matrix */ 
         for (k=1; k<=vSize; k++){
            if (cov.var[k]> 1.0/minv[k]) {
               cov.var[k] = 1.0/minv[k];
               nFloorVar++;
               mixFloored = TRUE;
            }
         }
         FixInvDiagGConst(hss.mp); 
         break;
      case FULLC: /* full covariance matrix */ 
         CovInvert(cov.inv,cov.inv);
         for (k=1; k<=vSize; k++){
            if (cov.inv[k][k]<minv[k]) {
               cov.inv[k][k] = minv[k];
               nFloorVar++;
               mixFloored = TRUE;
            }
         }
         FixFullGConst(hss.mp, -CovInvert(cov.inv,cov.inv) ); 
         break;
      case LLTC:
      case XFORMC:   
      default:
         HError(7023,"ApplyVFLoor: CovKind not (yet) supported");
      }

      if (mixFloored == TRUE) nFloorVarMix++;
   } while (GoNextMix(&hss,FALSE));
   EndHMMScan(&hss);

   /* summary */
   if (trace&T_TOP) {
      if (nFloorVar > 0)
         printf("ApplyVFloor: Total %d floored variance elements in %d out of %d mixture components\n",
                nFloorVar,nFloorVarMix,nMix);
      fflush(stdout);
   }
}


/* ---------------------- Print Profile Routines ---------------------- */

/* EXPORT->PrintHMMProfile: print out a profile of the given HMM */
void PrintHMMProfile(FILE *f, HLink hmm)
{
   int i,N,s,S;
   char buf[64];
   
   N = hmm->numStates; S = hmm->owner->swidth[0];
   fprintf(f," States   : "); 
   for (i=2;i<N;i++) fprintf(f,"%3d",i); fprintf(f," (width)\n");
   for (s=1;s<=S;s++){
      fprintf(f," Mixes  s%d: ",s); 
      for (i=2;i<N;i++) fprintf(f,"%3d",hmm->svec[i].info->pdf[s].nMix);
      fprintf(f," ( %2d  )\n",hmm->owner->swidth[s]);
      fprintf(f," Num Using: "); 
      for (i=2;i<N;i++) fprintf(f,"%3d",hmm->svec[i].info->nUse); 
      fprintf(f,"\n");
   }
   fprintf(f," Parm Kind:  %s\n",ParmKind2Str(hmm->owner->pkind,buf));
   fprintf(f," Number of owners = %d\n",hmm->nUse);
}

/* GetHashCounts: used for both mtab and pmap */
static void GetHashCounts(MLink * mtab,int size, int *min, int *max,
                          int *numz, int *mean, int *var)
{
   int h,i;
   float sum,sumsq,meen;
   MLink m;
   
   sum = sumsq = 0.0;
   *min = INT_MAX; *max = *numz = 0;
   for (h=0; h<size; h++){
      for (i=0,m=mtab[h]; m!=NULL; m=m->next,i++);
      if (i>*max) *max = i;
      if (i<*min) *min = i;
      if (i==0) ++(*numz);
      sum += i; sumsq += i*i;
   }
   meen = sum/(float)size;
   *mean = (int) meen;
   *var = (int) sqrt(sumsq/(float)size - meen*meen);
}

/* PrintHashUsage: print out usage of given hash table */
static void PrintHashUsage(FILE *f, HMMSet *hset)
{
   int min,max,numz,mean,var;
   
   fprintf(f," HashTab: numBins numZero  maxBin  minBin  meanBin  stdBin\n");
   GetHashCounts(hset->mtab,MACHASHSIZE,&min,&max,&numz,&mean,&var);
   fprintf(f,"   mtab : %6d  %6d  %6d  %6d  %6d  %6d\n",
           MACHASHSIZE,numz,max,min,mean,var);
   if (hset->pmap != NULL){
      GetHashCounts((MLink *)hset->pmap,PTRHASHSIZE,&min,&max,&numz,&mean,&var);
      fprintf(f,"   pmap : %6d  %6d  %6d  %6d  %6d  %6d\n",
              PTRHASHSIZE,numz,max,min,mean,var);
   }
}

/* EXPORT->PrintHSetProfile: Print a profile to f of given HMM set. */
void PrintHSetProfile(FILE *f, HMMSet *hset)
{
   char buf[20];
   MILink p;
   CovKind ck;
   int i,S;
   
   fprintf(f," Memory: "); PrintHeapStats(hset->hmem);
   fprintf(f," Files:");
   if (hset->numFiles > 0) 
      for (p=hset->mmfNames; p!=NULL; p=p->next)
         fprintf(f," %s",p->fName);
   else
      fprintf(f," None");
   fprintf(f,"\n");
   fprintf(f," Macros = %d;  Logical HMMs = %d;  Physical HMMs = %d\n",
           hset->numMacros,hset->numLogHMM,hset->numPhyHMM);
   if (HasMacros(hset,buf))
      fprintf(f," Macro types: %s\n",buf);
   S = hset->swidth[0];
   fprintf(f," VecSize = %d; Streams = %d [",hset->vecSize,S);
   for (i=1; i<=S; i++) 
      fprintf(f,"%d%c",hset->swidth[i],i==S?']':'|'); 
   fprintf(f,"\n Max mix comps = %d; Max states = %d\n",
           MaxMixInSet(hset),MaxStatesInSet(hset)); 
   fprintf(f," Parameter Kind = %s\n",ParmKind2Str(hset->pkind,buf)),
      fprintf(f," Duration  Kind = %s\n",DurKind2Str(hset->dkind,buf));
   fprintf(f," HMM Set   Kind = ");
   switch(hset->hsKind){
   case PLAINHS:  fprintf(f,"PLAIN\n"); break;
   case SHAREDHS: fprintf(f,"SHAREDHS\n"); break;
   case TIEDHS:   fprintf(f,"TIEDHS\n"); break;
   case DISCRETEHS:   fprintf(f,"DISCRETEHS\n"); break;
   }

   fprintf (f, " Global CovKind = %s\n", CovKind2Str(hset->ckind, buf));
   for (ck = 0; ck < NUMCKIND; ck++) 
      if (hset->ckUsage[ck] > 0 ) 
	 fprintf (f, "   CK %-8s = %d mixcomps\n", CovKind2Str (ck, buf), hset->ckUsage[ck]);

   PrintHashUsage(f,hset);
}


/* ------------------- HMM/Macro Load Routines -------------------- */

/* LoadAllMacros: loads macros from MMF file fname */
static ReturnStatus LoadAllMacros(HMMSet *hset, char *fname, short fidx)
{
   Source src;
   Token tok;
   Ptr structure;
   MLink m;
   char type,buf[MAXSTRLEN];
   LabId id;
   HLink dhmm;
   HMMSet dset;
   int nState=0;

   if (trace&T_MAC)
      printf("HModel: getting Macros from %s\n",fname);
   if(InitScanner(fname,&src,&tok,hset)<SUCCESS){
      HRError(7010,"LoadAllMacros: Can't open file");
      return(FAIL);
   }

   if(GetToken(&src,&tok)<SUCCESS){
      TermScanner(&src);
      HMError(&src,"LoadAllMacros: GetToken failed");
      return(FAIL);
   }
   while (tok.sym != EOFSYM){
      if (tok.sym != MACRO){
         TermScanner(&src);
         HMError(&src,"LoadAllMacros: Macro sym expected");
         return(FAIL);
      }
      type = tok.macroType;
      if (type == 'o') {    /* load an options macro */
         if(GetOptions(hset,&src,&tok,&nState)==FAIL){
            TermScanner(&src);
            HMError(&src,"LoadAllMacros: GetOptions Failed");
            return(FAIL);
         }
      }
      else {   
         if (!ReadString(&src,buf)){
            TermScanner(&src);
            HMError(&src,"LoadAllMacros: Macro name expected");
            return(FAIL);
         }
         id = GetLabId(buf,TRUE);
         if (type == 'h'){    /* load a HMM definition */
            m = FindMacroName(hset,'h',id);
            if (m == NULL) {
               if (!allowOthers){
                  TermScanner(&src);
                  HRError(7030,"LoadAllMacros: phys HMM %s unexpected in %s",
                          id->name,fname);
                  return(FAIL);
               }
               dset=*hset;
               dset.hmem=&gstack;
               dhmm = (HLink) New(&gstack,sizeof(HMMDef));
               dhmm->owner=NULL; dhmm->numStates=0; dhmm->nUse=0; dhmm->hook=NULL;
               if (trace&T_MAC)
                  printf("HModel: skipping HMM Def from macro %s\n",id->name);
               if(GetToken(&src,&tok)<SUCCESS){
                  TermScanner(&src);
                  Dispose(&gstack,dhmm);
                  HMError(&src,"LoadAllMacros: GetToken failed");       
                  return(FAIL);
               }               
               if(GetHMMDef(&dset,&src,&tok,dhmm,nState)<SUCCESS){
                  TermScanner(&src);
                  Dispose(&gstack,dhmm);
                  HMError(&src,"LoadAllMacros: GetHMMDef failed");       
                  return(FAIL);
               }               
               Dispose(&gstack,dhmm);
            }
            else {
               if (trace&T_MAC)
                  printf("HModel: getting HMM Def from macro %s\n",m->id->name);
               if(GetToken(&src,&tok)<SUCCESS){
                  TermScanner(&src);
                  HMError(&src,"LoadAllMacros: GetToken failed");
                  return(FAIL);
               }
               if(GetHMMDef(hset,&src,&tok,(HLink)m->structure,nState)<SUCCESS){
                  TermScanner(&src);
                  HMError(&src,"LoadAllMacros: GetHMMDef failed");       
                  return(FAIL);
               }               
               m->fidx = fidx;
            }
         } else {             /* load a shared structure */
            /* input transforms are store prior to the options !! */
            if ((type != 'j') && (CheckOptions(hset)<SUCCESS)) {
               TermScanner(&src);
               HMError(&src,"LoadAllMacros: CheckOptions failed");
               return(FAIL);
            }
            if(GetToken(&src,&tok)<SUCCESS){
               TermScanner(&src);
               HMError(&src,"LoadAllMacros: GetToken failed");
               return(FAIL);
            }
            switch(type){
            case 's': structure = GetStateInfo(hset,&src,&tok); break;
            case 'm': structure = GetMixPDF(hset,&src,&tok);    break;
            case 'u': structure = GetMean(hset,&src,&tok);      break;
            case 'v': structure = GetVariance(hset,&src,&tok);  break;
            case 'i': structure = GetCovar(hset,&src,&tok);     break;
            case 'c': structure = GetCovar(hset,&src,&tok);     break;
            case 'x': structure = GetTransform(hset,&src,&tok); break;
            case 'w': structure = GetSWeights(hset,&src,&tok);  break;
            case 't': structure = GetTransMat(hset,&src,&tok);  break;
            case 'd': structure = GetDuration(hset,&src,&tok);  break;
            case 'r': structure = GetRegTree(hset,&src,&tok);   break;
               /* code for input transform support */
            case 'f': structure = GetLinXForm(hset,&src,&tok);  break;
            case 'y': structure = GetBias(hset,&src,&tok);   break;
            case 'j': 
               structure = GetInputXForm(hset,&src,&tok);   
               ((InputXForm *)structure)->xformName = CopyString(hset->hmem,id->name);
               break;
            default :
               TermScanner(&src);
               HRError(7037,"LoadAllMacros: bad macro type in MMF %s",fname);
               return(FAIL);
            }
            if(structure==NULL){
               TermScanner(&src);
               HRError(7035,"LoadAllMacros: Get macro data failed in MMF %s",fname);
               return(FAIL);
            }
            NewMacro(hset,fidx,type,id,structure);
            if (trace&T_MAC)
               printf("HModel: storing macro ~%c %s -> %p\n",type,id->name,structure);
         }
      }
   }
   TermScanner(&src);
   return(SUCCESS);
}

/* LoadMacroFiles: scan file list of hset and load any macro files */
static ReturnStatus LoadMacroFiles(HMMSet *hset)
{
   MILink mmf;
   int i=0;
   ReturnStatus result=SUCCESS;
   
   for (mmf=hset->mmfNames; mmf!=NULL; mmf=mmf->next)
      if (!mmf->isLoaded){
         if(LoadAllMacros(hset,mmf->fName,++i)<SUCCESS) result=FAIL;
         mmf->isLoaded = TRUE;
      }
   return result;
}

/* IsShared: decide if given HMM set is SHAREDHS.  Note that shared means,
      vars, invs, xforms dont make it a SHARED system since there is no
      way to exploit the sharing at this level. */
static Boolean IsShared(HMMSet *hset)
{
   char types[20];
   
   if (HasMacros(hset,types)) {
      if (strchr(types,'m') != NULL ||  strchr(types,'s') != NULL)
         return TRUE;
   }
   return FALSE;
}

/* ConcatFN: make up file name */
static char *ConcatFN(char *path, char *base, char *ext, char *fname)
{
   char *s;

   fname[0] = '\0';
   if (path!=NULL && *path!='\0') {
      if ((s=strrchr(path,PATHCHAR))!=NULL && *(s+1)=='\0')
         strcat(fname,path);
      else
         sprintf(fname,"%s%c",path,PATHCHAR);
   }
   strcat(fname,base);
   if (ext!=NULL && *ext!='\0') {
      strcat(fname,".");
      strcat(fname,ext);
   }
   return fname;
}


void SetIndexes(HMMSet *hset)
{
   HMMScanState hss;
   StateInfo *si;
   MixPDF *mp;
   MLink m;
   int h,nm,nsm,ns,nss,nt;
   
   /* Reset indexes */
   nt=0;
   NewHMMScan(hset,&hss);
   while(GoNextState(&hss,FALSE))
      hss.si->sIdx=-1;
   EndHMMScan(&hss);
   if (hset->hsKind == PLAINHS || hset->hsKind == SHAREDHS) {
      NewHMMScan(hset,&hss);
      while(GoNextMix(&hss,FALSE))
         hss.mp->mIdx=-1;
      EndHMMScan(&hss);
   }

   NewHMMScan(hset,&hss);
   do {
      if (!IsSeenV(hss.hmm->transP)) {
         SetHook(hss.hmm->transP,(Ptr)(++nt));
         TouchV(hss.hmm->transP);
      }
      hss.hmm->tIdx=(int)GetHook(hss.hmm->transP);
   }
   while(GoNextHMM(&hss));
   EndHMMScan(&hss);
   NewHMMScan(hset,&hss);
   do {
      SetHook(hss.hmm->transP,NULL);
      UntouchV(hss.hmm->transP); /* Not needed (done by EndHMMScan) */
   }
   while(GoNextHMM(&hss));
   EndHMMScan(&hss);

   nsm=nss=0;
   for (h=0; h<MACHASHSIZE; h++)
      for (m=hset->mtab[h]; m!=NULL; m=m->next) {
         if (m->type=='s') {
            si=(StateInfo *)m->structure;
            if (si->sIdx<0) si->sIdx=++nss;
         }
         if (m->type=='m') {
            mp=(MixPDF *)m->structure;
            if (mp->mIdx<0) mp->mIdx=++nsm;
         }
      }
   nm=nsm;ns=nss;
   NewHMMScan(hset,&hss);
   while(GoNextState(&hss,FALSE))
      if (hss.si->sIdx<0) hss.si->sIdx=++ns;
   EndHMMScan(&hss);
   if (hset->hsKind == PLAINHS || hset->hsKind == SHAREDHS) {
      NewHMMScan(hset,&hss);
      while(GoNextMix(&hss,FALSE))
         if (hss.mp->mIdx<0) hss.mp->mIdx=++nm;
      EndHMMScan(&hss);
   }
   hset->numStates=ns; hset->numSharedStates=nss;
   hset->numMix=nm; hset->numSharedMix=nsm;
   hset->numTransP=nt;
}


/* SetCovKindUsage

     Count number of mixtures for using each covKind.
     Set covKind at HMMSet level, if it is currently NULLC and it is the same for
     all mixtures.
*/
void SetCovKindUsage (HMMSet *hset)
{
   HMMScanState hss;
   CovKind ck;

   for (ck = 0; ck < NUMCKIND; ck++)
      hset->ckUsage[ck] =0;

   if (hset->hsKind == PLAINHS || hset->hsKind == SHAREDHS) {
      NewHMMScan(hset,&hss);
      while(GoNextMix(&hss,FALSE))
         ++hset->ckUsage[hss.mp->ckind];
      EndHMMScan(&hss);

      /* set global covKind if currently unset (i.e. NULLC) */
      if (hset->ckind == NULLC) {
         CovKind lastSeenCK;
         int nCK = 0;

         for (ck = 0; ck < NUMCKIND; ck++)
            if (hset->ckUsage[ck] > 0) {
               ++nCK;
               lastSeenCK = ck;
            }
         if (nCK == 1)
            hset->ckind = lastSeenCK;
      }
   }
}

/*ResetHMMSet. New function ResetHMMSet - both hash tables removed and values of numLogHMM etc. set to zero. Dispose() done assuming MSTAK */
void ResetHMMSet(HMMSet *hset)
{
   int i;

   for(i=0; i<MACHASHSIZE; i++)
      hset->mtab[i]=NULL;

   for(i=0; hset->pmap!=NULL && i<PTRHASHSIZE; i++)
      hset->pmap[i]=NULL;
   hset->numLogHMM=0;
   hset->numPhyHMM=0;
   hset->numMacros=0;
   hset->numFiles=0;
   hset->mmfNames=NULL;
   Dispose(hset->hmem, hset->firstElem);
}

   
/* EXPORT->LoadHMMSet: Load all definition files for given hset */
ReturnStatus LoadHMMSet(HMMSet *hset, char *hmmDir, char *hmmExt)
{
   int h,nState=0;
   MLink p;
   char fname[MAXSTRLEN],buf[MAXSTRLEN];
   HLink hmm;
   Source src;
   Token tok;
   
   hset->hsKind = PLAINHS; /* default assumption */
   if(LoadMacroFiles(hset)<SUCCESS){
      HRError(7050,"LoadHMMSet: Macro name expected");         
      ResetHMMSet(hset);
      return(FAIL);
   }
   for (h=0; h<MACHASHSIZE; h++)
      for (p=hset->mtab[h]; p!=NULL; p=p->next) 
         if (p->type == 'h'){
            hmm = (HLink)p->structure;
            if (hmm->numStates == 0 ) {
               ConcatFN(hmmDir,p->id->name,hmmExt,fname); 
               if(InitScanner(fname,&src,&tok,hset)<SUCCESS){
                  HRError(7010,"LoadHMMSet: Can't find file");         
                  ResetHMMSet(hset);
                  return(FAIL);
               }
               if (trace&T_MAC)
                  printf("HModel: getting HMM Def from %s\n",fname);
               if(GetToken(&src,&tok)<SUCCESS){
                  TermScanner(&src);
                  ResetHMMSet(hset);
                  HMError(&src,"LoadHMMSet: GetToken failed");
                  return(FAIL);
               }
               while (tok.sym == MACRO)
                  switch (tok.macroType){
                  case 'o':
                     if(GetOptions(hset,&src,&tok, &nState)==FAIL){
                        TermScanner(&src);
                        ResetHMMSet(hset);
                        HMError(&src,"LoadHMMSet: GetOptions failed");
                        return(FAIL);
                     }
                     break;
                  case 'h':
                     if (!ReadString(&src,buf)){
                        TermScanner(&src);
                        ResetHMMSet(hset);
                        HMError(&src,"LoadHMMSet: Macro name failed");
                        return(FAIL);
                     }
                     if (GetLabId(buf,FALSE) != p->id){
                        TermScanner(&src);
                        ResetHMMSet(hset);
                        HMError(&src,"LoadHMMSet: Inconsistent HMM macro name");
                        return(FAIL);
                     }
                     if(GetToken(&src,&tok)<SUCCESS){
                        TermScanner(&src);
                        ResetHMMSet(hset);
                        HMError(&src,"LoadAllMacros: GetToken failed");
                        return(FAIL);
                     }
                     break;
                  default:
                     TermScanner(&src);
                     ResetHMMSet(hset);
                     HMError(&src,"LoadHMMSet: Unexpected macro in HMM def file");
                     return(FAIL);
                     break;
                  }              
               if(GetHMMDef(hset,&src,&tok,hmm,nState)<SUCCESS){
                  TermScanner(&src);
                  ResetHMMSet(hset);
                  HRError(7032,"LoadHMMSet: GetHMMDef failed");
                  return(FAIL);
               }
               TermScanner(&src);
            }
         }
   if (forceHSKind) {
      if (hset->hsKind == DISCRETEHS && cfHSKind != DISCRETEHS){
         HRError(7032,"LoadHMMSet: cannot change DISCRETE HMM Kind");
         ResetHMMSet(hset);
         return(FAIL);
      }
      if (hset->hsKind == TIEDHS && cfHSKind != TIEDHS){
         HRError(7032,"LoadHMMSet: cannot change TIED HMM Kind");
         ResetHMMSet(hset);
         return(FAIL);
      }
      hset->hsKind = cfHSKind;
   }else{
      if (hset->hsKind == PLAINHS && IsShared(hset))
         hset->hsKind = SHAREDHS;
   }
   if (checking) 
      if (CheckHSet(hset)<SUCCESS){
         ResetHMMSet(hset);
         HRError(7031,"LoadHMMSet: Invalid HMM data");
         return(FAIL);
      }
   SetIndexes(hset);
   SetCovKindUsage(hset);
   return(SUCCESS);
}

/* ------------------------ HMM Set Creation ----------------------- */

/* EXPORT->CreateHMMSet: create the basic HMMSet structure */
void CreateHMMSet(HMMSet *hset, MemHeap *heap, Boolean allowTMods)
{
   int s;

   /* set default values in hset structure */
   hset->hmem = heap;
   hset->hmmSetId = NULL;
   hset->mmfNames = NULL; hset->numFiles = 0;
   hset->allowTMods = allowTMods;  hset->optSet = FALSE;
   hset->vecSize = 0; hset->swidth[0] = 0;
   hset->dkind = NULLD; hset->ckind = NULLC; hset->pkind = 0;
   hset->numPhyHMM = hset->numLogHMM = hset->numMacros = 0;
   hset->xf = NULL;
   for (s=1; s<SMAX; s++) {
      hset->tmRecs[s].nMix = 0; hset->tmRecs[s].mixId = NULL;
      hset->tmRecs[s].probs = NULL; hset->tmRecs[s].mixes = NULL;
   }
   /* initialise the hash tables */
   hset->mtab = (MLink *)MakeHashTab(hset,MACHASHSIZE);
   hset->pmap = NULL;
}

/* CreateHMM: create logical macro. If pId is unknown, create macro for
              it too, along with HMMDef.  */
static ReturnStatus CreateHMM(HMMSet *hset, LabId lId, LabId pId)
{
   MLink m;
   HLink hmm;
   Boolean newMacro=FALSE; /* for memory clear up*/

   m = FindMacroName(hset,'l',lId);
   if (m != NULL){
      HRError(7036,"CreateHMM: multiple use of logical HMM name %s",lId->name);
      return(FAIL);
   }

   m = FindMacroName(hset,'h',pId);
   if (m == NULL) {  /* need to create new phys macro and HMMDef */
      hmm = (HLink)New(hset->hmem,sizeof(HMMDef));
      hmm->owner = hset; hmm->nUse = 1; hmm->hook = NULL;
      hmm->numStates = 0;  /* indicates HMMDef not yet defined */
      if((m=NewMacro(hset,0,'h',pId,hmm))==NULL){
         HRError(7091,"CreateHMM: NewMacro (Physical) failed"); /*will never happen*/
         return(FAIL);
      } 
      newMacro=TRUE;
      if (pId != lId) ++hmm->nUse;
   } 
   else {
      if (m->type != 'h'){
         HRError(7091,"CreateHMM: %s is not a physical HMM",pId->name);
         return(FAIL);
      }
      hmm = (HLink)m->structure;
      ++hmm->nUse;
   }
   if(NewMacro(hset,0,'l',lId,hmm)==NULL){ 
      HRError(7091,"CreateHMM: NewMacro (Logical) failed"); /*will never happen*/
      return(FAIL);
   }    
   return(SUCCESS);
}


/* InitHMMSet: Init a HMM set by reading the HMM list in fname.
               If isSingle, then fname is the name of a single HMM */
static ReturnStatus InitHMMSet(HMMSet *hset, char *fname, Boolean isSingle)
{
   Source src;
   char buf[MAXSTRLEN];
   LabId lId, pId;
    
   /* sets first element on heap to allow disposing of memory */
   hset->firstElem = (Boolean *) New(hset->hmem, sizeof(Boolean));

   if (isSingle){
      /* fname is a single HMM file to load as a singleton set */
      lId = GetLabId(fname,TRUE);
      pId = lId;
      if(CreateHMM(hset,lId,pId)<SUCCESS){
         HRError(7060,"InitHMMSet: Error in CreateHMM", hset->numLogHMM);
         return(FAIL);   
      }
   } 
   else {
      /* read the HMM list file and build the logical list */
      if(InitSource(fname,&src,HMMListFilter)<SUCCESS){
         HRError(7010,"InitHMMSet: Can't open list file %s", fname);           
         return(FAIL);
      }
      SkipWhiteSpace(&src);
      while(ReadString(&src,buf)) {
         lId = GetLabId(buf,TRUE);
         SkipWhiteSpace(&src);
         if (!src.wasNewline) {
            if (!ReadString(&src,buf)){
               CloseSource(&src);
               HRError(7060,"InitHMMSet: Expected a physical name for %d'th HMM", hset->numLogHMM);            
               return(FAIL);
            }
            pId = GetLabId(buf,TRUE);
            SkipWhiteSpace(&src);
         } else               /* physical name same as logical name */
            pId = lId;
         if(CreateHMM(hset,lId,pId)<SUCCESS){
            CloseSource(&src);
            HRError(7060,"InitHMMSet: Error in CreateHMM", hset->numLogHMM);
            return(FAIL);        
         }
         if (!src.wasNewline){
            CloseSource(&src);
            HRError(7060,"InitHMMSet: Expected newline after %d'th HMM",
                    hset->numLogHMM);
            return(FAIL);
         }
      }
      CloseSource(&src);
   }
   return(SUCCESS);
}


/* EXPORT->MakeHMMSet: Make a HMM set by reading the HMM list in fname */
ReturnStatus MakeHMMSet(HMMSet *hset, char *fname)
{
   if(InitHMMSet(hset, fname, FALSE)<SUCCESS){
      ResetHMMSet(hset);
      return(FAIL);
   }
   return(SUCCESS);
}

/* EXPORT->MakeOneHMM: Create a singleton for the HMM hname */
ReturnStatus MakeOneHMM(HMMSet *hset, char *hname)
{
   if(InitHMMSet(hset, hname, TRUE)<SUCCESS){
      ResetHMMSet(hset);
      return(FAIL);
   }
   return(SUCCESS);
}

/* -------------------- HMM/Macro Save Routines -------------------- */

/* SaveMacros: save all shared macro structures associated with fidx */
static void SaveMacros(FILE *f, HMMSet *hset, short fidx, Boolean binary)
{
   MLink m;
   int h;
      
   /* 
      First thing is to store the InputXForm if necessary. Rules are:
      1) If read from Macro file store using related file
      2) If Loaded using a config variable/Created  store in the first macrofile
   */
   if (saveInputXForm) { 
      for (h=0; h<MACHASHSIZE; h++)
         for (m=hset->mtab[h]; m!=NULL; m=m->next)
            if ((m->type == 'j') && (((InputXForm *)m->structure)->nUse>0) &&
		(((fidx == 1) && ((m->fidx == LOADFIDX) || (m->fidx == CREATEFIDX))) ||
		 (m->fidx == fidx)))
               PutInputXForm(hset,f,m,(InputXForm *)m->structure,TRUE,binary);     
   }
   fprintf(f,"~o\n");
   PutOptions(hset,f,binary);
   for (h=0; h<MACHASHSIZE; h++)
      for (m=hset->mtab[h]; m!=NULL; m=m->next)
         if (m->fidx == fidx)
            switch(m->type){            /* atomic macros first */
            case 'u':
               PutMean(hset,f,m,(SVector)m->structure,TRUE,binary);     
               break;
            case 'v':
               PutVariance(hset,f,m,(SVector)m->structure,TRUE,binary); 
               break;
            case 'i':
               PutCovar(hset,f,m,(SMatrix)m->structure,INVCOVAR,TRUE,binary); 
               break;
            case 'c':
               PutCovar(hset,f,m,(SMatrix)m->structure,LLTCOVAR,TRUE,binary); 
               break;
            case 'x':
               PutTransform(hset,f,m,(SMatrix)m->structure,TRUE,binary);
               break;
            case 'w':
               PutSWeights(hset,f,m,(SVector)m->structure,TRUE,binary); 
               break;
            case 'd':
               PutDuration(hset,f,m,(SVector)m->structure,TRUE,binary); 
               break;
            case 't':
               PutTransMat(hset,f,m,(SMatrix)m->structure,TRUE,binary);
               break;
            case 'r':
               if (saveRegClass)
                  PutRegTree(hset,f,m,(RegTree *)m->structure,TRUE,binary);
               break;
            default: break;
            }
   for (h=0; h<MACHASHSIZE; h++)    /* then higher levels */
      for (m=hset->mtab[h]; m!=NULL; m=m->next)
         if (m->fidx == fidx && m->type == 'm')
            PutMixPDF(hset,f,m,(MixPDF *)m->structure,TRUE,binary);
   for (h=0; h<MACHASHSIZE; h++)    /* then higher levels */
      for (m=hset->mtab[h]; m!=NULL; m=m->next)
         if (m->fidx == fidx && m->type == 's')
            PutStateInfo(hset,f,m,(StateInfo *)m->structure,TRUE,binary);
   for (h=0; h<MACHASHSIZE; h++)    /* then top level */
      for (m=hset->mtab[h]; m!=NULL; m=m->next)
         if (m->fidx == fidx && m->type == 'h')
            PutHMMDef(hset,f,m,TRUE,binary); 
}

static ReturnStatus GetXFormMacros(HMMSet *hset, Source *src, Token *tok, int fidx)
{
   char type, buf[MAXSTRLEN];
   LabId id;
   Ptr structure;

   while (tok->sym == MACRO){
      type = tok->macroType;
      if (!ReadString(src,buf)){
         TermScanner(src);
         HRError(999,"GetXFormMacros: Macro name expected in macro file %s",src->name);
         return(FAIL);
      }
      id = GetLabId(buf,TRUE);
      if(GetToken(src,tok)<SUCCESS){
         TermScanner(src);
         HRError(999,"GetXFormMacros: Macro name expected in macro file %s",src->name);
         return(FAIL);
      }
      switch(type){
      case 'f': structure = GetLinXForm(hset,src,tok);  break;
      case 'x': structure = GetTransform(hset,src,tok); break;
      case 'y': structure = GetBias(hset,src,tok);   break;
      case 'j': 
         structure = GetInputXForm(hset,src,tok);
         ((InputXForm *)structure)->xformName = CopyString(hset->hmem,id->name);
         break;
      default :
         TermScanner(src);
         HRError(7037,"GetXFormMacros: bad macro type in xform file %s",src->name);
         return(FAIL);
      }
      if(structure==NULL){
         TermScanner(src);
         HRError(7035,"GetXFormMacros: Get macro data failed in %s",src->name);
         return(FAIL);
      }
      NewMacro(hset,fidx,type,id,structure);
      if (trace&T_MAC)
         printf("HModel: storing macro ~%c %s -> %p\n",type,id->name,structure);
   }
   return(SUCCESS);
}

/* EXPORT->SaveInOneFile: ignore source files and store in fname */
void SaveInOneFile(HMMSet *hset, char *fname)
{
   int h;
   MLink m; 
   MILink p;

   /* unset all isLoaded flags to preserve original files */
   for (p=hset->mmfNames; p!=NULL; p=p->next)
      p->isLoaded = FALSE;
   /* Add new file and put its fidx in all loaded macros */
   p = AddMMF(hset,fname); p->isLoaded = TRUE;
   for (h=0; h<MACHASHSIZE; h++)
      for (m=hset->mtab[h]; m!=NULL; m=m->next)
         if (m->type != 'l')
            m->fidx = hset->numFiles;
}

/* FixOrphanMacros: make sure orphan macros get output somewhere. 
   This is a real hack and cannot guarantee that resulting file
   set can ever be reloaded - only alternative is for user to 
   store as single MMF. */
void FixOrphanMacros(HMMSet *hset)
{
   int h,numMacs=0,numHMMs = 0;
   short firstx,firsth,firsts,firstm;
   MLink m;
   MILink p;

   /* find first occ of major macro classes */
   firsth=firsts=firstm=firstx=hset->numFiles+1;
   for (h=0; h<MACHASHSIZE; h++)
      for (m=hset->mtab[h]; m!=NULL; m=m->next)
         if (m->fidx == 0)
            switch(m->type){
            case 'h':  ++numHMMs; break;
            case 'l':
            case 'o': break;
            default: ++numMacs; break;
            }   
         else
            switch (m->type){
            case 'h': 
               if (m->fidx<firsth) firsth=m->fidx;
               break;
            case 's':
               if (m->fidx<firsts) firsts=m->fidx;
               break;
            case 'm':
               if (m->fidx<firstm) firstm=m->fidx;
               break;
            }
   if ((numHMMs>1 && !keepDistinct) || numMacs>0){
      if (trace&T_ORP)
         printf("FixOrphanMacros: %d mac and %d hmm orphans to home\n",
                numMacs,numHMMs);
      if (hset->numFiles==0){
         p = AddMMF(hset,orphanMacFile); p->isLoaded = TRUE;
         if (trace&T_ORP)
            printf("   creating MMF %s\n",orphanMacFile);
      }
      /* now the horrid heuristic to */
      if (firsth>hset->numFiles) firsth=hset->numFiles;
      if (firsts>firsth) {
         firsts=firsth;
         if (firsts>2) --firsts;
      }
      if (firstm>firsts) firstm = firsts;
      if (firstx>firstm) firstx = firstm;
      if (trace&T_ORP)
         printf("  files=%d h=%d,s=%d,m=%d,x=%d\n",
                hset->numFiles,firsth,firsts,firstm,firstx);
      /* finally fix the fidx's */
      for (h=0; h<MACHASHSIZE; h++)
         for (m=hset->mtab[h]; m!=NULL; m=m->next)
            if (m->fidx==0){
               switch (m->type){
               case 'l':
               case '*': 
                  break;
               case 'h': 
                  if (!keepDistinct) m->fidx = firsth; 
                  break;
               case 's':
                  m->fidx = firsts; break;
               case 'm':
                  m->fidx = firstm; break;
               case 'o':
                  m->fidx = 1; break;
               default:
                  m->fidx = firstx; break;
               }
               if (trace&T_ORP && m->type != 'l' && m->type != '*')
                  printf("   %s[%c] to %d\n",m->id->name,m->type,m->fidx);
            }
   }
}

/* EXPORT->SaveHMMSet: save the given HMM set */
ReturnStatus SaveHMMSet(HMMSet *hset, char *hmmDir, char *hmmExt, Boolean binary)
{
   FILE *f;
   MILink p;
   char fname[256];
   int h,i;
   MLink m;
   Boolean isPipe;
   
   FixOrphanMacros(hset);
   binary = binary || saveBinary;
   /* First output to all named MMF files */
   for (p=hset->mmfNames,i=1; p!=NULL; p=p->next,i++) 
      if (p->isLoaded) {
         MakeFN(p->fName,hmmDir,NULL,fname);
         if ((f=FOpen(fname,HMMDefOFilter,&isPipe)) == NULL){
            HRError(7011,"SaveHMMSet: Cannot create MMF file %s",fname);
            return(FAIL);
         }
         if (trace&T_MAC)
            printf("HModel: saving Macros to %s\n",fname);
         SaveMacros(f,hset,i,binary);
         FClose(f,isPipe);
      }

   /* Second output all individual HMMs */
   for (h=0; h<MACHASHSIZE; h++)
      for (m=hset->mtab[h]; m!=NULL; m=m->next)
         if (m->type == 'h' && m->fidx == 0) {
            MakeFN(m->id->name,hmmDir,hmmExt,fname);
            if ((f=FOpen(fname,HMMDefOFilter,&isPipe)) == NULL){
               HRError(7011,"SaveHMMSet: Cannot create HMM file %s",fname);
               return(FAIL);
            }
            if (trace&T_MAC)
               printf("HModel: saving HMM Def to %s\n",fname);
            if (saveGlobOpts) {
               fprintf(f,"~o\n");
               PutOptions(hset,f,binary);
            }
            PutHMMDef(hset,f,m,TRUE,binary);
            FClose(f,isPipe);
         }
   return(SUCCESS);
}

/* EXPORT->SaveHMMList: Save a HMM list in fname describing given HMM set */
ReturnStatus SaveHMMList(HMMSet *hset, char *fname)
{
   int h;
   MLink m,p;
   HLink hmm;
   FILE *f;
   Boolean isPipe;

   if ((f=FOpen(fname,HMMListOFilter,&isPipe)) == NULL){
      HRError(7011,"SaveHMMList: Cannot create HMM list file %s",fname);
      return(FAIL);
   }
   for (h=0; h<MACHASHSIZE; h++)
      for (m=hset->mtab[h]; m!=NULL; m=m->next)
         if (m->type == 'l'){
            fprintf(f,"%s",ReWriteString(m->id->name,NULL,ESCAPE_CHAR));
            hmm = (HLink)m->structure;
            if (hmm->nUse > 0) {
               p = FindMacroStruct(hset,'h',hmm);
               if (p==NULL){
                  HRError(7020,"SaveHMMList: no phys hmm for %s",m->id->name);
                  return(FAIL);
               }
               if (p->id != m->id)
                  fprintf(f," %s",ReWriteString(p->id->name,NULL,ESCAPE_CHAR));
            }
            fprintf(f,"\n");
         }
   FClose(f,isPipe);
   return(SUCCESS);
}

/* -------------- Shared Structure "Seen" Flags -------------- */

/* EXPORT->IsSeen: return true if flag is "set" */
Boolean IsSeen(int flag)
{
   return (flag<0);
}

/* EXPORT->Touch: set given flag */
void Touch(int *flag)
{
   if (*flag==0) 
      *flag = INT_MIN; 
   else if (*flag>0)
      *flag = - *flag;  
}

/* EXPORT->Untouch: set given flag */
void Untouch(int *flag)
{
   if (*flag==INT_MIN) 
      *flag = 0; 
   else if (*flag<0)
      *flag = - *flag;
}

/* EXPORT->ClearStreams: clear flags in given stream */
void ClearStreams(HMMSet *hset, StreamElem *ste, ClearDepth depth)
{
   MixPDF *mp;
   int m;
   
   Untouch(&ste->nMix);
   if (depth==CLR_ALL && (hset->hsKind==PLAINHS || hset->hsKind==SHAREDHS)) {
      for (m=1; m<=ste->nMix; m++){
         mp = ste->spdf.cpdf[m].mpdf;
         Untouch(&mp->nUse);
         UntouchV(mp->mean);
         UntouchV(mp->cov.var);  
      }
   }
}

/* EXPORT->ClearSeenFlags: clear all seen flags downto given depth */
void ClearSeenFlags(HMMSet *hset, ClearDepth depth)
{
   int S,h,s,i;
   MLink m;
   HLink hmm;
   StateElem *se;
   StreamElem *ste;
   StateInfo *si;
   
   S = hset->swidth[0];
   for (h=0; h<MACHASHSIZE; h++)
      for (m=hset->mtab[h]; m!=NULL; m=m->next)
         if (m->type == 'h') {
            hmm = (HLink)m->structure;
            if (hmm->dur != NULL) UntouchV(hmm->dur);
            if (hmm->transP != NULL) UntouchV(hmm->transP);
            if (depth>=CLR_STATES)
               for (i=2,se=hmm->svec+2; i<hmm->numStates; i++,se++){
                  si = se->info;
                  Untouch(&si->nUse);
                  if (si->weights != NULL) UntouchV(si->weights);
                  if (si->dur != NULL) UntouchV(si->dur);
                  if (depth>=CLR_STREAMS)
                     for (s=1,ste=si->pdf+1; s<=S; s++,ste++)
                        ClearStreams(hset,ste,depth);
               }
         }
}

/* ----------------- General HMM Operations --------------------- */

/* EXPORT->MaxMixtures: return max num mixtures in given hmm */
int MaxMixtures(HLink hmm)
{
   int i,s,S,maxM = 0;
   StreamElem *se;
   
   S = hmm->owner->swidth[0];
   for (i=2; i<hmm->numStates; i++){
      se = hmm->svec[i].info->pdf +1;
      for (s=1; s<=S; s++,se++)
         if (se->nMix > maxM) maxM = se->nMix;
   }
   return maxM;
}

/* EXPORT->MaxMixInS: return max num mixes in given stream of given hmm */
int MaxMixInS(HLink hmm, int s)
{
   int i,maxM = 0;
   StreamElem *se;
   
   for (i=2; i<hmm->numStates; i++){
      se = hmm->svec[i].info->pdf +s;
      if (se->nMix > maxM) maxM = se->nMix;
   }
   return maxM;
}

/* EXPORT->MaxMixInSetS: rtn max num mixes in given stream of given hmmset */
int MaxMixInSetS(HMMSet *hset, int s)
{
   int max=0,m,h;
   HLink hmm;
   MLink q;
   Boolean found;

   switch(hset->hsKind){
   case PLAINHS:
   case SHAREDHS:
      for (h=0; h<MACHASHSIZE; h++)
         for (q=hset->mtab[h]; q!=NULL; q=q->next)
            if (q->type == 'h'){
               hmm = (HLink)q->structure;
               m = MaxMixInS(hmm, s);
               if (m>max) max = m;
            }
      break;
   case TIEDHS:
      max=hset->tmRecs[s].nMix;
      break;
   case DISCRETEHS:
      found=FALSE; h=0;
      while((h<MACHASHSIZE) && !found){
         q=hset->mtab[h];
         while((q!=NULL) && !found){
            if (q->type == 'h'){
               hmm = (HLink)q->structure;
               max = MaxMixInS(hmm, s);
               found=TRUE;
            }
            q=q->next;
         }
         h++;
      }
      break;
   }
   return max;     
}


/* EXPORT->MaxMixInSet: return max num mixtures in given set */
int MaxMixInSet(HMMSet *hset)
{
   int max=0,m,h,s;
   HLink hmm;
   MLink q;
   Boolean found;

   switch(hset->hsKind){
   case PLAINHS:
   case SHAREDHS:
      for (h=0; h<MACHASHSIZE; h++)
         for (q=hset->mtab[h]; q!=NULL; q=q->next)
            if (q->type == 'h'){
               hmm = (HLink)q->structure;
               m = MaxMixtures(hmm);
               if (m>max) max = m;
            }
      break;
   case TIEDHS:
      for (s=1;s<=hset->swidth[0];s++){
         m=hset->tmRecs[s].nMix;
         if (m>max) max = m;
      }
      break;
   case DISCRETEHS:
      found=FALSE; h=0;
      while((h<MACHASHSIZE) && !found){
         q=hset->mtab[h];
         while((q!=NULL) && !found){
            if (q->type == 'h'){
               hmm = (HLink)q->structure;
               max = MaxMixtures(hmm);
               found=TRUE;
            }
            q=q->next;
         }
         h++;
      }
      break;
   }
   return max;     
}

/* EXPORT->MaxStatesInSet: return max num states in given set */
int MaxStatesInSet(HMMSet *hset)
{
   int max=0,n,h;
   HLink hmm;
   MLink q;
   
   for (h=0; h<MACHASHSIZE; h++)
      for (q=hset->mtab[h]; q!=NULL; q=q->next)
         if (q->type == 'h'){
            hmm = (HLink)q->structure;
            n = hmm->numStates;
            if (n>max) max = n;
         }
   return max;
}

/* -------------------- Input Transform Operations ---------------- */

/* EXPORT->LoadInputXForm: loads, or returns, the specified transform */
InputXForm *LoadInputXForm(HMMSet *hset, char* macroname, char *fname)
{
   Source src;
   char *fn;
   LabId id;
   Token tok;
   InputXForm *xf;
   Ptr structure;
   MLink m;
   char buf[MAXSTRLEN], type;
   int fidx = LOADFIDX; /* indicates that these are loaded xform macros */

   if (hset == NULL) { /* read a transform with no model set */
      fn = InitXFormScanner(hset, macroname, fname, &src, &tok);
      SkipWhiteSpace(&src);
      if(GetToken(&src,&tok)<SUCCESS){
         TermScanner(&src);
         return(NULL);
      }
      /* handle the case where the macro header is included */
      if (tok.sym == MACRO) {
         type = tok.macroType;
         if (type != 'j') {
            HRError(999,"LoadInputXForm: Only Input Transform can be specified with no model set %s",src.name);
            return(NULL);
         }
         if (!ReadString(&src,buf)){
            TermScanner(&src);
            HRError(999,"LoadInputXForm: Input XForm macro name expected in macro file %s",src.name);
            return(NULL);
         }
         if (strcmp(buf,macroname)) {
            HRError(999,"LoadInputXForm: Inconsistent macro names  %s and %s in file %s",macroname,buf,src.name);
            return(NULL);
         }
         if(GetToken(&src,&tok)<SUCCESS){
            TermScanner(&src);
            return(NULL);
         }
      }
      xf = GetInputXForm(hset,&src,&tok);
      xf->xformName = CopyString(&xformStack,macroname);
      TermScanner(&src);
   } else { /* First see whether the macro exists */
      id = GetLabId(macroname,FALSE);
      if (id == NULL) { /* macro doesn't exist go find it */
         fn = InitXFormScanner(hset, macroname, fname, &src, &tok);
         SkipWhiteSpace(&src);
         if(GetToken(&src,&tok)<SUCCESS){
            TermScanner(&src);
            return(NULL);
         }
         if (GetXFormMacros(hset,&src,&tok,fidx)<SUCCESS)
            HError(999,"Error in XForm macro file %s",fn);
         /* All xform macros have been read. The speaker macro may be one of them */
         id = GetLabId(macroname,FALSE);
         if (id==NULL) { /* not stored in macro format */
            id = GetLabId(macroname,TRUE);
            structure = xf = GetInputXForm(hset,&src,&tok);
            if (xf->xformName == NULL) /* may have been stored without the macro header */
               xf->xformName = CopyString(hset->hmem,macroname);
            NewMacro(hset,fidx,'j',id,structure);
         } else {
            m = FindMacroName(hset,'j',id);
            xf = (InputXForm *)m->structure;
         }
         TermScanner(&src);
      } else { /* macro already exists so just return it */
         m = FindMacroName(hset,'j',id);
         xf = (InputXForm *)m->structure;
         xf->nUse++;
      }
   }
   if (trace&T_XFM)
      printf("  Using input xform macro \"%s\" from file %s\n",macroname,xf->fname);
   return xf;
}

/* EXPORT->SaveInputXForm: outputs an individual transform */
void SaveInputXForm(HMMSet *hset, InputXForm *xf, char *fname, Boolean binary)
{
   FILE *f;
   Boolean isPipe;

   binary = binary||saveBinary;
   if ((f=FOpen(fname,HMMDefOFilter,&isPipe)) == NULL){
      HError(7011,"SaveInputXForm: Cannot create output file %s",fname);
   }
   if (trace&T_MAC)
      printf("HModel: saving XForm to %s\n",fname);
   /* 
      store the macro header information without explicitly generating
      the macro 
   */
   fprintf(f,"~a %s",ReWriteString(xf->xformName,NULL,DBL_QUOTE));
   if (!binary) fprintf(f,"\n");
   PutInputXForm(hset,f,NULL,xf,FALSE,binary);
   FClose(f,isPipe);  
}


/* ----------------- Output Probability Calculations ---------------------- */


/* CmpTM: return sign(t2->prob - t1->prob) */
static int CmpTM(const void *t1, const void *t2)
{
   if (((TMProb *)t2)->prob < ((TMProb *)t1)->prob)
      return -1;
   else if (((TMProb *)t2)->prob > ((TMProb *)t1)->prob)
      return +1;
   return 0;
}

/* EXPORT->PrecomputeTMix: set up the tmRec[s].probs arrays */
void PrecomputeTMix(HMMSet *hset, Observation *x, float tmThresh, int topM)
{
   TMixRec *tr;
   int m,s,curM;
   LogFloat p,maxP,minP;
   
   for (s=1; s<=hset->swidth[0]; s++) {
      tr = hset->tmRecs+s;
      if (tr->nMix == 0 || tr->mixes == NULL || tr->probs == NULL)
         HError(7092,"PrecomputeTMix: badly formed TMixRec in stream %d",s);
      maxP = LZERO;
      for (m=1; m<=tr->nMix; m++){
         p = MOutP(x->fv[s],tr->mixes[m]);
         if (p>maxP) maxP = p;
         tr->probs[m].prob = p; tr->probs[m].index = m;
      }
      tr->maxP = maxP;
      qsort(tr->probs+1,tr->nMix,sizeof(TMProb),CmpTM);
      if (topM>0) {     /* use fixed number of mixtures */
         if (topM > tr->nMix) curM = tr->nMix;
         else curM = topM;
         for (m=1; m<=curM; m++){
            p = tr->probs[m].prob - maxP;
            tr->probs[m].prob = (p<MINEARG)?0.0:exp(p);
         }
         tr->topM = curM;
      } else {          /* use variable beam of mixtures */
         minP = maxP-tmThresh;
         for (m=1; m<=tr->nMix; m++){
            if (tr->probs[m].prob<minP) break;
            p = tr->probs[m].prob - maxP;
            tr->probs[m].prob = (p<MINEARG)?0.0:exp(p);
         }
         tr->topM = m-1;
      }
   }
}

/* DOutP: Log prob of x in given mixture - Diagonal Case */
static LogFloat DOutP(Vector x, int vecSize, MixPDF *mp)
{
   int i;
   float sum,xmm;

   sum = mp->gConst;
   for (i=1;i<=vecSize;i++) {
      xmm=x[i] - mp->mean[i];
      sum += xmm*xmm/mp->cov.var[i];
   }
   return -0.5*sum;
}

/* FOutP: Log prob of x in given mixture - Full Covariance Case */
static LogFloat FOutP(Vector x, int vecSize, MixPDF *mp)
{
   float sum;
   int i,j;
   Vector xmm;
   TriMat m = mp->cov.inv;
   
   xmm = CreateVector(&gstack,vecSize);
   for (i=1;i<=vecSize;i++)
      xmm[i] = x[i] - mp->mean[i];
   sum = 0.0;
   for (j=1;j<vecSize;j++)
      for (i=j+1;i<=vecSize;i++)
         sum += xmm[i]*xmm[j]*m[i][j];
   sum *= 2;
   sum += mp->gConst;
   for (i=1;i<=vecSize;i++)
      sum += xmm[i] * xmm[i] * m[i][i];
   FreeVector(&gstack,xmm);
   return -0.5*sum;
}


/* COutP: Log prob of x in given mixture - LLT (Choleski) Cov Case */
static LogFloat COutP(Vector x, int vecSize, MixPDF *mp)
{
   HError(7001,"COutP: Choleski storage not yet implemented");
   return 0.0;
}

/* XOutP: Log prob of x in given mixture - XForm Case */
static LogFloat XOutP(Vector x, int vecSize, MixPDF *mp)
{
   Vector xmm,trans_xmm;
   int i,j;
   int numrows;
   Vector xrow;
   LogFloat sum;

   xmm = CreateVector(&gstack,vecSize);
   trans_xmm = CreateVector(&gstack,vecSize);
   for (j=1;j<=vecSize;j++)
      xmm[j] = x[j] - mp->mean[j];
   numrows=NumRows(mp->cov.xform);
   for (i=1;i<=numrows;i++) {
      trans_xmm[i] = 0.0;
      xrow = mp->cov.xform[i];
      for (j=1;j<=vecSize;j++)
         trans_xmm[i] += xrow[j]*xmm[j];
   }        
   sum = 0.0;
   for (i=1;i<=numrows;i++)
      sum += trans_xmm[i]*trans_xmm[i];
   sum += mp->gConst;
   FreeVector(&gstack,xmm);
   return -0.5*sum;
}

/* IDOutP: Log prob of x in given mixture - Inverse Diagonal Case */
static LogFloat IDOutP(Vector x, int vecSize, MixPDF *mp)
{
   int i;
   float sum,xmm;

   sum = mp->gConst;
   for (i=1;i<=vecSize;i++) {
      xmm=x[i] - mp->mean[i];
      sum += xmm*xmm*mp->cov.var[i];
   }
   return -0.5*sum;
}




/* EXPORT-> MOutP: returns log prob of vector x for given mixture */
LogFloat MOutP(Vector x, MixPDF *mp)
{
   int vSize;
   LogFloat px;
   
   vSize=VectorSize(x);
   switch (mp->ckind) {
   case DIAGC:    px=DOutP(x,vSize,mp); break;
   case INVDIAGC: px=IDOutP(x,vSize,mp); break;
   case FULLC:    px=FOutP(x,vSize,mp); break;
   case LLTC:     px=COutP(x,vSize,mp); break;
   case XFORMC:   px=XOutP(x,vSize,mp); break;
   default:       px = LZERO;
   }
   return px;
}


/* EXPORT-> SOutP: returns log prob of stream s of observation x */
LogFloat SOutP(HMMSet *hset, int s, Observation *x, StreamElem *se)
{
   int m,vSize;
   LogDouble bx,px;
   double sum;
   MixtureElem *me;
   MixPDF *mp;
   TMixRec *tr;
   TMProb *tm;
   ShortVec uv;
   Vector v,tv;
   LogFloat wt;
   int ix;

   switch (hset->hsKind){
   case PLAINHS:
   case SHAREDHS:
      v = x->fv[s];
      vSize = VectorSize(v);
      if (vSize != hset->swidth[s])
         HError(7071,"SOutP: incompatible stream widths %d vs %d",
                vSize,hset->swidth[s]);
      me = se->spdf.cpdf+1;
      if (se->nMix == 1){     /* Single Mixture Case */
         mp = me->mpdf; 
         switch (mp->ckind) {
         case DIAGC:    px=DOutP(v,vSize,mp); break;
         case INVDIAGC: px=IDOutP(v,vSize,mp); break;
         case FULLC:    px=FOutP(v,vSize,mp); break;
         case LLTC:     px=COutP(v,vSize,mp); break;
         case XFORMC:   px=XOutP(v,vSize,mp); break;
         default:       px=LZERO;
         }
         return px;
      } else {
         bx = LZERO;                   /* Multi Mixture Case */
         for (m=1; m<=se->nMix; m++,me++) {
            wt=me->weight; wt=MixLogWeight(hset,wt);
            if (wt>LMINMIX) {  
               mp = me->mpdf; 
               switch (mp->ckind) {
               case DIAGC:    px=DOutP(v,vSize,mp); break;
               case INVDIAGC: px=IDOutP(v,vSize,mp); break;
               case FULLC:    px=FOutP(v,vSize,mp); break;
               case LLTC:     px=COutP(v,vSize,mp); break;
               case XFORMC:   px=XOutP(v,vSize,mp); break;
               default:       px = LZERO;
               }
               bx = LAdd(bx,wt+px);
            }
         }
      }
      return bx;
   case TIEDHS:
      v = x->fv[s];
      vSize = VectorSize(v);
      if (vSize != hset->swidth[s])
         HError(7071,"SOutP: incompatible stream widths %d vs %d",
                vSize,hset->swidth[s]);
      sum = 0.0; tr = hset->tmRecs+s;
      tm = tr->probs+1; tv = se->spdf.tpdf;
      for (m=1; m<=tr->topM; m++,tm++)
         sum += tm->prob * tv[tm->index];
      return (sum>=MINLARG)?log(sum)+tr->maxP:LZERO;
   case DISCRETEHS:
      uv = se->spdf.dpdf; m = x->vq[s];
      ix = uv[m];
      if (discreteLZero && ix == DLOGZERO) 
         return LZERO;
      return (float)ix/DLOGSCALE;
   default: HError(7071,"SOutP: bad hsKind %d\n",hset->hsKind);
   }
   return LZERO; /* to keep compiler happy */
}


/* EXPORT-> POutP: returns log prob of streams x for given StateInfo */
LogFloat POutP(HMMSet *hset,Observation *x, StateInfo *si)
{
   LogFloat bx;
   StreamElem *se;
   Vector w;
   int s,S = x->swidth[0];
   
   if (S==1 && si->weights==NULL)
      return SOutP(hset,1,x,si->pdf+1);
   bx=0.0; se=si->pdf+1; w = si->weights;
   for (s=1;s<=S;s++,se++)
      bx += w[s]*SOutP(hset,s,x,se);
   return bx;
}


/* EXPORT-> OutP: Returns probability (log) of observation x for given state */
LogFloat OutP(Observation *x, HLink hmm, int state)
{
   StateInfo *si;
   LogFloat bx;
   StreamElem *se;
   Vector w;
   int s,S = x->swidth[0];
   
   si = (hmm->svec+state)->info;
   if (S==1 && si->weights==NULL)
      return SOutP(hmm->owner,1,x,si->pdf+1);
   bx=0.0; se=si->pdf+1; w = si->weights;
   for (s=1;s<=S;s++,se++)
      bx += w[s]*SOutP(hmm->owner,s,x,se);
   return bx;
}

         
/* EXPORT->DProb2Short: convert prob to scaled log form */
short DProb2Short(float p)
{
   if (p <= MINDLOGP) return 32767;
   return (short) (log(p)*DLOGSCALE);
}

/* EXPORT->Short2DProb: convert scaled log form to prob */
LogFloat Short2DProb(short s)
{
   if (s==32767) return(LZERO);
   else return(exp(((float)s)/DLOGSCALE));
}

/* ----------------------- Fix GConsts ----------------------------- */

/* EXPORT->FixDiagGConst: Sets gConst for given MixPDF in DIAGC case */
void FixDiagGConst(MixPDF *mp)
{
   float sum;
   int i,n;
   LogFloat z;
   Vector v;
   
   v = mp->cov.var; n=VectorSize(v); sum = n*log(TPI);
   for (i=1; i<=n; i++){
      z = (v[i]<=MINLARG)?LZERO:log(v[i]);
      sum += z;
   }
   mp->gConst = sum;
}

/* EXPORT->FixInvDiagGConst: Sets gConst for given MixPDF in INVDIAGC case */
void FixInvDiagGConst(MixPDF *mp)
{
   float sum;
   int i,n;
   LogFloat z;
   Vector v;

   v = mp->cov.var; n=VectorSize(v); sum = n*log(TPI);
   for (i=1; i<=n; i++){
      z = (v[i]<=0.0)?-LZERO:-log(v[i]);
      sum += z;
   }
   mp->gConst = sum;
}

/* EXPORT->FixLLTGConst: Sets gConst for given MixPDF in INVDIAGC case */
void FixLLTGConst(MixPDF *mp)
{
   mp->gConst = 0.0;
}

/* EXPORT->FixFullGConst: Sets gConst for given MixPDF in FULLC case */
/*    note that ldet is passed as a parameter rather than computed  */
/*    here since it will often be known as a byproduct of inverting */
/*    the covariance matrix   */
void FixFullGConst(MixPDF *mp, LogFloat ldet)
{
   mp->gConst = TriMatSize(mp->cov.inv)*log(TPI) + ldet;
}

/* EXPORT->FixGConsts: Sets all gConsts in hmm (must be PLAINHS/SHAREDHS) */
void FixGConsts(HLink hmm)
{
   int n,m,s,S;
   StateElem *ste;
   MixtureElem *me;
   StreamElem *se;
   MixPDF *mp;
   
   ste = hmm->svec+2; S = hmm->owner->swidth[0];
   for (n=2; n<hmm->numStates; n++,ste++) {
      se = ste->info->pdf+1;
      for (s=1; s<=S; s++,se++){
         me = se->spdf.cpdf+1;
         for (m=1; m<=se->nMix; m++,me++)
            if (me->weight > MixFloor(hmm->owner)){
               mp = me->mpdf;
               switch (mp->ckind) {
               case DIAGC:    FixDiagGConst(mp); break;
               case INVDIAGC: FixInvDiagGConst(mp); break;
               case LLTC:     FixLLTGConst(mp); break;
               case FULLC:    FixFullGConst(mp,-CovDet(mp->cov.inv)); break;
               case XFORMC:   break;
               }
            }
      }
   }
}


/* FixTiedGConsts: Sets all gConsts in tied mix sets */
void FixTiedGConsts(HMMSet *hset)
{
   int m,s;
   MixPDF *mp;
   TMixRec *tr;
   
   for (s=1; s<=hset->swidth[0]; s++){
      tr = hset->tmRecs+s;
      for (m=1; m<=tr->nMix; m++){
         mp = tr->mixes[m];
         switch (mp->ckind) {
         case DIAGC:    FixDiagGConst(mp); break;
         case INVDIAGC: FixInvDiagGConst(mp); break;
         case FULLC:    FixFullGConst(mp,-CovDet(mp->cov.inv)); break;
         case XFORMC:   break;
         }
      }
   }
}

/* EXPORT->FixAllGConsts: fix all GConsts in HMM set */
void FixAllGConsts(HMMSet *hset)
{
   int h;
   HLink hmm;
   MLink q;
   
   if (hset->hsKind == PLAINHS || hset->hsKind == SHAREDHS) {
      for (h=0; h<MACHASHSIZE; h++)
         for (q=hset->mtab[h]; q!=NULL; q=q->next)
            if (q->type == 'h'){
               hmm = (HLink)q->structure;
               FixGConsts(hmm);
            }
   }
   else if (hset->hsKind == TIEDHS)
      FixTiedGConsts(hset);
}

/* ----------------------- Enum Conversions ------------------------------ */

/* EXPORT-> DurKind2Str: Return string representation of enum DurKind */
char *DurKind2Str(DurKind dkind, char *buf)
{
   static char *durmap[] = {"NULLD","POISSOND","GAMMAD","RELD","GEND"};
   return strcpy(buf,durmap[dkind]);
}


/* EXPORT-> CovKind2Str: Return string representation of enum CovKind */
char *CovKind2Str(CovKind ckind, char *buf)
{
   static char *covmap[] = {"DIAGC","INVDIAGC","FULLC","XFORMC"};
   return strcpy(buf,covmap[ckind]);
}

/* ------------------------- End of HModel.c --------------------------- */
