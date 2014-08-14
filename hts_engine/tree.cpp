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
/*   tree.c: decision trees handling functions                       */
/*  ---------------------------------------------------------------  */

/* hts_engine libraries */
#include "misc.hpp"
#include "tree.hpp"

/* ------ Pattern Matching functions ------ */
/* DPMatch: recursive matching */
HTS_Boolean DPMatch (char *str, char *pat, const int pos, const int max)
{
   if (pos>max) return 0;
   if (*str=='\0' && *pat=='\0') return 1;
   
   if (*pat=='*') {
      if (DPMatch(str+1, pat, pos+1, max)==1)
         return 1;
      else
         return DPMatch(str+1, pat+1, pos+1, max);
   }
   if (*str==*pat || *pat=='?') {
      if (DPMatch(str+1, pat+1, pos+1, max+1)==1)
         return 1;
      else {
         if (*(pat+1)=='*')
            return DPMatch(str+1, pat+2, pos+1, max+1);
      }
   }
   
   return 0;
}

/* PMatch: pattern matching function */
HTS_Boolean PMatch (char *str, char *pat)
{
   int i, max=0, nstar=0, nques=0;
   char buf[HTS_MAXBUFLEN];
   
   for (i=0; i<(int)strlen(pat); i++) {
      switch (pat[i]) {
         case '*': nstar++; break;
         case '?': nques++; max++; break;
         default:  max++;
      }
   }

   if (nstar==2 && nques==0 && pat[0]=='*' && pat[i-1]=='*') {
      /* only string matching is required */ 
      strncpy(buf, &pat[1], i-2); buf[i-2] = '\0';
      if (strstr(str, buf)!=NULL)
         return 1;
      else
         return 0;
   }
   else 
      return DPMatch(str, pat, 0, (int)(strlen(str)-max));
}

/* QMatch: check given string match given question */
HTS_Boolean QMatch (char *str, Question *q)
{
   HTS_Boolean flag = 0;
   Pattern *p;
  
   for (p=q->phead; p!=q->ptail; p=p->next) {
      flag = PMatch(str, p->pat);
      if (flag)
         return 1;
   }
   
   return 0;
}


/* ------ decition tree handling functions ------*/
/* SearchTree: tree search */
int SearchTree (char *str, Node *node)
{
   while (node!=NULL) {
      if (QMatch(str, node->quest)) {
         if (node->yes->pdf>0)
            return (node->yes->pdf);
         node = node->yes;
      }
      else {
         if (node->no->pdf>0)
             return(node->no->pdf);
         node = node->no;
      }
   }
   return -1;
}

/* LoadQuestions: Load questions from file */
void LoadQuestions (FILE *fp, Question *q, const HTS_Mtype type)
{
   char buf[HTS_MAXBUFLEN];

   HTS_GetToken(fp, buf);
   q->qName = HTS_Strdup(buf);
   q->phead = q->ptail = (Pattern *) HTS_Calloc(1, sizeof(Pattern));

   HTS_GetToken(fp,buf);
   if (strcmp(buf, "{")==0) {
      while (strcmp(buf,"}")!=0) {
          HTS_GetToken(fp, buf);
          q->ptail->pat = HTS_Strdup(buf);
          q->ptail->next = (Pattern *) HTS_Calloc(1, sizeof(Pattern));
          q->ptail = q->ptail->next;
          HTS_GetToken(fp, buf);
      }
   }
   
   return;
}

/* ParseTreePat: parse pattern specified for each tree */
void ParseTreePat (Tree *tree, char *buf) 
{
   char *s, *l, *r;
   
   /* parse tree pattern */
   s = buf;
   if ((l=strchr(s, '{'))!=NULL) {  /* pattern is specified */
      tree->phead = tree->ptail = (Pattern *) HTS_Calloc(1, sizeof(Pattern));
      
      /* manimupate pattern */
      s = l+1;
      if (*s=='(') ++s;
      
      r = strrchr(s,'}');
      if (l<r && *(r-1)==')')
         --r;
      *r = ',';
      
      /* parse pattern */
      while ((l=strchr(s,','))!=NULL) {
         *l = '\0';
         tree->ptail->pat = HTS_Strdup(s);
         tree->ptail->next = (Pattern *) HTS_Calloc(1, sizeof(Pattern));
         tree->ptail = tree->ptail->next;
            
         s = l+1;
      }
   }
   
   return;
}

/* IsTree: check given buffer is tree or not */
HTS_Boolean IsTree (Tree *tree, char *buf)
{
   char *s, *l, *r;

   s = buf;
   if (((l=strchr(s, '['))==NULL) || ((r=strrchr(s, ']'))==NULL)) {
      return 0;
   }
   else {
      *r = '\0';
      s = l+1;
      tree->state = atoi(s);
      ParseTreePat(tree,buf);
   }
   
   return 1;
}

/* IsNum: check given buffer is number or not */
HTS_Boolean IsNum (char *buf)
{
   int i;

   for (i=0; i<(int)strlen(buf); i++)
      if (!(isdigit(buf[i]) || (buf[i]=='-'))) 
         return 0;
      
   return 1;
}

/* FindQuestion: find question from question list */
Question *FindQuestion(TreeSet *ts, const HTS_Mtype type, char *buf)
{
   Question *q;
   
   for (q=ts->qhead[type]; q!=ts->qtail[type]; q=q->next)
      if (strcmp(buf, q->qName)==0)
         break;
         
   if (q==ts->qtail[type])
      HTS_Error(1, "FindQuestion: cannot find question %s.\n", buf);
   
   return q;
}

/* name2num: convert name of node to node index number */
int name2num(char *buf)
{
   return(atoi(strrchr(buf,'_')+1));
}

/* FindNode: find node for given node index */
Node *FindNode (Node *node, const int num)
{
   while(node){
      if (node->idx==num) return node;
      else node=node->next;
   }
   return NULL;
}
         
/* LoadTree: Load trees */
void LoadTree (TreeSet *ts, FILE *fp, Tree *tree, const HTS_Mtype type)
{
   char buf[HTS_MAXBUFLEN];
   Node *node;
   
   HTS_GetToken(fp, buf);
   node = (Node *) HTS_Calloc(1, sizeof(Node));
   tree->root = tree->leaf = node;
   
   if (strcmp(buf,"{")==0) {
      while (HTS_GetToken(fp,buf),strcmp(buf,"}")!= 0) {
         node = FindNode(tree->leaf, atoi(buf));
         HTS_GetToken(fp, buf);     /* load question at this node */
         
         node->quest = FindQuestion(ts, type, buf);
         node->yes = (Node *) HTS_Calloc(1, sizeof(Node));
         node->no  = (Node *) HTS_Calloc(1, sizeof(Node));

         HTS_GetToken(fp, buf);
         if (IsNum(buf)) {
            node->no->idx = atoi(buf);
         }
         else {
            node->no->pdf = name2num(buf);
         }
         node->no->next = tree->leaf;
         tree->leaf = node->no;
         
         HTS_GetToken(fp, buf);
         if (IsNum(buf)) {
            node->yes->idx = atoi(buf);
         }
         else {
            node->yes->pdf = name2num(buf);
         }
         node->yes->next = tree->leaf;
         tree->leaf = node->yes;
      }
   }
   else {
      node->pdf = name2num(buf);
   }
   
   return;
}
   
/* LoadTreeSet: Load decision tree information */
void LoadTreeSet (TreeSet *ts, const HTS_Mtype type)
{
   char buf[HTS_MAXBUFLEN];
   Question *q;
   Tree *t;
   FILE *fp = ts->fp[type];
   
   q = (Question *) HTS_Calloc(1, sizeof(Question));
   ts->qhead[type] = q;  
   ts->qtail[type] = NULL;

   t = (Tree *) HTS_Calloc(1, sizeof(Tree));
   t->next = NULL; t->root = t->leaf = NULL;
   t->phead = t->ptail = NULL;
   t->state=0;
   
   ts->thead[type] = t;  ts->ttail[type] = NULL;
   ts->nTrees[type] = 0;
   
   /* parse tree files */
   while (!feof(fp)) {
      HTS_GetToken(fp, buf);
      /* parse questions */
      if (strcmp(buf, "QS")==0) {
         LoadQuestions(fp, q, type);
         q->next = (Question *) HTS_Calloc(1, sizeof(Question));
         q = ts->qtail[type] = q->next;
         q->next = NULL;
      }
      /* parse trees */
      if (IsTree(t, buf)) {
         LoadTree(ts, fp, t, type);
         t->next = (Tree *) HTS_Calloc(1, sizeof(Tree));
         t = ts->ttail[type] = t->next;
         t->next = NULL; t->root = t->leaf = NULL;
         t->phead = t->ptail = NULL;
         t->state=0;
   
         /* increment # of trees */
         ts->nTrees[type]++;
      }
   }
   
   /* no tree information in tree file */
   if (ts->thead[type]->next==NULL) {
      if (type==DUR)
         HTS_Error(1, "LoadTreesFile: no trees for duration are loaded.\n");
      else if (type==LF0)
         HTS_Error(1, "LoadTreesFile: no trees for log f0 are loaded.\n");
      else
         HTS_Error(1, "LoadTreesFile: no trees for mel-cepstrum are loaded.\n");
   }
   
   return;
}

/* ClearNode: recursive function to free Node */
void ClearNode (Node *node)
{
   if (node->yes!=NULL) 
      ClearNode(node->yes);
   if (node->no!=NULL)
      ClearNode(node->no);

   HTS_Free(node);
   
   return;
}

/* ClearTree: clear given tree */
void ClearTree (Tree *tree)
{
   ClearNode(tree->root);

   return;
}

/* ClearQuestion: clear loaded question */
void ClearQuestion (Question *q)
{
   Pattern *p,*next=NULL;

   for (p=q->phead; p!=q->ptail; p=next) {
      HTS_Free(p->pat);
      next = p->next;
      HTS_Free(p);
   }
   HTS_Free(next);
   HTS_Free(q->qName);
   
   return;
}

/* ClearTreeSet: clear decision trees */
void ClearTreeSet (TreeSet *ts, const HTS_Mtype type)
{
   Question *q,*qnext=NULL;
   Tree *t,*tnext=NULL;
   
   /* close */
   fclose(ts->fp[type]);

   /* free questions */
   for (q=ts->qhead[type]; q!=ts->qtail[type]; q=qnext) {
      ClearQuestion(q);
      qnext = q->next;
      HTS_Free(q);
   }
   HTS_Free(qnext);
   
   /* free trees */
   for (t=ts->thead[type]; t!=ts->ttail[type]; t=tnext) {
      ClearTree(t);
      tnext = t->next;
      HTS_Free(t);
   }
   HTS_Free(tnext);
   
   return;
}      
   
/* InitTreeSet: Initialise TreeSet */
void InitTreeSet (TreeSet *ts) 
{
   ts->fp[DUR] = NULL;
   ts->fp[LF0] = NULL;
   ts->fp[MCP] = NULL;
   
   return; 
} 

/* -------------------- End of "tree.cc" -------------------- */
