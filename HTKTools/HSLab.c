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
/*              2001  Cambridge University                     */
/*                    Engineering Department                   */
/*                                                             */
/*   Use of this software is governed by a License Agreement   */
/*    ** See the file License for the Conditions of Use  **    */
/*    **     This banner notice must not be removed      **    */
/*                                                             */
/* ----------------------------------------------------------- */
/*         File: HSLab.c:   The Speech Label Editor            */
/* ----------------------------------------------------------- */

/*  *** THIS IS A MODIFIED VERSION OF HTK ***                        */
/* ----------------------------------------------------------------- */
/*           The HMM-Based Speech Synthesis System (HTS)             */
/*           developed by HTS Working Group                          */
/*           http://hts.sp.nitech.ac.jp/                             */
/* ----------------------------------------------------------------- */
/*                                                                   */
/*  Copyright (c) 2001-2008  Nagoya Institute of Technology          */
/*                           Department of Computer Science          */
/*                                                                   */
/*                2001-2008  Tokyo Institute of Technology           */
/*                           Interdisciplinary Graduate School of    */
/*                           Science and Engineering                 */
/*                                                                   */
/* All rights reserved.                                              */
/*                                                                   */
/* Redistribution and use in source and binary forms, with or        */
/* without modification, are permitted provided that the following   */
/* conditions are met:                                               */
/*                                                                   */
/* - Redistributions of source code must retain the above copyright  */
/*   notice, this list of conditions and the following disclaimer.   */
/* - Redistributions in binary form must reproduce the above         */
/*   copyright notice, this list of conditions and the following     */
/*   disclaimer in the documentation and/or other materials provided */
/*   with the distribution.                                          */
/* - Neither the name of the HTS working group nor the names of its  */
/*   contributors may be used to endorse or promote products derived */
/*   from this software without specific prior written permission.   */
/*                                                                   */
/* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND            */
/* CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,       */
/* INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF          */
/* MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE          */
/* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS */
/* BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,          */
/* EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED   */
/* TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,     */
/* DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON */
/* ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,   */
/* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY    */
/* OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE           */
/* POSSIBILITY OF SUCH DAMAGE.                                       */
/* ----------------------------------------------------------------- */

char *hslab_version = "!HVER!HSLab:   3.4 [CUED 25/04/06]";
char *hslab_vc_id = "$Id: HSLab.c,v 1.6 2008/06/24 03:19:04 zen Exp $";

/* 
   --------------------------------------------------------------
   
   In order to compile this program, you may need to use different 
   compiler options.  For example, under HP-UX 8.07 you would add 
   
   -D_HPUX_SOURCE
   
   to the compiler flags.  Also the environment variable 
   LPATH will need to be set to "/lib:/usr/lib:/usr/lib/X11R4"

   --------------------------------------------------------------
*/

/*
  Enable use of stat calls to reload files changed by another process - e.g. 
  label files generated by a recogniser

  #define USE_STAT
*/

#define USE_STAT


#include "HShell.h"
#include "HMem.h"
#include "HMath.h"
#include "HWave.h"
#include "HAudio.h"
#include "HLabel.h"
#include "HGraf.h"

#define WINNAME      "HSLab"
#define HSLAB_INFO "HSLab V2.0 by Valtcho Valtchev and Steve Young"

#define SLEN  256       /* string length for various buffers */
#define MAX_ZOOM 100       /* maximum level of zooming */

#define INIT_XPOS 10       /* initial x-position of the HSLab window */
#define INIT_YPOS 10       /* initial y-position of the HSLab window */
#define WIDTH    700       /* window width */
#define HEIGHT   500       /* window height */

#define MAX_LAB_LEN 128   /* maximum length of a string representing a label */
#define LAB_BUF_LEN 256   /* string length for buffers used to store labels */


/* define 16-bit resolution, assume that shorts are 16-bit */
#define MIN_AMPL   SHRT_MIN 
#define MAX_AMPL   SHRT_MAX
#define AMPL_RANGE USHRT_MAX

/* button configuration relative to the main window size */
#define BTN_AREA_W   1.00      /* area width as fraction of main window width */
#define BTN_AREA_H   0.25      /* area height as fraction of main window height */
#define BTN_AREA_X   0.0       /* top left corner x of button area as fraction of main win width */
#define BTN_AREA_Y   0.75      /* top left corner y of button area as fraction of main win height */
#define BTN_PER_ROW 10.0       /* how many buttons in a single row in button area */
#define BTN_PER_COL  3.0       /* how many buttons in a column in button area */
#define BTN_H_SPC    0.20      /* vertical spacing between rows of buttons as fraction of btn_width */
#define BTN_V_SPC    0.33      /* horizontal spacing between cols of buttons as fraction of btn_height */

#define SCROLL_PT       20
#define NUM_OF_SCALES    6     /* the number of different amplitude scales */
#define VOLUME_STEPS    10     /* volume control button steps */
#define TIMES_LABELLED   2     /* maximum number of times a frame can be labelled */

/* -------------------------- Global Types  ------------------------ */

typedef enum {
   WAVE_WIN, LAB_WIN, IO_WIN, NO_WIN 
} WinKind;

typedef enum {
   BT_NONE, BT_LOAD, BT_SAVE, BT_SCRLL, BT_SCRLR, 
   BT_ZOOM_IN, BT_ZOOM_OUT, BT_RESTORE, BT_PLAY, 
   BT_PLAY_VOL, BT_SCALE, BT_MARK, BT_UNMARK, 
   BT_LABSTR, BT_LABEL, BT_LABELAS, BT_LDELETE, 
   BT_LEDIT, BT_LSELECT, BT_ABOUT, BT_ADJUST, 
   BT_LABSET, BT_NEWSET,
   BT_UNDO, BT_SPCL, BT_REC, BT_PAUSE, BT_STOP, BT_QUIT
} BtnId;

typedef struct {           /* used to record different levels of zooming */
   long st, en;            /* start and end sample of the previous zoom level */
   int   regA, regB;       /* marked regions on the screen (if any) */
} ZoomRecord;

typedef enum {
   CREATE_OP, DELETE_OP, CHANGE_OP
} OpType;

typedef struct {        /* used to keep track of tha last operation carried out */
   OpType op;           /* the operation carried out */
   int lset;            /* label set */
   Label lab;           /* the label prior to change */
   LLink ptr;           /* pointer to the label last modified (if any) */
} UndoRecord; 

typedef struct{           /* the structure describing the various windows on the main window */
   int x, y, w, h, bw;    /* x-pos, y-pos, width, height, fg and bg colours */
   HColour fg, bg;
} RectWin;

/* ----------------------- Global Variables ---------------------------- */

#define STACKSIZE 100000     /* assume ~100K wave files */
static MemHeap tmpStack;     /* temporary storage with reset */
static MemHeap tmpCHeap;     /* storage allocated once */
static MemHeap wavStack;     /* storage for waveforms */
static MemHeap labStack;     /* storage for label objects */
static MemHeap audStack;

static short *data;          /* the data samples */
static long nSamples;        /* number of samples */
static HTime sampPeriod;     /* the sample period */
static Wave wave;            /* the input waveform */
static Boolean newData=FALSE;

static FileFormat ifmt=UNDEFF;      /* Label input file format */
static FileFormat ofmt=UNDEFF;      /* Label output file format */
static FileFormat wfmt=UNDEFF;      /* wave file format */

static HButton *btnList;     /* the list of buttons displayed on the HSLab window */
static RectWin waveWin;      /* the waveform window */
static RectWin fileWin;      /* the window for file names and other info */
static RectWin labWin1;      /* the label window */
static RectWin io_Win;       /* the I/O window used for messages and inputing strings */

static int trace = 0;                /* trace level */
static int nParm = 0;                /* total num params */
static ConfParam *cParm[MAXGLOBS];   /* configuration parameters */

static char labstr[LAB_BUF_LEN] = "Speech";
static HButton *labstr_btn;

static char volStr[10];      /* the string displayed in the volume button */
static short playVol = 0;    /* the volume for playing samples */
static HButton *vol_btn;     /* pointer to the volume button */

static HButton *lev_btn;     /* the level button */
static char levStr[10];      /* the level button string */

static int btn_v_spc;        /* calculated button vertical spacing */
static int btn_h_spc;        /* calculated button horizontal spacing */
static int btn_w;            /* calculated button width */
static int btn_h;            /* calculated button height */

static int *plotBuf;         /* buffer to store waveform samples to plot */

static int zoomLev;                    /* the current zoom level - ptr in the zoom record */
static ZoomRecord zoomRec[MAX_ZOOM];   /* the zoom record */

static float samplesPt;      /* samples per point of the graphics screen */
static long  T;              /* number of samples */
 
/* waveform scaling */
static float scValues[NUM_OF_SCALES] = {1.0, 2.0, 4.0, 8.0, 16.0, 32.0};
static char *scString[NUM_OF_SCALES] = {"x1", "x2", "x4", "x8", "x16", "x32"};
static short scIndex = 0;    /* initial scaling index == 1.0*/
static HButton *scale_btn;   /* pointer to the amplitude scale button */

static Boolean wavePtrOn = FALSE;   /* keeps track of the state of the waveform pointer */
static int thisWpos, lastWpos;      /* the positions of the waveform pointer */

static char labfn[SLEN];                   /* the label file name */
static char spfn[SLEN] = "noname.wav";     /* the speech (waveform) file name */
static char *ospfn = "noname.wav";         /* the speech file name on the command line */
static char *labDir = NULL;                /* directory for label files */
static char *labExt = "lab";               /* label file extension */
static long sStart, sEnd;                  /* the start and end sample currently visible on the screen */

static Transcription *trans;    /* the transcriptions */
static LabList *llist;          /* the current label list edited */
static int labSet;              /* the number of the current label list */
static int numSet;              /* number of alternative transcriptions */

static Boolean labsModified = FALSE;   /* tracks any changes made to the labels */
static Boolean newLabFile   = FALSE;   /* forcing the creation of a new (empty) label file */
static Boolean incLab       = FALSE;   /* increment global label */

static Boolean regnMarked = FALSE;     /* region marked */
static int markA, markB;               /* boundaries of a marked region */

static UndoRecord undo;                /* the undo record variable */
static Boolean undoEmpty;              /* shows the status of the record */

static HButton *spcl_btn;                 /* pointer to special button */
static char spcl_str[SLEN] = "Command";   /* special button string */

/* ---------------- Process Command Line ------------------------- */

void ReportUsage(void)
{
   printf("\nUSAGE: HSLab [options] waveformFile\n\n");
   printf("\nModified for HTS\n");
   printf(" Option                                       Default\n\n");
   printf(" -a      auto-increment global label          off\n");
   printf(" -i s    Output transcriptions to MLF s       off\n");  
   printf(" -n      Create new transcription             off\n");
   printf(" -s str  Set special button label             command\n");
   PrintStdOpts("FGILPX");
}

int main(int argc, char *argv[])
{
   char *s;

   void Initialise(void);
   void LoadFiles(void);
   void hRedrawWindow(void);
   void DecodeCommands(void);
   void CloseAudio(void);

   if(InitShell(argc,argv, hslab_version, hslab_vc_id)<SUCCESS)
      HError(1500,"HSLab: InitShell failed");

   InitMem(); InitMath(); 
   InitWave(); InitLabel();
   InitAudio(); InitGraf();

   if (!InfoPrinted() && NumArgs() == 0)
      ReportUsage();
   if (NumArgs() == 0) Exit(0);

   while (NextArg() == SWITCHARG) {
      s = GetSwtArg();
      if (strlen(s)!=1) 
         HError(1519,"HSLab: Bad switch %s; must be single letter",s);
      switch(s[0]){
      case 'a':
         incLab = TRUE;
         break;
      case 'i':
         if (NextArg()!=STRINGARG)
            HError(1519,"HSLab: Output MLF file name expected");
         if(SaveToMasterfile(GetStrArg())<SUCCESS)
            HError(1514,"HCopy: Cannot write to MLF");
         break;
      case 'n':
         newLabFile = TRUE;
         break;
      case 's':
         if (NextArg() != STRINGARG)
            HError(1519,"HSLab: Button label expected");
         strcpy(spcl_str, GetStrArg());
         break;
      case 'F':
         if (NextArg() != STRINGARG)
            HError(1519,"HSLab: HSLab: Data File format expected");
         if((wfmt = Str2Format(GetStrArg())) == ALIEN)
            HError(-1589,"HSLab: Warning ALIEN Data file format set");
         break;
      case 'G':
         if (NextArg() != STRINGARG)
            HError(1519,"HSLab: Label File format expected");
         if((ifmt = Str2Format(GetStrArg())) == ALIEN)
            HError(-1589,"HSLab: Warning ALIEN Label input file format set");
         break;
      case 'I':
         if (NextArg() != STRINGARG)
            HError(1519,"HSLab: Input MLF file name expected");
         LoadMasterFile(GetStrArg());
         break;
      case 'L':
         if (NextArg()!=STRINGARG)
            HError(1519,"HSLab: HSLab: Label file directory expected");
         labDir = GetStrArg(); break;
      case 'P':
         if (NextArg() != STRINGARG)
            HError(1519,"HSLab: Label File format expected");
         if((ofmt = Str2Format(GetStrArg())) == ALIEN)
            HError(-1589,"Warning ALIEN Label output file format set");
         break;
      case 'T':
         trace = GetChkedInt(0,65535,s); break;
      case 'X':
         if (NextArg()!=STRINGARG)
            HError(1519,"HSLab: Label file extension expected");
         labExt = GetStrArg(); break;
      default:
         HError(1519,"HSLab: Unknown switch %s",s);
      }
   }
   if (NextArg()==STRINGARG)
      strcpy(spfn, ospfn=GetStrArg());
   Initialise(); LoadFiles(); hRedrawWindow(); 
   DecodeCommands(); TermHGraf();
   
   ResetGraf();
   ResetAudio();
   ResetLabel();
   ResetWave();
   ResetMath();
   ResetMem();
   ResetShell();
   
   Exit(0);
   return (0);          /* never reached -- make compiler happy */
}


/* ------------------------- Bar handling ------------------------- */

#define BAR_WIDTH   200
#define BAR_HEIGHT  15
#define BAR_BORDER  2

typedef struct {
   int x, y;
   HColour fg;
   float step;
} BarType;

static void InitBar(BarType *bar, int x_ofs, HColour fg, int range, char *str)
{
   bar->x = io_Win.x + x_ofs + HTextWidth(str) + 10;
   bar->y = io_Win.y + (io_Win.h - BAR_HEIGHT)/2;
   bar->fg = fg;
   bar->step = (float) (BAR_WIDTH - 2*BAR_BORDER)/(float) range;
   HSetColour(BLACK);
   HPrintf(io_Win.x+x_ofs, io_Win.y+(io_Win.h/2)+HTextHeight(str)/2, "%s", str);
   HSetGrey(30);
   HFillRectangle(bar->x, bar->y, bar->x+BAR_WIDTH, bar->y+BAR_HEIGHT);
}

static void StepBar(BarType *bar, int pos)
{
   int x0, y0, x1, y1;

   HSetGrey(30);
   HFillRectangle(bar->x, bar->y, bar->x+BAR_WIDTH, bar->y+BAR_HEIGHT);
   x0 = bar->x + BAR_BORDER;  y0 = bar->y + BAR_BORDER;
   x1 = x0 + ((int) bar->step*pos);   y1 = bar->y + BAR_HEIGHT-BAR_BORDER; 
   HSetColour(bar->fg);
   HFillRectangle(x0, y0, x1, y1); 
}

/* ------------------------- Audio Drivers -------------------------- */

#define REC_BUF_SIZE 1920000
typedef signed short SampleType;

BtnId DoPause(void);
short ClipSample(short sample, float sampleScale);

/* Record: record signal from audio device */
Wave Record(long *nSamples, HTime *sampPeriod)
{
   Wave w;
   AudioIn ain;
   float bufDur;
   HButton *btn;
   Boolean done;
   HEventRec hev;
   BarType tm, vu;
   char sbuf[MAXSTRLEN];
   long i, chunk;
   long nWaiting, nToRecord;
   SampleType *buf, smin, smax;
   long nEmpty, nTotal, nToFlush;
   SampleType *recbuf;

   ResetHeap(&tmpStack);
   recbuf=(SampleType *) New(&tmpStack,sizeof(SampleType)*REC_BUF_SIZE);
   /* initialise parameters */
   *sampPeriod = 0;
   ain = OpenAudioInput(&tmpStack, sampPeriod, 0.0, 0.0);
   i = (long) (1.0E+07 / (*sampPeriod));
   chunk = (long) (i / 20.0);
   bufDur = (float) REC_BUF_SIZE / (double) i;
   sprintf(sbuf, "rec buffer (0 - %.2fsec)", bufDur);
   InitBar(&tm, 10, RED, REC_BUF_SIZE, sbuf);
   InitBar(&vu, 400, DARK_GREEN, 65536, "VU:");
   buf = recbuf; done = FALSE;
   nEmpty = REC_BUF_SIZE; nTotal =0;
   /* start audio input */
   StartAudioInput(ain, NULLSIG);
   while (nEmpty > 0){
      nWaiting = SamplesInAudio(ain);
      nToRecord = (nWaiting > chunk) ? nWaiting : chunk;
      if (nEmpty < nToRecord)
         nToRecord = nEmpty;
      GetRawAudio(ain, nToRecord, buf);
      /* get min/max samples and update display*/
      smin = 32767; smax = -smin;         
      for (i=0; i<nToRecord; i++){
         if (buf[i]<smin) smin=buf[i];
         if (buf[i]>smax) smax=buf[i];
      }
      StepBar(&tm, nTotal);
      StepBar(&vu, smax-smin);
      buf += nToRecord;
      nEmpty -= nToRecord;
      nTotal += nToRecord;
      if(HEventsPending() > 0){       /* break loop if mouse  button pressed  */
         hev = HGetEvent(FALSE,NULL);
         if(hev.event == HMOUSEDOWN){
            btn = CheckButtonList(btnList, hev.x, hev.y);
            if(btn != NULL){
               switch(btn->id){
               case BT_REC: 
                  done = TRUE; 
                  break;
               case BT_PAUSE: 
                  if(DoPause() == BT_REC)                 /* wait for next button press  */
                     done = TRUE;
                  else {
                     nWaiting = SamplesInAudio(ain);      /* flush accumulated samples */
                     do {
                        nToFlush = (nWaiting > nEmpty) ? nEmpty : nWaiting;
                        GetRawAudio(ain, nToFlush, buf);
                        nWaiting -= nToFlush;
                     } while (nWaiting > 0);
                  }     
                  break;
               }
            }
         }
      }
      if( nEmpty==0 || done) break;
   }
   /* shutdown audio and tranfer data into wave structure */
   StopAudioInput(ain);
   CloseAudioInput(ain);
   *nSamples = nTotal;
   w = OpenWaveOutput(&tmpStack, sampPeriod, *nSamples);
   PutWaveSample(w, *nSamples, recbuf);
   return w;
}

/* Playback: playback part of DataFile with given volume and scale */
void Playback(Wave w, long sampA, long sampB, int vol, int scale, Boolean *newData)
{
   long i;
   AudioOut aout;
   short *data, *p;
   HTime sampP = sampPeriod;
   static long nSamples;
   static short *playBuf;
   static int c_scale = -1;
   int  volmap[10] = {0,10,20,30,40,50,60,70,85,100};

   ResetHeap(&tmpStack);
   /* Reload buffer if new data or scale changed */
   if ((*newData) || (scale!=c_scale)) {
      data = GetWaveDirect(w, &nSamples);
      ResetHeap(&audStack);
      playBuf = (short *) New(&audStack, nSamples*sizeof(short));
      for (p=data, i=0; i<nSamples; i++, p++)
         playBuf[i] = ClipSample(*p, scale);
      *newData = FALSE; c_scale = scale;
   }
   aout = OpenAudioOutput(&tmpStack, &sampP);
   SetVolume(aout, volmap[vol]);
   /* Set up play buffer */
   if (sampA < 0) sampA = 0;
   if (sampB > nSamples) sampB = nSamples;
   i = sampB - sampA;
   if (i>0)   /* play segment of sound */
      StartAudioOutput(aout, i, playBuf+sampA);
   while(SamplesToPlay(aout)>0);
   CloseAudioOutput(aout);
}

/* ---------------------- Command Processing ------------------------ */

static char *command_map[] = {
   "HSLABCMD", "HSLABRUN"
};

typedef enum {
   HSLabCmd, HSLabRun
} HSlabCmdType;

/* CommandSet: returns true and puts cmd in s if envvar set */
Boolean CommandSet(HSlabCmdType hcmd, char *s)
{
   char *env;
 
   if ((env = getenv(command_map[hcmd])) != NULL){
      strcpy(s,env);
      return TRUE;
   } else
      return FALSE;
}

/* ----------------------- Rectangular Windows ----------------------- */

/* InitRectWin: create and initialise a RectWin structure */
void InitRectWin(RectWin *win, float x, float y, float w, float h, 
                 int bw, HColour fg, HColour bg)
{
   win->w = (int) (WIDTH * w);
   win->h = (int) (HEIGHT *h);
   win->x = (int) (WIDTH * x);
   win->y = (int) (HEIGHT * y);
   win->bw = bw; win->fg = fg; win->bg = bg;
}

/* DrawRectWin: draw a rect window on the screen */
void DrawRectWin(RectWin *win)
{
   HSetColour(win->bg);
   HFillRectangle(win->x, win->y, win->x + win->w, win->y + win->h);
   HSetColour(win->fg);
   HSetLineWidth(win->bw);
   HDrawRectangle(win->x, win->y, win->x + win->w, win->y + win->h);
}

/* IsInRectWin: check if (x, y) is in w */
Boolean IsInRectWin(RectWin *w, int x, int y)
{
   return IsInRect(x, y, w->x, w->y, w->x + w->w, w->y + w->h);
}

/* IsInWin: check if (x, y) is in w */
Boolean IsInWin(RectWin *w, int x, int y)
{
   return IsInRect(x, y, w->x, w->y, w->x + w->w, w->y + w->h);
}

/* GetWinKind: return window kind associated with given event */
WinKind GetWinKind(HEventRec hev)
{
   if (IsInWin(&waveWin,hev.x,hev.y)) return WAVE_WIN;
   if (IsInWin(&labWin1,hev.x,hev.y)) return LAB_WIN;
   if (IsInWin(&io_Win,hev.x,hev.y)) return IO_WIN;
   return NO_WIN;
}  

/* Point2Sample: converts a point from the wave form window into a sample number */
long Point2Sample(RectWin *win, int pt)
{
   int t;
   t=(int) (sStart + (float) (pt - (win->x))*samplesPt);
   if (t<sStart) t=sStart;if (t>sEnd) t=sEnd;
   return t;
}

/* ------------------------- Wave Ptr Handling --------------------------- */

/* PlotWaveWinPtr: plot the waveform window pointer in inverted colour */
void PlotWaveWinPtr(int pos)
{
   HSetXMode(GINVERT);
   HDrawLine(pos, waveWin.y + 1, pos, waveWin.y + waveWin.h);
   HDrawLine(pos, labWin1.y + 1, pos, labWin1.y + labWin1.h);
   HSetXMode(GCOPY);
}

/* WPtrOff: turn the waveform window pointer off */
void WPtrOff(void)
{
   if (wavePtrOn){
      PlotWaveWinPtr(thisWpos);
      wavePtrOn = FALSE;
   }
}

/* WPtrOn: turn the waveform window pointer on */
void WPtrOn(void)
{
   if (!wavePtrOn){
      PlotWaveWinPtr(thisWpos);
      wavePtrOn = TRUE;
   }
}

/* ----------------------- Misc Plot Routines ---------------------------- */

/* PlotGStripes: create and draw a grey scale image for about box */
void PlotGStripes(int x, int y, int width, int height)
{
   int i, j;
   unsigned char *pix;
   static unsigned char *gimage = NULL;

   if (gimage==NULL){
      gimage = (unsigned char *) New(&tmpCHeap,width*height*sizeof(unsigned char));
      pix = gimage;
      for (j = 0; j < height; j++)
         for (i = 0; i < width; i++)
            *pix++ = (unsigned char) (i % MAX_GREYS);
   }
   HDrawImage(gimage, x, y, width, height);
}

/* PrintMsg: print a message in the io_window */
void PrintMsg(RectWin *win, char *msg)
{
   int sx, sy, pos, pad = 4;
   char sbuf[MAXSTRLEN];

   HSetXMode(GCOPY);
   DrawRectWin(win);
   if (msg==NULL)
      return;
   strcpy(sbuf, msg); pos = strlen(sbuf);
   while(HTextWidth(sbuf) > (win->w - 2*pad)) 
      sbuf[--pos]='\0';
   HSetColour(win->fg);
   sx = win->x + pad; sy = CentreY(win->y + win->h/2, sbuf);
   HPrintf(sx, sy, "%s", sbuf);
}

/* ----------------------- LabList Processing ------------------------------ */

/* CreateLabObj: create and insert a Label into a sorted LabList */
Label *CreateLabObj(LabList *ll, LabId labid, long st, long en)
{
   LLink p, q;
   
   q = CreateLabel(&labStack, 10);
   q->start = st;    q->end = en;
   q->labid = labid; q->score = 0.0;
   q->succ = NULL;
   q->pred = NULL;
   
   for (p=ll->head->succ; p->succ!=NULL; p=p->succ)
      if (p->start > st)
         break;
   q->succ = p;
   q->pred = p->pred;
   q->pred->succ = q;
   p->pred = q;

   return q;
}

/* TimesLabelled: returns the number of times frame pos is labelled */
int TimesLabelled(LabList *ll, long pos)
{
   int c;
   LLink p;

   c = 0;
   for (p=ll->head->succ; p->succ!=NULL; p=p->succ)
      if (pos >= p->start && pos <= p->end)
         c++;
   return c;
}

/* MaxTimesLabelled: find frame with max number of labels and returns that number */
int MaxTimesLabelled(LabList *ll, long st, long en)
{
   long t;
   int a, max;

   for (max=0, t=st; t<en; t++){
      a = TimesLabelled(ll, t);
      if (a > max) max = a;
   }
   return max;
}

/* DeleteLabObj: deletes q from LabList ll, q must be in lablst */ 
void DeleteLabObj(LabList *ll, LLink q)
{
   LLink p;

   if (ll->head==q || ll->tail==q)
      HError(1590, "DeleteLabObj: attempt to delete sentinel");
   for (p=ll->head; p->succ!=NULL; p=p->succ)
      if (p==q)
         break;
   if (p!=NULL) {
      p->pred->succ = p->succ;
      p->succ->pred = p->pred;
   }
}

/* GetLabDistance: return a pointer to a label whose start or end point are within
   minimum distance of t. The code takes care of exact boundaries 
   i.e. the end of the current label matches the start of it's
   successor by taking (start + eps) and (end - eps) where eps is a 
   small number. */
LLink GetLabDistance(LabList *ll, long t, long st, long en, Boolean *isStart)
{
   long pos;
   LLink p, bp;
   float eps = 0.005;
   float d, min_d, fpos;

   bp=NULL; min_d=INT_MAX;
   for (p=ll->head->succ; p->succ!=NULL; p=p->succ){
      fpos = p->start; pos = (long) fpos; 
      if ((pos >= st) && (pos < en)){
         fpos += eps;
         d = fabs((float) t - fpos); 
         if (d < min_d) { 
            bp = p; min_d = d; *isStart = TRUE;
         }
      }
      fpos = p->end; pos = (long) fpos; 
      if ((pos >= st) && (pos < en)){
         fpos -= eps;
         d = fabs((float) t - fpos); 
         if (d < min_d) { 
            bp = p; min_d = d; *isStart = FALSE;
         }
      }
   }
   return bp;
}

/* GetLabT: returns a pointer to the first Label in the ll which contains t */
LLink GetLabT(LabList *ll, long t)
{
   LLink p;
   
   for (p=ll->head->succ; p->succ!=NULL; p=p->succ)
      if ((t > p->start) && (t < p->end))
         return p;
   return NULL;
}

/* Intersect: returns TRUE if the regions (a,b):{a<=b} and (a1,b1):{a1<=b1} intersect */
Boolean Intersect(long a, long b, long a1, long b1)
{
   return ( (!(((a1 < a) && (b1 < a)) || ((a1 >= b) && (b1 >= b))) ) ? TRUE:FALSE);
}


/* PlotLabels: draws the labels on the screen. This routine assumes that no frame is labelled 
               more than twice. The display is split into two levels for improved clarity.   
               If a label text is too large to fit in the corrsponding box on the screen the  
               code attempts to print a single dot, if there is still no room nothing is 
               printed. Also `..` are added to the start or the end of the label text if the
               labelled region continues off the screen.  */
void PlotLabels(RectWin *win, LabList *ll, long sStart, long sEnd)
{
   LLink p;
   float spt;
   long nSamples;
   HTime st, en, p_en;
   char sbuf[SLEN];
   short shift, p_shift;
   Boolean no_st_bar, no_en_bar;
   char *dot = ".", *dot2 = "..", *str;
   int x, y, h, w, yy, b1, b2, bb, h2, nPoints;

   nPoints  = win->w - 1;
   nSamples = sEnd - sStart;
   spt = (float)nSamples/(float)nPoints;
   x = win->x; y = win->y; w = win->w; h = win->h; 
   h2 = h/2; shift = 0; p_en = 0; p_shift = 0;
   HDrawLine(x, y+h2, x+w, y+h2); 
   /* traverse the list and plot visible labels on the screen */
   for (p=ll->head->succ; p->succ!=NULL; p=p->succ){
      st = p->start;  en = p->end;
      if (!Intersect(sStart, sEnd, (long) st, (long) en))
         continue;
      no_st_bar = no_en_bar = FALSE;   
      str = p->labid->name;
      if (st < sStart) { 
         st = sStart; no_st_bar = TRUE; 
      } 
      b1 =(int) (x + (float)(st - sStart)/spt);
      if (en >= sEnd){
         en = sEnd-1; no_en_bar = TRUE; 
      }
      b2 = (int) (x + (float)(en - sStart)/spt);
      sbuf[0]='\0';
      shift = (st < p_en) ? (p_shift + 1) % 2 : 0;
      yy = y + shift*h2;
      if (no_st_bar) 
         strcat(sbuf, dot2);
      else 
         HDrawLine(b1, yy+1, b1, yy+h2);
      strcat(sbuf, str);
      if (no_en_bar) 
         strcat(sbuf, dot2); 
      else 
         HDrawLine(b2, yy+1, b2, yy+h2);
      bb = (b2 - b1);
      if (HTextWidth(sbuf) < bb)
         HPrintf(CentreX(b1 + bb/2, sbuf), CentreY(yy + h2/2, sbuf), "%s", sbuf);
      else if (HTextWidth(dot) < bb)
         HPrintf(CentreX(b1 + bb/2, dot), CentreY(yy + h2/2, dot), "%", dot);
      if (en >= p_en){ 
         p_en = en; p_shift = shift; 
      }
   }
}

/* PlotLabWin: plot the label window */ 
void PlotLabWin(void)
{ 
   WPtrOff();
   HSetXMode(GCOPY);
   DrawRectWin(&labWin1);
   HSetColour(labWin1.fg);
   HSetLineWidth(0);
   PlotLabels(&labWin1, llist, sStart, sEnd);
}

/* ------------------------- WaveForm Processing --------------------------- */

/* ClipSample: scale a sample, apply clipping if necessary */
short ClipSample(short sample, float sampleScale)
{
   int intSample;

   intSample = (int) (sample * sampleScale);
   if (intSample > MAX_AMPL) intSample = MAX_AMPL;
   if (intSample < MIN_AMPL) intSample = MIN_AMPL;
   return (short) intSample;
}

/* AdjustBounds: adjust the sStart & sEnd such that there are at least 1.0 samples per pixel */
void AdjustBounds(RectWin *win)
{
   long   np, ns, st,   en;
   float p;

   np   = win->w - 1;
   ns = sEnd - sStart;
   p = (float)ns/(float)np;
   if (p > 1.0)
      return;
   en = sStart + np;
   if (en <= T) {
      sEnd = en;
      return;
   }
   st = sEnd - np;
   if (st >= 0) {
      sStart = st;
      return;
   }
}

/* CreatePlotBuf: create the waveform-plot buffer */
void CreatePlotBuf(int size)
{
}

/* PreparePlot: fill the plot buffer with st -> en samples from the data */
void PreparePlot(RectWin *win, short *data, int st, int en)
{
   int i, cent, nPoints;
   float s_st, s_en;
   long nSamples, t;
   short a, b, min, max, sample;
   float WaveScale;
   Boolean rising=FALSE;

   nPoints   = win->w - 1; cent = win->y + win->h/2; nSamples = en - st;
   WaveScale = (float) win->h / (float) 65535.0;
   samplesPt = (float) nSamples / (float)nPoints;
   s_st = (float) st;
   for (i=0; i<nPoints; i++){
      s_en = s_st + samplesPt;   min = MAX_AMPL;   max = MIN_AMPL;
      for (t=(int)s_st; t<(int)s_en; t++){
         sample = data[t];
         sample = ClipSample(sample, scValues[scIndex]);
         if (sample < min) { min = sample; rising = FALSE; }
         if (sample > max) { max = sample; rising = TRUE; }
      }
      if (rising){ a = min; b = max; } else { a = max; b = min; }
      plotBuf[2*i] = (int) (cent - a*WaveScale);
      plotBuf[2*i+1] = (int) (cent - b*WaveScale); s_st = s_en;
   }
}

/* PlotWaveForm: plot the waveform in the waveWin */
void PlotWaveForm(RectWin *win, int st_pt, int en_pt)
{
   int x, i, s;

   HSetColour(win->fg);
   HSetLineWidth(0);
   x = win->x;   s = (st_pt-1)*2;
   for (i=st_pt; i <= en_pt; i++){
      if (i > 1)
         HDrawLine(x+i-1, plotBuf[s-1], x+i, plotBuf[s]);
      HDrawLine(x+i, plotBuf[s], x+i, plotBuf[s+1]);
      s += 2;
   }
}

/* PlotWaveWin: plot the waveform window */ 
void PlotWaveWin(void)
{
   WPtrOff();
   HSetXMode(GCOPY);
   DrawRectWin(&waveWin);
   HSetColour(waveWin.fg);
   HSetLineWidth(0);
   PlotWaveForm(&waveWin, 1, waveWin.w - 1);
   regnMarked = FALSE;
}

/* -------------------------- Miscellaneous Processing --------------------- */

/* TrackWPtr: update the waveform window pointer position */
void TrackWPtr(void)
{
   int x, y;

   if (!HMousePos(&x, &y) || !IsInWin(&waveWin, x, y)){
      WPtrOff();
      return;
   }
   lastWpos = thisWpos; thisWpos = x;
   if (thisWpos==lastWpos)
      return;
   if (wavePtrOn){
      PlotWaveWinPtr(lastWpos);
   }
   PlotWaveWinPtr(thisWpos);
   wavePtrOn = TRUE;
}

/* InvertRegion: invert the window region a - b, taking care of the WavePtr */
void InvertRegion(RectWin *win, int a, int b)
{
   Boolean wasOn=FALSE;

   if (wavePtrOn) { wasOn = TRUE; WPtrOff(); }
   HSetXMode(GINVERT);
   HFillRectangle(a, win->y+1, b, win->y + win->h - 1);
   HSetXMode(GCOPY);
   regnMarked = (!regnMarked) ? TRUE:FALSE;
   if (wasOn) WPtrOn();
}

/* PlotFileWin: plot the file name window */ 
void PlotFileWin(void)
{
   char sbuf[SLEN], nullchar = '\0';
   int sx, sy, pos, pad = 4;
   double rate;

   DrawRectWin(&fileWin);
   HSetColour(fileWin.fg);
   rate = 1.0E+04/(double)sampPeriod;
   sprintf(sbuf, "Waveform: %s, Label: %s, Num samples %ld, HTK sampling rate: %.3fKHz", 
           spfn, labfn, T, rate);
   pos = strlen(sbuf); while(HTextWidth(sbuf) > (fileWin.w - 2*pad)) sbuf[--pos]=nullchar;
   sx = fileWin.x + pad; sy = CentreY(fileWin.y + fileWin.h/2, sbuf);
   HPrintf(sx, sy, "%s", sbuf);
}

/* hRedrawWindow: redraw the main window */
void hRedrawWindow(void)
{
   Boolean wasMarked;

   wasMarked = regnMarked;
   PlotWaveWin();
   PlotLabWin();
   PlotFileWin();
   DrawRectWin(&io_Win);
   RedrawHButtonList(btnList);
   if (wasMarked){
      regnMarked=FALSE;
      InvertRegion(&waveWin, markA, markB);
   }
}

/* GetWavePtrPos: returns the waveform window pointer position after a button press */
Boolean GetWavePtrPos(void)
{
   HEventRec hev;

   do {
      hev = HGetEvent(FALSE, NULL);
      switch (hev.event) {
      case HMOUSEDOWN:
         if (!IsInWin(&waveWin, hev.x, hev.y))
            return FALSE;
         do {
            hev = HGetEvent(TRUE, NULL);
            if (hev.event==HMOUSEMOVE)
               TrackWPtr();
         } while (hev.event!=HMOUSEUP);
         if (IsInWin(&waveWin, hev.x, hev.y))
            return TRUE;
         else
            return FALSE;
         break;
      case HREDRAW:
         hRedrawWindow();
         break;
      case HMOUSEMOVE:
         TrackWPtr(); 
         break;
      default:
         break;
      }
   } while (TRUE);
   return FALSE;
}

/* CalcEnergy: some dude wrote this */
float CalcEnergy(int st,int en)
{
   int i,sample;
   double sq, av;
   
   st=Point2Sample(&waveWin,st);
   en=Point2Sample(&waveWin,en);
   if (en==st) 
      return(0.0);
   else 
      if (st>en) {
         i=st; st=en; en=i;
      }
   sq=av=0.0;
   for (i=st;i<=en;i++) {
      sample = data[i];
      sq+=sample*sample;
      av+=sample;
   }
   av=av/(en-st+1);
   sq=sqrt((sq/(en-st+1))-av*av);
   return(sq);
}

/* GetRegionAB: prompts the user for 2 points on the wave form window */
Boolean GetRegionAB(int *pa, int *pb)
{
   int a;
   char buf[80];

   if (regnMarked) {
      /* if region marked already then use markA and markB */
      *pa = markA;
      *pb = markB;
      InvertRegion(&waveWin, *pa, *pb);
   } else {
      /* otherwisw prompt for a and b */
      PrintMsg(&io_Win, "Select point A...");
      if (!GetWavePtrPos())
         return FALSE;
      *pa = thisWpos;
      PlotWaveWinPtr(*pa);
      PrintMsg(&io_Win, "Select point B...");
      if (!GetWavePtrPos()){
         PlotWaveWinPtr(*pa);
         return FALSE;
      }
      PlotWaveWinPtr(*pa);
      *pb = thisWpos;
      if (*pa > *pb){
         a = *pb; *pb = *pa; *pa = a;
      }
   }
   sprintf(buf,"Region Energy %6.1f (%d..%d)",
           CalcEnergy(markA, markB),
           (int)Point2Sample(&waveWin,markA),
           (int)Point2Sample(&waveWin,markB));
   PrintMsg(&io_Win, buf);
   return TRUE;
}

static int mmWPos;

/* Drag: drag the mouse from MarkA */
void Drag(void)
{
   int x,y;
   
   HMousePos(&x,&y);
   if (mmWPos >=0) PlotWaveWinPtr(mmWPos);
   mmWPos = x;
   PlotWaveWinPtr(mmWPos); 
}

/* MouseMark: mark a region directly using the mouse */
void MouseMark(int x, int *markA, int *markB)
{
   HEventRec hev;
   Boolean done;
   char buf[80];
   
   /* check if region marked, if so unmark */
   if (regnMarked)
      InvertRegion(&waveWin, *markA, *markB);
   *markA = x;
   /* wait for mouse release */
   mmWPos = -1;
   do {
      hev = HGetEvent(TRUE, Drag);
      done = (hev.event==HMOUSEUP) ? TRUE:FALSE;
   } while (!done);
   if (mmWPos >=0) PlotWaveWinPtr(mmWPos);
   *markB = hev.x;
   if (*markA != *markB) {
      InvertRegion(&waveWin, *markA, *markB);
   }
   sprintf(buf,"Region Energy %6.1f (%d..%d)",
           CalcEnergy(*markA, *markB),
           (int)Point2Sample(&waveWin,*markA),
           (int)Point2Sample(&waveWin,*markB));
   PrintMsg(&io_Win, buf);
}

/* GetString:  prompt user to type in string of min length minlen and max length maxlen
               pressing escape returns the original string str, the routine returns 
               TRUE if the length of the typed string is greater than minlen, FALSE 
               otherwise.  The minlen concept is introduced in order to allow for 
               prompts which will not be deleted by the back space character.  If 
               the caller supplies a string whose length is less than minlen 
               backspacing will be allowed until the length of the current string 
               exceeds minlen, from then on backspace will not delete any characters 
               which will make the length < minlen.  */
Boolean GetString(RectWin *win, char *str, short minlen, short maxlen)
{
   HEventRec hev;
   int x, y, sx, sy, cptr, ofs;
   char cursor    = '_';
   char nullchar = '\0';
   char buf[SLEN];
   int pad = 4;
   Boolean done, hitesc;

#define STR_FG   { HSetColour(win->fg); HPrintf(sx, sy, "%s", buf+ofs); }
#define STR_BG   { HSetColour(win->bg); HPrintf(sx, sy, "%s", buf+ofs); }
#define ADJ_OFS {ofs = 0; while (HTextWidth(buf+ofs) > (win->w - 2*pad)) ofs++;}

   DrawRectWin(win);
   strcpy(buf, str);
   cptr = strlen(buf); buf[cptr]=cursor;   buf[cptr+1]=nullchar;
   sx = win->x + pad; sy = CentreY(win->y + win->h/2, buf);
   ADJ_OFS   STR_FG   done = hitesc = FALSE;
   do {
      hev = HGetEvent(FALSE, NULL);
      switch (hev.event) {
      case HKEYPRESS:
         if (!HMousePos(&x, &y))
            break;
         /* if (!IsInWin(win, x, y))
            break; */
         switch (hev.ktype) {
         case ENTERKEY:
            buf[cptr]=nullchar; strcpy(str, buf);
            done = TRUE;
            break;
         case ESCKEY:
            done = hitesc = TRUE;
            break;
         case NORMALKEY:
            if (cptr==maxlen)
               break;
            STR_BG
               buf[cptr++]=hev.c; buf[cptr]=cursor; 
            buf[cptr+1]=nullchar; 
            ADJ_OFS   STR_FG
               break;
         case DELKEY:
            if (cptr==minlen || cptr==0)
               break;
            STR_BG
               buf[--cptr]=cursor; buf[cptr+1]=nullchar; 
            ADJ_OFS   STR_FG
               break;
         default:
            break;
         }
         break;
      case HREDRAW :
         hRedrawWindow();
         STR_FG
            break;
      case HMOUSEMOVE:
         TrackWPtr(); 
         break;
      default:
         break;
      }
   } while (!done);
   DrawRectWin(win);
   return ( ((!hitesc) && (strlen(str) > minlen) ) ? TRUE:FALSE); 
}

/* FileExists: check to see if a file exsists */
Boolean FileExists(char *fn, char *fmode)
{
   Boolean isEXF;               /* File name is extended */
   char actfile[MAXFNAMELEN];   /* actual file name */
   long stindex,enindex;        /* segment indices */
   FILE *f;
   
   strcpy (actfile, fn);
   isEXF = GetFileNameExt (fn, actfile, &stindex, &enindex);

   if ((f=fopen(actfile, fmode)) == NULL)
      return FALSE;
   fclose(f);
   return TRUE;
}

/* LoadLabs: load transcription */
void LoadLabs(void)
{
   LabList *ll;
   LLink p;

   ResetHeap(&labStack);
   MakeFN(spfn, labDir, labExt, labfn);
   /* read labels if label file exists */
   if ((!newLabFile) && (NumMLFFiles()>0 || FileExists(labfn, "r"))){
      trans = LOpen(&labStack,labfn,ifmt);
      llist = trans->head; numSet = trans->numLists; 
      for (ll=trans->head; ll!=NULL; ll=ll->next)
         for (p=ll->head->succ; p->succ!=NULL; p=p->succ) {
            p->start /= sampPeriod; p->end /= sampPeriod;
         }  
   } else {
      trans = CreateTranscription(&labStack);
      llist = CreateLabelList(&labStack, 10);
      AddLabelList(llist, trans); numSet = 1;
   }
   labsModified = FALSE;
   zoomLev = 0; undoEmpty = TRUE; labSet = 0;
}

/* LoadData: read data file from disk */ 
void LoadData(void)
{
   int i;
   short sample0 = 0;

   ResetHeap(&wavStack);
   /* load existing file or create null waveform */
   if (!FileExists(spfn, "r")) {
      nSamples = 10000; sampPeriod = 625;
      wave = OpenWaveOutput(&wavStack, &sampPeriod, nSamples);
      for (i=0; i<nSamples; i++)
         PutWaveSample(wave,1,&sample0);
   } else {
      sampPeriod = 0.0;
      if((wave = OpenWaveInput(&wavStack, spfn, wfmt, 0, 0, &sampPeriod))==NULL)
         HError(1514,"LoadData: OpenWaveInput failed");
   }
   data = GetWaveDirect(wave, &nSamples);
   T = nSamples; sStart = 0; sEnd = T;
   PreparePlot(&waveWin, data, sStart, sEnd);
   zoomLev = 0; undoEmpty = TRUE;   
   newData = TRUE;
}

/* LoadFiles: load the speech & label files */
void LoadFiles(void)
{
   LoadData(); LoadLabs();
}   

/* ------------------------- Undo handling ---------------------------- */

/* RecordOp: record details of the current label operation */
void RecordOp(OpType op, LLink p)
{
   switch (op) {
   case CREATE_OP: 
      undo.ptr = p;
      break;
   case DELETE_OP:
      undo.lab = *p;
      break;
   case CHANGE_OP:
      undo.ptr = p;
      undo.lab = *p;
      break;
   default:
      break;
   }
   undo.op = op; undo.lset = labSet;
   undoEmpty = FALSE;
}

/* UndoOp: use the information stored in undo to undo the last operation */
void UndoOp(void)
{
   Label p;
   char sbuf[MAXSTRLEN];

   if (undoEmpty)
      return;
   if (undo.lset!=labSet) {
      sprintf(sbuf, "Undo operation applies to label set [%d]", undo.lset);
      PrintMsg(&io_Win, sbuf);
      return;
   }
   switch (undo.op) {
   case CREATE_OP:
      undo.op = DELETE_OP; undo.lab = *(undo.ptr);
      DeleteLabObj(llist, undo.ptr);
      break;
   case DELETE_OP:
      undo.op = CREATE_OP; 
      undo.ptr = CreateLabObj(llist, undo.lab.labid, (long) undo.lab.start, (long) undo.lab.end);
      break;
   case CHANGE_OP:
      p = *(undo.ptr);
      *(undo.ptr) = undo.lab;
      undo.lab = p;
      break;
   default:
      break;
   }
   labsModified = TRUE;
   PlotLabWin();
}

/* ---------------------------- Button Commands -------------------------- */

#define FIXLABSTR sprintf(levStr, "Set [%d]", labSet);

/* DoZoomOut: zoom out of the current view */
Boolean DoZoomOut(void)
{
   if (zoomLev==0){ /* nothing to zoom out from */
      return FALSE;
   }
   /* decrement zoomLev, restore previous view (not affected by any scrolls) */
   zoomLev--;
   sStart = zoomRec[zoomLev].st;
   sEnd    = zoomRec[zoomLev].en;
   PreparePlot(&waveWin, data, sStart, sEnd);
   markA = zoomRec[zoomLev].regA;
   markB = zoomRec[zoomLev].regB;
   /* plot the wave form */
   PlotWaveWin();
   PlotLabWin();
   /* re-display the marked region */
   regnMarked = FALSE;
   InvertRegion(&waveWin, markA, markB);
   return TRUE;
}

/* DoZoomIn: zoom into a region */
Boolean DoZoomIn(void)
{
   int zoomA, zoomB;
   int s_st, s_en;

   if (zoomLev==MAX_ZOOM || samplesPt <= 1.0){ /* too many zooms so exit */
      PrintMsg(&io_Win, "Maximum zoom level reached.");
      return FALSE;
   }
   if (!GetRegionAB(&zoomA, &zoomB)){
      PrintMsg(&io_Win, "Zoom In aborted.");
      return FALSE;
   }
   /* record the current zoom level start and end samples */ 
   zoomRec[zoomLev].st = sStart;
   zoomRec[zoomLev].en = sEnd;
   zoomRec[zoomLev].regA = zoomA;
   zoomRec[zoomLev].regB = zoomB;
   zoomLev++;
   /* compute the new start and end samples */
   s_st = Point2Sample(&waveWin, zoomA);
   s_en = Point2Sample(&waveWin, zoomB);
   sStart = s_st; sEnd = s_en;   
   AdjustBounds(&waveWin);
   PreparePlot(&waveWin, data, sStart, sEnd); 
   /* plot the wave form */
   PlotWaveWin();
   PlotLabWin();
   return TRUE;
}

/* DoMark: mark a region */
Boolean DoMark(int *markA, int *markB)
{
   /* check if region marked, if so unmark */
   if (regnMarked)
      InvertRegion(&waveWin, *markA, *markB);
   if (!GetRegionAB(markA, markB)){
      PrintMsg(&io_Win, "Mark region aborted.");
      return FALSE;
   }
   InvertRegion(&waveWin, *markA, *markB);
   return TRUE;
}

/* DoUnMark: unmark marked region */
void DoUnMark(int markA, int markB)
{
   if (regnMarked)
      InvertRegion(&waveWin, markA, markB);
}

/* DoRestore: restore to full view, unmark any regions */
void DoRestore(void)
{
   DrawRectWin(&waveWin);
   wavePtrOn=FALSE;
   HSetColour(BLACK);
   sStart = 0; sEnd = T;   zoomLev = 0;
   PreparePlot(&waveWin, data, sStart, sEnd);
   PlotWaveWin();
   PlotLabWin();
   regnMarked = FALSE;
}

/* DoScrollRight: scroll right thru the waveform */
void DoScrollRight(void)
{
   int x, y, w, h;
   int sShift, ScrollPt;
   
   /* check if already at the end sample, if so exit */
   if (sEnd==T){
      PrintMsg(&io_Win, "Last sample reached");
      return;
   }
   
   /* unmark any marked regions */
   if (regnMarked) 
      InvertRegion(&waveWin, markA, markB);

   sShift = (int) (SCROLL_PT*samplesPt);
   /* adjust the new start and end points */
   if (sEnd + sShift > T){
      sShift = T - sEnd;
      ScrollPt = (int) (sShift / samplesPt);
   } else
      ScrollPt = SCROLL_PT;
   WPtrOff(); 
   x = waveWin.x + 1;   y = waveWin.y + 1;   
   w = waveWin.w - 1;   h = waveWin.h - 1;
   HCopyArea(x + ScrollPt, y, w - ScrollPt, h, x, y);
   HSetColour(waveWin.bg);
   HFillRectangle(x + w - ScrollPt, y, x + w, y + h);
   HSetColour(waveWin.fg);
   sStart += sShift; sEnd += sShift;   
   PreparePlot(&waveWin, data, sStart, sEnd);
   wavePtrOn=FALSE;
   HSetColour(BLACK);
   PlotWaveForm(&waveWin, w - ScrollPt, w); 
   PlotLabWin();
}

/* DoScrollLeft: scroll left thru the waveform */
void DoScrollLeft(void)
{
   int x, y, w, h;
   int sShift, ScrollPt;

   /* check if already at the first sample, if so exit */
   if (sStart==0){
      PrintMsg(&io_Win, "First sample reached");
      return;
   }

   /* unmark region */
   if (regnMarked)
      InvertRegion(&waveWin, markA, markB);

   sShift = (int) (SCROLL_PT*samplesPt);
   if (sStart - sShift < 0){
      sShift = sStart;
      ScrollPt = (int) (sShift / samplesPt);
   } else
      ScrollPt = SCROLL_PT;
   WPtrOff(); 
   x = waveWin.x + 1;   y = waveWin.y + 1;   w = waveWin.w - 1;   h = waveWin.h - 1;
   HCopyArea(x, y, w - ScrollPt, h, x + ScrollPt, y);
   HSetColour(waveWin.bg);
   HFillRectangle(x, y, x + ScrollPt, y + h);
   HSetColour(waveWin.fg);
   sStart -= sShift; sEnd -= sShift;
   PreparePlot(&waveWin, data, sStart, sEnd);
   wavePtrOn=FALSE;
   HSetColour(BLACK);
   PlotWaveForm(&waveWin, 1, ScrollPt+1);
   PlotLabWin();
}

/* DoIncScale: increment the waveform scale factor */
void DoIncScale(void)
{
   Boolean was;
   
   was = regnMarked;
   scIndex = (scIndex + 1) % NUM_OF_SCALES;
   scale_btn->str = scString[scIndex];
   RedrawHButton(scale_btn);
   PreparePlot(&waveWin, data, sStart, sEnd);
   PlotWaveWin();
   if (was){
      regnMarked = FALSE;
      InvertRegion(&waveWin, markA, markB);
   }
   PlotLabWin();
}

/* DoIncVolume: increment the playback volume */
void DoIncVolume(void)
{
   playVol = (playVol + 1) % VOLUME_STEPS;
   sprintf(volStr, "Vol %d", playVol);
   vol_btn->str = volStr;
   RedrawHButton(vol_btn);
}

/* DoIncLabSet: select the next alternative transcription */
void DoIncLabSet(void)
{
   if (numSet==1)
      return;
   if (++labSet < numSet) {
      llist = llist->next;
   } else {
      labSet = 0;
      llist = trans->head;
   }
   if (llist==NULL)
      HError(1590,"DoIncLabSet: Label set/count mismatch error");
   FIXLABSTR;
   RedrawHButton(lev_btn);
   PlotLabWin();
}

/* DoNewLabSet: create a new alternative transcription */
void DoNewLabSet(void)
{
   LabList *ll;
   int c;

   for (c=0, ll=trans->head; ll!=NULL; ll=ll->next, c++)
      if (CountLabs(ll)==0)
         break;
   if (ll==NULL) {
      ll = CreateLabelList(&labStack, 10);
      AddLabelList(ll, trans);
      numSet = trans->numLists;
   }
   llist = ll; labSet=c;
   FIXLABSTR;
   RedrawHButton(lev_btn);
   PlotLabWin();
}

/* DoPlay: play a marked region or the waveform visible */
void DoPlay(void)
{
   long playA, playB;

   if (regnMarked) {
      /* if region marked already then use markA and markB */
      playA = (long) (sStart + (markA - waveWin.x)*samplesPt);
      playB = (long) (sStart + (markB - waveWin.x)*samplesPt);
   } else {
      playA = sStart;
      playB = sEnd;
   }
   Playback(wave, playA, playB, playVol, (int )scValues[scIndex], &newData);
}

/* DoSpecial: execute external command if environment variable set */ 
void DoSpecial(void)
{
   char strbuf[MAXSTRLEN];
   char cmdstr[MAXSTRLEN];
   char cmdfn[MAXSTRLEN];
   
   if (CommandSet(HSLabCmd, cmdstr)){
      strcpy(cmdfn, PathOf(spfn, strbuf));
      strcat(cmdfn, BaseOf(spfn, strbuf));
      while (strchr(cmdstr, '$')!=NULL)
         SubstFName(cmdfn, cmdstr);
      PrintMsg(&io_Win, cmdstr);
      system(cmdstr);
   }
}

void CheckForSave(void);

/* DoRecord: record from audio device */
void DoRecord(short *data)
{
   HButton *btn;
   static int fnkey = 0;   /* key appended to the file name */
   char bfn[SLEN], xfn[SLEN];

   /* change button state */
   SetActive(btnList, FALSE);
   btn = FindButton(btnList, BT_REC);
   btn->str = "stop"; btn->active = TRUE;
   btn = FindButton(btnList, BT_PAUSE);
   btn->active = TRUE;
   RedrawHButtonList(btnList);
   /* wipe out any labels */
   CheckForSave();
   /* record from audio device */
   wave = Record(&nSamples, &sampPeriod);
   T = nSamples;
   /* restore original button state */
   SetActive(btnList, TRUE);
   btn = FindButton(btnList, BT_REC);
   btn->str = "rec";
   btn = FindButton(btnList, BT_PAUSE);
   btn->active = FALSE;
   /* redraw buttons, save data file and reload */
   PrintMsg(&io_Win, "Saving data...");
   RedrawHButtonList(btnList);
   BaseOf(ospfn, bfn); ExtnOf(ospfn, xfn);
   sprintf(spfn, "%s_%d.%s", bfn, fnkey, xfn);
   fnkey = (fnkey+1) % 2;
   MakeFN(spfn, labDir, labExt, labfn);
   if(CloseWaveOutput(wave, wfmt, spfn)<SUCCESS)
      HError(1514,"DoRecord: CloseWaveOutput failed");
   LoadData();
   hRedrawWindow();
}


/* DoPause: pause recording */
BtnId DoPause(void)
{
   HEventRec hev;  
   BtnId btn_id;
   HButton *btn_pause;
   
   btn_pause = FindButton(btnList, BT_PAUSE);
   SetButtonLit(btn_pause, TRUE);
   do {
      do 
         hev = HGetEvent(FALSE, NULL);
      while (hev.event!=HMOUSEDOWN);
      btn_id = (BtnId) TrackButtons(btnList, hev);
   } while ((btn_id!=BT_PAUSE) && (btn_id!=BT_REC));
   SetButtonLit(btn_pause, FALSE);
   return btn_id;
}


/* DoChangeLabStr: change the global string used for labelling */
void DoChangeLabStr(void)
{
   char sbuf[LAB_BUF_LEN], *prompt = "New label string: ";

   strcpy(sbuf, prompt); strcat(sbuf, labstr);
   GetString(&io_Win, sbuf, strlen(prompt), MAX_LAB_LEN);
   strcpy(labstr, sbuf + strlen(prompt));
   RedrawHButton(labstr_btn);
}

/* IncLabStr: increment numeric component of the global string */
void IncLabStr(void)
{
   char sbuf[LAB_BUF_LEN], nbuf[LAB_BUF_LEN],*p,*q;
   Boolean isNum = FALSE;
   int num;

   strcpy(sbuf,labstr);
   for (p=sbuf; *p != '\0'; p++)
      if (isdigit((int) *p)){
         isNum = TRUE; break;
      }
   if (isNum){    
      for (q=p+1; *q != '\0'; q++)
         if (!isdigit(*q)) break;
      num = atoi(p) +1; sprintf(nbuf,"%d",num);
      *p = '\0'; 
      strcpy(labstr,sbuf); strcat(labstr,nbuf); strcat(labstr,q);
      RedrawHButton(labstr_btn);
   }
}

/* DoDelLab: delete a label */
Boolean DoDelLab(void)
{
   LLink p; 
   HTime t;
   
   if (llist==NULL)
      return FALSE;
   PrintMsg(&io_Win, "Select label to delete...");
   if (!GetWavePtrPos()){
      PrintMsg(&io_Win, "Delete label aborted.");
      return FALSE;
   }
   PrintMsg(&io_Win, NULL);
   t = Point2Sample(&waveWin, thisWpos);
   if ((p = GetLabT(llist, (long) t)) == NULL)
      return FALSE;
   RecordOp(DELETE_OP, p);
   DeleteLabObj(llist, p);
   PlotLabWin();
   labsModified = TRUE;
   return TRUE;
}

/* DoEditLabel: change the label string */
Boolean DoEditLabel(void)
{
   HTime t;
   LLink p;
   LabId labid;
   char sbuf[LAB_BUF_LEN], *prompt = "New label: ";

   if (llist==NULL)
      return FALSE;
   PrintMsg(&io_Win, "Select label to edit...");
   if (!GetWavePtrPos()){
      PrintMsg(&io_Win, "Edit label aborted.");
      return FALSE;
   }
   t = Point2Sample(&waveWin, thisWpos);
   if ((p = GetLabT(llist, (long) t))==NULL){
      PrintMsg(&io_Win, "Selected frame not labelled");
      return FALSE;
   }
   labid = p->labid;   
   strcpy(sbuf, prompt); strcat(sbuf, labid->name);
   if (!GetString(&io_Win, sbuf, strlen(prompt), MAX_LAB_LEN)){
      PrintMsg(&io_Win, "Edit label aborted.");
      return FALSE;
   }
   RecordOp(CHANGE_OP, p);
   p->labid = GetLabId(sbuf + strlen(prompt), TRUE);
   PlotLabWin();
   PrintMsg(&io_Win, NULL);
   labsModified = TRUE;
   return TRUE;
}

/* DoSelectlabel: use the boundaries of a chosen label to create a marked region */
Boolean DoSelectLabel(void)
{
   LLink p;
   HTime t, st, en;
   
   if (llist==NULL)
      return FALSE;
   PrintMsg(&io_Win, "Choose label to select as marked region...");
   if (!GetWavePtrPos()){
      PrintMsg(&io_Win, "Select label aborted.");
      return FALSE;
   }
   t = Point2Sample(&waveWin, thisWpos);
   if ((p = GetLabT(llist, (long) t))==NULL){
      PrintMsg(&io_Win, "Selected frame not labelled.");
      return FALSE;
   }
   /* unmark any marked regions */
   if (regnMarked)   
      InvertRegion(&waveWin, markA, markB);
   st = p->start;  if (st < sStart) st = sStart;
   en = p->end;    if (en >= sEnd)  en = sEnd - 1;
   markA = (int) (waveWin.x + (st - sStart)/samplesPt);
   markB = (int) (waveWin.x + (en - sStart)/samplesPt);
   InvertRegion(&waveWin, markA, markB);
   PrintMsg(&io_Win, NULL);
   return TRUE;
}

/* DoLabel: label a marked region or the waveform visible on the screen */
Boolean DoLabel(Boolean useLabStr)
{
   HTime st, en;
   char sbuf[LAB_BUF_LEN], *prompt = "Lab name: ";
   LabId labid;
   
   strcpy(sbuf, prompt);
   if (regnMarked) {
      /* if region marked already then use markA and markB */
      st = Point2Sample(&waveWin, markA);
      en = Point2Sample(&waveWin, markB);
   } else {
      st = sStart;
      en = sEnd - 1;
   }
   if (MaxTimesLabelled(llist, (long) st, (long) en)==TIMES_LABELLED){
      PrintMsg(&io_Win, "Area contains frames already labelled twice.");
      return FALSE;
   }
   if (useLabStr)
      labid = GetLabId(labstr, TRUE);
   else {
      if (!GetString(&io_Win, sbuf, strlen(prompt), MAX_LAB_LEN)){
         PrintMsg(&io_Win, "Labelling aborted.");
         return FALSE;
      }
      labid = GetLabId(sbuf + strlen(prompt), TRUE);
   }
   RecordOp(CREATE_OP, CreateLabObj(llist, labid, (long) st, (long) en));
   PlotLabWin();
   labsModified = TRUE;
   if (incLab) IncLabStr();
   return TRUE;
}

/* DoAdjustLabel: adjust the boundary of a selected label */
Boolean DoAdjustLabel(void)
{
   LLink p;
   Label q, qq;
   Boolean isStart;
   HTime t, nt, *bp;
   char sbuf[SLEN];

   if (llist==NULL)
      return FALSE;
   PrintMsg(&io_Win, "Select boundary to adjust...");
   if (!GetWavePtrPos()){
      PrintMsg(&io_Win, "Adjust boundary aborted");
      return FALSE;
   }
   t = Point2Sample(&waveWin, thisWpos);
   if ((p = GetLabDistance(llist, (long) t, sStart, sEnd, &isStart))==NULL){
      PrintMsg(&io_Win, "Selected frame not labelled or boundary not on-screen.");
      return FALSE;
   }
   q = qq = *p;
   strcpy(sbuf, "Select new boundary position for ["); 
   strcat(sbuf, q.labid->name);
   if (isStart){
      bp = &(q.start);
      strcat(sbuf, " : start]"); 
   } else {
      bp = &(q.end);
      strcat(sbuf, " : end]"); 
   }
   /* get the new boundary position */
   PrintMsg(&io_Win, sbuf);
   if (!GetWavePtrPos()){
      PrintMsg(&io_Win, "Adjust boundary aborted");
      return FALSE;
   }
   nt = Point2Sample(&waveWin, thisWpos);
   if (nt >= T) nt = T-1;
   if (nt <  0) nt = 0;
   *bp = nt; 
   if (MaxTimesLabelled(llist, sStart, sEnd) > TIMES_LABELLED){ 
      PrintMsg(&io_Win, "Boundary change will make some frames labelled more than twice.");
      *p = qq;
      return FALSE;
   }
   RecordOp(CHANGE_OP, p);
   *p = q; PlotLabWin();
   PrintMsg(&io_Win, NULL);
   labsModified = TRUE;
   return TRUE;
}

/* DoAbout: display information about author, version, etc. */
void DoAbout(void)
{
   int width = MAX_GREYS;

   PrintMsg(&io_Win, HSLAB_INFO);
   PlotGStripes(io_Win.x + io_Win.w - width, io_Win.y + 1, width, io_Win.h - 1);
}

/* DoSave: save the changes into a HTK-format label file */
void DoSave(void)
{
   char sbuf[SLEN], fnbuf[SLEN], *prompt = "Save label file: ";
   LabList *ll;
   LLink p;

   if (!labsModified){
      PrintMsg(&io_Win, "No changes need to be saved.");
      return;
   }
   for (ll=trans->head; ll!=NULL; ll=ll->next)
      if (CountLabs(ll)!=0)
         break;
   if (ll==NULL) {
      PrintMsg(&io_Win, "No labels found.");
      return;
   }
   strcpy(sbuf, prompt); strcat(sbuf, labfn);
   if (!GetString(&io_Win, sbuf, strlen(prompt), MAX_LAB_LEN)){
      PrintMsg(&io_Win, "Save file aborted.");
      return;
   }
   strcpy(fnbuf, sbuf + strlen(prompt));
   if (!FileExists(fnbuf, "w")){
      PrintMsg(&io_Win, "File I/O error.");
      return;
   }
   strcpy(labfn, fnbuf);
   for (ll=trans->head; ll!=NULL; ll=ll->next)
      for (p=ll->head->succ; p->succ!=NULL; p=p->succ) {
         p->start *= sampPeriod; p->end *= sampPeriod;
      }  
   if(LSave(labfn, trans, ofmt)<SUCCESS)
      HError(1514, "DoSave: Could not save label file %s", labfn);
   PrintMsg(&io_Win, "Label file saved.");
   labsModified = FALSE;
   PlotFileWin();
}


/* CheckForSave: check to see if changes have to be saved */
void CheckForSave(void)
{
   char sbuf[MAXSTRLEN], *prompt = "Save label file (Y/N): ", c;
   Boolean is0;

   if (!labsModified)
      return;
   do { 
      strcpy(sbuf, prompt);
      do 
         is0 = (!GetString(&io_Win, sbuf, strlen(prompt), MAX_LAB_LEN)) ? TRUE:FALSE; 
      while (is0);
      c = sbuf[strlen(prompt)];
   } while ((c!='y') && (c!='Y') && (c!='n') && (c!='N'));
   if ((c=='y') || (c=='Y'))
      DoSave();
}
   
/* DoLoad: load a new waveform/label files */
Boolean DoLoad(void)
{
   char sbuf[SLEN], fnbuf[SLEN], *prompt = "Load waveform file: ";

   CheckForSave();
   strcpy(sbuf, prompt); 
   if (NextArg()==STRINGARG) strcat(sbuf, GetStrArg());
   else strcat(sbuf, "st");
   if (!GetString(&io_Win, sbuf, strlen(prompt), MAX_LAB_LEN)){
      PrintMsg(&io_Win, "Load file aborted.");
      return FALSE;
   }
   strcpy(fnbuf, sbuf + strlen(prompt));
   if (!FileExists(fnbuf, "r")){
      PrintMsg(&io_Win, "File not found.");
      return FALSE;
   }
   strcpy(spfn, fnbuf);
   LoadFiles(); hRedrawWindow();
   return TRUE;
}

/* PlayLabel: play region for clicked on label */
void PlayLabel(int x)
{
   HTime t;
   LLink p;

   t = Point2Sample(&waveWin, x);
   if ((p = GetLabT(llist, (long) t))==NULL) 
      return;
   Playback(wave, (long) p->start, (long) p->end, playVol, (int) scValues[scIndex], &newData);
}

/* ------------------------- Initialisation ------------------------------- */

/* CreateButtons: create and initialise all buttons */
void CreateButtons(void)
{
   int btn_area_w;
   int btn_area_h;
   int btn_area_x;
   int btn_area_y;
   int x, y, w;
   char cmdstr[MAXSTRLEN];
   HButton *btn;

   btn_area_w = (int) (WIDTH *BTN_AREA_W);
   btn_area_h = (int) (HEIGHT*BTN_AREA_H);
   btn_area_x = (int) (0.0 + (float) WIDTH *BTN_AREA_X);
   btn_area_y = (int) (0.0 + (float) HEIGHT*BTN_AREA_Y);
   btn_w = (int) ((float) btn_area_w/(BTN_PER_ROW + (BTN_PER_ROW + 1.0)*BTN_H_SPC));
   btn_h = (int) ((float) btn_area_h/(BTN_PER_COL + (BTN_PER_COL + 1.0)*BTN_V_SPC));
   btn_h_spc = (int) (btn_w*BTN_H_SPC);
   btn_v_spc = (int) (btn_h*BTN_V_SPC);

   x = btn_area_x + btn_h_spc;   y = btn_area_y + btn_v_spc;
   btnList = CreateHButton(NULL,BT_LOAD,x,y,btn_w,btn_h,"Load",BLACK,LIGHT_GREY,NULL);
   x += btn_w + btn_h_spc;
   CreateHButton(btnList,BT_SAVE,x,y,btn_w,btn_h,"Save",BLACK,LIGHT_GREY,NULL);
   x += btn_w + btn_h_spc;
   CreateHButton(btnList,BT_ABOUT,x,y,btn_w,btn_h,"About",BLACK,LIGHT_GREY,NULL);
   x += btn_w + btn_h_spc;
   CreateHButton(btnList,BT_QUIT,x,y,btn_w,btn_h,"Quit",BLACK,LIGHT_GREY,NULL);
   x += 3*(btn_w + btn_h_spc); w = 2*btn_w + btn_h_spc;
   spcl_btn=CreateHButton(btnList,BT_SPCL,x,y,w,btn_h,spcl_str,BLACK,LIGHT_GREY,NULL);
   if (!CommandSet(HSLabCmd, cmdstr))
      spcl_btn->active = FALSE;
   x += w + btn_h_spc;
   vol_btn = CreateHButton(btnList,BT_PLAY_VOL,x,y,btn_w,btn_h,volStr,BLACK,LIGHT_GREY,NULL);
   x += btn_w + btn_h_spc;
   scale_btn = CreateHButton(btnList,BT_SCALE,x,y,btn_w,btn_h,scString[scIndex],BLACK,LIGHT_GREY,NULL);
   x += btn_w + btn_h_spc;

   /* start the second row of buttons */
   x = btn_area_x + btn_h_spc;   y += btn_h + btn_v_spc;
   CreateHButton(btnList,BT_MARK,x,y,btn_w,btn_h,"Mark",BLACK,LIGHT_GREY,NULL);
   x += btn_w + btn_h_spc;
   CreateHButton(btnList,BT_UNMARK,x,y,btn_w,btn_h,"Unmark",BLACK,LIGHT_GREY,NULL);
   x += btn_w + btn_h_spc;
   CreateHButton(btnList,BT_SCRLL,x,y,btn_w,btn_h,"<--",BLACK,LIGHT_GREY,DoScrollLeft);
   x += btn_w + btn_h_spc;
   CreateHButton(btnList,BT_SCRLR,x,y,btn_w,btn_h,"-->",BLACK,LIGHT_GREY,DoScrollRight);
   x += btn_w + btn_h_spc;
   CreateHButton(btnList,BT_ZOOM_IN,x,y,btn_w,btn_h,"Z.In",BLACK,LIGHT_GREY,NULL);
   x += btn_w + btn_h_spc;
   CreateHButton(btnList,BT_ZOOM_OUT,x,y,btn_w,btn_h,"Z.Out",BLACK,LIGHT_GREY,NULL);
   x += btn_w + btn_h_spc;
   CreateHButton(btnList,BT_RESTORE,x,y,btn_w,btn_h,"Restore",BLACK,LIGHT_GREY,NULL);
   x += btn_w + btn_h_spc;
   CreateHButton(btnList,BT_PLAY,x,y,btn_w,btn_h,"play",BLACK,LIGHT_GREY,NULL);
   x += btn_w + btn_h_spc;
   CreateHButton(btnList,BT_REC,x,y,btn_w,btn_h,"rec",BLACK,LIGHT_GREY,NULL);
   x += btn_w + btn_h_spc;
   btn = CreateHButton(btnList,BT_PAUSE,x,y,btn_w,btn_h,"pause",BLACK,LIGHT_GREY,NULL);
   btn->active = FALSE;

   /* start the third row of buttons - labelling functions */
   x = btn_area_x + btn_h_spc; y += btn_h + btn_v_spc; 
   CreateHButton(btnList,BT_LABEL,x,y,btn_w,btn_h,"Label",BLACK,LIGHT_GREY,NULL);
   x += btn_w + btn_h_spc;
   CreateHButton(btnList,BT_LABELAS,x,y,btn_w,btn_h,"Labelas",BLACK,LIGHT_GREY,NULL); 
   x += btn_w + btn_h_spc;
   CreateHButton(btnList,BT_LDELETE,x,y,btn_w,btn_h,"Delete",BLACK,LIGHT_GREY,NULL);
   x += btn_w + btn_h_spc;
   CreateHButton(btnList,BT_LEDIT,x,y,btn_w,btn_h,"Edit",BLACK,LIGHT_GREY,NULL);
   x += btn_w + btn_h_spc;
   CreateHButton(btnList,BT_LSELECT,x,y,btn_w,btn_h,"Select",BLACK,LIGHT_GREY,NULL);
   x += btn_w + btn_h_spc;
   CreateHButton(btnList,BT_ADJUST,x,y,btn_w,btn_h,"Adjust",BLACK,LIGHT_GREY,NULL);
   x += btn_w + btn_h_spc;
   lev_btn = CreateHButton(btnList,BT_LABSET,x,y,btn_w,btn_h,levStr,BLACK,LIGHT_GREY,NULL);
   x += btn_w + btn_h_spc;
   CreateHButton(btnList,BT_NEWSET,x,y,btn_w,btn_h,"New",BLACK,LIGHT_GREY,NULL);
   x += btn_w + btn_h_spc;
   CreateHButton(btnList,BT_UNDO,x,y,btn_w,btn_h,"Undo",BLACK,LIGHT_GREY,NULL);
   x += btn_w + btn_h_spc; w = 3*btn_w + 2*btn_h_spc;
   labstr_btn = CreateHButton(btnList,BT_LABSTR,x,y,btn_w,btn_h,labstr,BLACK,LIGHT_GREY,NULL);
}

/* CreateSubWindows: Create the various sub-windows */
void CreateSubWindows(void)
{
   CreateButtons();
   InitRectWin(&fileWin, 0.01, 0.01, 0.98, 0.05, 1, BLACK, WHITE);
   InitRectWin(&waveWin, 0.01, 0.06, 0.98, 0.50, 1, BLACK, WHITE);
   InitRectWin(&labWin1, 0.01, 0.56, 0.98, 0.10, 1, BLACK, LIGHT_GREY);
   InitRectWin(&io_Win,   0.01, 0.68, 0.98, 0.05, 1, BLACK, LIGHT_GREY);
}

/* SetConfParms: set conf parms relevant to HSLab */
void SetConfParms(void)
{
   int i;
   
   nParm = GetConfig("HSLAB", TRUE, cParm, MAXGLOBS);
   if (nParm>0) {
      if (GetConfInt(cParm,nParm,"TRACE",&i)) trace = i;
   }
}

/* Initialise: initialise everything */
void Initialise(void)
{
   SetConfParms();
   CreateHeap(&tmpStack, "tmpStack", MSTAK, 1, 1.2, 512, 4096);
   CreateHeap(&wavStack, "wavStack", MSTAK, 1, 1.2, 512, 4096);
   CreateHeap(&tmpCHeap, "tmpCHeap", CHEAP, 1, 0.0, 0, 0);
   CreateHeap(&labStack, "labStack", MSTAK, 1, 1.2, 4096, 8192);
   CreateHeap(&audStack, "audStack", MSTAK, 1, 4, REC_BUF_SIZE, LONG_MAX);
   MakeXGraf(WINNAME, INIT_XPOS, INIT_YPOS, WIDTH, HEIGHT, 4);
   /* init volume setting */
   playVol = (VOLUME_STEPS*2 / 3);
   sprintf(volStr, "Vol %d", playVol);
   /* set alternative transcription */
   labSet = 0; FIXLABSTR;
   CreateSubWindows();
   plotBuf = (int *) New(&tmpCHeap, 2*waveWin.w*sizeof(int));
}


/* -------------------- Main command loop -------------------- */

/* DecodeCommands: decode the button commands */
void DecodeCommands(void)
{
   HEventRec hev;
   BtnId btn_id;
   WinKind wkind;

   do {
      hev = HGetEvent(FALSE, NULL);
      switch (hev.event) {
      case HMOUSEDOWN :
         wkind = GetWinKind(hev);
         PrintMsg(&io_Win, NULL);
         switch (wkind) {
         case WAVE_WIN: 
            MouseMark(hev.x, &markA, &markB);
            break;
         case LAB_WIN: 
            PlayLabel(hev.x);
            break;
         case IO_WIN: break;
         default:       
            btn_id = (BtnId) TrackButtons(btnList, hev);
            switch(btn_id) {
            case BT_ZOOM_IN   : DoZoomIn();   break;
            case BT_ZOOM_OUT  : DoZoomOut(); break;
            case BT_RESTORE   : DoRestore();    break;
            case BT_SCRLL     : DoScrollLeft(); break;
            case BT_SCRLR     : DoScrollRight();break;
            case BT_PLAY      : DoPlay(); break;
            case BT_REC       : DoRecord(data); break;
            case BT_PLAY_VOL  : DoIncVolume(); break;
            case BT_SCALE     : DoIncScale(); break;
            case BT_MARK      : DoMark(&markA, &markB); break;
            case BT_UNMARK    : DoUnMark(markA, markB); break;
            case BT_LABSTR    : DoChangeLabStr(); break;
            case BT_LDELETE   : DoDelLab(); break;
            case BT_LEDIT     : DoEditLabel(); break;
            case BT_LSELECT   : DoSelectLabel(); break;
            case BT_LABELAS   : DoLabel(FALSE); break;
            case BT_LABEL     : DoLabel(TRUE); break;
            case BT_LABSET    : DoIncLabSet(); break;
            case BT_NEWSET    : DoNewLabSet(); break;
            case BT_ABOUT     : DoAbout(); break;
            case BT_ADJUST    : DoAdjustLabel(); break;
            case BT_SAVE      : DoSave(); break;
            case BT_LOAD      : DoLoad(); break;
            case BT_QUIT      : CheckForSave(); return;
            case BT_UNDO      : UndoOp(); break;
            case BT_SPCL      : DoSpecial(); break;
            default: break;
            }
            break;
         }
         break;
      case HREDRAW :
         hRedrawWindow();
         break;
      case HMOUSEMOVE:
         TrackWPtr(); 
         break;
      default:
         break;
      }
   } while(TRUE);
}

/* ----------------------------------------------------------- */
/*                      END:  HSLab.c                          */
/* ----------------------------------------------------------- */
