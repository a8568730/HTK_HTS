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
/*         File: HGraf.h:   Minimal Graphics Interface         */
/* ----------------------------------------------------------- */
/* Win32 port: Peter Silsbee                                   */

/*  *** THIS IS A MODIFIED VERSION OF HTK ***                        */
/* ----------------------------------------------------------------- */
/*           The HMM-Based Speech Synthesis System (HTS)             */
/*           developed by HTS Working Group                          */
/*           http://hts.sp.nitech.ac.jp/                             */
/* ----------------------------------------------------------------- */
/*                                                                   */
/*  Copyright (c) 2001-2010  Nagoya Institute of Technology          */
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

/* !HVER!HGraf:   3.4.1 [CUED 12/03/09] */

/*
   This module provides a minimal graphics facility.  It provides a 
   single fixed-size window onto which the graphics and text can 
   be painted.  A restricted event handler allows simple interactive
   graphics and window repair.
   
   The coordinate system used has (0,0) at the top left corner with
   x increasing left and y increasing down.  All draw operations use
   the current colour, line width and transfer mode.
*/

/*  
    WIN32: The window is resizable and can be minimized to the taskbar, 
    but the drawing area has a fixed size. The window is automatically
    repainted when it is uncovered, so most applications can ignore
    HREDRAW events, although they are still generated.
*/

#ifndef _HGRAF_H_
#define _HGRAF_H_

#ifdef __cplusplus
extern "C" {
#endif

#define MAX_GREYS   64      /* implementations may quantise further */
#define MAX_COLOURS 16

enum _HColour { WHITE, YELLOW, ORANGE, RED, MAUVE, PURPLE, DARK_BLUE, 
                LIGHT_BLUE, DARK_GREEN, LIGHT_GREEN, DARK_BROWN, LIGHT_BROWN, 
                LIGHT_GREY, GREY, DARK_GREY, BLACK};
typedef enum _HColour HColour;  /* implementations may map these onto grey */


enum _XferMode {GCOPY, GOR, GXOR, GINVERT};
typedef enum _XferMode XferMode;

/* The following types define the small subset of events which HGraf
   assumes can be generated by the graphics engine.  Where multi-button
   mouses are available - only the left button is used.  It is assumed
   that modifier keys such as the shift key generate both keypress
   and keyrelease events.  The redraw event indicates that part of the
   window is damaged and should be redrawn.  There is no facility
   to specify that only a subpart of the window has been damaged.
   At least one redraw event is guaranteed immediately after calling
   InitGraf
*/

enum _HEvent {HMOUSEDOWN,  /* (left) mouse button pressed */
              HMOUSEUP,    /* (left) mouse button released */
              HMOUSEMOVE,  /* mouse has moved */
              HKEYPRESS,   /* key pressed */
              HKEYRELEASE, /* key released */
              HREDRAW      /* window damaged */
};
enum _KeyType {NORMALKEY, SHIFTKEY, COMMANDKEY, CONTROLKEY, 
               ENTERKEY, DELKEY, ESCKEY};

typedef enum _HEvent HEvent;
typedef enum _KeyType KeyType;

typedef struct {
   HEvent event;     /* type of event */
   int x,y;          /* position of mouse */
   unsigned char c;  /* keypress char */
   KeyType ktype;    /* type of key pressed */
} HEventRec;

typedef struct {     
  short x, y;
} HPoint;

/* ------------------ Initialisation and Termination -----------------------*/

void InitGraf(void);
/*
   Initialise the module
*/

void ResetGraf(void);
/*
   Reset the module
*/

void MakeXGraf(char *wname, int x, int y, int w, int h, int bw);
/* 
   Create a window of width w and height h, with top left corner at
   x,y and border width bw.
*/

void TermHGraf(void);
/*
   Close down graphics window and clean up.  This call is superfluous
   in vanilla ANSI C implementations which implement at_exit.
*/


/* ------------------------ Event Handling --------------------------- */

HEventRec HGetEvent(Boolean anyEvent, void (*action)(void));
/* 
   Return next event in event queue.  If an action routine is
   supplied and there are no events pending, then that routine 
   is called repeatedly until an event occurs.  If anyEvent is 
   FALSE then only events for the HGraf window are returned,
   otherwise, all events are returned.

   WIN32: in practice, anyEvent is effectively ignored. If it is TRUE,
   then all events for the current thread are returned; but since only
   one window is allowed, all events will be associated with the HGraf
   window.
*/

int HEventsPending(void);
/* 
   Return the number of events pending in the event queue

   WIN32: returns 1 if (number` of events >= 1)
*/

Boolean HMousePos(int *x, int *y);
/* 
   Return mouse pos in x, y, returns TRUE if the pointer 
   is on the window 

   WIN32: Mouse position is not updated when the cursor is outside
   the window, UNLESS the button was pressed within the window and is
   still pressed. This is the only time the application should need
   to know the mouse position.
*/

Boolean IsInRect(int x, int y, int x0, int y0, int x1, int y1);
/*
   Return TRUE if the point (x,y) lies within the rectangle
   with top-left at (x0,y0) and bottom right at (x1,y1)
*/


/* ------------------------ Drawing Primitives ----------------------- */

void HDrawLine(int x0, int y0, int x1, int y1);
/* 
   Draw a line between x0, y0 and x1, y1 
*/

void HDrawLines(HPoint *points, int n);
/* 
   Draw the sequence of lines which connect the n points stored
   in the array points.  If the last point is the same as the first
   then a polygon will be drawn, otherwise an open curve will be 
   drawn
*/

void HFillPolygon(HPoint *points, int n);
/* 
   Draw/fill the convex polygon whose vertices are stored in
   the n points stored in the array points
*/

void HDrawRectangle(int x0, int y0, int x1, int y1);
void HFillRectangle(int x0, int y0, int x1, int y1);
/* 
   Draw/fill a rectangle with top left corner at x0,y0
   and bottom right corner at x1,y1
*/

void HDrawArc(int x0, int y0, int x1, int y1, int stAngle, int arcAngle);
void HFillArc(int x0, int y0, int x1, int y1, int stAngle, int arcAngle);
/*
   Draw/fill an arc from stAngle degrees thru arcAngle degrees such
   that if arcAngle=360, the closed contour would just fit inside the
   rect with top-left at (x0,y0) and bottom right at (x1,y1).  Positive
   angles are anticlockwise and 0 degrees is along the positive x
   axis.
*/

void HPrintf(int x, int y, char *format, ...);
/* 
   Print the text defined by format such that the left side of the baseline
   of the first character is at x,y.  The format and following
   argument list work as for the standard C printf.
*/

void HCopyArea(int srcx, int srcy, int width, int height, int destx, int desty);
/* 
   Copy a rectangular area of the given width and height from (srcx,srcy)
   to (destx,desty).   This routine is mainly used for scrolling.
*/

void HPlotVector(int x0,int y0,int x1,int y1, Vector v, int st, 
                 int en, float ymax, float ymin);
/*
   Plot the values from v[st]..v[en] within the rectangle defined 
   by x0,y0,x1,y1 such that peak values of ymax and ymin would just
   touch the top and bottom of the rectangle.  The plot is drawn
   the conventional way round ie positive values in the -ve y direction
*/

/* ----------------------------- Settings --------------------------------- */

void HSetColour(HColour c);
/*
   Set the current colour to c 
*/

void HSetGrey(int g);
/* 
   Set current colour to grey level g in range 0 to 63 
*/

void HSetFontSize(int size);
/* 
   Set the current font size in points, 0 selects the default font.
   Maximum range supported is 8 to 24 but implementations may offer
   less than this.  Actual font size set is closest to that offered
   by implementation.
*/

void HSetLineWidth(int w);
/* 
   Set the line width to w pixels 
*/

void HSetXMode(XferMode m);
/* 
   Set current transfer mode to m
*/

int CentreX(int x, char *str);
/* 
   Return the position at which the the h-center of str will be at x 
*/

int CentreY(int y, char *str);
/* 
   Return the position at which the the v-center of str will be at y 
*/

int HTextWidth(char *str);
/* 
   Return the width of string str in pixels 
*/

int HTextHeight(char *str);
/* 
   Return the height of string str in pixels 
*/


/* --------------------------- Misc Routines -----------------------------*/

void HDrawImage(unsigned char *p, int x, int y, int width, int height);
/*
   Draw a grey scale image width by height pixels with top left
   corner at (x,y).  Each pixel is a single byte stored in
   row order , p points to the first (top-left) pixel

  WIN32: p should be on a word boundary.
*/

void HFlush(void);
/* 
   For buffered graphics systems, this flushes pending draw requests.
   For non-buffered systems it has no effect.
*/

void HSpoolGraf(char *fname);
/* 
   Start saving an image of window in the file fname.  This routine may
   not be implemented.  It is intended for graphics systems which have
   an object oriented dump format such as PICT. 
   
	 WIN32: This call must be balanced by a call to HEndSpoolGraf.
*/

void HEndSpoolGraf(void);
/*
	WIN32: Close the Metafile opened in HSpoolGraf. 
	Should be followed by redraw of the window.
*/

void HDumpGraf(char *fname);
/* 
   Dump a pixel image of current display into fname.  This routine
   may not be implemented.  The format of the dumped image is 
   implementation dependent.

   WIN32: Routine is implemented and creates a BMP file. The size of
   the file depends on the display color resolution in effect when 
   the application is launched; full-color displays store 32 bits for
   each pixel!
*/

/* --------------------- Button facility --------------------- */

typedef struct _HButton *BtnLink;
typedef short ButtonId;

typedef struct _HButton { 
   int x, y, w, h;            /* size of button rectangle */
   HColour fg, bg;            /* colors */
   Boolean lit;               /* if true, invert colors */
   Boolean active;            /* if false, stipple gray */
   Boolean toggle;            /* if true, clicking toggles state */
   char *str;                 /* string in button */
   ButtonId id; 
   BtnLink next;
   void (*action)(void);      /* ptr to function to call whilst button is down */
} HButton;

HButton *CreateHButton(HButton *btnlst, ButtonId btnid, int x, int y, int w, 
                       int h, char *str, HColour fg, HColour bg, void (*action)(void));
/*
   CreateHButton: create a button object at position (x,y) of width/height (w,h),
   foreground/background colours (fg,bg) and str displayed inside, (*acton) will
   be executed when the button is de-pressed.
*/

void RedrawHButton(HButton *btn);
/*
    RedrawHButton: readraw the button on the screen.
*/

void RedrawHButtonList(HButton *btnlst);
/*
   RedrawHButtonList: redraw each item in the list of buttons.
*/

HButton *FindButton(HButton *btnlst, ButtonId key);
/*
   FindButton: return a pointer to the button object with ButtonId.
*/

void SetActive(HButton *btnlst, Boolean active);
/*
   SetActive: set the active field in each button object to active.
*/

HButton *CheckButtonList(HButton *btnlst, int x, int y);
/*
   CheckButtonList: return a pointer to the button object which contains (x,y).
*/

void SetButtonLit(HButton *btn, Boolean lit);
/*
   SetButtonLit: redraw the button object in inverse colours. 
*/

ButtonId TrackButtons(HButton *btnlist, HEventRec hev);
/*
   TrackButtons: track the mouse pointer until the mouse button released and
   return the ButtonId of the object which contains the mouse pointer, 0 if 
   none.
*/

#ifdef __cplusplus
}
#endif

#endif  /* _HGRAF_H_ */

/* --------------------------- End of HGraf.h ---------------------------- */
