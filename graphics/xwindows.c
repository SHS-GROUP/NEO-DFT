/*  15 Apr 2002

      CALCOMP-to-Xwindows conversion code.

      Written in March 1990 by Mike Schmidt

      Modified in October 1993 by Greg Atchity, Hideaki Shimizu,
      and Nikita Matsunaga to provide better interaction with the
      window manager, and to use pixmaps for instantaneous drawing.

      The following CALCOMP calls are (partially) implemented:
         traditional      herein
         ----------       ------
          PLOTS           OPENCC
            -             FLSHCC
          PLOT            CLOSCC and PLOT
          FACTOR          FACTOR (as a do nothing)
          WHERE           WHERE
          SYMBOL          SYMBOL
          NUMBER          NUMBER
          NEWPEN          NEWPEN
            -             ERASE
            -             ZAPBOX
          SCALE             -
          LINE              -
          AXIS              -
       For more information on each routine, and its usage from FORTRAN,
       see the comments preceeding each.  Note that not all of the
       functionality of each routine has been provided in every case.
*/
#ifdef UNIX
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#endif

#ifdef VMS
#include <decw$include/Xlib.h>
#include <decw$include/Xutil.h>
#endif
 
#include <string.h>
#include <stdio.h>
#include <math.h>

#define TRUE  1
#define FALSE 0
/*
     FONT1 is used by most machines,
     FONT2 is loaded for MACs running MACX at ISU.
     another font you might try is "courr10" or "helvr10"
*/
#define Font1 "-ADOBE-NEW CENTURY SCHOOLBOOK-MEDIUM-R-NORMAL--*-140-*-*-P-*"
#define Font2 "-ADOBE-NEW CENTURY SCHOOLBOOK-MEDIUM-R-NORMAL--*-100-*-*-P-*"

/*
    globally known Xwindows "structures"
*/
Display              *display;
Window               window;
GC                   gc;
Visual               *visual;
Screen               *screen;
Font                 font;
XColor               exactC,screenC;
Colormap             cmap;
XSetWindowAttributes xswa;
XGCValues            xgcv;
XEvent               event;
int                  screen_num;

#ifdef DIRECTDRAW
#define  off_screen  window
#else
Pixmap               off_screen;
#endif

/*
    globally known variables
*/
static char  *colors[16];
static int   height, width, monochrome;
static int   Height, Width;
static int   flip, xsave, ysave, xoff, yoff, ccopen=FALSE;
static float scale,factr,xmx,ymx;

/*
    Convert subroutine names to upper case or underscores
*/

#ifdef CAPS
#define opencc OPENCC
#define flshcc FLSHCC
#define closcc CLOSCC
#define plot   PLOT
#define erase  ERASE
#define zapbox ZAPBOX
#define factor FACTOR
#define newpen NEWPEN
#define where  WHERE
#define symbol SYMBOL
#define number NUMBER
#define xwinsz XWINSZ
#endif /* ifdef CAPS */

#ifdef UNDERSCORE
#define opencc opencc_
#define flshcc flshcc_
#define closcc closcc_
#define plot   plot_
#define erase  erase_
#define zapbox zapbox_
#define factor factor_
#define newpen newpen_
#define where  where_
#define symbol symbol_
#define number number_
#define xwinsz xwinsz_
#endif /* ifdef UNDERSCORE */

/* all functions in this file */
extern void opencc();
extern void killCC();
extern void flshcc();
extern void closcc();
extern void plot();
extern void erase();
extern void zapbox();
extern void factor();
extern void newpen();
extern void where();
extern void symbol();
extern void number();
extern int nonblank();
extern void xwinsz();

#ifdef ANSI
char *getenv(const char *name );
void exit(int status);
#else
char *getenv( );
#endif

/*  ----------------------------- OPENCC -------------------------------
             REAL XMAX,YMAX
             INTEGER ROTATE,WNAME(4)
             DATA WNAME/4HMy p,4Hrogr,4Ham  ,0/
             CALL OPENCC(XMAX,YMAX,ROTATE,WNAME)

    Open the CALCOMP-to-Xwindows emulator.
    See also FLSHCC and CLOSCC.

    XMAX and YMAX are the maximum value in "inches" in which your
    program will draw.  This code automatically converts these to
    the pixels in which X thinks.
    ROTATE should be 0 to select the normal orientation of CALCOMP,
    or 1 to select a 90 degree rotated version.   For more information
    about axis orientation, see the comments just above the end of
    this routine.
    WNAME is an integer array containing the name to be used for the
    window.  The last element of the array must be zero to terminate
    the string.  Integers seem to be more portable than passing a
    FORTRAN CHARACTER variable into this quiche.
*/
void opencc(xmax, ymax, rotate, wname)
float *xmax, *ymax;
int   *rotate;
char  *wname;
{
int   xcorner, ycorner, border;
int   depth, rootw, class, vmask;
int   x, y, win_x, win_y;
float xsize, ysize;
unsigned int mask;
char  macx[32], *icon_name;
Window child, root, focus_window;

   flip= *rotate;
   xmx = *xmax;
   ymx = *ymax;

/*   
     At ISU, we have some MacIntosh computers which simulate
     X-windows by running a program called MACX.  If the ISU user
     is on a MAC, the environment variable MACX will contain
     the string "true".  Since the MAC screen is ridiculously 
     small, we must decrease the window size, and the font below.
     Anything but a MACX screen gets treated normally.
*/
   if (getenv("MACX") != NULL) {
      strcpy(macx, getenv("MACX"));
   } else {
      strcpy(macx, "false");
   }
   if (strcmp(macx,"true") == 0) {
      width  =400;
      height =400;
      xcorner=100;
      ycorner=100;
      border =  3;
   } else {
      width  =640;
      height =640;
      xcorner=300;
      ycorner=100;
      border =  3;
   }

   factr  = 1.0;

   colors[0] = "white";
   colors[1] = "black";
   colors[2] = "red";
   colors[3] = "green";
   colors[4] = "blue";
   colors[5] = "cyan";
   colors[6] = "magenta";
   colors[7] = "yellow";
   colors[8] = "plum";
   colors[9] = "orange";
   colors[10] = "gray";
   colors[11] = "navy";
   colors[12] = "violet";
   colors[13] = "brown";
   colors[14] = "pink";
   colors[15] = "sea green";

   if (flip == TRUE) {
      xsize = width  / ymx;
      ysize = height / xmx;
      scale = xsize;
      if (ysize < scale) scale = ysize;
      width = scale * ymx;
      height= scale * xmx;
   } else {
      xsize = width  / xmx;
      ysize = height / ymx;
      scale = xsize;
      if (ysize < scale) scale = ysize;
      width = scale * xmx;
      height= scale * ymx;
   }
/* 
   keep real size of window to Height and Width
   height and width will be different from the real size of window
   because it determined to take a maximum region with a specific
   aspect in the window after window will be resized
*/
   Height = height;
   Width  = width;

/*
       Open the display.
       We do not "synchronize" here, because this slows down the
       plotting.  However, it means the calling program must be
       sure to CALL FLSHCC when the entire plot has been drawn.
       Otherwise, you will not see the last part of the plot.
       If you are debugging, synchronization ought to be on, and
       in addition, use -DDIRECTDRAW to avoid use of the pixmap.
*/
   display = XOpenDisplay(0);
   if (display == 0) {
      printf("OPENCC is unable to open a display.\n");
      printf("Most likely, you are not running Xwindows.\n");
      exit(-1);
   }
   XSynchronize(display, FALSE);

   screen = XDefaultScreenOfDisplay(display);
   visual = XDefaultVisualOfScreen(screen);
   depth  = XDefaultDepthOfScreen(screen);
   rootw  = XRootWindowOfScreen(screen);
   cmap   = XDefaultColormapOfScreen(screen);
   class  = InputOutput;
/*
       Locate window which currently has the input focus, 
       for this is the parent text window.
*/
   XQueryPointer(display, RootWindow(display,screen_num), 
                 &root, &child, &x, &y, &win_x, &win_y, &mask);
   if (child == None) focus_window = RootWindow(display,screen_num);
                 else focus_window = child;
/*
       Establish a window in the display.
       Set override_redirect to FALSE/TRUE if you want/do not want
       the window manager to interfere with this window.  The current
       attempt to leave input focus in the parent text window does
       not work so very well using FALSE.
*/
   xswa.event_mask        = ExposureMask;
   xswa.background_pixel  = WhitePixelOfScreen(screen);
   xswa.border_pixel      = BlackPixelOfScreen(screen);
   xswa.backing_store     = WhenMapped;
   xswa.override_redirect = TRUE;
   vmask  = (CWEventMask | CWBackPixel | CWBorderPixel | 
             CWOverrideRedirect | CWBackingStore);
   window = XCreateWindow(display, rootw, xcorner, ycorner, width, height, 
                          border, depth, class, visual, vmask, &xswa);

/*
     Is this a monochrome screen?

This is very hard to understand.  If you execute the following

      printf("There are %d color planes.\n",depth);
      switch (visual->class) {
         case PseudoColor: printf("case PseudoColor\n");  break;
         case StaticColor: printf("case Static Color\n"); break;
         case DirectColor: printf("case DirectColor\n");  break;
         case TrueColor:   printf("case TrueColor\n");    break;
         case GrayScale:   printf("case GrayScale\n");    break;
         case StaticGray:  printf("case StaticGray\n");   break;
      }

the following results are printed.  "color" is a decision made by
a human being standing in front of the screen:
                      color planes visual class
     DECstation 3100    no    1    StaticGray
     DECstation 5000    no    8    StaticGray
     VAXstation 3200    no    4    PseudoColor
     DECstation 3100   yes    8    PseudoColor
     IBM RS/6000 5081  yes    8    PseudoColor
     SGI Iris 4D/35    yes    8    PseudoColor
     XPac board in PC  yes    4    StaticColor
     Ardent TITAN      yes   24    DirectColor
so it appears that neither the class or the number of planes is
directly related to whether your computer has a color screen or not.

This matters, because many colors are invisible when you attempt
to use colors on a monochrome screen.  We want to draw strictly in
black if the screen does not do colors.

The following code attempts to decide correctly if the screen does
not support color, and works correctly for all known cases (so far).
*/
   monochrome = FALSE;
   if (visual->class == StaticGray) monochrome=TRUE;
   if (depth <= 4 && visual->class != StaticColor) monochrome=TRUE;
/*
        Establish the graphics context of the window
*/
   vmask  = (GCForeground | GCBackground);
   xgcv.foreground = BlackPixelOfScreen(screen);
   xgcv.background = WhitePixelOfScreen(screen);
   gc     = XCreateGC(display,window,vmask,&xgcv);
  
/*
        Establish fonts, and name the window
*/
   if (strcmp(macx,"true") == 0) {
      font = XLoadFont(display,Font2);
   } else {
      font = XLoadFont(display,Font1);
   }
   XSetFont(display,gc,font);
   XStoreName(display,window,wname);
/*
        Finally!  Put the danged window on the screen
*/
   XMapWindow(display,window);
   XClearWindow(display,window);
   XFlush(display);
#ifndef DIRECTDRAW
   off_screen = XCreatePixmap(display, window, width, height, depth);
#endif
   ccopen = TRUE;
   erase();
/*
        In the normal CALCOMP orientation, the origin is in the
        lower left corner, x runs across, y runs up.  This makes
        the window seem like a piece of graph paper.

        In the rotated extension to CALCOMP, the origin is in
        the upper left corner, x runs down, y runs across.  This
        is like a piece of graph paper rotated 90 degrees clockwise.

        In X, the origin is in the upper left corner, with x
        running across, and y running down.  The code like the
        following throughout this file converts to this illogical
        orientation when computing pixels.
*/
   if (flip == TRUE) {
      xsave = 0;
      ysave = 0;
      xoff  = 0;
      yoff  = 0;
   } else {
      xsave = 0;
      ysave = height;
      xoff  = 0;
      yoff  = 0;
   }
/*
     Now we wait until the window manager puts the window
     up on the screen, which is known as the "expose event".
     We'd like to return the input focus to the parent text window,
     so that FORTRAN I/O statements will read from that window.
*/
   for ( ; ; ) {
      XNextEvent(display,&event);
      switch (event.type) {
         case Expose: 
            if (xswa.override_redirect == FALSE) {
               XSetInputFocus(display, focus_window, 
                              RevertToParent, CurrentTime);
            }
            return; 
            break;
         default: break;
      }
   }
}

/*  ----------------------------- KILLCC -------------------------------
    This routine is not callable from FORTRAN.
*/
void killCC()
{
   printf("You must call OPENCC to open the Calcomp emulator!\n");
   printf("\n");
   printf("A typical program calls OPENCC to open a window,\n");
   printf("then calls PLOT, SYMBOL, NUMBER, ... to draw the plot,");
   printf("and finally calls FLSHCC and CLOSCC to quit.\n");
   printf("\n");
   printf("Your program dies now!\n");
   exit(-1);
}

/*  ----------------------------- FLSHCC -------------------------------
             CALL FLSHCC
    Flush the Xwindows output buffer so you can see all the plot.
    It is necessary to call this routine once your plot is completed
    so that the Xserver is forced to finish displaying the drawing.
*/
void flshcc()
{

   XCopyArea(display, off_screen, window, gc, 0, 0, Width, Height, 0, 0);
   XFlush(display);
   return;
}

/*  ----------------------------- CLOSCC -------------------------------
             CALL CLOSCC
    Close CALCOMP emulation and destroy the window.
*/
void closcc()
{
   if (ccopen == FALSE) killCC();

#ifndef DIRECTDRAW
   XFreePixmap(display, off_screen);
#endif
   XUnmapWindow(display, window);
   XDestroyWindow(display, window);
   XCloseDisplay(display);
   ccopen=FALSE;
   return;
}

/*  ----------------------------- PLOT ---------------------------------
             REAL X,Y
             INTEGER IPEN
             CALL PLOT(X,Y,IPEN)

    Moves pen from its previous position to the new point X,Y.
    IPEN = +/- 2 means move with pen down (drawing a line).
    IPEN = +/- 3 means move with pen up (not drawing a line).

    If IPEN is negative, the origin is redefined to be at X,Y.

    All other values of IPEN are illegal.  Note that the widespread
    but not universal use of IPEN=999 to close CALCOMP is instead
    implemented herein by CLOSCC.
*/
void plot(x,y,ipen)
float *x,*y;
int   *ipen;
{
int   pen,xnew,ynew;

   if (ccopen == FALSE) killCC();
   pen = *ipen;
   if ((abs(pen) != 2) & (abs(pen) != 3)) {
      printf("An illegal value of %d was given to IPEN in PLOT.\n",*ipen);
      exit(-1);
   }

   if (flip == TRUE) {
      xnew =          (int) (scale * *y + 0.5) + xoff;
      ynew =          (int) (scale * *x + 0.5) + yoff;
   }
   else {
      xnew =          (int) (scale * *x + 0.5) + xoff;
      ynew = height - (int) (scale * *y + 0.5) - yoff;
   }
   if (abs(pen) == 2)
           XDrawLine(display, off_screen, gc, xsave, ysave, xnew, ynew);
   xsave = xnew;
   ysave = ynew;
   if (pen > 0) return;
/*
        establish a new origin at the current point
*/
   if (flip == TRUE) {
      xoff = xsave;
      yoff = ysave;
   }
   else {
      xoff =          xsave;
      yoff = height - ysave;
   }
   return;
}

/*  ----------------------------- ERASE --------------------------------
             CALL ERASE
    Blank out the entire window, restore origin to original position.
*/
void erase()
{
   if (ccopen == FALSE) killCC();
   /* to clear a pixmap, reverse foreground and background ... */
   XSetForeground(display, gc, WhitePixelOfScreen(screen));
   /* ...and fill rectangle the size of the pixmap */
   XFillRectangle(display, off_screen, gc, 0, 0, Width, Height);
   /* don't foreget to reset
      here, the foreground assumed to be set other place using
      xgcv. If it isn't, the foreground will be lost */
   XSetForeground(display, gc, xgcv.foreground);

   if (flip == TRUE) {
      xsave = 0;
      ysave = 0;
      xoff  = 0;
      yoff  = 0;
   }
   else {
      xsave = 0;
      ysave = height;
      xoff  = 0;
      yoff  = 0;
   }
   return;
}

/*  ----------------------------- ZAPBOX -------------------------------
             REAL X1,Y1,X2,Y2
             CALL ZAPBOX(X1,Y1,X2,Y2)
    Blank out the box whose lower left corner is (X1,Y1) and
    whose upper right corner is (X2,Y2).
*/
void zapbox(x1,y1,x2,y2)
float *x1,*y1,*x2,*y2;
{
float xx1,yy1,xx2,yy2;
int xnew,ynew,w,h;

   if (ccopen == FALSE) killCC();
   xx1 = *x1;
   yy1 = *y1;
   xx2 = *x2;
   yy2 = *y2;

   if (flip == TRUE) {
      xnew =          (int) (scale * yy1 + 0.5) + xoff;
      ynew =          (int) (scale * xx1 + 0.5) + yoff;
      w = (int) (scale * (yy2-yy1));
      h = (int) (scale * (xx2-xx1));
   }
   else {
      xnew =          (int) (scale * xx1 + 0.5) + xoff;
      ynew = height - (int) (scale * yy2 + 0.5) - yoff;
      w = (int) (scale * (xx2-xx1));
      h = (int) (scale * (yy2-yy1));
   }

   /* to clear a pixmap, reverse foreground and background ... */
   XSetForeground(display, gc, WhitePixelOfScreen(screen));
   /* ...and fill rectangle the size of the pixmap */
   XFillRectangle(display, off_screen, gc, xnew, ynew, w, h);
   /* don't foreget to reset
      here, the foreground assumed to be set other place using
      xgcv. If it isn't, the foreground will be lost */
   XSetForeground(display, gc, xgcv.foreground);

   return;
}

/*  ----------------------------- FACTOR -------------------------------
             REAL FAC
             CALL FACTOR(FAC)
    Scale the entire plot by a factor FAC.
    This option is not implemented, except to store the value.
*/
void factor(fac)
float *fac;
{
   if (ccopen == FALSE) killCC();
   factr = *fac;
   return;
}

/*  ----------------------------- NEWPEN -------------------------------
             INTEGER IPEN
             CALL NEWPEN(IPEN)
    Change the pen color to IPEN.  See the routine OPENCC for the
    mapping of integers IPEN to color names.
    We extend Calcomp such that IPEN=0 draws in the background color.
*/
void newpen(ipen)
int   *ipen;
{
int pen,vmask;

   if (ccopen == FALSE) killCC();
   pen = *ipen;
   if (pen < 0) return;
   if (pen > 15) return;

/*             code for monochrome display     */
   if (monochrome == TRUE) {
      if (pen > 1) return;
      vmask  = GCForeground;
      if (pen == 0) xgcv.foreground = WhitePixelOfScreen(screen);
      if (pen == 1) xgcv.foreground = BlackPixelOfScreen(screen);
      XChangeGC(display, gc, vmask, &xgcv);
      return;
   }

/*             code for color displays         */
   XAllocNamedColor(display, cmap, colors[pen], &screenC, &exactC);
   vmask  = GCForeground;
   xgcv.foreground = screenC.pixel;
   XChangeGC(display, gc, vmask, &xgcv);
   return;
}

/*  ----------------------------- WHERE --------------------------------
             REAL X,Y,FAC
             CALL WHERE(X,Y,FAC)
    Return current pen position X,Y and scale factor FAC.
*/
void where(x, y, fac)
float *x,*y,*fac;
{
   if (ccopen == FALSE) killCC();
   if (flip == TRUE) {
      *x =           (float) ysave / scale;
      *y =           (float) xsave / scale;
   }
   else {
      *x =           (float) xsave / scale;
      *y = (float) (height - ysave) / scale;
   }
   *fac = factr;
   return;
}

/*  ----------------------------- SYMBOL -------------------------------
             REAL X,Y,HITE,THETA
             INTEGER NCHAR,INFO(?)
             CALL SYMBOL(X,Y,HITE,INFO,THETA,NCHAR)
    For positive NCHAR:
    Draw the first NCHAR letters of the character string INFO at X,Y.
    The character string must be passed packed into an integer array!
    The HITE and angle THETA at which to draw are not implemented.

    For negative NCHAR:
    Draw the special marker INFO at X,Y.
    INFO equals 1,2,3,... means square, circle, triangle, ...
    If NCHAR= -1, just draw the marker.
    If NCHAR= -2, draw a line from the current pen location
                  to X,Y before drawing the marker there.
    Note that negative NCHAR's are not yet implemented!
*/
void symbol(x, y, hite, info, theta, nchar)
float *x, *y, *hite, *theta;
char  *info;
int   *nchar;
{
int horiz, xnew, ynew, nch;

   if (ccopen == FALSE) killCC();
   nch = *nchar;
   if (nch == 0) return;
/*
      Code to draw a character string stored in an integer array.
      The X function below draws horizontal strings at a constant height.
*/
   if (nch > 0) {
   horiz = FALSE;
   if (flip == TRUE) {
      if(fabs(*theta-90.0) < 5.0) horiz=TRUE;
      xnew =          (int) (scale * *y + 0.5) + xoff;
      ynew =          (int) (scale * *x + 0.5) + yoff;
   }
   else {
      if(fabs(*theta) < 5.0) horiz=TRUE;
      xnew =          (int) (scale * *x + 0.5) + xoff;
      ynew = height - (int) (scale * *y + 0.5) - yoff;
   }
   if (horiz == FALSE) return;
   XDrawString(display, off_screen, gc, xnew, ynew, info, nch);
   return;
}
/*
      Code to draw a special symbol.
      Note that this is not implemented.
*/
   if (nch < -2) return;
   return;
}

/*  ----------------------------- NUMBER -------------------------------
             REAL X,Y,HITE,FLOTE,THETA
             INTEGER NCHAR
             CALL NUMBER(X,Y,HITE,FLOTE,THETA,NCHAR)
    Draw the floating point value FLOTE with NCHAR digits after
    the decimal at X,Y, at angle THETA and height HITE.

    This implementation treats negative NCHAR, and NCHAR=0 as
    if the user wants to draw an integer.
*/
void number(x, y, hite, flote, theta, nchar)
float *x, *y, *hite, *flote, *theta;
int   *nchar;
{
int  len, nch, xnew, ynew, horiz;
char info[30], fmt[30];

   if (ccopen == FALSE) killCC();
   nch = *nchar;
   if (nch < 0) nch=0;
   sprintf(fmt,"%%-30.%df",nch);
   sprintf(info,fmt,*flote);
   len = nonblank(info);

   horiz = FALSE;
   if (flip == TRUE) {
      if(fabs(*theta-90.0) < 5.0) horiz=TRUE;
      xnew =          (int) (scale * *y + 0.5) + xoff;
      ynew =          (int) (scale * *x + 0.5) + yoff;
   }
   else {
      if(fabs(*theta) < 5.0) horiz=TRUE;
      xnew =          (int) (scale * *x + 0.5) + xoff;
      ynew = height - (int) (scale * *y + 0.5) - yoff;
   }
/*
      the X function below draws horizontal strings at a constant height
*/
   if (horiz == FALSE) return;
   XDrawString(display, off_screen, gc, xnew, ynew, info, len);
   return;
}
/*  ------------------------- nonblank --------------------------------
    Returns the number of nonblank characters in the first word in
    string S.  Stolen from K&R, page 98.  Not callable from FORTRAN.
*/
int nonblank(s)
char *s;
{
char *p;

    p=s;
    while ((*p != ' ') & (*p != '\0')) p++;
    return (p-s);
}
/*  ------------------------- xwinsz --------------------------------
*/
void xwinsz(hite,wid)
int *hite, *wid;
{
   XWindowAttributes  winatt;
   Status ierr;
   float xsize, ysize;
   int test;
   ierr = XGetWindowAttributes(display, window, &winatt);
   *hite = winatt.height;
   *wid = winatt.width;
/* compare with the size of pixmap (window),
   if different, change configuration of window */
   if ( (Height != winatt.height) || (Width != winatt.width) ) {
      Height = height = *hite;
      Width  = width  = *wid;
      if (flip == TRUE) {
         xsize = width  / ymx;
         ysize = height / xmx;
         scale = xsize;
         if (ysize < scale) scale = ysize;
         width = scale * ymx;
         height= scale * xmx;
      } else {
         xsize = width  / xmx;
         ysize = height / ymx;
         scale = xsize;
         if (ysize < scale) scale = ysize;
         width = scale * xmx;
         height= scale * ymx;
      }
/*    also keep aspect of the window
      XResizeWindow(display, window, width, height); */
      Width  = width;
      Height = height;
#ifndef DIRECTDRAW
      XFreePixmap(display, off_screen);
      off_screen = XCreatePixmap(display, window, Width, Height,
                             XDefaultDepthOfScreen(screen));
#endif
   }
   return;
}
