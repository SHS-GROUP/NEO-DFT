/*  15 Apr 02 - MWS

      CALCOMP-to-POSTSCRIPT emulation code.

      Written in April and May 1990 by
         Francis F. MUGUET
         SPQR Laboratory
         Department of Chemistry
         Texas Tech University
         Lubbock, TX 79409
         AGGWI@TTACS.BITNET
         aggwi@ttacs1.ttu.edu
      starting from the CALCOMP-to-Xwindows conversion code
      previously written by Mike Schmidt.  Modified by MWS
      just a little in Jan 1992 to bring them online.

      Further modified by Jerry Boatz in Jan 1992 at Edwards
      Air Force Base Phillips Lab to make use of color (only
      for color postscript printers that recognize the CMYK
      color space.)

      The following CALCOMP calls are (partially) implemented:
         traditional      herein
         ----------       ------
          PLOTS           OPENPS
            -             EJCTPS
          PLOT            PLOT       and CLOSPS
          FACTOR          FACTOR
          WHERE           WHERE
          SYMBOL          SYMBOL
          NUMBER          NUMBER
          NEWPEN          NEWPEN
            -             ZAPBOX
          SCALE             -
          LINE              -
          AXIS              -

       At the bottom of this file, there are two C language routines,
       BLOTPS, BALLPS, and the commented out FORTRAN DRWMOL.  These
       are Francis' ideas for improved contours, and for a MOLPLT
       that can do XW and PS from a single executable.  Good ideas,
       and someday they should get some more attention, so we keep them
       here for safekeeping.  They are not actually used at present.

       For more information on each routine, and its usage from FORTRAN,
       see the comments preceeding each.  Note that not all of the
       functionality of each routine has been provided in every case.

       VAX/VMS FORTRAN-to-C requires capital letters, with no trailing
       underscores for subroutine names in this file (cc -DCAPS).

*/
#define TRUE 1
#define FALSE 0

#include  <stdio.h>
#include  <math.h>

/* globally known variables */

static int   height,width,flip;
static int   xsave,ysave,xoff,yoff,psopen=FALSE;
static int   npages,eps;
static float scale,factr;

/*  pointer to the PostScript file  */
FILE *fp;

/*
    Convert function names to upper case or underscore format
*/

#ifdef CAPS
#define openps OPENPS
#define ejctps EJCTPS
#define nextps NEXTPS
#define closps CLOSPS
#define plot   PLOT
#define zapbox ZAPBOX
#define factor FACTOR
#define newpen NEWPEN
#define where  WHERE
#define symbol SYMBOL
#define number NUMBER
#endif /* ifdef CAPS */

#ifdef UNDERSCORE
#define openps openps_
#define ejctps ejctps_
#define nextps nextps_
#define closps closps_
#define plot   plot_
#define zapbox zapbox_
#define factor factor_
#define newpen newpen_
#define where  where_
#define symbol symbol_
#define number number_
#endif /* ifdef UNDERSCORE */

#ifdef ANSI
char *getenv(const char *name );
void exit(int status);
#else
char *getenv( );
#endif

/*  ----------------------------- OPENPS -----------------------------
             REAL XMAX,YMAX
             INTEGER ROTATE,WNAME(4)
             DATA WNAME/4HMy p,4Hrogr,4Ham  ,0/
             CALL OPENPS(XMAX,YMAX,ROTATE,WNAME)

    Open a Postscript output file, write down general specifications.

    XMAX and YMAX are the maximum value in "inches" in which your
    program will draw.

    WNAME is an integer array containing the name to be used for the
    window.  The last element of the array must be zero to terminate
    the string.  Integers seem to be more portable than passing a
    FORTRAN CHARACTER variable.

    ROTATE should be 0 to select the normal orientation of CALCOMP,
    or 1 to select a 90 degree rotated version.

    ------------------------------------------------------------------
        In PostScript, the origin is the lower left corner,
        x runs across  ( 0 to 612 points)   (0 to 8.5 inches)   width
        y runs up      ( 0 to 792 points)   (0 to 11  inches)   height

        A point = 1/72 inch = 1/18 Pitch
        The Pitch is the standard typographical unit = 1/4 inch
    ------------------------------------------------------------------

*/
void openps(xmax,ymax,rotate,wname)
float *xmax,*ymax;
int   *rotate;
char  *wname;
{
float xsize,ysize,xmx,ymx;
char *filename;

   flip= *rotate;
   xmx = *xmax;
   ymx = *ymax;

/*  write all coordinates in points   */
   width  =600;
   height =780;

   factr  = 1.0;

   if (flip == TRUE) {
      xsize =  width/ymx;
      ysize = height/xmx;
                               /*  scale = min (xsize,ysize)   */
      scale = xsize;
      if (ysize < scale) scale = ysize;
   } else {
      xsize =  width/xmx;
      ysize = height/ymx;
      scale = xsize;
      if (ysize < scale) scale = ysize;
   }
/*
     Open the output PostScript file.
     The name of the file is defined by the environment variable PSTNAM
*/
   filename=getenv("PSTNAM");
   fp=fopen(filename,"w");
/*
     Write the comments according to
     the Adobe 2.0 structuring conventions.
     These conventions are really tedious and somewhat useless,
     so just try.
*/
   eps = FALSE;   /*   if Encapsulated PostScript, change this line  */
   if (eps==FALSE)   {
     fprintf(fp,"%%! PS-Adobe-2.0 \n");
   } else {
     fprintf(fp,"%%! PS-Adobe-2.0  EPS \n");
   }
   fprintf(fp,"%%%%Title: %s ",filename);
   fprintf(fp,"from GAMESS graphics utility %s \n",wname);
   fprintf(fp,"%%%%Creator: Francis F. Muguet \n");
   fprintf(fp,"%%%%CreationDate:  April-May 90 \n");
   fprintf(fp,"%%%%For: any GAMESS user \n");
   fprintf(fp,"%%%%BoundingBox:  \n");
   fprintf(fp,"%%%%Pages: (atend) \n");
   fprintf(fp,"%%%%DocumentFonts : Times-Roman\n");
   fprintf(fp,"%%%%EndComments\n");
   fprintf(fp,"%%%%EndProlog\n");
   fprintf(fp," \n");
/*      eps always false for now */
   if (eps == FALSE)   {
     fprintf(fp,"%%  reset graphics state\n") ;
     fprintf(fp,"initgraphics\n") ;
     fprintf(fp,"erasepage\n") ;
/*
    This set of statments checks to see if the postscript device
    understands the "setcmykcolor" operation.  If not, it redefines
    "setcmykcolor" as a no-op.
*/
     fprintf(fp,"/setcmykcolor where \n");
     fprintf(fp,"{pop} \n");
     fprintf(fp,"{/setcmykcolor {pop pop pop pop} def \n");
     fprintf(fp,"} ifelse \n");

/*
    This statement uses the CMYK color space and black as
    the initial color
*/
     fprintf(fp," 0.0 0.0 0.0 1.0 setcmykcolor \n") ;
   }
   fprintf(fp,"/inch {72 mul} def \n");
/*
      Define short-hand operator to reduce the output file's size
*/
   fprintf(fp,"%% short-hand operators \n");
   fprintf(fp,"/C {curveto} def \n");
   fprintf(fp,"/M {moveto}  def \n");
   fprintf(fp,"/L {lineto  stroke newpath}  def \n");
   fprintf(fp,"/circle {   0 360 arc} def \n");
/*
        Define line types and  set line parameters
*/
   fprintf(fp,"%% define line types \n");
   fprintf(fp,"/solid   { []        0 setdash } def \n");
   fprintf(fp,"/dot     { [2]       1 setdash } def \n");
   fprintf(fp,"/dash    { [4]       0 setdash } def \n");
   fprintf(fp,"/dashdot { [0 4 3 4] 0 setdash } def \n");
   fprintf(fp,"%% set the line width \n");
   fprintf(fp,"3.0 setlinewidth \n");
   fprintf(fp,"%% set line end type  0 butt 1 round 2 proj. square \n");
   fprintf(fp,"1 setlinecap \n");
   fprintf(fp,"%% set line join  0 miter 1 round 2 bevel join \n");
   fprintf(fp,"1 setlinejoin \n");
/*
       Print the graph name in the upper corner, in italics
*/
/*
   fprintf(fp,"/Souvenir-LightItalic  findfont \n");
   fprintf(fp,"15 scalefont \n");
   fprintf(fp,"setfont \n");
   fprintf(fp," 15 750 moveto\n");
   fprintf(fp," ( GAMESS graphics program : ");
   fprintf(fp,wname);
   fprintf(fp,"  ) show newpath \n");
*/

   npages = 1;
   fprintf(fp,"%%\n");
   fprintf(fp,"%%Page: %d\n",npages);
   fprintf(fp,"%%\n");
/*
       Translate just a little bit to be within the paper sheet
       If call requests it, rotate entire graph.
*/
   fprintf(fp," 15 15 translate\n");
   if (flip == TRUE) {
      fprintf(fp,"  0 730 translate\n");
      fprintf(fp,"  -90 rotate \n");
   }
/*
       All done
*/
   xsave = 0;
   ysave = 0;
   xoff  = 0;
   yoff  = 0;
   psopen = TRUE;
   }

/*  ----------------------------- KILLPS -------------------------------
    This routine is not callable from FORTRAN.
*/
void killPS()
{
printf("You must call OPENPS to open the Calcomp to Postscript emulator! \n");
   printf("\n");
   printf("A typical program calls OPENPS, then calls PLOT, SYMBOL,\n");
   printf("NUMBER, ..., and finally calls CLOSPS to quit.\n");
   printf("\n");
   printf("Your program dies now!\n");
   exit(-1);
}

/*  ----------------------------- EJCTPS --------------------------------
         CALL EJCTPS
    Eject the current plot, and restore origin to original position.
    You must call this once per page of output.

*/
void ejctps()
{
   if (psopen == FALSE) killPS();
   fprintf(fp," showpage \n");
   xsave = 0;
   ysave = 0;
   xoff  = 0;
   yoff  = 0;
   return;
}

/*  ----------------------------- NEXTPS --------------------------------
         CALL NEXTPS
    Set up additional plots...
    Translate just a little bit to be within the paper sheet.
    If opening call requests it, rotate entire graph.
    If a scale factor is in effect, use it on this page.
*/
void nextps()
{
   npages = npages + 1;
   fprintf(fp,"%%\n");
   fprintf(fp,"%%Page: %d\n",npages);
   fprintf(fp,"%%\n");
   fprintf(fp,"15 15 translate\n");
   if (flip == TRUE) {
      fprintf(fp,"0 730 translate\n");
      fprintf(fp,"-90 rotate \n");
   }
   fprintf(fp,"3.0 setlinewidth \n");
   if (factr != 1.0) fprintf(fp,"%5.2f   %5.2f  scale \n",factr,factr);
   return;
}

/*  ---------------------------- CLOSPS -------------------------------
          CALL CLOSPS
     Write trailer comments, close the PostScript output file.
     Note that this call does not print the last page,
     your calling program must call EJCTPS before CLOSPS!
*/
void closps()
    {int i;
        fprintf(fp,"%%%%Trailer \n") ;
        fprintf(fp,"%%%%Pages: %d \n",npages);
        i=fclose(fp);
        if (i) printf("Error closing PS file");
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
    implemented herein by a call to CLOSPS.

*/
void plot(x,y,ipen)
     float *x,*y;
           int   *ipen;
{
int   pen,xnew,ynew;

   if (psopen == FALSE) killPS();

   pen = *ipen;
   if ((abs(pen) != 2) & (abs(pen) != 3)) {
      printf("An illegal value of %d was given to IPEN in PLOT.\n",*ipen);
      exit(-1);
   }

   xnew = (int) (scale * (float)*x + 0.5) + xoff;
   ynew = (int) (scale * (float)*y + 0.5) + yoff;

   if (abs(pen) == 2)  {
      fprintf(fp,"  %d   %d  M ",xsave,ysave);
      fprintf(fp,"  %d   %d  L \n",xnew,ynew);
   }
   xsave = xnew;
   ysave = ynew;
   if (pen > 0) return;
/*
        establish a new origin at the current point
*/
   xoff = xsave;
   yoff = ysave;
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

   if (psopen == FALSE) killPS();
/*  first print the drawing  to be blanked out !  */
   fprintf(fp," copypage\n");
/*  translate a little bit to have the next frame within the paper sheet */
/*   fprintf(fp," 15 200 translate\n"); */
   xx1 = *x1;
   yy1 = *y1;
   xx2 = *x2;
   yy2 = *y2;

      xnew = (int) (scale * (float) xx1 + 0.5) + xoff;
      ynew = (int) (scale * (float) yy2 + 0.5) + yoff;
      w = (int) (scale * (xx2-xx1));
      h = (int) (scale * (yy2-yy1));
/*  xnew, ynew coordinates of the upper left corner of the box  */
   fprintf(fp,"%% print a white rectangular box, a little larger \n");
   fprintf(fp," newpath   %d 3 sub  %d  moveto \n",xnew,ynew);
   fprintf(fp," %d 3 add  0         rlineto \n",w);
   fprintf(fp," 0         %d 3 sub  rlineto \n",-h);
   fprintf(fp," %d 3 sub  0         rlineto \n",-w);
   fprintf(fp," closepath \n");

   fprintf(fp," 0.0 0.0 0.0 1.0 setcmykcolor fill \n");

/*  set default back to black                      */

   fprintf(fp,"newpath 0.0 0.0 0.0 1.0 setcmykcolor \n");

   return;
}

/*  ----------------------------- FACTOR -------------------------------
             REAL FAC
             CALL FACTOR(FAC)
    Scale the entire plot by a factor FAC.

*/
void factor(fac)
float *fac;
{
   if (psopen == FALSE) killPS();
   factr = *fac;
   fprintf(fp,"%5.2f   %5.2f  scale \n",factr,factr);
   return;
}

/*  ----------------------------- NEWPEN -------------------------------
             INTEGER IPEN
             CALL NEWPEN(IPEN)
    Change the pen color to IPEN.

*/
void newpen(ipen)
int   *ipen;
{
int pen;

   if (psopen == FALSE) killPS();
   pen = *ipen;
   if (pen <  0) return;
   if (pen > 15) return;
   switch (pen) {
                         /*             white                     */
     case  0:  fprintf(fp," 0.0 0.0 0.0 0.0 setcmykcolor \n"); break;
                         /*             black                     */
     case  1:  fprintf(fp," 0.0 0.0 0.0 1.0 setcmykcolor \n"); break;
                         /*             red                       */
     case  2:  fprintf(fp," 0.0 1.0 1.0 0.0 setcmykcolor \n"); break;
                         /*             green                     */
     case  3:  fprintf(fp," 1.0 0.0 1.0 0.0 setcmykcolor \n"); break;
                         /*             blue                      */
     case  4:  fprintf(fp," 1.0 1.0 0.0 0.0 setcmykcolor \n"); break;
                         /*             cyan                      */
     case  5:  fprintf(fp," 1.0 0.0 0.0 0.0 setcmykcolor \n"); break;
                         /*             magenta                   */
     case  6:  fprintf(fp," 0.0 1.0 0.0 0.0 setcmykcolor \n"); break;
                         /*             yellow                    */
     case  7:  fprintf(fp," 0.0 0.0 1.0 0.0 setcmykcolor \n"); break;
                         /*             plum                      */
     case  8:  fprintf(fp," 0.3 0.0 0.7 0.0 setcmykcolor \n"); break;
                         /*             orange                    */
     case  9:  fprintf(fp," 0.0 0.7 1.0 0.0 setcmykcolor \n"); break;
                         /*             gray                      */
     case 10:  fprintf(fp," 0.0 0.0 0.0 0.5 setcmykcolor \n"); break;
                         /*             navy                      */
     case 11:  fprintf(fp," 1.0 1.0 0.3 0.0 setcmykcolor \n"); break;
                         /*             violet                    */
     case 12:  fprintf(fp," 0.3 0.3 0.0 0.0 setcmykcolor \n"); break;
                         /*             brown                     */
     case 13:  fprintf(fp," 0.3 0.7 1.0 0.0 setcmykcolor \n"); break;
                         /*             pink                      */
     case 14:  fprintf(fp," 0.0 0.7 0.0 0.0 setcmykcolor \n"); break;
                         /*             sea green                 */
     case 15:  fprintf(fp," 1.0 0.0 0.6 0.0 setcmykcolor \n"); break;
   }
   return;
}

/*  ----------------------------- WHERE --------------------------------
             REAL X,Y,FAC
             CALL WHERE(X,Y,FAC)
    Return current pen position X,Y and scale factor FAC.

*/
void where(x,y,fac)
float *x,*y,*fac;
{
   if (psopen == FALSE) killPS();
   *x = (float)xsave / scale;
   *y = (float)ysave / scale;
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

    HITE seems to me the font size, the height of the letters
    The angle THETA, the rotation angle.

 For negative NCHAR:
      Draw the special marker INFO at X,Y.
      INFO equals 1,2,3,... means square, circle, triangle, ...
      If NCHAR= -1, just draw the marker.
      If NCHAR= -2, draw a line from the current pen location
                    to X,Y before drawing the marker there.
      Note that negative NCHAR's are not yet implemented!
      I haven't seen any call in the code with NCHAR < 0,
      so just ignore this,for now

*/
void symbol(x,y,hite,info,theta,nchar)
float *x,*y,*hite,*theta;
char *info;
int *nchar;
{
int xnew,ynew,nch,fontsi,i;
float angle,fonth;
   if (psopen == FALSE) killPS();
   angle = *theta;
   fonth = *hite;
   nch = *nchar;
   if (nch == 0) return;
   if (nch > 0) {
      xnew = (int) (scale * (float)*x + 0.5) + xoff;
      ynew = (int) (scale * (float)*y + 0.5) + yoff;
/* First save current graphic state and
         put in position and then
         rotate the reference frame
      fprintf(fp,"%%  SYMBOL call \n");  */
      fprintf(fp," gsave \n");
      fprintf(fp,"   %d   %d  moveto \n",xnew,ynew);
      fprintf(fp,"  %.2g   rotate \n",angle);

/* Second set up the font size, in the FORTRAN calls
   the HITE value varies from  0.125 to 0.30 i.e a ratio 1 to 2.5,
   a variation a little too wide, take the sqrt of this ratio
   8 points seems minimal for legibility.                       */
      fonth = 8.0 * sqrt(fonth/0.125) ;
      fontsi = (int) (fonth +0.5) ;
      fprintf(fp," /Times-Roman findfont \n");
      fprintf(fp," %d scalefont setfont \n",fontsi);
      fprintf(fp," (  ");
      for (i=0; i < nch; i++)  fprintf(fp,"%c",info[i]) ;
      fprintf(fp,"  ) show newpath \n");

/* restore graphic state  */
      fprintf(fp," grestore \n");
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
void number(x,y,hite,flote,theta,nchar)
float *x,*y,*hite,*flote,*theta;
int   *nchar;
{
int  nch,xnew,ynew,fontsi;
float angle,fonth,picnb;
char fmt[30];

   picnb = *flote;
   if (psopen == FALSE) killPS();
   angle = *theta;
   fonth = *hite;
   nch = *nchar;
   if (nch < 0) nch=0;
/* compose the format string =  "(   %.Nf  ) show newpath \n "
   for further printing
*/
      sprintf(fmt,"( %%.%df ) show newpath \n ",nch);

      xnew = (int) (scale * (float)*x + 0.5) + xoff;
      ynew = (int) (scale * (float)*y + 0.5) + yoff;

/* First save current graphic state , position and
         rotate the reference frame
      fprintf(fp,"%% NUMBER call \n");   */
      fprintf(fp," gsave \n");
      fprintf(fp,"   %d    %d    moveto  \n",xnew,ynew);
      fprintf(fp,"  %.2g   rotate \n",angle);

/* Second set up the font size, in the FORTRAN calls
   the HITE value varies from  0.125 to 0.30 i.e a ratio 1 to 2.5,
   a variation a little too wide, take the sqrt of this ratio
   8 points seems minimal for legibility.                       */
      fonth = 8.0 * sqrt(fonth/0.125) ;
      fontsi = (int) (fonth +0.5) ;
      fprintf(fp," /Times-Roman findfont \n");
      fprintf(fp," %d scalefont setfont \n",fontsi);
      fprintf(fp,fmt,picnb);

/* restore graphic state   */
      fprintf(fp," grestore \n");
   return;
}
/*  ----------------------------- BLOTPS ---------------------------------
             DOUBLE PRECISION X,Y
             DIMENSION X(N),Y(N)
             INTEGER N,IDASH
             CALL BLOT(X,Y,N,IDASH)   in KONTRS.CODE

    This subroutine replaces directly the FORTRAN subroutine BLOT,
    of KONTRS, therefore it
    has to be linked before KONTRS as to replace it.

  This replacement is not necessary, but the plotting quality is improved
  WYSIWYG is not true here, What you get is better than you see !

  If IDASH > 0  continuous solid line
           = 0  dash-dot line     Zero contour line
           < 0  dot-dot line

   maximum size of the array in KONTRS is 1600

    use of current point and curveto

*/
void blotps(x,y,in,idash)
/*   double x[],y[];  */
   double *x,*y;
          int *in,*idash;
{
int   dash,nleft,i,j,n;
int xn[1600],yn[1600] ;
double xx,yy;

   if (psopen == FALSE) killPS();
   dash = *idash;
   n = *in;
   fprintf(fp,"%% BLOT call N=  %d   IDASH = %d  \n",n,dash);
/* if the contour has less than 4 points ! don't bother with it    */
   if (n < 4) return;

   fprintf(fp," newpath \n");
   if  (dash > 0)     {   fprintf(fp," solid \n"); }
   if  (dash == 0)    {   fprintf(fp," dashdot \n"); }
   if  (dash < 0)     {   fprintf(fp," dot \n"); }

   for (i=0;i<n;i++)  {
/*  keep a pointer programming style instead of an array style
    equivalent to : */
/*    xx = (double) x[i];  yy = (double) y[i] ; */
/*  Within the pointer arithmetic, the pointer will be incremented
    by i*sizeof(double), it is what the compiler does internally anyway
    when dealing with arrays     */
      xx = (double) *(x+i) ;  yy = (double) *(y+i) ;
      xn[i] = (int) (scale * xx + 0.5) + xoff;
      yn[i] = (int) (scale * yy + 0.5) + yoff;  }

/* set up first position and current point  */
      fprintf(fp,"  %d   %d  moveto \n",xn[0],yn[0]);
/*  draw a Bezier cubic curve between the current point
    and next 3 specified points, i.e 4 points, or 3 intervals  */

    i = 1;  nleft=3;
  while (nleft >= 3) {
    fprintf(fp," %d %d ",xn[i],yn[i]);
    fprintf(fp," %d %d ",xn[i+1],yn[i+1]);
    fprintf(fp," %d %d C \n",xn[i+2],yn[i+2]);
/* the last point becomes the current point  */
    i = i + 3;
    nleft = n - i;  }

/* if there not enough of three points to finish the contour
    either 2 or 1  then  do straight lines   */

  for (j=i;j<n; j++)  {
    fprintf(fp,"  %d  %d lineto \n",xn[j],yn[j]);  }

  fprintf(fp," stroke ");
  fprintf(fp," newpath \n");

/*  reset  default line type */
   fprintf(fp," solid \n") ;
   return;
}


/*  ----------------------------- BALLPS  ---------------------------------
             INTEGER K      atom number in GAMESS
             REAL X,Y       center of atom
             REAL RADIUS    radius
             INTEGER COLOR  color code
             INTEGER IOPT   sphere rendering option
             CALL BALLPS(K,X,Y,RADIUS,COLOR,IOPT)

    This subroutine is called from the modified subroutine
    DRWMOL of MOLPLT

    option = 1    as X
             2    X plus atom number printed on spheres
             3    halft tone volume rendering of spheres

*/
void ballps(nb,x,y,radius,color,iopt)
       int *nb;
        float *x,*y,*radius;
                       int *color,*iopt;
{
int  label,couleur,xnew,ynew,rayon,i;
float xx,yy;
int  option;
/*   graylevel=1 white   graylevel=0.0  black  */
float  graylevel;
int rspot;

   if (psopen == FALSE) killPS();

/*      option =3;  */
      option = *iopt;
      label = *nb;
      couleur = *color;  /* not used, so far */

      xx = (float) *x;
      yy = (float) *y;
      xnew = (int) (scale * xx + 0.5) + xoff;
      ynew = (int) (scale * yy + 0.5) + yoff;
      rayon = (int) (scale * (float)*radius + 0.5 );

      fprintf(fp,"%% BALLPS call  Option =  %d \n",option);
      fprintf(fp," gsave \n");
      fprintf(fp," newpath \n");

 if (option == 1 || option == 2)
   {
      fprintf(fp," %d %d   %d circle \n",xnew,ynew,rayon);
/*  fill with white color  */
      fprintf(fp," 1 setgray fill \n");
/*  enlarge a little the linewidth   */
      fprintf(fp," 3.0 setlinewidth \n");
      fprintf(fp," %d %d   %d circle \n",xnew,ynew,rayon);
      fprintf(fp," 0 setgray stroke newpath  \n");
/* optionnally print the atom number  */
      if (option == 2)
      {
/* adjust a little bit where to print the label  */
      fprintf(fp," %d 3  sub  %d 2 sub moveto \n",xnew,ynew);
/* set the font  */
      fprintf(fp,"/Souvenir-Demi findfont \n");
      fprintf(fp,"9 scalefont  setfont \n");
      fprintf(fp,"(%d) show newpath \n",label);
      }
   }
/*  personnal artistic touch of mine
    no optometrics in it !  */
if (option ==3)
 {
/*  draw disk 1 point by 1 point, increasing the gray scale by 0.02 */
/*  determine size of the bright spot  */
    rspot = (int)  (  (float) rayon/ 2.5 );
/* determine the darkest gray tone  */
    graylevel = 0.85 - (float) ( rayon - rspot ) * 0.04 ;
/* print disks from the larger to the smaller */
    for (i=rayon; i > rspot; i--)
    {
    graylevel = graylevel + 0.04;
    if (graylevel < 0 )  {
      fprintf(fp,"  0  setgray \n");   }
    else {
      fprintf(fp," %.2f   setgray \n",graylevel); }
      fprintf(fp," %d %d  %d circle fill \n",xnew,ynew,i);
    }
/* print the center almost white spot  */
  fprintf(fp, " 0.90  setgray \n ");
  fprintf(fp, " %d  %d  %d circle fill newpath \n",xnew,ynew,rspot);

/* print the perimeter in black */
  fprintf(fp," 0 setgray \n");
  fprintf(fp," %d %d  %d circle  stroke \n",xnew,ynew,rayon);
 }
  fprintf(fp," grestore newpath \n");
   return;
}

/*

C     --------------------------------------
      SUBROUTINE DRWMOL(ERROR,IBD,JBD,NBMAX)
C     --------------------------------------
      IMPLICIT REAL(A-H,O-Z)
      DIMENSION IBD(NBMAX),JBD(NBMAX)
      LOGICAL ERROR
C*FM contains the depth ordered Z and Label=atom number
      REAL ZORD(1:50)
      INTEGER LABEL(1:50)
      CHARACTER*2 ENVVAR
C*FM
      COMMON /MPDATA/ X0(50),Y0(50),Z0(50),DX0(50),DY0(50),DZ0(50),
     *                X(50),Y(50),Z(50),DX(50),DY(50),DZ(50),
     *                KIND(50),BSIZ(50),SQMASS(50),
     *                IZAT(20),SIZE(20),KOLORS(20),
     *                SIZEP,VIEW,SCALEM,BNDWID,BNDLEN,
     *                NATOMS,NKINDS,NBONDS,KOLOR,MODE,NOBOX
      PARAMETER (ZERO=0.0E+00)
C
C     ----- DRAW THE ATOMIC SPHERES, BONDS, AND NORMAL MODE -----
C
C              ESTABLISH ORIGIN IN MIDDLE OF 14X14 PLOT AREA
C
      ERROR=.FALSE.
      CALL PLOT(7.0,7.0,-3)
      CALL CENTER(BIGGST)
C
C        SCALE MOLECULE TO FIT AVAILABLE AREA, ADDING PERSPECTIVE.
C
C        IF A NORMAL MODE IS TO BE DRAWN, GENERATE POINT (DX,DY,DZ)
C        USING THE DISPLACEMENT VECTOR (DX0,DY0,DZ0) FROM THE ATOM
C        CENTER (X0,Y0,Z0).  TO ENHANCE VIEWING OF SMALL DISPLACEMENTS,
C        THE DISPLACEMENT IS MASS-WEIGHTED, AND STARTS AT THE ATOM
C        SPHERE'S SURFACE, RATHER THAN THE CENTER.
C
      SCALE=4.5/BIGGST
      VDIST = VIEW*SCALE
      DO 200 I=1,NATOMS
         BSIZ(I)=SIZE(KIND(I))*SCALE*VDIST/(VDIST-Z(I))
         Z(I)=Z0(I)*SCALE
         X(I)=X0(I)*SCALE*VDIST/(VDIST-Z(I))
         Y(I)=Y0(I)*SCALE*VDIST/(VDIST-Z(I))
         IF(MODE.NE.0) THEN
               SC = SCALEM*SQMASS(I)
               SIZ = BSIZ(I)
               IF(SCALEM.LT.ZERO) SIZ = -SIZ
               DX(I) = X0(I) + SIZ*DX0(I) + SC*DX0(I)
               DY(I) = Y0(I) + SIZ*DY0(I) + SC*DY0(I)
               DZ(I) = Z0(I) + SIZ*DZ0(I) + SC*DZ0(I)
               DZ(I)=DZ(I)*SCALE
               DX(I)=DX(I)*SCALE*VDIST/(VDIST-DZ(I))
               DY(I)=DY(I)*SCALE*VDIST/(VDIST-DZ(I))
            END IF
  200 CONTINUE
C*FM with laser printer and display ,
C*FM it is easier to draw spheres filled
C*FM with color or blank, erasing part of the spheres drawn earlier
C*FM from the subroutine BALL we can infer that if Z(K)> Z(I)
C*FM atom K is not hidden by atom I , therefore we draw atoms
C*FM with from lower Z first,  sort the array from lower to higher Z
       DO I=1,NATOMS
           ZORD(I) = Z(I)
           LABEL(I) = I
       END DO

 210  ISWAP=0
      DO 250 I=1,NATOMS-1
C*FM test
         IF (ZORD(I).LE.ZORD(I+1)) GOTO 250
C*FM if not, swap the values and set ISWAP
         ZSWAP=ZORD(I)
         LSWAP=LABEL(I)
         ZORD(I)= ZORD(I+1)
         LABEL(I)=LABEL(I+1)
         ZORD(I+1)=ZSWAP
         LABEL(I+1)=LSWAP
           ISWAP=1
 250  CONTINUE
C*FM test if interchange, if yes, go for another pass
      IF (ISWAP.EQ.1) GOTO 210
C*FM end of sort, LABEL() contains the atoms nb sorted
C*FM get the rendering option
       CALL GETENV('OPTION',ENVVAR)
       READ(UNIT=ENVVAR,FMT='(I2)') IOPT
C
C              DRAW A SPHERE FOR EACH ATOM
C
      DO 300 I=1,NATOMS
         KOL = KOLORS(KIND(LABEL(I)))
C*FM
         CALL BALLPS(LABEL(I),X(LABEL(I)),Y(LABEL(I)),
     *               BSIZ(LABEL(I)),KOL,IOPT)
  300 CONTINUE
C
C             FIND ALL DISTANCES SHORTER THAN BNDLEN, EXCEPT H-H
C
      IF(NBONDS.LE.0) THEN
         NBONDS = 0
         DO 410 I=2,NATOMS
            DO 400 J=1,I-1
               IF(IZAT(KIND(I))*IZAT(KIND(J)).EQ.1) GO TO 400
               IF(DSTNCE(I,J).GT.BNDLEN) GO TO 400
               NBONDS=NBONDS+1
               IF(NBONDS.GT.NBMAX) GO TO 890
               IBD(NBONDS)=I
               JBD(NBONDS)=J
  400       CONTINUE
  410    CONTINUE
      END IF
C
C            DRAW BONDS WHERE APPROPRIATE
C            SKIP BONDS IF THE WIDTH IS ZERO.
C
      CALL NEWPEN(KOLOR)
      IF(BNDWID.EQ.ZERO) GO TO 600
C
C            A BOND CONSISTS OF NLINES STRAIGHT LINES,
C            HILDERBRANDT'S VERSION HAD NLINES=20
C
      NLINES=7
      DO 580 IJ=1,NBONDS
         I=IBD(IJ)
         J=JBD(IJ)
         DELTX=X(J)-X(I)
         DELTY=Y(J)-Y(I)
         DELTZ=Z(J)-Z(I)
         DD=SQRT(DELTX*DELTX+DELTY*DELTY)
         IF(DD.LT.BSIZ(I)) GO TO 580
         IF(DD.LT.BSIZ(J)) GO TO 580
         DELTXYZ=SQRT(DD*DD+DELTZ*DELTZ)
         IF(DELTZ.LT.ZERO) THEN
               DELTXI=BSIZ(I)*DELTX/DD
               DELTYI=BSIZ(I)*DELTY/DD
               DELTXJ=BSIZ(J)*DELTX/DELTXYZ
               DELTYJ=BSIZ(J)*DELTY/DELTXYZ
            ELSE
               DELTXI=BSIZ(I)*DELTX/DELTXYZ
               DELTYI=BSIZ(I)*DELTY/DELTXYZ
               DELTXJ=BSIZ(J)*DELTX/DD
               DELTYJ=BSIZ(J)*DELTY/DD
            END IF
         DO 560 LL=1,NLINES
            BDSZ=BNDWID*FLOAT(LL-1)/FLOAT(NLINES-1)
            DDI=BDSZ*SCALE*VDIST/(VDIST-Z(I))
            DDJ=BDSZ*SCALE*VDIST/(VDIST-Z(J))
            X1=X(I)+DELTXI-DDI*DELTY/DD
            Y1=Y(I)+DELTYI+DDI*DELTX/DD
            X2=X(J)-DELTXJ-DDJ*DELTY/DD
            Y2=Y(J)-DELTYJ+DDJ*DELTX/DD
            CALL STICK(X1,Y1,X2,Y2,I,J)
            X3=X(J)-DELTXJ+DDJ*DELTY/DD
            Y3=Y(J)-DELTYJ-DDJ*DELTX/DD
            CALL STICK(X2,Y2,X3,Y3,J,J)
            X4=X(I)+DELTXI+DDI*DELTY/DD
            Y4=Y(I)+DELTYI-DDI*DELTX/DD
            CALL STICK(X3,Y3,X4,Y4,J,I)
            X5=X(I)+DELTXI-DDI*DELTY/DD
            Y5=Y(I)+DELTYI+DDI*DELTX/DD
            CALL STICK(X4,Y4,X5,Y5,I,I)
  560    CONTINUE
  580 CONTINUE
C
C        SHOW NORMAL MODE DISPLACEMENTS
C
  600 CONTINUE
      IF(MODE.EQ.0) GO TO 700
      NLINES=3
      DO 680 I=1,NATOMS
         DELTX=DX(I)-X(I)
         DELTY=DY(I)-Y(I)
         DELTZ=DZ(I)-Z(I)
         DD=SQRT(DELTX*DELTX+DELTY*DELTY)
         IF(DD.LT.BSIZ(I)) GO TO 680
         DELTXYZ=SQRT(DD*DD+DELTZ*DELTZ)
         IF(DELTZ.LT.ZERO) THEN
               DELTXI=BSIZ(I)*DELTX/DD
               DELTYI=BSIZ(I)*DELTY/DD
            ELSE
               DELTXI=BSIZ(I)*DELTX/DELTXYZ
               DELTYI=BSIZ(I)*DELTY/DELTXYZ
            END IF
         DO 660 LL=1,NLINES
            BDSZ=BNDWID*FLOAT(LL-1)/FLOAT(NLINES-1)
            DDI=BDSZ*SCALE*VDIST/(VDIST-Z(I))
            X1=X(I)+DELTXI-DDI*DELTY/DD
            Y1=Y(I)+DELTYI+DDI*DELTX/DD
            X2=DX(I)
            Y2=DY(I)
            CALL STICK(X1,Y1,X2,Y2,I,I)
            X3=X2
            Y3=Y2
            X4=X(I)+DELTXI+DDI*DELTY/DD
            Y4=Y(I)+DELTYI-DDI*DELTX/DD
            CALL STICK(X3,Y3,X4,Y4,I,I)
            X5=X(I)+DELTXI-DDI*DELTY/DD
            Y5=Y(I)+DELTYI+DDI*DELTX/DD
            CALL STICK(X4,Y4,X5,Y5,I,I)
  660    CONTINUE
  680 CONTINUE
C
C        PUT ORIGIN BACK INTO LOWER RIGHT CORNER OF PLOT AREA
C
  700 CONTINUE
      CALL PLOT(-7.0,-7.0,-3)
      RETURN
C
C        EXCEEDED DIMENSION OF BOND ARRAYS
C
  890 CONTINUE
      WRITE(IW,900) I,J,NBMAX
      ERROR=.TRUE.
      RETURN
  900 FORMAT(1X,'Too many bonds, I,J=',2I6,' redimension NBMAX=',2I6)
      END
*/
