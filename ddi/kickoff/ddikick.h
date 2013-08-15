/* ------------------------------------------------------------------------ *\
 *    GENERALIZED DISTRIBUTED DATA INTERFACE - 'MPMD' DATA SERVER MODEL     *
 *           Kickoff Program Prototypes and Data Structures                 *
 *     2 Aug 13 - DGF - increase number of TCP/IP connection attempts.
\* ------------------------------------------------------------------------ */
 # include "mysystem.h"
 # include "common.h"

 # include "debug.h"

/*   increased from 2 to 100 by DGF for large scale job kickoff (8/2013). */
/*   JAB could launch onto 50 nodes iff this was 100, but not if 10.      */
/*   The value's usage is as the number of attempts, not a time value.    */
 # ifndef DDIKICK_TIMEOUT
 # define DDIKICK_TIMEOUT 100
 # endif
 
 # ifdef USE_PBS
 # include "tm.h"

 # define TASK_UNBORN  0
 # define TASK_ALIVE   1
 # define TASK_DEAD    2
 # define RESC_STR_LEN 256

   struct ddi_task {
      tm_event_t event;
      tm_task_id id;
      int        obitval;
      int        state;
   };
 # endif

/* ---------------------------------------------------------- *\
   Macro used for global variable -- Ensures unique namespace
\* ---------------------------------------------------------- */
 # define gv(a) __ddikick_global__ ## a ## __

   typedef struct {
      int np;
      int nn;
   } DDI_Comm;

/* ---------------- *\
   Global variables
\* ---------------- */
   extern DDI_Comm gv(ddi_base_comm);
   extern int   gv(master);
   extern int   gv(nprocs);
   extern int   gv(nnodes)[1];
   extern int  *gv(sockets);
   extern int   gv(ppn);
   extern char *gv(netext);

/* -------------------------------------- *\
   Structure for command-line information
\* -------------------------------------- */
   typedef struct {
     int argc;
     char **argv;
     int nprocs;
     int nnodes;
     int ddiarg;
     int nodearg;
     int kickoffsock;
     int kickoffport;
     char *kickoffhost;
     int *ports;
     int *sockets;
   } Cmdline_info;


/* ------------------- *\
   Kickoff subroutines
\* ------------------- */
   void Kickoff_Local(int,int,const Node_info*,const Cmdline_info*);
   void Kickoff_Remote(int*,int*,const Node_info*,const Cmdline_info*,char*);
 # ifdef USE_PBS
   void Kickoff_PBS(const Node_info*,const Cmdline_info*);
 # endif


/* ----------------------------- *\
   DDI Task handling subroutines
\* ----------------------------- */
   void *Accept_tasks(void *);
   void  Finalize_tasks();


/* -------------------------- *\
   Error handling subroutines
\* -------------------------- */
   void Filicide();
   const char *Info();
   const char *DDI_Id();


