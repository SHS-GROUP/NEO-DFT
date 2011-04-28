/* ----------------------------------------------------------------- *\
   GENERALIZED DISTRIBUTED DATA INTERFACE - 'MPMD' DATA SERVER MODEL
   =================================================================

   data_server.x program
   =====================

   DDI data_server program. Used to facilitate access to shared data 
   via point-to-point communications with remote compute processes.

   This is just a call to DDI_Init which branches to ddi_server 
   when the process rank (from ddikick) is >= nproc.

   Author:  Graham D. Fletcher  12.06.07

\* ----------------------------------------------------------------- */
 # include "mysystem.h"

  void DDI_Init(int, char **);

   int main(int argc, char *argv[]) { 
   DDI_Init(argc,argv);
   }
