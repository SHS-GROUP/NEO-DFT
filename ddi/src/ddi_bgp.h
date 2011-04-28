#ifndef DDI_BGP_H
#define DDI_BGP_H

#include "ddi_base.h"
#include <stdio.h>

/*   GDF 10.06.09  DDI_BGL_File_IO_rank and DDI_BGL_File_IONode_rank were not used on L?  */

#define DDI_RUNTIME_PRINT_(s)    DDI_BGP_Runtime_print(s)

/** print DDI Blue Gene/P runtime */
void DDI_BGP_Runtime_print(FILE *stream);

#endif
