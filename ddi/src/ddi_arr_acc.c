#include "ddi_base.h"

/* Subtle difference between DDI_Acc and DDI_ARR_Acc is scale factor alpha */
void DDI_ARR_acc_(DDI_Patch *dAPatch, double alpha, double *buf) {
  DDI_AccP(dAPatch->handle, dAPatch, alpha, buf);
}

void DDI_ARR_acc_server(DDI_Patch *dAPatch, int rank) {
  DDI_Acc_server(dAPatch, rank);
}
