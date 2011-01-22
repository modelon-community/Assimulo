
#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _REG_ROUTINES_H
#define _REG_ROUTINES_H


#include "slu_ddefs.h"

  /*
   * Routines used in regularisation
   */

SuperMatrix* regSparseMatrix( SuperMatrix *jac, double h);

SuperMatrix* getRegRHS( SuperMatrix *jac, SuperMatrix *B);

double getRegParam(SuperMatrix *jac, SuperMatrix *B);
#endif

#ifdef __cplusplus
}
#endif
