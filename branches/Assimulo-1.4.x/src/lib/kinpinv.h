
#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _PINV_H
#define _PINV_H

#include "kinsol_jmod.h"
#include <sundials/sundials_dense.h>

/*
 * -----------------------------------------------------------------
 * Function : KINpinv
 * -----------------------------------------------------------------
 * A call to the KINpinv function links the main solver with the
 * pinv linear solver. Its arguments are as follows:
 *
 * kinmem - pointer to an internal memory block allocated during a
 *          prior call to KINCreate
 *
 * N      - problem size
 *
 * The return value of KINpinv is one of:
 *    0                         if successful
 *    int different from zero   otherwise
 * -----------------------------------------------------------------
 */

  SUNDIALS_EXPORT int KINPinv(void *kinmem, int N);

#endif

#ifdef __cplusplus
}
#endif
