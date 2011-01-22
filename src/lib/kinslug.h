
#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _KINSLU_H
#define _KINSLU_H

#include "kinsol_jmod_wSLU.h"
#include <sundials/sundials_dense.h>

/*
 * -----------------------------------------------------------------
 * Function : KINSLUG
 * -----------------------------------------------------------------
 * A call to the KINSLUG function links the main solver with the
 * SLUG (SuperLU reGularization implementation)  linear solver.
 * Its arguments are as follows:
 *
 * kinmem - pointer to an internal memory block allocated during a
 *          prior call to KINCreate
 *
 * N      - problem size
 *
 * The return value of KINSLUG is one of:
 *    0                         if successful
 *    int different from zero   otherwise
 * -----------------------------------------------------------------
 */

  SUNDIALS_EXPORT int KINSLUG(void *kinmem, int N);

#endif

#ifdef __cplusplus
}
#endif
