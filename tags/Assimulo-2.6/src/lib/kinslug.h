/*
 * Copyright (C) 2010 Modelon AB / Copyright (c) 2002, The Regents of the University of California.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3 of the License.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.

 * You should have received a copy of the GNU Lesser General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * This file is a modification of the file kinsol_dense.h , revision 1.5,
 * from the SUNDIALS suite. The file is modified by:
 *
 * Johan Ylikiiskil� - johan.ylikiiskila@gmail.com
 *
 */

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
