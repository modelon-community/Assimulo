/*
    Copyright (C) 2022 Modelon AB, all rights reserved.
*/

#ifndef _ODE_EVENT_LOCATOR_H
#define _ODE_EVENT_LOCATOR_H

#define ID_NO_EVENT 0
#define ID_EVENT    2

typedef int (*FP_event_func)(int, int, double, double*, double*, void*);
typedef int (*FP_interpolation)(int, double, double*, void*);

int f_event_locator(int n_y, int n_g, double TOL, double t_low, double *t_high,
                    double *y_high, double *g_old, double *g_mid, double *g_high, 
                    FP_event_func f_event, void *f_event_EXT,
                    FP_interpolation f_interp, void *f_interp_EXT,
                    int *f_event_cnt);

#endif /* _ODE_EVENT_LOCATOR_H */
