/*
    Copyright (C) 2022 Modelon AB, all rights reserved.
*/

#include "ode_event_locator.h"

#define _abs(x) ((x) >= 0 ? (x) : -(x))
#define _max(a,b) ((a) >= (b) ? (a) : (b))
#define _min(a,b) ((a) < (b) ? (a) : (b))

/* Reference: "Discontinuities handled with events in Assimulo, a practical approach" by Emil Fredriksson */

/* Checks if an event occurs in [t_low, t_high], if that is the case event */
/* localization is started. Event localization finds the earliest small interval  */
/* that contains a change in domain. The right endpoint of this interval is then */
/* returned as the time to restart the integration at. */
int f_event_locator(int n_y, int n_g, double TOL_inp, double t_low, double *t_high,
                  double *y_high, double *g_low, double *g_mid, double *g_high, 
                  FP_event_func f_event, void *f_event_EXT,
                  FP_interpolation f_interp, void *f_interp_EXT,
                  int *f_event_cnt)
{
/* PARAMETERS */
/* n_y          Length of state vector */
/* n_g          Length of event indicator vector */
/* TOL_inp      Input tolerance used for event detection */
/* t_low        Lower bound for time-interval */
/* t_high       INPUT: Upper bound for time-interval */
/*              OUTPUT: Time of event */
/* y_high       IN & OUTPUT: (Length n_y) solution at *t_high */
/* g_low        (Length n_g) Event indicators at t_low */
/* g_mid        (Length n_g) Working array */
/* g_high       OUPUT: (Length n_g) output event indicators at t_high */
/* f_event      Function to obtain event indicators */
/* f_event_EXT  Additional input to f_event */
/* f_interp     Interpolation function */
/* f_interp_EXT Additional input to f_interp */
/* f_event_cnt  IN & OUTPUT: Counter for f_event evaluations */

/* RETURN, also see parameters marked with OUTPUT */
/* ID_NO_EVENT || ID_EVENT */
/* negative returns from any evaluations of f_event of f_interp are immediately returned. */
    int i, ret, side, sideprev, imax, event_left;
    double alpha, maxfrac, gfrac, t_mid, fracint;
    double TOL = _max(_abs(t_low), _abs(*t_high))*TOL_inp;

    /* evaluate event indicators at end-point of time-interval */
    ret = f_event(n_y, n_g, *t_high, y_high, g_high, f_event_EXT);
    *f_event_cnt += 1;
    if (ret < 0){ return ret;}

    /* Check for events in [t_low, t_high] */
    ret = ID_NO_EVENT;
    for (i = 0; i < n_g; i++){
        if ((g_low[i] > 0) != (g_high[i] > 0)){
            ret = ID_EVENT;
            break;
        }
    }
    
    /* no event */
    if (ret == ID_NO_EVENT){
        return ret;
    }

    /* event detected; localize */
    side = 0;
    sideprev = -1;
    event_left = -1;

    alpha = 1;
    while (_abs(*t_high - t_low) > TOL){
        /* Adjust alpha if the same side is choosen more than once in a row. */
        if (side == sideprev){
            if (side == 2){
                alpha *= 2.;
            }else{
                alpha /= 2.;
            }
        }else{ /* Otherwise alpha = 1 and the secant rule is used. */
            alpha = 1.;
        }
        /* Decide which event function to iterate with. */
        maxfrac = 0.;
        imax = 0;
        for (i = 0; i < n_g; i++){
            if ((g_low[i] > 0) != (g_high[i] > 0)){
                gfrac = _abs(g_high[i]/(g_low[i] - g_high[i]));
                if (gfrac >= maxfrac){
                    maxfrac = gfrac;
                    imax = i;
                }
            }
        }

        /* Hack for solving the slow converging case when g is zero for a large part of [t_low, t_high] */
        if ((g_high[imax] == 0) || (g_low[imax] == 0)){
            t_mid = (t_low + *t_high)/2.;
        }else{
            t_mid = *t_high - (*t_high - t_low)*g_high[imax]/(g_high[imax] - alpha*g_low[imax]);
        }

        /* Check if t_mid is to close to current brackets and adjust inwards if so is the case. */
        if (_abs(t_mid - t_low) < TOL/2.){
            fracint = _abs(t_low - *t_high)/TOL;
            t_mid = t_low + (*t_high - t_low)/(2.*_min(fracint, 5.));
        }

        if (_abs(t_mid - *t_high) < TOL/2.){
            fracint = _abs(t_low - *t_high)/TOL;
            t_mid = *t_high - (*t_high - t_low)/(2.*_min(fracint, 5.));
        }

        /* Calculate g at t_mid and check for events in [t_low, t_mid]. */
        ret = f_interp(n_y, t_mid, y_high, f_interp_EXT);
        if (ret < 0){ return ret;}
        ret = f_event(n_y, n_g, t_mid, y_high, g_mid, f_event_EXT);
        *f_event_cnt += 1;
        if (ret < 0){ return ret;}

        sideprev = side;
        event_left = 0; 
        for(i = 0; i < n_g; i++){
            if ((g_low[i] > 0) != (g_mid[i] > 0)){
                event_left = 1;
                break;
            }
        }
        if (event_left){ /* previous loop broken, event in [t_low, t_mid] */
            *t_high = t_mid;
            for(i = 0; i < n_g; i++){
                g_high[i] = g_mid[i];
            }
            side = 1;
        }else{ /* previous loop not broken */
            /* If there are no events in [t_low, t_mid] there must be some event in [t_mid, t_high]. */
            t_low = t_mid;
            for(i = 0; i < n_g; i++){
                g_low[i] = g_mid[i];
            }
            side = 2;
        }
    } /* while (_abs(*t_high - t_low) > TOL) */

    ret = f_interp(n_y, *t_high, y_high, f_interp_EXT); 
    if (ret < 0){ return ret;}

    return ID_EVENT;
}
        