#ifndef RADAU_DECSOL_H
#define RADAU_DECSOL_H

int radau5_c(integer *n, U_fp fcn, doublereal *x, doublereal *
            y, doublereal *xend, doublereal *h__, doublereal *rtol, doublereal *
            atol, integer *itol, U_fp jac, integer *ijac, integer *mljac, integer 
            *mujac, U_fp mas, integer *imas, integer *mlmas, integer *mumas, U_fp 
            solout, integer *iout, doublereal *work, integer *lwork, integer *
            iwork, integer *liwork, doublereal *rpar, integer *ipar, integer *
            idid);

doublereal contr5_c(integer *i__, doublereal *x, doublereal *cont, integer *
                   lrc);
                   
#endif