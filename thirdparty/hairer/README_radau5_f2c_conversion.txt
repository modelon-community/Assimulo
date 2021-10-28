Instructions for conversion of Radau5 (radau5_decsol.f) from Fortran to C via f2c:

Running f2c on radau5_decsol.f gives a similar issue as described here:

http://computer-programming-forum.com/49-fortran/1ac16746aa2d7d96.htm
https://stat.ethz.ch/pipermail/r-devel/2002-February/023967.html

The culprint is the "WERR" variable in the "RADCOR" subroutine. This can be fixed (to be tested/confirmed) by passing WERR as an additional argument into the "RADCOR" subroutine. This should enable the conversion from .f to .c code.

Afterwards, in the .c file: 

Remove the resulting extra function parameter of radcor_ in the resulting .c file and fix the function calls of radcor_ accordingly. 
In line 980 ish, replace

--werr;

by

doublereal *werr = (doublereal*) malloc(*n * sizeof(doublereal));

(This also requires including stdlib.h)

Finally, rename the following functions:

radau5_ -> radau5_c
contr5_ -> contr5_c

(This is meant to avoid name conflicts with the corresponding Fortran functions, when imported via Python. This should be fixable by other means on the Python side as well?)