Instructions for conversion of Radau5 (radau5_decsol.f) from Fortran to C via f2c:

Running f2c on radau5_decsol.f runs an issue that requires minor modification in the .f file. 
The culprint is the "WERR" variable in the "RADCOR" subroutine. The problem can be fixed by passing WERR as an additional argument into the "RADCOR" subroutine (in the .f file), to enable the f2c conversion.

Afterwards, in the .c file: 

Remove the resulting extra function parameter of radcor_ in the resulting .c file and fix the function calls of radcor_ accordingly. 

In line 980 ish, insert

doublereal *werr = (doublereal*) malloc(*n * sizeof(doublereal));
(This requires including stdlib.h)

Make sure the line "--werr;" line happens after this memory allocation.