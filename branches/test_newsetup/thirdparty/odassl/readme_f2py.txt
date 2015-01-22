F2py command

f2py -c odassl.pyf -m odassl *.f

* The file odassl.pyf is manually generated.

* Changes made to original ODASSL

  * RES signature changed in odacor and odajac
  * dfloat.f removed (obsolete)
