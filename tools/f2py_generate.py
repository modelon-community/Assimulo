"""Run f2py on a .pyf signature file for meson and make sure every declared output exists.

    python tools/f2py_generate.py <signature.pyf> <output dir> <module name>

f2py writes <module>module.c and, only when the wrapped routines need them,
<module>-f2pywrappers.f / <module>-f2pywrappers2.f90; meson's custom_target needs
a fixed list of outputs, so missing wrapper files are created empty (an empty
Fortran source compiles to an empty object).
"""
import os, subprocess, sys

pyf, outdir, module = sys.argv[1:4]
os.makedirs(outdir, exist_ok=True)
subprocess.check_call([sys.executable, "-m", "numpy.f2py", "--build-dir", outdir, pyf])
for name in ("%smodule.c" % module, "%s-f2pywrappers.f" % module, "%s-f2pywrappers2.f90" % module):
    p = os.path.join(outdir, name)
    if not os.path.exists(p):
        open(p, "w").close()
