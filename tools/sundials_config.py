"""Print SUNDIALS build facts for meson: version, index size, SuperLU_MT and rtol-vector support.

    python tools/sundials_config.py <include dir>

Output (one line): major.minor.patch index_bits superlu(0|1) rtolvec(0|1)
Mirrors what setup.py's check_SUNDIALS() reads from sundials/sundials_config.h.
"""
import os, sys

inc = sys.argv[1]
path = os.path.join(inc, "sundials", "sundials_config.h")
if not os.path.exists(path):
    sys.exit("no sundials_config.h under %s" % inc)
version, bits, superlu, rtolvec = None, "32", 0, 0
with open(path) as f:
    for line in f:
        if line.startswith("#define") and ("SUNDIALS_PACKAGE_VERSION" in line or "SUNDIALS_VERSION " in line):
            version = line.split()[-1].strip('"').split("-dev")[0]
        elif line.startswith("#define") and "SUNDIALS_INT64_T" in line:
            bits = "64"
        elif line.startswith("#define") and "SUNDIALS_INT32_T" in line:
            bits = "32"
        elif line.startswith("#define") and "SUNDIALS_SUPERLUMT" in line:
            superlu = 1
        elif line.startswith("#define") and "SUNDIALS_CVODE_RTOL_VEC" in line:
            rtolvec = 1
if version is None:
    sys.exit("no SUNDIALS version in %s" % path)
print(version, bits, superlu, rtolvec)
