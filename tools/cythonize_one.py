"""Cythonize one .pyx for meson, with a compile-time environment (the SUNDIALS facts) that
cannot be passed on cython's command line (tuples split at commas).

    python tools/cythonize_one.py <in.pyx> <out.c> -m <full.module.name> [-I dir ...] [-E NAME=PYTHON_LITERAL ...]
"""
import ast, os, sys

from Cython.Compiler import Options
from Cython.Compiler.Main import CompilationOptions, compile_single, default_options

args = sys.argv[1:]
pyx, out = args[0], args[1]
include_path, env, module = [], {}, None
i = 2
while i < len(args):
    if args[i] == "-m":
        module = args[i + 1]; i += 2
    elif args[i] == "-I":
        include_path.append(args[i + 1]); i += 2
    elif args[i] == "-E":
        name, value = args[i + 1].split("=", 1)
        env[name] = ast.literal_eval(value); i += 2
    else:
        sys.exit("unknown argument %s" % args[i])

options = CompilationOptions(default_options)
options.include_path = include_path
options.output_file = out
options.compile_time_env = env
options.language_level = "3str"
options.compiler_directives = {"language_level": "3str"}
options.emit_linenums = False
result = compile_single(pyx, options, full_module_name=module)
if result.num_errors:
    sys.exit(1)
