"""Run every Assimulo example and dump solver trajectories + statistics, so that two builds of
Assimulo (e.g. the numpy.distutils build and the meson-python build, or two compilers) can be
compared bit for bit.

    python tools/dump_examples.py OUT.npz          # dump with the currently installed assimulo
    python tools/dump_examples.py --compare A.npz B.npz

tests/test_build_equivalence.py runs the same comparison against the file named by the
environment variable ASSIMULO_REFERENCE_RESULTS.

Records per example: t_sol, y_sol (and yd_sol for implicit problems), the statistics dict, and
for Kinsol the solution vector. Examples that need a missing solver are skipped and listed.
"""
import importlib, io, json, os, sys, contextlib, traceback
import numpy as np



def run_all():
    """Runs the examples of the *installed* assimulo package (assimulo.examples)."""
    import matplotlib
    matplotlib.use("Agg")
    import assimulo.examples as ex
    examples_dir = os.path.dirname(ex.__file__)
    out, skipped = {}, {}
    for f in sorted(os.listdir(examples_dir)):
        if not f.endswith(".py") or f == "__init__.py":
            continue
        name = f[:-3]
        try:
            mod = importlib.import_module("assimulo.examples." + name)
            if not hasattr(mod, "run_example"):
                continue
            with contextlib.redirect_stdout(io.StringIO()):
                try:
                    res = mod.run_example(with_plots=False)
                except TypeError:
                    res = mod.run_example()
        except Exception as e:
            skipped[name] = "%s: %s" % (type(e).__name__, str(e)[:120])
            continue
        objs = res if isinstance(res, tuple) else (res,)
        rec = {}
        for o in objs:
            for attr in ("t_sol", "y_sol", "yd_sol"):
                if hasattr(o, attr) and getattr(o, attr) is not None and len(getattr(o, attr)):
                    rec[attr] = np.asarray(getattr(o, attr), dtype=float)
            if hasattr(o, "statistics"):
                st = o.statistics
                try:
                    rec["statistics"] = {k: float(st[k]) for k in st.statistics.keys()}
                except Exception:
                    pass
            if hasattr(o, "y") and "y_sol" not in rec:      # Kinsol: solution vector
                try:
                    rec["y"] = np.asarray(o.y, dtype=float)
                except Exception:
                    pass
        out[name] = rec
    return out, skipped


def save(path, out, skipped):
    flat = {}
    for name, rec in out.items():
        for k, v in rec.items():
            if k == "statistics":
                flat[name + "/statistics"] = np.array(json.dumps(v, sort_keys=True))
            else:
                flat[name + "/" + k] = v
    flat["__skipped__"] = np.array(json.dumps(skipped, sort_keys=True))
    import assimulo
    flat["__assimulo_file__"] = np.array(assimulo.__file__)
    np.savez_compressed(path, **flat)


def compare(a_path, b_path):
    A, B = np.load(a_path), np.load(b_path)
    print("A:", str(A["__assimulo_file__"]), "\nB:", str(B["__assimulo_file__"]))
    keys = sorted(set(A.files) | set(B.files))
    ident, differ, missing = 0, [], []
    for k in keys:
        if k.startswith("__"):
            continue
        if k not in A.files or k not in B.files:
            missing.append(k); continue
        a, b = A[k], B[k]
        if a.dtype.kind in "US":
            same = str(a) == str(b)
            if not same:
                differ.append((k, "statistics differ: %s vs %s" % (str(a)[:80], str(b)[:80])))
            else:
                ident += 1
            continue
        if a.shape != b.shape:
            differ.append((k, "shape %s vs %s" % (a.shape, b.shape))); continue
        if np.array_equal(a, b):
            ident += 1
        else:
            d = np.abs(a - b); rel = d.max() / (np.abs(a).max() + 1e-300)
            differ.append((k, "max abs diff %.3e, rel %.3e, %d of %d entries" % (d.max(), rel, int((d > 0).sum()), a.size)))
    print("identical arrays/statistics: %d; differing: %d; missing on one side: %d" % (ident, len(differ), len(missing)))
    for k, msg in differ:
        print("  DIFF %-45s %s" % (k, msg))
    for k in missing:
        print("  MISSING %s" % k)
    print("skipped A:", str(A["__skipped__"])); print("skipped B:", str(B["__skipped__"]))
    return len(differ) == 0 and len(missing) == 0


if __name__ == "__main__":
    if sys.argv[1] == "--compare":
        sys.exit(0 if compare(sys.argv[2], sys.argv[3]) else 1)
    out, skipped = run_all()
    save(sys.argv[1], out, skipped)
    print("dumped %d examples, skipped %d: %s" % (len(out), len(skipped), json.dumps(skipped, indent=1) if skipped else ""))
