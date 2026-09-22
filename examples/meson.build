# The examples ship inside the package (assimulo.examples); tests/test_examples.py runs them.
example_py = run_command(py, ['-c', 'import glob, os, sys; print("\\n".join(sorted(os.path.basename(f) for f in glob.glob(os.path.join(sys.argv[1], "*.py")))))',
                              meson.current_source_dir()], check: true).stdout().strip().split('\n')
example_files = []
foreach f : example_py
  example_files += files(f)
endforeach
py.install_sources(example_files, subdir: pkg / 'examples')
install_data('kinsol_ors_matrix.mtx', install_dir: py.get_install_dir() / pkg / 'examples')
