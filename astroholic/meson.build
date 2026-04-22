py_subdir = 'astroholic'

py3.install_sources([
  '__init__.py',
  'lic.py',
  ], subdir: py_subdir)

fortran_module_name = '_lic'

# Check for OpenMP (optional dependency), then get OpenMP flags if found
openmp_dep = dependency('openmp', required : false)
openmp_flags = []
if openmp_dep.found()
  openmp_flags = ['-fopenmp']
endif

fortran_cwrapper = custom_target(fortran_module_name,
  input: ['lic.f90'],
  output: [fortran_module_name+'module.c', fortran_module_name + '-f2pywrappers2.f90'],
  command: [py3, '-m', 'numpy.f2py', '--lower', '-m', fortran_module_name,
  '@INPUT@', '--build-dir', '@OUTDIR@']
)
py3.extension_module(fortran_module_name,
  'lic.f90',
  fortran_cwrapper,
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies: [py_dep, openmp_dep],
  install: true,
  fortran_args: openmp_flags,
  subdir: py_subdir
)
