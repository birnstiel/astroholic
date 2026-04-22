py_subdir = 'astroholic'

py3.install_sources([
  '__init__.py',
  'lic.py',
  ], subdir: py_subdir)

fortran_module_name = '_lic'

# Check for OpenMP (optional dependency), then get OpenMP flags if found
fc = meson.get_compiler('fortran')
openmp_flags = []
openmp_link_flags = []
if fc.has_argument('-fopenmp')
  openmp_flags = ['-fopenmp']
  openmp_link_flags = ['-lgomp']
  message('OpenMP: YES (-fopenmp supported by ' + fc.get_id() + ')')
else
  message('OpenMP: NO (-fopenmp not supported)')
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
  c_args: np_inc_flags,
  dependencies: py_dep,
  install: true,
  fortran_args: openmp_flags + np_inc_flags,
  link_args: openmp_link_flags,
  subdir: py_subdir
)
