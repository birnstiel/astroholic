project('astroholic', 'c')
add_languages('fortran')

py_mod = import('python')
py3 = py_mod.find_installation(pure: false)
py_dep = py3.dependency()

incdir_numpy = run_command(py3,
  ['-c', 'import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py3,
    ['-c', 'import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

np_inc_flags = ['-I' + incdir_numpy, '-I' + incdir_f2py]

subdir('astroholic')
