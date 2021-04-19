project(
  'testhiggs', 'cpp',
  version : '1.0.0',
  default_options : ['cpp_std=c++14', 'warning_level=1'],
  license : 'MIT'
)

cxx = meson.get_compiler('cpp')

if cxx.get_id() == 'clang'
  extra_args = ['-fclang-flag']
else
  extra_args = []
endif

glib_dep = [dependency('lhapdf'), dependency('gsl')]
cuba_lib = cxx.find_library('cuba', required : true)

inc = include_directories(['../include'])

include_directories('../include')

install_headers('catch.hpp')

sources = [
  '../src/splittings.cpp',
  '../src/luminosity.cpp',
  '../src/utils.cpp',
  '../src/higgspt.cpp',
  '../src/higgsptpartonic.cpp',
  '../src/hadronic.cpp',
  '../src/partonic.cpp',
  '../src/integration.cpp'
]

testlib = static_library(
  'higgsfo',
  sources,
  include_directories : inc,
  dependencies : [glib_dep, cuba_lib],
  install: true
)

higgsexe = executable(
  'testhiggs', 'TestHiggs.cpp',
  include_directories : inc,
  dependencies : [glib_dep],
  link_with : testlib,
  c_args : extra_args
)
