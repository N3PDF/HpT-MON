project(
  'higgsfo', 'cpp',
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
yaml_lib = cxx.find_library('yaml-cpp', required : true)

inc = include_directories(
  ['include']
)

subdir('include')
subdir('src')

higgsexe = executable(
  'higgs', 'main.cpp',
  include_directories : inc,
  dependencies : [glib_dep, yaml_lib],
  link_with : higgslib,
  c_args : extra_args
)

pkg_mod = import('pkgconfig')
pkg_mod.generate(
  libraries : higgslib,
  version : '1.0.0',
  name : 'higgsfo',
  description : 'A library for higgs FO calculations.'
)
