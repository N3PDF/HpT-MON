sources = [
  'splittings.cpp',
  'luminosity.cpp',
  'utils.cpp',
  'higgspt.cpp',
  'higgsptpartonic.cpp',
  'hadronic.cpp',
  'partonic.cpp',
  'integration.cpp'
]

higgslib = static_library(
  'higgsfo',
  sources,
  include_directories : inc,
  dependencies : [glib_dep, cuba_lib, yaml_lib],
  install: true
)
