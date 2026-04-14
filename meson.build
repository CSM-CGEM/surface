project(
    'meson-fortran-example',
    'cpp',
    version: '2026.0.0',
    default_options: [
        'cpp_std=c++20',
        'b_ndebug=if-release',
        'b_vscrt=from_buildtype',
        'default_library=shared',
    ],
    meson_version: '>=1.3',
)

inc_dir = include_directories('src')

surface_lib = library(
    'surface',
    'src/surface.cpp',
    include_directories:inc_dir,
    cpp_static_args: '-DSURFACE_STATIC',
    cpp_shared_args: '-DSURFACE_BUILD_DYNAMIC',
    gnu_symbol_visibility: 'hidden',
)

dep_args = []
if get_option('default_library') == 'static'
    dep_args += []
endif

surface_dep = declare_dependency(
    link_with: surface_lib,
    include_directories:inc_dir,
    compile_args:dep_args
)