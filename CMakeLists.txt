cmake_minimum_required(VERSION 3.1)

project(MathematicalOptimization CXX)

add_library(cmdline INTERFACE)
target_include_directories(cmdline SYSTEM INTERFACE
    "external/cmdline"
)

add_subdirectory(external/eigen)


add_library(report_interface INTERFACE)
target_compile_features(report_interface INTERFACE
    cxx_std_17
)
target_compile_options(report_interface INTERFACE
    $<$<CONFIG:Release>:-O3>
    $<$<CONFIG:Debug>:-Og>
    $<$<CONFIG:Debug>:-g>
    $<$<CONFIG:RelWithDebgInfo>:-O3>
    $<$<CONFIG:RelWithDebgInfo>:-g>

    "-ftemplate-backtrace-limit=0"
    "-fmax-errors=1"

    -Wall
    -Wextra
    -Wconversion
    -Wswitch-default
    -Wdisabled-optimization
    -Wformat
    -Winit-self
    -Woverloaded-virtual
    # -Wfloat-equal
    -Wno-old-style-cast
    -Wno-pragmas
    -Wno-ignored-attributes
    -Wreturn-type
    -Wno-cast-function-type
)
target_compile_definitions(report_interface INTERFACE
    $<$<CONFIG:Release>:EIGEN_NO_DEBUG>
)


add_subdirectory(report1)
