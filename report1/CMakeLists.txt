cmake_minimum_required(VERSION 3.0)


add_executable(q1
    1/main.cpp
)
target_link_libraries(q1 PUBLIC
    report_interface
    cmdline
    eigen
)


add_executable(q2-iteration
    2/iteration.cpp
)
target_link_libraries(q2-iteration PUBLIC
    report_interface
    cmdline
    eigen
)

add_executable(q2-performance
    2/performance.cpp
)
target_link_libraries(q2-performance PUBLIC
    report_interface
    eigen
)


add_executable(q3-iteration
    3/iteration.cpp
)
target_link_libraries(q3-iteration PUBLIC
    report_interface
    cmdline
    eigen
)

add_executable(q3-performance
    3/performance.cpp
)
target_link_libraries(q3-performance PUBLIC
    report_interface
    eigen
)
