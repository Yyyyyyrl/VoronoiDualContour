cmake_minimum_required(VERSION 3.1)
project(VDC)

set(CMAKE_BUILD_TYPE Release)

find_package(CGAL)

link_directories("/usr/lib/libteem.so.1.12.0")

add_executable(vdc vdc.cpp vdc_utilities.cpp vdc_io.cpp vdc_debug.cpp vdc_commandline.cpp vdc_cube.cpp vdc_delaunay.cpp vdc_globalvar.cpp vdc_grid.cpp vdc_voronoi.cpp)
target_link_libraries(vdc CGAL::CGAL teem)