cmake_minimum_required(VERSION 3.1)
project(DMR)

set(CMAKE_BUILD_TYPE Release)

find_package(CGAL)

link_directories("/usr/lib/libteem.so.1.12.0")

add_executable(dmr dmr.cpp dmr_utilities.cpp dmr_io.cpp dmr_debug.cpp dmr_commandline.cpp dmr_cube.cpp dmr_delaunay.cpp dmr_globalvar.cpp dmr_grid.cpp dmr_voronoi.cpp)
target_link_libraries(dmr CGAL::CGAL teem)