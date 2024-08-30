cmake_minimum_required(VERSION 3.1)
project(DMR)

set(CMAKE_BUILD_TYPE Release)

find_package(CGAL)

link_directories("/usr/lib/libteem.so.1.12.0")

add_executable(dmr dmr.cpp utilities.cpp io.cpp debug.cpp)
target_link_libraries(dmr CGAL::CGAL teem)