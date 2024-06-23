cmake_minimum_required(VERSION 3.1)
project(MyCGALProject)

set(CMAKE_BUILD_TYPE Release)

find_package(CGAL)

add_executable(dmr dmr.cpp)
target_link_libraries(dmr CGAL::CGAL)
