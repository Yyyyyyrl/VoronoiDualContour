cmake_minimum_required(VERSION 3.12)
project(CGAL_Comparison
    VERSION 1.0
    DESCRIPTION "CGAL Surface Mesher comparison tool for VDC"
    LANGUAGES CXX
)

# ─── Compiler settings ────────────────────────────────────────────────────────
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

# ─── Find dependencies ────────────────────────────────────────────────────────
find_package(CGAL REQUIRED COMPONENTS Core ImageIO)
find_package(ZLIB REQUIRED)

# Optional: Intel TBB for parallel meshing
find_package(TBB QUIET)
if(TBB_FOUND)
  message(STATUS "Found Intel TBB - parallel meshing available")
else()
  message(STATUS "Intel TBB not found - using sequential meshing")
endif()

# Teem doesn't ship a CMake config, so try:
#  1) an optional TEEM_ROOT hint
#  2) fallback to find_path/find_library
option(TEEM_ROOT "Root directory of your Teem install" "")

if(TEEM_ROOT)
  set(TEEM_INCLUDE_DIR  "${TEEM_ROOT}/include")
  set(TEEM_LIBRARY      "${TEEM_ROOT}/lib/libteem${CMAKE_SHARED_LIBRARY_SUFFIX}")
else()
  find_path(TEEM_INCLUDE_DIR NAMES teem/nrrd.h
            DOC "Where to find Teem headers (e.g. nrrd.h)")
  find_library(TEEM_LIBRARY NAMES teem
               DOC "Path to libteem")
endif()

if(NOT TEEM_INCLUDE_DIR OR NOT TEEM_LIBRARY)
  message(FATAL_ERROR
    "Could not locate Teem!
     Please install Teem or set -DTEEM_ROOT=/path/to/teem-prefix")
endif()

message(STATUS "Found Teem include: ${TEEM_INCLUDE_DIR}")
message(STATUS "Found Teem lib:     ${TEEM_LIBRARY}")

# ─── Include directories ──────────────────────────────────────────────────────
include_directories(${TEEM_INCLUDE_DIR})

# ─── Main executable ──────────────────────────────────────────────────────────
add_executable(cgal_mesher cgal_mesher.cpp)

# ─── Link libraries ───────────────────────────────────────────────────────────
target_link_libraries(cgal_mesher
    PRIVATE
      CGAL::CGAL
      CGAL::CGAL_ImageIO
      ZLIB::ZLIB
      ${TEEM_LIBRARY}
)

# Link TBB if available for parallel meshing
if(TBB_FOUND)
  target_link_libraries(cgal_mesher PRIVATE TBB::tbb)
  target_compile_definitions(cgal_mesher PRIVATE CGAL_CONCURRENT_MESH_3)
endif()

# ─── Installation (optional) ──────────────────────────────────────────────────
install(TARGETS cgal_mesher
        RUNTIME DESTINATION bin
)

# ─── Usage info ───────────────────────────────────────────────────────────────
message(STATUS "=================================================================")
message(STATUS "CGAL Comparison Tool Build Configuration")
message(STATUS "=================================================================")
message(STATUS "To build:")
message(STATUS "  cd cgal_comparison")
message(STATUS "  mkdir build && cd build")
message(STATUS "  cmake ..")
message(STATUS "  make")
message(STATUS "")
message(STATUS "To run:")
message(STATUS "  ./cgal_mesher <isovalue> <nrrd_file>")
message(STATUS "  ./cgal_mesher --help    (for all options)")
message(STATUS "=================================================================")
