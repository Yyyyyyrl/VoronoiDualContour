cmake_minimum_required(VERSION 3.12)
project(VDC
    VERSION 1.0
    DESCRIPTION "Dual‐contouring / Voronoi‐based mesh reconstruction"
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
find_package(CGAL REQUIRED COMPONENTS Core)
find_package(ZLIB REQUIRED)

# Teem doesn’t ship a CMake config, so try:
#  1) an optional TEEM_ROOT hint
#  2) fallback to find_path/find_library
option(TEEM_ROOT "Root directory of your Teem install" "")

if(TEEM_ROOT)
  set(TEEM_INCLUDE_DIR  "${TEEM_ROOT}/include")
  set(TEEM_LIBRARY      "${TEEM_ROOT}/lib/libteem${CMAKE_SHARED_LIBRARY_SUFFIX}")
else()
  find_path(TEEM_INCLUDE_DIR NAMES teem/ten.h
            DOC "Where to find teemplate headers (e.g. ten.h)")
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

# ─── Include & link directories ──────────────────────────────────────────────
# (CGAL::CGAL already brings in its include dirs)
include_directories(${TEEM_INCLUDE_DIR})

# ─── Common source files (all vdc*.cpp except vdc.cpp) ───────────────────────
file(GLOB_RECURSE COMMON_SOURCES
     "${CMAKE_CURRENT_SOURCE_DIR}/vdc_*.cpp"  # Matches vdc_* but not vdc.cpp
)

# ─── Main executable ─────────────────────────────────────────────────────────
add_executable(vdc vdc.cpp ${COMMON_SOURCES})

# ─── Test executable ─────────────────────────────────────────────────────────
add_executable(test_vor test_vor.cpp ${COMMON_SOURCES})

# ─── Link libraries ──────────────────────────────────────────────────────────
target_link_libraries(vdc
    PRIVATE
      CGAL::CGAL
      ZLIB::ZLIB
      ${TEEM_LIBRARY}
)

target_link_libraries(test_vor
    PRIVATE
      CGAL::CGAL
      ZLIB::ZLIB
      ${TEEM_LIBRARY}
)

# ─── Installation (optional) ─────────────────────────────────────────────────
install(TARGETS vdc test_vor
        RUNTIME DESTINATION bin
)

# ─── Usage info ──────────────────────────────────────────────────────────────
message(STATUS "To build:")
message(STATUS "  mkdir build && cd build")
message(STATUS "  cmake ..")
message(STATUS "  make")
message(STATUS "  ./vdc [your options]")
message(STATUS "  ./test_vor [options]")