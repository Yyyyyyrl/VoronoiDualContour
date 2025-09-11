# Voronoi-Based Dual Contouring (VDC) - 3D Isosurface Extraction

## Overview

VDC is a C++ implementation of a Voronoi diagram-based dual contouring method for high-quality isosurface extraction from 3D scalar fields. The project leverages CGAL for computational geometry operations


## Installation

### Prerequisites
- CMake (3.12+)
- CGAL (5.x) with Core component
- Zlib
- Teem (for NRRD/NHDR IO)
- Python 3 (for visualization scripts only)

### Building from Source
- Ensure dependencies are installed and discoverable (e.g., via your package manager).
- Teem does not provide a CMake config; if it is not in a default prefix, pass `-DTEEM_ROOT=/path/to/teem-prefix` where `include/` and `lib/` live.

```bash
git clone https://github.com/Yyyyyyrl/voronoi-dual-contouring.git
cd voronoi-dual-contouring
mkdir build && cd build
cmake -DTEEM_ROOT=/absolute/path/to/teem ..   # optional if Teem in default path
make -j
```

## Usage

### Basic Command Line Syntax
```bash
./vdc [OPTIONS] <isovalue> <(nhdr/nrrd) raw data file path>
```

### Options
- -o <file>: output basename (default derived from input)
- -off: write mesh in OFF format (default)
- -ply: write mesh in PLY format
- -out_csv <file>: dump Voronoi diagram to CSV
- -sep_isov: pick a subset of non-adjacent active cubes before triangulation
- -supersample <factor>: supersample the input grid by factor
- -multi_isov: enable multi iso-vertices mode
- -single_isov: use single iso-vertex mode (default)
- -conv_H: use CGAL convex hull for Voronoi cell construction
- -mod_cyc: modify-cycles pass after initial cycles (advanced)
- --help: print help

Advanced/debug options (subject to change):
- -bound_cells: add bounding cells around the domain
- --debug: enable verbose debug logs ([DEBUG]/[ISO]/[ISO-MATCH])

### Examples
- Basic run (OFF output):
  `./vdc 0.0 ./data/sphere-32.nrrd`
- PLY output with supersampling x2:
  `./vdc -ply -supersample 2 0.0 ./data/sphere-32.nrrd`


## File Structure

```
.
├── CMakeLists.txt            # CMake build configuration
├── README.md                 # Project documentation
├── cmake/                    # CMake helpers (if any)
├── include/                  # Header files
│   ├── core/                 # Core types, utilities, CLI
│   │   ├── vdc.h
│   │   ├── vdc_type.h
│   │   ├── vdc_utilities.h
│   │   ├── vdc_debug.h
│   │   └── vdc_commandline.h
│   ├── grid/                 # Grid data structures
│   │   └── vdc_grid.h
│   ├── io/                   # IO helpers (NRRD)
│   │   └── vdc_io.h
│   ├── delaunay/             # Delaunay types/utilities
│   │   └── vdc_delaunay.h
│   ├── voronoi/              # Voronoi structures and ops
│   │   └── vdc_voronoi.h
│   ├── algo/                 # High-level algorithm declarations
│   │   └── vdc_func.h
│   └── test/                 # Test helpers
│       └── test_vor.h
├── src/                      # Source files
│   ├── app/                  # Application entry points
│   │   └── vdc.cpp
│   ├── core/                 # Core utilities and CLI
│   │   ├── vdc_commandline.cpp
│   │   ├── vdc_debug.cpp
│   │   └── vdc_utilities.cpp
│   ├── grid/
│   │   └── vdc_grid.cpp
│   ├── io/
│   │   └── vdc_io.cpp
│   ├── algo/                 # Algorithm implementations
│   │   ├── vdc_func_delaunay.cpp
│   │   ├── vdc_func_helper.cpp
│   │   ├── vdc_func_isosurface.cpp
│   │   ├── vdc_func_modcyc.cpp
│   │   └── vdc_func_voronoi.cpp
│   ├── voronoi/              # Voronoi core and utilities
│   │   ├── vdc_voronoi_collapse.cpp
│   │   ├── vdc_voronoi_core.cpp
│   │   └── vdc_voronoi_matching.cpp
│   └── tests/                # Test drivers
│       ├── test_vor.cpp
│       └── test_modcyc.cpp
├── tools/
│   ├── plot_modcyc.py        # Visualize modify-cycles test case
│   └── compExec.py           # Compare execution output of different configurations
├── data/                     # Sample volumes (NRRD/NHDR/RAW) and the outputs
```

## API Documentation

The project provides comprehensive API documentation through Doxygen. To generate:
```bash
doxygen Doxyfile
```
