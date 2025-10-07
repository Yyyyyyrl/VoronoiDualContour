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
- `-o {output_filename}`: Specify output filename (default: derived from input filename)
- `-off`: Generate output in .off format (default)
- `-ply`: Generate output in .ply format
- `-out_csv {output_csv_name}`: Write the Voronoi diagram to a CSV file
- `-sep_isov_1`: Separation method I: Greedy cube-level (26-connectivity)
- `-sep_isov_3`: Separation method III: 3×3×3 subgrid-based separation
- `-supersample {factor}`: Supersample the input data by the given factor
- `-collapse_eps {eps}`: Set absolute collapse threshold in world units (default: 1% of grid spacing)
- `-multi_isov`: Use multi iso-vertices mode
- `-single_isov`: Use single iso-vertices mode (default)
- `-conv_H`: Use the Convex_Hull_3 from CGAL in voronoi cell construction
- `-mod_cyc`: After initial cycles, try facet rematching and recompute cycles
- `--summary_stats`: Print summary statistics after the run
- `--debug`: Enable debug logging ([DEBUG]/[ISO]/[ISO-MATCH]/[CYC-MOD])
- `--help`: Print help message

Advanced/debug options (subject to change):
- `-bound_cells`: Add bounding cells around the domain
- `--test_vor`: Flag for testing the Voronoi diagram construction

### Examples
- Basic run (OFF output):
  ```bash
  ./vdc 0.0 ./data/sphere-32.nrrd
  ```
- PLY output with supersampling x2:
  ```bash
  ./vdc -ply -supersample 2 0.0 ./data/sphere-32.nrrd
  ```
- Multi iso-vertices with separation method I:
  ```bash
  ./vdc -multi_isov -sep_isov_1 0.0 ./data/sphere-32.nrrd
  ```
- Export Voronoi diagram to CSV with custom output name:
  ```bash
  ./vdc -o sphere_output -out_csv voronoi_data.csv 0.0 ./data/sphere-32.nrrd
  ```
- Using modify-cycles with summary statistics:
  ```bash
  ./vdc -mod_cyc --summary_stats 0.0 ./data/sphere-32.nrrd
  ```


## File Structure

```
.
├── CMakeLists.txt            # CMake build configuration
├── README.md                 # Project documentation
├── cmake/                    # CMake helpers
├── include/                  # Header files
│   ├── core/                 # Core types, utilities, CLI
│   │   ├── vdc.h
│   │   ├── vdc_type.h
│   │   ├── vdc_utilities.h
│   │   ├── vdc_debug.h
│   │   └── vdc_commandline.h
│   ├── grid/                 # Grid data structures and separation methods
│   │   ├── vdc_grid.h
│   │   └── vdc_sep_isov.h
│   ├── io/                   # IO helpers (NRRD)
│   │   └── vdc_io.h
│   ├── delaunay/             # Delaunay types/utilities
│   │   └── vdc_delaunay.h
│   ├── voronoi/              # Voronoi structures and operations
│   │   └── vdc_voronoi.h
│   ├── algo/                 # High-level algorithm declarations
│   │   └── vdc_func.h
│   └── test/                 # Test helpers
│       └── test_vor.h
├── src/                      # Source files
│   ├── app/                  # Application entry point
│   │   └── vdc.cpp
│   ├── core/                 # Core utilities and CLI
│   │   ├── vdc_commandline.cpp
│   │   ├── vdc_debug.cpp
│   │   └── vdc_utilities.cpp
│   ├── grid/                 # Grid operations and separation methods
│   │   ├── vdc_grid.cpp
│   │   ├── vdc_sep_isov.cpp
│   │   └── vdc_sep_isov_debug.cpp
│   ├── io/                   # IO implementations
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
├── tools/                    # Utility scripts
│   ├── plot_modcyc.py        # Visualize modify-cycles test case
│   ├── plot_cells_and_edge.py # Visualize cells and edges
│   ├── run_vdc_batch.py      # Generate output for different isovalues
│   ├── compExec.py           # Compare execution output of different
```

## API Documentation

The project provides comprehensive API documentation through Doxygen. To generate:
```bash
doxygen Doxyfile
```
