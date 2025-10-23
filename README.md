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
- `-sep_isov_3_wide`: Testing variant of method III with 5×5×5 clearance in the 3× subgrid
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
├── include/                  # Header files (14 total)
│   ├── core/                 # Core infrastructure (7 files)
│   │   ├── vdc.h             # Main umbrella header
│   │   ├── vdc_type.h        # Type definitions and CGAL types
│   │   ├── vdc_utilities.h   # Utility functions and helpers
│   │   ├── vdc_debug.h       # Debug output and logging
│   │   ├── vdc_commandline.h # CLI parsing
│   │   ├── vdc_timing.h      # Performance timing
│   │   └── vdc_stats.h       # Statistics collection and reporting
│   ├── processing/           # All algorithm headers (5 files)
│   │   ├── vdc_func.h        # High-level algorithm orchestration
│   │   ├── vdc_delaunay.h    # Delaunay triangulation structures
│   │   ├── vdc_voronoi.h     # Voronoi diagram structures and operations
│   │   ├── vdc_grid.h        # Grid data structures
│   │   └── vdc_sep_isov.h    # Isosurface vertex separation methods
│   ├── vdc_io.h              # I/O operations (NRRD, OFF, PLY)
│   └── test_vor.h            # Test utilities
├── src/                      # Source files (20 total)
│   ├── core/                 # Core implementations (5 files)
│   │   ├── vdc_commandline.cpp
│   │   ├── vdc_debug.cpp
│   │   ├── vdc_utilities.cpp
│   │   ├── vdc_timing.cpp
│   │   └── vdc_stats.cpp     # Statistics implementation
│   ├── processing/           # All algorithmic implementations (11 files)
│   │   ├── vdc_func_delaunay.cpp      # Delaunay triangulation construction
│   │   ├── vdc_func_voronoi.cpp       # Voronoi diagram construction
│   │   ├── vdc_func_isosurface.cpp    # Isosurface extraction
│   │   ├── vdc_func_modcyc.cpp        # Modify-cycles algorithm
│   │   ├── vdc_func_helper.cpp        # Algorithm helper functions
│   │   ├── vdc_voronoi_core.cpp       # Core Voronoi operations
│   │   ├── vdc_voronoi_collapse.cpp   # Edge collapse operations
│   │   ├── vdc_voronoi_matching.cpp   # Bipolar edge matching
│   │   ├── vdc_grid.cpp               # Grid operations
│   │   ├── vdc_sep_isov.cpp           # Separation implementations
│   │   └── vdc_sep_isov_debug.cpp     # Separation debug utilities
│   ├── testing/              # Test drivers (2 files)
│   │   ├── test_vor.cpp      # Voronoi test driver
│   │   └── test_modcyc.cpp   # Modify-cycles test driver
│   ├── vdc_main.cpp          # Main application entry point
│   └── vdc_io.cpp            # I/O implementations
├── tools/                    # Utility scripts
│   ├── plot_modcyc.py        # Visualize modify-cycles test case
│   ├── plot_cells_and_edge.py # Visualize cells and edges
│   ├── run_vdc_batch.py      # Generate output for different isovalues
│   └── compExec.py           # Compare execution outputs
```

**Directory Organization Philosophy:**
- **`core/`**: Infrastructure code (types, utilities, I/O, CLI, debugging, timing, stats)
- **`processing/`**: All computational geometry algorithms (Delaunay, Voronoi, grid, isosurface extraction)
- **`testing/`**: Test drivers and utilities
- **Root level**: Main application and I/O implementations

## API Documentation

The project provides comprehensive API documentation through Doxygen. To generate:
```bash
doxygen Doxyfile
```
