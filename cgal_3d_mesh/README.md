# CGAL 3D Mesh Generator 

## Overview

The `cgal_mesher` program implements CGAL's modern **Mesh_3 package** (3D Mesh Generation) to extract isosurface meshes from NRRD volumetric data.

## Building the Program

### Prerequisites

- **CGAL** (with Core and ImageIO components)
- **Teem** (for NRRD I/O)
- **CMake** 3.12+
- **C++17** compiler
### Build Instructions

```bash
cd cgal_comparison
mkdir build
cd build
cmake ..
make
```



## Usage

### Basic Syntax

```bash
./cgal_mesher [OPTIONS] <isovalue> <nrrd_file>
```

### Arguments

- **`isovalue`**: The scalar value defining the isosurface (e.g., 128 for mid-range in 8-bit data)
- **`nrrd_file`**: Path to input NRRD volumetric data file

### Options

| Option | Description | Default |
|--------|-------------|---------|
| `-o <filename>` | Output mesh filename | `cgal_output.off` |
| `-facet_angle <degrees>` | Minimum facet angle (≤30° for guaranteed termination) | `20.0` |
| `-facet_size <value>` | Maximum facet size (Delaunay ball radius) | Auto (5% of max extent) |
| `-facet_distance <value>` | Maximum distance from circumcenter to surface | Auto (1% of max extent) |
| `-cell_ratio <value>` | Cell radius-edge ratio (≥2 for termination) | `3.0` |
| `-cell_size <value>` | Maximum cell size (tetrahedron circumradius) | Auto (8% of max extent) |
| `-surface_only` | Generate surface mesh only (no volume tetrahedralization) | Off (generates volume) |
| `-h, --help` | Display help message | - |

## Quality Parameters

### Facet Criteria (Surface Quality)

#### Facet Angle
- Controls the **minimum angle** in surface mesh triangles
- Must be ≤30° for theoretical termination guarantee
- Default: 25°, Range: typically 20-30°

#### Facet Size
- Controls the **maximum size** of mesh facets
- Upper bound on the radius of surface Delaunay balls
- Smaller values → finer mesh with more triangles
- Measured in world units (matches NRRD spacing)

#### Facet Distance
- Controls **surface approximation accuracy**
- Maximum distance between facet circumcenter and the actual surface
- Smaller values → closer fit to true isosurface, more triangles
- Measured in world units

### Cell Criteria (Volume Quality, if generating volume mesh)

#### Cell Radius-Edge Ratio
- Controls **tetrahedron shape quality**
- Ratio of circumradius to shortest edge length
- Must be >2 for termination guarantee
- Default: 3.0

#### Cell Size
- Upper bound on tetrahedron **circumradius**

### Surface-Only Mode
- Use `-surface_only` flag to skip volume mesh generation
- Faster execution, outputs only the boundary surface


## References

- [CGAL 3D Mesh Generation Documentation](https://doc.cgal.org/latest/Mesh_3/index.html)
- [CGAL Mesh_3 User Manual](https://doc.cgal.org/latest/Mesh_3/index.html#Chapter_3D_Mesh_Generation)
- [Teem NRRD Format](http://teem.sourceforge.net/nrrd/)
