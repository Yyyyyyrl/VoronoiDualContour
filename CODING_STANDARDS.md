# Coding Standards for DelaunayContour

This document outlines the coding conventions and standards used in the DelaunayContour project.

## Naming Conventions

### Types (Classes, Structs, Enums)
- **Format**: PascalCase
- **Examples**:
  - `VoronoiVertex`, `VoronoiEdge`, `VoronoiCell`
  - `DelaunayTriangle`, `UnifiedGrid`
  - `VertexInfo`, `CellInfo`, `DelaunayFacetInfo`
  - `GridFacets`, `VdcParam`

### Functions and Methods
- **Format**: snake_case
- **Examples**:
  - `construct_voronoi_vertices()`
  - `compute_dual_triangles()`
  - `is_bipolar()`, `is_cube_active()`
  - `write_off_single()`, `write_ply_multi()`
  - `read_points_from_file()`

### Member Variables
- **Format**: camelCase (to distinguish from function names)
- **Examples**:
  - `cellIndex`, `edgeIndex`
  - `verticesIndices`, `facetIndices`
  - `delaunayVertex`, `voronoiCellIndex`
  - `cellEdgeIndices`, `voronoiEdgeIndices`

### Local Variables and Parameters
- **Format**: snake_case
- **Examples**:
  - `cell_index`, `vertex_count`
  - `iso_surface`, `voronoi_diagram`
  - `grid_facets`, `active_cubes`

### Constants and Macros
- **Format**: SCREAMING_SNAKE_CASE
- **Examples**:
  - `DIM3`, `MAX_SIZE`
  - `ISO_DBG_ENABLED`

### Enum Values
- **Format**: SCREAMING_SNAKE_CASE
- **Examples**:
  - `SEP_POS`, `SEP_NEG`
  - `UNCONSTRAINED_MATCH`, `UNDEFINED_MATCH_TYPE`

## Documentation Style

### Doxygen Comments
All public functions, classes, and structs should have Doxygen-style documentation.

#### File Headers
```cpp
//! @file filename.h
//! @brief Brief description of the file's purpose
```

#### Function Documentation
```cpp
//! @brief Brief one-line description
/*!
 * Detailed description (optional, for complex functions)
 *
 * @param param1 Description of param1
 * @param param2 Description of param2
 * @return Description of return value
 */
return_type function_name(type param1, type param2);
```

#### Struct/Class Documentation
```cpp
//! @brief Brief description
/*!
 * Detailed description of the struct/class
 */
struct MyStruct
{
    int memberVar;  //!< Brief description of member
};
```

## Code Organization

### Header Files
- Use include guards: `#ifndef HEADER_NAME_H` / `#define HEADER_NAME_H` / `#endif`
- Order includes: System headers → Third-party libraries → Project headers
- Keep header files minimal - declare only what's necessary

### Source Files
- One class/major functionality per file
- Keep functions focused and concise
- Use anonymous namespaces for file-local helper functions

## Best Practices

### General
- Prefer const correctness
- Use `auto` for complex iterator types, but be explicit for readability where appropriate
- Avoid magic numbers - use named constants
- Initialize variables at declaration when possible

### Comments
- Use `//` for single-line comments
- Use `/* */` for multi-line comments
- Write self-documenting code - comment "why", not "what"
- Keep comments up to date with code changes

### Error Handling
- Check return values and file operations
- Use `std::cerr` for error messages
- Provide meaningful error messages

### Formatting
- Indentation: 4 spaces (no tabs)
- Maximum line length: 120 characters (flexible)
- Braces: Opening brace on same line for functions/control structures
- Space after keywords (if, for, while), not after function names

### Example
```cpp
//! @brief Checks if two values are bipolar with respect to an isovalue
/*!
 * @param val1 First scalar value
 * @param val2 Second scalar value  
 * @param isovalue The threshold value (default is 0)
 * @return true if values are bipolar, false otherwise
 */
bool is_bipolar(float val1, float val2, float isovalue = 0)
{
    return ((val1 < isovalue) && (val2 >= isovalue)) ||
           ((val1 >= isovalue) && (val2 < isovalue));
}
```

## Git Commit Messages
- Use present tense ("Add feature" not "Added feature")
- Start with a capital letter
- Keep first line under 72 characters
- Provide detailed description in body for complex changes
- Reference issue numbers when applicable

## Building and Testing
- Always build and test before committing
- Run all test suites: `vdc`, `test_vor`, `test_modcyc`
- Ensure no compiler warnings with `-Wall -Wextra`

## Dependencies
- CGAL for computational geometry
- Teem for NRRD file handling
- Boost (headers only)
- Standard C++17 or later

---

*Last updated: 2025-11-13*
*This document reflects the standardization completed in the November 2025 codebase cleanup.*
