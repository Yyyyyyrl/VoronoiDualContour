//! @file vdc_type.h
//! @brief Type definitions and includes for Voronoi Diagram and Delaunay Triangulation operations.

#ifndef VDC_TYPE_H
#define VDC_TYPE_H

// Standard library headers for various utilities.
#include <iostream>    // Input and output stream.
#include <fstream>     // File stream for reading and writing files.
#include <vector>      // STL vector container.
#include <array>       // STL array container.
#include <cmath>       // Math functions like abs(), pow(), etc.
#include <limits>      // To work with numeric limits of data types.
#include <algorithm>   // Common algorithms like sort(), min(), max().
#include <cstddef>     // Definitions for size_t, ptrdiff_t, etc.
#include <cstring>     // C-style string operations.
#include <teem/nrrd.h> // Teem library for handling NRRD data files.
#include <iomanip>     // For formatted output (e.g., precision control).

// CGAL headers for computational geometry operations.
#include <CGAL/bounding_box.h>                         // Compute bounding boxes for geometric objects.
#include <CGAL/intersections.h>                       // Perform geometric intersection tests.
#include <CGAL/Vector_3.h>                            // Represents a vector in 3D space.
#include <CGAL/Polyhedron_3.h>                        // Represents 3D polyhedral surfaces.
#include <CGAL/convex_hull_3.h>                       // Compute convex hulls in 3D.
#include <CGAL/centroid.h>                            // Compute centroids of geometric objects.
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h> // Kernel for exact predicates and inexact constructions.
#include <CGAL/Delaunay_triangulation_3.h>            // Delaunay triangulation in 3D.
#include <CGAL/Triangulation_vertex_base_with_info_3.h> // Vertex base with additional user data.
#include <CGAL/Delaunay_triangulation_cell_base_with_circumcenter_3.h> // Cell base with circumcenter support.

//! @brief CGAL Kernel.
/*!
 * Provides geometric primitives like points, vectors, lines, and planes with
 * exact predicates and inexact constructions.
 */
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

//! @brief Vertex base for triangulations with additional information.
/*!
 * The boolean value in the vertex base can be used to store metadata
 * (e.g., whether a vertex is marked or visited).
 */
typedef CGAL::Triangulation_vertex_base_with_info_3<bool, K> Vb;

//! @brief Cell base for Delaunay triangulations.
/*!
 * This base provides support for calculating circumcenters of cells.
 */
typedef CGAL::Delaunay_triangulation_cell_base_3<K> Cb;

//! @brief Data structure for triangulations.
/*!
 * Combines the vertex and cell bases into a complete triangulation data structure.
 */
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;

//! @brief 3D Delaunay triangulation.
/*!
 * A Delaunay triangulation is a geometric structure that divides 3D space
 * into tetrahedra such that no vertex lies inside the circumsphere of any tetrahedron.
 */
typedef CGAL::Delaunay_triangulation_3<K, Tds> Delaunay;

//! @brief 3D polyhedral surface representation.
/*!
 * A polyhedron is used to represent Voronoi cells and other geometric structures.
 */
typedef CGAL::Polyhedron_3<K> Polyhedron_3;

//! @brief 3D point type used in Delaunay triangulations.
/*!
 * Points are the fundamental objects used for constructing the triangulation.
 */
typedef Delaunay::Point Point;

//! @brief Handle to a cell in the Delaunay triangulation.
/*!
 * A cell handle allows access to the tetrahedra in the triangulation.
 */
typedef Delaunay::Cell_handle Cell_handle;

//! @brief Handle to a vertex in the Delaunay triangulation.
/*!
 * A vertex handle allows access to a specific vertex in the triangulation.
 */
typedef Delaunay::Vertex_handle Vertex_handle;

//! @brief Circulator for iterating through cells around a vertex.
/*!
 * This circulator iterates through all cells sharing a specific vertex.
 */
typedef Delaunay::Cell_circulator Cell_circulator;

//! @brief CGAL object wrapper for geometric primitives.
/*!
 * Used to store generic geometric objects like points, segments, or lines.
 */
typedef CGAL::Object Object;

//! @brief Represents an edge in the Delaunay triangulation.
typedef Delaunay::Edge Edge;

//! @brief Represents a facet (triangle face) in the Delaunay triangulation.
typedef Delaunay::Facet Facet;

//! @brief Represents a triangle in 3D space.
typedef Delaunay::Triangle Triangle;

//! @brief Represents a line segment in 3D space.
typedef K::Segment_3 Segment3;

//! @brief Represents a ray (semi-infinite line) in 3D space.
typedef K::Ray_3 Ray3;

//! @brief Represents a line in 3D space.
typedef K::Line_3 Line3;

//! @brief Represents a point in 3D space.
typedef K::Point_3 Point3;

//! @brief Represents a vector in 3D space.
/*!
 * Vectors are used for translations and other vector operations.
 */
typedef K::Vector_3 Vector3;

//! @brief Represents a plane in 3D space.
/*!
 * A plane is defined by a point and a normal vector or by three non-collinear points.
 */
typedef K::Plane_3 Plane_3;

#endif // VDC_TYPE_H
