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
#include <set>
#include <stack>
#include <iterator> // for std::back_inserter
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
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_3.h>            // Compute Halfspace Intersections in 3D
#include <CGAL/centroid.h>                            // Compute centroids of geometric objects.
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h> // Kernel for exact predicates and inexact constructions.
#include <CGAL/Delaunay_triangulation_3.h>            // Delaunay triangulation in 3D.
#include <CGAL/Triangulation_vertex_base_with_info_3.h> // Vertex base with additional user data.
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_with_circumcenter_3.h> // Cell base with circumcenter support.

//! @brief CGAL Kernel.
/*!
 * Provides geometric primitives like points, vectors, lines, and planes with
 * exact predicates and inexact constructions.
 */
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

struct VERTEX_INFO {
    bool is_dummy;                          // indicator of whether this delaunay vertex is added for bounding ( dummy )
    int voronoiCellIndex;                   // Index of the voronoi Cell that is dual to the vertex
    int index;                              // Index of the delaunay vertex in the triangulation
};

struct CELL_INFO {
    int dualVoronoiVertexIndex;             // Index of the voronoi vertex that is dual to this cell
    int index;                              // Index of this delaunay cell in the triangulation
};

//! @brief Vertex base for triangulations with additional information.
/*!
 * The boolean value in the vertex base can be used to store metadata
 * (e.g., whether a vertex is marked or visited).
 */
typedef CGAL::Triangulation_vertex_base_with_info_3<VERTEX_INFO, K> Vb;


//! @brief Cell base for Delaunay triangulations.
/*!
* This base provides support for calculating circumcenters of cells.
*/
typedef CGAL::Delaunay_triangulation_cell_base_with_circumcenter_3<K> Cb2;
 
typedef CGAL::Triangulation_cell_base_with_info_3<CELL_INFO,K,Cb2> Cb;

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


//
// Custom Output << Operators
//

//! @brief Output operator for Cell_handle
/*! 
 * This operator allows printing of Cell_handle objects to output streams.
 * It displays the cell's vertices, circumcenter, and additional info.
 *
 * @param os The output stream to write to
 * @param ch The Cell_handle to output
 * @return Reference to the output stream
 */
template <typename OSTREAM_TYPE>
OSTREAM_TYPE& operator<<(OSTREAM_TYPE& os, const Cell_handle& ch) {
    if (ch == nullptr) {
        os << "Cell_handle(nullptr)";
        return os;
    }

    // Print cell info
    os << "Cell_handle[" << &(*ch) << "]:\n";
    os << "  Info: { index: " << ch->info().index
       << ", dualVoronoiVertexIndex: " << ch->info().dualVoronoiVertexIndex << " }\n";
    
    // Print vertices
    os << "  Vertices:\n";
    for (int i = 0; i < 4; ++i) {
        Vertex_handle vh = ch->vertex(i);
        os << "    [" << i << "]: (" 
           << vh->point().x() << ", " 
           << vh->point().y() << ", " 
           << vh->point().z() << ")";
        
        if (vh->info().is_dummy) {
            os << " (dummy)";
        }
        os << " index: " << vh->info().index;
        os << "\n";
    }
    
    // Print circumcenter
    try {
        Point cc = ch->circumcenter();
        os << "  Circumcenter: (" 
           << cc.x() << ", " 
           << cc.y() << ", " 
           << cc.z() << ")\n";
    } catch (...) {
        os << "  Circumcenter: (calculation failed)\n";
    }
    
    return os;
}

//! @brief Output operator for Cell_circulator
/*!
 * This operator allows printing of Cell_circulator objects to output streams.
 * It displays the current cell the circulator points to.
 *
 * @param os The output stream to write to
 * @param cc The Cell_circulator to output
 * @return Reference to the output stream
 */
template <typename OSTREAM_TYPE>
OSTREAM_TYPE& operator<<(OSTREAM_TYPE& os, const Cell_circulator& cc) {
    if (cc == Cell_circulator()) {
        os << "Cell_circulator(default)";
        return os;
    }
    
    // Get the current cell and output it
    Cell_handle ch = cc;
    os << "Cell_circulator -> " << ch;
    
    return os;
}

//! @brief Output operator for Delaunay triangulation
/*!
 * This operator allows printing of Delaunay triangulation objects to output streams.
 * It displays basic statistics about the triangulation including:
 * - Number of vertices
 * - Number of finite cells
 * - Number of infinite cells
 * - Dimension of the triangulation
 *
 * @param os The output stream to write to
 * @param dt The Delaunay triangulation to output
 * @return Reference to the output stream
 */
template <typename OSTREAM_TYPE>
OSTREAM_TYPE& operator<<(OSTREAM_TYPE& os, const Delaunay& dt) {
    // Basic statistics
    os << "Delaunay Triangulation Statistics:\n"
       << "  Dimension: " << dt.dimension() << "\n"
       << "  Is valid: " << (dt.is_valid() ? "yes" : "no") << "\n\n";
    
    // Vertex statistics
    int total_vertices = dt.number_of_vertices();
    int dummy_vertices = 0;
    int regular_vertices = 0;
    
    for(auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
        if(vit->info().is_dummy) {
            dummy_vertices++;
        } else {
            regular_vertices++;
        }
    }
    
    os << "Vertex Statistics:\n"
       << "  Total vertices: " << total_vertices << "\n"
       << "  Regular vertices: " << regular_vertices << "\n"
       << "  Dummy vertices: " << dummy_vertices << "\n";
    
    // Cell statistics
    int finite_cells = dt.number_of_finite_cells();
    int infinite_cells = dt.number_of_cells() - finite_cells;
    
    os << "\nCell Statistics:\n"
       << "  Total cells: " << dt.number_of_cells() << "\n"
       << "  Finite cells: " << finite_cells << "\n"
       << "  Infinite cells: " << infinite_cells << "\n";
       
    // Output detailed information for finite cells
    os << "\nCell Details (limited to first 5 finite cells):\n";
    int cell_count = 0;
    for(auto cit = dt.finite_cells_begin(); cit != dt.finite_cells_end(); ++cit) {
        Cell_handle ch = cit;
        os << ch; 
        
        // Limit the number of cells to display to avoid excessive output
        if(++cell_count >= 5) {
            os << "... and " << (finite_cells - 5) << " more finite cells\n";
            break;
        }
        
        // Add a separator line
        os << "  -----------------------------------------\n";
    }
    
    // Vertex degree distribution (number of adjacent cells per vertex)
    std::map<int, int> degree_distribution;
    for(auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
        int degree = dt.degree(vit);
        degree_distribution[degree]++;
    }
    
    os << "\nVertex Degree Distribution:\n"
       << "  (degree: number of vertices)\n";
    for(const auto& pair : degree_distribution) {
        os << "  " << pair.first << " cells: " << pair.second << " vertices\n";
    }
    
    // Vertex detailed information
    os << "\nVertex Details:\n";
    int vertex_count = 0;
    for(auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
        os << "  [" << vit->info().index << "] ("
           << std::fixed << std::setprecision(3)
           << vit->point().x() << ", "
           << vit->point().y() << ", "
           << vit->point().z() << ")"
           << (vit->info().is_dummy ? " (dummy)" : "")
           << " voronoi_cell: " << vit->info().voronoiCellIndex
           << " degree: " << dt.degree(vit)
           << "\n";
        
        // Limit the number of vertices to display to avoid excessive output
        if(++vertex_count >= 10) {
            os << "  ... and " << (total_vertices - 10) << " more vertices\n";
            break;
        }
    }
    
    return os;
}

//! @brief Helper function to print all cells in a Cell_circulator
/*!
 * This function iterates through all cells in a Cell_circulator and prints them.
 * It's useful for debugging purposes to see all cells around a vertex.
 *
 * @param os The output stream to write to
 * @param start The Cell_circulator to start from
 * @param vertex_index The index of the vertex being circulated around (for information)
 */
template <typename OSTREAM_TYPE>
void print_cell_circuit(OSTREAM_TYPE& os, Cell_circulator start, int vertex_index) {
    if (start == Cell_circulator()) {
        os << "Empty Cell_circulator\n";
        return;
    }
    
    os << "Cells around vertex " << vertex_index << ":\n";
    
    Cell_circulator current = start;
    int count = 0;
    do {
        os << "Cell " << count << ":\n" << current << "\n";
        ++current;
        ++count;
    } while (current != start && count < 100); // Safety limit to prevent infinite loops
    
    if (count >= 100) {
        os << "Warning: Circuit may be infinite, stopped after 100 cells\n";
    }
    
    os << "Total: " << count << " cells\n";
}

#endif // VDC_TYPE_H
