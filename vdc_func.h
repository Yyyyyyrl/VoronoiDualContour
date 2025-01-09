//! @file vdc_func.h
//! @brief Header file for Voronoi Diagram and Isosurface computation functions.
#ifndef VDC_FUNC_H
#define VDC_FUNC_H

#include "vdc_type.h"
#include "vdc_utilities.h"
#include "vdc_io.h"
#include "vdc_globalvar.h"

/*Struct*/

//! @brief Structure for comparing points for approximate equality.
/*!
 * This structure defines a custom comparator for points, allowing comparison
 * with a small epsilon tolerance.
 */
struct PointApproxEqual
{

    //! @brief Compares two points for approximate equality.
    /*!
     * @param p1 First point.
     * @param p2 Second point.
     * @return `true` if the points are approximately equal, `false` otherwise.
     */
    bool operator()(const Point &p1, const Point &p2) const
    {
        const double EPSILON = 1e-6;
        return (CGAL::squared_distance(p1, p2) < EPSILON * EPSILON);
    }
};

/*
Functions for both single/multi isov that computes the vertices and faces of the final mesh
*/

//! @brief Computes the dual triangles for the final mesh in single iso vertex case.
/*!
 * @param voronoi_edges Vector of Voronoi edges.
 * @param vertexValueMap Map of Voronoi vertices to scalar values.
 * @param bbox Bounding box of the domain.
 * @param delaunay_facet_to_voronoi_edge_map Map linking Delaunay facets to Voronoi edges.
 * @param dt Delaunay triangulation.
 * @param grid Scalar grid containing scalar values.
 * @return A vector of Delaunay triangles.
 */
std::vector<DelaunayTriangle> computeDualTriangles(std::vector<CGAL::Object> &voronoi_edges, std::map<Point, float> &vertexValueMap, CGAL::Epick::Iso_cuboid_3 &bbox, std::map<Object, std::vector<Facet>, ObjectComparator> &delaunay_facet_to_voronoi_edge_map, Delaunay &dt, ScalarGrid &grid);

//! @brief Computes the dual triangles for the final mesh in the multi iso vertex case.
/*!
 * @param voronoiDiagram The Voronoi diagram to compute from.
 * @param bbox Bounding box of the domain.
 * @param delaunay_facet_to_voronoi_edge_map Map linking Delaunay facets to Voronoi edges.
 * @param grid Scalar grid containing scalar values.
 * @param isovalue The isovalue for mesh computation.
 */
void computeDualTrianglesMulti(
    VoronoiDiagram &voronoiDiagram,
    CGAL::Epick::Iso_cuboid_3 &bbox,
    std::map<CGAL::Object, std::vector<Facet>, ObjectComparator> &delaunay_facet_to_voronoi_edge_map,
    ScalarGrid &grid,
    float isovalue);

//! @brief Computes isosurface vertices for the multi-isovertex case.
/*!
 * @param voronoiDiagram The Voronoi diagram to compute vertices for.
 * @param isovalue The isovalue to use for vertex computation.
 */
void Compute_Isosurface_Vertices_Multi(VoronoiDiagram &voronoiDiagram, float isovalue);

//! @brief Computes isosurface vertices for the single-isovertex case.
/*!
 * @param voronoiDiagram The Voronoi diagram to compute vertices for.
 * @param grid The scalar grid containing scalar values.
 */
void Compute_Isosurface_Vertices_Single(VoronoiDiagram &voronoiDiagram, ScalarGrid &grid);

/*
Setting up the delaunay triangulation
*/

//! @brief Constructs a Delaunay triangulation from a grid and grid facets.
/*!
 * @param grid The grid containing scalar values.
 * @param grid_facets The grid facets to use in constructing the triangulation.
 */
void construct_delaunay_triangulation(Grid &grid, const std::vector<std::vector<GRID_FACETS>> &grid_facets);

//! @brief Adds dummy points from a facet for Voronoi diagram bounding.
/*!
 * @param facet The facet to extract dummy points from.
 * @param data_grid The grid containing data for point computation.
 * @return A vector of points added from the facet.
 */
std::vector<Point> add_dummy_from_facet(const GRID_FACETS &facet, const Grid &grid);


/*
Functions that relates to building the voronoi diagram
*/

//! @brief Constructs Voronoi vertices for the given diagram.
/*!
 * @param voronoiDiagram The Voronoi diagram to construct vertices for.
 */
void construct_voronoi_vertices(VoronoiDiagram &voronoiDiagram);

//! @brief Constructs Voronoi cells from the Delaunay triangulation.
/*!
 * @param voronoiDiagram The Voronoi diagram to populate with cells.
 */
void construct_voronoi_cells(VoronoiDiagram &voronoiDiagram);

//! @brief Computes Voronoi vertex values using scalar grid interpolation.
/*!
 * @param voronoiDiagram The Voronoi diagram to compute values for.
 * @param grid The scalar grid containing data.
 */
void compute_voronoi_values(VoronoiDiagram &voronoiDiagram, ScalarGrid &grid);


//! @brief Constructs Voronoi edges from Delaunay facets.
/*!
 * @param voronoiDiagram The Voronoi diagram to populate with edges.
 * @param delaunay_facet_to_voronoi_edge_map Map linking Delaunay facets to Voronoi edges.
 */
void construct_voronoi_edges(
    VoronoiDiagram &voronoiDiagram,
    std::map<CGAL::Object, std::vector<Facet>, ObjectComparator> &delaunay_facet_to_voronoi_edge_map);

//! @brief Handles output mesh generation.
/*!
 * @param retFlag Reference to a flag indicating success or failure.
 * @param vd The Voronoi diagram containing mesh data.
 * @return An integer representing the exit status.
 */
int handle_output_mesh(bool &retFlag, VoronoiDiagram &vd);

#endif