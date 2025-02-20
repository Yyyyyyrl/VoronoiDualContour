//! @file vdc_func.h
//! @brief Header file for Voronoi Diagram and Isosurface computation functions.
#ifndef VDC_FUNC_H
#define VDC_FUNC_H

#include "vdc_type.h"
#include "vdc_utilities.h"
#include "vdc_io.h"
#include "vdc_commandline.h"
#include "vdc_voronoi.h"

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

//! @brief Computes the dual triangles for the final mesh in the single isovertex case.
/*!
 * This function calculates the Delaunay triangles dual to bipolar edges in
 * the Voronoi diagram for a single isovalue.
 * 
 * @param iso_surface Instance of IsoSurface
 * @param voronoi_edges Vector of Voronoi edges.
 * @param vertexValueMap Map of Voronoi vertices to scalar values.
 * @param bbox Bounding box of the computational domain.
 * @param delaunay_facet_to_voronoi_edge_map Map linking Delaunay facets to Voronoi edges.
 * @param dt Delaunay triangulation structure.
 * @param grid Scalar grid containing scalar values.
 * @param isovalue the isovalue used for computing.
 * @param point_index_map the map between the points in the delaunay triangulation to its index
 */
void computeDualTriangles(IsoSurface &iso_surface, std::vector<CGAL::Object> &voronoi_edges, std::map<Point, float> &vertexValueMap, CGAL::Epick::Iso_cuboid_3 &bbox, std::map<Object, std::vector<Facet>, ObjectComparator> &delaunay_facet_to_voronoi_edge_map, Delaunay &dt, ScalarGrid &grid, float isovalue, std::map<Point, int> &point_index_map);

//! @brief Computes the dual triangles for the final mesh in the multi-isovertex case.
/*!
 * This function calculates the Delaunay triangles dual to bipolar edges in
 * the Voronoi diagram for multiple isovalues.
 * 
 * @param voronoiDiagram The Voronoi diagram to compute from.
 * @param bbox Bounding box of the computational domain.
 * @param delaunay_facet_to_voronoi_edge_map Map linking Delaunay facets to Voronoi edges.
 * @param grid Scalar grid containing scalar values.
 * @param isovalue The isovalue for mesh computation.
 * @param iso_surface Instance of IsoSurface contains the isosurface vertices and faces
 */
void computeDualTrianglesMulti(
    VoronoiDiagram &voronoiDiagram,
    CGAL::Epick::Iso_cuboid_3 &bbox,
    std::map<CGAL::Object, std::vector<Facet>, ObjectComparator> &delaunay_facet_to_voronoi_edge_map,
    ScalarGrid &grid,
    float isovalue,
    IsoSurface &iso_surface);

//! @brief Computes isosurface vertices for the multi-isovertex case.
/*!
 * @param voronoiDiagram The Voronoi diagram to compute vertices for.
 * @param isovalue The isovalue to use for vertex computation.
 * @param iso_surface Instance of IsoSurface contains the isosurface vertices and faces
 */
void Compute_Isosurface_Vertices_Multi(VoronoiDiagram &voronoiDiagram, float isovalue, IsoSurface &iso_surface);

//! @brief Computes isosurface vertices for the single-isovertex case.
/*!
 * @param grid The scalar grid containing scalar values.
 * @param isovalue The iso value to use for calculation
 * @param iso_surface Instance of IsoSurface contains the isosurface vertices and faces
 * @param data_grid The grid containing input data
 * @param activeCubeCenters The list of center points of active cubes
 */
void Compute_Isosurface_Vertices_Single(ScalarGrid &grid, float isovalue, IsoSurface &iso_surface, Grid &data_grid, std::vector<Point> &activeCubeCenters);

/*
Setting up the delaunay triangulation
*/

//! @brief Constructs a Delaunay triangulation from a grid and grid facets.
/*!
 * This function constructs a 3D Delaunay triangulation using the grid's scalar values
 * and the facets of the active cubes.
 * 
 * @param dt The Deluanay Triangulation instance
 * @param grid The grid containing scalar values.
 * @param grid_facets The grid facets to use in constructing the triangulation.
 * @param vdc_param The VDC_PARAM instance that holds the commandline options the user input
 * @param activeCubeCenters The list of center points of active cubes
 * @param point_index_map the map between the points in the delaunay triangulation to its index (Used in Single Iso-V Case ONLY)
 */
void construct_delaunay_triangulation(Delaunay &dt, Grid &grid, const std::vector<std::vector<GRID_FACETS>> &grid_facets, VDC_PARAM &vdc_param, std::vector<Point> &activeCubeCenters, std::map<Point, int> &point_index_map);

//! @brief Adds dummy points from a facet for Voronoi diagram bounding.
/*!
 * This function adds additional dummy points to ensure that the Voronoi diagram is
 * bounded within the computational domain.
 * 
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
 * Generates the vertices of the Voronoi diagram based on the Delaunay triangulation.
 * 
 * @param voronoiDiagram The Voronoi diagram to construct vertices for.
 * @param dt The Deluanay Triangulation corresponding(dual) to the voronoi diagram
 */
void construct_voronoi_vertices(VoronoiDiagram &voronoiDiagram, Delaunay &dt);

//! @brief Constructs Voronoi cells from the Delaunay triangulation.
/*!
 * Populates the Voronoi diagram with polyhedral cells derived from the Delaunay triangulation.
 * 
 * @param voronoiDiagram The Voronoi diagram to populate with cells.
 * @param dt The Deluanay Triangulation corresponding(dual) to the voronoi diagram
 */
void construct_voronoi_cells(VoronoiDiagram &voronoiDiagram, Delaunay &dt);

//! @brief (in dev) Construct the voronoi cells routine that doesn't use Convex_Hull_3
/*!
 * Populates the Voronoi diagram with polyhedral cells derived from the Delaunay triangulation without using Convex_Hull_3
 * @param dt The Deluanay Triangulation corresponding(dual) to the voronoi diagram
 */
void construct_voronoi_cells_non_convex_hull(VoronoiDiagram &voronoiDiagram, Delaunay &dt);

//! @brief (in dev) Construct the voronoi cells routine as the intersection of halfspaces.
/*!
 *
 * @param dt The Deluanay Triangulation corresponding(dual) to the voronoi diagram
 */
void construct_voronoi_cells_halfspace(VoronoiDiagram &voronoiDiagram, Delaunay &dt);

//! @brief Computes Voronoi vertex values using scalar grid interpolation.
/*!
 * Interpolates scalar values from the grid to each vertex of the Voronoi diagram.
 * 
 * @param voronoiDiagram The Voronoi diagram to compute values for.
 * @param grid The scalar grid containing data.
 * @param VertexValueMap The map between vertices in the voronoi diagram to its scalar values
 */
void compute_voronoi_values(VoronoiDiagram &voronoiDiagram, ScalarGrid &grid, std::map<Point, float> &vertexValueMap);


//! @brief Constructs Voronoi edges from Delaunay facets.
/*!
 * Derives the edges of the Voronoi diagram by processing the facets of the Delaunay triangulation.
 * 
 * @param voronoiDiagram The Voronoi diagram to populate with edges.
 * @param delaunay_facet_to_voronoi_edge_map Map linking Delaunay facets to Voronoi edges.
 * @param dt The Deluanay Triangulation corresponding(dual) to the voronoi diagram
 */
void construct_voronoi_edges(
    VoronoiDiagram &voronoiDiagram,
    std::map<CGAL::Object, std::vector<Facet>, ObjectComparator> &delaunay_facet_to_voronoi_edge_map,
    Delaunay &dt);


//! @brief Constructs the VoronoiCelledges in the VoronoiDiagram and link them
/*!
 *
 * @param voronoiDiagram The Voronoi diagram to populate with edges.
 * @param delaunay_facet_to_voronoi_edge_map Map linking Delaunay facets to Voronoi edges.
 * @param bbox The bounding box used for clipping Ray and Line voronoi Edges
 * @param dt The Deluanay Triangulation corresponding(dual) to the voronoi diagram
 */
void construct_voronoi_cell_edges(VoronoiDiagram &voronoiDiagram,
    std::map<CGAL::Object, std::vector<Facet>, ObjectComparator> &voronoi_edge_to_delaunay_facet_map,
    CGAL::Epick::Iso_cuboid_3 &bbox,
    Delaunay &dt);
    
/* Output Handling */

//! @brief Handles output mesh generation.
/*!
 * Generates the final output mesh based on the Voronoi diagram and handles file output.
 * 
 * @param retFlag Reference to a flag indicating success or failure.
 * @param vd The Voronoi diagram containing mesh data.
 * @param vdc_param The VDC_PARAM instance that contains the user input options
 * @param iso_surface The Instance of IsoSurface containing the vertices and faces of the isosurface
 * @param point_index_map the map between the points in the delaunay triangulation to its index (Used in Single Iso-V Case ONLY)
 * @return An integer representing the exit status.
 */
int handle_output_mesh(bool &retFlag, VoronoiDiagram &vd, VDC_PARAM &vdc_param, IsoSurface &iso_surface, std::map<Point, int> &point_index_map);


//! @brief Wrap up function of all process of building the voronoi diagram from the Delaunay Triangulation
/*!
 *
 */
void construct_voronoi_diagram(VoronoiDiagram &vd, VDC_PARAM &vdc_param, std::map<CGAL::Object, std::vector<Facet>, ObjectComparator> &voronoi_edge_to_delaunay_facet_map, ScalarGrid &grid, std::map<Point, float> &vertexValueMap, CGAL::Epick::Iso_cuboid_3 &bbox, Delaunay &dt);

//! @brief Wrap up function for all process of building the isosurface (vertices and faces) from the Voronoi Diagram / Delaunay Triangulation
/*!
 *
 */
void construct_iso_surface(Delaunay &dt, VoronoiDiagram &vd, VDC_PARAM &vdc_param, IsoSurface &iso_surface, ScalarGrid &grid, Grid &data_grid, std::vector<Point> &activeCubeCenters,std::map<CGAL::Object, std::vector<Facet>, ObjectComparator> &voronoi_edge_to_delaunay_facet_map, std::map<Point, float> &vertexValueMap, CGAL::Epick::Iso_cuboid_3 &bbox, std::map<Point, int> &point_index_map);
#endif