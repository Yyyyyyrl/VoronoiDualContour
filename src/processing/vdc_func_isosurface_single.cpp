//! @file vdc_func_isosurface_single.cpp
//! @brief Single-isovertex isosurface construction functions.
/*!
 * This file contains all functions specific to the single-isovertex mode
 * where each active cube generates exactly one isovertex. These functions
 * are simpler than their multi-isovertex counterparts and do not require
 * cycle detection or conflict resolution.
 */

#include "processing/vdc_func.h"
#include "processing/vdc_grid.h"
#include "core/vdc_utilities.h"
#include "core/vdc_timing.h"

//! @brief Processes a segment edge for dual triangle computation.
/*!
 * Checks if the segment is bipolar and generates triangles for each associated
 * Delaunay facet.
 *
 * @param edge The CGAL object representing the edge.
 * @param vd The Voronoi Diagram.
 * @param isovalue The isovalue for bipolarity check.
 * @param dt The Delaunay triangulation.
 * @param dualTriangles Vector to store generated triangles.
 */
static void process_segment_edge(
    VoronoiEdge &edge,
    VoronoiDiagram &vd,
    float isovalue,
    Delaunay &dt,
    std::vector<DelaunayTriangle> &dualTriangles)
{
    Point v1 = vd.vertices[edge.vertex1].coord;
    Point v2 = vd.vertices[edge.vertex2].coord;
    float v1_val = vd.vertices[edge.vertex1].value;
    float v2_val = vd.vertices[edge.vertex2].value;

    if (is_bipolar(v1_val, v2_val, isovalue))
    {
        for (const auto &facet : edge.delaunayFacets)
        {
            int iFacet = facet.second;
            Cell_handle c = facet.first;
            int d1 = (iFacet + 1) % 4;
            int d2 = (iFacet + 2) % 4;
            int d3 = (iFacet + 3) % 4;

            Vertex_handle p1 = c->vertex(d1);
            Vertex_handle p2 = c->vertex(d2);
            Vertex_handle p3 = c->vertex(d3);

            int iOrient = get_orientation(iFacet, v1, v2, v1_val, v2_val);
            generate_triangle(p1, p2, p3, iOrient, dt.is_infinite(c), dualTriangles);
        }
    }
}

//! @brief Processes a ray edge for dual triangle computation.
/*!
 * Intersects the ray with the bounding box, checks bipolarity, and generates
 * triangles for associated Delaunay facets.
 *
 * @param ray The ray edge to process.
 * @param edge The CGAL object representing the edge.
 * @param bbox The bounding box for intersection.
 * @param grid The scalar grid for interpolation.
 * @param isovalue The isovalue for bipolarity check.
 * @param dt The Delaunay triangulation.
 * @param dualTriangles Vector to store generated triangles.
 */
static void process_ray_edge(
    VoronoiEdge &edge,
    VoronoiDiagram &vd,
    CGAL::Epick::Iso_cuboid_3 &bbox,
    UnifiedGrid &grid,
    float isovalue,
    Delaunay &dt,
    std::vector<DelaunayTriangle> &dualTriangles)
{
    Ray3 ray;
    CGAL::assign(ray, edge.edgeObject);
    CGAL::Object intersectObj = CGAL::intersection(bbox, ray);
    Segment3 iseg;
    if (CGAL::assign(iseg, intersectObj))
    {
        Point v1 = iseg.source();
        Point v2 = iseg.target();
        int idx_v1 = edge.vertex1;
        float v1_val = vd.vertices[idx_v1].value;
        float iPt_value = trilinear_interpolate(adjust_outside_bound_points(v2, grid, v1, v2), grid);

        if (is_bipolar(v1_val, iPt_value, isovalue))
        {
            Point positive = (v1_val >= iPt_value) ? v1 : v2;

            for (const auto &facet : edge.delaunayFacets)
            {
                Facet mirror_f = dt.mirror_facet(facet);
                Object e = dt.dual(facet);

                int iFacet = facet.second;
                Cell_handle c = facet.first;
                int d1 = (iFacet + 1) % 4;
                int d2 = (iFacet + 2) % 4;
                int d3 = (iFacet + 3) % 4;

                Vertex_handle p1 = c->vertex(d1);
                Vertex_handle p2 = c->vertex(d2);
                Vertex_handle p3 = c->vertex(d3);

                int iOrient = get_orientation(iFacet, v1, v2, v1_val, iPt_value);
                generate_triangle(p1, p2, p3, iOrient, dt.is_infinite(c), dualTriangles);
            }
        }
    }
}

//! @brief Processes a line edge for dual triangle computation.
/*!
 * Intersects the line with the bounding box, checks bipolarity, and generates
 * triangles for associated Delaunay facets.
 *
 * @param line The line edge to process.
 * @param edge The CGAL object representing the edge.
 * @param grid The scalar grid for interpolation.
 * @param isovalue The isovalue for bipolarity check.
 * @param bbox The bounding box for intersection.
 * @param dt The Delaunay triangulation.
 * @param dualTriangles Vector to store generated triangles.
 */
static void process_line_edge(
    const Line3 &line,
    VoronoiEdge &edge,
    UnifiedGrid &grid,
    float isovalue,
    CGAL::Epick::Iso_cuboid_3 &bbox,
    Delaunay &dt,
    std::vector<DelaunayTriangle> &dualTriangles)
{
    CGAL::Object intersectObj = CGAL::intersection(bbox, line);
    Segment3 iseg;
    if (CGAL::assign(iseg, intersectObj))
    {
        Point intersection1 = iseg.source();
        Point intersection2 = iseg.target();

        float iPt1_val = trilinear_interpolate(adjust_outside_bound_points(intersection1, grid, intersection1, intersection2), grid);
        float iPt2_val = trilinear_interpolate(adjust_outside_bound_points(intersection2, grid, intersection1, intersection2), grid);

        if (is_bipolar(iPt1_val, iPt2_val, isovalue))
        {
            Point positive = (iPt1_val >= iPt2_val) ? intersection1 : intersection2;

            for (const auto &facet : edge.delaunayFacets)
            {
                int iFacet = facet.second;
                Cell_handle c = facet.first;
                int d1 = (iFacet + 1) % 4;
                int d2 = (iFacet + 2) % 4;
                int d3 = (iFacet + 3) % 4;

                Vertex_handle p1 = c->vertex(d1);
                Vertex_handle p2 = c->vertex(d2);
                Vertex_handle p3 = c->vertex(d3);

                int iOrient = get_orientation(iFacet, intersection1, intersection2, iPt1_val, iPt2_val);
                generate_triangle(p1, p2, p3, iOrient, dt.is_infinite(c), dualTriangles);
            }
        }
    }
}

//! @brief Computes the dual triangles for the final mesh in the single isovertex case.
/*!
 * Iterates over Voronoi edges, processes segments, rays, and lines to generate
 * Delaunay triangles dual to bipolar edges.
 *
 * @param iso_surface Instance of IsoSurface to store triangles.
 * @param voronoi_edges Vector of Voronoi edges.
 * @param bbox Bounding box of the computational domain.
 * @param dt Delaunay triangulation structure.
 * @param grid Scalar grid containing scalar values.
 * @param isovalue The isovalue used for computing.
 */
void compute_dual_triangles(
    IsoSurface &iso_surface,
    VoronoiDiagram &vd,
    CGAL::Epick::Iso_cuboid_3 &bbox,
    Delaunay &dt,
    UnifiedGrid &grid,
    float isovalue)
{
    std::vector<DelaunayTriangle> dualTriangles;

    auto vertex_is_valid = [&](int idx) -> bool
    {
        return idx >= 0 && idx < static_cast<int>(vd.vertices.size());
    };

    for (size_t edgeIdx = 0; edgeIdx < vd.edges.size(); ++edgeIdx)
    {
        auto &edge = vd.edges[edgeIdx];
        Segment3 seg;
        Ray3 ray;
        Line3 line;

        if (edge.type == 0)
        {
            if (!vertex_is_valid(edge.vertex1) ||
                !vertex_is_valid(edge.vertex2) ||
                edge.vertex1 == edge.vertex2)
            {
                if (debug)
                {
                    std::cerr << "[ISO] Skipping segment edge " << edgeIdx
                              << " with invalid vertices (" << edge.vertex1
                              << ", " << edge.vertex2 << ")\n";
                }
                continue;
            }
            process_segment_edge(edge, vd, isovalue, dt, dualTriangles);
        }
        else if (edge.type == 1)
        {
            if (!vertex_is_valid(edge.vertex1))
            {
                if (debug)
                {
                    std::cerr << "[ISO] Skipping ray edge " << edgeIdx
                              << " with invalid vertex " << edge.vertex1 << "\n";
                }
                continue;
            }
            process_ray_edge(edge, vd, bbox, grid, isovalue, dt, dualTriangles);
        }
        else if (edge.type == 2)
        {
            process_line_edge(line, edge, grid, isovalue, bbox, dt, dualTriangles);
        }
    }

    iso_surface.isosurfaceTrianglesSingle = dualTriangles;
}

//! @brief Computes isosurface vertices for the single-isovertex case.
/*!
 * Directly uses pre-computed accurate iso-crossing points as isosurface vertices.
 * The accurate iso-crossings were already calculated during the separation phase.
 */
void compute_isosurface_vertices_single(UnifiedGrid &grid, float isovalue, IsoSurface &iso_surface, std::vector<Point> &activeCubeIsoCrossingPoints)
{
    // Directly use the accurate iso-crossing points as isosurface vertices
    // These were pre-computed during the separation phase as the centroid of edge-isovalue intersections
    iso_surface.isosurfaceVertices = activeCubeIsoCrossingPoints;
}
