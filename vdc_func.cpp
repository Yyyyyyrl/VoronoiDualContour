//! @file vdc_func.cpp
//! @brief Implementation of functions for Voronoi Diagram and Isosurface computation.

#include "vdc_func.h"

// Helper for positive mod
static int positive_mod(int val, int mod)
{
    int res = val % mod;
    return (res < 0) ? res + mod : res;
}

//! @brief Generates a Delaunay triangle based on orientation and cell finiteness.
/*!
 * Adds a triangle to the dualTriangles vector, adjusting vertex order based on
 * the orientation and whether the cell is infinite.
 *
 * @param p1 First vertex of the triangle.
 * @param p2 Second vertex of the triangle.
 * @param p3 Third vertex of the triangle.
 * @param iOrient Orientation value determining vertex order.
 * @param isInfinite Flag indicating if the cell is infinite.
 * @param dualTriangles Vector to store the generated triangles.
 */
static void generate_triangle(
    const Vertex_handle &p1, const Vertex_handle &p2, const Vertex_handle &p3,
    int iOrient, bool isInfinite,
    std::vector<DelaunayTriangle> &dualTriangles)
{
    if (isInfinite)
    {
        if (iOrient < 0)
        {
            dualTriangles.push_back(DelaunayTriangle(p1, p2, p3));
        }
        else
        {
            dualTriangles.push_back(DelaunayTriangle(p1, p3, p2));
        }
    }
    else
    {
        if (iOrient >= 0)
        {
            dualTriangles.push_back(DelaunayTriangle(p1, p2, p3));
        }
        else
        {
            dualTriangles.push_back(DelaunayTriangle(p1, p3, p2));
        }
    }
}

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

    for (auto &edge : vd.edges)
    {
        Segment3 seg;
        Ray3 ray;
        Line3 line;

        if (edge.type == 0)
        {
            process_segment_edge(edge, vd, isovalue, dt, dualTriangles);
        }
        else if (edge.type == 1)
        {
            process_ray_edge(edge, vd, bbox, grid, isovalue, dt, dualTriangles);
        }
        else if (edge.type == 2)
        {
            process_line_edge(line, edge, grid, isovalue, bbox, dt, dualTriangles);
        }
    }

    iso_surface.isosurfaceTrianglesSingle = dualTriangles;
}

static inline int select_isovertex_from_cell_edge(
    const VoronoiDiagram &voronoiDiagram,
    int cellIndex, int globalEdgeIndex)
{
    // Lookup the VoronoiCellEdge
    auto it = voronoiDiagram.cellEdgeLookup.find(std::make_pair(cellIndex, globalEdgeIndex));
    if (it == voronoiDiagram.cellEdgeLookup.end())
    {
        // No such cell-edge found, pass
        std::cout << "didn't find cell-edge for edge " << globalEdgeIndex << std::endl;
        return -1;
    }
    int ceIdx = it->second;
    int starting = ceIdx;
    VoronoiCellEdge cellEdge = voronoiDiagram.cellEdges[ceIdx];

    while (cellEdge.cycleIndices.empty())
    {
        if (cellEdge.nextCellEdge == starting)
        {
            return -1;
        }
        else
        {
            cellEdge = voronoiDiagram.cellEdges[cellEdge.nextCellEdge];
        }
    }

    const VoronoiCell &vc = voronoiDiagram.cells[cellIndex];
    // The final isosurface vertex index in iso_surface.isosurfaceVertices:
    int isoVtxIndex = vc.isoVertexStartIndex + cellEdge.cycleIndices[0];
    return isoVtxIndex;
}

//! @brief Generates a triangle for the multi-isovertex case.
/*!
 * Adds a triangle to the isosurface's triangle list, adjusting vertex order based
 * on orientation, and logs problematic triangles if invalid.
 *
 * @param iso_surface The isosurface to store the triangle.
 * @param idx1 First vertex index.
 * @param idx2 Second vertex index.
 * @param idx3 Third vertex index.
 * @param iOrient Orientation value determining vertex order.
 * @param isValid Flag indicating if the triangle is valid.
 */
static void generate_triangle_multi(
    IsoSurface &iso_surface,
    int idx1, int idx2, int idx3,
    int iOrient,
    bool isValid)
{
    if (isValid)
    {
        if (iOrient < 0)
        {
            iso_surface.isosurfaceTrianglesMulti.emplace_back(idx1, idx2, idx3);
        }
        else
        {
            iso_surface.isosurfaceTrianglesMulti.emplace_back(idx1, idx3, idx2);
        }
    }
}

//! @brief Selects isovertices for a Delaunay facet.
/*!
 * Retrieves isovertex indices for the three vertices of a facet, ensuring they
 * are valid and not dummy vertices.
 *
 * @param voronoiDiagram The Voronoi diagram containing cell data.
 * @param facet The Delaunay facet to process.
 * @param globalEdgeIndex The global edge index for isovertex selection.
 * @param idx1 Output parameter for the first vertex index.
 * @param idx2 Output parameter for the second vertex index.
 * @param idx3 Output parameter for the third vertex index.
 * @param cellIndex1 Output parameter for the first cell index.
 * @param cellIndex2 Output parameter for the second cell index.
 * @param cellIndex3 Output parameter for the third cell index.
 * @return True if the vertices are valid, false otherwise.
 */
static bool select_isovertices(
    const VoronoiDiagram &voronoiDiagram,
    const Facet &facet,
    int globalEdgeIndex,
    int &idx1, int &idx2, int &idx3,
    int &cellIndex1, int &cellIndex2, int &cellIndex3)
{
    int iFacet = facet.second;
    Cell_handle c = facet.first;
    int d1 = (iFacet + 1) % 4;
    int d2 = (iFacet + 2) % 4;
    int d3 = (iFacet + 3) % 4;

    Vertex_handle delaunay_vertex1 = c->vertex(d1);
    Vertex_handle delaunay_vertex2 = c->vertex(d2);
    Vertex_handle delaunay_vertex3 = c->vertex(d3);

    if (delaunay_vertex1->info().is_dummy || delaunay_vertex2->info().is_dummy || delaunay_vertex3->info().is_dummy)
    {
        return false;
    }

    cellIndex1 = delaunay_vertex1->info().voronoiCellIndex;
    cellIndex2 = delaunay_vertex2->info().voronoiCellIndex;
    cellIndex3 = delaunay_vertex3->info().voronoiCellIndex;

    idx1 = select_isovertex_from_cell_edge(voronoiDiagram, cellIndex1, globalEdgeIndex);
    idx2 = select_isovertex_from_cell_edge(voronoiDiagram, cellIndex2, globalEdgeIndex);
    idx3 = select_isovertex_from_cell_edge(voronoiDiagram, cellIndex3, globalEdgeIndex);

    return (idx1 != idx2 && idx2 != idx3 && idx1 != idx3 && idx1 >= 0 && idx2 >= 0 && idx3 >= 0);
}

//! @brief Processes a segment edge for multi-isovertex triangle computation.
/*!
 * Checks if the segment is bipolar, retrieves the global edge index, and generates
 * triangles for associated Delaunay facets.
 *
 * @param seg The segment edge to process.
 * @param edge The CGAL object representing the edge.
 * @param voronoiDiagram The Voronoi diagram containing edge and cell data.
 * @param isovalue The isovalue for bipolarity check.
 * @param iso_surface The isosurface to store triangles.
 */
static void process_segment_edge_multi(
    VoronoiEdge edge,
    VoronoiDiagram &voronoiDiagram,
    float isovalue,
    IsoSurface &iso_surface)
{

    int idx_v1 = edge.vertex1;
    int idx_v2 = edge.vertex2;
    Point v1 = voronoiDiagram.vertices[idx_v1].coord;
    Point v2 = voronoiDiagram.vertices[idx_v2].coord;
    float val1 = voronoiDiagram.vertices[idx_v1].value;
    float val2 = voronoiDiagram.vertices[idx_v2].value;

    if (is_bipolar(val1, val2, isovalue))
    {
        if (idx_v1 > idx_v2)
            std::swap(idx_v1, idx_v2);
        auto itEdge = voronoiDiagram.segmentVertexPairToEdgeIndex.find(std::make_pair(idx_v1, idx_v2));
        if (itEdge == voronoiDiagram.segmentVertexPairToEdgeIndex.end())
            return;

        int globalEdgeIndex = itEdge->second;

        for (const auto &facet : edge.delaunayFacets)
        {
            int idx1, idx2, idx3, cellIndex1, cellIndex2, cellIndex3;
            bool isValid = select_isovertices(voronoiDiagram, facet, globalEdgeIndex, idx1, idx2, idx3, cellIndex1, cellIndex2, cellIndex3);
            int iOrient = get_orientation(facet.second, v1, v2, val1, val2);
            generate_triangle_multi(iso_surface, idx1, idx2, idx3, iOrient, isValid);
        }
    }
}

//! @brief Processes a ray edge for multi-isovertex triangle computation.
/*!
 * Intersects the ray with the bounding box, checks bipolarity, and generates
 * triangles using the first isovertex from each cell.
 *
 * @param ray The ray edge to process.
 * @param edge The CGAL object representing the edge.
 * @param voronoiDiagram The Voronoi diagram containing edge and cell data.
 * @param grid The scalar grid for interpolation.
 * @param isovalue The isovalue for bipolarity check.
 * @param bbox The bounding box for intersection.
 * @param iso_surface The isosurface to store triangles.
 */
static void process_ray_edge_multi(
    int source_pt,
    const Ray3 &ray,
    std::vector<Facet> dualDelaunayFacets,
    VoronoiDiagram &voronoiDiagram,
    UnifiedGrid &grid,
    float isovalue,
    CGAL::Epick::Iso_cuboid_3 &bbox,
    IsoSurface &iso_surface)
{
    CGAL::Object intersectObj = CGAL::intersection(bbox, ray);
    Segment3 iseg;
    if (CGAL::assign(iseg, intersectObj))
    {
        Point v1 = ray.source();
        Point v2 = iseg.target();
        int idx_v1 = source_pt;
        float val1 = voronoiDiagram.vertices[idx_v1].value;
        float val2 = trilinear_interpolate(v2, grid);

        if (is_bipolar(val1, val2, isovalue))
        {

            for (const auto &facet : dualDelaunayFacets)
            {
                int iFacet = facet.second;
                Cell_handle c = facet.first;
                int d1 = (iFacet + 1) % 4;
                int d2 = (iFacet + 2) % 4;
                int d3 = (iFacet + 3) % 4;

                Vertex_handle delaunay_vertex1 = c->vertex(d1);
                Vertex_handle delaunay_vertex2 = c->vertex(d2);
                Vertex_handle delaunay_vertex3 = c->vertex(d3);

                if (delaunay_vertex1->info().is_dummy || delaunay_vertex2->info().is_dummy || delaunay_vertex3->info().is_dummy)
                    continue;

                int cellIndex1 = delaunay_vertex1->info().voronoiCellIndex;
                int cellIndex2 = delaunay_vertex2->info().voronoiCellIndex;
                int cellIndex3 = delaunay_vertex3->info().voronoiCellIndex;

                VoronoiCell &vc1 = voronoiDiagram.cells[cellIndex1];
                VoronoiCell &vc2 = voronoiDiagram.cells[cellIndex2];
                VoronoiCell &vc3 = voronoiDiagram.cells[cellIndex3];

                int idx1 = vc1.isoVertexStartIndex;
                int idx2 = vc2.isoVertexStartIndex;
                int idx3 = vc3.isoVertexStartIndex;

                int iOrient = get_orientation(iFacet, v1, v2, val1, val2);
                bool isValid = (idx1 != idx2 && idx2 != idx3 && idx1 != idx3);
                generate_triangle_multi(iso_surface, idx1, idx2, idx3, iOrient, isValid);
            }
        }
    }
}

//! @brief Processes a line edge for multi-isovertex triangle computation.
/*!
 * Intersects the line with the bounding box, checks bipolarity, and generates
 * triangles using the first isovertex from each cell.
 *
 * @param line The line edge to process.
 * @param edge The CGAL object representing the edge.
 * @param voronoiDiagram The Voronoi diagram containing edge and cell data.
 * @param grid The scalar grid for interpolation.
 * @param isovalue The isovalue for bipolarity check.
 * @param bbox The bounding box for intersection.
 * @param iso_surface The isosurface to store triangles.
 */
static void process_line_edge_multi(
    const Line3 &line,
    std::vector<Facet> dualDelaunayFacets,
    VoronoiDiagram &voronoiDiagram,
    UnifiedGrid &grid,
    float isovalue,
    CGAL::Epick::Iso_cuboid_3 &bbox,
    IsoSurface &iso_surface)
{
    CGAL::Object intersectObj = CGAL::intersection(bbox, line);
    Segment3 iseg;
    if (CGAL::assign(iseg, intersectObj))
    {
        Point v1 = iseg.source();
        Point v2 = iseg.target();
        float val1 = trilinear_interpolate(v1, grid);
        float val2 = trilinear_interpolate(v2, grid);

        if (is_bipolar(val1, val2, isovalue))
        {

            for (const auto &facet : dualDelaunayFacets)
            {
                int iFacet = facet.second;
                Cell_handle c = facet.first;
                int d1 = (iFacet + 1) % 4;
                int d2 = (iFacet + 2) % 4;
                int d3 = (iFacet + 3) % 4;

                Vertex_handle delaunay_vertex1 = c->vertex(d1);
                Vertex_handle delaunay_vertex2 = c->vertex(d2);
                Vertex_handle delaunay_vertex3 = c->vertex(d3);

                if (delaunay_vertex1->info().is_dummy || delaunay_vertex2->info().is_dummy || delaunay_vertex3->info().is_dummy)
                    continue;

                int cellIndex1 = delaunay_vertex1->info().voronoiCellIndex;
                int cellIndex2 = delaunay_vertex2->info().voronoiCellIndex;
                int cellIndex3 = delaunay_vertex3->info().voronoiCellIndex;

                VoronoiCell &vc1 = voronoiDiagram.cells[cellIndex1];
                VoronoiCell &vc2 = voronoiDiagram.cells[cellIndex2];
                VoronoiCell &vc3 = voronoiDiagram.cells[cellIndex3];

                int idx1 = vc1.isoVertexStartIndex;
                int idx2 = vc2.isoVertexStartIndex;
                int idx3 = vc3.isoVertexStartIndex;

                int iOrient = get_orientation(iFacet, v1, v2, val1, val2);
                bool isValid = (idx1 != idx2 && idx2 != idx3 && idx1 != idx3);
                generate_triangle_multi(iso_surface, idx1, idx2, idx3, iOrient, isValid);
            }
        }
    }
}

//! @brief Computes the dual triangles for the final mesh in the multi-isovertex case.
/*!
 * Iterates over Voronoi edges, processes segments, rays, and lines to generate
 * Delaunay triangles dual to bipolar edges for multiple isovalues.
 *
 * @param voronoiDiagram The Voronoi diagram to compute from.
 * @param bbox Bounding box of the computational domain.
 * @param grid Scalar grid containing scalar values.
 * @param isovalue The isovalue for mesh computation.
 * @param iso_surface Instance of IsoSurface containing the isosurface vertices and faces.
 */
void compute_dual_triangles_multi(
    VoronoiDiagram &voronoiDiagram,
    CGAL::Epick::Iso_cuboid_3 &bbox,
    UnifiedGrid &grid,
    float isovalue,
    IsoSurface &iso_surface)
{
    for (const auto &edge : voronoiDiagram.edges)
    {
        Segment3 seg;
        Ray3 ray;
        Line3 line;
        std::vector<Facet> dualDelaunayFacets = edge.delaunayFacets;

        // std::cout << "[DEBUG] processing edge (" << edge.vertex1 << ", " << edge.vertex2 << "), type = " << edge.type << std::endl;
        if (edge.type == 0)
        {
            process_segment_edge_multi(edge, voronoiDiagram, isovalue, iso_surface);
        }
        else if (edge.type == 1)
        {
            CGAL::assign(ray, edge.edgeObject);
            int source = edge.vertex1;
            process_ray_edge_multi(source, ray, dualDelaunayFacets, voronoiDiagram, grid, isovalue, bbox, iso_surface);
        }
        else if (edge.type == 2)
        {
            CGAL::assign(line, edge.edgeObject);
            process_line_edge_multi(line, dualDelaunayFacets, voronoiDiagram, grid, isovalue, bbox, iso_surface);
        }
    }
}

//! @brief Computes isosurface vertices for the single-isovertex case.
void compute_isosurface_vertices_single(UnifiedGrid &grid, float isovalue, IsoSurface &iso_surface, std::vector<Point> &activeCubeCenters)
{
    const int cubeVertices[8][3] = {
        {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    const int cubeEdges[12][2] = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6}, {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};

    int vertexIndex = 0;
    for (const auto &center : activeCubeCenters)
    {
        std::vector<Point> intersectionPoints;
        std::array<float, 8> scalarValues;

        for (int i = 0; i < 8; i++)
        {
            Point vertex(
                center.x() + (cubeVertices[i][0] - 0.5f) * grid.dx,
                center.y() + (cubeVertices[i][1] - 0.5f) * grid.dy,
                center.z() + (cubeVertices[i][2] - 0.5f) * grid.dz);
            scalarValues[i] = grid.get_scalar_value_at_point(vertex);
        }

        for (const auto &edge : cubeEdges)
        {
            int idx1 = edge[0];
            int idx2 = edge[1];
            float val1 = scalarValues[idx1];
            float val2 = scalarValues[idx2];

            if (is_bipolar(val1, val2, isovalue))
            {
                Point p1(
                    center.x() + (cubeVertices[idx1][0] - 0.5f) * grid.dx,
                    center.y() + (cubeVertices[idx1][1] - 0.5f) * grid.dy,
                    center.z() + (cubeVertices[idx1][2] - 0.5f) * grid.dz);
                Point p2(
                    center.x() + (cubeVertices[idx2][0] - 0.5f) * grid.dx,
                    center.y() + (cubeVertices[idx2][1] - 0.5f) * grid.dy,
                    center.z() + (cubeVertices[idx2][2] - 0.5f) * grid.dz);

                Point intersect = interpolate(p1, p2, val1, val2, isovalue, grid);
                intersectionPoints.push_back(intersect);
            }
        }

        if (!intersectionPoints.empty())
        {
            Point centroid = compute_centroid(intersectionPoints);
            iso_surface.isosurfaceVertices.push_back(centroid);
        }
        else
            std::cerr << "[WARNING] No intersection points for cube at center: " << center << std::endl;
    }
}

//! @brief Collects midpoints for bipolar edges in a Voronoi cell's facets.
/*!
 * Identifies bipolar edges, computes their midpoints, and stores them along with
 * facet associations.
 *
 * @param vc The Voronoi cell to process.
 * @param voronoiDiagram The Voronoi diagram containing facet and vertex data.
 * @param isovalue The isovalue for bipolarity check.
 * @param midpoints Vector to store computed midpoints.
 * @param edge_to_midpoint_index Map linking edge keys to midpoint indices.
 * @param facet_midpoint_indices Vector storing midpoint indices for each facet.
 */
static void collect_midpoints(
    VoronoiCell &vc,
    VoronoiDiagram &voronoiDiagram,
    float isovalue,
    std::vector<MidpointNode> &midpoints,
    std::map<std::pair<int, int>, int> &edge_to_midpoint_index,
    std::vector<std::vector<int>> &facet_midpoint_indices)
{
    for (size_t i = 0; i < vc.facet_indices.size(); ++i)
    {
        int facet_index = vc.facet_indices[i];
        VoronoiCellFacet &facet = voronoiDiagram.facets[facet_index];
        size_t num_vertices = facet.vertices_indices.size();

        std::vector<int> current_facet_midpoints;

        for (size_t j = 0; j < num_vertices; ++j)
        {
            size_t idx1 = j;
            size_t idx2 = (j + 1) % num_vertices;

            float val1 = voronoiDiagram.vertices[facet.vertices_indices[idx1]].value;
            float val2 = voronoiDiagram.vertices[facet.vertices_indices[idx2]].value;

            if (is_bipolar(val1, val2, isovalue))
            {
                int vertex_index1 = facet.vertices_indices[idx1];
                int vertex_index2 = facet.vertices_indices[idx2];

                Point p1 = voronoiDiagram.vertices[vertex_index1].coord;
                Point p2 = voronoiDiagram.vertices[vertex_index2].coord;

                double t = (isovalue - val1) / (val2 - val1);
                Point midpoint = p1 + (p2 - p1) * t;

                auto edge_key = std::make_pair(std::min(vertex_index1, vertex_index2),
                                               std::max(vertex_index1, vertex_index2));

                if (edge_to_midpoint_index.find(edge_key) == edge_to_midpoint_index.end())
                {
                    int globalEdgeIndex = -1;
                    auto iter_glob = voronoiDiagram.segmentVertexPairToEdgeIndex.find(edge_key);
                    if (iter_glob != voronoiDiagram.segmentVertexPairToEdgeIndex.end())
                    {
                        globalEdgeIndex = iter_glob->second;
                    }

                    MidpointNode node;
                    node.point = midpoint;
                    node.facet_index = facet_index;
                    node.cycle_index = -1;
                    node.global_edge_index = globalEdgeIndex;

                    midpoints.push_back(node);
                    int midpoint_index = midpoints.size() - 1;
                    edge_to_midpoint_index[edge_key] = midpoint_index;
                    current_facet_midpoints.push_back(midpoint_index);
                }
                else
                {
                    int midpoint_index = edge_to_midpoint_index[edge_key];
                    current_facet_midpoints.push_back(midpoint_index);
                }
            }
        }

        facet_midpoint_indices.push_back(current_facet_midpoints);
    }
}

// Helper: get undirected edge key from a global facet boundary slot
static inline std::pair<int, int>
facet_slot_edge_key(const VoronoiFacet &gf, int slot)
{
    int m = (int)gf.vertices_indices.size();
    int a = gf.vertices_indices[slot];
    int b = gf.vertices_indices[(slot + 1) % m];
    return (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
}

// It reuses the canonical pairing from global_facets.
static void connect_midpoints_via_global_matches(
    const VoronoiDiagram &vd,
    const VoronoiCell &vc,
    const std::map<std::pair<int, int>, int> &edge_to_midpoint_index,
    std::vector<MidpointNode> &midpoints)
{
    for (int cf : vc.facet_indices)
    {
        if (cf < 0 || cf >= (int)vd.facets.size())
            continue;

        const auto &cellFacet = vd.facets[cf];
        int vfi = cellFacet.voronoi_facet_index;
        if (vfi < 0 || vfi >= (int)vd.global_facets.size())
            continue;

        const auto &gf = vd.global_facets[vfi];
        for (const auto &pr : gf.bipolar_matches)
        {
            int sA = pr.first; // boundary slot in global facet order
            int sB = pr.second;

            auto ekA = facet_slot_edge_key(gf, sA);
            auto ekB = facet_slot_edge_key(gf, sB);

            auto itA = edge_to_midpoint_index.find(ekA);
            auto itB = edge_to_midpoint_index.find(ekB);

            if (itA == edge_to_midpoint_index.end() || itB == edge_to_midpoint_index.end()){
                
                continue; // this cell doesn't have both midpoints (e.g., not bipolar here), skip
            }

            int iA = itA->second, iB = itB->second;
            midpoints[iA].connected_to.push_back(iB);
            midpoints[iB].connected_to.push_back(iA);
        }
    }
}

//! @brief Extracts cycles from the midpoint connectivity graph.
/*!
 * Uses depth-first search to identify closed cycles in the midpoint graph.
 *
 * @param midpoints Vector of midpoints with connectivity information.
 * @param cycles Vector to store the extracted cycles as lists of midpoint indices.
 */
static void extract_cycles(
    const std::vector<MidpointNode> &midpoints,
    std::vector<std::vector<int>> &cycles)
{
    const int n = (int)midpoints.size();
    std::vector<char> used(n, 0);

    // build unique adjacency
    std::vector<std::vector<int>> adj(n);
    for (int i = 0; i < n; ++i) {
        std::unordered_set<int> uniq(midpoints[i].connected_to.begin(), midpoints[i].connected_to.end());
        adj[i].assign(uniq.begin(), uniq.end());
    }

    //alert if degree != 2
    for (int i = 0; i < n; ++i) {
        if (!adj[i].empty() && adj[i].size() != 2) {
            std::cerr << "[warn] midpoint " << i << " has degree " << adj[i].size() << " (expected 2)\n";
        }
    }

    for (int s = 0; s < n; ++s) {
        if (used[s] || adj[s].empty()) continue;
        int prev = -1, cur = s;
        std::vector<int> cyc;
        while (true) {
            used[cur] = 1;
            cyc.push_back(cur);

            if (adj[cur].empty()) break; // broken

            int nxt = (adj[cur].size()==1) ? adj[cur][0] : (adj[cur][0]==prev ? adj[cur][1] : adj[cur][0]);

            if (nxt == s) { cycles.push_back(cyc); break; }        // closed
            if (nxt < 0 || nxt >= n || used[nxt]) {                 // broken/self-intersecting
                cycles.push_back(cyc); break;
            }

            prev = cur;
            cur = nxt;
        }
    }
}

// TODO: Rename the routine to show its true functionality, as it computes the isovertices AND the cycle centroids.
//! @brief Computes centroids for cycles and updates the isosurface.
/*!
 * Calculates the centroid for each cycle, updates the Voronoi cell's cycles,
 * and adds the centroid to the isosurface vertices.
 *
 * @param vc The Voronoi cell to update.
 * @param voronoiDiagram The Voronoi diagram for edge lookup.
 * @param midpoints Vector of midpoints used for centroid computation.
 * @param cycles Vector of cycles as lists of midpoint indices.
 * @param iso_surface The isosurface to store vertices.
 */
static void compute_cycle_centroids(
    VoronoiCell &vc,
    VoronoiDiagram &voronoiDiagram,
    std::vector<MidpointNode> &midpoints,
    const std::vector<std::vector<int>> &cycles,
    IsoSurface &iso_surface)
{
    vc.isoVertexStartIndex = iso_surface.isosurfaceVertices.size();
    vc.numIsoVertices = cycles.size();

    int cycIdx = 0;
    for (const auto &single_cycle : cycles)
    {
        Cycle cycle;
        cycle.voronoi_cell_index = vc.cellIndex;
        cycle.midpoint_indices = single_cycle;

        for (size_t i = 0; i < single_cycle.size(); ++i)
        {
            int idx1 = single_cycle[i];
            int idx2 = single_cycle[(i + 1) % single_cycle.size()];
            cycle.edges.emplace_back(idx1, idx2);
        }

        cycle.compute_centroid(midpoints);

        for (int ptIdx : single_cycle)
        {
            midpoints[ptIdx].cycle_index = cycIdx;

            int globalEdgeIdx = midpoints[ptIdx].global_edge_index;
            if (globalEdgeIdx >= 0)
            {
                std::pair<int, int> key = std::make_pair(vc.cellIndex, globalEdgeIdx);
                auto iter_cEdge = voronoiDiagram.cellEdgeLookup.find(key);
                if (iter_cEdge != voronoiDiagram.cellEdgeLookup.end())
                {
                    int cEdgeIdx = iter_cEdge->second;
                    auto &cyclesVec = voronoiDiagram.cellEdges[cEdgeIdx].cycleIndices;

                    if (std::find(cyclesVec.begin(), cyclesVec.end(), cycIdx) == cyclesVec.end())
                    {
                        cyclesVec.push_back(cycIdx);
                    }
                }
            }
        }

        vc.cycles.push_back(cycle);
        iso_surface.isosurfaceVertices.push_back(cycle.isovertex);
        cycIdx++;
    }
}

//! @brief Computes isosurface vertices for the multi-isovertex case.
/*!
 * Processes each Voronoi cell to identify bipolar edges, form cycles, and compute
 * centroids as isosurface vertices.
 *
 * @param voronoiDiagram The Voronoi diagram to compute vertices for.
 * @param isovalue The isovalue to use for vertex computation.
 * @param iso_surface Instance of IsoSurface containing the isosurface vertices and faces.
 */
void compute_isosurface_vertices_multi(VoronoiDiagram &voronoiDiagram, float isovalue, IsoSurface &iso_surface)
{
    voronoiDiagram.compute_bipolar_matches(isovalue);
    for (auto &vc : voronoiDiagram.cells)
    {
        std::vector<MidpointNode> midpoints;
        std::map<std::pair<int, int>, int> edge_to_midpoint_index;
        std::vector<std::vector<int>> facet_midpoint_indices;

        collect_midpoints(vc, voronoiDiagram, isovalue, midpoints, edge_to_midpoint_index, facet_midpoint_indices);
        // connect_midpoints(facet_midpoint_indices, midpoints);
        connect_midpoints_via_global_matches(voronoiDiagram, vc, edge_to_midpoint_index, midpoints);

        size_t num_edges_added = 0;
        for (auto &n : midpoints)
            num_edges_added += n.connected_to.size();
        if (num_edges_added == 0)
        {
            std::cerr << "[DBG] cell " << vc.cellIndex
                      << " had " << midpoints.size()
                      << " midpoints but 0 connections; check bipolar_matches and keys.\n";
        }
        std::vector<std::vector<int>> cycles;
        extract_cycles(midpoints, cycles);

        compute_cycle_centroids(vc, voronoiDiagram, midpoints, cycles, iso_surface);
    }
}

//! @brief Adds dummy points from a facet for Voronoi diagram bounding.
std::vector<Point> add_dummy_from_facet(const GRID_FACETS &facet, const UnifiedGrid &data_grid)
{
    std::vector<Point> points;

    // 2D slice dimension
    int dim0 = facet.axis_size[0];
    int dim1 = facet.axis_size[1];

    // For convenience
    int d = facet.orth_dir;
    int d1 = facet.axis_dir[0];
    int d2 = facet.axis_dir[1];

    // and localSize[] = (maxIndex[i] - minIndex[i] + 1)
    // The grid spacing in each dimension
    double dx[3] = {data_grid.dx, data_grid.dy, data_grid.dz};

    // Loop over the 2D slice
    for (int coord1 = 0; coord1 < dim1; coord1++)
    {
        for (int coord0 = 0; coord0 < dim0; coord0++)
        {
            if (!facet.CubeFlag(coord0, coord1))
                continue;

            // localX[d1] = coord0, localX[d2] = coord1
            int localX[3] = {0, 0, 0};
            localX[d1] = coord0;
            localX[d2] = coord1;

            // side=0 => localX[d] = 0, side=1 => localX[d] = localSize[d]-1
            localX[d] = (facet.side == 0) ? 0 : (facet.localSize[d] - 1);

            // Convert localX -> global indices
            int g[3];
            for (int i = 0; i < 3; i++)
            {
                g[i] = localX[i] + facet.minIndex[i];
            }

            // Compute center in real-world coordinates
            double cx = (g[0] + 0.5) * dx[0];
            double cy = (g[1] + 0.5) * dx[1];
            double cz = (g[2] + 0.5) * dx[2];

            // Offset by +/- dx[d]
            double offset = (facet.side == 0) ? -dx[d] : dx[d];
            if (d == 0)
                cx += offset;
            else if (d == 1)
                cy += offset;
            else
                cz += offset;

            points.emplace_back(cx, cy, cz);
        }
    }

    return points;
}

//! @brief Collects points for the Delaunay triangulation.
/*!
 * Gathers original points and dummy points from grid facets for multi-isovertex mode,
 * or uses only active cube centers for single-isovertex mode.
 *
 * @param grid The grid containing data.
 * @param grid_facets The grid facets for dummy point generation.
 * @param activeCubeCenters The list of center points of active cubes.
 * @param vdc_param The VDC_PARAM instance containing user input options.
 * @param delaunay_points Output vector for all points (original + dummy).
 */
static int collect_delaunay_points(UnifiedGrid &grid,
                                   const std::vector<std::vector<GRID_FACETS>> &grid_facets,
                                   const std::vector<Point> &activeCubeCenters,
                                   VDC_PARAM &vdc_param,
                                   std::vector<Point> &delaunay_points)
{
    delaunay_points = activeCubeCenters;
    int first_dummy_index = delaunay_points.size(); // Dummies start here

    if (vdc_param.multi_isov)
    {
        for (int d = 0; d < 3; ++d)
        { // Assuming 3 dimensions
            for (const auto &f : grid_facets[d])
            {
                auto pointsf = add_dummy_from_facet(f, grid);
                delaunay_points.insert(delaunay_points.end(), pointsf.begin(), pointsf.end());
            }
        }
    }
    return first_dummy_index;
}

//! @brief Inserts points into the Delaunay triangulation.
/*!
 * Inserts the collected points into the triangulation and writes debug output if enabled.
 *
 * @param dt The Delaunay triangulation to insert points into.
 * @param delaunay_points The points to insert.
 * @param activeCubeCenters The list of center points of active cubes.
 * @param vdc_param The VDC_PARAM instance containing user input options.
 */
static Vertex_handle insert_point_into_delaunay_triangulation(Delaunay &dt,
                                                              const Point &p,
                                                              int index,
                                                              bool is_dummy)
{
    // Insert point and retrieve handle
    Vertex_handle vh = dt.insert(p);
    // Immediately assign index and dummy status
    vh->info().index = index;
    vh->info().is_dummy = is_dummy;
    return vh;
}

//! @brief Constructs a Delaunay triangulation from a grid and grid facets.
/*!
 * Constructs a 3D Delaunay triangulation using the grid's scalar values
 * and the facets of the active cubes.
 *
 * @param dt The Delaunay triangulation instance.
 * @param grid The grid containing scalar values.
 * @param grid_facets The grid facets to use in constructing the triangulation.
 * @param vdc_param The VDC_PARAM instance holding user input options.
 * @param activeCubeCenters The list of center points of active cubes.
 */
void construct_delaunay_triangulation(Delaunay &dt, UnifiedGrid &grid, const std::vector<std::vector<GRID_FACETS>> &grid_facets, VDC_PARAM &vdc_param, std::vector<Point> &activeCubeCenters)
{
    std::clock_t start = std::clock();

    std::vector<Point> delaunay_points;
    size_t first_dummy_index = collect_delaunay_points(grid, grid_facets, activeCubeCenters, vdc_param, delaunay_points);

    std::clock_t after_collect = std::clock();
    double collect_time = static_cast<double>(after_collect - start) / CLOCKS_PER_SEC;
    std::cout << "[INFO] Time to build point list: " << collect_time << " seconds" << std::endl;

    std::cout << "[DEBUG] Number of vertices: " << delaunay_points.size() << std::endl;

    dt.clear();

    // Batch insert all points
    dt.insert(delaunay_points.begin(), delaunay_points.end());

    std::clock_t after_insert = std::clock();
    double insert_time = static_cast<double>(after_insert - after_collect) / CLOCKS_PER_SEC;
    std::cout << "[INFO] Vertices insert time: " << insert_time << " seconds" << std::endl;

    // Create map from point to original index
    std::map<Point, size_t> point_to_index;
    for (size_t i = 0; i < delaunay_points.size(); ++i)
    {
        point_to_index[delaunay_points[i]] = i;
    }

    // Assign info to vertices
    for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit)
    {
        const Point &p = vit->point();
        auto it = point_to_index.find(p);
        if (it == point_to_index.end())
        {
            std::cerr << "[ERROR] Vertex point not found in original points!" << std::endl;
            continue;
        }
        size_t original_index = it->second;
        vit->info().index = original_index;
        vit->info().is_dummy = (original_index >= first_dummy_index);
        vit->info().voronoiCellIndex = -1; // Initialize if needed
    }

    std::clock_t after_assign = std::clock();
    double assign_time = static_cast<double>(after_assign - after_insert) / CLOCKS_PER_SEC;
    std::cout << "[INFO] Time to assign vertex info: " << assign_time << " seconds" << std::endl;
}

//! @brief Constructs Voronoi vertices for the given voronoi Diagram instance.
void construct_voronoi_vertices(VoronoiDiagram &voronoiDiagram, Delaunay &dt)
{
    int vertexIndex = 0;
    voronoiDiagram.vertices.reserve(dt.number_of_finite_cells());
    for (Delaunay::Finite_cells_iterator cit = dt.finite_cells_begin(); cit != dt.finite_cells_end(); ++cit)
    {
        // Skip degenerate cells to avoid invalid circumcenters (e.g., NaN coordinates)
        if (is_degenerate(cit))
        {
            if (debug)
            {
                std::cerr << "[DEBUG] Skipping degenerate cell with vertices: "
                          << cit->vertex(0)->point() << ", "
                          << cit->vertex(1)->point() << ", "
                          << cit->vertex(2)->point() << ", "
                          << cit->vertex(3)->point() << "\n";
            }
            continue;
        }

        Point circumcenter = dt.dual(cit);

        if (std::isnan(circumcenter.x()) || std::isnan(circumcenter.y()) || std::isnan(circumcenter.z()))
        {
            if (debug)
            {
                std::cerr << "[DEBUG] Skipping cell with NaN circumcenter: vertices "
                          << cit->vertex(0)->point() << ", "
                          << cit->vertex(1)->point() << ", "
                          << cit->vertex(2)->point() << ", "
                          << cit->vertex(3)->point() << "\n";
            }
            continue;
        }
        VoronoiVertex vv(circumcenter);
        vv.index = vertexIndex;
        voronoiDiagram.vertices.push_back(vv);
        cit->info().dualVoronoiVertexIndex = vertexIndex;
        vertexIndex++;
    }
}

//! @brief Computes Voronoi Vertex values using scalar grid interpolation
void compute_voronoi_values(VoronoiDiagram &voronoiDiagram, UnifiedGrid &grid)
{
    for (size_t i = 0; i < voronoiDiagram.vertices.size(); ++i)
    {
        Point vertex = voronoiDiagram.vertices[i].coord;
        voronoiDiagram.vertices[i].value = trilinear_interpolate(vertex, grid);
    }
}

//! @brief Constructs Voronoi cells from the Delaunay triangulation.
void construct_voronoi_cells_as_convex_hull(VoronoiDiagram &voronoiDiagram, Delaunay &dt)
{
    int index = 0;
    for (auto delaunay_vertex = dt.finite_vertices_begin(); delaunay_vertex != dt.finite_vertices_end(); ++delaunay_vertex)
    {
        if (delaunay_vertex->info().is_dummy)
        {
            continue;
        }
        VoronoiCell vc(delaunay_vertex);
        vc.cellIndex = index;

        std::vector<Cell_handle> incident_cells;
        dt.finite_incident_cells(delaunay_vertex, std::back_inserter(incident_cells));

        // Collect vertex indices combinatorially
        std::set<int> unique_vertex_indices_set;
        for (Cell_handle ch : incident_cells)
        {
            if (dt.is_infinite(ch))
            {
                continue; // Skip infinite cells
            }
            // Use direct index instead of dual point + search
            int vertex_index = ch->info().dualVoronoiVertexIndex;
            unique_vertex_indices_set.insert(vertex_index);
        }

        // Copy unique indices to vector
        vc.vertices_indices.assign(unique_vertex_indices_set.begin(), unique_vertex_indices_set.end());

        // Build vertex_points and vector for lookup (allow duplicates by using first idx for matching points)
        std::vector<Point> vertex_points;
        std::vector<std::pair<Point, int>> point_index_pairs;
        for (int idx : vc.vertices_indices)
        {
            Point p = voronoiDiagram.vertices[idx].coord;
            vertex_points.push_back(p);
            // Check if point already added (approx equal), if not, add pair
            bool found = false;
            for (const auto &pair : point_index_pairs)
            {
                if (PointApproxEqual()(pair.first, p))
                {
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                point_index_pairs.emplace_back(p, idx); // Use first idx for this point
            }
        }

        // No remove duplicates: pass all to hull, it will handle

        CGAL::convex_hull_3(vertex_points.begin(), vertex_points.end(), vc.polyhedron);

        // Extract facets from polyhedron
        for (auto facet_it = vc.polyhedron.facets_begin();
             facet_it != vc.polyhedron.facets_end(); ++facet_it)
        {
            VoronoiCellFacet vf;
            auto h = facet_it->facet_begin();
            do
            {
                Point p = h->vertex()->point();
                // Linear lookup in point_index_pairs
                bool found = false;
                for (const auto &pair : point_index_pairs)
                {
                    if (PointApproxEqual()(pair.first, p))
                    {
                        vf.vertices_indices.push_back(pair.second);
                        found = true;
                        break;
                    }
                }
                if (!found)
                {
                    std::cerr << "[WARNING] Point not found during facet extraction: " << p << "\n";
                }
                ++h;
            } while (h != facet_it->facet_begin());

            // Skip degenerate facets
            if (vf.vertices_indices.size() < 3)
                continue;

            int facet_index = voronoiDiagram.facets.size();
            voronoiDiagram.facets.push_back(vf);
            vc.facet_indices.push_back(facet_index);
        }

        voronoiDiagram.cells.push_back(vc);
        delaunay_vertex->info().voronoiCellIndex = index;
        index++;
    }
}

//! @brief Creates a Voronoi cell for a Delaunay vertex.
/*!
 * Initializes a Voronoi cell with the given cell index and Delaunay vertex handle.
 *
 * @param delaunay_vertex The Delaunay vertex to create the cell for.
 * @param cellIndex The index to assign to the cell.
 * @return The initialized Voronoi cell.
 */
static VoronoiCell create_voronoi_cell(Vertex_handle delaunay_vertex, int cellIndex)
{
    VoronoiCell vc(delaunay_vertex);
    vc.cellIndex = cellIndex;
    return vc;
}

//! @brief Collects unique vertex indices from incident cells.
/*!
 * Retrieves the Voronoi vertex indices from cells incident to a Delaunay vertex,
 * applying the old-to-new vertex index mapping.
 *
 * @param dt The Delaunay triangulation.
 * @param delaunay_vertex The Delaunay vertex to process.
 * @param voronoiDiagram The Voronoi diagram containing vertex mappings.
 * @param vertices_indices Vector to store the collected vertex indices.
 */
static void collcet_cell_vertices(
    Delaunay &dt,
    Vertex_handle delaunay_vertex,
    VoronoiDiagram &voronoiDiagram,
    std::vector<int> &vertices_indices)
{
    std::vector<Cell_handle> incidentCells;
    dt.finite_incident_cells(delaunay_vertex, std::back_inserter(incidentCells));

    std::set<int> uniqueVertexIndices;
    for (Cell_handle c : incidentCells)
    {
        int vertex_index = c->info().dualVoronoiVertexIndex;
        uniqueVertexIndices.insert(vertex_index);
    }
    vertices_indices.assign(uniqueVertexIndices.begin(), uniqueVertexIndices.end());
}

//! @brief Builds a facet from an incident edge using cell circulators.
/*!
 * Constructs a Voronoi facet by collecting vertices around an incident edge,
 * ordering them cyclically, and assigning scalar values.
 *
 * @param dt The Delaunay triangulation.
 * @param ed The incident edge to process.
 * @param delaunay_vertex The Delaunay vertex associated with the cell.
 * @param voronoiDiagram The Voronoi diagram containing vertex and value data.
 * @param facet_indices Vector to store the facet index.
 * @return The constructed Voronoi facet, or an empty facet if invalid.
 */
static VoronoiCellFacet build_facet_from_edge(
    Delaunay &dt,
    const Edge &ed,
    Vertex_handle delaunay_vertex,
    VoronoiDiagram &voronoiDiagram,
    std::vector<int> &facet_indices,
    std::map<std::pair<int, int>, std::vector<int>> &edge_to_facets,
    int vcIdx)
{
    Cell_handle cell_ed = ed.first;
    int i = ed.second;
    int j = ed.third;
    Vertex_handle v1 = cell_ed->vertex(i);
    Vertex_handle v2 = cell_ed->vertex(j);

    Delaunay::Cell_circulator cc = dt.incident_cells(ed);
    Delaunay::Cell_circulator start = cc;
    std::vector<int> facetVertexIndices;
    int finite_cell_count = 0;
    bool involves_infinite = false;
    do
    {
        if (dt.is_infinite(cc))
        {
            involves_infinite = true;
        }
        else
        {
            int vertex_index = cc->info().dualVoronoiVertexIndex;
            facetVertexIndices.push_back(vertex_index);
            finite_cell_count++;
        }
        ++cc;
    } while (cc != start);

    std::set<int> unique_vertices(facetVertexIndices.begin(), facetVertexIndices.end());

    if (unique_vertices.size() >= 3)
    {
        // Clean duplicates
        std::vector<int> cleaned;
        for (size_t k = 0; k < facetVertexIndices.size(); ++k)
        {
            if (k == 0 || facetVertexIndices[k] != facetVertexIndices[k - 1])
            {
                cleaned.push_back(facetVertexIndices[k]);
            }
        }
        if (cleaned.size() > 1 && cleaned.front() == cleaned.back())
        {
            cleaned.pop_back();
        }
        if (cleaned.size() < 3)
        {
            std::cout << "[DEBUG] Degenerate after cleaning: " << cleaned.size() << " verts\n";
            return VoronoiCellFacet();
        }
        std::vector<int> orderedFacetVertices = std::move(cleaned);

        // Determine orientation of the facet using CGAL::orientation ( taking determinant )
        Point P1 = voronoiDiagram.vertices[orderedFacetVertices[0]].coord;
        Point P2 = voronoiDiagram.vertices[orderedFacetVertices[1]].coord;
        Point P3 = voronoiDiagram.vertices[orderedFacetVertices[2]].coord;
        Point site = delaunay_vertex->point();
        CGAL::Orientation orient = CGAL::orientation(P1, P2, P3, site);
        if (orient == CGAL::POSITIVE)
        {
            std::reverse(orderedFacetVertices.begin(), orderedFacetVertices.end());
        }
        else if (orient == CGAL::ZERO)
        {
        }

        VoronoiCellFacet facet;
        facet.vertices_indices = orderedFacetVertices;

        int facetIndex = voronoiDiagram.facets.size();
        voronoiDiagram.facets.push_back(facet);
        facet_indices.push_back(facetIndex);
        std::pair<int, int> edge_key = std::make_pair(
            std::min(v1->info().index, v2->info().index),
            std::max(v1->info().index, v2->info().index));
        edge_to_facets[edge_key].push_back(facetIndex);

        // Map vertex loop → cell-edge loop for this cell facet
        facet.cell_edge_indices.clear();
        const int n = (int)facet.vertices_indices.size();

        for (int i = 0; i < n; ++i)
        {
            const int a = facet.vertices_indices[i];
            const int b = facet.vertices_indices[(i + 1) % n];
            const int vmin = std::min(a, b), vmax = std::max(a, b);

            // 1) find global segment edge
            const auto eit = voronoiDiagram.segmentVertexPairToEdgeIndex.find({vmin, vmax});
            if (eit == voronoiDiagram.segmentVertexPairToEdgeIndex.end())
            {
                facet.cell_edge_indices.push_back(-1);
                continue;
            }
            const int globalEdge = eit->second;

            // 2) find the cell-edge for (this cell, globalEdge)
            const int thisCellIndex = vcIdx;
            const auto ceit = voronoiDiagram.cellEdgeLookup.find({thisCellIndex, globalEdge});
            facet.cell_edge_indices.push_back(ceit == voronoiDiagram.cellEdgeLookup.end() ? -1 : ceit->second);
        }
        return facet;
    }
    else
    {
        std::cout << "[DEBUG] Degenerate facet for edge with " << finite_cell_count << " finite cells\n";
        return VoronoiCellFacet();
    }
}

//! @brief Processes incident edges to build facets for a Voronoi cell.
/*!
 * Iterates over incident edges to construct facets and add them to the cell.
 *
 * @param dt The Delaunay triangulation.
 * @param delaunay_vertex The Delaunay vertex to process.
 * @param voronoiDiagram The Voronoi diagram to update.
 * @param vc The Voronoi cell to populate with facets.
 */
static void process_incident_edges(
    Delaunay &dt,
    Vertex_handle delaunay_vertex,
    VoronoiDiagram &voronoiDiagram,
    VoronoiCell &vc,
    std::map<std::pair<int, int>, std::vector<int>> &edge_to_facets)
{
    std::vector<Edge> incidentEdges;
    dt.incident_edges(delaunay_vertex, std::back_inserter(incidentEdges));

    for (const Edge &ed : incidentEdges)
    {
        // Count finite incident cells
        int finite_cell_count = 0;
        Delaunay::Cell_circulator cc = dt.incident_cells(ed);
        Delaunay::Cell_circulator start = cc;
        do
        {
            if (!dt.is_infinite(cc))
            {
                finite_cell_count++;
            }
            ++cc;
        } while (cc != start);

        if (finite_cell_count < 3)
        {
            // std::cout << "[INFO] Skipping edge with " << finite_cell_count << " finite incident cells (insufficient for interior facet)\n";
            continue;
        }

        // Build facet only if edge has 3+ finite cells
        VoronoiCellFacet facet = build_facet_from_edge(dt, ed, delaunay_vertex, voronoiDiagram, vc.facet_indices, edge_to_facets, vc.cellIndex);
        if (facet.vertices_indices.empty())
        {
            std::cout << "[WARNING] Facet construction failed for edge with " << finite_cell_count << " finite cells\n";
            continue;
        }

        // Verify facet validity
        if (facet.vertices_indices.size() < 3)
        {
            std::cout << "[ERROR] Degenerate interior facet detected with " << facet.vertices_indices.size() << " vertices despite " << finite_cell_count << " finite cells\n";
            continue;
        }

        int facetIndex = voronoiDiagram.facets.size() - 1; // Assuming facet was just added
        if (std::find(vc.facet_indices.begin(), vc.facet_indices.end(), facetIndex) == vc.facet_indices.end())
        {
            vc.facet_indices.push_back(facetIndex);
        }
    }
}

//! @brief Retrieves the edge key for a pair of vertices.
static std::pair<int, int> get_edge_key(int u, int v)
{
    return {std::min(u, v), std::max(u, v)};
}

//! @brief Validates that each Facet in the Voronoi Diagram has the correct
/// orientation and normal vector.
/*!
 *
 */
void validate_facet_orientations_and_normals(VoronoiDiagram &voronoiDiagram)
{
    // propagation per cell
    for (auto &cell : voronoiDiagram.cells)
    {
        if (cell.facet_indices.empty())
            continue;

        // Build facet adjacency: vector of maps: index in facet_indices -> {adj_index in facet_indices: shared_edge_key}
        size_t num_facets = cell.facet_indices.size();
        std::vector<std::map<size_t, std::pair<int, int>>> facet_adj(num_facets);
        for (size_t i = 0; i < num_facets; ++i)
        {
            int f1 = cell.facet_indices[i];
            const auto &verts1 = voronoiDiagram.facets[f1].vertices_indices;
            std::map<std::pair<int, int>, size_t> edge_to_pos1;
            for (size_t j = 0; j < verts1.size(); ++j)
            {
                int u = verts1[j];
                int v = verts1[(j + 1) % verts1.size()];
                edge_to_pos1[get_edge_key(u, v)] = j;
            }

            for (size_t k = i + 1; k < num_facets; ++k)
            {
                int f2 = cell.facet_indices[k];
                const auto &verts2 = voronoiDiagram.facets[f2].vertices_indices;
                std::pair<int, int> shared_edge = {-1, -1};
                int shared_count = 0;
                for (size_t j = 0; j < verts2.size(); ++j)
                {
                    int u = verts2[j];
                    int v = verts2[(j + 1) % verts2.size()];
                    auto key = get_edge_key(u, v);
                    if (edge_to_pos1.count(key))
                    {
                        shared_edge = key;
                        shared_count++;
                        if (shared_count > 1)
                            break; // Not adjacent if >1 edge shared
                    }
                }
                if (shared_count == 1)
                {
                    facet_adj[i][k] = shared_edge;
                    facet_adj[k][i] = shared_edge;
                }
            }
        }

        // BFS to propagate
        std::vector<bool> visited(num_facets, false);
        std::queue<size_t> q;
        q.push(0);
        visited[0] = true;

        while (!q.empty())
        {
            size_t curr = q.front();
            q.pop();
            int curr_f = cell.facet_indices[curr];
            auto &curr_verts = voronoiDiagram.facets[curr_f].vertices_indices;

            for (const auto &kv : facet_adj[curr])
            {
                size_t adj = kv.first;
                if (visited[adj])
                    continue;
                visited[adj] = true;
                q.push(adj);

                auto shared_edge = kv.second;
                int adj_f = cell.facet_indices[adj];
                auto &adj_verts = voronoiDiagram.facets[adj_f].vertices_indices;

                // Find direction in curr: true if u to v (min to max)
                bool curr_dir_uv = false;
                for (size_t j = 0; j < curr_verts.size(); ++j)
                {
                    int a = curr_verts[j];
                    int b = curr_verts[(j + 1) % curr_verts.size()];
                    if (get_edge_key(a, b) == shared_edge)
                    {
                        curr_dir_uv = (a == shared_edge.first && b == shared_edge.second);
                        break;
                    }
                }

                // In adj
                bool adj_dir_uv = false;
                for (size_t j = 0; j < adj_verts.size(); ++j)
                {
                    int a = adj_verts[j];
                    int b = adj_verts[(j + 1) % adj_verts.size()];
                    if (get_edge_key(a, b) == shared_edge)
                    {
                        adj_dir_uv = (a == shared_edge.first && b == shared_edge.second);
                        break;
                    }
                }

                // If same direction, reverse adj
                if (curr_dir_uv == adj_dir_uv)
                {
                    std::reverse(adj_verts.begin(), adj_verts.end());
                    // std::cout << "[INFO] Reversed intra-cell facet " << adj_f << " in cell " << cell.cellIndex << " to match opposite edge {" << shared_edge.first << "," << shared_edge.second << "} with facet " << curr_f << "\n";
                }
            }
        }

        int outward_count = 0;
        int total_non_deg = 0;
        for (size_t f = 0; f < num_facets; ++f)
        {
            int fi = cell.facet_indices[f];
            auto &V = voronoiDiagram.facets[fi].vertices_indices;
            if (V.size() < 3)
                continue;

            Point centroid(0, 0, 0);
            for (int idx : V)
            {
                centroid = centroid + (voronoiDiagram.vertices[idx].coord - CGAL::ORIGIN) / V.size();
            }

            Vector3 normal(0, 0, 0);
            size_t n = V.size();
            for (size_t k = 0; k < n; ++k)
            {
                const Point &p1 = voronoiDiagram.vertices[V[k]].coord;
                const Point &p2 = voronoiDiagram.vertices[V[(k + 1) % n]].coord;
                normal = normal + Vector3(
                                      (p1.y() - p2.y()) * (p1.z() + p2.z()),
                                      (p1.z() - p2.z()) * (p1.x() + p2.x()),
                                      (p1.x() - p2.x()) * (p1.y() + p2.y()));
            }
            normal = normal / 2.0;

            double sq_norm = normal.squared_length();
            if (sq_norm > 1e-10)
            { // Non-degenerate
                Vector3 v = cell.delaunay_vertex->point() - centroid;
                double dot = CGAL::scalar_product(normal, v);
                if (dot <= 0)
                    outward_count++;
                total_non_deg++;
            }
        }

        if (total_non_deg > 0 && outward_count < total_non_deg / 2)
        {
            // Majority inward, reverse all facets in the cell
            // std::cout << "[INFO] Reversing all facets in cell " << cell.cellIndex << " to make majority outward (outward_count: " << outward_count << " / " << total_non_deg << ")\n";
            for (size_t f = 0; f < num_facets; ++f)
            {
                int fi = cell.facet_indices[f];
                auto &V = voronoiDiagram.facets[fi].vertices_indices;
                std::reverse(V.begin(), V.end());
            }
        }
    }
}

//! @brief Constructs Voronoi cells without using Convex_Hull_3 (in development).
/*!
 * Populates the Voronoi diagram with polyhedral cells derived from the Delaunay
 * triangulation by processing incident edges and cell circulators.
 *
 * @param voronoiDiagram The Voronoi diagram to populate with cells.
 * @param dt The Delaunay triangulation corresponding (dual) to the Voronoi diagram.
 */
void construct_voronoi_cells_from_delaunay_triangulation(VoronoiDiagram &voronoiDiagram, Delaunay &dt)
{
    std::map<std::pair<int, int>, std::vector<int>> edge_to_facets;
    int cellIndex = 0;

    for (Vertex_handle v : dt.finite_vertex_handles())
    {
        if (v->info().is_dummy)
            continue;

        VoronoiCell vc = create_voronoi_cell(v, cellIndex);
        collcet_cell_vertices(dt, v, voronoiDiagram, vc.vertices_indices);
        process_incident_edges(dt, v, voronoiDiagram, vc, edge_to_facets);

        if (vc.facet_indices.size() < 4)
        {
            std::cout << "[WARNING] Cell " << cellIndex << " has only " << vc.facet_indices.size() << " facets, skipping\n";
        }
        else
        {
            voronoiDiagram.cells.push_back(vc);
            v->info().voronoiCellIndex = cellIndex;
            cellIndex++;
        }
    }

    for (const auto &kv : edge_to_facets)
    {
        const std::vector<int> &dfacets = kv.second;
        if (dfacets.size() == 2)
        {
            int f1 = dfacets[0];
            int f2 = dfacets[1];
            voronoiDiagram.facets[f1].mirror_facet_index = f2;
            voronoiDiagram.facets[f2].mirror_facet_index = f1;

            auto &A = voronoiDiagram.facets[f1].vertices_indices;
            auto &B = voronoiDiagram.facets[f2].vertices_indices;
        }
        else if (dfacets.size() == 1)
        {
            int f = dfacets[0];
            voronoiDiagram.facets[f].mirror_facet_index = -1;
        }
    }
}

// Helper function to check if two directions are approximately equal
bool directions_equal(const Vector3 &d1, const Vector3 &d2, double epsilon)
{
    Vector3 n1 = d1 / std::sqrt(d1.squared_length());      // Normalize d1
    Vector3 n2 = d2 / std::sqrt(d2.squared_length());      // Normalize d2
    return (n1 - n2).squared_length() < epsilon * epsilon; // Compare squared distance of normalized vectors
}

//! @brief Constructs Voronoi edges from Delaunay facets.
void construct_voronoi_edges(VoronoiDiagram &voronoiDiagram, Delaunay &dt)
{
    voronoiDiagram.edges.clear();
    std::map<std::pair<int, int>, int> segmentMap;              // Maps sorted vertex pairs to edge indices
    std::map<int, std::vector<std::pair<Vector3, int>>> rayMap; // Maps vertex index to (direction, edgeIdx) pairs
    const double EPSILON = 1e-6;

    for (Delaunay::Finite_facets_iterator fit = dt.finite_facets_begin(); fit != dt.finite_facets_end(); ++fit)
    {
        Facet facet = *fit;
        CGAL::Object edgeobj = dt.dual(facet);
        Segment3 seg;
        Ray3 ray;

        if (CGAL::assign(seg, edgeobj))
        {
            Cell_handle c1 = facet.first;
            Cell_handle c2 = c1->neighbor(facet.second);
            int idx1 = dt.is_infinite(c1) ? -1 : c1->info().dualVoronoiVertexIndex;
            int idx2 = dt.is_infinite(c2) ? -1 : c2->info().dualVoronoiVertexIndex;

            if (idx1 != -1 && idx2 != -1 && idx1 != idx2)
            { // Finite segment
                int v1 = std::min(idx1, idx2);
                int v2 = std::max(idx1, idx2);
                auto it = segmentMap.find({v1, v2});

                VoronoiEdge vEdge(edgeobj);
                vEdge.type = 0;
                vEdge.vertex1 = v1;
                vEdge.vertex2 = v2;
                int edgeIdx = voronoiDiagram.edges.size();
                vEdge.delaunayFacets.push_back(facet);
                voronoiDiagram.edges.push_back(vEdge);
                segmentMap[{v1, v2}] = edgeIdx;
                voronoiDiagram.segmentVertexPairToEdgeIndex[{v1, v2}] = edgeIdx;
            }
        }
        else if (CGAL::assign(ray, edgeobj))
        {
            Cell_handle c1 = facet.first;
            Cell_handle c2 = c1->neighbor(facet.second);
            int vertex1 = -1;
            if (!dt.is_infinite(c1))
            {
                vertex1 = c1->info().dualVoronoiVertexIndex;
            }
            else if (!dt.is_infinite(c2))
            {
                vertex1 = c2->info().dualVoronoiVertexIndex;
            }
            if (vertex1 != -1)
            {
                Vector3 dir = ray.direction().vector();
                bool found = false;
                auto it = rayMap.find(vertex1);
                if (it != rayMap.end())
                {
                    for (const auto &pair : it->second)
                    {
                        if (directions_equal(pair.first, dir, EPSILON))
                        {
                            voronoiDiagram.edges[pair.second].delaunayFacets.push_back(facet);
                            found = true;
                            break;
                        }
                    }
                }
                if (!found)
                {
                    VoronoiEdge vEdge(edgeobj);
                    vEdge.type = 1;
                    vEdge.vertex1 = vertex1;
                    vEdge.vertex2 = -1;
                    vEdge.source = ray.source();
                    vEdge.direction = dir;
                    int edgeIdx = voronoiDiagram.edges.size();
                    vEdge.delaunayFacets.push_back(facet);
                    voronoiDiagram.edges.push_back(vEdge);
                    rayMap[vertex1].push_back({dir, edgeIdx});
                }
            }
        }
        // Lines are not expected for finite facets
    }
}

//! @brief Builds Voronoi cell edges for each edge in the diagram.
/*!
 * Creates VoronoiCellEdge entries for cells sharing each edge, collecting cell indices
 * from associated Delaunay facets.
 *
 * @param voronoiDiagram The Voronoi diagram to populate with cell edges.
 * @param dt The Delaunay triangulation.
 */
static void build_cell_edges(
    VoronoiDiagram &voronoiDiagram,
    Delaunay &dt)
{
    voronoiDiagram.cellEdges.clear();
    std::vector<std::unordered_set<int>> cellIndicesPerEdge(voronoiDiagram.edges.size());

#pragma omp parallel for
    for (int edgeIdx = 0; edgeIdx < voronoiDiagram.edges.size(); ++edgeIdx)
    {

        const std::vector<Facet> &sharedFacets = voronoiDiagram.edges[edgeIdx].delaunayFacets;

        std::unordered_set<int> cellIndices;

        for (const Facet &f : sharedFacets)
        {
            Cell_handle c = f.first;
            if (dt.is_infinite(c))
                continue;
            int opp = f.second; // Opposite vertex index
            for (int corner = 0; corner < 4; ++corner)
            {
                if (corner == opp)
                    continue; // Skip opposite, add only facet's 3 vertices
                Vertex_handle delaunay_vertex = c->vertex(corner);
                if (!delaunay_vertex->info().is_dummy)
                {
                    int cellIdx = delaunay_vertex->info().voronoiCellIndex;
                    cellIndices.insert(cellIdx);
                }
            }
        }
        cellIndicesPerEdge[edgeIdx] = std::move(cellIndices);
    }

    // Collect all cell edges in a single pass
    for (int edgeIdx = 0; edgeIdx < voronoiDiagram.edges.size(); ++edgeIdx)
    {
        for (int cIdx : cellIndicesPerEdge[edgeIdx])
        {
            VoronoiCellEdge cellEdge;
            cellEdge.cellIndex = cIdx;
            cellEdge.edgeIndex = edgeIdx;
            cellEdge.cycleIndices = {};
            cellEdge.nextCellEdge = -1;
            voronoiDiagram.cellEdges.push_back(cellEdge);
        }
    }
}

//! @brief Links Voronoi cell edges in a circular ring.
/*!
 * Connects cell edges sharing the same edge index using the nextCellEdge field
 * to form a closed loop.
 *
 * @param voronoiDiagram The Voronoi diagram containing cell edges to link.
 */
static void link_cell_edges(
    VoronoiDiagram &voronoiDiagram)
{
    std::unordered_map<int, std::vector<int>> edgeIdx_to_cellEdges;
    for (int ceIdx = 0; ceIdx < (int)voronoiDiagram.cellEdges.size(); ++ceIdx)
    {
        const VoronoiCellEdge &ce = voronoiDiagram.cellEdges[ceIdx];
        edgeIdx_to_cellEdges[ce.edgeIndex].push_back(ceIdx);
    }

    for (auto &kv : edgeIdx_to_cellEdges)
    {
        auto &cellEdgeIndices = kv.second;
        int N = (int)cellEdgeIndices.size();
        for (int i = 0; i < N; i++)
        {
            int ceIdx = cellEdgeIndices[i];
            int nextIdx = cellEdgeIndices[(i + 1) % N];
            voronoiDiagram.cellEdges[ceIdx].nextCellEdge = nextIdx;
        }
    }
}

//! @brief Processes edge mapping for a single Voronoi edge.
/*!
 * Updates the segmentVertexPairToEdgeIndex map for segments, rays, and lines
 * after intersecting with the bounding box.
 *
 * @param voronoiDiagram The Voronoi diagram to update.
 * @param edgeObj The CGAL object representing the edge.
 * @param edgeIdx The index of the edge in the diagram.
 * @param bbox The bounding box for intersection.
 */
static void process_edge_mapping(VoronoiDiagram &voronoiDiagram, VoronoiEdge &edge, int edgeIdx, CGAL::Epick::Iso_cuboid_3 &bbox)
{
    if (edge.type == 0) // Only process finite segments combinatorially
    {
        int idx1 = edge.vertex1;
        int idx2 = edge.vertex2;
        if (idx1 != -1 && idx2 != -1)
        {
            int v1 = std::min(idx1, idx2);
            int v2 = std::max(idx1, idx2);
            voronoiDiagram.segmentVertexPairToEdgeIndex[{v1, v2}] = edgeIdx;
        }
    }
    // Rays and lines are skipped; no mapping for infinite edges
}

//! @brief Updates edge mappings for all Voronoi edges.
/*!
 * Processes all edges to update segmentVertexPairToEdgeIndex and cellEdgeLookup maps.
 *
 * @param voronoiDiagram The Voronoi diagram to update.
 * @param bbox The bounding box for intersection.
 */
static void update_edge_mapping(VoronoiDiagram &voronoiDiagram, CGAL::Epick::Iso_cuboid_3 &bbox)
{
    for (int edgeIdx = 0; edgeIdx < static_cast<int>(voronoiDiagram.edges.size()); ++edgeIdx)
    {
        VoronoiEdge &edge = voronoiDiagram.edges[edgeIdx];
        process_edge_mapping(voronoiDiagram, edge, edgeIdx, bbox);
    }

    voronoiDiagram.cellEdgeLookup.clear();
    for (int ceIdx = 0; ceIdx < static_cast<int>(voronoiDiagram.cellEdges.size()); ++ceIdx)
    {
        const VoronoiCellEdge &ce = voronoiDiagram.cellEdges[ceIdx];
        std::pair<int, int> key = std::make_pair(ce.cellIndex, ce.edgeIndex);
        voronoiDiagram.cellEdgeLookup[key] = ceIdx;
    }
}

//! @brief Constructs the Voronoi cell edges in the Voronoi diagram and links them.
/*!
 * Builds cell edges for each Voronoi edge, links them in a circular ring, and
 * updates edge mappings.
 *
 * @param voronoiDiagram The Voronoi diagram to populate with edges.
 * @param bbox The bounding box used for clipping rays and lines.
 * @param dt The Delaunay triangulation.
 */
void construct_voronoi_cell_edges(
    VoronoiDiagram &voronoiDiagram,
    CGAL::Epick::Iso_cuboid_3 &bbox,
    Delaunay &dt)
{
    voronoiDiagram.cellEdges.clear();

    std::clock_t start = std::clock();

    build_cell_edges(voronoiDiagram, dt);
    std::clock_t check1 = std::clock();
    double duration1 = static_cast<double>(check1 - start) / CLOCKS_PER_SEC;

    std::cout << "build cell edge Execution time: " << duration1 << " seconds" << std::endl;

    link_cell_edges(voronoiDiagram);
    std::clock_t check2 = std::clock();
    double duration2 = static_cast<double>(check2 - check1) / CLOCKS_PER_SEC;

    std::cout << "link cell edge Execution time: " << duration2 << " seconds" << std::endl;

    update_edge_mapping(voronoiDiagram, bbox);

    std::clock_t check3 = std::clock();

    double duration3 = static_cast<double>(check3 - check2) / CLOCKS_PER_SEC;

    std::cout << "update edge mapping Execution time: " << duration3 << " seconds" << std::endl;
}

//! @brief Wrap up function of constructing voronoi diagram
void construct_voronoi_diagram(VoronoiDiagram &vd, VDC_PARAM &vdc_param, UnifiedGrid &grid, CGAL::Epick::Iso_cuboid_3 &bbox, Delaunay &dt)
{
    std::clock_t start = std::clock();
    std::cout << "[INFO] Start constructing Voronoi vertices and edges..." << std::endl;
    construct_voronoi_vertices(vd, dt);
    std::clock_t check0 = std::clock();
    double duration0 = static_cast<double>(check0 - start) / CLOCKS_PER_SEC;
    std::cout << "construct vertices Execution time: " << duration0 << " seconds" << std::endl;

    construct_voronoi_edges(vd, dt);
    std::clock_t check1 = std::clock();
    double duration1 = static_cast<double>(check1 - start) / CLOCKS_PER_SEC;
    std::cout << "construct edges Execution time: " << duration1 << " seconds" << std::endl;

    compute_voronoi_values(vd, grid);
    std::clock_t check2 = std::clock();
    double duration2 = static_cast<double>(check2 - check1) / CLOCKS_PER_SEC;

    if (vdc_param.multi_isov)
    {
        if (vdc_param.convex_hull)
        {
            construct_voronoi_cells_as_convex_hull(vd, dt);
        }
        else
        {
            construct_voronoi_cells_from_delaunay_triangulation(vd, dt);
            std::clock_t check3 = std::clock();
            double duration3 = static_cast<double>(check3 - check2) / CLOCKS_PER_SEC;
            std::cout << "construct cells Execution time: " << duration3 << " seconds" << std::endl;
            validate_facet_orientations_and_normals(vd);
            std::clock_t check4 = std::clock();
            double duration4 = static_cast<double>(check4 - check3) / CLOCKS_PER_SEC;
            std::cout << "validate facet orientations and normals in cells Execution time: " << duration4 << " seconds" << std::endl;
        }

        std::clock_t t = std::clock();

        construct_voronoi_cell_edges(vd, bbox, dt);
        std::clock_t check5 = std::clock();
        double duration5 = static_cast<double>(check5 - t) / CLOCKS_PER_SEC;
        std::cout << "construct cell edges Execution time: " << duration5 << " seconds" << std::endl;

        vd.create_global_facets();
    }

    std::clock_t t2 = std::clock();
    vd.check(false);
    std::clock_t check6 = std::clock();
    double duration6 = static_cast<double>(check6 - t2) / CLOCKS_PER_SEC;
    std::cout << "vd.check() Execution time: " << duration6 << " seconds" << std::endl;
}

// ！@brief Wrap up function for constructing iso surface
void construct_iso_surface(Delaunay &dt, VoronoiDiagram &vd, VDC_PARAM &vdc_param, IsoSurface &iso_surface, UnifiedGrid &grid, std::vector<Point> &activeCubeCenters, CGAL::Epick::Iso_cuboid_3 &bbox)
{
    std::clock_t start_time = std::clock();
    if (vdc_param.multi_isov)
    {
        compute_isosurface_vertices_multi(vd, vdc_param.isovalue, iso_surface);
    }
    else
    {
        compute_isosurface_vertices_single(grid, vdc_param.isovalue, iso_surface, activeCubeCenters);
    }

    std::clock_t check_time = std::clock();
    double duration = static_cast<double>(check_time - start_time) / CLOCKS_PER_SEC;
    std::cout << "Compute isosurface vertices Execution time: " << duration << " seconds" << std::endl;

    if (vdc_param.multi_isov)
    {
        compute_dual_triangles_multi(vd, bbox, grid, vdc_param.isovalue, iso_surface);
    }
    else
    {
        compute_dual_triangles(iso_surface, vd, bbox, dt, grid, vdc_param.isovalue);
    }

    std::clock_t check2_time = std::clock();
    double duration2 = static_cast<double>(check2_time - check_time) / CLOCKS_PER_SEC;
    std::cout << "Compute isosurface facets Execution time: " << duration2 << " seconds" << std::endl;
}

//! @brief Handles output mesh generation.
/*!
 * Writes the final isosurface mesh to file in the specified format (OFF or PLY).
 * Supports both single and multi-isovertex modes.
 *
 * @param retFlag Output parameter indicating success/failure
 * @param vd The Voronoi diagram
 * @param vdc_param Configuration parameters
 * @param iso_surface The isosurface to output
 * @return EXIT_SUCCESS on success, EXIT_FAILURE on error
 */
int handle_output_mesh(bool &retFlag, VoronoiDiagram &vd, VDC_PARAM &vdc_param, IsoSurface &iso_surface)
{
    retFlag = true;

    std::cout << "Result file at: " << vdc_param.output_filename << std::endl;

    // Multi-isovertex mode output
    if (vdc_param.multi_isov)
    {
        if (vdc_param.output_format == "off")
        {
            writeOFFMulti(vdc_param.output_filename, vd, iso_surface.isosurfaceTrianglesMulti, iso_surface);
        }
        else if (vdc_param.output_format == "ply")
        {
            writePLYMulti(vdc_param.output_filename, vd, iso_surface.isosurfaceTrianglesMulti, iso_surface);
        }
        else
        {
            std::cerr << "Unsupported output format: " << vdc_param.output_format << std::endl;
            return EXIT_FAILURE;
        }
    }
    // Single-isovertex mode output
    else
    {
        if (vdc_param.output_format == "off")
        {
            writeOFFSingle(vdc_param.output_filename, iso_surface.isosurfaceVertices, iso_surface.isosurfaceTrianglesSingle);
        }
        else if (vdc_param.output_format == "ply")
        {
            writePLYSingle(vdc_param.output_filename, iso_surface.isosurfaceVertices, iso_surface.isosurfaceTrianglesSingle);
        }
        else
        {
            std::cerr << "Unsupported output format: " << vdc_param.output_format << std::endl;
            return EXIT_FAILURE;
        }
    }

    retFlag = false;
    return EXIT_SUCCESS;
}
