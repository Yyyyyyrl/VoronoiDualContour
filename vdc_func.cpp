//! @file vdc_func.cpp
//! @brief Implementation of functions for Voronoi Diagram and Isosurface computation.

#include "vdc_func.h"

//! @brief Helper function - returns the index of the vertex matching p, or -1 if not found.

//! @brief Implementation of collapse_small_voronoi_edges function
void collapse_small_voronoi_edges(VoronoiDiagram &vd, double D, CGAL::Epick::Iso_cuboid_3 &bbox)
{
    // Simply call the member function
    vd.collapseSmallEdges(D, bbox);
}
int find_vertex_index(const VoronoiDiagram &vd, const Point &p)
{
    for (const auto &vVertex : vd.vertices)
    {
        // Use an appropriate comparison.
        if (PointApproxEqual()(vVertex.vertex, p))
            return vVertex.index;
    }
    return -1; // Not found (should not happen if all points are valid)
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
static void generateTriangle(
    const Point &p1, const Point &p2, const Point &p3,
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
 * @param seg The segment edge to process.
 * @param edge The CGAL object representing the edge.
 * @param vertexValueMap Map of Voronoi vertices to scalar values.
 * @param isovalue The isovalue for bipolarity check.
 * @param voronoi_edge_to_delaunay_facet_map Map linking edges to Delaunay facets.
 * @param dt The Delaunay triangulation.
 * @param dualTriangles Vector to store generated triangles.
 */
static void processSegmentEdge(
    const Segment3 &seg,
    const CGAL::Object &edge,
    std::map<Point, float> &vertexValueMap,
    float isovalue,
    const std::map<Object, std::vector<Facet>, ObjectComparator> &voronoi_edge_to_delaunay_facet_map,
    Delaunay &dt,
    std::vector<DelaunayTriangle> &dualTriangles)
{
    Point v1 = seg.source();
    Point v2 = seg.target();

    if (is_bipolar(vertexValueMap[v1], vertexValueMap[v2], isovalue))
    {
        Point positive = (vertexValueMap[v1] >= vertexValueMap[v2]) ? v1 : v2;

        for (const auto &facet : voronoi_edge_to_delaunay_facet_map.at(edge))
        {
            int iFacet = facet.second;
            Cell_handle c = facet.first;
            int d1 = (iFacet + 1) % 4;
            int d2 = (iFacet + 2) % 4;
            int d3 = (iFacet + 3) % 4;

            Point p1 = c->vertex(d1)->point();
            Point p2 = c->vertex(d2)->point();
            Point p3 = c->vertex(d3)->point();

            int iOrient = get_orientation(iFacet, v1, v2, vertexValueMap[v1], vertexValueMap[v2]);
            generateTriangle(p1, p2, p3, iOrient, dt.is_infinite(c), dualTriangles);
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
 * @param vertexValueMap Map of Voronoi vertices to scalar values.
 * @param bbox The bounding box for intersection.
 * @param grid The scalar grid for interpolation.
 * @param isovalue The isovalue for bipolarity check.
 * @param voronoi_edge_to_delaunay_facet_map Map linking edges to Delaunay facets.
 * @param dt The Delaunay triangulation.
 * @param dualTriangles Vector to store generated triangles.
 */
static void processRayEdge(
    const Ray3 &ray,
    const CGAL::Object &edge,
    std::map<Point, float> &vertexValueMap,
    CGAL::Epick::Iso_cuboid_3 &bbox,
    ScalarGrid &grid,
    float isovalue,
    const std::map<Object, std::vector<Facet>, ObjectComparator> &voronoi_edge_to_delaunay_facet_map,
    Delaunay &dt,
    std::vector<DelaunayTriangle> &dualTriangles)
{
    CGAL::Object intersectObj = CGAL::intersection(bbox, ray);
    Segment3 iseg;
    if (CGAL::assign(iseg, intersectObj))
    {
        Point v1 = iseg.source();
        Point v2 = iseg.target();
        float iPt_value = trilinear_interpolate(adjust_outside_bound_points(v2, grid, v1, v2), grid);

        if (is_bipolar(vertexValueMap[v1], iPt_value, isovalue))
        {
            Point positive = (vertexValueMap[v1] >= iPt_value) ? v1 : v2;

            for (const auto &facet : voronoi_edge_to_delaunay_facet_map.at(edge))
            {
                Facet mirror_f = dt.mirror_facet(facet);
                Object e = dt.dual(facet);

                int iFacet = facet.second;
                Cell_handle c = facet.first;
                int d1 = (iFacet + 1) % 4;
                int d2 = (iFacet + 2) % 4;
                int d3 = (iFacet + 3) % 4;

                Point p1 = c->vertex(d1)->point();
                Point p2 = c->vertex(d2)->point();
                Point p3 = c->vertex(d3)->point();

                int iOrient = get_orientation(iFacet, v1, v2, vertexValueMap[v1], iPt_value);
                generateTriangle(p1, p2, p3, iOrient, dt.is_infinite(c), dualTriangles);
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
 * @param voronoi_edge_to_delaunay_facet_map Map linking edges to Delaunay facets.
 * @param dt The Delaunay triangulation.
 * @param dualTriangles Vector to store generated triangles.
 */
static void processLineEdge(
    const Line3 &line,
    const CGAL::Object &edge,
    ScalarGrid &grid,
    float isovalue,
    CGAL::Epick::Iso_cuboid_3 &bbox,
    const std::map<Object, std::vector<Facet>, ObjectComparator> &voronoi_edge_to_delaunay_facet_map,
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

            for (const auto &facet : voronoi_edge_to_delaunay_facet_map.at(edge))
            {
                int iFacet = facet.second;
                Cell_handle c = facet.first;
                int d1 = (iFacet + 1) % 4;
                int d2 = (iFacet + 2) % 4;
                int d3 = (iFacet + 3) % 4;

                Point p1 = c->vertex(d1)->point();
                Point p2 = c->vertex(d2)->point();
                Point p3 = c->vertex(d3)->point();

                int iOrient = get_orientation(iFacet, intersection1, intersection2, iPt1_val, iPt2_val);
                generateTriangle(p1, p2, p3, iOrient, dt.is_infinite(c), dualTriangles);
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
 * @param vertexValueMap Map of Voronoi vertices to scalar values.
 * @param bbox Bounding box of the computational domain.
 * @param delaunay_facet_to_voronoi_edge_map Map linking Delaunay facets to Voronoi edges.
 * @param dt Delaunay triangulation structure.
 * @param grid Scalar grid containing scalar values.
 * @param isovalue The isovalue used for computing.
 * @param point_index_map The map between points in the Delaunay triangulation to their indices.
 */
void computeDualTriangles(
    IsoSurface &iso_surface,
    std::vector<CGAL::Object> &voronoi_edges,
    std::map<Point, float> &vertexValueMap,
    CGAL::Epick::Iso_cuboid_3 &bbox,
    std::map<Object, std::vector<Facet>, ObjectComparator> &delaunay_facet_to_voronoi_edge_map,
    Delaunay &dt,
    ScalarGrid &grid,
    float isovalue,
    std::map<Point, int> &pointToIndexMap)
{
    std::vector<DelaunayTriangle> dualTriangles;

    for (const auto &edge : voronoi_edges)
    {
        Segment3 seg;
        Ray3 ray;
        Line3 line;

        if (CGAL::assign(seg, edge))
        {
            processSegmentEdge(seg, edge, vertexValueMap, isovalue, delaunay_facet_to_voronoi_edge_map, dt, dualTriangles);
        }
        else if (CGAL::assign(ray, edge))
        {
            processRayEdge(ray, edge, vertexValueMap, bbox, grid, isovalue, delaunay_facet_to_voronoi_edge_map, dt, dualTriangles);
        }
        else if (CGAL::assign(line, edge))
        {
            processLineEdge(line, edge, grid, isovalue, bbox, delaunay_facet_to_voronoi_edge_map, dt, dualTriangles);
        }
    }

    iso_surface.isosurfaceTrianglesSingle = dualTriangles;
}

static inline int selectIsovertexFromCellEdge(
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

    // If no cycleIndex was assigned, access its next cellEdge or return -1 if already traversed through all 3 cells surrounding the vEdge
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

    // The cell in question:
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
static void generateTriangleMulti(
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
    else
    {
        std::cout << "Problematic triangle" << std::endl;
        std::cout << "Vertex 1: " << idx1 << std::endl;
        std::cout << "Vertex 2: " << idx2 << std::endl;
        std::cout << "Vertex 3: " << idx3 << std::endl;
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
static bool selectIsovertices(
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

    idx1 = selectIsovertexFromCellEdge(voronoiDiagram, cellIndex1, globalEdgeIndex);
    idx2 = selectIsovertexFromCellEdge(voronoiDiagram, cellIndex2, globalEdgeIndex);
    idx3 = selectIsovertexFromCellEdge(voronoiDiagram, cellIndex3, globalEdgeIndex);

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
 * @param voronoi_edge_to_delaunay_facet_map Map linking edges to Delaunay facets.
 * @param iso_surface The isosurface to store triangles.
 */
static void processSegmentEdgeMulti(
    const Segment3 &seg,
    const CGAL::Object &edge,
    VoronoiDiagram &voronoiDiagram,
    float isovalue,
    const std::map<CGAL::Object, std::vector<Facet>, ObjectComparator> &voronoi_edge_to_delaunay_facet_map,
    IsoSurface &iso_surface)
{
    Point v1 = seg.source();
    Point v2 = seg.target();
    int idx_v1 = find_vertex_index(voronoiDiagram, v1);
    int idx_v2 = find_vertex_index(voronoiDiagram, v2);
    float val1 = voronoiDiagram.vertexValues[idx_v1];
    float val2 = voronoiDiagram.vertexValues[idx_v2];

    if (is_bipolar(val1, val2, isovalue))
    {
        if (idx_v1 > idx_v2)
            std::swap(idx_v1, idx_v2);
        auto itEdge = voronoiDiagram.segmentVertexPairToEdgeIndex.find(std::make_pair(idx_v1, idx_v2));
        if (itEdge == voronoiDiagram.segmentVertexPairToEdgeIndex.end())
            return;

        int globalEdgeIndex = itEdge->second;
        auto it = voronoi_edge_to_delaunay_facet_map.find(edge);
        if (it != voronoi_edge_to_delaunay_facet_map.end())
        {
            for (const auto &facet : it->second)
            {
                int idx1, idx2, idx3, cellIndex1, cellIndex2, cellIndex3;
                bool isValid = selectIsovertices(voronoiDiagram, facet, globalEdgeIndex, idx1, idx2, idx3, cellIndex1, cellIndex2, cellIndex3);
                int iOrient = get_orientation(facet.second, v1, v2, val1, val2);
                generateTriangleMulti(iso_surface, idx1, idx2, idx3, iOrient, isValid);
            }
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
 * @param voronoi_edge_to_delaunay_facet_map Map linking edges to Delaunay facets.
 * @param iso_surface The isosurface to store triangles.
 */
static void processRayEdgeMulti(
    const Ray3 &ray,
    const CGAL::Object &edge,
    VoronoiDiagram &voronoiDiagram,
    ScalarGrid &grid,
    float isovalue,
    CGAL::Epick::Iso_cuboid_3 &bbox,
    const std::map<CGAL::Object, std::vector<Facet>, ObjectComparator> &voronoi_edge_to_delaunay_facet_map,
    IsoSurface &iso_surface)
{
    CGAL::Object intersectObj = CGAL::intersection(bbox, ray);
    Segment3 iseg;
    if (CGAL::assign(iseg, intersectObj))
    {
        Point v1 = ray.source();
        Point v2 = iseg.target();
        int idx_v1 = find_vertex_index(voronoiDiagram, v1);
        float val1 = voronoiDiagram.vertexValues[idx_v1];
        float val2 = trilinear_interpolate(v2, grid);

        if (is_bipolar(val1, val2, isovalue))
        {
            auto it = voronoi_edge_to_delaunay_facet_map.find(edge);
            if (it != voronoi_edge_to_delaunay_facet_map.end())
            {
                for (const auto &facet : it->second)
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
                    generateTriangleMulti(iso_surface, idx1, idx2, idx3, iOrient, isValid);
                }
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
 * @param voronoi_edge_to_delaunay_facet_map Map linking edges to Delaunay facets.
 * @param iso_surface The isosurface to store triangles.
 */
static void processLineEdgeMulti(
    const Line3 &line,
    const CGAL::Object &edge,
    VoronoiDiagram &voronoiDiagram,
    ScalarGrid &grid,
    float isovalue,
    CGAL::Epick::Iso_cuboid_3 &bbox,
    const std::map<CGAL::Object, std::vector<Facet>, ObjectComparator> &voronoi_edge_to_delaunay_facet_map,
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
            auto it = voronoi_edge_to_delaunay_facet_map.find(edge);
            if (it != voronoi_edge_to_delaunay_facet_map.end())
            {
                for (const auto &facet : it->second)
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
                    generateTriangleMulti(iso_surface, idx1, idx2, idx3, iOrient, isValid);
                }
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
 * @param delaunay_facet_to_voronoi_edge_map Map linking Delaunay facets to Voronoi edges.
 * @param grid Scalar grid containing scalar values.
 * @param isovalue The isovalue for mesh computation.
 * @param iso_surface Instance of IsoSurface containing the isosurface vertices and faces.
 */
void computeDualTrianglesMulti(
    VoronoiDiagram &voronoiDiagram,
    CGAL::Epick::Iso_cuboid_3 &bbox,
    std::map<CGAL::Object, std::vector<Facet>, ObjectComparator> &delaunay_facet_to_voronoi_edge_map,
    ScalarGrid &grid,
    float isovalue,
    IsoSurface &iso_surface)
{
    for (const auto &edge : voronoiDiagram.edges)
    {
        Segment3 seg;
        Ray3 ray;
        Line3 line;

        if (CGAL::assign(seg, edge))
        {
            processSegmentEdgeMulti(seg, edge, voronoiDiagram, isovalue, delaunay_facet_to_voronoi_edge_map, iso_surface);
        }
        else if (CGAL::assign(ray, edge))
        {
            processRayEdgeMulti(ray, edge, voronoiDiagram, grid, isovalue, bbox, delaunay_facet_to_voronoi_edge_map, iso_surface);
        }
        else if (CGAL::assign(line, edge))
        {
            processLineEdgeMulti(line, edge, voronoiDiagram, grid, isovalue, bbox, delaunay_facet_to_voronoi_edge_map, iso_surface);
        }
    }
}

//! @brief Computes isosurface vertices for the single-isovertex case.
void Compute_Isosurface_Vertices_Single(ScalarGrid &grid, float isovalue, IsoSurface &iso_surface, Grid &data_grid, std::vector<Point> &activeCubeCenters)
{
    const int cubeVertices[8][3] = {
        {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};

    const int cubeEdges[12][2] = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6}, {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};

    for (const auto &center : activeCubeCenters)
    {
        std::vector<Point> intersectionPoints;
        float cubeSize = grid.dx; // Assuming grid.dx is the cube size
        std::array<float, 8> scalarValues;

        // Compute scalar values at cube vertices
        for (int i = 0; i < 8; i++)
        {
            Point vertex(
                center.x() + (cubeVertices[i][0] - 0.5f) * cubeSize,
                center.y() + (cubeVertices[i][1] - 0.5f) * cubeSize,
                center.z() + (cubeVertices[i][2] - 0.5f) * cubeSize);
            scalarValues[i] = grid.get_scalar_value_at_point(vertex);
        }

        // Check each edge for intersection with the isovalue
        for (const auto &edge : cubeEdges)
        {
            int idx1 = edge[0];
            int idx2 = edge[1];
            float val1 = scalarValues[idx1];
            float val2 = scalarValues[idx2];

            if (is_bipolar(val1, val2, isovalue))
            {
                Point p1(
                    center.x() + (cubeVertices[idx1][0] - 0.5f) * cubeSize,
                    center.y() + (cubeVertices[idx1][1] - 0.5f) * cubeSize,
                    center.z() + (cubeVertices[idx1][2] - 0.5f) * cubeSize);

                Point p2(
                    center.x() + (cubeVertices[idx2][0] - 0.5f) * cubeSize,
                    center.y() + (cubeVertices[idx2][1] - 0.5f) * cubeSize,
                    center.z() + (cubeVertices[idx2][2] - 0.5f) * cubeSize);

                Point intersect = interpolate(p1, p2, val1, val2, isovalue, data_grid);
                intersectionPoints.push_back(intersect);
            }
        }

        // Compute the centroid of the intersection points
        if (!intersectionPoints.empty())
        {
            Point centroid = compute_centroid(intersectionPoints);
            iso_surface.isosurfaceVertices.push_back(centroid);
        }
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
static void collectMidpoints(
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
        VoronoiFacet &facet = voronoiDiagram.facets[facet_index];
        size_t num_vertices = facet.vertices_indices.size();

        std::vector<int> current_facet_midpoints;

        for (size_t j = 0; j < num_vertices; ++j)
        {
            size_t idx1 = j;
            size_t idx2 = (j + 1) % num_vertices;

            float val1 = facet.vertex_values[idx1];
            float val2 = facet.vertex_values[idx2];

            if (is_bipolar(val1, val2, isovalue))
            {
                int vertex_index1 = facet.vertices_indices[idx1];
                int vertex_index2 = facet.vertices_indices[idx2];

                Point p1 = voronoiDiagram.vertices[vertex_index1].vertex;
                Point p2 = voronoiDiagram.vertices[vertex_index2].vertex;

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

//! @brief Connects midpoints within each facet to form a graph.
/*!
 * Links pairs of midpoints in each facet to establish connectivity for cycle detection.
 *
 * @param facet_midpoint_indices Vector storing midpoint indices for each facet.
 * @param midpoints Vector of midpoints to update with connectivity.
 */
static void connectMidpoints(
    const std::vector<std::vector<int>> &facet_midpoint_indices,
    std::vector<MidpointNode> &midpoints)
{
    for (const auto &facet_midpoints : facet_midpoint_indices)
    {
        size_t num_midpoints = facet_midpoints.size();
        for (size_t k = 0; k + 1 < num_midpoints; k += 2)
        {
            int idx1 = facet_midpoints[k];
            int idx2 = facet_midpoints[k + 1];
            midpoints[idx1].connected_to.push_back(idx2);
            midpoints[idx2].connected_to.push_back(idx1);
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
static void extractCycles(
    const std::vector<MidpointNode> &midpoints,
    std::vector<std::vector<int>> &cycles)
{
    std::set<int> visited;

    for (size_t i = 0; i < midpoints.size(); ++i)
    {
        if (visited.find(i) == visited.end())
        {
            std::vector<int> cycle;
            std::stack<int> stack;
            stack.push(i);

            while (!stack.empty())
            {
                int current = stack.top();
                stack.pop();

                if (visited.find(current) != visited.end())
                {
                    continue;
                }

                visited.insert(current);
                cycle.push_back(current);

                for (int neighbor : midpoints[current].connected_to)
                {
                    if (visited.find(neighbor) == visited.end())
                    {
                        stack.push(neighbor);
                    }
                }
            }

            if (!cycle.empty())
            {
                cycles.push_back(cycle);
            }
        }
    }
}

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
static void computeCycleCentroids(
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
void Compute_Isosurface_Vertices_Multi(VoronoiDiagram &voronoiDiagram, float isovalue, IsoSurface &iso_surface)
{
    for (auto &vc : voronoiDiagram.cells)
    {
        std::vector<MidpointNode> midpoints;
        std::map<std::pair<int, int>, int> edge_to_midpoint_index;
        std::vector<std::vector<int>> facet_midpoint_indices;

        collectMidpoints(vc, voronoiDiagram, isovalue, midpoints, edge_to_midpoint_index, facet_midpoint_indices);
        connectMidpoints(facet_midpoint_indices, midpoints);

        std::vector<std::vector<int>> cycles;
        extractCycles(midpoints, cycles);

        computeCycleCentroids(vc, voronoiDiagram, midpoints, cycles, iso_surface);
    }
}

//! @brief Adds dummy points from a facet for Voronoi diagram bounding.
std::vector<Point> add_dummy_from_facet(const GRID_FACETS &facet, const Grid &data_grid)
{
    std::vector<Point> points;

    // 2D slice dimension
    int dim0 = facet.axis_size[0];
    int dim1 = facet.axis_size[1];

    // For convenience
    int d = facet.orth_dir;
    int d1 = facet.axis_dir[0];
    int d2 = facet.axis_dir[1];

    // We have bounding-box in facet.minIndex[], facet.maxIndex[],
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
 * @param dummy_points Output vector for dummy points.
 */
static void collectDelaunayPoints(
    Grid &grid,
    const std::vector<std::vector<GRID_FACETS>> &grid_facets,
    const std::vector<Point> &activeCubeCenters,
    VDC_PARAM &vdc_param,
    std::vector<Point> &delaunay_points,
    std::vector<Point> &dummy_points)
{
    if (vdc_param.multi_isov)
    {
        for (const auto &p : activeCubeCenters)
        {
            delaunay_points.push_back(p);
        }

        for (int d = 0; d < 3; d++)
        {
            for (const auto &f : grid_facets[d])
            {
                std::vector<Point> pointsf = add_dummy_from_facet(f, grid);
                dummy_points.insert(dummy_points.end(), pointsf.begin(), pointsf.end());
            }
        }

        delaunay_points.insert(delaunay_points.end(), dummy_points.begin(), dummy_points.end());
    }
    else
    {
        delaunay_points = activeCubeCenters;
    }
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
static void insertPointsIntoTriangulation(
    Delaunay &dt,
    std::vector<Point> &delaunay_points,
    std::vector<Point> &activeCubeCenters,
    VDC_PARAM &vdc_param)
{
    dt.insert(delaunay_points.begin(), delaunay_points.end());

    if (vdc_param.multi_isov)
    {
        write_triangulation(dt, activeCubeCenters, vdc_param.output_filename);
    }
}

//! @brief Sets vertex info for the Delaunay triangulation.
/*!
 * Marks vertices as dummy or regular based on their presence in the dummy points list.
 *
 * @param dt The Delaunay triangulation to update.
 * @param dummy_points The list of dummy points.
 */
static void setVertexInfo(Delaunay &dt, const std::vector<Point> &dummy_points)
{
    for (auto delaunay_vertex = dt.finite_vertices_begin(); delaunay_vertex != dt.finite_vertices_end(); ++delaunay_vertex)
    {
        Point p = delaunay_vertex->point();
        delaunay_vertex->info().is_dummy = (std::find(dummy_points.begin(), dummy_points.end(), p) != dummy_points.end());
    }
}

//! @brief Assigns indices to points in the point index map.
/*!
 * Builds the point index map for multi-isovertex or single-isovertex mode.
 *
 * @param dt The Delaunay triangulation.
 * @param activeCubeCenters The list of center points of active cubes.
 * @param vdc_param The VDC_PARAM instance containing user input options.
 * @param point_index_map The map to store point-to-index mappings.
 */
static void assignPointIndices(
    Delaunay &dt,
    const std::vector<Point> &activeCubeCenters,
    VDC_PARAM &vdc_param,
    std::map<Point, int> &point_index_map)
{
    int i = 0;
    if (vdc_param.multi_isov)
    {
        for (auto delaunay_vertex = dt.finite_vertices_begin(); delaunay_vertex != dt.finite_vertices_end(); ++delaunay_vertex)
        {
            if (delaunay_vertex->info().is_dummy)
            {
                continue;
            }
            Point p = delaunay_vertex->point();
            point_index_map[p] = i++;
        }
    }
    else
    {
        for (const auto &pt : activeCubeCenters)
        {
            point_index_map[pt] = i++;
        }
    }
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
 * @param point_index_map The map between points in the Delaunay triangulation to their indices (used in single isovertex case only).
 */
void construct_delaunay_triangulation(
    Delaunay &dt,
    Grid &grid,
    const std::vector<std::vector<GRID_FACETS>> &grid_facets,
    VDC_PARAM &vdc_param,
    std::vector<Point> &activeCubeCenters,
    std::map<Point, int> &point_index_map)
{
    std::vector<Point> delaunay_points;
    std::vector<Point> dummy_points;

    collectDelaunayPoints(grid, grid_facets, activeCubeCenters, vdc_param, delaunay_points, dummy_points);

    if (vdc_param.multi_isov && debug)
    {
        write_dummy_points(grid, dummy_points);
    }

    insertPointsIntoTriangulation(dt, delaunay_points, activeCubeCenters, vdc_param);

    if (vdc_param.multi_isov)
    {
        setVertexInfo(dt, dummy_points);
    }

    assign_indices_dt(dt, point_index_map);
    assignPointIndices(dt, activeCubeCenters, vdc_param, point_index_map);
}

/**
 * @brief Assigns unique indices to vertices and tetrahedra in a Delaunay Triangulation (DT)
 *
 * This function assigns unique indices to all vertices and tetrahedra in the given
 * Delaunay triangulation and maps these indices to the specified mapping.
 *
 * @param dt The Delaunay triangulation object
 * @param pointToIndexMap Mapping from points to indices, where keys are points and values are their corresponding indices
 */
void assign_indices_dt(Delaunay &dt, std::map<Point, int> &pointToIndexMap)
{
    // First, assign unique indices to all vertices
    int vertex_counter = 0;
    for (auto vertex = dt.finite_vertices_begin(); vertex != dt.finite_vertices_end(); ++vertex)
    {
        vertex->info().index = vertex_counter++;
        vertex->info().voronoiCellIndex = -1; // Initialize voronoi cell index
        pointToIndexMap[vertex->point()] = vertex->info().index;
    }

    // Then, assign unique indices to all cells (tetrahedra)
    int cell_counter = 0;
    for (auto cell = dt.finite_cells_begin(); cell != dt.finite_cells_end(); ++cell)
    {
        cell->info().index = cell_counter++;
        cell->info().dualVoronoiVertexIndex = -1; // Initialize dual voronoi vertex index
    }

    std::cout << " Stat of DT: # of vertices: " << vertex_counter << " # of cells: " << cell_counter << std::endl;
}

//! @brief Constructs Voronoi vertices for the given voronoi Diagram instance.
void construct_voronoi_vertices(VoronoiDiagram &voronoiDiagram, Delaunay &dt)
{
    voronoiDiagram.vertices.clear();
    int vertex_index = 0;
    for (Delaunay::Finite_cells_iterator cit = dt.finite_cells_begin();
         cit != dt.finite_cells_end(); ++cit)
    {
        Point voronoi_vertex = dt.dual(cit);
        VoronoiVertex vVertex(voronoi_vertex);
        vVertex.index = vertex_index;
        voronoiDiagram.vertices.push_back(vVertex);
        cit->info().dualVoronoiVertexIndex = vertex_index;
        vertex_index++;
    }
}

//! @brief Computes Voronoi Vertex values using scalar grid interpolation
void compute_voronoi_values(VoronoiDiagram &voronoiDiagram, ScalarGrid &grid, std::map<Point, float> &vertexValueMap)
{
    voronoiDiagram.vertexValues.resize(voronoiDiagram.vertices.size());
    for (size_t i = 0; i < voronoiDiagram.vertices.size(); ++i)
    {
        Point vertex = voronoiDiagram.vertices[i].vertex;
        float value = trilinear_interpolate(vertex, grid);
        voronoiDiagram.vertexValues[i] = value;
        vertexValueMap[vertex] = value;
    }
}

//! @brief Constructs Voronoi cells from the Delaunay triangulation.
void construct_voronoi_cells(VoronoiDiagram &voronoiDiagram, Delaunay &dt)
{
    int index = 0;
    for (auto delaunay_vertex = dt.finite_vertices_begin(); delaunay_vertex != dt.finite_vertices_end(); ++delaunay_vertex)
    {
        if (delaunay_vertex->info().is_dummy)
        {
            // std::cout << "Dummy Point excluded: " << delaunay_vertex->point() << std::endl;
            continue;
        }
        VoronoiCell vc(delaunay_vertex);
        vc.cellIndex = index;

        std::vector<Cell_handle> incident_cells;
        dt.finite_incident_cells(delaunay_vertex, std::back_inserter(incident_cells));

        // Collect vertex indices, ensuring uniqueness
        std::set<int> unique_vertex_indices_set;
        for (Cell_handle ch : incident_cells)
        {
            if (dt.is_infinite(ch))
            {
                // Through an error, should not be happening after checking dummy vertices
                continue; // Skip infinite cells
            }
            Point voronoi_vertex = dt.dual(ch);

            // Check if voronoi_vertex is within domain and exclude dummy points in the dt
            int vertex_index = find_vertex_index(voronoiDiagram, voronoi_vertex);
            unique_vertex_indices_set.insert(vertex_index);
        }

        // Copy unique indices to vector
        vc.vertices_indices.assign(unique_vertex_indices_set.begin(), unique_vertex_indices_set.end());

        // Build convex hull and extract facets
        std::vector<Point> vertex_points;
        for (int idx : vc.vertices_indices)
        {
            vertex_points.push_back(voronoiDiagram.vertices[idx].vertex);
        }

        // Remove duplicate points
        std::sort(vertex_points.begin(), vertex_points.end(), [](const Point &a, const Point &b)
                  { return a.x() < b.x() || (a.x() == b.x() && (a.y() < b.y() || (a.y() == b.y() && a.z() < b.z()))); });
        vertex_points.erase(std::unique(vertex_points.begin(), vertex_points.end(), PointApproxEqual()), vertex_points.end());

        CGAL::convex_hull_3(vertex_points.begin(), vertex_points.end(), vc.polyhedron);

        // Extract facets from polyhedron
        for (auto facet_it = vc.polyhedron.facets_begin();
             facet_it != vc.polyhedron.facets_end(); ++facet_it)
        {
            VoronoiFacet vf;
            auto h = facet_it->facet_begin();
            do
            {
                Point p = h->vertex()->point();
                int vertex_index = find_vertex_index(voronoiDiagram, p);
                vf.vertices_indices.push_back(vertex_index);
                float value = voronoiDiagram.vertexValues[vertex_index];
                vf.vertex_values.push_back(value);
                ++h;
            } while (h != facet_it->facet_begin());

            int facet_index = voronoiDiagram.facets.size();
            voronoiDiagram.facets.push_back(vf);
            vc.facet_indices.push_back(facet_index);
        }

        voronoiDiagram.cells.push_back(vc);
        delaunay_vertex->info().voronoiCellIndex = index;
        index++;
    }
}

//
// Helper function: Order a set of circumcenters (given by indices) in cyclic order,
// for the facet dual to the Delaunay edge between p0 and p1.
//
void orderFacetVertices(std::vector<int> &indices,
                        const Point &p0,
                        const Point &p1,
                        const std::vector<VoronoiVertex> &vertices)
{
    // Compute the Delaunay edge direction.
    Vector3 edgeDir = p1 - p0;
    double norm = std::sqrt(edgeDir.squared_length());
    if (norm < 1e-10)
        return; // Degenerate edge.
    Vector3 edgeDirUnit = edgeDir / norm;

    // Choose an arbitrary vector not parallel to edgeDirUnit.
    Vector3 arbitrary;
    if (std::abs(CGAL::scalar_product(edgeDirUnit, Vector3(1, 0, 0))) < 0.9)
        arbitrary = Vector3(1, 0, 0);
    else
        arbitrary = Vector3(0, 1, 0);

    // First basis vector: perpendicular to edgeDirUnit.
    Vector3 vRef = CGAL::cross_product(edgeDirUnit, arbitrary);
    double vRefNorm = std::sqrt(vRef.squared_length());
    if (vRefNorm < 1e-10)
    {
        arbitrary = Vector3(0, 0, 1);
        vRef = CGAL::cross_product(edgeDirUnit, arbitrary);
        vRefNorm = std::sqrt(vRef.squared_length());
        if (vRefNorm < 1e-10)
            return;
    }
    vRef = vRef / vRefNorm;

    // Second basis vector.
    Vector3 vRef2 = CGAL::cross_product(edgeDirUnit, vRef);
    double vRef2Norm = std::sqrt(vRef2.squared_length());
    if (vRef2Norm >= 1e-10)
        vRef2 = vRef2 / vRef2Norm;

    // Compute the centroid of the circumcenters.
    double sumX = 0, sumY = 0, sumZ = 0;
    for (int idx : indices)
    {
        const Point &pt = vertices[idx].vertex;
        sumX += pt.x();
        sumY += pt.y();
        sumZ += pt.z();
    }
    double count = static_cast<double>(indices.size());
    Point center(sumX / count, sumY / count, sumZ / count);

    // Sort indices by their angle (using atan2) in the (vRef, vRef2) coordinate system.
    std::sort(indices.begin(), indices.end(), [&](int a, int b)
              {
        const Point &pa = vertices[a].vertex;
        const Point &pb = vertices[b].vertex;

        Vector3 va = pa - center;
        Vector3 vb = pb - center;

        // Remove the component along the edge direction.
        double dotA = CGAL::scalar_product(va, edgeDirUnit);
        Vector3 va_proj = va - dotA * edgeDirUnit;
        double dotB = CGAL::scalar_product(vb, edgeDirUnit);
        Vector3 vb_proj = vb - dotB * edgeDirUnit;

        double angleA = std::atan2(CGAL::scalar_product(va_proj, vRef2),
                                   CGAL::scalar_product(va_proj, vRef));
        double angleB = std::atan2(CGAL::scalar_product(vb_proj, vRef2),
                                   CGAL::scalar_product(vb_proj, vRef));
        return angleA < angleB; });
}

//! @brief Creates a Voronoi cell for a Delaunay vertex.
/*!
 * Initializes a Voronoi cell with the given cell index and Delaunay vertex handle.
 *
 * @param delaunay_vertex The Delaunay vertex to create the cell for.
 * @param cellIndex The index to assign to the cell.
 * @return The initialized Voronoi cell.
 */
static VoronoiCell createVoronoiCell(Vertex_handle delaunay_vertex, int cellIndex)
{
    VoronoiCell vc(delaunay_vertex);
    vc.cellIndex = cellIndex;
    std::cout << "[DEBUG] Created VoronoiCell with index " << cellIndex << "\n";
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
static void collectCellVertices(
    Delaunay &dt,
    Vertex_handle delaunay_vertex,
    VoronoiDiagram &voronoiDiagram,
    std::vector<int> &vertices_indices)
{
    std::vector<Cell_handle> incidentCells;
    dt.finite_incident_cells(delaunay_vertex, std::back_inserter(incidentCells));
    std::cout << "[DEBUG] Found " << incidentCells.size() << " incident cells\n";

    std::set<int> uniqueVertexIndices;
    for (Cell_handle c : incidentCells)
    {
        int oldIdx = c->info().dualVoronoiVertexIndex;
        if (oldIdx >= 0 && oldIdx < voronoiDiagram.oldToNewVertexIndex.size())
        {
            int newIdx = voronoiDiagram.oldToNewVertexIndex[oldIdx];
            if (newIdx >= 0 && newIdx < voronoiDiagram.vertices.size())
            {
                uniqueVertexIndices.insert(newIdx);
            }
            else
            {
                std::cerr << "[ERROR] Invalid new vertex index " << newIdx
                          << " (vertices size: " << voronoiDiagram.vertices.size() << ")\n";
            }
        }
        else
        {
            std::cerr << "[ERROR] Invalid old vertex index " << oldIdx
                      << " (oldToNewVertexIndex size: " << voronoiDiagram.oldToNewVertexIndex.size() << ")\n";
        }
    }
    vertices_indices.assign(uniqueVertexIndices.begin(), uniqueVertexIndices.end());
    std::cout << "[DEBUG] Collected " << vertices_indices.size() << " unique vertex indices\n";
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
static VoronoiFacet buildFacetFromEdge(
    Delaunay &dt,
    const Edge &ed,
    Vertex_handle delaunay_vertex,
    VoronoiDiagram &voronoiDiagram,
    std::vector<int> &facet_indices,
    std::map<std::pair<int, int>, std::vector<int>> &edge_to_facets)
{
    Cell_handle cell_ed = ed.first;
    int i = ed.second;
    int j = ed.third;
    Vertex_handle v1 = cell_ed->vertex(i);
    Vertex_handle v2 = cell_ed->vertex(j);
    if (v1 != delaunay_vertex)
        std::swap(v1, v2); // Ensure v1 is the current vertex

    Delaunay::Cell_circulator cc = dt.incident_cells(ed);
    Delaunay::Cell_circulator start = cc;
    std::vector<int> facetVertexIndices;
    int cell_count = 0;
    do
    {
        if (!dt.is_infinite(cc))
        {
            int oldIdx = cc->info().dualVoronoiVertexIndex;
            int newIdx = voronoiDiagram.oldToNewVertexIndex[oldIdx];
            facetVertexIndices.push_back(newIdx);
            cell_count++;
        }
        ++cc;
    } while (cc != start);

    std::cerr << "[DEBUG] Edge has " << cell_count << " finite incident cells\n";
    std::set<int> unique_vertices(facetVertexIndices.begin(), facetVertexIndices.end());
    std::cerr << "[DEBUG] Facet has " << unique_vertices.size() << " unique vertices\n";

    // Remove consecutive duplicates while preserving order
    std::vector<int> uniqueFacetVertices;
    for (size_t k = 0; k < facetVertexIndices.size(); ++k)
    {
        if (k == 0 || facetVertexIndices[k] != facetVertexIndices[k - 1])
        {
            uniqueFacetVertices.push_back(facetVertexIndices[k]);
        }
    }
    if (uniqueFacetVertices.size() > 1 && uniqueFacetVertices.front() == uniqueFacetVertices.back())
    {
        uniqueFacetVertices.pop_back();
    }

    if (uniqueFacetVertices.size() >= 3)
    {
        // Compute normal and adjust orientation
        Point p0 = voronoiDiagram.vertices[uniqueFacetVertices[0]].vertex;
        Point p1 = voronoiDiagram.vertices[uniqueFacetVertices[1]].vertex;
        Point p2 = voronoiDiagram.vertices[uniqueFacetVertices[2]].vertex;
        Vector3 normal = CGAL::cross_product(p1 - p0, p2 - p0);
        Point centroid(0, 0, 0);
        for (int idx : uniqueFacetVertices)
        {
            centroid = centroid + (voronoiDiagram.vertices[idx].vertex - CGAL::ORIGIN);
        }
        centroid = CGAL::ORIGIN + (centroid - CGAL::ORIGIN) / uniqueFacetVertices.size();
        Point cell_center = delaunay_vertex->point();
        Vector3 v = centroid - cell_center;
        if (CGAL::scalar_product(normal, v) < 0)
        {
            std::reverse(uniqueFacetVertices.begin(), uniqueFacetVertices.end());
        }

        VoronoiFacet facet;
        facet.vertices_indices = uniqueFacetVertices;
        for (int idx : uniqueFacetVertices)
        {
            facet.vertex_values.push_back(voronoiDiagram.vertexValues[idx]);
        }
        int facetIndex = voronoiDiagram.facets.size();
        voronoiDiagram.facets.push_back(facet);
        facet_indices.push_back(facetIndex);
        std::pair<int, int> edge_key = std::make_pair(
            std::min(v1->info().index, v2->info().index),
            std::max(v1->info().index, v2->info().index));
        edge_to_facets[edge_key].push_back(facetIndex);
        return facet;
    }
    else
    {
        std::cerr << "[WARNING] Facet has only " << uniqueFacetVertices.size() << " unique vertices, skipping\n";
        return VoronoiFacet();
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
static void processIncidentEdges(
    Delaunay &dt,
    Vertex_handle delaunay_vertex,
    VoronoiDiagram &voronoiDiagram,
    VoronoiCell &vc,
    std::map<std::pair<int, int>, std::vector<int>> &edge_to_facets)
{
    std::vector<Edge> incidentEdges;
    dt.incident_edges(delaunay_vertex, std::back_inserter(incidentEdges));
    std::cerr << "[DEBUG] Processing " << incidentEdges.size() << " incident edges for cell " << vc.cellIndex << "\n";

    for (const Edge &ed : incidentEdges)
    {
        VoronoiFacet facet = buildFacetFromEdge(dt, ed, delaunay_vertex, voronoiDiagram, vc.facet_indices, edge_to_facets);
        if (!facet.vertices_indices.empty())
        {
            int facetIndex = voronoiDiagram.facets.size() - 1; // Assuming facet was just added
            if (std::find(vc.facet_indices.begin(), vc.facet_indices.end(), facetIndex) == vc.facet_indices.end())
            {
                vc.facet_indices.push_back(facetIndex);
                std::cerr << "[DEBUG] Added facet " << facetIndex << " with " << facet.vertices_indices.size() << " vertices to cell " << vc.cellIndex << "\n";
            }
        }
    }
    std::cerr << "[DEBUG] Cell " << vc.cellIndex << " has " << vc.facet_indices.size() << " facets\n";
}

//! @brief Constructs Voronoi cells without using Convex_Hull_3 (in development).
/*!
 * Populates the Voronoi diagram with polyhedral cells derived from the Delaunay
 * triangulation by processing incident edges and cell circulators.
 *
 * @param voronoiDiagram The Voronoi diagram to populate with cells.
 * @param dt The Delaunay triangulation corresponding (dual) to the Voronoi diagram.
 */
void construct_voronoi_cells_non_convex_hull(VoronoiDiagram &voronoiDiagram, Delaunay &dt)
{
    std::map<std::pair<int, int>, std::vector<int>> edge_to_facets;
    int cellIndex = 0;

    for (Vertex_handle v : dt.finite_vertex_handles())
    {
        if (v->info().is_dummy)
            continue; // Skip dummy vertices

        VoronoiCell vc = createVoronoiCell(v, cellIndex);
        collectCellVertices(dt, v, voronoiDiagram, vc.vertices_indices);
        processIncidentEdges(dt, v, voronoiDiagram, vc, edge_to_facets);

        // Validate the number of facets
        if (vc.facet_indices.size() < 4)
        {
            std::cerr << "[WARNING] Cell " << cellIndex << " has only " << vc.facet_indices.size() << " facets, skipping\n";
            std::cerr << "[DEBUG] Facet indices: ";
            for (int fi : vc.facet_indices)
            {
                std::cerr << fi << " (" << voronoiDiagram.facets[fi].vertices_indices.size() << " vertices) ";
            }
            std::cerr << "\n";
        }
        else
        {
            voronoiDiagram.cells.push_back(vc);
            v->info().voronoiCellIndex = cellIndex;
            std::cerr << "[INFO] Cell " << cellIndex << " constructed with " << vc.facet_indices.size() << " facets\n";
            cellIndex++;
        }
    }

    // Link mirror facets
    for (const auto &kv : edge_to_facets)
    {
        const std::vector<int> &facets = kv.second;
        if (facets.size() == 2)
        {
            int f1 = facets[0];
            int f2 = facets[1];
            voronoiDiagram.facets[f1].mirror_facet_index = f2;
            voronoiDiagram.facets[f2].mirror_facet_index = f1;

            // Verify opposite orientations
            const auto &v1 = voronoiDiagram.facets[f1].vertices_indices;
            const auto &v2 = voronoiDiagram.facets[f2].vertices_indices;
            std::vector<int> v2_rev(v2.rbegin(), v2.rend());
            if (v1 != v2_rev)
            {
                std::cerr << "[WARNING] Facets " << f1 << " and " << f2 << " do not have opposite orientations\n";
            }
        }
        else if (facets.size() == 1)
        {
            int f = facets[0];
            voronoiDiagram.facets[f].mirror_facet_index = -1; // Boundary facet
            std::cerr << "[INFO] Facet " << f << " identified as a boundary facet\n";
        }
        else
        {
            std::cerr << "[ERROR] Edge has " << facets.size() << " facets, expected 1 or 2\n";
        }
    }

    std::cout << "[INFO] Constructed " << voronoiDiagram.cells.size() << " finite Voronoi cells\n";
}

//! @brief Constructs Voronoi cells from the Delaunay triangulation without using convex hull computation.
void construct_voronoi_cells_halfspace(VoronoiDiagram &voronoiDiagram, Delaunay &dt)
{
    int index = 0;
    // Iterate through all finite Delaunay vertices.
    for (auto delaunay_vertex = dt.finite_vertices_begin(); delaunay_vertex != dt.finite_vertices_end(); ++delaunay_vertex)
    {
        // Exclude dummy points.
        if (delaunay_vertex->info().is_dummy)
        {
            continue;
        }
        VoronoiCell vc(delaunay_vertex);
        vc.cellIndex = index;

        // For each Delaunay vertex, the Voronoi cell is the intersection of
        // the half-spaces defined by the perpendicular bisectors between delaunay_vertex and each neighbor.
        std::vector<Vertex_handle> neighbors;
        dt.finite_adjacent_vertices(delaunay_vertex, std::back_inserter(neighbors));

        std::vector<Plane_3> halfspaces;
        for (auto neighbor : neighbors)
        {
            // Get the two Delaunay points, their midpoint and normal
            Point p = delaunay_vertex->point();
            Point q = neighbor->point();
            Point mid = CGAL::midpoint(p, q);
            Vector3 normal = q - p;

            // Construct the bisector plane.
            // With this construction, delaunay_vertex->point() lies in the negative half-space:
            //   (q - p) · (delaunay_vertex->point() - mid) < 0,
            // so the desired half-space is: (q - p) · (x - mid) <= 0.
            Plane_3 bisector(mid, normal);
            halfspaces.push_back(bisector);
        }

        // Use delaunay_vertex->point() as an interior point of the Voronoi cell.
        Point interior = delaunay_vertex->point();
        Polyhedron_3 poly;
        // Compute the intersection of the halfspaces.
        CGAL::halfspace_intersection_3(halfspaces.begin(), halfspaces.end(), poly, interior);
        vc.polyhedron = poly;

        // Extract the vertex indices from the computed polyhedron.
        vc.vertices_indices.clear();
        for (auto vit = vc.polyhedron.vertices_begin(); vit != vc.polyhedron.vertices_end(); ++vit)
        {
            Point p = vit->point();
            int vertex_index = find_vertex_index(voronoiDiagram, p);
            vc.vertices_indices.push_back(vertex_index);
        }

        // Extract the facets from the polyhedron. For each facet, record the list of vertex indices (and associated scalar values) that form its boundary.
        for (auto facet_it = vc.polyhedron.facets_begin(); facet_it != vc.polyhedron.facets_end(); ++facet_it)
        {
            VoronoiFacet vf;
            auto h = facet_it->facet_begin();
            do
            {
                Point p = h->vertex()->point();
                int vertex_index = find_vertex_index(voronoiDiagram, p);
                vf.vertices_indices.push_back(vertex_index);
                float value = voronoiDiagram.vertexValues[vertex_index];
                vf.vertex_values.push_back(value);
                ++h;
            } while (h != facet_it->facet_begin());
            int facet_index = voronoiDiagram.facets.size();
            voronoiDiagram.facets.push_back(vf);
            vc.facet_indices.push_back(facet_index);
        }

        // Add the newly constructed cell to the diagram and update the mapping.
        voronoiDiagram.cells.push_back(vc);
        delaunay_vertex->info().voronoiCellIndex = index;
        index++;
    }
}

//! @brief Constructs Voronoi edges from Delaunay facets.
void construct_voronoi_edges(
    VoronoiDiagram &voronoiDiagram,
    std::map<CGAL::Object, std::vector<Facet>, ObjectComparator> &voronoi_edge_to_delaunay_facet_map,
    Delaunay &dt)
{
    voronoiDiagram.edges.clear();
    for (Delaunay::Finite_facets_iterator fit = dt.finite_facets_begin();
         fit != dt.finite_facets_end(); ++fit)
    {
        Facet facet = *fit;
        Cell_handle c1 = facet.first;
        int i = facet.second;
        Cell_handle c2 = c1->neighbor(i);
        CGAL::Object vEdge = dt.dual(facet);
        /*         if (isDegenerate(vEdge))
                {
                    continue; // Skip edges where source == target, though rare with distinct indices
                } */
        int idx1 = -1;
        int idx2 = -1;
        if (!dt.is_infinite(c1))
            idx1 = c1->info().dualVoronoiVertexIndex;
        if (!dt.is_infinite(c2))
            idx2 = c2->info().dualVoronoiVertexIndex;

        // TODO: Duplicate check
        if (idx1 != -1 && idx2 != -1)
        {
            voronoiDiagram.edgeVertexIndices.push_back(std::make_pair(idx1, idx2));
        }
        else if (idx1 != -1)
        {
            voronoiDiagram.edgeVertexIndices.push_back(std::make_pair(idx1, -1));
        }
        else if (idx2 != -1)
        {
            voronoiDiagram.edgeVertexIndices.push_back(std::make_pair(idx2, -1));
        }
        else
        {
            std::cerr << "Error: finite facet with both adjacent cells infinite.\n";
        }

        voronoiDiagram.edges.push_back(vEdge);
        voronoi_edge_to_delaunay_facet_map[vEdge].push_back(facet);
    }
}

//! @brief Builds Voronoi cell edges for each edge in the diagram.
/*!
 * Creates VoronoiCellEdge entries for cells sharing each edge, collecting cell indices
 * from associated Delaunay facets.
 *
 * @param voronoiDiagram The Voronoi diagram to populate with cell edges.
 * @param voronoi_edge_to_delaunay_facet_map Map linking Voronoi edges to Delaunay facets.
 * @param dt The Delaunay triangulation.
 */
static void buildCellEdges(
    VoronoiDiagram &voronoiDiagram,
    const std::map<CGAL::Object, std::vector<Facet>, ObjectComparator> &voronoi_edge_to_delaunay_facet_map,
    Delaunay &dt)
{
    voronoiDiagram.cellEdges.clear();
    std::vector<std::unordered_set<int>> cellIndicesPerEdge(voronoiDiagram.edges.size());

    // Parallelize edge processing
    #pragma omp parallel for
    for (int edgeIdx = 0; edgeIdx < voronoiDiagram.edges.size(); ++edgeIdx)
    {
        const CGAL::Object &edgeObj = voronoiDiagram.edges[edgeIdx];
        auto it = voronoi_edge_to_delaunay_facet_map.find(edgeObj);
        if (it == voronoi_edge_to_delaunay_facet_map.end())
            continue;

        const std::vector<Facet> &sharedFacets = it->second;
        std::unordered_set<int> cellIndices;

        for (const Facet &f : sharedFacets)
        {
            Cell_handle c = f.first;
            if (dt.is_infinite(c))
                continue;

            for (int corner = 0; corner < 4; ++corner)
            {
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
static void linkCellEdges(
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
static void processEdgeMapping(
    VoronoiDiagram &voronoiDiagram,
    const CGAL::Object &edgeObj,
    int edgeIdx,
    CGAL::Epick::Iso_cuboid_3 &bbox)
{
    Segment3 seg;
    Ray3 ray;
    Line3 line;
    Point p1, p2;
    bool isSegment = false;

    if (CGAL::assign(seg, edgeObj))
    {
        p1 = seg.source();
        p2 = seg.target();
        isSegment = true;
    }
    else if (CGAL::assign(ray, edgeObj))
    {
        CGAL::Object clippedObj = CGAL::intersection(bbox, ray);
        Segment3 clippedSeg;
        if (CGAL::assign(clippedSeg, clippedObj))
        {
            p1 = clippedSeg.source();
            p2 = clippedSeg.target();
            isSegment = true;
        }
    }
    else if (CGAL::assign(line, edgeObj))
    {
        CGAL::Object clippedObj = CGAL::intersection(bbox, line);
        Segment3 clippedSeg;
        if (CGAL::assign(clippedSeg, clippedObj))
        {
            p1 = clippedSeg.source();
            p2 = clippedSeg.target();
            isSegment = true;
        }
    }

    if (isSegment)
    {
        int idx1 = find_vertex_index(voronoiDiagram, p1);
        int idx2 = find_vertex_index(voronoiDiagram, p2);
        if (idx1 != -1 && idx2 != -1)
        {
            int v1 = std::min(idx1, idx2);
            int v2 = std::max(idx1, idx2);
            voronoiDiagram.segmentVertexPairToEdgeIndex[{v1, v2}] = edgeIdx;
        }
    }
}

//! @brief Updates edge mappings for all Voronoi edges.
/*!
 * Processes all edges to update segmentVertexPairToEdgeIndex and cellEdgeLookup maps.
 *
 * @param voronoiDiagram The Voronoi diagram to update.
 * @param bbox The bounding box for intersection.
 */
static void updateEdgeMappings(
    VoronoiDiagram &voronoiDiagram,
    CGAL::Epick::Iso_cuboid_3 &bbox)
{
    for (int edgeIdx = 0; edgeIdx < (int)voronoiDiagram.edges.size(); ++edgeIdx)
    {
        const CGAL::Object &edgeObj = voronoiDiagram.edges[edgeIdx];
        processEdgeMapping(voronoiDiagram, edgeObj, edgeIdx, bbox);
    }

    voronoiDiagram.cellEdgeLookup.clear();
    for (int ceIdx = 0; ceIdx < (int)voronoiDiagram.cellEdges.size(); ++ceIdx)
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
 * @param voronoi_edge_to_delaunay_facet_map Map linking Voronoi edges to Delaunay facets.
 * @param bbox The bounding box used for clipping rays and lines.
 * @param dt The Delaunay triangulation.
 */
void construct_voronoi_cell_edges(
    VoronoiDiagram &voronoiDiagram,
    std::map<CGAL::Object, std::vector<Facet>, ObjectComparator> &voronoi_edge_to_delaunay_facet_map,
    CGAL::Epick::Iso_cuboid_3 &bbox,
    Delaunay &dt)
{
    voronoiDiagram.cellEdges.clear();

    std::clock_t start = std::clock();

    buildCellEdges(voronoiDiagram, voronoi_edge_to_delaunay_facet_map, dt);
    std::clock_t check1 = std::clock();
    double duration1 = static_cast<double>(check1 - start) / CLOCKS_PER_SEC;

    std::cout << "build cell edge Execution time: " << duration1 << " seconds" << std::endl;


    linkCellEdges(voronoiDiagram);
    std::clock_t check2 = std::clock();
    double duration2 = static_cast<double>(check2 - check1) / CLOCKS_PER_SEC;

    std::cout << "link cell edge Execution time: " << duration2 << " seconds" << std::endl;


    updateEdgeMappings(voronoiDiagram, bbox);

    std::clock_t check3 = std::clock();

    double duration3 = static_cast<double>(check3 - check2) / CLOCKS_PER_SEC;

    std::cout << "update edge mapping Execution time: " << duration3 << " seconds" << std::endl;
}

//! @brief Wrap up function of constructing voronoi diagram
void construct_voronoi_diagram(VoronoiDiagram &vd, VDC_PARAM &vdc_param, std::map<CGAL::Object, std::vector<Facet>, ObjectComparator> &voronoi_edge_to_delaunay_facet_map, ScalarGrid &grid, std::map<Point, float> &vertexValueMap, CGAL::Epick::Iso_cuboid_3 &bbox, Delaunay &dt)
{
    construct_voronoi_vertices(vd, dt);
    construct_voronoi_edges(vd, voronoi_edge_to_delaunay_facet_map, dt);
    vd.collapseSmallEdges(0.001, bbox);
    compute_voronoi_values(vd, grid, vertexValueMap);
    if (vdc_param.multi_isov)
    {
        if (vdc_param.convex_hull)
        {
            construct_voronoi_cells(vd, dt);
        }
        else
        {
            construct_voronoi_cells_non_convex_hull(vd, dt);
        }
        construct_voronoi_cell_edges(vd, voronoi_edge_to_delaunay_facet_map, bbox, dt);
    }
    vd.check();
}

// ！@brief Wrap up function for constructing iso surface
void construct_iso_surface(Delaunay &dt, VoronoiDiagram &vd, VDC_PARAM &vdc_param, IsoSurface &iso_surface, ScalarGrid &grid, Grid &data_grid, std::vector<Point> &activeCubeCenters, std::map<CGAL::Object, std::vector<Facet>, ObjectComparator> &voronoi_edge_to_delaunay_facet_map, std::map<Point, float> &vertexValueMap, CGAL::Epick::Iso_cuboid_3 &bbox, std::map<Point, int> &pointToIndexMap)
{
    if (vdc_param.multi_isov)
    {
        Compute_Isosurface_Vertices_Multi(vd, vdc_param.isovalue, iso_surface);
    }
    else
    {
        Compute_Isosurface_Vertices_Single(grid, vdc_param.isovalue, iso_surface, data_grid, activeCubeCenters);
    }

    if (vdc_param.multi_isov)
    {
        computeDualTrianglesMulti(vd, bbox, voronoi_edge_to_delaunay_facet_map, grid, vdc_param.isovalue, iso_surface);
    }
    else
    {
        computeDualTriangles(iso_surface, vd.edges, vertexValueMap, bbox, voronoi_edge_to_delaunay_facet_map, dt, grid, vdc_param.isovalue, pointToIndexMap);
    }
}
//! @brief Handles output mesh generation.
int handle_output_mesh(bool &retFlag, VoronoiDiagram &vd, VDC_PARAM &vdc_param, IsoSurface &iso_surface, std::map<Point, int> &pointToIndexMap)
{
    retFlag = true;

    std::cout << "Result file at: " << vdc_param.output_filename << std::endl;
    // Use locations of isosurface vertices as vertices of Delaunay triangles and write the output mesh
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
    else
    {
        if (vdc_param.output_format == "off")
        {
            writeOFFSingle(vdc_param.output_filename, iso_surface.isosurfaceVertices, iso_surface.isosurfaceTrianglesSingle, pointToIndexMap);
        }
        else if (vdc_param.output_format == "ply")
        {
            writePLYSingle(vdc_param.output_filename, iso_surface.isosurfaceVertices, iso_surface.isosurfaceTrianglesSingle, pointToIndexMap);
        }
        else
        {
            std::cerr << "Unsupported output format: " << vdc_param.output_format << std::endl;
            return EXIT_FAILURE;
        }
    }
    retFlag = false;
    return {};
}
