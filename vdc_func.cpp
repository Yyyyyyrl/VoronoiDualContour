//! @file vdc_func.cpp
//! @brief Implementation of functions for Voronoi Diagram and Isosurface computation.

#include "vdc_func.h"

//! @brief Helper function - returns the index of the vertex matching p, or -1 if not found.


int find_vertex_index(const VoronoiDiagram &vd, const Point &p)
{
    return vd.find_vertex(p);
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
static void processSegmentEdge(
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
 * @param bbox The bounding box for intersection.
 * @param grid The scalar grid for interpolation.
 * @param isovalue The isovalue for bipolarity check.
 * @param dt The Delaunay triangulation.
 * @param dualTriangles Vector to store generated triangles.
 */
static void processRayEdge(
    const Ray3 &ray,
    VoronoiEdge &edge,
    VoronoiDiagram &vd,
    CGAL::Epick::Iso_cuboid_3 &bbox,
    ScalarGrid &grid,
    float isovalue,
    Delaunay &dt,
    std::vector<DelaunayTriangle> &dualTriangles)
{
    CGAL::Object intersectObj = CGAL::intersection(bbox, ray);
    Segment3 iseg;
    if (CGAL::assign(iseg, intersectObj))
    {
        Point v1 = iseg.source();
        Point v2 = iseg.target();
        int idx_v1 = find_vertex_index(vd, v1);
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
 * @param dt The Delaunay triangulation.
 * @param dualTriangles Vector to store generated triangles.
 */
static void processLineEdge(
    const Line3 &line,
    VoronoiEdge &edge,
    ScalarGrid &grid,
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
 * @param bbox Bounding box of the computational domain.
 * @param dt Delaunay triangulation structure.
 * @param grid Scalar grid containing scalar values.
 * @param isovalue The isovalue used for computing.
 */
void computeDualTriangles(
    IsoSurface &iso_surface,
    VoronoiDiagram &vd,
    CGAL::Epick::Iso_cuboid_3 &bbox,
    Delaunay &dt,
    ScalarGrid &grid,
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
            processSegmentEdge( edge, vd, isovalue, dt, dualTriangles);
        }
        else if (edge.type == 1)
        {
            CGAL::assign(ray, edge.edgeObject);
            processRayEdge(ray, edge, vd, bbox, grid, isovalue, dt, dualTriangles);
        }
        else if (edge.type == 2)
        {
            processLineEdge(line, edge, grid, isovalue, bbox, dt, dualTriangles);
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
 * @param iso_surface The isosurface to store triangles.
 */
static void processSegmentEdgeMulti(
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
                bool isValid = selectIsovertices(voronoiDiagram, facet, globalEdgeIndex, idx1, idx2, idx3, cellIndex1, cellIndex2, cellIndex3);
                int iOrient = get_orientation(facet.second, v1, v2, val1, val2);
                generateTriangleMulti(iso_surface, idx1, idx2, idx3, iOrient, isValid);
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
static void processRayEdgeMulti(
    const Ray3 &ray,
    std::vector<Facet> dualDelaunayFacets,
    VoronoiDiagram &voronoiDiagram,
    ScalarGrid &grid,
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
        int idx_v1 = find_vertex_index(voronoiDiagram, v1);
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
                    generateTriangleMulti(iso_surface, idx1, idx2, idx3, iOrient, isValid);
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
static void processLineEdgeMulti(
    const Line3 &line,
    std::vector<Facet> dualDelaunayFacets,
    VoronoiDiagram &voronoiDiagram,
    ScalarGrid &grid,
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
                    generateTriangleMulti(iso_surface, idx1, idx2, idx3, iOrient, isValid);
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
void computeDualTrianglesMulti(
    VoronoiDiagram &voronoiDiagram,
    CGAL::Epick::Iso_cuboid_3 &bbox,
    ScalarGrid &grid,
    float isovalue,
    IsoSurface &iso_surface)
{
    for (const auto &edge : voronoiDiagram.edges)
    {
        Segment3 seg;
        Ray3 ray;
        Line3 line;
        std::vector<Facet> dualDelaunayFacets = edge.delaunayFacets;

        if (edge.type == 0)
        {
            processSegmentEdgeMulti(edge, voronoiDiagram, isovalue, iso_surface);
        }
        else if (edge.type == 1)
        {
            CGAL::assign(ray, edge.edgeObject);
            processRayEdgeMulti(ray, dualDelaunayFacets, voronoiDiagram, grid, isovalue, bbox, iso_surface);
        }
        else if (edge.type == 2)
        {
            CGAL::assign(line, edge.edgeObject);
            processLineEdgeMulti(line, dualDelaunayFacets, voronoiDiagram, grid, isovalue, bbox, iso_surface);
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
static void collectDelaunayPoints(Grid &grid,
                                  const std::vector<std::vector<GRID_FACETS>> &grid_facets,
                                  const std::vector<Point> &activeCubeCenters,
                                  VDC_PARAM &vdc_param,
                                  std::vector<Point> &delaunay_points,
                                  std::vector<int> &dummy_point_indices)
{
    // Start with active cube centers
    delaunay_points = activeCubeCenters;
    dummy_point_indices.clear();

    if (vdc_param.multi_isov) {
        // For each facet, generate dummy points and record their indices
        for (size_t d = 0; d < grid_facets.size(); ++d) {
            for (const auto &f : grid_facets[d]) {
                auto pointsf = add_dummy_from_facet(f, grid);
                for (const auto &p : pointsf) {
                    delaunay_points.push_back(p);
                    dummy_point_indices.push_back(static_cast<int>(delaunay_points.size()) - 1);
                }
            }
        }
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
static Vertex_handle insertPointIntoTriangulation(Delaunay &dt,
                                                  const Point &p,
                                                  int index,
                                                  bool is_dummy)
{
    // Insert point and retrieve handle
    Vertex_handle vh = dt.insert(p);
    // Immediately assign index and dummy status
    vh->info().index    = index;
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
void construct_delaunay_triangulation(Delaunay &dt,
                                      Grid &grid,
                                      const std::vector<std::vector<GRID_FACETS>> &grid_facets,
                                      VDC_PARAM &vdc_param,
                                      std::vector<Point> &activeCubeCenters)
{
    // Build point list and dummy indices
    std::vector<Point> delaunay_points;
    std::vector<int>   dummy_point_indices;
    collectDelaunayPoints(grid, grid_facets, activeCubeCenters,
                           vdc_param, delaunay_points, dummy_point_indices);


    // Clear existing triangulation
    dt.clear();

    // Insert each point using the helper
    for (int i = 0; i < static_cast<int>(delaunay_points.size()); ++i) {
        bool is_dummy = (std::find(dummy_point_indices.begin(), dummy_point_indices.end(), i) != dummy_point_indices.end());
        insertPointIntoTriangulation(dt,delaunay_points[i],i,is_dummy);
    }

}


//! @brief Constructs Voronoi vertices for the given voronoi Diagram instance.
void construct_voronoi_vertices(VoronoiDiagram &voronoiDiagram, Delaunay &dt) {
    voronoiDiagram.vertices.clear();
    const double EPSILON = 1e-6;
    const double SCALE_FACTOR = 1e6;

    for (Delaunay::Finite_cells_iterator cit = dt.finite_cells_begin(); cit != dt.finite_cells_end(); ++cit) {
        Point P = dt.dual(cit);
        int ix = static_cast<int>(std::round(P.x() * SCALE_FACTOR));
        int iy = static_cast<int>(std::round(P.y() * SCALE_FACTOR));
        int iz = static_cast<int>(std::round(P.z() * SCALE_FACTOR));
        std::tuple<int, int, int> key(ix, iy, iz);

        auto it = voronoiDiagram.vertexMap.find(key);
        int vertex_index = -1;
        if (it != voronoiDiagram.vertexMap.end()) {
            for (int idx : it->second) {
                if (CGAL::squared_distance(P, voronoiDiagram.vertices[idx].coord) < EPSILON * EPSILON) {
                    vertex_index = idx;
                    break;
                }
            }
        }
        if (vertex_index == -1) {
            vertex_index = voronoiDiagram.vertices.size();
            VoronoiVertex vVertex(P);
            vVertex.index = vertex_index;
            voronoiDiagram.vertices.push_back(vVertex);
            voronoiDiagram.vertexMap[key].push_back(vertex_index);
        }
        cit->info().dualVoronoiVertexIndex = vertex_index;
    }
}

//! @brief Computes Voronoi Vertex values using scalar grid interpolation
void compute_voronoi_values(VoronoiDiagram &voronoiDiagram, ScalarGrid &grid)
{
    for (size_t i = 0; i < voronoiDiagram.vertices.size(); ++i)
    {
        Point vertex = voronoiDiagram.vertices[i].coord;
        voronoiDiagram.vertices[i].value = trilinear_interpolate(vertex, grid);
        float value = trilinear_interpolate(vertex, grid);
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
            vertex_points.push_back(voronoiDiagram.vertices[idx].coord);
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
        const Point &pt = vertices[idx].coord;
        sumX += pt.x();
        sumY += pt.y();
        sumZ += pt.z();
    }
    double count = static_cast<double>(indices.size());
    Point center(sumX / count, sumY / count, sumZ / count);

    // Sort indices by their angle (using atan2) in the (vRef, vRef2) coordinate system.
    std::sort(indices.begin(), indices.end(), [&](int a, int b)
              {
        const Point &pa = vertices[a].coord;
        const Point &pb = vertices[b].coord;

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

    std::set<int> uniqueVertexIndices;
    for (Cell_handle c : incidentCells)
    {
        Point vor_vertex = dt.dual(c);
        int vertex_index = voronoiDiagram.find_vertex(vor_vertex);
        if (vertex_index >= 0 && vertex_index < voronoiDiagram.vertices.size()) {
            uniqueVertexIndices.insert(vertex_index);
        } else {
            std::cerr << "[ERROR] Vertex not found for point " << vor_vertex << "\n";
        }

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
            Point vor_vertex = dt.dual(cc);
            int newIdx = voronoiDiagram.find_vertex(vor_vertex);
            if (newIdx >= 0 && newIdx < voronoiDiagram.vertices.size()) {
                facetVertexIndices.push_back(newIdx);
            } else {
                std::cerr << "[ERROR] Vertex not found for point " << vor_vertex << "\n";
            }
            cell_count++;
        }
        ++cc;
    } while (cc != start);

    std::set<int> unique_vertices(facetVertexIndices.begin(), facetVertexIndices.end());

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
        Point p0 = voronoiDiagram.vertices[uniqueFacetVertices[0]].coord;
        Point p1 = voronoiDiagram.vertices[uniqueFacetVertices[1]].coord;
        Point p2 = voronoiDiagram.vertices[uniqueFacetVertices[2]].coord;
        Vector3 normal = CGAL::cross_product(p1 - p0, p2 - p0);
        Point centroid(0, 0, 0);
        for (int idx : uniqueFacetVertices)
        {
            centroid = centroid + (voronoiDiagram.vertices[idx].coord - CGAL::ORIGIN);
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

    for (const Edge &ed : incidentEdges)
    {
        VoronoiFacet facet = buildFacetFromEdge(dt, ed, delaunay_vertex, voronoiDiagram, vc.facet_indices, edge_to_facets);
        if (!facet.vertices_indices.empty())
        {
            int facetIndex = voronoiDiagram.facets.size() - 1; // Assuming facet was just added
            if (std::find(vc.facet_indices.begin(), vc.facet_indices.end(), facetIndex) == vc.facet_indices.end())
            {
                vc.facet_indices.push_back(facetIndex);
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
            continue; // Skip dummy vertices

        VoronoiCell vc = createVoronoiCell(v, cellIndex);
        collectCellVertices(dt, v, voronoiDiagram, vc.vertices_indices);
        processIncidentEdges(dt, v, voronoiDiagram, vc, edge_to_facets);

        // Validate the number of facets
        if (vc.facet_indices.size() < 4)
        {
            std::cerr << "[WARNING] Cell " << cellIndex << " has only " << vc.facet_indices.size() << " facets, skipping\n";
        }
        else
        {
            voronoiDiagram.cells.push_back(vc);
            v->info().voronoiCellIndex = cellIndex;
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
        }
    }

}

// Helper function to check if two directions are approximately equal
bool directionsEqual(const Vector3& d1, const Vector3& d2, double epsilon) {
    Vector3 n1 = d1 / std::sqrt(d1.squared_length()); // Normalize d1
    Vector3 n2 = d2 / std::sqrt(d2.squared_length()); // Normalize d2
    return (n1 - n2).squared_length() < epsilon * epsilon; // Compare squared distance of normalized vectors
}

//! @brief Constructs Voronoi edges from Delaunay facets.
void construct_voronoi_edges(VoronoiDiagram &voronoiDiagram, Delaunay &dt) {
    voronoiDiagram.edges.clear();
    std::map<std::pair<int, int>, int> segmentMap; // Maps sorted vertex pairs to edge indices
    std::map<int, std::vector<std::pair<Vector3, int>>> rayMap; // Maps vertex index to (direction, edgeIdx) pairs
    const double EPSILON = 1e-6;

    for (Delaunay::Finite_facets_iterator fit = dt.finite_facets_begin(); fit != dt.finite_facets_end(); ++fit) {
        Facet facet = *fit;
        CGAL::Object edgeobj = dt.dual(facet);
        Segment3 seg;
        Ray3 ray;

        if (CGAL::assign(seg, edgeobj)) {
            Cell_handle c1 = facet.first;
            Cell_handle c2 = c1->neighbor(facet.second);
            int idx1 = dt.is_infinite(c1) ? -1 : c1->info().dualVoronoiVertexIndex;
            int idx2 = dt.is_infinite(c2) ? -1 : c2->info().dualVoronoiVertexIndex;

            if (idx1 != -1 && idx2 != -1 && idx1 != idx2) { // Finite segment
                int v1 = std::min(idx1, idx2);
                int v2 = std::max(idx1, idx2);
                auto it = segmentMap.find({v1, v2});
                if (it != segmentMap.end()) {
                    voronoiDiagram.edges[it->second].delaunayFacets.push_back(facet);
                } else {
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
        } else if (CGAL::assign(ray, edgeobj)) {
            Cell_handle c1 = facet.first;
            Cell_handle c2 = c1->neighbor(facet.second);
            int vertex1 = -1;
            if (!dt.is_infinite(c1)) {
                vertex1 = c1->info().dualVoronoiVertexIndex;
            } else if (!dt.is_infinite(c2)) {
                vertex1 = c2->info().dualVoronoiVertexIndex;
            }
            if (vertex1 != -1) {
                Vector3 dir = ray.direction().vector();
                bool found = false;
                auto it = rayMap.find(vertex1);
                if (it != rayMap.end()) {
                    for (const auto& pair : it->second) {
                        if (directionsEqual(pair.first, dir, EPSILON)) {
                            voronoiDiagram.edges[pair.second].delaunayFacets.push_back(facet);
                            found = true;
                            break;
                        }
                    }
                }
                if (!found) {
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
        // Lines are not expected for finite facets, so we skip them
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
static void buildCellEdges(
    VoronoiDiagram &voronoiDiagram,
    Delaunay &dt)
{
    voronoiDiagram.cellEdges.clear();
    std::vector<std::unordered_set<int>> cellIndicesPerEdge(voronoiDiagram.edges.size());

    // Parallelize edge processing
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
        std::cout << "[DEBUG] Processing edge " << edgeIdx << " with vertices: (" << idx1 << ", " << idx2 << ") with edge index: " << edgeIdx << std::endl;
        if (idx1 == -1) {
            std::cout << "[WARNING] Failed to find vertex index for point: " << p1 << std::endl;
        } else if (idx2 == -1 ){
            std::cout << "[WARNING] Failed to find vertex index for point: " << p2 << std::endl;
        } else
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
        const CGAL::Object &edgeObj = voronoiDiagram.edges[edgeIdx].edgeObject;
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

    buildCellEdges(voronoiDiagram, dt);
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
void construct_voronoi_diagram(VoronoiDiagram &vd, VDC_PARAM &vdc_param, ScalarGrid &grid, CGAL::Epick::Iso_cuboid_3 &bbox, Delaunay &dt)
{
    construct_voronoi_vertices(vd, dt);
    construct_voronoi_edges(vd, dt);
    VoronoiDiagram vd2 = collapseSmallEdges(vd, 0.001, bbox);
    compute_voronoi_values(vd2, grid);
    if (vdc_param.multi_isov)
    {
        if (vdc_param.convex_hull)
        {
            construct_voronoi_cells_as_convex_hull(vd2, dt);
        }
        else
        {
            construct_voronoi_cells_from_delaunay_triangulation(vd2, dt);
        }
        construct_voronoi_cell_edges(vd2, bbox, dt);
    }
    vd2.check();
    vd = std::move(vd2);

    if (debug) {
        Point hole_center(52, 34, 35);
        double threshold = 2.0;
        std::cout << "===========\n[DEBUG] Voronoi vertices near hole:\n";
        for (size_t i = 0; i < vd.vertices.size(); ++i) {
            const Point& p = vd.vertices[i].coord;
            if (CGAL::squared_distance(p, hole_center) < threshold * threshold) {
                std::cout << "[DEBUG] Voronoi vertex " << i << ": " << p << "\n";
            }
        }
        std::cout << "===========\n";
    }
}

// ！@brief Wrap up function for constructing iso surface
void construct_iso_surface(Delaunay &dt, VoronoiDiagram &vd, VDC_PARAM &vdc_param, IsoSurface &iso_surface, ScalarGrid &grid, Grid &data_grid, std::vector<Point> &activeCubeCenters, CGAL::Epick::Iso_cuboid_3 &bbox)
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
        computeDualTrianglesMulti(vd, bbox, grid, vdc_param.isovalue, iso_surface);
    }
    else
    {
        computeDualTriangles(iso_surface, vd, bbox, dt, grid, vdc_param.isovalue);
    }

    if (debug) {
        std::set<int> problem_vertices = {1730, 1554, 1731, 1553};
        std::cout << "[DEBUG] Cell assignments for isovertices:\n";
        for (size_t cell_idx = 0; cell_idx < vd.cells.size(); ++cell_idx) {
            const VoronoiCell& vc = vd.cells[cell_idx];
            int start = vc.isoVertexStartIndex;
            int end = start + vc.numIsoVertices;
            for (int idx = start; idx < end; ++idx) {
                if (problem_vertices.count(idx) > 0) {
                    int cycle_idx = idx - start;
                    std::cout << "[DEBUG] Problematic Isovertex " << idx << " in Cell " << cell_idx 
                            << ", Cycle " << cycle_idx << "\n";
                    if (cycle_idx < vc.cycles.size()) {
                        const Cycle& cycle = vc.cycles[cycle_idx];
                        std::cout << "[DEBUG] Cycle midpoints: ";
                        for (int mid_idx : cycle.midpoint_indices) {
                            std::cout << mid_idx << " ";
                        }
                        std::cout << "\n";
                    }
                }
            }
        }
    }

    if (debug) {
    Point hole_center(52, 34, 35);
    double threshold = 2.0;
    std::cout << "[DEBUG] Isovertices near hole region:\n";
    for (size_t i = 0; i < iso_surface.isosurfaceVertices.size(); ++i) {
        const Point& p = iso_surface.isosurfaceVertices[i];
        if (CGAL::squared_distance(p, hole_center) < threshold * threshold) {
            std::cout << "[DEBUG] Isovertex " << i << ": " << p << "\n";
        }
    }

}
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
