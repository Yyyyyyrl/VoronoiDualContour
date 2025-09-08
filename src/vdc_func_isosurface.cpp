#include "vdc_func.h"

// ===== Debug instrumentation for isosurface pipeline =====
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <unordered_set>

static int ISO_DBG_FOCUS_CELL  = -1;
static int ISO_DBG_FOCUS_GFACET= -1;
static int ISO_DBG_FOCUS_EDGE  = -1;
static bool ISO_DBG_ONLY_ERRORS= false;
static bool ISO_DBG_ENABLED    = true;  // set false to silence all

static int iso_getenv_int(const char* k, int defv){ if(const char* v=std::getenv(k)){ try{ return std::stoi(v);}catch(...){}} return defv; }
static bool iso_getenv_bool(const char* k, bool defv){ if(const char* v=std::getenv(k)){ std::string s(v); for(char &c: s) c=std::tolower(c); if(s=="1"||s=="true"||s=="yes"||s=="on") return true; if(s=="0"||s=="false"||s=="no"||s=="off") return false;} return defv; }

static void ISO_DBG_LOAD_ENV(){
    ISO_DBG_FOCUS_CELL   = iso_getenv_int("ISO_DEBUG_CELL", -1);
    ISO_DBG_FOCUS_GFACET = iso_getenv_int("ISO_DEBUG_GLOBAL_FACET", -1);
    ISO_DBG_FOCUS_EDGE   = iso_getenv_int("ISO_DEBUG_EDGE", -1);
    ISO_DBG_ONLY_ERRORS  = iso_getenv_bool("ISO_DEBUG_ONLY_ERRORS", false);
    ISO_DBG_ENABLED      = !iso_getenv_bool("ISO_DEBUG_OFF", false);
}

static inline bool iso_dbg_cell_ok(int c){ return ISO_DBG_FOCUS_CELL < 0 || ISO_DBG_FOCUS_CELL == c; }
static inline bool iso_dbg_gfacet_ok(int g){ return ISO_DBG_FOCUS_GFACET < 0 || ISO_DBG_FOCUS_GFACET == g; }
static inline bool iso_dbg_edge_ok(int e){ return ISO_DBG_FOCUS_EDGE < 0 || ISO_DBG_FOCUS_EDGE == e; }

struct IsoStats {
    size_t cells_seen=0, cells_with_midpts=0, cells_with_zero_connect=0, cycles_total=0;
    size_t edges_seg=0, edges_ray=0, edges_line=0;
    size_t seg_bip=0, seg_skip=0, ray_bip=0, ray_skip=0, line_bip=0, line_skip=0;
    size_t tri_ok=0, tri_bad=0, tri_dupverts=0, sel_fail=0, sel_fallback=0;
    void dump_summary() const {
        if(!ISO_DBG_ENABLED) return;
        std::cerr << "[ISO] ===== Isosurface Summary =====\n"
                  << "[ISO] Cells: seen=" << cells_seen
                  << " with_midpoints=" << cells_with_midpts
                  << " zero_connections=" << cells_with_zero_connect
                  << " cycles=" << cycles_total << "\n"
                  << "[ISO] Edges: seg=" << edges_seg << " ray=" << edges_ray << " line=" << edges_line << "\n"
                  << "[ISO] Bipolar pass: seg=" << seg_bip << "/" << (seg_bip+seg_skip)
                  << " ray=" << ray_bip << "/" << (ray_bip+ray_skip)
                  << " line=" << line_bip << "/" << (line_bip+line_skip) << "\n"
                  << "[ISO] Triangles: ok=" << tri_ok
                  << " bad=" << tri_bad
                  << " dupVerts=" << tri_dupverts
                  << " isoSelectFail=" << sel_fail
                  << " isoSelectFallback=" << sel_fallback << "\n"
                  << "[ISO] Filters: CELL=" << ISO_DBG_FOCUS_CELL
                  << " GFACET=" << ISO_DBG_FOCUS_GFACET
                  << " EDGE=" << ISO_DBG_FOCUS_EDGE
                  << " ONLY_ERRORS=" << (ISO_DBG_ONLY_ERRORS? "1":"0") << "\n"
                  << "[ISO] =================================\n";
    }
};
static IsoStats ISO_STATS;
// ===== End debug instrumentation header =====

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

// Pick an isovertex index for a given (cell, globalEdge).
// Strategy: try exact (cell,edge)→cellEdge→cycle mapping; if none,
// fall back to the first iso-vertex in that cell; else fail (-1).
static inline int select_isovertex_from_cell_edge(
    const VoronoiDiagram &vd,
    int cellIndex, int globalEdgeIndex)
{
    // Guard: valid cell index
    if (cellIndex < 0 || cellIndex >= static_cast<int>(vd.cells.size())) {
        ISO_STATS.sel_fail++;
        return -1;
    }

    // 1) Exact (cell,edge) lookup → walk ring to a cellEdge carrying cycles
    auto it = vd.cellEdgeLookup.find(std::make_pair(cellIndex, globalEdgeIndex));
    if (it != vd.cellEdgeLookup.end()) {
        int ceIdx = it->second;
        const int start = ceIdx;

        // Walk nextCellEdge ring until a non-empty cycleIndices is found
        while (ceIdx >= 0 &&
               ceIdx < static_cast<int>(vd.cellEdges.size()) &&
               vd.cellEdges[ceIdx].cycleIndices.empty())
        {
            const int nxt = vd.cellEdges[ceIdx].nextCellEdge;
            if (nxt < 0 || nxt == start) break;
            ceIdx = nxt;
        }

        if (ceIdx >= 0 &&
            ceIdx < static_cast<int>(vd.cellEdges.size()) &&
            !vd.cellEdges[ceIdx].cycleIndices.empty())
        {
            const int cycLocal = vd.cellEdges[ceIdx].cycleIndices[0];
            const VoronoiCell &vc = vd.cells[cellIndex];

            if (cycLocal >= 0 && cycLocal < vc.numIsoVertices) {
                if (ISO_DBG_ENABLED && iso_dbg_cell_ok(cellIndex) && iso_dbg_edge_ok(globalEdgeIndex)) {
                    std::cerr << "[ISO] pick cell " << cellIndex
                              << " edge " << globalEdgeIndex
                              << " via cellEdge#" << ceIdx
                              << " cycleLocal=" << cycLocal << "\n";
                }
                return vc.isoVertexStartIndex + cycLocal;
            }
        }
    }

    // 2) Fallback: use the first iso-vertex in this cell if any exist
    const VoronoiCell &vc = vd.cells[cellIndex];
    if (vc.numIsoVertices > 0) {
        if (ISO_DBG_ENABLED && iso_dbg_cell_ok(cellIndex) && iso_dbg_edge_ok(globalEdgeIndex)) {
            std::cerr << "[ISO] pick cell " << cellIndex
                      << " edge " << globalEdgeIndex
                      << " FALLBACK first isoVert\n";
        }
        ISO_STATS.sel_fallback++;
        return vc.isoVertexStartIndex; // first cycle’s isovertex in this cell
    }

    // 3) No iso-vertices in this cell at all
    if (ISO_DBG_ENABLED && iso_dbg_cell_ok(cellIndex) && iso_dbg_edge_ok(globalEdgeIndex)) {
        std::cerr << "[ISO] pick cell " << cellIndex
                  << " edge " << globalEdgeIndex
                  << " FAIL (no iso-verts)\n";
    }
    ISO_STATS.sel_fail++;
    return -1;
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
// Emit a triangle (by isovertex indices) if valid; count & log otherwise.
// iOrient has the same meaning as in single-isov mode: >=0 keeps (1,2,3), <0 swaps 2<->3.
static void generate_triangle_multi(
    IsoSurface &iso_surface,
    int idx1, int idx2, int idx3,
    int iOrient,
    bool isValid)
{
    if (!isValid) {
        ISO_STATS.tri_bad++;
        if (ISO_DBG_ENABLED && !ISO_DBG_ONLY_ERRORS) {
            std::cerr << "[ISO] TRI SKIP invalid indices (" << idx1 << "," << idx2 << "," << idx3 << ")\n";
        }
        return;
    }
    if (idx1 == idx2 || idx2 == idx3 || idx1 == idx3) {
        ISO_STATS.tri_dupverts++;
        if (ISO_DBG_ENABLED && !ISO_DBG_ONLY_ERRORS) {
            std::cerr << "[ISO] TRI SKIP duplicate indices (" << idx1 << "," << idx2 << "," << idx3 << ")\n";
        }
        return;
    }

    ISO_STATS.tri_ok++;
    if (ISO_DBG_ENABLED && !ISO_DBG_ONLY_ERRORS) {
        std::cerr << "[ISO] TRI EMIT (" << idx1 << "," << idx2 << "," << idx3 << ")\n";
    }

    if (iOrient >= 0) {
        iso_surface.isosurfaceTrianglesMulti.emplace_back(idx1, idx2, idx3);
    } else {
        iso_surface.isosurfaceTrianglesMulti.emplace_back(idx1, idx3, idx2);
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
static inline bool select_isovertices(
    const VoronoiDiagram &voronoiDiagram,
    const Facet &facet,
    int globalEdgeIndex,
    int &idx1, int &idx2, int &idx3,
    int &cellIndex1, int &cellIndex2, int &cellIndex3)
{
    idx1 = idx2 = idx3 = -1;
    cellIndex1 = cellIndex2 = cellIndex3 = -1;

    const int iFacet = facet.second;
    Cell_handle c = facet.first;

    const int d1 = (iFacet + 1) % 4;
    const int d2 = (iFacet + 2) % 4;
    const int d3 = (iFacet + 3) % 4;

    Vertex_handle v1 = c->vertex(d1);
    Vertex_handle v2 = c->vertex(d2);
    Vertex_handle v3 = c->vertex(d3);

    // Skip facets involving dummy vertices
    if (v1->info().is_dummy || v2->info().is_dummy || v3->info().is_dummy) {
        return false;
    }

    cellIndex1 = v1->info().voronoiCellIndex;
    cellIndex2 = v2->info().voronoiCellIndex;
    cellIndex3 = v3->info().voronoiCellIndex;

    idx1 = select_isovertex_from_cell_edge(voronoiDiagram, cellIndex1, globalEdgeIndex);
    idx2 = select_isovertex_from_cell_edge(voronoiDiagram, cellIndex2, globalEdgeIndex);
    idx3 = select_isovertex_from_cell_edge(voronoiDiagram, cellIndex3, globalEdgeIndex);

    if (ISO_DBG_ENABLED && iso_dbg_edge_ok(globalEdgeIndex)) {
        std::cerr << "[ISO] facet pick edge=" << globalEdgeIndex
                  << " cells=(" << cellIndex1 << "," << cellIndex2 << "," << cellIndex3 << ")"
                  << " iso=(" << idx1 << "," << idx2 << "," << idx3 << ")";
        if (idx1 < 0 || idx2 < 0 || idx3 < 0) std::cerr << " [MISS]";
        else if (idx1 == idx2 || idx2 == idx3 || idx1 == idx3) std::cerr << " [DUP]";
        else std::cerr << " [OK]";
        std::cerr << "\n";
    }

    return (idx1 >= 0 && idx2 >= 0 && idx3 >= 0 &&
            idx1 != idx2 && idx2 != idx3 && idx1 != idx3);
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
    ISO_STATS.edges_seg++;

    int idx_v1 = edge.vertex1;
    int idx_v2 = edge.vertex2;
    Point v1 = voronoiDiagram.vertices[idx_v1].coord;
    Point v2 = voronoiDiagram.vertices[idx_v2].coord;
    float val1 = voronoiDiagram.vertices[idx_v1].value;
    float val2 = voronoiDiagram.vertices[idx_v2].value;

    if (!is_bipolar(val1, val2, isovalue))
    {
        ISO_STATS.seg_skip++;
        if(ISO_DBG_ENABLED && !ISO_DBG_ONLY_ERRORS){ std::cerr << "[ISO] SEG non-bipolar v=("<<val1<<","<<val2<<") iso="<<isovalue<<"\n"; }
        return;
    }
    ISO_STATS.seg_bip++;
    {
        if (idx_v1 > idx_v2)
            std::swap(idx_v1, idx_v2);
        auto itEdge = voronoiDiagram.segmentVertexPairToEdgeIndex.find(std::make_pair(idx_v1, idx_v2));
        if (itEdge == voronoiDiagram.segmentVertexPairToEdgeIndex.end())
            return;

        int globalEdgeIndex = itEdge->second;
        if(ISO_DBG_ENABLED && iso_dbg_edge_ok(globalEdgeIndex)) std::cerr << "[ISO] SEG bipolar edge → globalEdge="<<globalEdgeIndex<<" dualFacets="<<edge.delaunayFacets.size()<<"\n";

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
 * @param globalEdgeIndex The global edge index.
 * @param source_pt The index of the source point.
 * @param ray The ray edge to process.
 * @param edge The CGAL object representing the edge.
 * @param voronoiDiagram The Voronoi diagram containing edge and cell data.
 * @param grid The scalar grid for interpolation.
 * @param isovalue The isovalue for bipolarity check.
 * @param bbox The bounding box for intersection.
 * @param iso_surface The isosurface to store triangles.
 */
static void process_ray_edge_multi(
    int globalEdgeIndex,
    int source_pt,
    const Ray3 &ray,
    std::vector<Facet> dualDelaunayFacets,
    VoronoiDiagram &voronoiDiagram,
    UnifiedGrid &grid,
    float isovalue,
    CGAL::Epick::Iso_cuboid_3 &bbox,
    IsoSurface &iso_surface)
{
    ISO_STATS.edges_ray++;
    CGAL::Object intersectObj = CGAL::intersection(bbox, ray);
    Segment3 iseg;
    if (CGAL::assign(iseg, intersectObj))
    {
        Point v1 = ray.source();
        Point v2 = iseg.target();
        float val1 = voronoiDiagram.vertices[source_pt].value;
        float val2 = trilinear_interpolate(adjust_outside_bound_points(v2, grid, v1, v2), grid);

        if (!is_bipolar(val1, val2, isovalue)){
            ISO_STATS.ray_skip++;
            if(ISO_DBG_ENABLED && iso_dbg_edge_ok(globalEdgeIndex) && !ISO_DBG_ONLY_ERRORS){ std::cerr << "[ISO] RAY non-bipolar edge="<<globalEdgeIndex<<" v=("<<val1<<","<<val2<<")\n"; }
            return; }
        ISO_STATS.ray_bip++;
        if(ISO_DBG_ENABLED && iso_dbg_edge_ok(globalEdgeIndex)) std::cerr << "[ISO] RAY bipolar edge="<<globalEdgeIndex<<" dualFacets="<<dualDelaunayFacets.size()<<"\n";

        for (const auto &facet : dualDelaunayFacets)
        {
            int idx1, idx2, idx3, cellIndex1, cellIndex2, cellIndex3;
            // Use the same facet/edge-aware vertex selection used for finite segments.
            bool ok = select_isovertices(voronoiDiagram, facet, globalEdgeIndex,
                                         idx1, idx2, idx3,
                                         cellIndex1, cellIndex2, cellIndex3);
            int iOrient = get_orientation(facet.second, v1, v2, val1, val2);
            bool isValid = ok && (idx1 != idx2 && idx2 != idx3 && idx1 != idx3);
            generate_triangle_multi(iso_surface, idx1, idx2, idx3, iOrient, isValid);
        }
    }
}

//! @brief Processes a line edge for multi-isovertex triangle computation.
/*!
 * Intersects the line with the bounding box, checks bipolarity, and generates
 * triangles using the first isovertex from each cell.
 *
 * @param globalEdgeIndex The global index of the edge.
 * @param line The line edge to process.
 * @param edge The CGAL object representing the edge.
 * @param voronoiDiagram The Voronoi diagram containing edge and cell data.
 * @param grid The scalar grid for interpolation.
 * @param isovalue The isovalue for bipolarity check.
 * @param bbox The bounding box for intersection.
 * @param iso_surface The isosurface to store triangles.
 */
static void process_line_edge_multi(
    int globalEdgeIndex,
    const Line3 &line,
    std::vector<Facet> dualDelaunayFacets,
    VoronoiDiagram &voronoiDiagram,
    UnifiedGrid &grid,
    float isovalue,
    CGAL::Epick::Iso_cuboid_3 &bbox,
    IsoSurface &iso_surface)
{
    ISO_STATS.edges_line++;
    CGAL::Object intersectObj = CGAL::intersection(bbox, line);
    Segment3 iseg;
    if (CGAL::assign(iseg, intersectObj))
    {
        Point v1 = iseg.source();
        Point v2 = iseg.target();
        // Interpolate values at the clipped segment endpoints
        float val1 = trilinear_interpolate(adjust_outside_bound_points(v1, grid, v1, v2), grid);
        float val2 = trilinear_interpolate(adjust_outside_bound_points(v2, grid, v1, v2), grid);

        if (!is_bipolar(val1, val2, isovalue)){
            ISO_STATS.line_skip++;
            if(ISO_DBG_ENABLED && iso_dbg_edge_ok(globalEdgeIndex) && !ISO_DBG_ONLY_ERRORS){ std::cerr << "[ISO] LINE non-bipolar edge="<<globalEdgeIndex<<" v=("<<val1<<","<<val2<<")\n"; }
            return; }
        ISO_STATS.line_bip++;
        if(ISO_DBG_ENABLED && iso_dbg_edge_ok(globalEdgeIndex)) std::cerr << "[ISO] LINE bipolar edge="<<globalEdgeIndex<<" dualFacets="<<dualDelaunayFacets.size()<<"\n";

        for (const auto &facet : dualDelaunayFacets)
        {
            int idx1, idx2, idx3, cellIndex1, cellIndex2, cellIndex3;
            bool ok = select_isovertices(voronoiDiagram, facet, globalEdgeIndex,
                                         idx1, idx2, idx3,
                                         cellIndex1, cellIndex2, cellIndex3);
            int iOrient = get_orientation(facet.second, v1, v2, val1, val2);
            bool isValid = ok && (idx1 != idx2 && idx2 != idx3 && idx1 != idx3);
            generate_triangle_multi(iso_surface, idx1, idx2, idx3, iOrient, isValid);
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
    // Iterate with a stable edge index so we can pass the edge id to downstream selectors.
    for (int edgeIdx = 0; edgeIdx < static_cast<int>(voronoiDiagram.edges.size()); ++edgeIdx)
    {
        const auto &edge = voronoiDiagram.edges[edgeIdx];
        Segment3 seg;
        Ray3 ray;
        Line3 line;
        const std::vector<Facet> &dualDelaunayFacets = edge.delaunayFacets;

        if (edge.type == 0)
        {
            // Finite segment: existing path already uses select_isovertices via the segment map.
            process_segment_edge_multi(edge, voronoiDiagram, isovalue, iso_surface);
        }
        else if (edge.type == 1)
        {
            // Ray: use edgeIdx to bind facet -> cell-edge -> local cycle -> isovertex
            CGAL::assign(ray, edge.edgeObject);
            const int source = edge.vertex1; // Voronoi vertex index of the ray source
            process_ray_edge_multi(edgeIdx, source, ray, dualDelaunayFacets,
                                   voronoiDiagram, grid, isovalue, bbox, iso_surface);
        }
        else if (edge.type == 2)
        {
            // Line: same unification as ray
            CGAL::assign(line, edge.edgeObject);
            process_line_edge_multi(edgeIdx, line, dualDelaunayFacets,
                                    voronoiDiagram, grid, isovalue, bbox, iso_surface);
        }
    }
    if(ISO_DBG_ENABLED){ ISO_STATS.dump_summary(); }
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
    if(ISO_DBG_ENABLED && iso_dbg_cell_ok(vc.cellIndex)) {
        std::cerr << "[ISO] cell " << vc.cellIndex << " collecting midpoints; facets=" << vc.facet_indices.size() << "\n";
    }
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

        if(ISO_DBG_ENABLED && iso_dbg_cell_ok(vc.cellIndex) && (!ISO_DBG_ONLY_ERRORS || !current_facet_midpoints.empty())){
            std::cerr << "  [ISO] cell " << vc.cellIndex << " facet#" << facet_index << " midpoints=" << current_facet_midpoints.size() << "\n";
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
    if(ISO_DBG_ENABLED && iso_dbg_cell_ok(vc.cellIndex)) {
        std::cerr << "[ISO] cell " << vc.cellIndex << " connect via global matches\n";
    }
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

    // degree histogram for diagnostics
    if(ISO_DBG_ENABLED){
        int d0=0,d1=0,d2=0,dg=0;
        for(int i=0;i<n;++i){int d=(int)adj[i].size(); if(d==0)++d0; else if(d==1)++d1; else if(d==2)++d2; else ++dg;}
        if(!ISO_DBG_ONLY_ERRORS || (d0||d1||dg)){
            std::cerr << "[ISO] adj degree histogram: deg0="<<d0<<" deg1="<<d1<<" deg2="<<d2<<" deg>=3="<<dg<<"\n";
        }
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

            if (nxt == s) { cycles.push_back(cyc); if(ISO_DBG_ENABLED){ std::cerr << "[ISO] cycle "<< (int)cycles.size()-1 << " len=" << (int)cyc.size() << "\n"; } break; }        // closed
            if (nxt < 0 || nxt >= n || used[nxt]) {                 // broken/self-intersecting
                cycles.push_back(cyc); if(ISO_DBG_ENABLED){ std::cerr << "[ISO] cycle "<< (int)cycles.size()-1 << " len=" << (int)cyc.size() << " (BROKEN/OPEN)\n"; } break;
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
    if(ISO_DBG_ENABLED && iso_dbg_cell_ok(vc.cellIndex)){
        std::cerr << "[ISO] cell " << vc.cellIndex << " cycles=" << cycles.size() << "\n";
    }
    ISO_STATS.cycles_total += cycles.size();
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
        if(ISO_DBG_ENABLED && iso_dbg_cell_ok(vc.cellIndex)){
            int isoIdx = (int)iso_surface.isosurfaceVertices.size()-1;
            std::cerr << "  [ISO] cell " << vc.cellIndex << " isoV@" << isoIdx
                      << " centroid=(" << std::fixed << std::setprecision(6)
                      << cycle.isovertex.x() << "," << cycle.isovertex.y() << "," << cycle.isovertex.z() << ")\n";
        }
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
    ISO_DBG_LOAD_ENV();
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
        ISO_STATS.cells_seen++;
        if (!midpoints.empty()) {
            ISO_STATS.cells_with_midpts++;
            if (num_edges_added == 0) {
                ISO_STATS.cells_with_zero_connect++;
                if(ISO_DBG_ENABLED && iso_dbg_cell_ok(vc.cellIndex)){
                    std::cerr << "[ISO] WARN cell " << vc.cellIndex << " midpoints=" << midpoints.size() << " but 0 connections (check gf matches)\n";
                }
            }
        }
        std::vector<std::vector<int>> cycles;
        extract_cycles(midpoints, cycles);

        compute_cycle_centroids(vc, voronoiDiagram, midpoints, cycles, iso_surface);
    }
}


// ！@brief Wrap up function for constructing iso surface
void construct_iso_surface(Delaunay &dt, VoronoiDiagram &vd, VDC_PARAM &vdc_param, IsoSurface &iso_surface, UnifiedGrid &grid, std::vector<Point> &activeCubeCenters, CGAL::Epick::Iso_cuboid_3 &bbox)
{
    ISO_DBG_LOAD_ENV();
    if(ISO_DBG_ENABLED){ std::cerr << "[ISO] Debug filters: CELL="<<ISO_DBG_FOCUS_CELL<<" GFACET="<<ISO_DBG_FOCUS_GFACET<<" EDGE="<<ISO_DBG_FOCUS_EDGE<<" ONLY_ERRORS="<<(ISO_DBG_ONLY_ERRORS? "1":"0")<<"\n"; }
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
