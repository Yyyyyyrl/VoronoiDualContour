#include "algo/vdc_func.h"

// ===== Debug instrumentation for isosurface pipeline =====
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <limits>
#include <unordered_set>
#include <unordered_map>
#include "core/vdc_debug.h"

static int ISO_DBG_FOCUS_CELL = -1;
static int ISO_DBG_FOCUS_GFACET = -1;
static int ISO_DBG_FOCUS_EDGE = -1;
static bool ISO_DBG_ONLY_ERRORS = false;
static bool ISO_DBG_ENABLED = true; // set via ISO_DBG_LOAD_ENV(), tied to global 'debug'

static int iso_getenv_int(const char *k, int defv)
{
    if (const char *v = std::getenv(k))
    {
        try
        {
            return std::stoi(v);
        }
        catch (...)
        {
        }
    }
    return defv;
}
static bool iso_getenv_bool(const char *k, bool defv)
{
    if (const char *v = std::getenv(k))
    {
        std::string s(v);
        for (char &c : s)
            c = std::tolower(c);
        if (s == "1" || s == "true" || s == "yes" || s == "on")
            return true;
        if (s == "0" || s == "false" || s == "no" || s == "off")
            return false;
    }
    return defv;
}

static void ISO_DBG_LOAD_ENV()
{
    ISO_DBG_FOCUS_CELL = iso_getenv_int("ISO_DEBUG_CELL", -1);
    ISO_DBG_FOCUS_GFACET = iso_getenv_int("ISO_DEBUG_GLOBAL_FACET", -1);
    ISO_DBG_FOCUS_EDGE = iso_getenv_int("ISO_DEBUG_EDGE", -1);
    ISO_DBG_ONLY_ERRORS = iso_getenv_bool("ISO_DEBUG_ONLY_ERRORS", false);
    // Respect global 'debug' flag; allow env var to force off
    ISO_DBG_ENABLED = debug && !iso_getenv_bool("ISO_DEBUG_OFF", false);
}

static inline bool iso_dbg_cell_ok(int c) { return ISO_DBG_FOCUS_CELL < 0 || ISO_DBG_FOCUS_CELL == c; }
static inline bool iso_dbg_gfacet_ok(int g) { return ISO_DBG_FOCUS_GFACET < 0 || ISO_DBG_FOCUS_GFACET == g; }
static inline bool iso_dbg_edge_ok(int e) { return ISO_DBG_FOCUS_EDGE < 0 || ISO_DBG_FOCUS_EDGE == e; }

struct IsoStats
{
    size_t cells_seen = 0, cells_with_midpts = 0, cells_with_zero_connect = 0, cycles_total = 0;
    size_t edges_seg = 0, edges_ray = 0, edges_line = 0;
    size_t seg_bip = 0, seg_skip = 0, ray_bip = 0, ray_skip = 0, line_bip = 0, line_skip = 0;
    size_t tri_ok = 0, tri_bad = 0, tri_dupverts = 0, sel_fail = 0;
    size_t isovertex_clipped = 0;
    double max_clip_distance = 0.0;
    void dump_summary() const
    {
        if (!ISO_DBG_ENABLED)
            return;
        std::cerr << "[ISO] ===== Isosurface Summary =====\n"
                  << "[ISO] Cells: seen=" << cells_seen
                  << " with_midpoints=" << cells_with_midpts
                  << " zero_connections=" << cells_with_zero_connect
                  << " cycles=" << cycles_total << "\n"
                  << "[ISO] Edges: seg=" << edges_seg << " ray=" << edges_ray << " line=" << edges_line << "\n"
                  << "[ISO] Bipolar pass: seg=" << seg_bip << "/" << (seg_bip + seg_skip)
                  << " ray=" << ray_bip << "/" << (ray_bip + ray_skip)
                  << " line=" << line_bip << "/" << (line_bip + line_skip) << "\n"
                  << "[ISO] Triangles: ok=" << tri_ok
                  << " bad=" << tri_bad
                  << " dupVerts=" << tri_dupverts
                  << " isoSelectFail=" << sel_fail << "\n"
                  << "[ISO] Clipping: clipped=" << isovertex_clipped
                  << " max_distance=" << std::fixed << std::setprecision(6) << max_clip_distance << std::defaultfloat << "\n"
                  << "[ISO] Filters: CELL=" << ISO_DBG_FOCUS_CELL
                  << " GFACET=" << ISO_DBG_FOCUS_GFACET
                  << " EDGE=" << ISO_DBG_FOCUS_EDGE
                  << " ONLY_ERRORS=" << (ISO_DBG_ONLY_ERRORS ? "1" : "0") << "\n"
                  << "[ISO] =================================\n";
    }
};
static IsoStats ISO_STATS;
using CycleKey = long long;
using CyclePair = std::pair<CycleKey, CycleKey>;
struct CyclePairHash
{
    size_t operator()(const CyclePair &p) const noexcept
    {
        return std::hash<CycleKey>{}(p.first) ^ (std::hash<CycleKey>{}(p.second) << 1);
    }
};
struct PairBindingInfo
{
    int edgeIndex;
    int count;
    bool foreignConflict;
};
using CyclePairBindingMap = std::unordered_map<CyclePair, PairBindingInfo, CyclePairHash>;

static inline CycleKey make_cycle_key(int cellIndex, int cycleLocal)
{
    if (cellIndex < 0 || cycleLocal < 0)
        return -1;
    return (static_cast<long long>(cellIndex) << 32) ^ static_cast<unsigned int>(cycleLocal);
}

struct EdgeKeyHash
{
    size_t operator()(const std::pair<int, int> &p) const noexcept
    {
        return std::hash<int>{}(p.first) ^ (std::hash<int>{}(p.second) << 1);
    }
};

static void prune_duplicate_edge_triangles(IsoSurface &iso_surface, const VoronoiDiagram &vd)
{
    static const bool modcycDebug = (std::getenv("MODCYC_DEBUG") != nullptr);
    const size_t triCount = iso_surface.isosurfaceTrianglesMulti.size();
    if (triCount == 0)
        return;

    std::vector<char> remove(triCount, 0);
    std::unordered_map<std::pair<int, int>, std::vector<int>, EdgeKeyHash> edgeToTriangles;
    edgeToTriangles.reserve(triCount * 3);

    std::unordered_map<std::pair<int, int>, int, EdgeKeyHash> edgeMultiplicity;
    edgeMultiplicity.reserve(triCount * 3);

    auto normalize_edge = [](int u, int v) -> std::pair<int, int>
    {
        return (u < v) ? std::make_pair(u, v) : std::make_pair(v, u);
    };

    for (size_t idx = 0; idx < triCount; ++idx)
    {
        const auto &tri = iso_surface.isosurfaceTrianglesMulti[idx];
        int a = std::get<0>(tri);
        int b = std::get<1>(tri);
        int c = std::get<2>(tri);

        std::pair<int, int> e1 = normalize_edge(a, b);
        std::pair<int, int> e2 = normalize_edge(b, c);
        std::pair<int, int> e3 = normalize_edge(c, a);

        edgeToTriangles[e1].push_back(static_cast<int>(idx));
        edgeToTriangles[e2].push_back(static_cast<int>(idx));
        edgeToTriangles[e3].push_back(static_cast<int>(idx));

        edgeMultiplicity[e1]++;
        edgeMultiplicity[e2]++;
        edgeMultiplicity[e3]++;
    }

    auto edge_priority = [&](int triIdx) -> std::tuple<int, int, int>
    {
        int edgeId = (triIdx < static_cast<int>(iso_surface.triangleSourceEdges.size())) ? iso_surface.triangleSourceEdges[triIdx] : std::numeric_limits<int>::max();
        int edgeType = 1;
        if (edgeId >= 0 && edgeId < static_cast<int>(vd.edges.size()))
            edgeType = vd.edges[edgeId].type;
        return std::make_tuple(edgeType, edgeId, triIdx);
    };

    for (auto &kv : edgeToTriangles)
    {
        auto &list = kv.second;
        if (list.size() <= 2)
            continue;

        std::vector<int> candidates = list;
        std::sort(candidates.begin(), candidates.end(), [&](int lhs, int rhs)
                  {
                      return edge_priority(lhs) < edge_priority(rhs);
                  });

        // Protect the two best triangles for this edge
        std::unordered_set<int> protectedSet;
        for (size_t i = 0; i < candidates.size() && i < 2; ++i)
            protectedSet.insert(candidates[i]);

        // Attempt to remove remaining triangles starting from worst priority
        for (auto it = candidates.rbegin(); it != candidates.rend(); ++it)
        {
            int triIdx = *it;
            if (protectedSet.count(triIdx) || remove[triIdx])
                continue;

            const auto &tri = iso_surface.isosurfaceTrianglesMulti[triIdx];
            std::pair<int, int> e1 = normalize_edge(std::get<0>(tri), std::get<1>(tri));
            std::pair<int, int> e2 = normalize_edge(std::get<1>(tri), std::get<2>(tri));
            std::pair<int, int> e3 = normalize_edge(std::get<2>(tri), std::get<0>(tri));

            if (edgeMultiplicity[e1] <= 2 || edgeMultiplicity[e2] <= 2 || edgeMultiplicity[e3] <= 2)
                continue;

            remove[triIdx] = 1;
            edgeMultiplicity[e1]--;
            edgeMultiplicity[e2]--;
            edgeMultiplicity[e3]--;
        }
    }

    // Verify no edge dropped below 2; if it did, abort pruning
    for (const auto &kv : edgeMultiplicity)
    {
        if (kv.second < 2)
        {
            std::fill(remove.begin(), remove.end(), 0);
            break;
        }
    }

    size_t removed = 0;
    std::vector<std::tuple<int, int, int>> newTriangles;
    newTriangles.reserve(triCount);
    std::vector<int> newSources;
    newSources.reserve(triCount);

    for (size_t idx = 0; idx < triCount; ++idx)
    {
        if (!remove[idx])
        {
            newTriangles.push_back(iso_surface.isosurfaceTrianglesMulti[idx]);
            if (idx < iso_surface.triangleSourceEdges.size())
                newSources.push_back(iso_surface.triangleSourceEdges[idx]);
        }
        else
        {
            ++removed;
        }
    }

    if (removed == 0)
        return;

    newTriangles.shrink_to_fit();
    newSources.shrink_to_fit();
    iso_surface.isosurfaceTrianglesMulti.swap(newTriangles);
    iso_surface.triangleSourceEdges.swap(newSources);

    if (ISO_DBG_ENABLED)
    {
        std::cerr << "[ISO] post-process removed " << removed << " duplicate triangles\n";
    }

    if (modcycDebug)
    {
        std::unordered_map<std::pair<int, int>, int, EdgeKeyHash> finalCounts;
        for (const auto &tri : iso_surface.isosurfaceTrianglesMulti)
        {
            int a = std::get<0>(tri);
            int b = std::get<1>(tri);
            int c = std::get<2>(tri);
            auto addEdge = [&](int u, int v)
            {
                if (u > v)
                    std::swap(u, v);
                finalCounts[{u, v}]++;
            };
            addEdge(a, b);
            addEdge(b, c);
            addEdge(c, a);
        }
        for (const auto &kvFinal : finalCounts)
        {
            if (kvFinal.second > 2)
            {
                std::cerr << "[MODCYC] remaining overfull edge {" << kvFinal.first.first << "," << kvFinal.first.second
                          << "} count=" << kvFinal.second << "\n";
            }
        }
    }
}

static bool adjust_conflicting_facets(VoronoiDiagram &vd,
                                      const IsoSurface &iso_surface,
                                      std::vector<int> *outTweakedFacets = nullptr)
{
    static const bool modcycDebug = (std::getenv("MODCYC_DEBUG") != nullptr);
    if (iso_surface.isosurfaceTrianglesMulti.empty())
        return false;

    std::vector<std::pair<int, int>> vertexToCellCycle(iso_surface.isosurfaceVertices.size(), {-1, -1});
    for (const auto &cell : vd.cells)
    {
        if (cell.isoVertexStartIndex < 0 || cell.numIsoVertices <= 0)
            continue;
        for (int i = 0; i < cell.numIsoVertices; ++i)
        {
            const int idx = cell.isoVertexStartIndex + i;
            if (idx >= 0 && idx < static_cast<int>(vertexToCellCycle.size()))
                vertexToCellCycle[idx] = {cell.cellIndex, i};
        }
    }

    std::unordered_map<std::pair<int, int>, std::vector<int>, EdgeKeyHash> edgeToTriangles;
    edgeToTriangles.reserve(iso_surface.isosurfaceTrianglesMulti.size() * 3);

    for (size_t triIdx = 0; triIdx < iso_surface.isosurfaceTrianglesMulti.size(); ++triIdx)
    {
        const auto &tri = iso_surface.isosurfaceTrianglesMulti[triIdx];
        int a = std::get<0>(tri);
        int b = std::get<1>(tri);
        int c = std::get<2>(tri);

        auto addEdge = [&](int u, int v)
        {
            if (u > v)
                std::swap(u, v);
            edgeToTriangles[{u, v}].push_back(static_cast<int>(triIdx));
        };

        addEdge(a, b);
        addEdge(b, c);
        addEdge(c, a);
    }

    bool adjusted = false;
    std::unordered_set<int> tweakedFacets;

    for (const auto &kv : edgeToTriangles)
    {
        if (kv.second.size() <= 2)
            continue;

        int vA = kv.first.first;
        int vB = kv.first.second;
        if (vA < 0 || vB < 0 || vA >= static_cast<int>(vertexToCellCycle.size()) || vB >= static_cast<int>(vertexToCellCycle.size()))
            continue;

        const auto &ccA = vertexToCellCycle[vA];
        const auto &ccB = vertexToCellCycle[vB];
        if (ccA.first < 0 || ccB.first < 0)
            continue;

        if (modcycDebug)
        {
            std::cerr << "[MODCYC] conflict edge {" << vA << "," << vB << "} tris=" << kv.second.size()
                      << " cellA=" << ccA.first << ":" << ccA.second
                      << " cellB=" << ccB.first << ":" << ccB.second << "\n";
            for (int triIdx : kv.second)
            {
                const auto &tri = iso_surface.isosurfaceTrianglesMulti[triIdx];
                int edgeId = (triIdx < static_cast<int>(iso_surface.triangleSourceEdges.size())) ? iso_surface.triangleSourceEdges[triIdx] : -1;
                int a = std::get<0>(tri);
                int b = std::get<1>(tri);
                int c = std::get<2>(tri);
                int third = (a != vA && a != vB) ? a : (b != vA && b != vB) ? b : c;
                auto ccThird = (third >= 0 && third < static_cast<int>(vertexToCellCycle.size())) ? vertexToCellCycle[third] : std::pair<int,int>{-1,-1};
                std::cerr << "  tri=" << triIdx << " edgeId=" << edgeId
                          << " third=" << third << " cell=" << ccThird.first << ":" << ccThird.second << "\n";
            }
        }

        for (size_t triIdxPos = 0; triIdxPos < kv.second.size(); ++triIdxPos)
        {
            int triIdx = kv.second[triIdxPos];
            (void)triIdx;
        }

        for (size_t vfi = 0; vfi < vd.global_facets.size(); ++vfi)
        {
            auto &vf = vd.global_facets[vfi];
            if (vf.iso_segments.empty())
                continue;

            for (const auto &seg : vf.iso_segments)
            {
                const int cell0 = vf.incident_cell_indices[0];
                const int cell1 = vf.incident_cell_indices[1];
                const int comp0 = seg.comp[0];
                const int comp1 = seg.comp[1];

                bool match = false;
                if (cell0 == ccA.first && comp0 == ccA.second && cell1 == ccB.first && comp1 == ccB.second)
                    match = true;
                else if (cell0 == ccB.first && comp0 == ccB.second && cell1 == ccA.first && comp1 == ccA.second)
                    match = true;

                if (!match)
                    continue;

                if (!tweakedFacets.insert(static_cast<int>(vfi)).second)
                    continue;

                auto &facet = vd.global_facets[vfi];
                if (modcycDebug)
                {
                    std::cerr << "  [MODCYC] facet " << vfi << " method=" << matchMethodToString(facet.bipolar_match_method)
                              << " slots(" << seg.slotA << ',' << seg.slotB << ")"
                              << " edges(" << seg.edgeA << ',' << seg.edgeB << ")\n";
                }
                BIPOLAR_MATCH_METHOD nextMethod;
                switch (facet.bipolar_match_method)
                {
                case BIPOLAR_MATCH_METHOD::SEP_POS:
                    nextMethod = BIPOLAR_MATCH_METHOD::SEP_NEG;
                    break;
                case BIPOLAR_MATCH_METHOD::SEP_NEG:
                    nextMethod = BIPOLAR_MATCH_METHOD::UNCONSTRAINED_MATCH;
                    break;
                default:
                    nextMethod = BIPOLAR_MATCH_METHOD::UNCONSTRAINED_MATCH;
                    break;
                }

                if (facet.bipolar_match_method != nextMethod)
                {
                    facet.bipolar_match_method = nextMethod;
                    adjusted = true;
                    if (modcycDebug)
                    {
                        std::cerr << "    [MODCYC] facet " << vfi << " switched to " << matchMethodToString(nextMethod) << "\n";
                    }
                    if (outTweakedFacets)
                    {
                        outTweakedFacets->push_back(static_cast<int>(vfi));
                    }
                }

                break;
            }
        }
    }

    return adjusted;
}

static inline bool bind_cycle_pair_to_edge(CyclePairBindingMap &binding,
                                           CycleKey a,
                                           CycleKey b,
                                           int edgeIndex,
                                           bool debugLog,
                                           bool &conflictFlag)
{
    if (a < 0 || b < 0)
        return false;
    if (a == b)
        return true;
    CyclePair key = (a < b) ? CyclePair{a, b} : CyclePair{b, a};
    auto it = binding.find(key);
    if (it == binding.end())
    {
        binding.emplace(key, PairBindingInfo{edgeIndex, 1, false});
        return true;
    }
    if (it->second.edgeIndex != edgeIndex)
    {
        it->second.foreignConflict = true;
        conflictFlag = true;
        if (debugLog)
        {
            std::cerr << "[MODCYC] pair conflict: edge " << edgeIndex
                      << " wants component pair already owned by edge " << it->second.edgeIndex << "\n";
        }
        return true;
    }
    ++it->second.count;
    return true;
}

static void dump_iso_vertex_cycles(const VoronoiDiagram &vd,
                                   const IsoSurface &iso_surface,
                                   const std::string &path)
{
    std::ofstream out(path);
    if (!out)
        return;
    out << "iso_index,cell_index,cycle_local\n";
    for (const auto &cell : vd.cells)
    {
        if (cell.isoVertexStartIndex < 0 || cell.numIsoVertices <= 0)
            continue;
        for (int local = 0; local < cell.numIsoVertices; ++local)
        {
            int idx = cell.isoVertexStartIndex + local;
            if (idx < 0 || idx >= static_cast<int>(iso_surface.isosurfaceVertices.size()))
                continue;
            out << idx << ',' << cell.cellIndex << ',' << local << '\n';
        }
    }
}
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

    auto vertex_is_valid = [&](int idx) -> bool {
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

// Pick an isovertex index for a given (cell, globalEdge).
// Strategy: try exact (cell,edge)→cellEdge→cycle mapping; if none,
// fall back to the first iso-vertex in that cell; else fail (-1).
static inline bool select_isovertex_from_cell_edge(
    const VoronoiDiagram &vd,
    int cellIndex,
    int globalEdgeIndex,
    int &isoIndexOut,
    int &cycleLocalOut)
{
    isoIndexOut = -1;
    cycleLocalOut = -1;

    // Guard: valid cell index
    if (cellIndex < 0 || cellIndex >= static_cast<int>(vd.cells.size()))
    {
        ISO_STATS.sel_fail++;
        return false;
    }

    // 1) Exact (cell,edge) lookup → walk ring to a cellEdge carrying cycles
    auto it = vd.cellEdgeLookup.find(std::make_pair(cellIndex, globalEdgeIndex));
    if (it != vd.cellEdgeLookup.end())
    {
        int ceIdx = it->second;
        const int start = ceIdx;

        while (ceIdx >= 0 &&
               ceIdx < static_cast<int>(vd.cellEdges.size()) &&
               vd.cellEdges[ceIdx].cycleIndices.empty())
        {
            const int nxt = vd.cellEdges[ceIdx].nextCellEdge;
            if (nxt < 0 || nxt == start)
                break;
            ceIdx = nxt;
        }

        if (ceIdx >= 0 &&
            ceIdx < static_cast<int>(vd.cellEdges.size()) &&
            !vd.cellEdges[ceIdx].cycleIndices.empty())
        {
            const int cycLocal = vd.cellEdges[ceIdx].cycleIndices[0];
            const VoronoiCell &vc = vd.cells[cellIndex];

            if (cycLocal >= 0 && cycLocal < vc.numIsoVertices)
            {
                if (ISO_DBG_ENABLED && iso_dbg_cell_ok(cellIndex) && iso_dbg_edge_ok(globalEdgeIndex))
                {
                    std::cerr << "[ISO] pick cell " << cellIndex
                              << " edge " << globalEdgeIndex
                              << " via cellEdge#" << ceIdx
                              << " cycleLocal=" << cycLocal << "\n";
                }
                isoIndexOut = vc.isoVertexStartIndex + cycLocal;
                cycleLocalOut = cycLocal;
                return true;
            }
        }
    }

    if (ISO_DBG_ENABLED && iso_dbg_cell_ok(cellIndex) && iso_dbg_edge_ok(globalEdgeIndex))
    {
        std::cerr << "[ISO] pick cell " << cellIndex
                  << " edge " << globalEdgeIndex
                  << " FAIL (no cycle mapped to edge)\n";
    }
    ISO_STATS.sel_fail++;
    return false;
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
    bool isValid,
    int sourceEdge)
{
    if (!isValid)
    {
        ISO_STATS.tri_bad++;
        if (ISO_DBG_ENABLED && !ISO_DBG_ONLY_ERRORS)
        {
            std::cerr << "[ISO] TRI SKIP invalid indices (" << idx1 << "," << idx2 << "," << idx3 << ")\n";
        }
        return;
    }
    if (idx1 == idx2 || idx2 == idx3 || idx1 == idx3)
    {
        ISO_STATS.tri_dupverts++;
        if (ISO_DBG_ENABLED && !ISO_DBG_ONLY_ERRORS)
        {
            std::cerr << "[ISO] TRI SKIP duplicate indices (" << idx1 << "," << idx2 << "," << idx3 << ")\n";
        }
        return;
    }

    ISO_STATS.tri_ok++;
    if (ISO_DBG_ENABLED && !ISO_DBG_ONLY_ERRORS)
    {
        std::cerr << "[ISO] TRI EMIT (" << idx1 << "," << idx2 << "," << idx3 << ")\n";
    }

    if (iOrient >= 0)
    {
        iso_surface.isosurfaceTrianglesMulti.emplace_back(idx1, idx2, idx3);
    }
    else
    {
        iso_surface.isosurfaceTrianglesMulti.emplace_back(idx1, idx3, idx2);
    }
    iso_surface.triangleSourceEdges.push_back(sourceEdge);
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
    int &cellIndex1, int &cellIndex2, int &cellIndex3,
    int &cycleLocal1, int &cycleLocal2, int &cycleLocal3)
{
    idx1 = idx2 = idx3 = -1;
    cellIndex1 = cellIndex2 = cellIndex3 = -1;
    cycleLocal1 = cycleLocal2 = cycleLocal3 = -1;

    const int iFacet = facet.second;
    Cell_handle c = facet.first;

    const int d1 = (iFacet + 1) % 4;
    const int d2 = (iFacet + 2) % 4;
    const int d3 = (iFacet + 3) % 4;

    Vertex_handle v1 = c->vertex(d1);
    Vertex_handle v2 = c->vertex(d2);
    Vertex_handle v3 = c->vertex(d3);
    
    int b1 = (v1->info().is_dummy) ? 1 : 0;
    int b2 = (v2->info().is_dummy) ? 1 : 0;
    int b3 = (v3->info().is_dummy) ? 1 : 0;

    // Skip facets involving any dummy vertex (cannot map to 3 valid Voronoi cells)
    if (b1 + b2 + b3 >= 1)
    {
        if (ISO_DBG_ENABLED)
        {
            std::cerr << "[ISO] SKIP FACET: contains dummy vertices\n";
            std::cerr << "v1: " << (v1->info().is_dummy ? "dummy" : "real") << ", coords: " << v1->point() << "\n";
            std::cerr << "v2: " << (v2->info().is_dummy ? "dummy" : "real") << ", coords: " << v2->point() << "\n";
            std::cerr << "v3: " << (v3->info().is_dummy ? "dummy" : "real") << ", coords: " << v3->point() << "\n";
        }
        return false;
    }

    cellIndex1 = v1->info().voronoiCellIndex;
    cellIndex2 = v2->info().voronoiCellIndex;
    cellIndex3 = v3->info().voronoiCellIndex;

    select_isovertex_from_cell_edge(voronoiDiagram, cellIndex1, globalEdgeIndex, idx1, cycleLocal1);
    select_isovertex_from_cell_edge(voronoiDiagram, cellIndex2, globalEdgeIndex, idx2, cycleLocal2);
    select_isovertex_from_cell_edge(voronoiDiagram, cellIndex3, globalEdgeIndex, idx3, cycleLocal3);

    if (ISO_DBG_ENABLED && iso_dbg_edge_ok(globalEdgeIndex))
    {
        std::cerr << "[ISO] facet pick edge=" << globalEdgeIndex
                  << " cells=(" << cellIndex1 << "," << cellIndex2 << "," << cellIndex3 << ")"
                  << " iso=(" << idx1 << "," << idx2 << "," << idx3 << ")";
        if (idx1 < 0 || idx2 < 0 || idx3 < 0)
            std::cerr << " [MISS]";
        else if (idx1 == idx2 || idx2 == idx3 || idx1 == idx3)
            std::cerr << " [DUP]";
        else
            std::cerr << " [OK]";
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
    IsoSurface &iso_surface,
    CyclePairBindingMap &pairBinding,
    bool &bindingConflict)
{
    ISO_STATS.edges_seg++;

    int idx_v1 = edge.vertex1;
    int idx_v2 = edge.vertex2;
    //TODO: Figure out why it's happening
    if (idx_v1 > idx_v2){
        std::swap(idx_v1, idx_v2);
    }

    Point v1 = voronoiDiagram.vertices[idx_v1].coord;
    Point v2 = voronoiDiagram.vertices[idx_v2].coord;
    float val1 = voronoiDiagram.vertices[idx_v1].value;
    float val2 = voronoiDiagram.vertices[idx_v2].value;

    if (!is_bipolar(val1, val2, isovalue))
    {
        ISO_STATS.seg_skip++;
        if (ISO_DBG_ENABLED && !ISO_DBG_ONLY_ERRORS)
        {
            std::cerr << "[ISO] SEG non-bipolar v=(" << val1 << "," << val2 << ") iso=" << isovalue << "\n";
        }
        return;
    }

    else
    {
        ISO_STATS.seg_bip++;

        auto itEdge = voronoiDiagram.segmentVertexPairToEdgeIndex.find(std::make_pair(idx_v1, idx_v2));
        if (itEdge == voronoiDiagram.segmentVertexPairToEdgeIndex.end())
            return;

        int globalEdgeIndex = itEdge->second;
        if (ISO_DBG_ENABLED && iso_dbg_edge_ok(globalEdgeIndex)) {
            std::cerr << "[ISO] SEG bipolar edge -> globalEdge=" << globalEdgeIndex << " dualFacets=" << edge.delaunayFacets.size() << "\n";
            std::cerr << voronoiDiagram.edges[globalEdgeIndex];}

        for (const auto &facet : edge.delaunayFacets)
        {
            int idx1, idx2, idx3;
            int cellIndex1, cellIndex2, cellIndex3;
            int cycleLocal1, cycleLocal2, cycleLocal3;
            bool isValid = select_isovertices(voronoiDiagram, facet, globalEdgeIndex,
                                             idx1, idx2, idx3,
                                             cellIndex1, cellIndex2, cellIndex3,
                                             cycleLocal1, cycleLocal2, cycleLocal3);

            // Combinatorial orientation using Delaunay facet index parity
            int iFacet = facet.second;
            int iOrient = get_orientation(iFacet, v1, v2, val1, val2);

            if (isValid)
            {
                CycleKey key1 = make_cycle_key(cellIndex1, cycleLocal1);
                CycleKey key2 = make_cycle_key(cellIndex2, cycleLocal2);
                CycleKey key3 = make_cycle_key(cellIndex3, cycleLocal3);
                bind_cycle_pair_to_edge(pairBinding, key1, key2, globalEdgeIndex, ISO_DBG_ENABLED && iso_dbg_edge_ok(globalEdgeIndex), bindingConflict);
                bind_cycle_pair_to_edge(pairBinding, key2, key3, globalEdgeIndex, ISO_DBG_ENABLED && iso_dbg_edge_ok(globalEdgeIndex), bindingConflict);
                bind_cycle_pair_to_edge(pairBinding, key3, key1, globalEdgeIndex, ISO_DBG_ENABLED && iso_dbg_edge_ok(globalEdgeIndex), bindingConflict);
            }
            generate_triangle_multi(iso_surface, idx1, idx2, idx3, iOrient, isValid, globalEdgeIndex);
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
    IsoSurface &iso_surface,
    CyclePairBindingMap &pairBinding,
    bool &bindingConflict)
{
    ISO_STATS.edges_ray++;
    CGAL::Object intersectObj = CGAL::intersection(bbox, ray);
    Segment3 iseg;
    if (CGAL::assign(iseg, intersectObj))
    {
        Point v1 = ray.source();
        Point v2 = iseg.target();
        // Rays may not carry a valid Voronoi vertex index (source_pt can be -1).
        // Always sample scalar values from the grid at the clipped endpoints.
        float val1 = trilinear_interpolate(adjust_outside_bound_points(v1, grid, v1, v2), grid);
        float val2 = trilinear_interpolate(adjust_outside_bound_points(v2, grid, v1, v2), grid);

        if (!is_bipolar(val1, val2, isovalue))
        {
            ISO_STATS.ray_skip++;
            if (ISO_DBG_ENABLED && iso_dbg_edge_ok(globalEdgeIndex) && !ISO_DBG_ONLY_ERRORS)
            {
                std::cerr << "[ISO] RAY non-bipolar edge=" << globalEdgeIndex << " v=(" << val1 << "," << val2 << ")\n";
            }
            return;
        }
        ISO_STATS.ray_bip++;
        if (ISO_DBG_ENABLED && iso_dbg_edge_ok(globalEdgeIndex)) {
            std::cerr << "[ISO] RAY bipolar edge=" << globalEdgeIndex << " dualFacets=" << dualDelaunayFacets.size() << "\n";
            std::cerr << voronoiDiagram.edges[globalEdgeIndex];}

        for (const auto &facet : dualDelaunayFacets)
        {
            int idx1, idx2, idx3;
            int cellIndex1, cellIndex2, cellIndex3;
            int cycleLocal1, cycleLocal2, cycleLocal3;
            bool ok = select_isovertices(voronoiDiagram, facet, globalEdgeIndex,
                                         idx1, idx2, idx3,
                                         cellIndex1, cellIndex2, cellIndex3,
                                         cycleLocal1, cycleLocal2, cycleLocal3);

            // Combinatorial orientation using Delaunay facet index parity
            int iFacet = facet.second;
            int iOrient = get_orientation(iFacet, v1, v2, val1, val2);
            bool isValid = ok && (idx1 != idx2 && idx2 != idx3 && idx1 != idx3);
            if (isValid)
            {
                CycleKey key1 = make_cycle_key(cellIndex1, cycleLocal1);
                CycleKey key2 = make_cycle_key(cellIndex2, cycleLocal2);
                CycleKey key3 = make_cycle_key(cellIndex3, cycleLocal3);
                bind_cycle_pair_to_edge(pairBinding, key1, key2, globalEdgeIndex, ISO_DBG_ENABLED && iso_dbg_edge_ok(globalEdgeIndex), bindingConflict);
                bind_cycle_pair_to_edge(pairBinding, key2, key3, globalEdgeIndex, ISO_DBG_ENABLED && iso_dbg_edge_ok(globalEdgeIndex), bindingConflict);
                bind_cycle_pair_to_edge(pairBinding, key3, key1, globalEdgeIndex, ISO_DBG_ENABLED && iso_dbg_edge_ok(globalEdgeIndex), bindingConflict);
            }
            generate_triangle_multi(iso_surface, idx1, idx2, idx3, iOrient, isValid, globalEdgeIndex);
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
    IsoSurface &iso_surface,
    CyclePairBindingMap &pairBinding,
    bool &bindingConflict)
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

        if (!is_bipolar(val1, val2, isovalue))
        {
            ISO_STATS.line_skip++;
            if (ISO_DBG_ENABLED && iso_dbg_edge_ok(globalEdgeIndex) && !ISO_DBG_ONLY_ERRORS)
            {
                std::cerr << "[ISO] LINE non-bipolar edge=" << globalEdgeIndex << " v=(" << val1 << "," << val2 << ")\n";
            }
            return;
        }
        ISO_STATS.line_bip++;
        if (ISO_DBG_ENABLED && iso_dbg_edge_ok(globalEdgeIndex))
            std::cerr << "[ISO] LINE bipolar edge=" << globalEdgeIndex << " dualFacets=" << dualDelaunayFacets.size() << "\n";

        for (const auto &facet : dualDelaunayFacets)
        {
            int idx1, idx2, idx3;
            int cellIndex1, cellIndex2, cellIndex3;
            int cycleLocal1, cycleLocal2, cycleLocal3;
            bool ok = select_isovertices(voronoiDiagram, facet, globalEdgeIndex,
                                         idx1, idx2, idx3,
                                         cellIndex1, cellIndex2, cellIndex3,
                                         cycleLocal1, cycleLocal2, cycleLocal3);

            // Combinatorial orientation using Delaunay facet index parity
            int iFacet = facet.second;
            int iOrient = get_orientation(iFacet, v1, v2, val1, val2);
            bool isValid = ok && (idx1 != idx2 && idx2 != idx3 && idx1 != idx3);
            if (isValid)
            {
                CycleKey key1 = make_cycle_key(cellIndex1, cycleLocal1);
                CycleKey key2 = make_cycle_key(cellIndex2, cycleLocal2);
                CycleKey key3 = make_cycle_key(cellIndex3, cycleLocal3);
                bind_cycle_pair_to_edge(pairBinding, key1, key2, globalEdgeIndex, ISO_DBG_ENABLED && iso_dbg_edge_ok(globalEdgeIndex), bindingConflict);
                bind_cycle_pair_to_edge(pairBinding, key2, key3, globalEdgeIndex, ISO_DBG_ENABLED && iso_dbg_edge_ok(globalEdgeIndex), bindingConflict);
                bind_cycle_pair_to_edge(pairBinding, key3, key1, globalEdgeIndex, ISO_DBG_ENABLED && iso_dbg_edge_ok(globalEdgeIndex), bindingConflict);
            }
            generate_triangle_multi(iso_surface, idx1, idx2, idx3, iOrient, isValid, globalEdgeIndex);
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
bool compute_dual_triangles_multi(
    VoronoiDiagram &voronoiDiagram,
    CGAL::Epick::Iso_cuboid_3 &bbox,
    UnifiedGrid &grid,
    float isovalue,
    IsoSurface &iso_surface)
{
    CyclePairBindingMap pairBinding;
    iso_surface.isosurfaceTrianglesMulti.clear();
    iso_surface.triangleSourceEdges.clear();
    bool bindingConflict = false;

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
            process_segment_edge_multi(edge, voronoiDiagram, isovalue, iso_surface, pairBinding, bindingConflict);
        }
        else if (edge.type == 1)
        {
            // Ray: use edgeIdx to bind facet -> cell-edge -> local cycle -> isovertex
            CGAL::assign(ray, edge.edgeObject);
            const int source = edge.vertex1; // Voronoi vertex index of the ray source
            process_ray_edge_multi(edgeIdx, source, ray, dualDelaunayFacets,
                                   voronoiDiagram, grid, isovalue, bbox, iso_surface, pairBinding, bindingConflict);
        }
        else if (edge.type == 2)
        {
            // Line: same unification as ray
            CGAL::assign(line, edge.edgeObject);
            process_line_edge_multi(edgeIdx, line, dualDelaunayFacets,
                                    voronoiDiagram, grid, isovalue, bbox, iso_surface, pairBinding, bindingConflict);
        }
    }
    if (ISO_DBG_ENABLED)
    {
        ISO_STATS.dump_summary();
    }

    prune_duplicate_edge_triangles(iso_surface, voronoiDiagram);
    return bindingConflict;
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
    if (ISO_DBG_ENABLED && iso_dbg_cell_ok(vc.cellIndex))
    {
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

        if (ISO_DBG_ENABLED && iso_dbg_cell_ok(vc.cellIndex) && (!ISO_DBG_ONLY_ERRORS || !current_facet_midpoints.empty()))
        {
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
    if (ISO_DBG_ENABLED && iso_dbg_cell_ok(vc.cellIndex))
    {
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

            if (itA == edge_to_midpoint_index.end() || itB == edge_to_midpoint_index.end())
            {

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
    for (int i = 0; i < n; ++i)
    {
        std::unordered_set<int> uniq(midpoints[i].connected_to.begin(), midpoints[i].connected_to.end());
        adj[i].assign(uniq.begin(), uniq.end());
    }

    // degree histogram for diagnostics
    if (ISO_DBG_ENABLED)
    {
        int d0 = 0, d1 = 0, d2 = 0, dg = 0;
        for (int i = 0; i < n; ++i)
        {
            int d = (int)adj[i].size();
            if (d == 0)
                ++d0;
            else if (d == 1)
                ++d1;
            else if (d == 2)
                ++d2;
            else
                ++dg;
        }
        if (!ISO_DBG_ONLY_ERRORS || (d0 || d1 || dg))
        {
            std::cerr << "[ISO] adj degree histogram: deg0=" << d0 << " deg1=" << d1 << " deg2=" << d2 << " deg>=3=" << dg << "\n";
        }
    }

    // alert if degree != 2
    for (int i = 0; i < n; ++i)
    {
        if (!adj[i].empty() && adj[i].size() != 2 && debug)
        {
            std::cerr << "[warn] midpoint " << i << " has degree " << adj[i].size() << " (expected 2)\n";
        }
    }

    for (int s = 0; s < n; ++s)
    {
        if (used[s] || adj[s].empty())
            continue;
        int prev = -1, cur = s;
        std::vector<int> cyc;
        while (true)
        {
            used[cur] = 1;
            cyc.push_back(cur);

            if (adj[cur].empty())
                break; // broken

            int nxt = (adj[cur].size() == 1) ? adj[cur][0] : (adj[cur][0] == prev ? adj[cur][1] : adj[cur][0]);

            if (nxt == s)
            {
                cycles.push_back(cyc);
                if (ISO_DBG_ENABLED)
                {
                    std::cerr << "[ISO] cycle " << (int)cycles.size() - 1 << " len=" << (int)cyc.size() << "\n";
                }
                break;
            } // closed
            if (nxt < 0 || nxt >= n || used[nxt])
            { // broken/self-intersecting
                cycles.push_back(cyc);
                if (ISO_DBG_ENABLED)
                {
                    std::cerr << "[ISO] cycle " << (int)cycles.size() - 1 << " len=" << (int)cyc.size() << " (BROKEN/OPEN)\n";
                }
                break;
            }

            prev = cur;
            cur = nxt;
        }
    }
}

//! @brief Clips an isovertex to the circumscribed sphere of a cube.
/*!
 * If the isovertex lies outside the circumscribed sphere centered at the cube center,
 * it is projected onto the sphere surface along the direction from center to isovertex.
 *
 * @param isovertex The isosurface vertex (centroid) to potentially clip.
 * @param cube_center The center of the active cube (Delaunay vertex).
 * @param cube_side_length The side length of the active cube.
 * @return The clipped isovertex (same as input if already inside sphere).
 */
static Point clip_isovertex_to_circumscribed_sphere(
    const Point &isovertex,
    const Point &cube_center,
    float cube_side_length)
{
    // Circumscribed sphere radius = (sqrt(3) / 2) * side_length
    const double circumscribed_radius = (std::sqrt(3.0) / 2.0) * cube_side_length;

    // Vector from cube center to isovertex
    Vector3 direction = isovertex - cube_center;
    double distance = std::sqrt(CGAL::to_double(direction.squared_length()));

    // If isovertex is outside the circumscribed sphere, project it onto the sphere
    if (distance > circumscribed_radius)
    {
        // Track clipping statistics
        ISO_STATS.isovertex_clipped++;
        double clip_distance = distance - circumscribed_radius;
        if (clip_distance > ISO_STATS.max_clip_distance)
        {
            ISO_STATS.max_clip_distance = clip_distance;
        }

        // Normalize direction and scale to sphere radius
        Vector3 normalized = direction / distance;
        Point clipped = cube_center + normalized * circumscribed_radius;

        if (ISO_DBG_ENABLED && debug)
        {
            std::cerr << "[ISO] Clipped isovertex: distance=" << distance
                      << " > radius=" << circumscribed_radius
                      << " (clipped by " << clip_distance << ")\n";
        }

        return clipped;
    }

    return isovertex;
}

// TODO: Rename the routine to show its true functionality, as it computes the isovertices AND the cycle centroids.
//! @brief Computes centroids for cycles and updates the isosurface.
/*!
 * Calculates the centroid for each cycle, updates the Voronoi cell's cycles,
 * and adds the centroid to the isosurface vertices. Optionally clips centroids
 * to the circumscribed sphere of the active cube.
 *
 * @param vc The Voronoi cell to update.
 * @param voronoiDiagram The Voronoi diagram for edge lookup.
 * @param midpoints Vector of midpoints used for centroid computation.
 * @param cycles Vector of cycles as lists of midpoint indices.
 * @param iso_surface The isosurface to store vertices.
 * @param cube_side_length The side length of the active cube (0 to disable clipping).
 */
static void compute_cycle_centroids(
    VoronoiCell &vc,
    VoronoiDiagram &voronoiDiagram,
    std::vector<MidpointNode> &midpoints,
    const std::vector<std::vector<int>> &cycles,
    IsoSurface &iso_surface,
    float cube_side_length = 0.0f)
{
    if (ISO_DBG_ENABLED && iso_dbg_cell_ok(vc.cellIndex))
    {
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

        // Clip the centroid to the circumscribed sphere if cube_side_length is specified
        if (cube_side_length > 0.0f)
        {
            Point cube_center = vc.delaunay_vertex->point();
            cycle.isovertex = clip_isovertex_to_circumscribed_sphere(
                cycle.isovertex, cube_center, cube_side_length);
        }

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
        if (ISO_DBG_ENABLED && iso_dbg_cell_ok(vc.cellIndex))
        {
            int isoIdx = (int)iso_surface.isosurfaceVertices.size() - 1;
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
 * centroids as isosurface vertices. Optionally clips centroids to circumscribed spheres.
 *
 * @param voronoiDiagram The Voronoi diagram to compute vertices for.
 * @param isovalue The isovalue to use for vertex computation.
 * @param iso_surface Instance of IsoSurface containing the isosurface vertices and faces.
 * @param grid The grid containing spacing information for clipping (optional).
 */
void compute_isosurface_vertices_multi(VoronoiDiagram &voronoiDiagram, float isovalue, IsoSurface &iso_surface, const UnifiedGrid *grid)
{
    ISO_DBG_LOAD_ENV();
    // Extract cube side length from grid if provided (for clipping)
    float cube_side_length = (grid != nullptr) ? grid->dx : 0.0f;

    // Expect callers to keep voronoiDiagram.global_facets[vfi].bipolar_matches in sync.
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
        if (!midpoints.empty())
        {
            ISO_STATS.cells_with_midpts++;
            if (num_edges_added == 0)
            {
                ISO_STATS.cells_with_zero_connect++;
                if (ISO_DBG_ENABLED && iso_dbg_cell_ok(vc.cellIndex))
                {
                    std::cerr << "[ISO] WARN cell " << vc.cellIndex << " midpoints=" << midpoints.size() << " but 0 connections (check gf matches)\n";
                }
            }
        }
        std::vector<std::vector<int>> cycles;
        extract_cycles(midpoints, cycles);

        compute_cycle_centroids(vc, voronoiDiagram, midpoints, cycles, iso_surface, cube_side_length);
    }
}

// ！@brief Wrap up function for constructing iso surface
void construct_iso_surface(Delaunay &dt, VoronoiDiagram &vd, VDC_PARAM &vdc_param, IsoSurface &iso_surface, UnifiedGrid &grid, std::vector<Point> &activeCubeCenters, CGAL::Epick::Iso_cuboid_3 &bbox, int *out_interior_flips, int *out_boundary_flips, std::size_t *out_clipped_count, double *out_max_clip_distance)
{
    ISO_DBG_LOAD_ENV();
    if (ISO_DBG_ENABLED)
    {
        std::cerr << "[ISO] Debug filters: CELL=" << ISO_DBG_FOCUS_CELL << " GFACET=" << ISO_DBG_FOCUS_GFACET << " EDGE=" << ISO_DBG_FOCUS_EDGE << " ONLY_ERRORS=" << (ISO_DBG_ONLY_ERRORS ? "1" : "0") << "\n";
    }
    std::clock_t start_time = std::clock();
    // Helper used before every new attempt of the multi-isov pipeline. Any facet flip
    // or cycle modification invalidates previously built iso vertices and triangle
    // bookkeeping, so we clear the shared buffers here to guarantee clean state.
    auto reset_iso_accumulators = [&]() {
        iso_surface.isosurfaceVertices.clear();
        iso_surface.isosurfaceTrianglesMulti.clear();
        iso_surface.triangleSourceEdges.clear();
        for (auto &ce : vd.cellEdges)
            ce.cycleIndices.clear();
        for (auto &cell : vd.cells)
        {
            cell.cycles.clear();
            cell.isoVertexStartIndex = -1;
            cell.numIsoVertices = 0;
        }
    };

    ModifyCyclesResult mod_cyc_result;
    bool dualBuilt = false;
    bool bindingConflict = false;

    if (vdc_param.multi_isov)
    {
        const int maxResolve = 8;
        // Track whether facet bipolar matches / iso segments are stale. The very first
        // attempt recomputes everything. Subsequent attempts only touch facets whose
        // match method changed, which keeps the retry cost bounded.
        bool matchesDirty = true;
        bool recomputeAllMatches = true;
        std::vector<int> facetWorklist;
        std::unordered_set<int> facetWorkset;

        // Mark a set of global facets as dirty so we can refresh their matches before
        // the next attempt. We preserve insertion order in facetWorklist to keep the
        // recomputation deterministic while avoiding duplicates with facetWorkset.
        auto mark_facets_dirty = [&](const std::vector<int> &facets) {
            if (facets.empty())
                return;
            matchesDirty = true;
            if (recomputeAllMatches)
                return;
            for (int vfi : facets)
            {
                if (vfi < 0 || vfi >= static_cast<int>(vd.global_facets.size()))
                    continue;
                if (facetWorkset.insert(vfi).second)
                    facetWorklist.push_back(vfi);
            }
        };

        for (int attempt = 0; attempt < maxResolve && !dualBuilt; ++attempt)
        {
            // Refresh bipolar matches and iso segments as required. Either we rebuild
            // the entire structure (first attempt) or only the facets that were flipped
            // during the previous iteration's conflict resolution.
            if (matchesDirty)
            {
                if (recomputeAllMatches)
                {
                    vd.compute_bipolar_matches(vdc_param.isovalue);
                    for (size_t vfi = 0; vfi < vd.global_facets.size(); ++vfi)
                    {
                        build_iso_segments_for_facet(vd, static_cast<int>(vfi), vdc_param.isovalue);
                    }
                }
                else
                {
                    for (int vfi : facetWorklist)
                    {
                        if (vfi < 0 || vfi >= static_cast<int>(vd.global_facets.size()))
                            continue;
                        recompute_bipolar_matches_for_facet(vd, vfi, vdc_param.isovalue);
                        build_iso_segments_for_facet(vd, vfi, vdc_param.isovalue);
                    }
                }
                matchesDirty = false;
                recomputeAllMatches = false;
                facetWorklist.clear();
                facetWorkset.clear();
            }

            if (attempt > 0)
                reset_iso_accumulators();

            // Build cycles and isovertex centroids using the current set of facet matches.
            compute_isosurface_vertices_multi(vd, vdc_param.isovalue, iso_surface, &grid);

            if (vdc_param.mod_cyc)
            {
                std::vector<int> tweakedFacets;
                // Some facets may still disagree after the fresh vertex build. When they
                // do we flip their match policy and loop to rebuild with the updated data.
                if (adjust_conflicting_facets(vd, iso_surface, &tweakedFacets))
                {
                    mark_facets_dirty(tweakedFacets);
                    continue;
                }

                // modify_cycles_pass already recomputes matches for the facets it flips,
                // so we can treat the facet cache as up to date once this returns.
                mod_cyc_result = modify_cycles_pass(vd, vdc_param.isovalue);

                matchesDirty = false; // modify_cycles_pass already recalculates matches locally.
                recomputeAllMatches = false;
                facetWorklist.clear();
                facetWorkset.clear();

                reset_iso_accumulators();
                // Rebuild iso vertices to capture any new cycle assignments produced by
                // modify_cycles_pass before we run the downstream conflict checks.
                compute_isosurface_vertices_multi(vd, vdc_param.isovalue, iso_surface, &grid);

                tweakedFacets.clear();
                if (adjust_conflicting_facets(vd, iso_surface, &tweakedFacets))
                {
                    mark_facets_dirty(tweakedFacets);
                    continue;
                }
            }

            bindingConflict = compute_dual_triangles_multi(vd, bbox, grid, vdc_param.isovalue, iso_surface);
            dualBuilt = true;

            if (bindingConflict && vdc_param.mod_cyc)
            {
                std::vector<int> tweakedFacets;
                // Triangle binding conflicts expose another set of facets whose method
                // must be flipped. We feed them back into the dirty queue and retry.
                if (adjust_conflicting_facets(vd, iso_surface, &tweakedFacets))
                {
                    mark_facets_dirty(tweakedFacets);
                    dualBuilt = false;
                    continue;
                }
            }
        }

        if (!dualBuilt && ISO_DBG_ENABLED)
        {
            std::cerr << "[MODCYC] warning: exceeded conflict resolution iterations; output may still contain duplicates.\n";
        }
    }
    else
    {
        compute_isosurface_vertices_single(grid, vdc_param.isovalue, iso_surface, activeCubeCenters);
        compute_dual_triangles(iso_surface, vd, bbox, dt, grid, vdc_param.isovalue);
        dualBuilt = true;
    }

    if (vdc_param.multi_isov && std::getenv("MODCYC_DEBUG"))
    {
        dump_iso_vertex_cycles(vd, iso_surface, "modcyc_cycles.csv");
    }

    std::clock_t check_time = std::clock();
    double duration = static_cast<double>(check_time - start_time) / CLOCKS_PER_SEC;
    std::cout << "Compute isosurface vertices Execution time: " << duration << " seconds" << std::endl;

    if (vdc_param.multi_isov && !dualBuilt)
    {
        bindingConflict = compute_dual_triangles_multi(vd, bbox, grid, vdc_param.isovalue, iso_surface);
    }
    else if (!vdc_param.multi_isov && !dualBuilt)
    {
        compute_dual_triangles(iso_surface, vd, bbox, dt, grid, vdc_param.isovalue);
    }

    std::clock_t check2_time = std::clock();
    double duration2 = static_cast<double>(check2_time - check_time) / CLOCKS_PER_SEC;
    std::cout << "Compute isosurface facets Execution time: " << duration2 << " seconds" << std::endl;

    if (out_interior_flips)
        *out_interior_flips = mod_cyc_result.interior_flips;
    if (out_boundary_flips)
        *out_boundary_flips = mod_cyc_result.boundary_flips;
    if (out_clipped_count)
        *out_clipped_count = ISO_STATS.isovertex_clipped;
    if (out_max_clip_distance)
        *out_max_clip_distance = ISO_STATS.max_clip_distance;
}
