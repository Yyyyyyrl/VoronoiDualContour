#include "vdc_voronoi.h"

static std::vector<int> collectFacetVoronoiEdges(const VoronoiDiagram &vd, const std::vector<int> &verts)
{
    std::vector<int> out;
    const int n = (int)verts.size();
    out.reserve(n);
    for (int i = 0; i < n; ++i)
    {
        int a = verts[i], b = verts[(i + 1) % n];
        int vmin = std::min(a, b), vmax = std::max(a, b);
        auto it = vd.segmentVertexPairToEdgeIndex.find({vmin, vmax});
        out.push_back(it == vd.segmentVertexPairToEdgeIndex.end() ? -1 : it->second);
    }
    return out;
}

void VoronoiDiagram::create_global_facets()
{
    // Group cell-facets by canonical full-vertex key (order-invariant, no collisions).
    std::map<std::vector<int>, std::vector<int>> keyToCellFacets;

    for (size_t fi = 0; fi < facets.size(); ++fi)
    {
        const auto &F = facets[fi].vertices_indices;
        auto key = getFacetHashKey(F);
        keyToCellFacets[key].push_back(static_cast<int>(fi));
    }

    // Build unique (global) facets from the groups
    for (const auto &kv : keyToCellFacets)
    {
        const auto &fvec = kv.second;
        if (fvec.size() > 2)
        {
            // With a correct key, this should never happen for a proper 3D Voronoi complex.
            // If it does, the input geometry/topology is degenerate or inconsistent.
            std::ostringstream oss;
            oss << "Facet shared by more than 2 cells. Key has " << fvec.size() << " facets; "
                << "vertices: ";
            for (int v : kv.first) oss << v << " ";
            throw std::runtime_error(oss.str());
        }

        VoronoiFacet vf;
        vf.index = static_cast<int>(global_facets.size());

        // Choose a primary representative facet from this group
        const int primary = fvec[0];
        vf.vertices_indices = facets[primary].vertices_indices;

        // Boundary Voronoi edges for this polygon (k -> edge(v[k], v[k+1]))
        vf.voronoi_edge_indices = collectFacetVoronoiEdges(*this, vf.vertices_indices);
        if (vf.voronoi_edge_indices.size() != vf.vertices_indices.size())
        {
            // Defensive: keep sizes consistent; mark missing with -1
            vf.voronoi_edge_indices.resize(vf.vertices_indices.size(), -1);
        }

        vf.primary_cell_facet_index = primary;
        vf.bipolar_match_method = BIPOLAR_MATCH_METHOD::SEP_POS; // default; can be changed later
        global_facets.push_back(vf);

        // Wire the primary back to this global facet
        facets[primary].voronoi_facet_index = vf.index;
        facets[primary].orientation = 1;

        // If there is exactly one mirror, it must have opposite winding
        if (fvec.size() == 2)
        {
            const int secondary = fvec[1];
            const bool opposite =
                haveOppositeOrientation(facets[primary].vertices_indices,
                                        facets[secondary].vertices_indices);
            if (!opposite)
            {
                throw std::runtime_error(
                    "Paired facets " + std::to_string(primary) + " and " +
                    std::to_string(secondary) + " do not have opposite orientations.");
            }
            facets[secondary].voronoi_facet_index = vf.index;
            facets[secondary].orientation = -1;
        }
    }
}

// --- helper: classify and pair facet bipolar edges -------------------------


// Pairs bipolar edges inside a VoronoiFacet according to the method:
//  - SEP_NEG: start from first (+,−) edge, then pair (start,next), (next2,next3), ...
//  - SEP_POS: start from first (−,+) edge, then pair likewise
//  - UNCONSTRAINED_MATCH: pair consecutive in the bipolar list
//
// Output:
//  - vf.bipolar_edge_indices: indices k (edges) around the facet boundary that are bipolar
//  - vf.bipolar_matches: pairs of indices into vf.bipolar_edge_indices (NOT global edge ids)
static void match_facet_bipolar_edges(const VoronoiDiagram &vd,
                                      VoronoiFacet &vf,
                                      float isovalue)
{
    vf.bipolar_edge_indices.clear();
    vf.bipolar_matches.clear();

    const int m = (int)vf.vertices_indices.size();
    if (m < 2) return;

    // Classify edges in the facet’s OWN boundary order
    // edgeType[k] = +1 for (− → +), -1 for (+ → −), 0 otherwise
    std::vector<int> edgeType(m, 0);
    for (int k = 0; k < m; ++k)
    {
        int i0 = vf.vertices_indices[k];
        int i1 = vf.vertices_indices[(k + 1) % m];

        // skip if not a finite segment edge in global map
        if (k < (int)vf.voronoi_edge_indices.size()) {
            int ei = vf.voronoi_edge_indices[k];
            if (ei < 0 || ei >= (int)vd.edges.size() || vd.edges[ei].type != 0) continue;
        }

        float v0 = vd.vertices[i0].value;
        float v1 = vd.vertices[i1].value;
        if (((v0 < isovalue) && (v1 >= isovalue)) || ((v0 >= isovalue) && (v1 < isovalue)))
        {
            edgeType[k] = (((v0 < isovalue) && (v1 >= isovalue)) || ((v0 >= isovalue) && (v1 < isovalue))) ? +1 : -1;
            vf.bipolar_edge_indices.push_back(k); // store boundary slot k
        }
    }

    const int nB = (int)vf.bipolar_edge_indices.size();
    if (nB == 0) return;

    // ensure nB is even
    if (nB % 2) std::cerr << "[warn] facet " << vf.index << " has odd # bipolar edges\n";

    const int want = (vf.bipolar_match_method == BIPOLAR_MATCH_METHOD::SEP_POS) ? +1 : -1;

    // Find anchor: first bipolar slot with desired sign; if missing, fall back to 0
    int startPos = -1;
    for (int i = 0; i < nB; ++i) {
        if (edgeType[vf.bipolar_edge_indices[i]] == want) { startPos = i; break; }
    }
    if (startPos < 0) startPos = 0;

    // Disjoint, circular pairing: (s, s+1), (s+2, s+3), ...
    for (int t = 0; t + 1 < nB; t += 2) {
        int a = vf.bipolar_edge_indices[(startPos + t) % nB];
        int b = vf.bipolar_edge_indices[(startPos + t + 1) % nB];
        vf.bipolar_matches.emplace_back(a, b);
    }
}

void VoronoiDiagram::compute_bipolar_matches(float isovalue)
{
    for (auto &vf : global_facets)
    {
        match_facet_bipolar_edges(*this, vf, isovalue);
    }
}

std::vector<int> VoronoiDiagram::get_vertices_for_facet(int cell_facet_index) const
{
    int vfi = facets[cell_facet_index].voronoi_facet_index;
    std::vector<int> vert = global_facets[vfi].vertices_indices;
    if (facets[cell_facet_index].orientation == -1)
    {
        std::reverse(vert.begin(), vert.end());
    }
    return vert;
}

std::string matchMethodToString(BIPOLAR_MATCH_METHOD method)
{
    switch (method)
    {
    case BIPOLAR_MATCH_METHOD::SEP_POS:
        return "SEP_POS";
    case BIPOLAR_MATCH_METHOD::SEP_NEG:
        return "SEP_NEG";
    case BIPOLAR_MATCH_METHOD::UNCONSTRAINED_MATCH:
        return "UNCONSTRAINED_MATCH";
    case BIPOLAR_MATCH_METHOD::UNDEFINED_MATCH_TYPE:
        return "UNDEFINED_MATCH_TYPE";
    default:
        return "UNKNOWN";
    }
}

namespace
{

    // --------- Small utilities -------------------------------------------------

    struct DSU
    {
        std::vector<int> p, r; // parent, rank
        explicit DSU(int n = 0) : p(n), r(n, 0) { std::iota(p.begin(), p.end(), 0); }
        void reset(int n)
        {
            p.resize(n);
            r.assign(n, 0);
            std::iota(p.begin(), p.end(), 0);
        }
        int find(int x) { return p[x] == x ? x : p[x] = find(p[x]); }
        void unite(int a, int b)
        {
            a = find(a);
            b = find(b);
            if (a == b)
                return;
            if (r[a] < r[b])
                std::swap(a, b);
            p[b] = a;
            if (r[a] == r[b])
                ++r[a];
        }
    };

    inline bool isFiniteSegmentEdge(const VoronoiEdge &e)
    {
        return e.type == 0 && e.vertex1 >= 0 && e.vertex2 >= 0;
    }

    inline double sqr(double x) { return x * x; }

    inline double squaredDist(const Point &a, const Point &b)
    {
        return sqr(a.x() - b.x()) + sqr(a.y() - b.y()) + sqr(a.z() - b.z());
    }

    // Keep representative vertex as the one with smallest original index.  This is
    // stable/deterministic and matches common collapse conventions.
    int chooseRepresentative(const std::vector<int> &group)
    {
        return *std::min_element(group.begin(), group.end());
    }

    // Deduplicate a sequence while keeping the first occurrence order.
    static inline std::vector<int> dedupKeepFirst(const std::vector<int> &seq)
    {
        std::vector<int> out;
        out.reserve(seq.size());
        std::unordered_set<int> seen;
        seen.reserve(seq.size() * 2 + 1);
        for (int v : seq)
            if (!seen.count(v))
            {
                seen.insert(v);
                out.push_back(v);
            }
        return out;
    }

}

//! @brief Adds a vertex to the Voronoi diagram with the given point and value.
/*!
 * Creates a new Voronoi vertex and updates both the vertex list and spatial hash map.
 * Maintains data structure consistency for efficient spatial queries.
 *
 * Implementation details:
 * - Vertex indices are assigned sequentially
 * - Initializes empty cell indices list
 *
 * Thread safety:
 * - Not thread-safe due to shared vertex list and map modifications
 * - External synchronization required for concurrent access
 *
 * @param p The 3D point coordinates for the new vertex (must be finite)
 * @param value The scalar value associated with this vertex (e.g., potential value)
 * @return The index of the newly added vertex (non-negative integer)
 * @throws May throw std::bad_alloc if memory allocation fails
 */
int VoronoiDiagram::AddVertex(const Point &p, float value)
{
    // Create new vertex
    int idx = vertices.size();
    VoronoiVertex v(p);
    v.index = idx;
    v.value = value;
    v.cellIndices = {}; // Initialize empty cell indices

    // Add to vertex list
    vertices.push_back(v);

    return idx;
}

//! @brief Adds a segment edge to the Voronoi diagram.
/*!
 * Manages undirected segment edges between Voronoi vertices with deduplication.
 * Uses sorted vertex indices to ensure consistent representation of undirected edges.
 *
 * Edge types:
 * - Type 0: Finite segment edge (this function)
 * - Type 1: Ray edge (infinite)
 * - Type 2: Line edge (infinite both ways)
 *
 * Data structures:
 * - Maintains segmentVertexPairToEdgeIndex map for O(1) duplicate checking
 * - Edge list stores all edges sequentially
 *
 * @param v1 Index of first vertex (must be valid, >= 0)
 * @param v2 Index of second vertex (must be valid, >= 0)
 * @param seg CGAL Segment_3 object defining the edge geometry
 * @return Index of the edge (existing if duplicate, new otherwise)
 * @throws std::out_of_range if vertex indices are invalid
 */
int VoronoiDiagram::AddSegmentEdge(int v1, int v2, const Segment3 &seg)
{
    // Use sorted vertex indices for undirected edge representation
    int minV = std::min(v1, v2);
    int maxV = std::max(v1, v2);

    // Check if edge already exists
    auto it = segmentVertexPairToEdgeIndex.find({minV, maxV});
    if (it != segmentVertexPairToEdgeIndex.end())
    {
        return it->second; // Return existing edge index
    }

    // Create new edge
    VoronoiEdge edge(CGAL::make_object(seg));
    edge.vertex1 = v1;
    edge.vertex2 = v2;
    edge.type = 0; // Type 0 = segment

    // Add to edge list and update mapping
    int edgeIdx = edges.size();
    edges.push_back(edge);
    segmentVertexPairToEdgeIndex[{minV, maxV}] = edgeIdx;

    return edgeIdx;
}

//! @brief Adds a ray edge to the Voronoi diagram.
/*!
 * Creates a semi-infinite ray edge in the Voronoi diagram.
 * Used for Voronoi edges that extend to infinity from a starting point.
 *
 * Characteristics:
 * - Type 1 edge (ray)
 * - vertex1 and vertex2 set to -1 (no finite vertices)
 * - Stores source point and direction vector
 *
 * Mathematical representation:
 * edge(t) = source + t * direction, t >= 0
 *
 * @param ray CGAL Ray_3 object defining the infinite edge
 * @return Index of the newly created edge
 * @note Ray edges don't participate in spatial hashing
 */
int VoronoiDiagram::AddRayEdge(const Ray3 &ray)
{
    // Create edge object
    VoronoiEdge edge(CGAL::make_object(ray));
    edge.vertex1 = -1; // No start vertex (infinite)
    edge.vertex2 = -1; // No end vertex (infinite)
    edge.type = 1;     // Type 1 = ray
    edge.source = ray.source();
    edge.direction = ray.direction().vector();

    // Add to edge list
    int edgeIdx = edges.size();
    edges.push_back(edge);

    return edgeIdx;
}

//! @brief Adds a line edge to the Voronoi diagram.
/*!
 * Creates a bi-infinite line edge in the Voronoi diagram.
 * Used for Voronoi edges that extend infinitely in both directions.
 *
 * Characteristics:
 * - Type 2 edge (line)
 * - vertex1 and vertex2 set to -1 (no finite vertices)
 * - Stores source point and direction vector
 *
 * Mathematical representation:
 * edge(t) = source + t * direction, t ∈ ℝ
 *
 * @param line CGAL Line_3 object defining the infinite edge
 * @return Index of the newly created edge
 * @note Line edges don't participate in spatial hashing
 */
int VoronoiDiagram::AddLineEdge(const Line3 &line)
{
    // Create edge object
    VoronoiEdge edge(CGAL::make_object(line));
    edge.vertex1 = -1; // No start vertex (infinite)
    edge.vertex2 = -1; // No end vertex (infinite)
    edge.type = 2;     // Type 2 = line
    edge.source = line.point(0);
    edge.direction = line.to_vector();

    // Add to edge list
    int edgeIdx = edges.size();
    edges.push_back(edge);

    return edgeIdx;
}

//! @brief Adds a facet to the Voronoi diagram.
/*!
 * Creates a new polygonal facet (face) in the Voronoi diagram.
 * Facets form the boundaries between Voronoi cells.
 *
 * Topological properties:
 * - Each facet is a simple polygon (no holes)
 * - Vertices must form a closed loop
 * - No duplicate vertices allowed
 * - Minimum 3 vertices required
 *
 * Data management:
 * - Initializes mirror_facet_index to -1 (no mirror)
 * - Facet indices are assigned sequentially
 *
 * @param vertices_indices Vector of vertex indices forming the facet (size >= 3)
 * @return Index of the newly created facet
 * @throws std::invalid_argument if vertices_indices has fewer than 3 elements
 */
int VoronoiDiagram::AddCellFacet(const std::vector<int> &vertices_indices)
{
    // Create new facet
    VoronoiCellFacet facet;
    facet.vertices_indices = vertices_indices;
    facet.mirror_facet_index = -1; // Initialize with no mirror facet

    // Add to facet list
    int facetIdx = facets.size();
    facets.push_back(facet);

    return facetIdx;
}

//! @brief Adds a cell to the Voronoi diagram.
/*!
 * Creates a new Voronoi cell associated with a Delaunay tetrahedralization vertex.
 * Voronoi cells are convex polyhedra in 3D space.
 *
 * Duality properties:
 * - Each Voronoi cell corresponds to exactly one Delaunay vertex
 * - Cell geometry is the dual of the Delaunay tetrahedralization
 *
 * Initial state:
 * - cellIndex set to current cells size
 * - vertices_indices and facet_indices empty
 * - delaunay_vertex stored for reference
 *
 * @param delaunay_vertex Handle to the associated Delaunay vertex
 * @return Index of the newly created cell
 * @note The cell geometry must be populated separately
 */
int VoronoiDiagram::AddCell(Vertex_handle delaunay_vertex)
{
    // Create new cell
    VoronoiCell cell(delaunay_vertex);
    cell.cellIndex = cells.size();

    // Add to cell list
    int cellIdx = cells.size();
    cells.push_back(cell);

    return cellIdx;
}

/*
 * Small Edge Collapsing Routines
 */

// Helper for union-find
static int find(std::vector<int> &mapto, int x)
{
    if (mapto[x] != x)
    {
        mapto[x] = find(mapto, mapto[x]);
    }
    return mapto[x];
}

//! @brief Computes the centroid of a cycle using the associated midpoints.
/*!
 * Calculates the geometric center (isovertex) of a cycle by averaging the positions
 * of its constituent midpoints using CGAL's centroid function.
 *
 * Mathematical properties:
 * - Computes the arithmetic mean of point coordinates
 * - Result is the geometric center of mass
 * - Only valid for non-empty midpoint sets
 *
 * @param midpoints Vector of all MidpointNode objects in the diagram
 * @note Skips computation if midpoint_indices is empty
 * @see Cycle::compute_centroid(const std::vector<VoronoiVertex>&)
 */
void Cycle::compute_centroid(const std::vector<MidpointNode> &midpoints)
{
    if (midpoint_indices.empty())
    {
        return; // No midpoints to compute centroid.
    }

    // Collect the points corresponding to the midpoint indices.
    std::vector<Point> points;
    for (int idx : midpoint_indices)
    {
        points.push_back(midpoints[idx].point);
    }

    // Compute the centroid using CGAL's centroid function.
    isovertex = CGAL::centroid(points.begin(), points.end());
}


bool approx_equal_points(const Point &p1, const Point &p2, double eps_sq = 1e-20)
{
    return CGAL::squared_distance(p1, p2) < eps_sq;
}

Vector3 normalize_dir(const Vector3 &v)
{
    double norm = std::sqrt(v.squared_length());
    if (norm < 1e-10)
        return v;
    Vector3 nv = v / norm;
    double components[3] = {CGAL::to_double(nv.x()), CGAL::to_double(nv.y()), CGAL::to_double(nv.z())};
    for (int i = 0; i < 3; ++i)
    {
        if (std::abs(components[i]) > 1e-10)
        {
            if (components[i] < 0)
            {
                return -nv;
            }
            break;
        }
    }
    return nv;
}

bool directions_approx_equal(const Vector3 &d1, const Vector3 &d2, double eps = 1e-10)
{
    Vector3 n1 = normalize_dir(d1);
    Vector3 n2 = normalize_dir(d2);
    double dot = CGAL::scalar_product(n1, n2);
    return std::abs(dot - 1.0) < eps || std::abs(dot + 1.0) < eps; // Parallel or anti-parallel
}

bool lines_approx_equal(const Line3 &l1, const Line3 &l2, double eps_sq = 1e-20)
{
    if (!directions_approx_equal(l1.to_vector(), l2.to_vector()))
        return false;
    // Check if a point on l1 is on l2 (distance to l2 ==0)
    return CGAL::squared_distance(l1.point(), l2) < eps_sq;
}

// Standalone function to collapse small edges
//! @brief Collapses small edges in a Voronoi diagram.
/*!
 * Processes a Voronoi diagram to merge vertices connected by edges shorter than D.
 * Uses union-find data structure to track vertex merges.
 *
 * @param input_vd Input Voronoi diagram
 * @param D Distance threshold for edge collapsing
 * @param bbox Bounding box of the diagram (unused)
 * @return New Voronoi diagram with small edges collapsed
 */
VoronoiDiagram collapseSmallEdges(const VoronoiDiagram &input_vd,
                                  double D,
                                  const CGAL::Epick::Iso_cuboid_3 & /*bbox*/,
                                  Delaunay & /*dt*/)
{
    // 0) Set up
    VoronoiDiagram out;      // fresh diagram — don't mutate input_vd
    const double D2 = D * D; // squared distances

    const int nV = static_cast<int>(input_vd.vertices.size());
    const int nE = static_cast<int>(input_vd.edges.size());
    const int nC = static_cast<int>(input_vd.cells.size());
    const int nF = static_cast<int>(input_vd.facets.size());

    // 1) Decide merges: union endpoints of every segment edge shorter than D.
    DSU dsu(nV);
    for (int ei = 0; ei < nE; ++ei)
    {
        const VoronoiEdge &e = input_vd.edges[ei];
        if (!isFiniteSegmentEdge(e))
            continue;
        const Point &a = input_vd.vertices[e.vertex1].coord;
        const Point &b = input_vd.vertices[e.vertex2].coord;
        if (squaredDist(a, b) < D2)
        {
            dsu.unite(e.vertex1, e.vertex2); // collapse this short edge
        }
    }

    // 2) Build groups and pick a representative per merged set.
    std::unordered_map<int, std::vector<int>> groups;
    groups.reserve(nV);
    for (int v = 0; v < nV; ++v)
        groups[dsu.find(v)].push_back(v);

    // 3) oldV -> representative oldV (root); then oldV -> newV index mapping
    std::vector<int> oldRoot(nV);
    for (int v = 0; v < nV; ++v)
        oldRoot[v] = dsu.find(v);

    // Determine insertion order for new vertices by ascending representative.
    // Also select a deterministic representative element within each group.
    std::vector<int> reps;
    reps.reserve(groups.size());
    for (auto &kv : groups)
        reps.push_back(chooseRepresentative(kv.second));
    std::sort(reps.begin(), reps.end());

    std::vector<int> oldToNewV(nV, -1);

    // 4) Insert merged vertices into `out`.
    //    represent a merged vertex by the *representative* original
    //    vertex’s coordinate/value, and union of cell membership.
    for (int rep : reps)
    {
        const auto &bucket = groups[dsu.find(rep)];

        // Representative: smallest original index in the bucket
        const int chosen = chooseRepresentative(bucket);
        const VoronoiVertex &origVV = input_vd.vertices[chosen];

        // Merge cell memberships (if used elsewhere) and only keep unique set.
        std::vector<int> mergedCells;
        {
            std::unordered_set<int> s;
            for (int vOld : bucket)
            {
                for (int ci : input_vd.vertices[vOld].cellIndices)
                    s.insert(ci);
            }
            mergedCells.reserve(s.size());
            for (int ci : s)
                mergedCells.push_back(ci);
            std::sort(mergedCells.begin(), mergedCells.end());
        }

        const int newIdx = out.AddVertex(origVV.coord, origVV.value);
        // Preserve back-references to cells
        out.vertices[newIdx].cellIndices = std::move(mergedCells);

        // Map all old vertices in the bucket to this new index
        for (int vOld : bucket)
            oldToNewV[vOld] = newIdx;
    }

    // 5) Rebuild edges — only keep UNcollapsed ones. Also build oldE->newE map.
    //    Preserve and MERGE edge.delaunayFacets across duplicates after collapse.
    std::vector<int> oldToNewE(nE, -1);

    // Helper to append unique Facets (Delaunay::Facet is usually a pair<Cell_handle,int>)
    auto appendUniqueFacets = [](std::vector<Facet> &dst, const std::vector<Facet> &src)
    {
        for (const auto &f : src)
        {
            if (std::find(dst.begin(), dst.end(), f) == dst.end())
                dst.push_back(f);
        }
    };

    // Local cache for segment edges keyed by normalized (vmin,vmax)
    std::map<std::pair<int, int>, int> localSegMap;

    for (int ei = 0; ei < nE; ++ei)
    {
        const VoronoiEdge &e = input_vd.edges[ei];

        if (e.type == 0)
        { // segment
            if (e.vertex1 < 0 || e.vertex2 < 0)
            {
                oldToNewE[ei] = -1;
                continue;
            }
            const int aNew = oldToNewV[e.vertex1];
            const int bNew = oldToNewV[e.vertex2];
            if (aNew == bNew || aNew < 0 || bNew < 0)
            { // collapsed or invalid
                oldToNewE[ei] = -1;
                continue;
            }

            const int vmin = std::min(aNew, bNew);
            const int vmax = std::max(aNew, bNew);
            const auto key = std::make_pair(vmin, vmax);

            // If already have this segment in the new diagram, merge facets
            auto it = localSegMap.find(key);
            if (it != localSegMap.end())
            {
                const int existing = it->second;
                oldToNewE[ei] = existing;
                appendUniqueFacets(out.edges[existing].delaunayFacets, e.delaunayFacets);
                continue;
            }

            // Otherwise create it and copy facets
            const Segment3 seg(out.vertices[aNew].coord, out.vertices[bNew].coord);
            const int ne = out.AddSegmentEdge(aNew, bNew, seg);
            oldToNewE[ei] = ne;

            // Preserve original edge's delaunayFacets
            out.edges[ne].delaunayFacets.clear();
            appendUniqueFacets(out.edges[ne].delaunayFacets, e.delaunayFacets);

            // Maintain both the local and public lookups
            localSegMap[key] = ne;
            out.segmentVertexPairToEdgeIndex[key] = ne;
        }
        else if (e.type == 1)
        { // ray
            Ray3 ray;
            if (!CGAL::assign(ray, e.edgeObject))
                ray = Ray3(e.source, e.direction);
            const int ne = out.AddRayEdge(ray);
            oldToNewE[ei] = ne;

            // Preserve facets for rays
            out.edges[ne].delaunayFacets = e.delaunayFacets;
        }
        else if (e.type == 2)
        { // line
            Line3 line;
            if (!CGAL::assign(line, e.edgeObject))
                line = Line3(e.source, e.direction);
            const int ne = out.AddLineEdge(line);
            oldToNewE[ei] = ne;

            // Preserve facets for lines
            out.edges[ne].delaunayFacets = e.delaunayFacets;
        }
        else
        {
            oldToNewE[ei] = -1; // unknown type
        }
    }

    // 6) Rebuild Cells & Facets with remapped vertex indices. Drop degenerate
    //    facets that end up with < 3 unique vertices after merging.
    std::vector<int> oldToNewCell(nC, -1);
    std::vector<int> oldToNewFacet(nF, -1);

    // clone cells in order to keep indices stable
    for (int ci = 0; ci < nC; ++ci)
    {
        const auto &oldCell = input_vd.cells[ci];
        const int nc = out.AddCell(oldCell.delaunay_vertex);
        oldToNewCell[ci] = nc;

        // Remap the cell’s vertex list ( preserves the order of vertices and also do deduplicate )
        std::vector<int> mappedVerts;
        mappedVerts.reserve(oldCell.vertices_indices.size());
        for (int ov : oldCell.vertices_indices)
        {
            if (ov < 0)
                continue;
            int nv = oldToNewV[ov];
            if (nv >= 0)
                mappedVerts.push_back(nv);
        }
        mappedVerts = dedupKeepFirst(mappedVerts);
        out.cells[nc].vertices_indices = std::move(mappedVerts);

        // Copy scalar/iso bookkeeping (if any)
        out.cells[nc].isoVertexStartIndex = oldCell.isoVertexStartIndex;
        out.cells[nc].numIsoVertices = oldCell.numIsoVertices;
    }

    // Then, rebuild facets in the same order so outside code can keep indices
    for (int fi = 0; fi < nF; ++fi)
    {
        const auto &oldFacet = input_vd.facets[fi];
        std::vector<int> mappedFacetVerts;
        mappedFacetVerts.reserve(oldFacet.vertices_indices.size());
        for (int ov : oldFacet.vertices_indices)
        {
            if (ov < 0)
                continue; // defensive
            const int nv = oldToNewV[ov];
            if (nv >= 0)
                mappedFacetVerts.push_back(nv);
        }
        mappedFacetVerts = dedupKeepFirst(mappedFacetVerts);

        if (mappedFacetVerts.size() < 3)
        {
            // degenerate after collapsing — drop it
            oldToNewFacet[fi] = -1;
            continue;
        }

        const int nf = out.AddCellFacet(mappedFacetVerts);
        oldToNewFacet[fi] = nf;

        // Carry auxiliary fields when present
        out.facets[nf].orientation = input_vd.facets[fi].orientation;
        out.facets[nf].mirror_facet_index = -1;  // will be repaired if needed elsewhere
        out.facets[nf].voronoi_facet_index = -1; // re-created later by create_global_facets()
    }

    // Patch each cell’s facet_indices with the new facet ids, skipping dropped
    // ones; preserve order of the remaining facets.
    for (int ci = 0; ci < nC; ++ci)
    {
        const auto &oldCell = input_vd.cells[ci];
        auto &newCell = out.cells[oldToNewCell[ci]];
        newCell.facet_indices.clear();
        newCell.facet_indices.reserve(oldCell.facet_indices.size());
        for (int of : oldCell.facet_indices)
        {
            if (of < 0 || of >= nF)
                continue;
            const int nf = oldToNewFacet[of];
            if (nf >= 0)
                newCell.facet_indices.push_back(nf);
        }
    }

    // 7) Rebuild VoronoiCellEdges and cellEdgeLookup by remapping & filtering
    //     collapsed edges. Also rebuild the nextCellEdge ring per edge.
    // If input cellEdges are not critical in your downstream pipeline, this
    // section can be simplified safely. Here attempt to preserve order.
    std::vector<int> oldToNewCE;
    oldToNewCE.reserve(input_vd.cellEdges.size());
    for (size_t i = 0; i < input_vd.cellEdges.size(); ++i)
        oldToNewCE.push_back(-1);

    // First pass: insert only those cell-edges whose Voronoi edge survived
    for (size_t cei = 0; cei < input_vd.cellEdges.size(); ++cei)
    {
        const auto &oce = input_vd.cellEdges[cei];
        if (oce.edgeIndex < 0 || oce.edgeIndex >= nE)
            continue;
        const int ne = oldToNewE[oce.edgeIndex];
        if (ne < 0)
            continue; // collapsed away -> skip this cell-edge

        const int nc = (oce.cellIndex >= 0 && oce.cellIndex < nC) ? oldToNewCell[oce.cellIndex] : -1;
        if (nc < 0)
            continue;

        VoronoiCellEdge nce;
        nce.cellIndex = nc;
        nce.edgeIndex = ne;
        nce.cycleIndices.clear(); // cycles remap is optional and pipeline-specific
        nce.nextCellEdge = -1;    // will be wired in the second pass

        const int newIdx = static_cast<int>(out.cellEdges.size());
        out.cellEdges.push_back(nce);
        out.cellEdgeLookup[{nc, ne}] = newIdx;
        oldToNewCE[cei] = newIdx;
    }

    // Second pass: rebuild the per-edge ring (nextCellEdge) keeping original
    // relative ordering as much as possible.
    {
        // Group new cell-edge indices by Voronoi edge index
        std::unordered_map<int, std::vector<int>> edgeToCEs;
        edgeToCEs.reserve(out.cellEdges.size());
        for (int idx = 0; idx < static_cast<int>(out.cellEdges.size()); ++idx)
        {
            edgeToCEs[out.cellEdges[idx].edgeIndex].push_back(idx);
        }
        for (auto &kv : edgeToCEs)
        {
            auto &ring = kv.second;
            // Preserve insertion order
            if (ring.size() >= 1)
            {
                const int m = static_cast<int>(ring.size());
                for (int i = 0; i < m; ++i)
                {
                    const int cur = ring[i];
                    const int nxt = ring[(i + 1) % m];
                    out.cellEdges[cur].nextCellEdge = nxt;
                }
            }
        }
    }

    // 8) Rebuild/refresh global facets & other derived structures.
    out.create_global_facets();
    //out.compute_bipolar_matches(/*isovalue*/ 0.5);

    // 9) Sanity checks
    out.check(false);
    out.checkAdvanced(false);

    return out;
}


//! @brief Checks internal consistency of the VoronoiDiagram.
/*!
 * Performs comprehensive validation of the Voronoi diagram data structure,
 * including topological, geometric and combinatorial properties.
 *
 * Validation stages:
 * 1. Cell-edge lookup consistency
 * 2. Next-cell-edge pointer validity
 * 3. Cell facet relationships
 * 4. Advanced geometric checks
 *
 * @throws std::runtime_error if any inconsistency is detected
 * @note This is an expensive operation (O(V+E+F) time complexity)
 */
void VoronoiDiagram::check(bool check_norm) const
{
    std::cout << "Running validity checks.\n";
    checkCellEdgeLookup();
    checkNextCellEdgeConsistency();
    checkCellFacets();
    checkAdvanced(check_norm);
    std::cout << "VoronoiDiagram::check() passed all advanced checks.\n";
}

//! @brief Verifies that `cellEdgeLookup` matches the data in `cellEdges`.
/*!
 * Validates the bidirectional mapping between (cell,edge) pairs and cellEdges indices.
 * Ensures the lookup table accurately reflects the cellEdges vector contents.
 *
 * Checks performed:
 * 1. All lookup table entries point to valid cellEdges indices
 * 2. The referenced cellEdges entries match the lookup key
 *
 * @throws std::runtime_error if any mapping inconsistency is found
 * @note O(K) time complexity where K is number of cell-edge pairs
 */
void VoronoiDiagram::checkCellEdgeLookup() const
{
    for (const auto &kv : cellEdgeLookup)
    {
        // kv.first is (ic, ie) and kv.second is the index in cellEdges.
        int ic = kv.first.first;  // cellIndex
        int ie = kv.first.second; // edgeIndex
        int cellEdgeIdx = kv.second;

        if (cellEdgeIdx < 0 || cellEdgeIdx >= static_cast<int>(cellEdges.size()))
        {
            throw std::runtime_error("cellEdgeLookup points to invalid VoronoiCellEdge index.");
        }

        const VoronoiCellEdge &ce = cellEdges[cellEdgeIdx];

        if (ce.cellIndex != ic || ce.edgeIndex != ie)
        {
            std::cerr << "ERROR: cellEdgeLookup mismatch!\n";
            std::cerr << "  Lookup says (cell=" << ic << ", edge=" << ie
                      << ") => cellEdgeIdx=" << cellEdgeIdx << "\n";
            std::cerr << "  But VoronoiCellEdge at cellEdgeIdx has (cellIndex="
                      << ce.cellIndex << ", edgeIndex=" << ce.edgeIndex << ")\n";
            throw std::runtime_error("Inconsistent cellEdgeLookup data.");
        }
    }
}

//! @brief Verifies validity of nextCellEdge pointers.
/*!
 * Validates the circular linked list structure of cell edges around each Voronoi edge.
 * Ensures proper connectivity in the cell-edge adjacency graph.
 *
 * Checks performed:
 * 1. All nextCellEdge indices are valid (-1 or within bounds)
 * 2. Linked edges share the same edgeIndex
 * 3. No null pointers in the middle of a cycle
 *
 * @throws std::runtime_error if any pointer inconsistency is found
 * @note Helps maintain the doubly-connected edge list (DCEL) invariants
 */
void VoronoiDiagram::checkNextCellEdgeValidity() const
{
    for (int ceIdx = 0; ceIdx < static_cast<int>(cellEdges.size()); ++ceIdx)
    {
        const VoronoiCellEdge &ce = cellEdges[ceIdx];
        int nxt = ce.nextCellEdge;
        if (nxt < 0)
        {
            continue;
        }
        if (nxt >= static_cast<int>(cellEdges.size()))
        {
            std::cerr << "ERROR: VoronoiCellEdge[" << ceIdx << "].nextCellEdge=" << nxt
                      << " is out of range.\n";
            throw std::runtime_error("Invalid nextCellEdge index.");
        }
        const VoronoiCellEdge &ceNext = cellEdges[nxt];
        if (ceNext.edgeIndex != ce.edgeIndex)
        {
            std::cerr << "ERROR: VoronoiCellEdge[" << ceIdx << "] -> edgeIndex="
                      << ce.edgeIndex << " but nextCellEdge=" << nxt
                      << " has edgeIndex=" << ceNext.edgeIndex << "\n";
            throw std::runtime_error("Inconsistent nextCellEdge edgeIndex.");
        }
    }
}

//! @brief Verifies edge cycles form complete rings.
/*!
 * Validates the topological structure of edge cycles in the Voronoi diagram.
 * Each cycle should form a complete ring around its corresponding edge.
 *
 * Cycle properties verified:
 * 1. All edges with same edgeIndex form a single cycle
 * 2. No premature termination (incomplete cycles)
 * 3. No sub-loops or disconnected components
 *
 * @throws std::runtime_error if cycle structure is invalid
 * @note Essential for maintaining proper cell adjacency relationships
 */
void VoronoiDiagram::checkEdgeCycles() const
{
    std::unordered_map<int, std::vector<int>> edgeIndexToCellEdges;
    for (int ceIdx = 0; ceIdx < static_cast<int>(cellEdges.size()); ++ceIdx)
    {
        int eIdx = cellEdges[ceIdx].edgeIndex;
        edgeIndexToCellEdges[eIdx].push_back(ceIdx);
    }

    for (const auto &kv : edgeIndexToCellEdges)
    {
        const int eIdx = kv.first;
        const auto &ceIndices = kv.second;
        if (ceIndices.empty())
            continue;

        int start = ceIndices[0];
        std::set<int> visited;
        visited.insert(start);

        int current = cellEdges[start].nextCellEdge;
        while (current != start)
        {
            if (current < 0)
            {
                std::cerr << "ERROR: The edges for edgeIndex=" << eIdx
                          << " do not form a complete cycle (nextCellEdge=-1 encountered).\n";
                throw std::runtime_error("Incomplete ring around an edge.");
            }
            if (visited.find(current) != visited.end())
            {
                std::cerr << "ERROR: The edges for edgeIndex=" << eIdx
                          << " contain a sub-loop. Edge " << current
                          << " was already visited.\n";
                throw std::runtime_error("Multiple loops or early cycle detected.");
            }
            visited.insert(current);
            const VoronoiCellEdge &ceNext = cellEdges[current];
            current = ceNext.nextCellEdge;
        }

        if (visited.size() != ceIndices.size())
        {
            std::cerr << "ERROR: For edgeIndex=" << eIdx
                      << ", visited " << visited.size()
                      << " edges, but expected " << ceIndices.size() << ".\n"
                      << "Implying there's a second disconnected cycle or missing edges.\n";
            throw std::runtime_error("Ring does not include all edges for edgeIndex.");
        }
    }
}

//! @brief Verifies consistency of nextCellEdge relationships.
/*!
 * Combines checks for nextCellEdge validity and edge cycles
 * to ensure the complete consistency of edge relationships.
 */
void VoronoiDiagram::checkNextCellEdgeConsistency() const
{
    checkNextCellEdgeValidity();
    checkEdgeCycles();
}

//! @brief Checks each VoronoiCell's facets to ensure that every facet's vertices are in the cell's vertex set.
/*!
 * Validates the vertex-facet relationships within each Voronoi cell.
 * Ensures all facet vertices belong to their containing cell's vertex set.
 *
 * Checks performed:
 * 1. All facet vertex indices are valid
 * 2. Each facet vertex appears in the cell's vertex set
 * 3. No dangling references
 *
 * @throws std::runtime_error if any vertex-facet inconsistency is found
 * @note O(F*V) time complexity in worst case
 */
void VoronoiDiagram::checkCellFacets() const
{
    for (int cIdx = 0; cIdx < static_cast<int>(cells.size()); ++cIdx)
    {
        const VoronoiCell &cell = cells[cIdx];
        // Build a set of the cell's vertex indices for quick membership testing.
        std::set<int> cellVertexSet(cell.vertices_indices.begin(), cell.vertices_indices.end());

        for (int fIdx : cell.facet_indices)
        {
            if (fIdx < 0 || fIdx >= static_cast<int>(facets.size()))
            {
                std::cerr << "ERROR: cell " << cIdx << " has invalid facet index " << fIdx << "\n";
                throw std::runtime_error("Facet index out of range.");
            }

            const VoronoiCellFacet &facet = facets[fIdx];
            for (int vIdx : facet.vertices_indices)
            {
                if (cellVertexSet.find(vIdx) == cellVertexSet.end())
                {
                    std::cerr << "ERROR: VoronoiFacet " << fIdx
                              << " has vertex " << vIdx
                              << " not present in cell[" << cIdx << "]'s vertices.\n";
                    throw std::runtime_error("Inconsistent facet vertex in cell.");
                }
            }
        }
    }
}

// ————————————————————————————————————————————————————————————————
// Helper implementations for checks
// ————————————————————————————————————————————————————————————————

//! @brief Generates a hash key for a facet based on its vertices.
/*!
 * Creates a canonical representation of a facet by selecting three
 * representative vertices in sorted order. Used for facet deduplication
 * and mirror facet identification.
 *
 * Key properties:
 * - Order-invariant (same key for any vertex ordering)
 * - Uses smallest three indices for consistency
 * - Suitable for use as map key
 *
 * @param verts Vector of vertex indices (must have >= 3 elements)
 * @return Tuple of three smallest indices in ascending order
 * @throws std::invalid_argument if verts has fewer than 3 elements
 */
std::vector<int> VoronoiDiagram::getFacetHashKey(const std::vector<int> &verts) const
{
    if (verts.size() < 3)
        throw std::invalid_argument("getFacetHashKey: facet must have at least 3 vertices.");

    // Canonicalize: sort + deduplicate the complete vertex set
    std::vector<int> key = verts;
    std::sort(key.begin(), key.end());
    key.erase(std::unique(key.begin(), key.end()), key.end());

    if (key.size() < 3)
        throw std::runtime_error("getFacetHashKey: facet has fewer than 3 unique vertices after dedup.");

    return key;
}

//! @brief Checks if two facets have the same vertex ordering.
/*!
 * Determines if two facets have identical vertex sequences (considering cyclic shifts).
 * Used to verify proper orientation of mirror facets.
 *
 * Algorithm:
 * - Checks all possible cyclic shifts of f2 against f1
 * - Early termination on first match
 *
 * @param f1 First facet's vertex indices
 * @param f2 Second facet's vertex indices
 * @return true if sequences match under any cyclic shift
 * @note O(n²) time complexity for n-vertex facets
 */
bool VoronoiDiagram::haveSameOrientation(const std::vector<int> &f1, const std::vector<int> &f2) const
{
    if (f1.size() != f2.size())
    {
        throw std::runtime_error("Facet vertex counts differ: " + std::to_string(f1.size()) + " vs " + std::to_string(f2.size()));
    }
    size_t n = f1.size();
    if (n <= 1)
    {
        return true; // Both empty or single element
    }

    // Find location of f1[0] in f2
    auto it = std::find(f2.begin(), f2.end(), f1[0]);
    if (it == f2.end())
    {
        throw std::runtime_error("f1[0] not found in f2 - facets have different vertex sets");
    }
    size_t iloc = std::distance(f2.begin(), it);

    // Check the rest sequentially with wrap-around
    for (size_t i = 1; i < n; ++i)
    {
        iloc = (iloc + 1) % n;
        if (f1[i] != f2[iloc])
        {
            return false; // Mismatch in sequence (implies different ordering or sets)
        }
    }
    return true;
}

//! @brief Checks if two facets have opposite vertex ordering.
bool VoronoiDiagram::haveOppositeOrientation(const std::vector<int> &f1,
                                             const std::vector<int> &f2) const
{
    if (f1.size() != f2.size())
        return false;
    std::vector<int> rev_f2 = f2;
    std::reverse(rev_f2.begin(), rev_f2.end());
    return haveSameOrientation(f1, rev_f2);
}

//! @brief Verifies all facets have at least 3 vertices.
/*!
 * Validates the minimum vertex count requirement for Voronoi facets.
 * Ensures all facets can form proper polygonal faces.
 *
 * Geometric requirements:
 * - Minimum 3 vertices to form a polygon
 * - No degenerate (collinear) vertex sets
 *
 * @throws std::runtime_error if any facet violates minimum vertex count
 * @note Part of the Eulerian polyhedron validation
 */
void VoronoiDiagram::checkFacetVertexCount() const
{
    for (size_t fi = 0; fi < facets.size(); ++fi)
    {
        auto const &F = facets[fi].vertices_indices;
        if (F.size() < 3)
            throw std::runtime_error("Facet " + std::to_string(fi) +
                                     " has fewer than 3 vertices.");
    }
}

//! @brief Verifies all cells have at least 4 facets.
/*!
 * Validates the minimum facet count requirement for Voronoi cells.
 * Ensures all cells can form proper 3D polyhedrons.
 *
 * Topological requirements:
 * - Minimum 4 facets to form a 3D polyhedron
 * - Follows Euler's formula (V-E+F=2)
 *
 * @throws std::runtime_error if any cell violates minimum facet count
 * @note Based on the tetrahedron being the simplest 3D polyhedron
 */
void VoronoiDiagram::checkCellFacetCount() const
{
    for (size_t ci = 0; ci < cells.size(); ++ci)
    {
        if (cells[ci].facet_indices.size() < 4)
            throw std::runtime_error("Cell " + std::to_string(ci) +
                                     " has fewer than 4 facets.");
    }
}

//! @brief Verifies facet sharing between cells.
/*!
 * Validates the proper sharing relationship of facets between adjacent cells.
 * Ensures the Voronoi diagram forms a valid cell complex.
 *
 * Sharing properties:
 * - Each facet shared by exactly 2 cells (in 3D)
 * - No unshared or multiply-shared facets
 *
 * @throws std::runtime_error if facet sharing is invalid
 * @note Uses spatial hashing for efficient facet lookup
 */
void VoronoiDiagram::checkFacetCellCount() const
{
    std::map<std::vector<int>, std::vector<int>> keyToFacets;

    for (size_t fi = 0; fi < facets.size(); ++fi)
    {
        const auto &F = facets[fi].vertices_indices;
        auto key = getFacetHashKey(F);
        auto &bucket = keyToFacets[key];
        bucket.push_back(static_cast<int>(fi));

        if (bucket.size() > 2)
        {
            std::ostringstream oss;
            oss << "Facet key appears in >2 cells (count=" << bucket.size() << "). "
                << "Vertices: ";
            for (int v : key) oss << v << " ";
            throw std::runtime_error(oss.str());
        }
    }
}

//! @brief Verifies edge-facet relationships.
/*!
 * Validates the boundary relationships between edges and facets.
 * Ensures proper manifold structure of the Voronoi diagram.
 *
 * Boundary conditions:
 * - Each edge appears exactly twice (once in each direction)
 * - Forms complete facet boundaries
 * - No dangling edges
 *
 * @throws std::runtime_error if edge-facet relationships are invalid
 * @note Essential for watertight mesh generation
 */
// Modified checkEdgeFacetCount in vdc_voronoi.cpp
// Explanation: Added logging when throwing to print the two facets sharing the edge and their vertex orders. This helps debug which facets have same direction for the edge.
void VoronoiDiagram::checkEdgeFacetCount() const
{
    for (const auto &cell : cells)
    {
        std::map<std::pair<int, int>, int> edge_count;
        for (int facet_idx : cell.facet_indices)
        {
            const auto &vertices = facets[facet_idx].vertices_indices;
            for (size_t i = 0; i < vertices.size(); ++i)
            {
                int u = vertices[i];
                int v = vertices[(i + 1) % vertices.size()];
                std::pair<int, int> edge = {std::min(u, v), std::max(u, v)};
                edge_count[edge]++;
            }
        }
        for (const auto &[edge, count] : edge_count)
        {
            if (count == 2)
            { // Internal edge in the cell (shared by two facets)
                int u = edge.first;
                int v = edge.second;
                int uv_count = 0, vu_count = 0;
                std::vector<int> sharing_facets;
                for (int facet_idx : cell.facet_indices)
                {
                    const auto &vertices = facets[facet_idx].vertices_indices;
                    for (size_t i = 0; i < vertices.size(); ++i)
                    {
                        if (vertices[i] == u && vertices[(i + 1) % vertices.size()] == v)
                            uv_count++;
                        if (vertices[i] == v && vertices[(i + 1) % vertices.size()] == u)
                            vu_count++;
                    }
                    // Collect sharing facets
                    for (size_t i = 0; i < vertices.size(); ++i)
                    {
                        if ((vertices[i] == u && vertices[(i + 1) % vertices.size()] == v) ||
                            (vertices[i] == v && vertices[(i + 1) % vertices.size()] == u))
                        {
                            sharing_facets.push_back(facet_idx);
                            break;
                        }
                    }
                }
                if (uv_count != 1 || vu_count != 1)
                {
                    std::cerr << "[DEBUG] Inconsistent orientations in cell " << cell.cellIndex
                              << " for edge {" << u << "," << v << "}: uv=" << uv_count << ", vu=" << vu_count << "\n";
                    if (sharing_facets.size() == 2)
                    {
                        std::cerr << "[DEBUG] Sharing facets: " << sharing_facets[0] << " and " << sharing_facets[1] << "\n";
                        std::cerr << "[DEBUG] Facet " << sharing_facets[0] << " vertices: ";
                        for (int vi : facets[sharing_facets[0]].vertices_indices)
                            std::cerr << vi << " ";
                        std::cerr << "\n[DEBUG] Facet " << sharing_facets[1] << " vertices: ";
                        for (int vi : facets[sharing_facets[1]].vertices_indices)
                            std::cerr << vi << " ";
                        std::cerr << "\n";
                    }
                    throw std::runtime_error("In cell " + std::to_string(cell.cellIndex) + " edge {" + std::to_string(u) + "," + std::to_string(v) + "} does not have opposite orientations.");
                }
            }
        }
    }
}

//! @brief Verifies facet normals point outward.
/*!
 * Validates the orientation consistency of facet normals relative to cell centers.
 * Ensures proper inside/outside determination for Voronoi cells.
 *
 * Orientation test:
 * - Uses CGAL's orientation predicate
 * - Checks normal points away from Delaunay vertex
 *
 * @throws std::runtime_error if any facet normal is inward-pointing
 * @note Critical for correct geometric computations
 */
void VoronoiDiagram::checkFacetNormals() const
{
    for (const auto &cell : cells)
    {
        auto site = cell.delaunay_vertex->point();
        for (int fi : cell.facet_indices)
        {
            auto const &V = facets[fi].vertices_indices;
            if (V.size() < 3)
                continue;

            // Centroid
            Point centroid(0, 0, 0);
            for (int idx : V)
            {
                centroid = centroid + (vertices[idx].coord - CGAL::ORIGIN) / V.size();
            }

            // calculate normal
            Vector3 normal(0, 0, 0);
            size_t n = V.size();
            for (size_t k = 0; k < n; ++k)
            {
                const Point &p1 = vertices[V[k]].coord;
                const Point &p2 = vertices[V[(k + 1) % n]].coord;
                normal = normal + Vector3(
                                      (p1.y() - p2.y()) * (p1.z() + p2.z()),
                                      (p1.z() - p2.z()) * (p1.x() + p2.x()),
                                      (p1.x() - p2.x()) * (p1.y() + p2.y()));
            }
            normal = normal / 2.0;

            Vector3 v = site - centroid;
            double dot = CGAL::scalar_product(normal, v);
            if (dot > 0 || (dot > -1e-10 && normal.squared_length() < 1e-20))
            { // Allow small negative or zero for degenerates
                std::cerr << "[DEBUG] Facet " << fi << " in cell " << cell.cellIndex << " has dot " << dot << " (should <0), sq_normal " << normal.squared_length() << "\n";
                std::cerr << "[DEBUG] Full vertices: ";
                for (int vi : V)
                    std::cerr << vi << " ";
                std::cerr << "\n[DEBUG] Coords: \n";
                for (int vi : V)
                {
                    const Point &p = vertices[vi].coord;
                    std::cerr << p << "\n";
                }
                std::cerr << "[DEBUG] Centroid: " << centroid << ", site: " << site << "\n";
                if (dot > 0)
                {
                    throw std::runtime_error("Facet " + std::to_string(fi) +
                                             " in cell " + std::to_string(cell.cellIndex) +
                                             " has inward‐pointing normal.");
                }
            }
        }
    }
}

//! @brief Verifies paired facets have opposite orientations.
/*!
 * Validates the orientation consistency of mirror facets.
 * Ensures proper surface normal continuity across cell boundaries.
 *
 * Mirror properties:
 * - Shared facets must have opposite winding orders
 * - Maintains consistent surface orientation
 *
 * @throws std::runtime_error if mirror facets have same orientation
 * @note Preserves manifold property of the diagram
 */
void VoronoiDiagram::checkPairedFacetOrientations() const
{
    std::map<std::vector<int>, std::vector<int>> keyToFacets;

    for (size_t fi = 0; fi < facets.size(); ++fi)
    {
        auto key = getFacetHashKey(facets[fi].vertices_indices);
        keyToFacets[key].push_back(static_cast<int>(fi));
    }

    for (const auto &kv : keyToFacets)
    {
        const auto &fvec = kv.second;

        if (fvec.size() == 2)
        {
            const auto &A = facets[fvec[0]].vertices_indices;
            const auto &B = facets[fvec[1]].vertices_indices;

            if (!haveOppositeOrientation(A, B))
            {
                std::ostringstream oss;
                oss << "Facet \"" << fvec[0] << "\" and \"" << fvec[1]
                    << "\" do not have opposite orientations.";
                throw std::runtime_error(oss.str());
            }
        }
        else if (fvec.size() > 2)
        {
            // Should have been caught by checkFacetCellCount; keep this here defensively.
            std::ostringstream oss;
            oss << "Facet key appears in >2 cells (count=" << fvec.size() << "). "
                << "Vertices: ";
            for (int v : kv.first) oss << v << " ";
            throw std::runtime_error(oss.str());
        }
        // fvec.size()==1 is fine (boundary facet)
    }
}

void VoronoiDiagram::checkAdvanced(bool check_norm) const
{
    checkFacetVertexCount();
    checkCellFacetCount();
    checkFacetCellCount();
    checkEdgeFacetCount(); /**/
    if (check_norm)
    {
        checkFacetNormals();            /**/
        checkPairedFacetOrientations(); /**/
    }
}

std::string get_directory(const std::string &path)
{
    size_t pos = path.find_last_of("/\\");
    if (pos == std::string::npos)
        return "";
    return path.substr(0, pos + 1);
}

std::string get_basename_without_ext(const std::string &path)
{
    std::string filename = path.substr(path.find_last_of("/\\") + 1);
    size_t dot_pos = filename.find_last_of('.');
    if (dot_pos == std::string::npos)
        return filename;
    return filename.substr(0, dot_pos);
}

void write_voronoiDiagram(VoronoiDiagram &vd, std::string &output_filename)
{
    std::string dir = get_directory(output_filename);
    std::string base = get_basename_without_ext(output_filename);
    std::string txt_filename = dir + "VoronoiDiagram_" + base + ".txt";

    std::ofstream file(txt_filename);
    if (!file)
    {
        std::cerr << "Error opening output file: " << txt_filename << "\n";
        exit(EXIT_FAILURE);
    }

    file << "VoronoiDiagram:\n";

    file << vd;

    file.close();
    std::cout << "voronoi diagram saved to " << txt_filename << "\n";
}