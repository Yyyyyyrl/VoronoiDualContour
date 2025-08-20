#include "vdc_voronoi.h"


void VoronoiDiagram::create_global_facets() {
    std::map<std::tuple<int, int, int>, std::vector<int>> hashToCellFacets;
    for (size_t fi = 0; fi < facets.size(); ++fi) {
        auto key = getFacetHashKey(facets[fi].vertices_indices);
        hashToCellFacets[key].push_back(static_cast<int>(fi));
    }

    for (const auto& kv : hashToCellFacets) {
        const auto& fvec = kv.second;
        if (fvec.size() > 2) {
            throw std::runtime_error("Facet shared by more than 2 cells.");
        }

        VoronoiFacet vf;
        vf.index = static_cast<int>(global_facets.size());
        int primary = fvec[0];
        vf.vertices_indices = facets[primary].vertices_indices;
        vf.primary_cell_facet_index = primary;
        global_facets.push_back(vf);

        facets[primary].voronoi_facet_index = vf.index;
        facets[primary].orientation = 1;

        if (fvec.size() == 2) {
            int secondary = fvec[1];
            bool opposite = haveOppositeOrientation(facets[primary].vertices_indices, facets[secondary].vertices_indices);
            if (!opposite) {
                throw std::runtime_error("Paired facets " + std::to_string(primary) + " and " + std::to_string(secondary) + " do not have opposite orientations.");
            }
            facets[secondary].voronoi_facet_index = vf.index;
            facets[secondary].orientation = -1;
        }
    }

}

void VoronoiDiagram::compute_bipolar_matches(float isovalue) {
    for (auto& vf : global_facets) {
        // Collect bipolar edges and their types
        std::vector<bool> types;
        std::vector<int> bipolar_edge_indices;
        const auto& vert_indices = vf.vertices_indices;
        int num = vert_indices.size();
        int num_sep_neg = 0;
        int num_sep_pos = 0;
        for (int i = 0; i < num; ++i) {
            int vi0 = vert_indices[i];
            int vi1 = vert_indices[(i + 1) % num];
            float val0 = vertices[vi0].value;
            float val1 = vertices[vi1].value;
            if (is_bipolar(val0, val1, isovalue)) {
                bipolar_edge_indices.push_back(i);  // Index in full facet edges
                if (val0 >= isovalue && val1 < isovalue) {
                    types.push_back(false);  // pos → neg
                    num_sep_neg++;
                } else {
                    types.push_back(true);  // neg → pos
                    num_sep_pos++;
                }
            }
        }
        vf.bipolar_edge_indices = bipolar_edge_indices;  // Store in vf

        // Check for non-alternating types (should always be equal in valid facets)
        if (num_sep_neg != num_sep_pos) {
            std::cerr << "[DEBUG] Facet " << vf.index << ": num_sep_neg=" << num_sep_neg << ", num_sep_pos=" << num_sep_pos << ", values=";
            for (int vi : vf.vertices_indices) std::cerr << " " << vertices[vi].value;
            std::cerr << "\nVertices: ";
            for (int vi : vf.vertices_indices) std::cerr << " " << vertices[vi].coord;
            std::cerr << "\n";  
            throw std::runtime_error("Non-alternating bipolar edge types on facet " + std::to_string(vf.index) + ". This indicates a potential error in facet construction or value interpolation.");
        }

        int num_bipolar = types.size();
        if (num_bipolar == 0 || num_bipolar % 2 != 0) {
            vf.bipolar_match_method = BIPOLAR_MATCH_METHOD::UNDEFINED_MATCH_TYPE;
            vf.bipolar_matches.clear();
            continue;
        }

        // Determine match method (set as member of vf)
        BIPOLAR_MATCH_METHOD match_method = vf.bipolar_match_method;

        // Compute matches using the new routine
        match_facet_bipolar_edges(vf, types, match_method, num_bipolar);
    }
}

void VoronoiDiagram::match_facet_bipolar_edges(VoronoiFacet &vf, const std::vector<bool> &types, BIPOLAR_MATCH_METHOD match_method, int num_bipolar) {

    std::vector<std::pair<int, int>> matches;

    if (match_method == BIPOLAR_MATCH_METHOD::UNCONSTRAINED_MATCH) {
        // Arbitrary consecutive pairing 
        for (int k = 0; k < num_bipolar; k += 2) {
            matches.emplace_back(k, k + 1);
        }
    } else {
        // Standard pairing based on method (connect to next in circular order)
        bool use_pos_neg_start = (match_method == BIPOLAR_MATCH_METHOD::SEP_POS); // Will be true for SEP_POS, false for SEP_NEG
        std::vector<int> starts;
        for (int k = 0; k < num_bipolar; ++k) {
            if (types[k] == use_pos_neg_start) {  
                starts.push_back(k);
            }
        }
        for (int s : starts) {
            int next = (s + 1) % num_bipolar;
            matches.emplace_back(s, next);
        }
    }

    vf.bipolar_matches = matches;
}

std::vector<int> VoronoiDiagram::get_vertices_for_facet(int cell_facet_index) const {
    int vfi = facets[cell_facet_index].voronoi_facet_index;
    std::vector<int> vert = global_facets[vfi].vertices_indices;
    if (facets[cell_facet_index].orientation == -1) {
        std::reverse(vert.begin(), vert.end());
    }
    return vert;
}

std::string matchMethodToString(BIPOLAR_MATCH_METHOD method) {
    switch (method) {
        case BIPOLAR_MATCH_METHOD::SEP_POS: return "SEP_POS";
        case BIPOLAR_MATCH_METHOD::SEP_NEG: return "SEP_NEG";
        case BIPOLAR_MATCH_METHOD::UNCONSTRAINED_MATCH: return "UNCONSTRAINED_MATCH";
        case BIPOLAR_MATCH_METHOD::UNDEFINED_MATCH_TYPE: return "UNDEFINED_MATCH_TYPE";
        default: return "UNKNOWN";
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

//! @brief Processes edges and marks vertices for merging based on distance threshold.
/*!
 * Implements the first phase of edge collapsing by identifying vertices closer than
 * threshold D using a union-find (disjoint-set) data structure.
 *
 * Algorithm phases:
 * 1. Iterate through all finite edges
 * 2. Compute Euclidean distance between vertices
 * 3. Union vertex sets when distance <= D
 *
 * Complexity:
 * - O(E) edge iterations
 * - O(α(V)) per union operation (inverse Ackermann)
 *
 * @param vd Voronoi diagram to process (read-only)
 * @param mapto Union-find structure (initially each vertex points to itself)
 * @param D Distance threshold for merging (must be >= 0)
 * @note Only processes finite edges (skips rays and lines)
 */
static void processEdges(const VoronoiDiagram &vd, std::vector<int> &mapto, double D)
{
    for (size_t edgeIdx = 0; edgeIdx < vd.edges.size(); ++edgeIdx)
    {

        const auto &edge = vd.edges[edgeIdx];

        // Get vertex indices for the current edge
        int v1 = edge.vertex1;
        int v2 = edge.vertex2;

        // Skip edges that are not finite segments (rays or lines)
        if (v1 < 0 || v2 < 0)
            continue;

        // Find the roots of the vertex sets in the union-find structure
        int root1 = v1;
        while (mapto[root1] != root1)
            root1 = mapto[root1];
        int root2 = v2;
        while (mapto[root2] != root2)
            root2 = mapto[root2];

        // If vertices are already in the same set, skip
        if (root1 == root2)
            continue;

        // Compute the distance between the vertex positions
        Point p1 = vd.vertices[v1].coord;
        Point p2 = vd.vertices[v2].coord;
        double dist = CGAL::sqrt(CGAL::squared_distance(p1, p2));

        // Merge vertices if the distance is less than or equal to threshold D
        if (dist <= D)
        {
            // Union by setting the parent; use smaller index as parent for consistency
            if (root1 < root2)
            {
                mapto[root2] = root1;
            }
            else
            {
                mapto[root1] = root2;
            }
        }
    }
}

//! @brief Performs path compression on the union-find structure.
/*!
 * Optimizes the union-find data structure by flattening the tree structure,
 * making future find operations nearly constant time.
 *
 * Algorithm details:
 * - Performs one-pass path compression
 * - Updates all nodes along path to point directly to root
 * - Preserves union-find invariants
 *
 * @param mapto Union-find structure to compress (modified in-place)
 * @note Should be called after all unions are complete
 * @see processEdges()
 */
static void compressMapping(std::vector<int> &mapto)
{
    for (size_t i = 0; i < mapto.size(); ++i)
    {
        int root = i;
        while (mapto[root] != root)
            root = mapto[root];
        mapto[i] = root;
    }
}

//! @brief Rebuilds vertex list after edge collapsing.
/*!
 * Creates a consolidated vertex list containing only representative vertices
 * after edge collapsing, and builds index remapping table.
 *
 * Postconditions:
 * - vd.vertices contains only root vertices
 * - oldToNewVertexIndex maps each original index to its new representative
 * - Vertex count <= original count
 *
 * @param vd Voronoi diagram to modify (vertices updated)
 * @param mapto Union-find structure after path compression
 * @param oldToNewVertexIndex Output mapping (resized and populated)
 * @note Clears and rebuilds the vertexMap spatial index
 */
static void rebuildVertices(VoronoiDiagram &vd, const std::vector<int> &mapto, std::vector<int> &oldToNewVertexIndex)
{
    std::map<int, int> oldToNewIndex;
    std::vector<VoronoiVertex> newVertices;
    for (size_t i = 0; i < vd.vertices.size(); ++i)
    {
        int root = mapto[i];
        if (oldToNewIndex.count(root) == 0)
        {
            VoronoiVertex v = vd.vertices[root];
            v.index = newVertices.size();
            oldToNewIndex[root] = v.index;
            newVertices.push_back(v);
        }
    }
    vd.vertices = std::move(newVertices);

    oldToNewVertexIndex.resize(mapto.size());
    for (size_t i = 0; i < mapto.size(); ++i)
    {
        oldToNewVertexIndex[i] = oldToNewIndex[mapto[i]];
    }
}

//! @brief Rebuilds edge list after vertex collapsing.
/*!
 * Consolidates and deduplicates edges after vertex merging, handling:
 * - Segment edges (finite)
 * - Ray edges (semi-infinite)
 * - Line edges (bi-infinite)
 *
 * Edge processing:
 * 1. Segments: deduplicate using vertex pairs
 * 2. Rays/Lines: deduplicate using geometric equality
 * 3. Updates all vertex references
 *
 * @param vd Voronoi diagram to modify (edges updated)
 * @param oldToNewVertexIndex Vertex index mapping from rebuildVertices()
 * @note Rebuilds segmentVertexPairToEdgeIndex lookup table
 */
static void rebuildEdges(VoronoiDiagram &vd, std::vector<int> &oldToNewVertexIndex)
{
    std::map<std::pair<int, int>, size_t> segmentMap;
    std::map<std::pair<Point, Vector3>, size_t, RayKeyComparator> rayMap;
    std::map<std::pair<Point, Vector3>, size_t, LineKeyComparator> lineMap;

    std::vector<VoronoiEdge> newEdges;

    for (const VoronoiEdge &oldEdge : vd.edges)
    {
        VoronoiEdge newEdge = oldEdge;

        if (oldEdge.type == 0)
        {
            int newV1 = oldToNewVertexIndex[oldEdge.vertex1];
            int newV2 = oldToNewVertexIndex[oldEdge.vertex2];
            if (newV1 == newV2)
                continue;
            int minV = std::min(newV1, newV2);
            int maxV = std::max(newV1, newV2);
            auto it = segmentMap.find({minV, maxV});
            if (it != segmentMap.end())
            {
                size_t keptIdx = it->second;
                newEdges[keptIdx].delaunayFacets.insert(newEdges[keptIdx].delaunayFacets.end(),
                                                        oldEdge.delaunayFacets.begin(), oldEdge.delaunayFacets.end());
            }
            else
            {
                newEdge.vertex1 = newV1;
                newEdge.vertex2 = newV2;
                newEdge.edgeObject = CGAL::make_object(Segment3(vd.vertices[newV1].coord, vd.vertices[newV2].coord));
                newEdges.push_back(newEdge);
                segmentMap[{minV, maxV}] = newEdges.size() - 1;
            }
        }
        else if (oldEdge.type == 1)
        {
            if (oldEdge.vertex1 >= 0)
            {
                int newV1 = oldToNewVertexIndex[oldEdge.vertex1];
                newEdge.vertex1 = newV1;
                newEdge.source = vd.vertices[newV1].coord;
            }
            std::pair<Point, Vector3> key = {newEdge.source, newEdge.direction};
            auto it = rayMap.find(key);
            if (it != rayMap.end())
            {
                size_t keptIdx = it->second;
                newEdges[keptIdx].delaunayFacets.insert(newEdges[keptIdx].delaunayFacets.end(),
                                                        oldEdge.delaunayFacets.begin(), oldEdge.delaunayFacets.end());
            }
            else
            {
                newEdge.edgeObject = CGAL::make_object(Ray3(newEdge.source, Direction3(newEdge.direction)));
                newEdges.push_back(newEdge);
                rayMap[key] = newEdges.size() - 1;
            }
        }
        else if (oldEdge.type == 2)
        {
            Vector3 dir = newEdge.direction;
            if (dir.x() < 0 || (dir.x() == 0 && dir.y() < 0) || (dir.x() == 0 && dir.y() == 0 && dir.z() < 0))
            {
                dir = -dir;
            }
            newEdge.direction = dir;
            std::pair<Point, Vector3> key = {newEdge.source, dir};
            auto it = lineMap.find(key);
            bool merged = false;
            if (it != lineMap.end())
            {
                size_t keptIdx = it->second;
                newEdges[keptIdx].delaunayFacets.insert(newEdges[keptIdx].delaunayFacets.end(),
                                                        oldEdge.delaunayFacets.begin(), oldEdge.delaunayFacets.end());
                merged = true;
            }
            else
            {
                // Check for flipped direction or collinear
                for (auto &existing_entry : lineMap)
                {
                    const auto &exist_key = existing_entry.first;
                    Vector3 exist_dir = exist_key.second;
                    if (!(exist_dir == dir || exist_dir == -dir))
                        continue;
                    Point exist_source = exist_key.first;
                    Vector3 to_new = newEdge.source - exist_source;
                    double proj = CGAL::scalar_product(to_new, exist_dir) / CGAL::scalar_product(exist_dir, exist_dir);
                    Vector3 proj_vec = proj * exist_dir;
                    if ((to_new - proj_vec).squared_length() < 1e-10)
                    { // collinear with tolerance
                        size_t keptIdx = existing_entry.second;
                        newEdges[keptIdx].delaunayFacets.insert(newEdges[keptIdx].delaunayFacets.end(),
                                                                oldEdge.delaunayFacets.begin(), oldEdge.delaunayFacets.end());
                        merged = true;
                        break;
                    }
                }
            }
            if (!merged)
            {
                newEdge.edgeObject = CGAL::make_object(Line3(newEdge.source, newEdge.source + newEdge.direction));
                newEdges.push_back(newEdge);
                lineMap[key] = newEdges.size() - 1;
            }
        }
    }

    vd.edges = std::move(newEdges);

    vd.segmentVertexPairToEdgeIndex.clear();
    for (size_t edgeIdx = 0; edgeIdx < vd.edges.size(); ++edgeIdx)
    {
        const auto &edge = vd.edges[edgeIdx];
        if (edge.type == 0 && edge.vertex1 >= 0 && edge.vertex2 >= 0)
        {
            int v1 = std::min(edge.vertex1, edge.vertex2);
            int v2 = std::max(edge.vertex1, edge.vertex2);
            vd.segmentVertexPairToEdgeIndex[{v1, v2}] = edgeIdx;
        }
    }
}

std::vector<int> update_facet_vertices(const std::vector<int> &old_vertices, const std::vector<int> &oldToNewVertexIndex)
{
    std::vector<int> new_vertices;
    for (int old_v : old_vertices)
    {
        int new_v = oldToNewVertexIndex[old_v];
        if (new_vertices.empty() || new_v != new_vertices.back())
        {
            new_vertices.push_back(new_v);
        }
    }
    if (new_vertices.size() > 1 && new_vertices.front() == new_vertices.back())
    {
        new_vertices.pop_back();
    }
    return new_vertices;
}

// Helper for union-find
static int find(std::vector<int> &mapto, int x)
{
    if (mapto[x] != x)
    {
        mapto[x] = find(mapto, mapto[x]);
    }
    return mapto[x];
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
VoronoiDiagram collapseSmallEdges(const VoronoiDiagram &input_vd, double D, const CGAL::Epick::Iso_cuboid_3 &bbox, Delaunay &dt)
{
    VoronoiDiagram vd = input_vd;
    std::cout << "[DEBUG] Starting collapseSmallEdges with D = " << D << "\n";
    std::cout << "[DEBUG] Initial vertex count: " << vd.vertices.size() << ", edge count: " << vd.edges.size() << "\n";
    std::vector<int> mapto(vd.vertices.size());
    for (size_t i = 0; i < mapto.size(); ++i)
        mapto[i] = i;
    processEdges(vd, mapto, D);
    std::cout << "[DEBUG] Performing path compression\n";
    compressMapping(mapto);
    std::vector<int> oldToNewVertexIndex;
    std::cout << "[DEBUG] Rebuilding vertices\n";
    rebuildVertices(vd, mapto, oldToNewVertexIndex);
    std::cout << "[DEBUG] New vertex count: " << vd.vertices.size() << "\n";
    std::cout << "[DEBUG] Rebuilding edges with duplicate removal\n";

    // Rebuild edges with merging and mapping old->new
    std::vector<VoronoiEdge> newEdges;
    std::map<int, int> mergedEdgeMap;                 // oldIdx -> newIdx, -1 if discarded
    std::map<std::pair<int, int>, int> newSegmentMap; // For new edges

    // For rays: list of <ray, newIdx> for loop check
    std::vector<std::pair<Ray3, int>> rayList;

    // For lines: list of <line, newIdx>
    std::vector<std::pair<Line3, int>> lineList;

    for (size_t oldEdgeIdx = 0; oldEdgeIdx < vd.edges.size(); ++oldEdgeIdx)
    {
        VoronoiEdge edge = vd.edges[oldEdgeIdx]; // Copy to modify
        Segment3 seg;
        Ray3 ray;
        Line3 line;

        if (CGAL::assign(seg, edge.edgeObject))
        {
            int oldV1 = edge.vertex1;
            int oldV2 = edge.vertex2;
            if (oldV1 < 0 || oldV2 < 0)
                continue; // Skip invalid
            int newV1 = oldToNewVertexIndex[oldV1];
            int newV2 = oldToNewVertexIndex[oldV2];
            if (newV1 == newV2)
            {
                mergedEdgeMap[oldEdgeIdx] = -1; // Discarded
                continue;
            }
            int minV = std::min(newV1, newV2);
            int maxV = std::max(newV1, newV2);
            auto it = newSegmentMap.find({minV, maxV});
            if (it != newSegmentMap.end())
            {
                // Merge to existing: append facets
                int existingIdx = it->second;
                newEdges[existingIdx].delaunayFacets.insert(
                    newEdges[existingIdx].delaunayFacets.end(),
                    edge.delaunayFacets.begin(), edge.delaunayFacets.end());
                mergedEdgeMap[oldEdgeIdx] = existingIdx;
            }
            else
            {
                // New edge
                edge.vertex1 = newV1;
                edge.vertex2 = newV2;
                int newIdx = newEdges.size();
                newEdges.push_back(edge);
                newSegmentMap[{minV, maxV}] = newIdx;
                mergedEdgeMap[oldEdgeIdx] = newIdx;
            }
        }
        else if (CGAL::assign(ray, edge.edgeObject))
        {
            // Check for duplicate ray
            bool found = false;
            Point src = ray.source();
            Vector3 dir = normalize_dir(ray.to_vector()); // Use normalized for compare
            for (const auto &pair : rayList)
            {
                const Ray3 &existingRay = pair.first;
                if (approx_equal_points(src, existingRay.source()) &&
                    directions_approx_equal(dir, existingRay.to_vector()))
                {
                    // Merge facets
                    int existingIdx = pair.second;
                    newEdges[existingIdx].delaunayFacets.insert(
                        newEdges[existingIdx].delaunayFacets.end(),
                        edge.delaunayFacets.begin(), edge.delaunayFacets.end());
                    mergedEdgeMap[oldEdgeIdx] = existingIdx;
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                int newIdx = newEdges.size();
                newEdges.push_back(edge);
                rayList.emplace_back(ray, newIdx);
                mergedEdgeMap[oldEdgeIdx] = newIdx;
            }
        }
        else if (CGAL::assign(line, edge.edgeObject))
        {
            // Check for duplicate line
            bool found = false;
            Vector3 dir = normalize_dir(line.to_vector());
            Point pt = line.point(0);
            for (const auto &pair : lineList)
            {
                const Line3 &existingLine = pair.first;
                if (directions_approx_equal(dir, existingLine.to_vector()) &&
                    CGAL::squared_distance(pt, existingLine) < 1e-20)
                { // Point on line
                    // Merge facets
                    int existingIdx = pair.second;
                    newEdges[existingIdx].delaunayFacets.insert(
                        newEdges[existingIdx].delaunayFacets.end(),
                        edge.delaunayFacets.begin(), edge.delaunayFacets.end());
                    mergedEdgeMap[oldEdgeIdx] = existingIdx;
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                int newIdx = newEdges.size();
                newEdges.push_back(edge);
                lineList.emplace_back(line, newIdx);
                mergedEdgeMap[oldEdgeIdx] = newIdx;
            }
        }
    }

    vd.edges = std::move(newEdges);
    vd.segmentVertexPairToEdgeIndex.clear();

    // Rebuild segmentVertexPairToEdgeIndex for new edges
    for (size_t edgeIdx = 0; edgeIdx < vd.edges.size(); ++edgeIdx)
    {
        const auto &edge = vd.edges[edgeIdx];
        if (edge.type == 0 && edge.vertex1 >= 0 && edge.vertex2 >= 0)
        {
            int v1 = std::min(edge.vertex1, edge.vertex2);
            int v2 = std::max(edge.vertex1, edge.vertex2);
            vd.segmentVertexPairToEdgeIndex[{v1, v2}] = edgeIdx;
        }
    }

    std::cout << "[DEBUG] New edge count: " << vd.edges.size() << "\n";
    std::cout << "[DEBUG] Cleared segmentVertexPairToEdgeIndex\n";

    // Rebuild facets (existing code, no change)
    std::vector<VoronoiCellFacet> newFacets;
    std::vector<int> oldToNewFacet(vd.facets.size(), -1);
    for (size_t f = 0; f < vd.facets.size(); ++f)
    {
        VoronoiCellFacet oldF = vd.facets[f];
        std::vector<int> newVertsList;
        for (int oldV : oldF.vertices_indices)
        {
            int newV = oldToNewVertexIndex[oldV];
            newVertsList.push_back(newV);
        }
        // Remove consecutive duplicates
        std::vector<int> cleaned;
        for (size_t k = 0; k < newVertsList.size(); ++k)
        {
            if (k == 0 || newVertsList[k] != newVertsList[k - 1])
            {
                cleaned.push_back(newVertsList[k]);
            }
        }
        // Remove closing duplicate if cyclic
        if (cleaned.size() > 1 && cleaned.front() == cleaned.back())
        {
            cleaned.pop_back();
        }
        if (cleaned.size() < 3)
            continue; // Degenerate: discard
        VoronoiCellFacet newF;
        newF.vertices_indices = std::move(cleaned);
        newF.mirror_facet_index = oldF.mirror_facet_index; // Temporary
        int newFIdx = newFacets.size();
        newFacets.push_back(newF);
        oldToNewFacet[f] = newFIdx;
    }
    // Remap mirror indices
    for (auto &newF : newFacets)
    {
        if (newF.mirror_facet_index >= 0)
        {
            newF.mirror_facet_index = oldToNewFacet[newF.mirror_facet_index];
            if (newF.mirror_facet_index < 0)
                newF.mirror_facet_index = -1; // Mirror was degenerate
        }
    }
    vd.facets = std::move(newFacets);

    // Rebuild cells
    for (auto &cell : vd.cells)
    {
        // Remap vertices (unique)
        std::set<int> newVerts;
        for (int oldV : cell.vertices_indices)
        {
            newVerts.insert(oldToNewVertexIndex[oldV]);
        }
        cell.vertices_indices.assign(newVerts.begin(), newVerts.end());
        // Remap facets
        std::vector<int> newFacetIndices;
        for (int oldF : cell.facet_indices)
        {
            int newF = oldToNewFacet[oldF];
            if (newF >= 0)
                newFacetIndices.push_back(newF);
        }
        cell.facet_indices = std::move(newFacetIndices);
        // Clear polyhedron (can rebuild if needed, but not used post-collapse)
        cell.polyhedron.clear();
    }

    // Update cell edges without full clear
    std::vector<VoronoiCellEdge> updatedCellEdges;
    for (const auto &ce : vd.cellEdges)
    {
        auto it = mergedEdgeMap.find(ce.edgeIndex);
        if (it != mergedEdgeMap.end() && it->second != -1)
        {
            VoronoiCellEdge newCE = ce;
            newCE.edgeIndex = it->second;
            newCE.nextCellEdge = -1; // Reset for re-linking
            updatedCellEdges.push_back(newCE);
        }
        // Discarded if edge collapsed or not mapped
    }

    // Re-group by new edgeIndex and re-link rings
    std::unordered_map<int, std::vector<size_t>> newEdgeToCEIndices;
    for (size_t newCEIdx = 0; newCEIdx < updatedCellEdges.size(); ++newCEIdx)
    {
        int newEdge = updatedCellEdges[newCEIdx].edgeIndex;
        newEdgeToCEIndices[newEdge].push_back(newCEIdx);
    }

    for (auto &kv : newEdgeToCEIndices)
    {
        auto &indices = kv.second;
        size_t M = indices.size();
        if (M < 1)
            continue;
        for (size_t j = 0; j < M; ++j)
        {
            size_t ceIdx = indices[j];
            size_t nextIdx = indices[(j + 1) % M];
            updatedCellEdges[ceIdx].nextCellEdge = nextIdx; // For M==1, next=self
        }
    }

    vd.cellEdges = std::move(updatedCellEdges);

    // Rebuild cellEdgeLookup
    vd.cellEdgeLookup.clear();
    for (size_t ceIdx = 0; ceIdx < vd.cellEdges.size(); ++ceIdx)
    {
        const VoronoiCellEdge &ce = vd.cellEdges[ceIdx];
        vd.cellEdgeLookup[{ce.cellIndex, ce.edgeIndex}] = ceIdx;
    }

    return vd;
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
                      << " edges, but we expected " << ceIndices.size() << ".\n"
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
std::tuple<int, int, int> VoronoiDiagram::getFacetHashKey(const std::vector<int> &verts) const
{
    // pick three smallest indices
    std::array<int, 3> a;
    std::vector<int> tmp = verts;
    std::sort(tmp.begin(), tmp.end());
    a = {tmp[0], tmp[1], tmp[2]};
    return std::make_tuple(a[0], a[1], a[2]);
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
bool VoronoiDiagram::haveSameOrientation(const std::vector<int> &f1, const std::vector<int> &f2) const {
    if (f1.size() != f2.size()) {
        throw std::runtime_error("Facet vertex counts differ: " + std::to_string(f1.size()) + " vs " + std::to_string(f2.size()));
    }
    size_t n = f1.size();
    if (n <= 1) {
        return true;  // Both empty or single element
    }
    
    // Find location of f1[0] in f2
    auto it = std::find(f2.begin(), f2.end(), f1[0]);
    if (it == f2.end()) {
        throw std::runtime_error("f1[0] not found in f2 - facets have different vertex sets");
    }
    size_t iloc = std::distance(f2.begin(), it);
    
    // Check the rest sequentially with wrap-around
    for (size_t i = 1; i < n; ++i) {
        iloc = (iloc + 1) % n;
        if (f1[i] != f2[iloc]) {
            return false;  // Mismatch in sequence (implies different ordering or sets)
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
    std::map<std::tuple<int, int, int>, std::vector<int>> hashToFacets;
    for (size_t fi = 0; fi < facets.size(); ++fi)
    {
        auto key = getFacetHashKey(facets[fi].vertices_indices);
        hashToFacets[key].push_back((int)fi);
        if (hashToFacets[key].size() > 2)
            throw std::runtime_error("Facet with hash “" + std::to_string(std::get<0>(key)) + "," +
                                     std::to_string(std::get<1>(key)) + "," + std::to_string(std::get<2>(key)) +
                                     "” appears in >2 cells.");
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
    std::map<std::tuple<int, int, int>, std::vector<int>> hashToFacets;
    for (size_t fi = 0; fi < facets.size(); ++fi)
    {
        auto key = getFacetHashKey(facets[fi].vertices_indices);
        hashToFacets[key].push_back((int)fi);
    }
    for (auto const &kv : hashToFacets)
    {
        auto const &fvec = kv.second;
        if (fvec.size() == 2)
        {
            auto &A = facets[fvec[0]].vertices_indices;
            auto &B = facets[fvec[1]].vertices_indices;
            if (!haveOppositeOrientation(A, B))
                throw std::runtime_error("Facet “" + std::to_string(fvec[0]) + "” and “" +
                                         std::to_string(fvec[1]) +
                                         "” do not have opposite orientations in their two cells.");
        }
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

std::string get_directory(const std::string& path) {
    size_t pos = path.find_last_of("/\\");
    if (pos == std::string::npos) return "";
    return path.substr(0, pos + 1);
}

std::string get_basename_without_ext(const std::string& path) {
    std::string filename = path.substr(path.find_last_of("/\\") + 1);
    size_t dot_pos = filename.find_last_of('.');
    if (dot_pos == std::string::npos) return filename;
    return filename.substr(0, dot_pos);
}

void write_voronoiDiagram(VoronoiDiagram &vd, std::string &output_filename) {
    std::string dir = get_directory(output_filename);
    std::string base = get_basename_without_ext(output_filename);
    std::string txt_filename = dir + "VoronoiDiagram_" + base + ".txt";

    std::ofstream file(txt_filename);
    if (!file) {
        std::cerr << "Error opening output file: " << txt_filename << "\n";
        exit(EXIT_FAILURE);
    }
    
    file << "VoronoiDiagram:\n";

    file << vd;

    file.close();
    std::cout << "voronoi diagram saved to " << txt_filename << "\n";
}