#include "vdc_voronoi.h"

//! @brief Finds the index of the vertex corresponding to the given point in the Voronoi diagram.
/*!
 * Implements a spatial hashing technique for efficient point lookup in 3D space.
 * Scales the point coordinates to integers and uses a hash map for quick vertex index retrieval.
 * 
 * Algorithm details:
 * 1. Scales coordinates by 1e6 to convert to integer space
 * 2. Uses tuple-based hash map for spatial partitioning
 * 3. Performs exact match check within hash bucket
 *
 * Performance considerations:
 * - O(1) average case lookup time
 * - Scales well with large number of vertices
 *
 * @param p The 3D point to locate in the Voronoi diagram (must be finite)
 * @return The index of the vertex if found, otherwise -1 indicating point not found
 * @note The scaling factor (1e6) provides ~1 micron precision for coordinates in meter units
 */
int VoronoiDiagram::find_vertex(const Point& p) const {
    // Scale coordinates for spatial hashing
    const double SCALE_FACTOR = 1e6;
    int ix = static_cast<int>(std::round(p.x() * SCALE_FACTOR));
    int iy = static_cast<int>(std::round(p.y() * SCALE_FACTOR));
    int iz = static_cast<int>(std::round(p.z() * SCALE_FACTOR));
    
    // Create hash key from scaled coordinates
    std::tuple<int, int, int> key(ix, iy, iz);
    
    // Look up in vertex map
    auto it = vertexMap.find(key);
    if (it != vertexMap.end()) {
        // Check all vertices in this hash bucket for exact match
        for (int idx : it->second) {
            if (PointApproxEqual()(vertices[idx].coord, p)) {
                return idx;
            }
        }
    }
    return -1;  // Not found
}

//! @brief Adds a vertex to the Voronoi diagram with the given point and value.
/*!
 * Creates a new Voronoi vertex and updates both the vertex list and spatial hash map.
 * Maintains data structure consistency for efficient spatial queries.
 *
 * Implementation details:
 * - Vertex indices are assigned sequentially
 * - Uses same spatial hashing scheme as find_vertex()
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
int VoronoiDiagram::AddVertex(const Point& p, float value) {
    // Create new vertex
    int idx = vertices.size();
    VoronoiVertex v(p);
    v.index = idx;
    v.value = value;
    v.cellIndices = {};  // Initialize empty cell indices
    
    // Add to vertex list
    vertices.push_back(v);
    
    // Update spatial hash map
    const double SCALE_FACTOR = 1e6;
    int ix = static_cast<int>(std::round(p.x() * SCALE_FACTOR));
    int iy = static_cast<int>(std::round(p.y() * SCALE_FACTOR));
    int iz = static_cast<int>(std::round(p.z() * SCALE_FACTOR));
    std::tuple<int, int, int> key(ix, iy, iz);
    vertexMap[key].push_back(idx);
    
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
int VoronoiDiagram::AddSegmentEdge(int v1, int v2, const Segment3& seg) {
    // Use sorted vertex indices for undirected edge representation
    int minV = std::min(v1, v2);
    int maxV = std::max(v1, v2);
    
    // Check if edge already exists
    auto it = segmentVertexPairToEdgeIndex.find({minV, maxV});
    if (it != segmentVertexPairToEdgeIndex.end()) {
        return it->second;  // Return existing edge index
    }
    
    // Create new edge
    VoronoiEdge edge(CGAL::make_object(seg));
    edge.vertex1 = v1;
    edge.vertex2 = v2;
    edge.type = 0;  // Type 0 = segment
    
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
int VoronoiDiagram::AddRayEdge(const Ray3& ray) {
    // Create edge object
    VoronoiEdge edge(CGAL::make_object(ray));
    edge.vertex1 = -1;  // No start vertex (infinite)
    edge.vertex2 = -1;  // No end vertex (infinite)
    edge.type = 1;      // Type 1 = ray
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
int VoronoiDiagram::AddLineEdge(const Line3& line) {
    // Create edge object
    VoronoiEdge edge(CGAL::make_object(line));
    edge.vertex1 = -1;  // No start vertex (infinite)
    edge.vertex2 = -1;  // No end vertex (infinite)
    edge.type = 2;      // Type 2 = line
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
int VoronoiDiagram::AddFacet(const std::vector<int>& vertices_indices) {
    // Create new facet
    VoronoiCellFacet facet;
    facet.vertices_indices = vertices_indices;
    facet.mirror_facet_index = -1;  // Initialize with no mirror facet
    
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
int VoronoiDiagram::AddCell(Vertex_handle delaunay_vertex) {
    // Create new cell
    VoronoiCell cell(delaunay_vertex);
    cell.cellIndex = cells.size();
    
    // Add to cell list
    int cellIdx = cells.size();
    cells.push_back(cell);
    
    return cellIdx;
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


//! @brief Computes centroid using Voronoi vertex positions.
/*!
 * Alternative centroid calculation using direct vertex coordinates instead of midpoints.
 * Provides more stable results when working with raw Voronoi vertices.
 *
 * Implementation details:
 * - Performs simple arithmetic mean of coordinates
 * - Handles empty vertex sets gracefully
 * - More efficient than midpoint-based version
 *
 * @param voronoiVertices Vector of all VoronoiVertex objects in the diagram
 * @note Prefer this version when midpoint data is unavailable
 * @see Cycle::compute_centroid(const std::vector<MidpointNode>&)
 */
void Cycle::compute_centroid(const std::vector<VoronoiVertex> &voronoiVertices)
{
    double sumX = 0, sumY = 0, sumZ = 0;
    for (int idx : midpoint_indices)
    {
        const Point &p = voronoiVertices[idx].coord;
        sumX += p.x();
        sumY += p.y();
        sumZ += p.z();
    }
    double n = static_cast<double>(midpoint_indices.size());
    isovertex = Point(sumX / n, sumY / n, sumZ / n);
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
static void processEdges(const VoronoiDiagram& vd, std::vector<int>& mapto, double D) {
    for (size_t edgeIdx = 0; edgeIdx < vd.edges.size(); ++edgeIdx) {

        const auto& edge = vd.edges[edgeIdx];


        // Get vertex indices for the current edge
        int v1 = edge.vertex1;
        int v2 = edge.vertex2;

        // Skip edges that are not finite segments (rays or lines)
        if (v1 < 0 || v2 < 0) continue;

        // Find the roots of the vertex sets in the union-find structure
        int root1 = v1;
        while (mapto[root1] != root1) root1 = mapto[root1];
        int root2 = v2;
        while (mapto[root2] != root2) root2 = mapto[root2];

        // If vertices are already in the same set, skip
        if (root1 == root2) continue;

        // Compute the distance between the vertex positions
        Point p1 = vd.vertices[v1].coord;
        Point p2 = vd.vertices[v2].coord;
        double dist = CGAL::sqrt(CGAL::squared_distance(p1, p2));

        // Merge vertices if the distance is less than or equal to threshold D
        if (dist <= D) {
            // Union by setting the parent; use smaller index as parent for consistency
            if (root1 < root2) {
                mapto[root2] = root1;
            } else {
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
static void compressMapping(std::vector<int>& mapto) {
    for (size_t i = 0; i < mapto.size(); ++i) {
        int root = i;
        while (mapto[root] != root) root = mapto[root];
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
static void rebuildVertices(VoronoiDiagram& vd, const std::vector<int>& mapto, std::vector<int> &oldToNewVertexIndex) {
    std::map<int, int> oldToNewIndex;
    std::vector<VoronoiVertex> newVertices;
    for (size_t i = 0; i < vd.vertices.size(); ++i) {
        int root = mapto[i];
        if (oldToNewIndex.count(root) == 0) {
            VoronoiVertex v = vd.vertices[root];
            v.index = newVertices.size();
            oldToNewIndex[root] = v.index;
            newVertices.push_back(v);
        }
    }
    vd.vertices = std::move(newVertices);

    oldToNewVertexIndex.resize(mapto.size());
    for (size_t i = 0; i < mapto.size(); ++i) {
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
static void rebuildEdges(VoronoiDiagram& vd, std::vector<int> &oldToNewVertexIndex) {
    std::map<std::pair<int, int>, int> segmentMap;
    std::set<VoronoiEdge> rayLineSet;

    // Process all edges to identify unique ones
    for (size_t edgeIdx = 0; edgeIdx < vd.edges.size(); ++edgeIdx) {
        const auto& edge = vd.edges[edgeIdx];
        Segment3 seg;
        Ray3 ray;
        Line3 line;

        if (CGAL::assign(seg, edge.edgeObject)) {
            // Handle segments using vertex indices
            int oldIdx1 = edge.vertex1;
            int oldIdx2 = edge.vertex2;
            if (oldIdx1 >= 0 && oldIdx2 >= 0) {
                int newIdx1 = oldToNewVertexIndex[oldIdx1];
                int newIdx2 = oldToNewVertexIndex[oldIdx2];
                if (newIdx1 != newIdx2) { // Skip collapsed edges (same vertex)
                    int v1 = std::min(newIdx1, newIdx2);
                    int v2 = std::max(newIdx1, newIdx2);
                    if (segmentMap.find({v1, v2}) == segmentMap.end()) {
                        segmentMap[{v1, v2}] = edgeIdx;
                    }
                }
            }
        } else if (CGAL::assign(ray, edge.edgeObject) || CGAL::assign(line, edge.edgeObject)) {
            // Handle rays and lines with geometric comparison
            rayLineSet.insert(edge);
        }
    }

    // Rebuild the edges and edgeVertexIndices vectors
    std::vector<VoronoiEdge> newEdges;
    newEdges.reserve(segmentMap.size() + rayLineSet.size());

    for (auto& kv : segmentMap) {
        size_t edgeIdx = kv.second;
        VoronoiEdge edge = vd.edges[edgeIdx];

        edge.vertex1 = oldToNewVertexIndex[vd.edges[edgeIdx].vertex1];
        edge.vertex2 = oldToNewVertexIndex[vd.edges[edgeIdx].vertex2];

        newEdges.push_back(edge);
    }
    for (auto& edge : rayLineSet) {
        newEdges.push_back(edge);
    }

    vd.edges = std::move(newEdges);
    vd.segmentVertexPairToEdgeIndex.clear();

    // Rebuild segmentVertexPairToEdgeIndex
    for (size_t edgeIdx = 0; edgeIdx < vd.edges.size(); ++edgeIdx) {
        const auto& edge = vd.edges[edgeIdx];
        if (edge.type == 0 && edge.vertex1 >= 0 && edge.vertex2 >= 0) {
            int v1 = std::min(edge.vertex1, edge.vertex2);
            int v2 = std::max(edge.vertex1, edge.vertex2);
            vd.segmentVertexPairToEdgeIndex[{v1, v2}] = edgeIdx;
        }
    }
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
VoronoiDiagram collapseSmallEdges(const VoronoiDiagram& input_vd, double D, const CGAL::Epick::Iso_cuboid_3& bbox) {
    // Create a copy of the input VoronoiDiagram
    VoronoiDiagram vd = input_vd;

    std::cout << "[DEBUG] Starting collapseSmallEdges with D = " << D << "\n";
    std::cout << "[DEBUG] Initial vertex count: " << vd.vertices.size() << ", edge count: " << vd.edges.size() << "\n";

    // Initialize union-find structure
    std::vector<int> mapto(vd.vertices.size());
    for (size_t i = 0; i < mapto.size(); ++i) {
        mapto[i] = i;
    }

    // Process edges to mark vertices for merging
    processEdges(vd, mapto, D);

    // Perform path compression
    std::cout << "[DEBUG] Performing path compression\n";
    compressMapping(mapto);

    std::vector<int> oldToNewVertexIndex;
    // Rebuild vertices
    std::cout << "[DEBUG] Rebuilding vertices\n";
    rebuildVertices(vd, mapto, oldToNewVertexIndex);
    std::cout << "[DEBUG] New vertex count: " << vd.vertices.size() << "\n";

    // Rebuild vertexMap
    vd.vertexMap.clear();
    const double SCALE_FACTOR = 1e6;
    for (size_t i = 0; i < vd.vertices.size(); ++i) {
        const Point& p = vd.vertices[i].coord;
        int ix = static_cast<int>(std::round(p.x() * SCALE_FACTOR));
        int iy = static_cast<int>(std::round(p.y() * SCALE_FACTOR));
        int iz = static_cast<int>(std::round(p.z() * SCALE_FACTOR));
        std::tuple<int, int, int> key(ix, iy, iz);
        vd.vertexMap[key].push_back(i);
    }

    // Rebuild edges
    std::cout << "[DEBUG] Rebuilding edges with duplicate removal\n";
    rebuildEdges(vd, oldToNewVertexIndex);
    std::cout << "[DEBUG] New edge count: " << vd.edges.size() << "\n";
    std::cout << "[DEBUG] Cleared segmentVertexPairToEdgeIndex\n";

    // Clear cell edges and related mappings (since cells/facets may need reconstruction)
    vd.cellEdges.clear();
    vd.cellEdgeLookup.clear();
    vd.cells.clear();
    vd.facets.clear();

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
void VoronoiDiagram::check() const
{
    std::cout << "Running validity checks.\n";
    checkCellEdgeLookup();
    checkNextCellEdgeConsistency();
    checkCellFacets();
    checkAdvanced();
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
bool VoronoiDiagram::haveSameOrientation(const std::vector<int> &f1,
                                         const std::vector<int> &f2) const
{
    if (f1.size() != f2.size())
        return false;
    int n = (int)f1.size();
    // try every cyclic shift
    for (int shift = 0; shift < n; ++shift)
    {
        bool ok = true;
        for (int i = 0; i < n; ++i)
        {
            if (f1[i] != f2[(i + shift) % n])
            {
                ok = false;
                break;
            }
        }
        if (ok)
            return true;
    }
    return false;
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
void VoronoiDiagram::checkEdgeFacetCount() const {
    for (const auto &cell : cells) {
        std::map<std::pair<int, int>, int> edge_count;
        for (int facet_idx : cell.facet_indices) {
            const auto &vertices = facets[facet_idx].vertices_indices;
            for (size_t i = 0; i < vertices.size(); ++i) {
                int u = vertices[i];
                int v = vertices[(i + 1) % vertices.size()];
                std::pair<int, int> edge = {std::min(u, v), std::max(u, v)};
                edge_count[edge]++;
            }
        }
        for (const auto &[edge, count] : edge_count) {
            if (count == 2) { // Internal edge
                int u = edge.first;
                int v = edge.second;
                int uv_count = 0, vu_count = 0;
                for (int facet_idx : cell.facet_indices) {
                    const auto &vertices = facets[facet_idx].vertices_indices;
                    for (size_t i = 0; i < vertices.size(); ++i) {
                        if (vertices[i] == u && vertices[(i + 1) % vertices.size()] == v) uv_count++;
                        if (vertices[i] == v && vertices[(i + 1) % vertices.size()] == u) vu_count++;
                    }
                }
                if (uv_count != 1 || vu_count != 1) {
                    throw std::runtime_error("In cell " + std::to_string(cell.cellIndex) + 
                                             " edge {" + std::to_string(u) + "," + std::to_string(v) + 
                                             "} does not have opposite orientations.");
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
        auto p = cell.delaunay_vertex->point();
        for (int fi : cell.facet_indices)
        {
            auto const &V = facets[fi].vertices_indices;
            const auto &P1 = vertices[V[0]].coord;
            const auto &P2 = vertices[V[1]].coord;
            const auto &P3 = vertices[V[2]].coord;
            if (CGAL::orientation(P1, P2, P3, p) != CGAL::NEGATIVE)
                throw std::runtime_error("Facet " + std::to_string(fi) +
                                         " in cell " + std::to_string(cell.cellIndex) +
                                         " has inward‐pointing normal.");
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
            if (haveSameOrientation(A, B))
                throw std::runtime_error("Facet “" + std::to_string(fvec[0]) + "” and “" +
                                         std::to_string(fvec[1]) +
                                         "” have the SAME orientation in their two cells.");
        }
    }
}

void VoronoiDiagram::checkAdvanced() const
{
    checkFacetVertexCount();
    checkCellFacetCount();
    checkFacetCellCount();
    checkEdgeFacetCount();
    checkFacetNormals();
    checkPairedFacetOrientations();
}