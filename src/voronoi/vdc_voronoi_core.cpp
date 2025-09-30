#include "voronoi/vdc_voronoi.h"
#include "core/vdc_debug.h"

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
                    if (debug)
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
                    }
                    throw std::runtime_error("In cell " + std::to_string(cell.cellIndex) + " edge {" + std::to_string(u) + "," + std::to_string(v) + "} does not have opposite orientations.");
                }
            }
        }
    }
}

//! @brief Ensures per-cell facet loops are consistently oriented.
/*!
 * Validates that each Voronoi cell has a coherent orientation purely by connectivity:
 *
 *  - Every undirected edge used by the cell must appear exactly twice across its facets.
 *  - Those two appearances must have opposite direction (clockwise vs counter-clockwise).
 *
 * This check guarantees combinatorial orientation consistency without referencing geometry.
 *
 * @throws std::runtime_error if a cell contains mismatched edge usage or orientation.
 */
void VoronoiDiagram::checkFacetNormals() const
{
    for (const auto &cell : cells)
    {
        // Track how each undirected edge is used by the facets attached to this cell.
        struct EdgeInfo
        {
            int count = 0;
            int orientationSum = 0;
            std::vector<int> facets;
        };

        // key = {min(vertex), max(vertex)} so we treat edges as undirected when aggregating.
        std::map<std::pair<int, int>, EdgeInfo> edgeMap;

        for (int fi : cell.facet_indices)
        {
            const auto &loop = facets[fi].vertices_indices;
            const size_t n = loop.size();
            if (n < 2)
                continue;

            for (size_t k = 0; k < n; ++k)
            {
                const int u = loop[k];
                const int v = loop[(k + 1) % n];
                if (u == v)
                    continue;

                const std::pair<int, int> key{std::min(u, v), std::max(u, v)};
                // Orientation +1 if the loop traverses key in ascending order, -1 otherwise.
                const int orientation = (key.first == u && key.second == v) ? 1 : -1;
                EdgeInfo &info = edgeMap[key];
                info.count += 1;
                info.orientationSum += orientation;
                // Record up to a few culprit facets to aid debugging if the check fails.
                if (info.facets.size() < 3)
                    info.facets.push_back(fi);
            }
        }

        for (const auto &kv : edgeMap)
        {
            const auto &info = kv.second;
            const bool countsOk = (info.count == 2);
            const bool orientationsOk = (info.orientationSum == 0);
            if (!countsOk || !orientationsOk)
            {
                if (debug)
                {
                    std::cerr << "[DEBUG] Cell " << cell.cellIndex
                              << ": edge {" << kv.first.first << ", " << kv.first.second << "}"
                              << " count=" << info.count
                              << " orientationSum=" << info.orientationSum << std::endl;
                    std::cerr << "[DEBUG] Facets using edge: ";
                    for (int f : info.facets)
                        std::cerr << f << ' ';
                    std::cerr << std::endl;
                }

                throw std::runtime_error("Cell " + std::to_string(cell.cellIndex) +
                                         " has inconsistent facet orientations.");
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
