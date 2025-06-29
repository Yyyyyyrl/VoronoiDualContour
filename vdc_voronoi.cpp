#include "vdc_voronoi.h"

//! @brief Finds the index of the vertex corresponding to the given point in the Voronoi diagram.
/*!
 * Scales the point coordinates and uses a hash map to quickly find the vertex index.
 *
 * @param p The point to find in the Voronoi diagram.
 * @return The index of the vertex if found, otherwise -1.
 */
int VoronoiDiagram::find_vertex(const Point& p) const {
    const double SCALE_FACTOR = 1e6;
    int ix = static_cast<int>(std::round(p.x() * SCALE_FACTOR));
    int iy = static_cast<int>(std::round(p.y() * SCALE_FACTOR));
    int iz = static_cast<int>(std::round(p.z() * SCALE_FACTOR));
    std::tuple<int, int, int> key(ix, iy, iz);
    auto it = vertexMap.find(key);
    if (it != vertexMap.end()) {
        for (int idx : it->second) {
            if (PointApproxEqual()(vertices[idx].coord, p)) return idx;
        }
    }
    return -1;
}

//! @brief Adds a vertex to the Voronoi diagram with the given point and value.
/*!
 * Scales the point coordinates and updates the vertex map and vertex list.
 *
 * @param p The point to add as a vertex.
 * @param value The associated scalar value at the vertex.
 * @return The index of the newly added vertex.
 */
int VoronoiDiagram::AddVertex(const Point& p, float value) {
    int idx = vertices.size();
    VoronoiVertex v(p);
    v.index = idx;
    v.value = value;
    v.cellIndices = {};
    vertices.push_back(v);
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
 * Checks for existing edges and adds a new edge if not found.
 *
 * @param v1 The index of the first vertex of the segment.
 * @param v2 The index of the second vertex of the segment.
 * @param seg The CGAL segment object representing the edge.
 * @return The index of the newly added edge.
 */
int VoronoiDiagram::AddSegmentEdge(int v1, int v2, const Segment3& seg) {
    int minV = std::min(v1, v2);
    int maxV = std::max(v1, v2);
    auto it = segmentVertexPairToEdgeIndex.find({minV, maxV});
    if (it != segmentVertexPairToEdgeIndex.end()) return it->second;
    VoronoiEdge edge(CGAL::make_object(seg));
    edge.vertex1 = v1;
    edge.vertex2 = v2;
    edge.type = 0;
    int edgeIdx = edges.size();
    edges.push_back(edge);
    segmentVertexPairToEdgeIndex[{minV, maxV}] = edgeIdx;
    return edgeIdx;
}

//! @brief Adds a ray edge to the Voronoi diagram.
/*!
 * Adds a new ray edge to the edge list.
 *
 * @param ray The CGAL ray object representing the edge.
 * @return The index of the newly added edge.
 */
int VoronoiDiagram::AddRayEdge(const Ray3& ray) {
    VoronoiEdge edge(CGAL::make_object(ray));
    edge.vertex1 = -1;
    edge.vertex2 = -1;
    edge.type = 1;
    edge.source = ray.source();
    edge.direction = ray.direction().vector();
    int edgeIdx = edges.size();
    edges.push_back(edge);
    return edgeIdx;
}

//! @brief Adds a line edge to the Voronoi diagram.
/*!
 * Adds a new line edge to the edge list.
 *
 * @param line The CGAL line object representing the edge.
 * @return The index of the newly added edge.
 */
int VoronoiDiagram::AddLineEdge(const Line3& line) {
    VoronoiEdge edge(CGAL::make_object(line));
    edge.vertex1 = -1;
    edge.vertex2 = -1;
    edge.type = 2;
    edge.source = line.point(0);
    edge.direction = line.to_vector();
    int edgeIdx = edges.size();
    edges.push_back(edge);
    return edgeIdx;
}

//! @brief Adds a facet to the Voronoi diagram.
/*!
 * Adds a new facet with the given vertex indices.
 *
 * @param vertices_indices The indices of the vertices comprising the facet.
 * @return The index of the newly added facet.
 */
int VoronoiDiagram::AddFacet(const std::vector<int>& vertices_indices) {
    VoronoiFacet facet;
    facet.vertices_indices = vertices_indices;
    facet.mirror_facet_index = -1;
    int facetIdx = facets.size();
    facets.push_back(facet);
    return facetIdx;
}

//! @brief Adds a cell to the Voronoi diagram.
/*!
 * Initializes a new cell with the given Delaunay vertex handle.
 *
 * @param delaunay_vertex The Delaunay vertex handle associated with the cell.
 * @return The index of the newly added cell.
 */
int VoronoiDiagram::AddCell(Vertex_handle delaunay_vertex) {
    VoronoiCell cell(delaunay_vertex);
    cell.cellIndex = cells.size();
    int cellIdx = cells.size();
    cells.push_back(cell);
    return cellIdx;
}
//! @brief Computes the centroid of a cycle using the associated midpoints.
/*!
 * The centroid is calculated using CGAL's `centroid` function, which averages
 * the positions of the midpoints in the cycle. If no midpoints exist, the centroid
 * computation is skipped.
 *
 * @param midpoints Vector of all midpoints in the Voronoi diagram.
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
 * Calculates the geometric center of a cycle by averaging
 * the coordinates of its constituent Voronoi vertices.
 * 
 * @param voronoiVertices List of all Voronoi vertices in the diagram
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
 * Iterates through all edges in the Voronoi diagram and marks vertices for merging
 * if their distance is less than or equal to the threshold D.
 *
 * @param vd The Voronoi diagram containing edges to process
 * @param mapto Union-find structure tracking vertex merges
 * @param D Distance threshold for merging vertices
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
 * Compresses paths in the union-find structure to make future lookups faster.
 * Updates each entry to point directly to its root ancestor.
 *
 * @param mapto Union-find structure to compress
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
 * Creates new vertex list by keeping only root vertices from the union-find structure.
 * Also builds mapping from old vertex indices to new ones.
 *
 * @param vd Voronoi diagram to modify
 * @param mapto Union-find structure tracking vertex merges
 * @param oldToNewVertexIndex Output mapping from old to new vertex indices
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
 * Processes edges to remove duplicates and update vertex references.
 * Handles segment edges, ray edges and line edges separately.
 *
 * @param vd Voronoi diagram to modify
 * @param oldToNewVertexIndex Mapping from old to new vertex indices
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
 * Checks that all nextCellEdge pointers in VoronoiCellEdge objects:
 * - Point to valid indices within cellEdges vector
 * - Reference edges with matching edgeIndex values
 *
 * @throws std::runtime_error if invalid pointers are found
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

// Helper function to check edge cycles
//! @brief Verifies edge cycles form complete rings.
/*!
 * Checks that edges sharing the same edgeIndex form complete cycles
 * through nextCellEdge pointers without breaks or sub-loops.
 *
 * @throws std::runtime_error if cycles are incomplete or malformed
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

            const VoronoiFacet &facet = facets[fIdx];
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
 * Creates a unique identifier for a facet by sorting and selecting
 * the three smallest vertex indices.
 *
 * @param verts Vector of vertex indices comprising the facet
 * @return Tuple containing the three smallest vertex indices
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
 * Compares two facets' vertex orderings by checking all cyclic permutations.
 *
 * @param f1 First facet's vertex indices
 * @param f2 Second facet's vertex indices
 * @return True if facets have identical vertex ordering, false otherwise
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

// Helper function to check facet vertex count
//! @brief Verifies all facets have at least 3 vertices.
/*!
 * Checks that no facet in the diagram has fewer than 3 vertices,
 * which would make it geometrically invalid.
 *
 * @throws std::runtime_error if any facet has <3 vertices
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

// Helper function to check cell facet count
//! @brief Verifies all cells have at least 4 facets.
/*!
 * Checks that no cell in the diagram has fewer than 4 facets,
 * which is the minimum for a 3D polyhedron.
 *
 * @throws std::runtime_error if any cell has <4 facets
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

// Helper function to check facet cell count
//! @brief Verifies facet sharing between cells.
/*!
 * Ensures that each facet is shared by at most 2 cells,
 * as expected in a proper cell complex.
 *
 * @throws std::runtime_error if any facet is shared by >2 cells
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

// Helper function to check edge facet count
//! @brief Verifies edge-facet relationships.
/*!
 * Checks that edges appear exactly twice in cell facets,
 * once in each direction, forming proper boundaries.
 *
 * @throws std::runtime_error if edge counts are incorrect
 */
void VoronoiDiagram::checkEdgeFacetCount() const
{
    for (const auto &cell : cells)
    {
        std::map<std::pair<int, int>, int> orientedCount;
        for (int fi : cell.facet_indices)
        {
            auto const &V = facets[fi].vertices_indices;
            int M = (int)V.size();
            for (int k = 0; k < M; ++k)
            {
                int u = V[k], v = V[(k + 1) % M];
                orientedCount[{u, v}]++;
            }
        }
        std::set<std::pair<int, int>> seen;
        for (auto const &kv : orientedCount)
        {
            auto uv = kv.first;
            auto cnt_uv = kv.second;
            auto vu = std::make_pair(uv.second, uv.first);
            if (seen.count(uv) || seen.count(vu))
                continue;
            int cnt_vu = orientedCount.count(vu) ? orientedCount.at(vu) : 0;
            if (cnt_uv + cnt_vu != 2 || cnt_uv != 1 || cnt_vu != 1)
                throw std::runtime_error("In cell " + std::to_string(cell.cellIndex) +
                                         " edge {" + std::to_string(uv.first) + "," +
                                         std::to_string(uv.second) +
                                         "} does not appear exactly twice with opposite orientations.");
            seen.insert(uv);
            seen.insert(vu);
        }
    }
}

// Helper function to check facet normals
//! @brief Verifies facet normals point outward.
/*!
 * Ensures all facet normals in cells point away from
 * the cell's Delaunay vertex, as required for proper orientation.
 *
 * @throws std::runtime_error if any normal points inward
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

// Helper function to check paired facet orientations
//! @brief Verifies paired facets have opposite orientations.
/*!
 * Checks that facets shared between two cells have
 * opposite vertex orderings (clockwise vs counter-clockwise).
 *
 * @throws std::runtime_error if orientations match
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