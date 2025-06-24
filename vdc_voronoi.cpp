#include "vdc_voronoi.h"

int VoronoiDiagram::find_vertex(const Point &p) const
{
    const double SCALE_FACTOR = 1e6;
    int ix = static_cast<int>(std::round(p.x() * SCALE_FACTOR));
    int iy = static_cast<int>(std::round(p.y() * SCALE_FACTOR));
    int iz = static_cast<int>(std::round(p.z() * SCALE_FACTOR));
    std::tuple<int, int, int> key(ix, iy, iz);
    auto it = vertexMap.find(key);
    if (it != vertexMap.end())
    {
        for (int idx : it->second)
        {
            if (PointApproxEqual()(vertices[idx].vertex, p))
            {
                return idx;
            }
        }
    }
    return -1;
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


// Compute centroid of the cycle using the positions in voronoiVertices.
void Cycle::compute_centroid(const std::vector<VoronoiVertex> &voronoiVertices)
{
    double sumX = 0, sumY = 0, sumZ = 0;
    for (int idx : midpoint_indices)
    {
        const Point &p = voronoiVertices[idx].vertex;
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

 // Helper function to process edges and mark vertices for merging
static void processEdges(const VoronoiDiagram& vd, std::vector<int>& mapto, double D) {
    for (size_t edgeIdx = 0; edgeIdx < vd.edges.size(); ++edgeIdx) {
        // Get vertex indices for the current edge
        int v1 = vd.edgeVertexIndices[edgeIdx].first;
        int v2 = vd.edgeVertexIndices[edgeIdx].second;

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
        Point p1 = vd.vertices[v1].vertex;
        Point p2 = vd.vertices[v2].vertex;
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

// Helper function to merge vertices within a small tolerance
static void mergeCloseVertices(const VoronoiDiagram& vd, double mergeTolerance, std::vector<int>& mapto) {
    for (size_t i = 0; i < vd.vertices.size(); ++i) {
        for (size_t j = i + 1; j < vd.vertices.size(); ++j) {
            double dist = CGAL::sqrt(CGAL::squared_distance(vd.vertices[i].vertex, vd.vertices[j].vertex));
            if (dist <= mergeTolerance) {
                int v1 = i, v2 = j;
                while (mapto[v1] != v1) v1 = mapto[v1];
                while (mapto[v2] != v2) v2 = mapto[v2];
                if (v1 < v2)
                    mapto[v2] = v1;
                else
                    mapto[v1] = v2;
            }
        }
    }
}

// Helper function to perform path compression
static void compressMapping(std::vector<int>& mapto) {
    for (size_t i = 0; i < mapto.size(); ++i) {
        int root = i;
        while (mapto[root] != root) root = mapto[root];
        mapto[i] = root;
    }
}

// Helper function to rebuild vertices
static void rebuildVertices(VoronoiDiagram& vd, const std::vector<int>& mapto) {
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

    vd.oldToNewVertexIndex.resize(mapto.size());
    for (size_t i = 0; i < mapto.size(); ++i) {
        vd.oldToNewVertexIndex[i] = oldToNewIndex[mapto[i]];
    }
}

// Helper function to rebuild edges with duplicate removal
static void rebuildEdges(VoronoiDiagram& vd) {
    std::map<std::pair<int, int>, int> segmentMap;
    std::set<CGAL::Object, ObjectComparator> rayLineSet;

    // Process all edges to identify unique ones
    for (size_t edgeIdx = 0; edgeIdx < vd.edges.size(); ++edgeIdx) {
        const auto& edge = vd.edges[edgeIdx];
        Segment3 seg;
        Ray3 ray;
        Line3 line;

        if (CGAL::assign(seg, edge)) {
            // Handle segments using vertex indices
            int oldIdx1 = vd.edgeVertexIndices[edgeIdx].first;
            int oldIdx2 = vd.edgeVertexIndices[edgeIdx].second;
            if (oldIdx1 >= 0 && oldIdx2 >= 0) {
                int newIdx1 = vd.oldToNewVertexIndex[oldIdx1];
                int newIdx2 = vd.oldToNewVertexIndex[oldIdx2];
                if (newIdx1 != newIdx2) { // Skip collapsed edges (same vertex)
                    int v1 = std::min(newIdx1, newIdx2);
                    int v2 = std::max(newIdx1, newIdx2);
                    if (segmentMap.find({v1, v2}) == segmentMap.end()) {
                        segmentMap[{v1, v2}] = edgeIdx;
                    }
                }
            }
        } else if (CGAL::assign(ray, edge) || CGAL::assign(line, edge)) {
            // Handle rays and lines with geometric comparison
            rayLineSet.insert(edge);
        }
    }

    // Rebuild the edges and edgeVertexIndices vectors
    std::vector<Object> newEdges;
    std::vector<std::pair<int, int>> newEdgeVertexIndices;
    newEdges.reserve(segmentMap.size() + rayLineSet.size());
    newEdgeVertexIndices.reserve(segmentMap.size());

    for (const auto& kv : segmentMap) {
        size_t edgeIdx = kv.second;
        newEdges.push_back(vd.edges[edgeIdx]);
        int oldIdx1 = vd.edgeVertexIndices[edgeIdx].first;
        int oldIdx2 = vd.edgeVertexIndices[edgeIdx].second;
        int newIdx1 = vd.oldToNewVertexIndex[oldIdx1];
        int newIdx2 = vd.oldToNewVertexIndex[oldIdx2];
        newEdgeVertexIndices.emplace_back(newIdx1, newIdx2);
    }
    for (const auto& edge : rayLineSet) {
        newEdges.push_back(edge);
        newEdgeVertexIndices.emplace_back(-1, -1); // Rays and lines have no vertex indices
    }

    vd.edges = std::move(newEdges);
    vd.edgeVertexIndices = std::move(newEdgeVertexIndices);
    vd.segmentVertexPairToEdgeIndex.clear();
}

// Standalone function to collapse small edges
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

    // Rebuild vertices
    std::cout << "[DEBUG] Rebuilding vertices\n";
    rebuildVertices(vd, mapto);
    std::cout << "[DEBUG] New vertex count: " << vd.vertices.size() << "\n";

    // Rebuild vertexMap
    vd.vertexMap.clear();
    const double SCALE_FACTOR = 1e6;
    for (size_t i = 0; i < vd.vertices.size(); ++i) {
        const Point& p = vd.vertices[i].vertex;
        int ix = static_cast<int>(std::round(p.x() * SCALE_FACTOR));
        int iy = static_cast<int>(std::round(p.y() * SCALE_FACTOR));
        int iz = static_cast<int>(std::round(p.z() * SCALE_FACTOR));
        std::tuple<int, int, int> key(ix, iy, iz);
        vd.vertexMap[key].push_back(i);
    }

    // Rebuild edges
    std::cout << "[DEBUG] Rebuilding edges with duplicate removal\n";
    rebuildEdges(vd);
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

std::tuple<int, int, int> VoronoiDiagram::getFacetHashKey(const std::vector<int> &verts) const
{
    // pick three smallest indices
    std::array<int, 3> a;
    std::vector<int> tmp = verts;
    std::sort(tmp.begin(), tmp.end());
    a = {tmp[0], tmp[1], tmp[2]};
    return std::make_tuple(a[0], a[1], a[2]);
}

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
void VoronoiDiagram::checkFacetNormals() const
{
    for (const auto &cell : cells)
    {
        auto p = cell.delaunay_vertex->point();
        for (int fi : cell.facet_indices)
        {
            auto const &V = facets[fi].vertices_indices;
            const auto &P1 = vertices[V[0]].vertex;
            const auto &P2 = vertices[V[1]].vertex;
            const auto &P3 = vertices[V[2]].vertex;
            if (CGAL::orientation(P1, P2, P3, p) != CGAL::NEGATIVE)
                throw std::runtime_error("Facet " + std::to_string(fi) +
                                         " in cell " + std::to_string(cell.cellIndex) +
                                         " has inward‐pointing normal.");
        }
    }
}

// Helper function to check paired facet orientations
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