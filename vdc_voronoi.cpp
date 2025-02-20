#include "vdc_voronoi.h"

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

//! @brief Checks internal consistency of the VoronoiDiagram.
void VoronoiDiagram::check() const
{
    checkCellEdgeLookup();
    checkNextCellEdgeConsistency();
    checkCellFacets();
    std::cout << "VoronoiDiagram::check() passed all checks.\n";
}

//! @brief Verifies that `cellEdgeLookup` matches the data in `VoronoiCellEdges`.
void VoronoiDiagram::checkCellEdgeLookup() const
{
    for (const auto &kv : cellEdgeLookup)
    {
        // kv.first is (ic, ie) and kv.second is the index in VoronoiCellEdges.
        int ic = kv.first.first;  // cellIndex
        int ie = kv.first.second; // edgeIndex
        int cellEdgeIdx = kv.second;

        if (cellEdgeIdx < 0 || cellEdgeIdx >= static_cast<int>(VoronoiCellEdges.size()))
        {
            throw std::runtime_error("cellEdgeLookup points to invalid VoronoiCellEdge index.");
        }

        const VoronoiCellEdge &ce = VoronoiCellEdges[cellEdgeIdx];

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

//! @brief Checks that each cellEdge's `nextCellEdge` points to another edge with the same `edgeIndex`.
void VoronoiDiagram::checkNextCellEdgeConsistency() const
{
    // 1. Check that each nextCellEdge is valid and has the same edgeIndex.
    for (int ceIdx = 0; ceIdx < static_cast<int>(VoronoiCellEdges.size()); ++ceIdx)
    {
        const VoronoiCellEdge &ce = VoronoiCellEdges[ceIdx];
        int nxt = ce.nextCellEdge;
        if (nxt < 0)
        {
            // -1 might be allowed for boundary conditions or partial references.
            continue;
        }
        if (nxt >= static_cast<int>(VoronoiCellEdges.size()))
        {
            std::cerr << "ERROR: VoronoiCellEdge[" << ceIdx << "].nextCellEdge=" << nxt
                      << " is out of range.\n";
            throw std::runtime_error("Invalid nextCellEdge index.");
        }
        const VoronoiCellEdge &ceNext = VoronoiCellEdges[nxt];
        if (ceNext.edgeIndex != ce.edgeIndex)
        {
            std::cerr << "ERROR: VoronoiCellEdge[" << ceIdx << "] -> edgeIndex="
                      << ce.edgeIndex << " but nextCellEdge=" << nxt
                      << " has edgeIndex=" << ceNext.edgeIndex << "\n";
            throw std::runtime_error("Inconsistent nextCellEdge edgeIndex.");
        }
    }

    // 2. Extended check: each set of cell edges with the same edgeIndex forms a single closed cycle.
    std::unordered_map<int, std::vector<int>> edgeIndexToCellEdges;
    for (int ceIdx = 0; ceIdx < static_cast<int>(VoronoiCellEdges.size()); ++ceIdx)
    {
        int eIdx = VoronoiCellEdges[ceIdx].edgeIndex;
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

        int current = VoronoiCellEdges[start].nextCellEdge;
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
            const VoronoiCellEdge &ceNext = VoronoiCellEdges[current];
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

//! @brief Checks each VoronoiCell's facets to ensure that every facet's vertices are in the cell's vertex set.
void VoronoiDiagram::checkCellFacets() const
{
    for (int cIdx = 0; cIdx < static_cast<int>(voronoiCells.size()); ++cIdx)
    {
        const VoronoiCell &cell = voronoiCells[cIdx];
        // Build a set of the cell's vertex indices for quick membership testing.
        std::set<int> cellVertexSet(cell.vertices_indices.begin(), cell.vertices_indices.end());

        for (int fIdx : cell.facet_indices)
        {
            if (fIdx < 0 || fIdx >= static_cast<int>(voronoiFacets.size()))
            {
                std::cerr << "ERROR: cell " << cIdx << " has invalid facet index " << fIdx << "\n";
                throw std::runtime_error("Facet index out of range.");
            }

            const VoronoiFacet &facet = voronoiFacets[fIdx];
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