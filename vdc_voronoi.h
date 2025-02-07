//! @file vdc_voronoi.h
//! @brief Header file defining structures and functions for Voronoi Diagram construction, cycles, facets, and isosurfaces.

#ifndef VDC_VORONOI_H
#define VDC_VORONOI_H

#include "vdc_type.h"
#include "vdc_delaunay.h"

//! @brief Represents a midpoint of an edge in a Voronoi diagram.
/*!
 * A midpoint lies on the edge of a Voronoi facet. This structure holds the
 * geometric coordinates of the midpoint, its connectivity, and metadata about
 * its relationship to facets and cycles.
 */
struct MidpointNode
{
    Point point;                   //!< Geometric coordinates of the midpoint.
    std::vector<int> connected_to; //!< Indices of midpoints connected to this one, forming graph edges.
    int facet_index;               //!< Index of the facet this midpoint lies on.
    int cycle_index;               //!< Index of the cycle this midpoint belongs to.
    int global_edge_index;
};

//! @brief Represents a facet in a Voronoi diagram.
/*!
 * A Voronoi facet is defined by a list of vertices and their scalar values.
 * Facets are part of Voronoi cells.
 */
struct VoronoiFacet
{
    std::vector<int> vertices_indices; //!< Indices of vertices forming this facet, ordered.
    std::vector<float> vertex_values;  //!< Scalar values at the facet's vertices.

    //! @brief Default constructor.
    VoronoiFacet() = default;
};

//! @brief Represents a closed cycle in a Voronoi cell formed by midpoints.
/*!
 * A cycle is a loop of midpoints connected by edges. Each cycle is associated
 * with a Voronoi cell and may have its centroid computed as the cycle's isovertex.
 */
struct Cycle
{
    std::vector<int> midpoint_indices;      //!< Indices of midpoints forming this cycle.
    std::vector<std::pair<int, int>> edges; //!< Edges forming the cycle, represented by pairs of midpoint indices.
    Point isovertex;                        //!< Geometric centroid of the cycle, representing the isovertex.
    int voronoi_cell_index;                 //!< Index of the Voronoi cell this cycle belongs to.

    // TODO: record it using the indices of the voronoi edges instead of the midpoints
    //! @brief Computes the centroid of the cycle.
    /*!
     * The centroid (isovertex) is calculated using the CGAL `centroid` function
     * based on the geometric positions of the midpoints.
     *
     * @param midpoints Vector of all midpoints in the diagram. Only the midpoints
     *                  associated with this cycle are used in the computation.
     */
    void compute_centroid(const std::vector<MidpointNode> &midpoints);
};

//! @brief Represents a vertex in a Voronoi diagram.
/*!
 * A Voronoi vertex is a point where multiple Voronoi edges meet. Each vertex
 * may belong to multiple Voronoi cells.
 */
struct VoronoiVertex
{
    Point vertex;                 //!< Geometric coordinates of the vertex.
    std::vector<int> cellIndices; //!< Indices of Voronoi cells that contain this vertex.

    //! @brief Constructor to initialize a Voronoi vertex.
    /*!
     * @param p Coordinates of the vertex.
     */
    VoronoiVertex(Point p) : vertex(p) {}
};

//! @brief Represents a Voronoi cell (polytope) in the Voronoi diagram.
/*!
 * A Voronoi cell is a polyhedral region associated with a single Delaunay vertex.
 * It contains facets, vertices, and cycles, and may also include isosurface data.
 */
struct VoronoiCell
{
    Vertex_handle delaunay_vertex;     //!< Handle to the corresponding Delaunay vertex.
    int cellIndex;                     //!< Index of this cell in the Voronoi diagram.
    std::vector<int> vertices_indices; //!< Indices of Voronoi vertices forming this cell.
    std::vector<int> facet_indices;    //!< Indices of Voronoi facets belonging to this cell.
    CGAL::Polyhedron_3<K> polyhedron;  //!< Geometric representation of the cell as a polyhedron.
    std::vector<Cycle> cycles;         //!< Cycles (loops) within this cell.
    int isoVertexStartIndex;           //!< Starting index of isosurface vertices associated with this cell.
    int numIsoVertices;                //!< Number of isosurface vertices in this cell.

    //! @brief Constructor to initialize a Voronoi cell.
    /*!
     * @param vh Handle to the corresponding Delaunay vertex.
     */
    VoronoiCell(Vertex_handle vh)
        : delaunay_vertex(vh), isoVertexStartIndex(-1), numIsoVertices(0) {}
};

struct VoronoiCellEdge
{
    int cellIndex;                 //!< Index of the VoronoiCell that contains this CellEdge
    int edgeIndex;                 //!< Index of this edge in the vector of voronoiEdges in the VoronoiDiagram instance that contains the cell this egde is in
    std::vector<int> cycleIndices; //!< Indices of cycles in this VoronoiCell that corresponding to this edge
    int nextCellEdge;              //!< Index of next cell edge around the Voronoi Edge ( VoronoiDiagram.voronoiEdges[edgeIndex])
};

//! @brief Represents the overall Voronoi diagram.
/*!
 * The Voronoi diagram consists of vertices, edges, cells, facets, and isosurface data.
 */
struct VoronoiDiagram
{
    std::vector<VoronoiVertex> voronoiVertices;                       //!< List of Voronoi vertices in the diagram.
    std::vector<Object> voronoiEdges;                                 //!< List of edges in the diagram (e.g., line segments).
    std::vector<VoronoiCellEdge> VoronoiCellEdges;                    //!< List of Cell Edges in the diagram
    std::vector<float> voronoiVertexValues;                           //!< Scalar values at the Voronoi vertices.
    std::vector<VoronoiCell> voronoiCells;                            //!< List of Voronoi cells in the diagram.
    std::vector<VoronoiFacet> voronoiFacets;                          //!< List of facets in the diagram.
    std::map<Cell_handle, int> delaunay_cell_to_voronoi_vertex_index;                  //!< Map from Voronoi cells to vertex indices.
    std::map<Point, int> point_to_vertex_index;                       //!< Map from Voronoi vertices to their indices.
    std::map<Vertex_handle, int> delaunayVertex_to_voronoiCell_index; //!< Map from Delaunay vertices to Voronoi cells.

    std::map<std::pair<int, int>, int> cellEdgeLookup;               //!< Maps (cellIndex, edgeIndex) -> index in VoronoiCellEdges
    std::map<std::pair<int, int>, int> segmentVertexPairToEdgeIndex; //!< a map from a pair of Voronoi vertex indices (v_1, v_2) (in ascending order) to the edgeIndex in voronoiDiagram

    //! @brief Checks internal consistency of the VoronoiDiagram.
    void check() const
    {
        checkCellEdgeLookup();
        checkNextCellEdgeConsistency();
        checkCellFacets();
        std::cout << "VoronoiDiagram::check() passed all checks.\n";
    }

private:
    //! @brief Verifies that `cellEdgeLookup` matches the data in `VoronoiCellEdges`.
    void checkCellEdgeLookup() const
    {
        for (const auto &kv : cellEdgeLookup)
        {
            // kv.first is (ic, ie)
            // kv.second is the index in VoronoiCellEdges
            int ic = kv.first.first;  // cellIndex
            int ie = kv.first.second; // edgeIndex
            int cellEdgeIdx = kv.second;

            if (cellEdgeIdx < 0 || cellEdgeIdx >= VoronoiCellEdges.size())
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
    void checkNextCellEdgeConsistency() const
{
    // 1. Original checks: each nextCellEdge is valid and has the same edgeIndex
    for (int ceIdx = 0; ceIdx < (int)VoronoiCellEdges.size(); ++ceIdx)
    {
        const VoronoiCellEdge &ce = VoronoiCellEdges[ceIdx];
        int nxt = ce.nextCellEdge;
        if (nxt < 0)
        {
            // -1 might be allowed for boundary conditions or partial references
            // If your logic demands a strict ring (always nextCellEdge != -1),
            // then throw or log an error here instead.
            continue;
        }
        if (nxt >= (int)VoronoiCellEdges.size())
        {
            std::cerr << "ERROR: VoronoiCellEdge[" << ceIdx << "].nextCellEdge=" << nxt
                      << " is out of range.\n";
            throw std::runtime_error("Invalid nextCellEdge index.");
        }
        // Check that next edge has the same edgeIndex
        const VoronoiCellEdge &ceNext = VoronoiCellEdges[nxt];
        if (ceNext.edgeIndex != ce.edgeIndex)
        {
            std::cerr << "ERROR: VoronoiCellEdge[" << ceIdx << "] -> edgeIndex="
                      << ce.edgeIndex << " but nextCellEdge=" << nxt
                      << " has edgeIndex=" << ceNext.edgeIndex << "\n";
            throw std::runtime_error("Inconsistent nextCellEdge edgeIndex.");
        }
    }

    // 2. Extended check: each set of cell edges with the same edgeIndex forms a single closed cycle
    //    (assuming your data structure is meant to form exactly one ring per edgeIndex).
    std::unordered_map<int, std::vector<int>> edgeIndexToCellEdges;
    for (int ceIdx = 0; ceIdx < (int)VoronoiCellEdges.size(); ++ceIdx)
    {
        int eIdx = VoronoiCellEdges[ceIdx].edgeIndex;
        edgeIndexToCellEdges[eIdx].push_back(ceIdx);
    }

    for (const auto &kv : edgeIndexToCellEdges)
    {
        const int eIdx = kv.first;
        const auto &ceIndices = kv.second;
        if (ceIndices.empty()) {
            continue; // no edges? skip
        }

        // Start from the first cell edge in this group
        int start = ceIndices[0];

        // Follow nextCellEdge pointers until we loop back
        std::set<int> visited;
        visited.insert(start);

        int current = VoronoiCellEdges[start].nextCellEdge;
        // If we allow -1 nextCellEdge, we might skip the ring-check
        // But let's assume every edge in a ring has nextCellEdge >= 0
        while (current != start)
        {
            if (current < 0)
            {
                // Means we reached a cellEdge with nextCellEdge = -1 => open chain
                std::cerr << "ERROR: The edges for edgeIndex=" << eIdx
                          << " do not form a complete cycle (nextCellEdge=-1 encountered).\n";
                throw std::runtime_error("Incomplete ring around an edge.");
            }

            if (visited.find(current) != visited.end())
            {
                // Means we looped early -> smaller cycle, or we have multiple cycles
                std::cerr << "ERROR: The edges for edgeIndex=" << eIdx
                          << " contain a sub-loop. Edge " << current
                          << " was already visited.\n";
                throw std::runtime_error("Multiple loops or early cycle detected.");
            }

            visited.insert(current);

            const VoronoiCellEdge &ceNext = VoronoiCellEdges[current];
            // Continue walking
            current = ceNext.nextCellEdge;
        }

        // Now we've come back to 'start'
        // Check if we've visited exactly all edges in ceIndices
        if (visited.size() != ceIndices.size())
        {
            std::cerr << "ERROR: For edgeIndex=" << eIdx
                      << ", visited " << visited.size()
                      << " edges, but we expected " << ceIndices.size() << ".\n"
                      << "Implying there's a second disconnected cycle or missing edges.\n";
            throw std::runtime_error("Ring does not include all edges for edgeIndex.");
        }

        // If get here, we have exactly one ring containing all edges in ceIndices
    }
}

    //! @brief Checks each VoronoiCell's facets to ensure that every facet's vertices are in the cell's vertex set.
    void checkCellFacets() const
    {
        for (int cIdx = 0; cIdx < (int)voronoiCells.size(); ++cIdx)
        {
            const VoronoiCell &cell = voronoiCells[cIdx];
            // Build a set of the cell's vertex indices for quick membership testing
            std::set<int> cellVertexSet(cell.vertices_indices.begin(), cell.vertices_indices.end());

            for (int fIdx : cell.facet_indices)
            {
                if (fIdx < 0 || fIdx >= (int)voronoiFacets.size())
                {
                    std::cerr << "ERROR: cell " << cIdx << " has invalid facet index " << fIdx << "\n";
                    throw std::runtime_error("Facet index out of range.");
                }

                const VoronoiFacet &facet = voronoiFacets[fIdx];
                for (int vIdx : facet.vertices_indices)
                {
                    // Check if vIdx is in the cell's vertex set
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
};

//! @brief Represents an isosurface in the domain.
/*!
 * An isosurface is a 3D surface represented by vertices and triangles,
 * extracted based on scalar values.
 */
struct IsoSurface
{
    std::vector<Point> isosurfaceVertices;                           //!< Vertices of the isosurface.
    std::vector<std::tuple<int, int, int>> isosurfaceTrianglesMulti; //!< Triangles forming the isosurface.
    std::vector<DelaunayTriangle> isosurfaceTrianglesSingle;
};

//! @brief Represents a midpoint on an edge, along with its facet information.
struct EdgeMidpoint
{
    Point midpoint;  //!< Geometric coordinates of the midpoint.
    int facet_index; //!< Index of the facet this midpoint belongs to.
};

//! @brief Represents a vertex in the Delaunay triangulation.
/*!
 * A Delaunay vertex can either be part of the actual triangulation or a dummy
 * point added to bound the Voronoi diagram.
 */
struct DelaunayVertex
{
    Point point;   //!< Coordinates of the vertex.
    bool is_dummy; //!< Flag indicating whether the vertex is a dummy point.
};

//! @brief Comparator for points with approximate equality.
/*!
 * This comparator allows points to be compared with a small tolerance
 * (epsilon) to account for floating-point inaccuracies.
 */
struct PointComparator
{
    //! @brief Compares two points for approximate equality.
    /*!
     * @param a First point.
     * @param b Second point.
     * @return `true` if the points are approximately equal, `false` otherwise.
     */
    bool operator()(const Point &a, const Point &b) const
    {
        const double epsilon = 1e-6;
        return (std::fabs(a.x() - b.x()) < epsilon) &&
               (std::fabs(a.y() - b.y()) < epsilon) &&
               (std::fabs(a.z() - b.z()) < epsilon);
    }
};

//! @brief Overloaded output operator for VoronoiFacet.
template <typename OSTREAM_TYPE>
OSTREAM_TYPE &operator<<(OSTREAM_TYPE &os, const VoronoiFacet &vf)
{
    os << "VoronoiFacet:\n";
    os << "  Vertices indices: ";
    for (const int idx : vf.vertices_indices)
        os << idx << " ";
    os << "\n";

    os << "  Vertex values: ";
    for (const float val : vf.vertex_values)
        os << val << " ";
    os << "\n";

    return os;
}

//! @brief Overloaded output operator for Cycle
template <typename OSTREAM_TYPE>
OSTREAM_TYPE &operator<<(OSTREAM_TYPE &os, const Cycle &cycle)
{
    os << "Cycle:\n";
    os << "  Voronoi cell index: " << cycle.voronoi_cell_index << "\n";
    os << "  Isovertex: " << cycle.isovertex << "\n";

    os << "  Midpoint indices: ";
    for (const int idx : cycle.midpoint_indices)
        os << idx << " ";
    os << "\n";

    os << "  Edges: ";
    for (const auto &edge : cycle.edges)
        os << "(" << edge.first << ", " << edge.second << ") ";
    os << "\n";

    return os;
}

//! @brief Overloaded output operator for VoronoiVertex
template <typename OSTREAM_TYPE>
OSTREAM_TYPE &operator<<(OSTREAM_TYPE &os, const VoronoiVertex &vv)
{
    os << "VoronoiVertex:\n";
    os << "  Vertex: " << vv.vertex << "\n";

    return os;
}

//! @brief Overloaded output operator for VoronoiCell
template <typename OSTREAM_TYPE>
OSTREAM_TYPE &operator<<(OSTREAM_TYPE &os, const VoronoiCell &vc)
{
    os << "VoronoiCell:\n";
    os << "  Cell index: " << vc.cellIndex << "\n";
    os << "  Delaunay vertex: " << vc.delaunay_vertex->point() << "\n";

    os << "  Voronoi Vertices indices: ";
    for (const int idx : vc.vertices_indices)
        os << idx << " ";
    os << "\n";

    os << "  Facet indices: ";
    for (const int idx : vc.facet_indices)
        os << idx << " ";
    os << "\n";

    os << "  IsoVertex start index: " << vc.isoVertexStartIndex << "\n";
    os << "  Number of isoVertices: " << vc.numIsoVertices << "\n";

    os << "  Cycles:\n";
    for (const auto &cycle : vc.cycles)
        os << cycle;

    return os;
}

//! @brief Overloaded output operator for VoronoiDiagram
template <typename OSTREAM_TYPE>
OSTREAM_TYPE &operator<<(OSTREAM_TYPE &os, const VoronoiDiagram &vd)
{
    os << "VoronoiDiagram:\n";

    // 1. Voronoi Vertices
    os << "\nVoronoiVertices:\n";
    for (size_t i = 0; i < vd.voronoiVertices.size(); ++i)
    {
        os << "Index " << i << ":\n";
        os << vd.voronoiVertices[i];
    }

    // 2. Voronoi Edges
    os << "\nVoronoiEdges:\n";
    for (const auto &edge : vd.voronoiEdges)
    {
        os << "  Edge: ";
        Segment3 segment;
        Line3 line;
        Ray3 ray;

        if (CGAL::assign(segment, edge))
        {
            os << "Segment(" << segment.source() << " - " << segment.target() << ")\n";
        }
        else if (CGAL::assign(line, edge))
        {
            os << "Line(" << line.point(0) << " - " << line.point(1) << ")\n";
        }
        else if (CGAL::assign(ray, edge))
        {
            os << "Ray(" << ray.source() << ", direction: " << ray.direction() << ")\n";
        }
        else
        {
            os << "Unknown edge type.\n";
        }
    }

    // 3. Voronoi Vertex Values
    os << "\nVoronoiVertexValues:\n";
    for (size_t i = 0; i < vd.voronoiVertexValues.size(); ++i)
    {
        os << "  Index " << i << ": " << vd.voronoiVertexValues[i] << "\n";
    }

    // 4. Voronoi Facets
    os << "\nVoronoiFacets:\n";
    for (size_t i = 0; i < vd.voronoiFacets.size(); ++i)
    {
        os << "Index " << i << ":\n";
        os << vd.voronoiFacets[i];
    }

    // 5. Voronoi Cells
    os << "\nVoronoiCells:\n";
    for (const auto &cell : vd.voronoiCells)
    {
        os << "\n" << cell;
    }

    // 6. Voronoi CellEdges
    os << "\nVoronoiCellEdges:\n";
    for (size_t i = 0; i < vd.VoronoiCellEdges.size(); ++i)
    {
        os << "Index " << i << ":\n";
        os << vd.VoronoiCellEdges[i];
    }

    // 7. Print the two new maps

    // 7a. cellEdgeLookup
    os << "\ncellEdgeLookup ( (cellIndex, edgeIndex) -> VoronoiCellEdges index ):\n";
    for (const auto &kv : vd.cellEdgeLookup)
    {
        int cellIndex = kv.first.first;
        int edgeIndex = kv.first.second;
        int ceIdx = kv.second;
        os << "  ( " << cellIndex << ", " << edgeIndex << " ) -> " << ceIdx << "\n";
    }

    // 7b. segmentVertexPairToEdgeIndex
    os << "\nsegmentVertexPairToEdgeIndex ( (v1, v2) -> edgeIndex ):\n";
    for (const auto &kv : vd.segmentVertexPairToEdgeIndex)
    {
        int v1 = kv.first.first;
        int v2 = kv.first.second;
        int edgeIndex = kv.second;
        os << "  ( " << v1 << ", " << v2 << " ) -> " << edgeIndex << "\n";
    }

    // -------------------------------------------------------------------------
    // 8. Print the three maps that were missing:
    // -------------------------------------------------------------------------

    // 8a. delaunay_cell_to_voronoi_vertex_index
    os << "\nDelaunayCell -> VoronoiVertexIndex:\n";
    for (const auto &kv : vd.delaunay_cell_to_voronoi_vertex_index)
    {
        // kv.first is a Cell_handle, kv.second is an int
        // We'll print the address of the cell handle (or any ID we want).
        os << "  Cell_handle@" << &(*kv.first) << " -> " << kv.second << "\n";
    }

    // 8b. point_to_vertex_index
    os << "\nPoint -> VoronoiVertexIndex:\n";
    for (const auto &kv : vd.point_to_vertex_index)
    {
        // kv.first is a Point, kv.second is an int
        // We can usually print a CGAL::Point_3 directly, or explicitly by coords.
        os << "  " << kv.first << " -> " << kv.second << "\n";
    }

    // 8c. delaunayVertex_to_voronoiCell_index
    os << "\nDelaunayVertex -> VoronoiCellIndex:\n";
    for (const auto &kv : vd.delaunayVertex_to_voronoiCell_index)
    {
        // kv.first is a Vertex_handle, kv.second is an int
        // Similarly, we just print the handle as a raw pointer or any ID we prefer.
        os << "  Vertex_handle@" << &(*kv.first) << " -> " << kv.second << "\n";
    }

    return os;
}

//! @brief Stream operator for VoronoiCellEdge
template <typename OSTREAM_TYPE>
OSTREAM_TYPE &operator<<(OSTREAM_TYPE &os, const VoronoiCellEdge &ce)
{
    os << "VoronoiCellEdge:\n"
       << "  cellIndex:    " << ce.cellIndex << "\n"
       << "  edgeIndex:    " << ce.edgeIndex << "\n"
       << "  nextCellEdge: " << ce.nextCellEdge << "\n"
       << "  cycleIndices: ";

    if (ce.cycleIndices.empty())
    {
        os << "(none)";
    }
    else
    {
        os << "{ ";
        for (size_t i = 0; i < ce.cycleIndices.size(); ++i)
        {
            os << ce.cycleIndices[i];
            if (i + 1 < ce.cycleIndices.size())
                os << ", ";
        }
        os << " }";
    }
    os << "\n";
    return os;
}

#endif