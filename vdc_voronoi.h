//! @file vdc_voronoi.h
//! @brief Header file defining structures and functions for Voronoi Diagram construction, cycles, facets, and isosurfaces.

#ifndef VDC_VORONOI_H
#define VDC_VORONOI_H

#include "vdc_type.h"
#include "vdc_delaunay.h"

//! @brief Represents a vertex in a Voronoi diagram.
/*!
 * A Voronoi vertex is a point where multiple Voronoi edges meet. Each vertex
 * may belong to multiple Voronoi cells.
 */
struct VoronoiVertex
{
    Point vertex;                 //!< Geometric coordinates of the vertex.
    int index;                    //!<
    std::vector<int> cellIndices; //!< Indices of Voronoi cells that contain this vertex.

    //! @brief Constructor to initialize a Voronoi vertex.
    /*!
     * @param p Coordinates of the vertex.
     */
    VoronoiVertex(Point p) : vertex(p) {}
};

struct VoronoiEdge
{
    CGAL::Object edgeObject;
    int type; //!< 0 for segment, 1 for ray and 2 for lines; -1 for unknown

    //! @brief constructor
    /*!
     * @param obj instance of the edge in CGAL::Object
     */
    VoronoiEdge(CGAL::Object obj) : edgeObject(obj) {}
};

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
    std::vector<float> vertex_values;       //!< Scalar values at the facet's vertices.
    std::vector<int> vertex_values_indices; //!< The locations in voronoiDiagram.vertexValues for vertices in the facet
    int mirror_facet_index = -1; //!< Index of the mirror facet in the adjacent cell
    
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
    int voronoi_cell_index;                 //!< Index of the Voronoi cell this cycle belongs to;

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
    // Compute centroid of the cycle using the positions in vertices.
    void compute_centroid(const std::vector<VoronoiVertex> &vertices);
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
    int edgeIndex;                 //!< Index of this edge in the vector of edges in the VoronoiDiagram instance that contains the cell this egde is in
    std::vector<int> cycleIndices; //!< Indices of cycles in this VoronoiCell that corresponding to this edge
    int nextCellEdge;              //!< Index of next cell edge around the Voronoi Edge ( VoronoiDiagram.edges[edgeIndex])
};


struct pair_hash {
    std::size_t operator()(const std::pair<int, int>& p) const {
        auto h1 = std::hash<int>{}(p.first);
        auto h2 = std::hash<int>{}(p.second);
        return h1 ^ (h2 << 1); // Shift to reduce collisions
    }
};

//! @brief Represents the overall Voronoi diagram.
/*!
 * The Voronoi diagram consists of vertices, edges, cells, facets, and isosurface data.
 */
struct VoronoiDiagram
{
    std::vector<VoronoiVertex> vertices;    //!< List of Voronoi vertices in the diagram.
    std::vector<Object> edges;              //!< List of edges in the diagram (e.g., line segments).
    std::vector<std::pair<int, int>> edgeVertexIndices; //!< where each pair corresponds to an edge in edges. For segments, both indices are valid (>=0); for rays, one is -1; for lines, both are -1.
    std::vector<VoronoiCellEdge> cellEdges; //!< List of Cell Edges in the diagram
    std::vector<float> vertexValues;        //!< Scalar values at the Voronoi vertices.
    std::vector<VoronoiCell> cells;         //!< List of Voronoi cells in the diagram.
    std::vector<VoronoiFacet> facets;       //!< List of facets in the diagram.
    std::vector<int> oldToNewVertexIndex;  //!< Mapping from old to new vertex indices after collapse

    std::map<std::pair<int, int>, int> cellEdgeLookup;               //!< Maps (cellIndex, edgeIndex) -> index in cellEdges
    std::map<std::pair<int, int>, int> segmentVertexPairToEdgeIndex; //!< a map from a pair of Voronoi vertex indices (v_1, v_2) (in ascending order) to the edgeIndex in voronoiDiagram

    //! @brief Checks internal consistency of the VoronoiDiagram.
    void check() const;

    //! @brief Comprehensive checker for Voronoi diagram consistency.
    void checkAdvanced() const;

    //! @brief Collapse all Voronoi vertices closer than D, rebuild cells/facets, and re‐check
    void collapseSmallEdges(double D, CGAL::Epick::Iso_cuboid_3& bbox);

    int find_vertex(const Point &p);

private:
    //! @brief Verifies that `cellEdgeLookup` matches the data in `cellEdges`.
    void checkCellEdgeLookup() const;

    //! @brief Checks that each cellEdge's `nextCellEdge` points to another edge with the same `edgeIndex`.
    void checkNextCellEdgeConsistency() const;

    //! @brief Checks each VoronoiCell's facets to ensure that every facet's vertices are in the cell's vertex set.
    void checkCellFacets() const;

    //! @brief Helper methods for collapseSmallEdges
    void processEdges(std::vector<int>& mapto, double D);
    void mergeCloseVertices(double mergeTolerance, std::vector<int> &mapto);
    void compressMapping(std::vector<int> &mapto);
    void rebuildVertices(const std::vector<int> &mapto);
    void rebuildEdges();

    //! @brief Helper methods for checkAdvanced
    void checkFacetVertexCount() const;
    void checkCellFacetCount() const;
    void checkFacetCellCount() const;
    void checkEdgeFacetCount() const;
    void checkFacetNormals() const;
    void checkPairedFacetOrientations() const;

    //! @brief Helper methods for checkNextCellEdgeConsistency
    void checkNextCellEdgeValidity() const;
    void checkEdgeCycles() const;

    /// Hash a facet by its three smallest vertex indices (for “at most two cells” check)
    std::tuple<int,int,int> getFacetHashKey(const std::vector<int>& verts) const;

    /// Do two facet‐vertex‐sequences represent the same cyclic orientation?
    bool haveSameOrientation(const std::vector<int>& f1,
                             const std::vector<int>& f2) const;
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

struct RayKeyComparator {
    bool operator()(const std::pair<CGAL::Point_3<CGAL::Epick>, CGAL::Vector_3<CGAL::Epick>>& a,
                    const std::pair<CGAL::Point_3<CGAL::Epick>, CGAL::Vector_3<CGAL::Epick>>& b) const {
        // Compare points first
        if (a.first < b.first) return true;
        if (b.first < a.first) return false;
        // If points are equal, compare vectors lexicographically
        auto a_begin = a.second.cartesian_begin();
        auto a_end = a.second.cartesian_end();
        auto b_begin = b.second.cartesian_begin();
        auto b_end = b.second.cartesian_end();
        return std::lexicographical_compare(a_begin, a_end, b_begin, b_end);
    }
};

struct LineKeyComparator {
    bool operator()(const std::pair<CGAL::Point_3<CGAL::Epick>, CGAL::Vector_3<CGAL::Epick>>& a,
                    const std::pair<CGAL::Point_3<CGAL::Epick>, CGAL::Vector_3<CGAL::Epick>>& b) const {
        // Compare points first
        if (a.first < b.first) return true;
        if (b.first < a.first) return false;
        // If points are equal, compare vectors lexicographically
        auto a_begin = a.second.cartesian_begin();
        auto a_end = a.second.cartesian_end();
        auto b_begin = b.second.cartesian_begin();
        auto b_end = b.second.cartesian_end();
        return std::lexicographical_compare(a_begin, a_end, b_begin, b_end);
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
    os << "\nvertices:\n";
    for (size_t i = 0; i < vd.vertices.size(); ++i)
    {
        os << "Index " << i << ":\n";
        os << vd.vertices[i];
    }

    // 2. Voronoi Edges
    os << "\nVoronoiEdges:\n";
    for (const auto &edge : vd.edges)
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
    os << "\nvertexValues:\n";
    for (size_t i = 0; i < vd.vertexValues.size(); ++i)
    {
        os << "  Index " << i << ": " << vd.vertexValues[i] << "\n";
    }

    // 4. Voronoi Facets
    os << "\nVoronoiFacets:\n";
    for (size_t i = 0; i < vd.facets.size(); ++i)
    {
        os << "Index " << i << ":\n";
        os << vd.facets[i];
    }

    // 5. Voronoi Cells
    os << "\ncells:\n";
    for (const auto &cell : vd.cells)
    {
        os << "\n"
           << cell;
    }

    // 6. Voronoi CellEdges
    os << "\ncellEdges:\n";
    for (size_t i = 0; i < vd.cellEdges.size(); ++i)
    {
        os << "Index " << i << ":\n";
        os << vd.cellEdges[i];
    }

    // 7. Print the two new maps

    // 7a. cellEdgeLookup
    os << "\ncellEdgeLookup ( (cellIndex, edgeIndex) -> cellEdges index ):\n";
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

// Helper function declarations (internal linkage)
std::tuple<int, int, int> getFacetHashKey(const std::vector<int> &vertices);
bool haveSameOrientation(const std::vector<int> &facet1, const std::vector<int> &facet2);
#endif