//! @file vdc_voronoi.h
//! @brief Header file defining structures and functions for Voronoi Diagram construction, cycles, facets, and isosurfaces.

#ifndef VDC_VORONOI_H
#define VDC_VORONOI_H

#include "vdc_type.h"
#include "vdc_delaunay.h"
#include "vdc_utilities.h"
#include <utility> // for std::make_tuple

//! @brief Represents a vertex in a Voronoi diagram.
/*!
 * A Voronoi vertex is a point where multiple Voronoi edges meet. Each vertex
 * may belong to multiple Voronoi cells.
 */
struct VoronoiVertex
{
    Point coord;                  //!< Geometric coordinates of the vertex.
    int index;                    //!< The index of the vertex in the Voronoi diagram
    float value;                  //!< The scalar value of the vertex, used for isosurface extraction
    std::vector<int> cellIndices; //!< Indices of Voronoi cells that contain this vertex.

    //! @brief Constructor to initialize a Voronoi vertex.
    /*!
     * @param p Coordinates of the vertex.
     */
    VoronoiVertex(Point p) : coord(p) {}
};

struct VoronoiEdge
{
    CGAL::Object edgeObject;
    int vertex1; //!< Index of the first vertex ( -1 if infinite )
    int vertex2; //!< Index of the second vertex ( -1 if infinite )
    int type;    //!< 0 for segment, 1 for ray and 2 for lines; -1 for unknown
    int index = -1;
    Point source;                      //!< Source point for rays and lines; empty for segments
    Vector3 direction;                 //!< Direction vector for rays and lines; empty for segments
    std::vector<Facet> delaunayFacets; //!< Facet indices in the Delaunay triangulation that correspond to this edge

    //! @brief constructor
    /*!
     * @param obj instance of the edge in CGAL::Object
     */
    VoronoiEdge(CGAL::Object obj) : edgeObject(obj), vertex1(-1), vertex2(-1), type(-1), source(), direction(), delaunayFacets()
    {
    }

    //! @brief operator< for sorting Voronoi edges
    bool operator<(VoronoiEdge const &other) const
    {
        if (type != other.type)
        {
            return type < other.type;
        }
        if (type == 0)
        {
            // compare segment by its endpoint-indices
            return std::tie(vertex1, vertex2) < std::tie(other.vertex1, other.vertex2);
        }
        // for rays/lines: lexicographic on (source, direction)
        auto t1 = std::make_tuple(source.x(), source.y(), source.z(),
                                  direction.x(), direction.y(), direction.z());
        auto t2 = std::make_tuple(other.source.x(), other.source.y(), other.source.z(),
                                  other.direction.x(), other.direction.y(), other.direction.z());
        return t1 < t2;
    }
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
    bool is_bipolar = false;
};

//! @brief Represents a facet in a Voronoi diagram.
/*!
 * A Voronoi facet is defined by a list of vertices and their scalar values.
 * Facets are part of Voronoi cells.
 */
struct VoronoiCellFacet
{
    std::vector<int> vertices_indices;
    int mirror_facet_index = -1;  //!< Index of the mirror cell-facet in the adjacent cell.
    int facet_index = -1;         //!< Index of this cell-facet in facets (for compatibility).
    int voronoi_facet_index = -1; //!< Index of the unique facet in voronoi_facets.
    int orientation = 1;          //!< 1 if using vertices_indices as-is, -1 if reversed.    //! @brief Default constructor.
    VoronoiCellFacet() = default;
};

//! @brief Represents a facet in the Voronoi diagram (unique, shared by cells).
struct VoronoiFacet
{
    std::vector<int> vertices_indices;                                                      //!< Ordered indices of vertices forming this facet.
    int index = -1;                                                                         //!< Index of the facet in voronoi_facets.
    int primary_cell_facet_index = -1;                                                      //!< Index in facets of the primary cell-facet.
    BIPOLAR_MATCH_METHOD bipolar_match_method = BIPOLAR_MATCH_METHOD::UNDEFINED_MATCH_TYPE; //!< Matching method for bipolar edges.
    std::vector<int> edge_indices;                                                          // edge_indices[i] is the global edge index for local edge i (between vertices i and i+1)
    std::vector<std::pair<int, int>> bipolar_matches;                                       //!< Pairs of edge indices for matched bipolar edges.
    std::vector<int> bipolar_edge_indices;                                                  // Full-facet edge indices (0 to num-1) where bipolarstd::vector<int> bipolar_voronoi_edge_indices;
};

//! @brief Represents a closed cycle in a Voronoi cell formed by midpoints.
/*!
 * A cycle is a loop of midpoints connected by edges. Each cycle is associated
 * with a Voronoi cell and may have its centroid computed as the cycle's isovertex.
 */
struct Cycle
{
    std::vector<int> midpoint_indices;
    Point isovertex;                               //!< The computed isovertex (centroid) for this cycle.
    int voronoi_cell_index;                        //!< Index of the Voronoi cell.
    int facet_index;                               //!< Index of the cell facet this cycle belongs to.
    std::vector<int> bipolar_voronoi_edge_indices; // Global Voronoi edge indices for the paired bipolar edges
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

//! @brief Hash function for std::tuple<int,int,int>
/*!
 * Provides a hash function for tuples of three integers,
 * used for vertex indexing in Voronoi diagrams.
 */
struct TupleHash
{
    std::size_t operator()(const std::tuple<int, int, int> &t) const
    {
        std::size_t seed = 0;
        seed ^= std::hash<int>{}(std::get<0>(t)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= std::hash<int>{}(std::get<1>(t)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= std::hash<int>{}(std::get<2>(t)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

//! @brief Represents the overall Voronoi diagram.
/*!
 * The Voronoi diagram consists of vertices, edges, cells, facets, and isosurface data.
 */
struct VoronoiDiagram
{
    std::vector<VoronoiVertex> vertices;     //!< List of Voronoi vertices in the diagram.
    std::vector<VoronoiEdge> edges;          //!< List of edges in the diagram
    std::vector<VoronoiCellEdge> cellEdges;  //!< List of Cell Edges in the diagram
    std::vector<VoronoiCell> cells;          //!< List of Voronoi cells in the diagram.
    std::vector<VoronoiCellFacet> facets;    //!< List of facets in the diagram.
    std::vector<VoronoiFacet> global_facets; //!< List of unique Voronoi facets.

    std::map<std::pair<int, int>, int> cellEdgeLookup;               //!< Maps (cellIndex, edgeIndex) -> index in cellEdges
    std::map<std::pair<int, int>, int> segmentVertexPairToEdgeIndex; //!< a map from a pair of Voronoi vertex indices (v_1, v_2) (in ascending order) to the edgeIndex in voronoiDiagram

    // Member Functions

    void compute_bipolar_matches(float isovalue);
    void create_global_facets();
     
    //! @brief Matches bipolar edges on a specific facet using the given matching method.
    /*!
     * This routine pairs bipolar edges on the facet based on the provided BIPOLAR_MATCH_METHOD.
     * - For SEP_NEG: Pairs pos→neg edges to the next bipolar edge (typically neg→pos if alternating).
     * - For SEP_POS: Pairs neg→pos edges to the next bipolar edge (typically pos→neg if alternating).
     * - For UNCONSTRAINED_MATCH: Pairs consecutive bipolar edges arbitrarily (0-1, 2-3, etc.).
     * Pairs are stored in vf.bipolar_matches as indices into vf.bipolar_edge_indices.
     * Assumes types and bipolar_edge_indices are pre-computed and valid (even count, alternating).
     *
     * @param vf The VoronoiFacet to update with matches.
     * @param types Vector of types for each bipolar edge (pos→neg for true or neg→pos for false).
     * @param match_method The method to use for pairing.
     * @param num_bipolar Number of bipolar edges (types.size()).
     */
    void match_facet_bipolar_edges(VoronoiFacet &vf, const std::vector<bool> &types, BIPOLAR_MATCH_METHOD match_method, int num_bipolar);

    // Helper to get oriented vertices for a cell facet
    std::vector<int> get_vertices_for_facet(int cell_facet_index) const;

    //! @brief Adds a vertex to the Voronoi diagram.
    /*!
     * @param p The point coordinates of the vertex
     * @param value The scalar value associated with the vertex
     * @return Index of the newly added vertex
     */
    int AddVertex(const Point &p, float value = 0.0f);

    //! @brief Adds a segment edge to the Voronoi diagram.
    /*!
     * @param v1 Index of the first vertex
     * @param v2 Index of the second vertex
     * @param seg The CGAL segment object
     * @return Index of the newly added edge
     */
    int AddSegmentEdge(int v1, int v2, const Segment3 &seg);

    //! @brief Adds a ray edge to the Voronoi diagram.
    /*!
     * @param ray The CGAL ray object
     * @return Index of the newly added edge
     */
    int AddRayEdge(const Ray3 &ray);

    //! @brief Adds a line edge to the Voronoi diagram.
    /*!
     * @param line The CGAL line object
     * @return Index of the newly added edge
     */
    int AddLineEdge(const Line3 &line);

    //! @brief Adds a facet to the Voronoi diagram.
    /*!
     * @param vertices_indices Indices of vertices forming the facet
     * @return Index of the newly added facet
     */
    int AddCellFacet(const std::vector<int> &vertices_indices);

    //! @brief Adds a cell to the Voronoi diagram.
    /*!
     * @param delaunay_vertex Handle to the corresponding Delaunay vertex
     * @return Index of the newly added cell
     */
    int AddCell(Vertex_handle delaunay_vertex);

    //! @brief Checks internal consistency of the VoronoiDiagram.
    void check(bool check_norm) const;

    //! @brief Comprehensive checker for Voronoi diagram consistency.
    void checkAdvanced(bool check_norm) const;

    /// Do two facet‐vertex‐sequences represent the same cyclic orientation?
    bool haveSameOrientation(const std::vector<int> &f1,
                             const std::vector<int> &f2) const;

    bool haveOppositeOrientation(const std::vector<int> &f1,
                                 const std::vector<int> &f2) const;

private:
    //! @brief Verifies that `cellEdgeLookup` matches the data in `cellEdges`.
    void checkCellEdgeLookup() const;

    //! @brief Checks that each cellEdge's `nextCellEdge` points to another edge with the same `edgeIndex`.
    void checkNextCellEdgeConsistency() const;

    //! @brief Checks each VoronoiCell's facets to ensure that every facet's vertices are in the cell's vertex set.
    void checkCellFacets() const;

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
    std::tuple<int, int, int> getFacetHashKey(const std::vector<int> &verts) const;
};

// Add standalone function declaration at the end of vdc_voronoi.h
VoronoiDiagram collapseSmallEdges(const VoronoiDiagram &vd, double D, const CGAL::Epick::Iso_cuboid_3 &bbox, Delaunay &dt);

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

//! @brief Comparator for approximate point equality
/*!
 * Compares points with a tolerance (epsilon) to account for
 * floating-point inaccuracies in geometric computations.
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

//! @brief Comparator for ray keys (point + direction)
/*!
 * Used to order rays lexicographically by their source point
 * and direction vector for efficient storage and retrieval.
 */
struct RayKeyComparator
{
    bool operator()(const std::pair<CGAL::Point_3<CGAL::Epick>, CGAL::Vector_3<CGAL::Epick>> &a,
                    const std::pair<CGAL::Point_3<CGAL::Epick>, CGAL::Vector_3<CGAL::Epick>> &b) const
    {
        // Compare points first
        if (a.first < b.first)
            return true;
        if (b.first < a.first)
            return false;
        // If points are equal, compare vectors lexicographically
        auto a_begin = a.second.cartesian_begin();
        auto a_end = a.second.cartesian_end();
        auto b_begin = b.second.cartesian_begin();
        auto b_end = b.second.cartesian_end();
        return std::lexicographical_compare(a_begin, a_end, b_begin, b_end);
    }
};

//! @brief Comparator for line keys (point + direction)
/*!
 * Used to order lines lexicographically by their point and
 * direction vector for efficient storage and retrieval.
 */
struct LineKeyComparator
{
    bool operator()(const std::pair<CGAL::Point_3<CGAL::Epick>, CGAL::Vector_3<CGAL::Epick>> &a,
                    const std::pair<CGAL::Point_3<CGAL::Epick>, CGAL::Vector_3<CGAL::Epick>> &b) const
    {
        // Compare points first
        if (a.first < b.first)
            return true;
        if (b.first < a.first)
            return false;
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
OSTREAM_TYPE &operator<<(OSTREAM_TYPE &os, const VoronoiCellFacet &vf)
{
    os << "VoronoiCellFacet:\n";
    os << " Index: " << vf.facet_index << "\n";
    os << " Global facet index: " << vf.voronoi_facet_index << "\n";
    os << " Orientation: " << vf.orientation << "\n";
    os << "\n";

    return os;
}

//! @brief Overloaded output operator for Cycle
template <typename OSTREAM_TYPE>
OSTREAM_TYPE &operator<<(OSTREAM_TYPE &os, const Cycle &cycle)
{
    os << "Cycle:\n";
    os << "  Isovertex: " << cycle.isovertex << "\n";
    os << "  Voronoi cell index: " << cycle.voronoi_cell_index << "\n";
    os << "  Facet index: " << cycle.facet_index << "\n";
    return os;
}

//! @brief Overloaded output operator for VoronoiVertex
template <typename OSTREAM_TYPE>
OSTREAM_TYPE &operator<<(OSTREAM_TYPE &os, const VoronoiVertex &vv)
{
    os << "VoronoiVertex:\n";
    os << "  Vertex: " << vv.coord << "\n";

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

// Helper to print BIPOLAR_MATCH_METHOD
std::string matchMethodToString(BIPOLAR_MATCH_METHOD method);

// Overloaded output operator for VoronoiFacet
template <typename OSTREAM_TYPE>
OSTREAM_TYPE &operator<<(OSTREAM_TYPE &os, const VoronoiFacet &vf)
{
    os << "VoronoiFacet:\n";
    os << "  Index: " << vf.index << "\n";
    os << "  Primary cell facet index: " << vf.primary_cell_facet_index << "\n";
    os << "  Bipolar match method: " << matchMethodToString(vf.bipolar_match_method) << "\n";

    os << "  Vertices indices: ";
    for (const int idx : vf.vertices_indices)
    {
        os << idx << " ";
    }
    os << "\n";

    os << "  Bipolar matches: ";
    if (vf.bipolar_matches.empty())
    {
        os << "(none)";
    }
    else
    {
        for (const auto &match : vf.bipolar_matches)
        {
            os << "(" << match.first << ", " << match.second << ") ";
        }
    }
    os << "\n";

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

        if (CGAL::assign(segment, edge.edgeObject))
        {
            os << "Segment(" << segment.source() << " - " << segment.target() << ")\n";
        }
        else if (CGAL::assign(line, edge.edgeObject))
        {
            os << "Line(" << line.point(0) << " - " << line.point(1) << ")\n";
        }
        else if (CGAL::assign(ray, edge.edgeObject))
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
    for (size_t i = 0; i < vd.vertices.size(); ++i)
    {
        os << "  Index " << i << ": " << vd.vertices[i].value << "\n";
    }

    // 4. Voronoi Cell Facets (existing)
    os << "\nVoronoiCellFacets:\n";
    for (size_t i = 0; i < vd.facets.size(); ++i)
    {
        os << "Index " << i << ":\n";
        os << vd.facets[i];
    }

    // 5. Voronoi Global Facets (new addition)
    os << "\nVoronoiGlobalFacets:\n";
    for (size_t i = 0; i < vd.global_facets.size(); ++i)
    {
        os << "Index " << i << ":\n";
        os << vd.global_facets[i];
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

inline std::ostream &operator<<(std::ostream &os, IsoSurface const &iso)
{
    os << "IsoSurface: "
       << iso.isosurfaceVertices.size() << " verts, "
       << iso.isosurfaceTrianglesMulti.size() << " tris\n";

    os << "  Vertices:\n";
    for (auto const &p : iso.isosurfaceVertices)
    {
        // CGAL::Point_3 already has operator<<, prints e.g. "(x,y,z)"
        os << "    " << p << "\n";
    }

    os << "  Multi-triangles:\n";
    for (auto const &t : iso.isosurfaceTrianglesMulti)
    {
        os << "    ("
           << std::get<0>(t) << ", "
           << std::get<1>(t) << ", "
           << std::get<2>(t) << ")\n";
    }

    return os;
}

// Helper function declarations (internal linkage)
std::tuple<int, int, int> getFacetHashKey(const std::vector<int> &vertices);
bool haveSameOrientation(const std::vector<int> &facet1, const std::vector<int> &facet2);

void write_voronoiDiagram(VoronoiDiagram &vd, std::string &output_filename);
#endif