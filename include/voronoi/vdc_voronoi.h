//! @file vdc_voronoi.h
//! @brief Header file defining structures and functions for Voronoi Diagram construction, cycles, facets, and isosurfaces.

#ifndef VDC_VORONOI_H
#define VDC_VORONOI_H

#include "core/vdc_type.h"
#include "delaunay/vdc_delaunay.h"
#include "core/vdc_utilities.h"
#include <array>
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

    //! @brief Print the Voronoi vertex information for debugging
    template <typename OSTREAM_TYPE>
    void Print(OSTREAM_TYPE & out) const {
        out << "VoronoiVertex:\n";
        out << "  Coordinate: " << coord << "\n";
        out << "  Index: " << index << "\n";
        out << "  Scalar value: " << value << "\n";
        out << "  Cell indices: [";
        for (size_t i = 0; i < cellIndices.size(); ++i) {
            out << cellIndices[i];
            if (i + 1 < cellIndices.size()) out << ", ";
        }
        out << "]\n";
    }
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
    std::vector<Facet> delaunayFacets; //!< Facet in the Delaunay triangulation that correspond to this edge

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

    //! @brief Print basic edge information
    template <typename OSTREAM_TYPE>
    void PrintBasic(OSTREAM_TYPE & out) const {
        out << "  Index: " << index << "\n";
        out << "  Type: " << type << " (0=segment, 1=ray, 2=line, -1=unknown)\n";
        out << "  Vertex1: " << vertex1 << "\n";
        out << "  Vertex2: " << vertex2 << "\n";
    }

    //! @brief Print geometric information
    template <typename OSTREAM_TYPE>
    void PrintGeometry(OSTREAM_TYPE & out) const {
        if (type == 0) {
            // Segment
            Segment3 segment;
            if (CGAL::assign(segment, edgeObject)) {
                out << "  Segment: " << segment.source() << " -> " << segment.target() << "\n";
            }
        } else if (type == 1) {
            // Ray
            out << "  Ray source: " << source << "\n";
            out << "  Ray direction: " << direction << "\n";
        } else if (type == 2) {
            // Line
            out << "  Line source: " << source << "\n";
            out << "  Line direction: " << direction << "\n";
        }
    }

    //! @brief Print Delaunay facets associated with this edge
    template <typename OSTREAM_TYPE>
    void PrintDelaunayFacets(OSTREAM_TYPE & out) const {
        out << "  Delaunay facets: " << delaunayFacets.size() << " facet(s)\n";
    }

    //! @brief Print complete edge information for debugging
    template <typename OSTREAM_TYPE>
    void Print(OSTREAM_TYPE & out) const {
        out << "VoronoiEdge:\n";
        PrintBasic(out);
        PrintGeometry(out);
        PrintDelaunayFacets(out);
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
    Point point; //!< Geometric coordinates of the midpoint.
    // TODO: Change this into a tuple ( should expect fixed size of 2)
    std::vector<int> connected_to; //!< Indices of midpoints connected to this one, forming graph edges.
    int facet_index;               //!< Index of the facet this midpoint lies on.
    int cycle_index;               //!< Index of the cycle this midpoint belongs to.
    int global_edge_index = -1;
    bool is_bipolar = false;

    //! @brief Print midpoint node information for debugging
    template <typename OSTREAM_TYPE>
    void Print(OSTREAM_TYPE & out) const {
        out << "MidpointNode:\n";
        out << "  Point: " << point << "\n";
        out << "  Connected to: [";
        for (size_t i = 0; i < connected_to.size(); ++i) {
            out << connected_to[i];
            if (i + 1 < connected_to.size()) out << ", ";
        }
        out << "]\n";
        out << "  Facet index: " << facet_index << "\n";
        out << "  Cycle index: " << cycle_index << "\n";
        out << "  Global edge index: " << global_edge_index << "\n";
        out << "  Is bipolar: " << (is_bipolar ? "true" : "false") << "\n";
    }
};

//! @brief Represents a facet in a Voronoi diagram.
/*!
 * A Voronoi facet is defined by a list of vertices and their scalar values.
 * Facets are part of Voronoi cells.
 */
struct VoronoiCellFacet
{
    std::vector<int> vertices_indices;
    std::vector<int> cell_edge_indices; // ordered edges (indices into VoronoiDiagram.cellEdges)

    int mirror_facet_index = -1;  //!< Index of the mirror cell-facet in the adjacent cell.
    int facet_index = -1;         //!< Index of this cell-facet in facets (for compatibility).
    int voronoi_facet_index = -1; //!< Index of the unique facet in voronoi_facets.
    int orientation = 1;          //!< 1 if using vertices_indices as-is, -1 if reversed.    //! @brief Default constructor.
    VoronoiCellFacet() = default;

    //! @brief Print vertex indices
    template <typename OSTREAM_TYPE>
    void PrintVertices(OSTREAM_TYPE & out) const {
        out << "  Vertex indices: [";
        for (size_t i = 0; i < vertices_indices.size(); ++i) {
            out << vertices_indices[i];
            if (i + 1 < vertices_indices.size()) out << ", ";
        }
        out << "]\n";
    }

    //! @brief Print cell edge indices
    template <typename OSTREAM_TYPE>
    void PrintCellEdgeIndices(OSTREAM_TYPE & out) const {
        out << "  Cell edge indices: [";
        for (size_t i = 0; i < cell_edge_indices.size(); ++i) {
            out << cell_edge_indices[i];
            if (i + 1 < cell_edge_indices.size()) out << ", ";
        }
        out << "]\n";
    }

    //! @brief Print metadata
    template <typename OSTREAM_TYPE>
    void PrintMetadata(OSTREAM_TYPE & out) const {
        out << "  Mirror facet index: " << mirror_facet_index << "\n";
        out << "  Facet index: " << facet_index << "\n";
        out << "  Voronoi facet index: " << voronoi_facet_index << "\n";
        out << "  Orientation: " << orientation << "\n";
    }

    //! @brief Print complete cell facet information for debugging
    template <typename OSTREAM_TYPE>
    void Print(OSTREAM_TYPE & out) const {
        out << "VoronoiCellFacet:\n";
        PrintVertices(out);
        PrintCellEdgeIndices(out);
        PrintMetadata(out);
    }
};

// One iso-segment per bipolar match inside a global facet.
// No geometry stored; this is bookkeeping only.
struct IsoSegment
{
    int global_facet_index = -1; // vfi
    int slotA = -1;              // boundary slot in global facet
    int slotB = -1;              // boundary slot in global facet
    int edgeA = -1;              // global Voronoi edge id for slotA (optional, for logging/asserts)
    int edgeB = -1;              // global Voronoi edge id for slotB
    int comp[2] = {-1, -1};      // component (=cycle id) in incident cells C0/C1; -1 if none/unknown

    //! @brief Print iso-segment information for debugging
    template <typename OSTREAM_TYPE>
    void Print(OSTREAM_TYPE & out) const {
        out << "IsoSegment:\n";
        out << "  Global facet index: " << global_facet_index << "\n";
        out << "  Slot A: " << slotA << "\n";
        out << "  Slot B: " << slotB << "\n";
        out << "  Edge A: " << edgeA << "\n";
        out << "  Edge B: " << edgeB << "\n";
        out << "  Components: [" << comp[0] << ", " << comp[1] << "]\n";
    }
};

//! @brief Represents a facet in the Voronoi diagram (unique, shared by cells).
struct VoronoiFacet
{
    std::vector<int> vertices_indices;     //!< Ordered indices of vertices forming this facet.
    std::vector<int> voronoi_edge_indices; //!< ordered edges along boundary; k -> edge(v[k], v[k+1])
    int index = -1;                        //!< Index of the facet in voronoi_facets.
    int primary_cell_facet_index = -1;     //!< Index in facets of the primary cell-facet.

    BIPOLAR_MATCH_METHOD bipolar_match_method = BIPOLAR_MATCH_METHOD::UNDEFINED_MATCH_TYPE; //!< Matching method for bipolar edges.

    std::vector<std::pair<int, int>> bipolar_matches; //!< Pairs of edge indices for matched bipolar edges.
    std::vector<int> bipolar_edge_indices;            //!< Boundary slot indices (0..m-1) around the facet that are bipolar

    // NEW
    std::array<int, 2> incident_cell_indices = {-1, -1};       // {C0, C1} or {-1, X}
    std::array<int, 2> incident_cell_facet_indices = {-1, -1}; // cell-facet indices for C0/C1
    std::vector<IsoSegment> iso_segments;

    //! @brief Print Voronoi facet information for debugging
    template <typename OSTREAM_TYPE>
    void Print(OSTREAM_TYPE & out) const {
        out << "VoronoiFacet:\n";
        out << "  Index: " << index << "\n";
        out << "  Primary cell facet index: " << primary_cell_facet_index << "\n";
        out << "  Vertex indices: [";
        for (size_t i = 0; i < vertices_indices.size(); ++i) {
            out << vertices_indices[i];
            if (i + 1 < vertices_indices.size()) out << ", ";
        }
        out << "]\n";
        out << "  Voronoi edge indices: [";
        for (size_t i = 0; i < voronoi_edge_indices.size(); ++i) {
            out << voronoi_edge_indices[i];
            if (i + 1 < voronoi_edge_indices.size()) out << ", ";
        }
        out << "]\n";
        out << "  Bipolar match method: " << static_cast<int>(bipolar_match_method) << "\n";
        out << "  Bipolar edge indices: [";
        for (size_t i = 0; i < bipolar_edge_indices.size(); ++i) {
            out << bipolar_edge_indices[i];
            if (i + 1 < bipolar_edge_indices.size()) out << ", ";
        }
        out << "]\n";
        out << "  Bipolar matches: [";
        for (size_t i = 0; i < bipolar_matches.size(); ++i) {
            out << "(" << bipolar_matches[i].first << "," << bipolar_matches[i].second << ")";
            if (i + 1 < bipolar_matches.size()) out << ", ";
        }
        out << "]\n";
        out << "  Incident cells: [" << incident_cell_indices[0] << ", " << incident_cell_indices[1] << "]\n";
        out << "  Incident cell facets: [" << incident_cell_facet_indices[0] << ", " << incident_cell_facet_indices[1] << "]\n";
        out << "  Iso segments count: " << iso_segments.size() << "\n";
    }
};

//! @brief Represents a closed cycle in a Voronoi cell formed by midpoints.
/*!
 * A cycle is a loop of midpoints connected by edges. Each cycle is associated
 * with a Voronoi cell and may have its centroid computed as the cycle's isovertex.
 */
struct Cycle
{
    std::vector<int> midpoint_indices;
    std::vector<std::pair<int, int>> edges;        //!< Edges forming the cycle, represented by pairs of midpoint indices.
    Point isovertex;                               //!< The computed isovertex (centroid) for this cycle.
    int voronoi_cell_index;                        //!< Index of the Voronoi cell.
    int facet_index;                               //!< Index of the cell facet this cycle belongs to.
    std::vector<int> bipolar_voronoi_edge_indices; // Global Voronoi edge indices for the paired bipolar edges

    void compute_centroid(const std::vector<MidpointNode> &midpoints);

    //! @brief Print cycle information for debugging
    template <typename OSTREAM_TYPE>
    void Print(OSTREAM_TYPE & out) const {
        out << "Cycle:\n";
        out << "  Isovertex: " << isovertex << "\n";
        out << "  Voronoi cell index: " << voronoi_cell_index << "\n";
        out << "  Facet index: " << facet_index << "\n";
        out << "  Midpoint indices: [";
        for (size_t i = 0; i < midpoint_indices.size(); ++i) {
            out << midpoint_indices[i];
            if (i + 1 < midpoint_indices.size()) out << ", ";
        }
        out << "]\n";
        out << "  Edges: [";
        for (size_t i = 0; i < edges.size(); ++i) {
            out << "(" << edges[i].first << "," << edges[i].second << ")";
            if (i + 1 < edges.size()) out << ", ";
        }
        out << "]\n";
        out << "  Bipolar voronoi edge indices: [";
        for (size_t i = 0; i < bipolar_voronoi_edge_indices.size(); ++i) {
            out << bipolar_voronoi_edge_indices[i];
            if (i + 1 < bipolar_voronoi_edge_indices.size()) out << ", ";
        }
        out << "]\n";
    }
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

    //! @brief Print Voronoi cell information for debugging
    template <typename OSTREAM_TYPE>
    void Print(OSTREAM_TYPE & out) const {
        out << "VoronoiCell:\n";
        out << "  Cell index: " << cellIndex << "\n";
        out << "  Delaunay vertex: " << delaunay_vertex->point() << "\n";
        out << "  Vertices indices: [";
        for (size_t i = 0; i < vertices_indices.size(); ++i) {
            out << vertices_indices[i];
            if (i + 1 < vertices_indices.size()) out << ", ";
        }
        out << "]\n";
        out << "  Facet indices: [";
        for (size_t i = 0; i < facet_indices.size(); ++i) {
            out << facet_indices[i];
            if (i + 1 < facet_indices.size()) out << ", ";
        }
        out << "]\n";
        out << "  IsoVertex start index: " << isoVertexStartIndex << "\n";
        out << "  Number of isoVertices: " << numIsoVertices << "\n";
        out << "  Cycles count: " << cycles.size() << "\n";
    }
};

struct VoronoiCellEdge
{
    int cellIndex;                 //!< Index of the VoronoiCell that contains this CellEdge
    int edgeIndex;                 //!< Index of this edge in the vector of edges in the VoronoiDiagram instance that contains the cell this egde is in
    std::vector<int> cycleIndices; //!< Indices of cycles in this VoronoiCell that corresponding to this edge
    int nextCellEdge;              //!< Index of next cell edge around the Voronoi Edge ( VoronoiDiagram.edges[edgeIndex])

    //! @brief Print Voronoi cell edge information for debugging
    template <typename OSTREAM_TYPE>
    void Print(OSTREAM_TYPE & out) const {
        out << "VoronoiCellEdge:\n";
        out << "  Cell index: " << cellIndex << "\n";
        out << "  Edge index: " << edgeIndex << "\n";
        out << "  Next cell edge: " << nextCellEdge << "\n";
        out << "  Cycle indices: [";
        for (size_t i = 0; i < cycleIndices.size(); ++i) {
            out << cycleIndices[i];
            if (i + 1 < cycleIndices.size()) out << ", ";
        }
        out << "]\n";
    }
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
    std::vector<VoronoiVertex> vertices;    //!< List of Voronoi vertices in the diagram.
    std::vector<VoronoiEdge> edges;         //!< List of edges in the diagram
    std::vector<VoronoiCellEdge> cellEdges; //!< List of Cell Edges in the diagram
    std::vector<VoronoiCell> cells;         //!< List of Voronoi cells in the diagram.
    // TODO: Rename CellFacet->cell_facets ; Rename global_facets->surface_facets
    std::vector<VoronoiCellFacet> facets;    //!< List of facets in the diagram.
    std::vector<VoronoiFacet> global_facets; //!< List of unique Voronoi facets.

    std::map<std::pair<int, int>, int> cellEdgeLookup;               //!< Maps (cellIndex, edgeIndex) -> index in cellEdges
    std::map<std::pair<int, int>, int> segmentVertexPairToEdgeIndex; //!< a map from a pair of Voronoi vertex indices (v_1, v_2) (in ascending order) to the edgeIndex in voronoiDiagram

    // Member Functions

    void compute_bipolar_matches(float isovalue);
    void create_global_facets();

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

    //! @brief Print Voronoi diagram information for debugging
    template <typename OSTREAM_TYPE>
    void Print(OSTREAM_TYPE & out) const {
        out << "VoronoiDiagram:\n";
        out << "  Vertices: " << vertices.size() << " vertex(es)\n";
        out << "  Edges: " << edges.size() << " edge(s)\n";
        out << "  Cell edges: " << cellEdges.size() << " cell edge(s)\n";
        out << "  Cells: " << cells.size() << " cell(s)\n";
        out << "  Cell facets: " << facets.size() << " cell facet(s)\n";
        out << "  Global facets: " << global_facets.size() << " global facet(s)\n";
        out << "  Cell edge lookup entries: " << cellEdgeLookup.size() << "\n";
        out << "  Segment vertex pair to edge index entries: " << segmentVertexPairToEdgeIndex.size() << "\n";
    }

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
    std::vector<int> getFacetHashKey(const std::vector<int> &verts) const;
};

//! @brief Classifies bipolar boundary edges within a global facet.
/*!
 * Scans a facet’s boundary edges and records which slots are bipolar at the
 * given isovalue, along with their sign orientation (−→+ vs +→−).
 *
 * @param vd Voronoi diagram providing vertex values
 * @param vf The global facet to analyze
 * @param isovalue Threshold for bipolarity
 * @param out_bipolar_edge_indices Output: facet‑local edge slot indices that are bipolar
 * @param out_is_pos_dir Output: per‑slot booleans; true for (−→+), false for (+→−)
 */
static void collect_facet_bipolars(
    const VoronoiDiagram &vd,
    const VoronoiFacet &vf,
    float isovalue,
    std::vector<int> &out_bipolar_edge_indices,
    std::vector<bool> &out_is_pos_dir);

//! @brief Pairs bipolar edges along a facet boundary according to a method.
/*!
 * Produces disjoint pairs of facet‑local edge slots describing how cycles should
 * connect across the facet.
 *
 * @param bipolar_edge_indices Indices of bipolar slots along facet boundary
 * @param is_pos_dir Sign for each slot; true for (−→+), false for (+→−)
 * @param method Pairing strategy (SEP_POS, SEP_NEG, UNCONSTRAINED_MATCH)
 * @param out_pairs Output pairs (indices into `bipolar_edge_indices`)
 */
static void match_facet_bipolar_edges(
    const std::vector<int> &bipolar_edge_indices,
    const std::vector<bool> &is_pos_dir,
    BIPOLAR_MATCH_METHOD method,
    std::vector<std::pair<int, int>> &out_pairs);

//! @brief Recomputes bipolar pairing for a single global facet.
/*!
 * Re-evaluates bipolar slots and updates the facet’s `bipolar_matches` using
 * its current `bipolar_match_method`.
 *
 * @param vd Voronoi diagram owning the facet
 * @param vfi Index of the global facet
 * @param isovalue Threshold for bipolarity
 */
void recompute_bipolar_matches_for_facet(VoronoiDiagram &vd, int vfi, float isovalue);

//! @brief Collapses edges shorter than a threshold and rebuilds a diagram.
/*!
 * Merges vertices connected by very short finite segment edges and reconstructs
 * the Voronoi diagram and auxiliary structures while preserving index stability
 * where possible.
 *
 * @param vd Input Voronoi diagram
 * @param D Length threshold (edges with length < D collapse)
 * @param bbox Clipping box (not used during collapse but kept for symmetry)
 * @param dt Delaunay triangulation handle (not modified)
 * @param out_vertex_mapping Output parameter: maps old vertex indices to new vertex indices after collapse
 * @return A new Voronoi diagram with small edges collapsed
 */
VoronoiDiagram collapseSmallEdges(const VoronoiDiagram &vd, double D, const CGAL::Epick::Iso_cuboid_3 &bbox, Delaunay &dt, std::vector<int> &out_vertex_mapping);

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
    std::vector<int> triangleSourceEdges;                            //!< Source Voronoi edge index for each multi-isov triangle (debug/fixup).
    std::array<double, 3> vertex_scale{1.0, 1.0, 1.0};               //!< Per-axis scale to convert from grid units to physical space.

    //! @brief Print isosurface information for debugging
    template <typename OSTREAM_TYPE>
    void Print(OSTREAM_TYPE & out) const {
        out << "IsoSurface:\n";
        out << "  Vertices: " << isosurfaceVertices.size() << " vertex(es)\n";
        out << "  Multi-isov triangles: " << isosurfaceTrianglesMulti.size() << " triangle(s)\n";
        out << "  Single-isov triangles: " << isosurfaceTrianglesSingle.size() << " triangle(s)\n";
        out << "  Triangle source edges: " << triangleSourceEdges.size() << " entry(es)\n";
        out << "  Vertex scale: [" << vertex_scale[0] << ", " << vertex_scale[1] << ", " << vertex_scale[2] << "]\n";
    }
};

//! @brief Represents a midpoint on an edge, along with its facet information.
struct EdgeMidpoint
{
    Point midpoint;  //!< Geometric coordinates of the midpoint.
    int facet_index; //!< Index of the facet this midpoint belongs to.

    //! @brief Print edge midpoint information for debugging
    template <typename OSTREAM_TYPE>
    void Print(OSTREAM_TYPE & out) const {
        out << "EdgeMidpoint:\n";
        out << "  Midpoint: " << midpoint << "\n";
        out << "  Facet index: " << facet_index << "\n";
    }
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

    //! @brief Print Delaunay vertex information for debugging
    template <typename OSTREAM_TYPE>
    void Print(OSTREAM_TYPE & out) const {
        out << "DelaunayVertex:\n";
        out << "  Point: " << point << "\n";
        out << "  Is dummy: " << (is_dummy ? "true" : "false") << "\n";
    }
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
    os << "  midpoints indices: ";
    for (int i = 0; i < cycle.midpoint_indices.size(); i++)
    {
        os << i << " ";
    }
    return os;
}

//! @brief Overloaded output operator for VoronoiVertex
template <typename OSTREAM_TYPE>
OSTREAM_TYPE &operator<<(OSTREAM_TYPE &os, const VoronoiVertex &vv)
{
    os << "VoronoiVertex:\n";
    os << " Vertex: " << vv.coord << "\n";
    os << " Scalar Value: " << vv.value << "\n";

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

//! @brief Human‑readable string for a bipolar matching method.
/*!
 * @param method Enum value
 * @return String name such as "SEP_POS", "SEP_NEG"
 */
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

    os << "  Edges: ";
    for (const auto &edge : vf.voronoi_edge_indices)
    {
        os << "(" << edge << ")";
    }
    os << "\n";

    os << "  Bipolar Edges: [ ";
    for (const auto &be : vf.bipolar_edge_indices)
    {
        os << be << " ";
    }
    os << "]\n  Bipolar Matches: [ ";
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
    os << "]\n";

    return os;
}

template <typename OSTREAM_TYPE>
OSTREAM_TYPE &operator<<(OSTREAM_TYPE &os, const VoronoiEdge &ve)
{
    Segment3 segment;
    Line3 line;
    Ray3 ray;

    if (CGAL::assign(segment, ve.edgeObject))
    {
        os << "Segment(" << segment.source() << " - " << segment.target() << ")\n";
    }
    else if (CGAL::assign(line, ve.edgeObject))
    {
        os << "Line(" << line.point(0) << " - " << line.point(1) << ")\n";
    }
    else if (CGAL::assign(ray, ve.edgeObject))
    {
        os << "Ray(" << ray.source() << ", direction: " << ray.direction() << ")\n";
    }
    else
    {
        os << "Unknown edge type.\n";
    }
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

//! @brief Populates incident cell references for each global facet.
/*!
 * Fills `incident_cell_indices` and `incident_cell_facet_indices` in
 * `VoronoiFacet` so later passes can map facet slots back to specific cells.
 *
 * @param vd Voronoi diagram to annotate
 */
void populate_incident_cells_for_global_facets(VoronoiDiagram &vd);

// Map a global facet slot index to the local slot index in a given cell facet,
// accounting for orientation (CW/CCW) of that cell facet w.r.t. the global facet.
//! @brief Maps a facet‑global slot index to a cell‑local slot index.
/*!
 * Accounts for whether the cell facet’s orientation matches or opposes the
 * global facet’s orientation.
 *
 * @param vf Global facet (provides canonical order)
 * @param cf Cell facet instance
 * @param slot_global Boundary slot index in the global facet
 * @return Corresponding local slot index in the cell facet
 */
int map_global_slot_to_cell(const VoronoiFacet &vf,
                            const VoronoiCellFacet &cf,
                            int slot_global);

// Return the cycle id for the bipolar edge at given local slot in cell facet cf
// (single-slot lookup via cf.cell_edge_indices[local_slot]; NO hashing).
//! @brief Looks up the cycle id touched by a bipolar edge slot in a cell facet.
/*!
 * Uses the `cellEdgeLookup` and per‑edge rings to find a cycle index within a
 * given cell that corresponds to the specified bipolar edge.
 *
 * @param vd Voronoi diagram
 * @param cellIndex Index of the incident cell
 * @param cellFacetIndex Index of the cell facet within `vd.facets`
 * @param slot_global Slot index on the global facet boundary
 * @return Cycle id local to the cell (−1 if none)
 */
int find_cycle_for_bipolar_edge(const VoronoiDiagram &vd,
                                int cellIndex,
                                int cellFacetIndex,
                                int slot_global);

// Build vf.iso_segments and fill comp[2] via single-slot lookup per side.
// Pre: bipolar_matches exist; cellEdges[*].cycleIndices filled.
//! @brief Builds iso‑segments per global facet from current bipolar matches.
/*!
 * Creates bookkeeping segments for each paired bipolar slot and records which
 * per‑cell cycle they touch on both sides when available.
 *
 * @param vd Voronoi diagram
 * @param vfi Global facet index
 * @param isovalue Current isovalue threshold
 */
void build_iso_segments_for_facet(VoronoiDiagram &vd,
                                  int vfi,
                                  float isovalue);

// Detect whether a facet has duplicate (comp0, comp1) among its iso-segments.
// Returns true if problematic; optionally returns one offending pair index.
//! @brief Detects duplicate (comp0,comp1) pairs among a facet’s iso‑segments.
/*!
 * Reports whether the iso‑segments imply incompatible cycle connectivity across
 * the facet and optionally returns one offending pair.
 *
 * @param vd Voronoi diagram
 * @param vfi Global facet index
 * @param offending_pair Optional output pair of indices into `iso_segments`
 * @return true if a problematic duplicate is detected
 */
bool facet_has_problematic_iso_segments(const VoronoiDiagram &vd,
                                        int vfi,
                                        std::pair<int, int> *offending_pair = nullptr);

// Flip matching method for a facet (SEP_POS <-> SEP_NEG).
//! @brief Toggles a facet’s matching method between SEP_POS and SEP_NEG.
/*!
 * Leaves UNCONSTRAINED_MATCH/UNDEFINED unchanged.
 * @param vf Global facet to update
 */
void flip_bipolar_match_method(VoronoiFacet &vf);

// Recompute cycles for a single cell after facet-local rematching.
// Clears tags in that cell, rebuilds midpoint graph using current global matches,
// extracts cycles, and writes cycle tags to vd.cellEdges[*].cycleIndices.
//! @brief Rebuilds midpoint graph and cycles for a single cell.
/*!
 * Clears old cycle tags and repopulates `cellEdges[*].cycleIndices` using the
 * current global facet matches.
 *
 * @param vd Voronoi diagram
 * @param cellIndex Cell to recompute
 * @param isovalue Isovalue used for bipolar classification
 */
void recompute_cell_cycles_for_matches_single_cell(VoronoiDiagram &vd,
                                                   int cellIndex,
                                                   float isovalue);

// Full “modify cycles” pass over all global facets at this isovalue.
//! @brief One pass of facet rematching and per‑cell cycle recomputation.
/*!
 * For each global facet, builds iso‑segments, flips matching if necessary, and
 * then recomputes affected cells’ cycles once.
 *
 * @param vd Voronoi diagram to modify
 * @param isovalue Isovalue used for bipolar classification
 */
struct ModifyCyclesResult {
    int interior_flips = 0;
    int boundary_flips = 0;
    int total_flips = 0;

    //! @brief Print modify cycles result for debugging
    template <typename OSTREAM_TYPE>
    void Print(OSTREAM_TYPE & out) const {
        out << "ModifyCyclesResult:\n";
        out << "  Interior flips: " << interior_flips << "\n";
        out << "  Boundary flips: " << boundary_flips << "\n";
        out << "  Total flips: " << total_flips << "\n";
    }
};

ModifyCyclesResult modify_cycles_pass(VoronoiDiagram &vd, float isovalue);

//! @brief Writes a text dump of the Voronoi diagram next to the mesh output.
/*!
 * Produces a human‑readable representation including vertices, edges, cells,
 * facets and connectivity useful for debugging.
 *
 * @param vd Voronoi diagram to serialize
 * @param output_filename Path of the mesh file; determines the dump file name
 */
void write_voronoiDiagram(VoronoiDiagram &vd, std::string &output_filename);
#endif
