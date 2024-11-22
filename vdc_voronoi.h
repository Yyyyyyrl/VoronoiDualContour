#ifndef VDC_VORONOI_H
#define VDC_VORONOI_H

#include "vdc_type.h"

struct MidpointNode {
    Point point;                   // Midpoint coordinates
    std::vector<int> connected_to; // Indices of connected midpoints
    int facet_index;               // Index of the facet the midpoint is on
    int cycle_index;               // Index of the cycle the midpoint belongs to
};


struct VoronoiFacet {
    std::vector<int> vertices_indices;    // Ordered indices of vertices in VoronoiDiagram.voronoiVertices
    std::vector<float> vertex_values;     // Scalar values at the vertices
    VoronoiFacet() {}
};

struct Cycle {
    std::vector<int> midpoint_indices;        // Indices of midpoints in the cycle
    std::vector<std::pair<int, int>> edges;   // Edges forming the cycle (indices into midpoints)
    Point isovertex;                          // Centroid of the cycle
    int voronoi_cell_index;                   // Index of the VoronoiCell it belongs to

    void compute_centroid(const std::vector<MidpointNode>& midpoints);
};

struct VoronoiVertex {
    Point vertex;
    std::vector<int> cellIndices; // Indices of Voronoi cells that contain this vertex
    VoronoiVertex(Point p) : vertex(p) {}
};

// VoronoiCell class to represent a Voronoi cell (polytope)
struct VoronoiCell {
    Vertex_handle delaunay_vertex;
    int cellIndex;
    std::vector<int> vertices_indices;       // Indices into VoronoiDiagram.voronoiVertices
    std::vector<int> facet_indices;          // Indices into VoronoiDiagram.voronoiFacets
    CGAL::Polyhedron_3<K> polyhedron;
    std::vector<Cycle> cycles;
    int isoVertexStartIndex;                 // Starting index of isoVertices in VoronoiDiagram.isosurfaceVertices
    int numIsoVertices;                      // Number of isoVertices in this cell
    VoronoiCell(Vertex_handle vh)
        : delaunay_vertex(vh), isoVertexStartIndex(-1), numIsoVertices(0) {}
};

// Construct a new struct store all the Voronoi Diagram together 
// -> Put all voronoi_Vertices in an array and in the Voronoi_Cell
// it should only store the index or reference to the vertices
// -> In facets also the vertices and vertex values should be references

//TODO: Class of isosurface Vertex:
// -> Need Point
// -> reference to the cells that contains it
// -> reference to the cycles that contains it


struct VoronoiDiagram {
    std::vector<VoronoiVertex> voronoiVertices;
    std::vector<Object> voronoiEdges;
    std::vector<float> voronoiVertexValues;
    std::vector<VoronoiCell> voronoiCells;
    std::vector<VoronoiFacet> voronoiFacets;
    std::vector<Point> isosurfaceVertices;
    std::map<Cell_handle, int> cell_to_vertex_index; // Add this line
    std::map<Point, int> point_to_vertex_index;
    std::map<Vertex_handle, int> delaunayVertex_to_voronoiCell_index;
};


// Define a structure to hold the midpoint and associated information
struct EdgeMidpoint {
    Point midpoint;
    int facet_index; // Index of the facet this edge belongs to
};

struct LabeledPoint {
    Point point;
    bool is_dummy;
};

struct PointComparator {
    bool operator()(const Point& a, const Point& b) const {
        const double epsilon = 1e-6;
        return (std::fabs(a.x() - b.x()) < epsilon) &&
               (std::fabs(a.y() - b.y()) < epsilon) &&
               (std::fabs(a.z() - b.z()) < epsilon);
    }
};


template <typename OSTREAM_TYPE>
OSTREAM_TYPE& operator<<(OSTREAM_TYPE& os, const VoronoiFacet& vf)
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

template <typename OSTREAM_TYPE>
OSTREAM_TYPE& operator<<(OSTREAM_TYPE& os, const Cycle& cycle)
{
    os << "Cycle:\n";
    os << "  Voronoi cell index: " << cycle.voronoi_cell_index << "\n";
    os << "  Isovertex: " << cycle.isovertex << "\n";

    os << "  Midpoint indices: ";
    for (const int idx : cycle.midpoint_indices)
        os << idx << " ";
    os << "\n";

    os << "  Edges: ";
    for (const auto& edge : cycle.edges)
        os << "(" << edge.first << ", " << edge.second << ") ";
    os << "\n";

    return os;
}

template <typename OSTREAM_TYPE>
OSTREAM_TYPE& operator<<(OSTREAM_TYPE& os, const VoronoiVertex& vv)
{
    os << "VoronoiVertex:\n";
    os << "  Vertex: " << vv.vertex << "\n";
    os << "  Cell indices: ";
    for (const int idx : vv.cellIndices)
        os << idx << " ";
    os << "\n";

    return os;
}

template <typename OSTREAM_TYPE>
OSTREAM_TYPE& operator<<(OSTREAM_TYPE& os, const VoronoiCell& vc)
{
    os << "VoronoiCell:\n";
    os << "  Cell index: " << vc.cellIndex << "\n";
    os << "  Delaunay vertex: " << vc.delaunay_vertex->point() << "\n";

    os << "  Vertices indices: ";
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
    for (const auto& cycle : vc.cycles)
        os << cycle;

    return os;
}


template <typename OSTREAM_TYPE>
OSTREAM_TYPE& operator<<(OSTREAM_TYPE& os, const VoronoiDiagram& vd)
{
    os << "VoronoiDiagram:\n";

    // Voronoi Vertices
    os << "\nVoronoiVertices:\n";
    for (size_t i = 0; i < vd.voronoiVertices.size(); ++i)
    {
        os << "Index " << i << ":\n";
        os << vd.voronoiVertices[i];
    }

    // Voronoi Edges
    os << "\nVoronoiEdges:\n";
    for (const auto& edge : vd.voronoiEdges)
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

    // Voronoi Vertex Values
    os << "\nVoronoiVertexValues:\n";
    for (size_t i = 0; i < vd.voronoiVertexValues.size(); ++i)
    {
        os << "  Index " << i << ": " << vd.voronoiVertexValues[i] << "\n";
    }

    // Voronoi Cells
    os << "\nVoronoiCells:\n";
    for (const auto& cell : vd.voronoiCells)
    {
        os << cell;
    }

    // Voronoi Facets
    os << "\nVoronoiFacets:\n";
    for (size_t i = 0; i < vd.voronoiFacets.size(); ++i)
    {
        os << "Index " << i << ":\n";
        os << vd.voronoiFacets[i];
    }

    // Isosurface Vertices
    os << "\nIsosurfaceVertices:\n";
    for (size_t i = 0; i < vd.isosurfaceVertices.size(); ++i)
    {
        os << "  Index " << i << ": " << vd.isosurfaceVertices[i] << "\n";
    }

    return os;
}




#endif