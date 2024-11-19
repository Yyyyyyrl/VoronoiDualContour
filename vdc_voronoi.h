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
#endif