#ifndef VDC_VORONOI_H
#define VDC_VORONOI_H

#include "vdc_type.h"

struct VoronoiFacet {
    std::vector<Point> vertices;          // Ordered vertices of the facet
    std::vector<float> vertex_values;     // Scalar values at the vertices

    VoronoiFacet(const std::vector<Point>& verts, const std::vector<float>& values)
        : vertices(verts), vertex_values(values) {}
};

struct Cycle {
    std::vector<Point> midpoints; // Midpoints forming the cycle \remove
    std::vector<std::pair<int, int>> edges; // Edges that forms the cycle where each edge is represented using indices of the two points
    Point isovertex;              // Centroid of the cycle
    int voronoi_cell_index;       // Index of the VoronoiCell it belongs to
    
    // Stores edges which mid points forms the cycle and reference to them
    void compute_centroid();

};

struct VoronoiVertex {
    Point vertex;
    std::vector<int> cellIndices; // Indices of voronoi cells that contains this vertex
    std::vector<int> cycleIndices; // Indices of cycles that contains this vertex

    VoronoiVertex(Point p) : vertex(p) {} 
};

// VoronoiCell class to represent a Voronoi cell (polytope)
struct VoronoiCell {
    Vertex_handle delaunay_vertex;                    // Corresponding Delaunay vertex
    std::vector<Point> voronoi_vertices;              // Voronoi vertices (coordinates)
    std::vector<float> vertex_values;                 // Scalar values at the vertices
    std::vector<VoronoiFacet> facets;                        // Facets of the Voronoi cell
    CGAL::Polyhedron_3<K> polyhedron;                 // Polyhedron representing the cell
    std::vector<Cycle> cycles;
    std::map<VoronoiFacet, int> facet_to_cycle_map; // Map from facet to cycle index

    VoronoiCell(Vertex_handle vh) : delaunay_vertex(vh) {}
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
    std::vector<Point> voronoiVertices;
    std::vector<Object> voronoiEdges;
    std::map<Point, float> voronoiVertex2ValuesMap;
    std::vector<Point> midpoints;
    std::map<VoronoiCell, int> voronoiCells2IndexMap;
    std::map<VoronoiFacet, int> voronoiFacets2IndexMap;  
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

struct MidpointNode {
    Point point;                   // Midpoint coordinates
    std::vector<int> connected_to; // Indices of connected midpoints
    int facet_index;               // Index of the facet the midpoint is on
    int cycle_index;               // Index of the cycle the midpoint belongs to
};

#endif