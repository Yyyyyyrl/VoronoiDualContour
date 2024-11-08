#ifndef VDC_DELAUNAY_H
#define VDC_DELAUNAY_H

#include "vdc_type.h"

struct DelaunayTriangle
{
    Point vertex1, vertex2, vertex3;
    DelaunayTriangle(Point v1, Point v2, Point v3) : vertex1(v1), vertex2(v2), vertex3(v3) {}
};

struct IsoTriangle {
    int vertex_indices[3]; // Indices into isosurfaceVertices vector
    IsoTriangle(int idx1, int idx2, int idx3) {
        vertex_indices[0] = idx1;
        vertex_indices[1] = idx2;
        vertex_indices[2] = idx3;
    }
};




#endif