#ifndef VDC_GLOBALVAR_H
#define VDC_GLOBALVAR_H

#include "vdc_type.h"
#include "vdc_delaunay.h"
#include "vdc_grid.h"
#include "vdc_voronoi.h"

//Global Variable

extern std::map<Point, float> vertexValueMap;
extern std::vector<Point> activeCubeCenters;
extern std::vector<Object> bipolar_voronoi_edges;
extern std::vector<Point> isosurfaceVertices;

extern K::Iso_cuboid_3 delaunayBbox;

extern std::vector<DelaunayVertex> delaunay_vertices;
extern std::vector<std::pair<Point, bool>> points_with_info;

extern std::map<Point, int> point_index_map;
extern std::map<Point, std::vector<int>> vertex_to_isovertex_indices; // Map Delaunay vertex to its isovertices' indices
extern std::vector<DelaunayTriangle> dualTriangles;
extern std::vector<std::tuple<int, int, int>> isoTriangles;

extern Grid data_grid;
extern Delaunay dt;

extern int isovertex_index; // Global index for isovertexes



#endif