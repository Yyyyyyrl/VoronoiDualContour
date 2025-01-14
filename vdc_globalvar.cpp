#include "vdc_globalvar.h"


std::map<Point, float> vertexValueMap;
std::vector<Point> activeCubeCenters;
std::vector<Object> bipolar_voronoi_edges;
std::vector<Point> isosurfaceVertices;

std::vector<DelaunayVertex> delaunay_vertices;
std::vector<std::pair<Point, bool>> points_with_info;

K::Iso_cuboid_3 delaunayBbox;

std::map<Point, int> point_index_map;
std::map<Point, std::vector<int>> vertex_to_isovertex_indices; // Map Delaunay vertex to its isovertices' indices
std::vector<DelaunayTriangle> dualTriangles;
std::vector<std::tuple<int, int, int>> isoTriangles;


Grid data_grid;
Delaunay dt;

int isovertex_index = 0; // Global index for isovertexes