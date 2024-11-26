#include "vdc_globalvar.h"

std::string file_path;
float isovalue;
std::string output_format;
std::string output_filename;
std::string out_csv_name;
bool out_csv = false;
bool sep_isov = false;
bool supersample = false;
bool multi_isov = true;
int supersample_r;

std::map<Point, float> vertexValueMap;
std::vector<Point> activeCubeCenters;
std::vector<Object> bipolar_voronoi_edges;
std::vector<Point> isosurfaceVertices;

std::vector<LabeledPoint> all_points;
std::vector<std::pair<Point, bool>> points_with_info;

K::Iso_cuboid_3 delaunayBbox;

std::map<Point, int> point_index_map;
std::map<Point, std::vector<int>> vertex_to_isovertex_indices; // Map Delaunay vertex to its isovertices' indices
std::vector<DelaunayTriangle> dualTriangles;
std::vector<std::tuple<int, int, int>> isoTriangles;


Grid data_grid;
Delaunay dt;

int isovertex_index = 0; // Global index for isovertexes
