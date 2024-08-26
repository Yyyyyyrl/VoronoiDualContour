#ifndef DMR_H
#define DMR_H

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>
#include <limits>
#include <algorithm>
#include <teem/nrrd.h>
#include <CGAL/bounding_box.h>
#include <CGAL/intersections.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_with_circumcenter_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K> Vb;
typedef CGAL::Delaunay_triangulation_cell_base_with_circumcenter_3<K> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds> Delaunay;
typedef Delaunay::Point Point;
typedef Delaunay::Cell_handle Cell_handle;
typedef Delaunay::Vertex_handle Vertex_handle;
typedef Delaunay::Cell_circulator Cell_circulator;
typedef CGAL::Object Object;
typedef Delaunay::Edge Edge;
typedef Delaunay::Facet Facet;
typedef Delaunay::Triangle Triangle;
typedef K::Segment_3 Segment3;
typedef K::Ray_3 Ray3;
typedef K::Line_3 Line3;
typedef K::Point_3 Point3;
typedef K::Vector_3 Vector3;

// Structures

struct Grid
{
    std::vector<float> data;
    int nx, ny, nz;
    float dx, dy, dz;
};

struct DelaunayTriangle
{
    Point vertex1, vertex2, vertex3;
    DelaunayTriangle(Point v1, Point v2, Point v3) : vertex1(v1), vertex2(v2), vertex3(v3) {}
};

struct ScalarGrid
{
    std::vector<std::vector<std::vector<float>>> data; // 3D vector to store scalar values
    int nx, ny, nz;                                    // Dimensions of the grid
    float dx, dy, dz;                                  // Voxel dimensions
    float min_x, min_y, min_z;                         // Minimum coordinates of the grid

    ScalarGrid(int nx, int ny, int nz, float dx, float dy, float dz, float min_x, float min_y, float min_z);

    // Method to get a value from the grid
    float get_value(int x, int y, int z) const;

    // Method to set a value in the grid
    void set_value(int x, int y, int z, float value);

    // Method to load grid data from an external source
    void load_from_source(const std::vector<std::vector<std::vector<float>>> &source);
};




// Function prototypes
void read_data(const std::string& filename, std::vector<Point>& points);
void compute_triangulation(const std::vector<Point>& points, Delaunay& triangulation);
void save_results(const Delaunay& triangulation, const std::string& filename);
void print_cell(Delaunay::Cell c);
void initialize_scalar_grid(ScalarGrid &grid, const Grid &nrrdGrid);
void write_dt_to_off(const Delaunay &dt, const std::string &filename);
void print_facet(Facet f);
void printActiveCubeCenters(std::vector<Point> activeCubeCenters);



Grid load_nrrd_data(const std::string &file_path);


bool is_cube_active(const Grid &grid, int x, int y, int z, float isovalue);
bool is_bipolar(float val1, float val2, float isovalue = 0);
bool isDegenerate(const Object &obj);


std::vector<Point> find_active_cubes(const Grid &grid, float isovalue);
std::vector<Point> load_grid_points(const Grid &grid);
std::vector<Point> process_active_cube(const std::array<Point, 8> &cube_corners, const std::array<float, 8> &scalar_values, float isovalue);

Point interpolate(const Point &p1, const Point &p2, float val1, float val2, float isovalue);
Point compute_centroid(const std::vector<Point> &points);

float trilinear_interpolate(const Point &p, const ScalarGrid &grid);
float get_scalar_value_at_point(const Point &point, const ScalarGrid &grid);

std::array<Point, 8> get_cube_corners(const Point &center, float side_length);
std::tuple<int, int, int> point_to_grid_index(const Point &point, const ScalarGrid &grid);

K::Vector_3 calculate_normal(const Point &p0, const Point &p1, const Point &p2);

void writeOFF(const std::string &filename, const std::vector<Point> &vertices, const std::vector<DelaunayTriangle> &triangles, std::map<Point, int> &pointIndexMap);
void writePLY(const std::string &filename, const std::vector<Point> &vertices, const std::vector<DelaunayTriangle> &triangles, std::map<Point, int> &pointIndexMap);

int get_orientation(const int iFacet, const Point v1, const Point v2, const float f1, const float f2);


std::string objectToString(const CGAL::Object &obj);


#endif // DMR_H
