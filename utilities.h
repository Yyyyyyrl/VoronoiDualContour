#ifndef UTILITIES_H
#define UTILITIES_H

#include "debug.h"

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

// Define the kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

// Define the vertex base with information
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K> Vb;
typedef CGAL::Delaunay_triangulation_cell_base_with_circumcenter_3<K> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;

// Define Delaunay Triangulation
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


/*
Structs
*/

struct Cube
{
    // The vertex stands for the lowest corner vertex of the cube, for example the cube defined by (0,0,0) and (1,1,1) will have repVertex (0,0,0)
    Point repVertex, center;
    int sidelength;
    // Default constructor
    Cube() : repVertex(Point3(0, 0, 0)), center(Point3(0, 0, 0)), sidelength(1) {}
    Cube(Point v, Point c, int len) : repVertex(v), sidelength(len), center(c) {}
};


struct Grid
{
    std::vector<float> data;
    int nx, ny, nz;
    float dx, dy, dz;
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

    float get_scalar_value_at_point(const Point &point);

    std::tuple<int, int, int> point_to_grid_index(const Point &point);


};

struct DelaunayTriangle
{
    Point vertex1, vertex2, vertex3;
    DelaunayTriangle(Point v1, Point v2, Point v3) : vertex1(v1), vertex2(v2), vertex3(v3) {}
};

std::string objectToString(const CGAL::Object &obj);

struct ObjectComparator
{
    bool operator()(const Object &obj1, const Object &obj2) const
    {
        return objectToString(obj1) < objectToString(obj2);
    }
};

// Functions for loading nrrd data
template <typename T>
std::vector<float> convert_to_float_vector(T *data_ptr, size_t total_size);


/*Preprocessing*/
void initialize_scalar_grid(ScalarGrid &grid, const Grid &nrrdGrid);
Grid load_nrrd_data(const std::string &file_path);
Grid supersample_grid(const Grid &grid, int n);


/*
Functions for Active Cube Centers
*/
bool is_cube_active(const Grid &grid, int x, int y, int z, float isovalue);
bool is_adjacent(const Cube &cubeA, const Cube &cubeB);
int get_cube_index(const Point &repVertex, int nx, int ny);
std::vector<int> find_neighbor_indices(const Point3& repVertex, int nx, int ny);
std::vector<Point> get_cube_centers(const std::vector<Cube> &cubes);
std::vector<Cube> separate_active_cubes_greedy(std::vector<Cube> &activeCubes, int nx, int ny, int nz);
std::vector<Cube> separate_active_cubes_graph(std::vector<Cube> &activeCubes);

std::vector<Cube> find_active_cubes(const Grid &grid, float isovalue);

std::vector<Point> load_grid_points(const Grid &grid);
bool is_bipolar(float val1, float val2, float isovalue = 0);
bool isDegenerate(const Object &obj);


/*
General Helper Functions
*/
Point compute_centroid(const std::vector<Point> &points);
Point interpolate(const Point &p1, const Point &p2, float val1, float val2, float isovalue, const Grid &data_grid);
float trilinear_interpolate(const Point &p, const ScalarGrid &grid);
float trilinear_interpolate(const Point &p, const Grid &grid);
std::array<Point, 8> get_cube_corners(const Point &center, float side_length);
int get_orientation(const int iFacet, const Point v1, const Point v2, const float f1, const float f2);

#endif