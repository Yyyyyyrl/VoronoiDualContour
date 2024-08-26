#ifndef UTILITIES_H
#define UTILITIES_H

#include "io.h"
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
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K> Vb;
typedef CGAL::Object Object;
typedef K::Segment_3 Segment3;
typedef K::Ray_3 Ray3;
typedef K::Line_3 Line3;
typedef K::Point_3 Point3;
typedef K::Vector_3 Vector3;

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

struct ObjectComparator
{
    bool operator()(const Object &obj1, const Object &obj2) const
    {
        return objectToString(obj1) < objectToString(obj2);
    }
};

/*Preprocessing*/
void initialize_scalar_grid(ScalarGrid &grid, const Grid &nrrdGrid);

/*
Functions for Active Cube Centers
*/
bool is_cube_active(const Grid &grid, int x, int y, int z, float isovalue);

std::vector<Point> find_active_cubes(const Grid &grid, float isovalue);

std::vector<Point> load_grid_points(const Grid &grid);
bool is_bipolar(float val1, float val2, float isovalue = 0);
bool isDegenerate(const Object &obj);

std::string objectToString(const CGAL::Object &obj);

/*
General Helper Functions
*/
Point compute_centroid(const std::vector<Point> &points);
Point interpolate(const Point &p1, const Point &p2, float val1, float val2, float isovalue);
float trilinear_interpolate(const Point &p, const ScalarGrid &grid);
std::array<Point, 8> get_cube_corners(const Point &center, float side_length);
int get_orientation(const int iFacet, const Point v1, const Point v2, const float f1, const float f2);


#endif UTILITIES_H