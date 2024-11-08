#ifndef VDC_CUBE_H
#define VDC_CUBE_H

#include "vdc_type.h"

struct Cube
{
    // The vertex stands for the lowest corner vertex of the cube, for example the cube defined by (0,0,0) and (1,1,1) will have repVertex (0,0,0)
    Point repVertex, center;
    int sidelength;
    // Default constructor
    Cube() : repVertex(Point3(0, 0, 0)), center(Point3(0, 0, 0)), sidelength(1) {}
    Cube(Point v, Point c, int len) : repVertex(v), sidelength(len), center(c) {}
};

bool is_adjacent(const Cube &cubeA, const Cube &cubeB);
int get_cube_index(const Point &repVertex, int nx, int ny);
std::vector<int> find_neighbor_indices(const Point3& repVertex, int nx, int ny);
std::vector<Point> get_cube_centers(const std::vector<Cube> &cubes);
std::vector<Cube> separate_active_cubes_greedy(std::vector<Cube> &activeCubes, int nx, int ny, int nz);
std::vector<Cube> separate_active_cubes_graph(std::vector<Cube> &activeCubes);
#endif