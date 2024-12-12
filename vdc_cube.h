#ifndef VDC_CUBE_H
#define VDC_CUBE_H

#include "vdc_type.h"

struct Cube {
    Point repVertex, center;
    int sidelength;
    int i, j, k; // store indices of the cube
    Cube() : repVertex(0,0,0), center(0,0,0), sidelength(1), i(0), j(0), k(0) {}
    Cube(Point v, Point c, int len, int ix, int iy, int iz) : repVertex(v), center(c), sidelength(len), i(ix), j(iy), k(iz) {}
};


bool is_adjacent(const Cube &cubeA, const Cube &cubeB);
int get_cube_index(const Point &repVertex, int nx, int ny);
std::vector<int> find_neighbor_indices(const Point3& repVertex, int nx, int ny);
std::vector<Point> get_cube_centers(const std::vector<Cube> &cubes);
std::vector<Cube> separate_active_cubes_greedy(std::vector<Cube> &activeCubes, int nx, int ny, int nz);
std::vector<Cube> separate_active_cubes_graph(std::vector<Cube> &activeCubes);
#endif


//TODO: 