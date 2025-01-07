#ifndef VDC_GRID_H
#define VDC_GRID_H

#include "vdc_cube.h"
#include "vdc_type.h"

// DIM = 3 for a 3D grid
static const int DIM3 = 3;

struct Grid
{
    std::vector<float> data;
    int nx, ny, nz;
    float dx, dy, dz;

    void print_grid();
};

struct ScalarGrid
{
    std::vector<std::vector<std::vector<float>>> data; // 3D vector to store scalar values
    int nx, ny, nz;                                    // Dimensions of the grid
    float dx, dy, dz;                                  // Voxel dimensions
    float min_x, min_y, min_z, max_x, max_y, max_z;                         // Minimum coordinates of the grid

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


/**
 * @brief A structure representing a 2D "facet" of a 3D grid,
 *        orthogonal to one axis (orth_dir).
 *
 * For example, if orth_dir = 0 (the x-axis), then this facet stores
 * data along the y- and z-axes.  axis_dir[0] = 1, axis_dir[1] = 2.
 */
struct GRID_FACETS 
{
    int orth_dir;     // Orthogonal direction: 0=x, 1=y, 2=z
    int side;         // Which side: 0=lower(min), 1=upper(max)

    // two axes define this facet (d+1) % 3 and (d+2) % 3
    int axis_dir[2];      

    // Size of this facet along those two axes
    // axis_size[0] = size along axis_dir[0], axis_size[1] = size along axis_dir[1]
    int axis_size[2];    

    //The local bounding-box offset (min/max Index) for each dimension, used for mapping local <-> global
    int minIndex[DIM3];
    int maxIndex[DIM3];
    int localSize[DIM3];

    // Boolean flags for each (coord0, coord1) in the facet
    std::vector<bool> cube_flag;

    /**
     * @param d        Orthogonal direction
     * @param s        Side (0 or 1)
     * @param minIdx   Global bounding-box min for this block of cubes
     * @param maxIdx   Global bounding-box max for this block of cubes
     */
    GRID_FACETS(int d, int s, const int minIdx[DIM3], const int maxIdx[DIM3]);

    // Default destructor is fine since we use std::vector<bool>
    ~GRID_FACETS() = default;

    // Set the flag for a particular (coord0, coord1) in this facet
    void SetFlag(int coord0, int coord1, bool flag);

    // Get the flag for a particular (coord0, coord1) in this facet
    bool CubeFlag(int coord0, int coord1) const;

private:
    // Convert (coord0, coord1) to linear index
    int index(int coord0, int coord1) const;
};



// Functions for loading nrrd data
template <typename T>
std::vector<float> convert_to_float_vector(T *data_ptr, size_t total_size);

Point adjust_outside_bound_points(const Point &p, const ScalarGrid &grid, const Point &v1, const Point &v2);

/*Preprocessing*/
void initialize_scalar_grid(ScalarGrid &grid, const Grid &nrrdGrid);
Grid load_nrrd_data(const std::string &file_path);
Grid supersample_grid(const Grid &grid, int n);

bool is_cube_active(const Grid &grid, int x, int y, int z, float isovalue);
void find_active_cubes(const Grid &grid, float isovalue, std::vector<Cube> &cubes);

std::vector<Point> load_grid_points(const Grid &grid);

bool is_point_inside_grid(const Point &p, const ScalarGrid &grid);

Point interpolate(const Point &p1, const Point &p2, float val1, float val2, float isovalue, const Grid &data_grid);
float trilinear_interpolate(const Point &p, const ScalarGrid &grid);
float trilinear_interpolate(const Point &p, const Grid &grid);

std::vector<std::vector<GRID_FACETS>> create_grid_facets(const std::vector<Cube> &activeCubes);
#endif