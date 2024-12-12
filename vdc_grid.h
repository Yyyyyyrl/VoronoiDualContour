#ifndef VDC_GRID_H
#define VDC_GRID_H

#include "vdc_cube.h"
#include "vdc_type.h"

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

#include <vector>

struct GRID_FACETS {
    int orth_dir;      // Which axis: 0=x, 1=y, 2=z
    int side;          // Which side: 0=lower(min), 1=upper(max)
    int axis_size[3];  // Number of cubes along each dimension: {Nx, Ny, Nz}
    std::vector<bool> cube_flag; // Use vector<bool> to store flags

    // Constructor
    GRID_FACETS(int OrthDir, int Side, int Nx, int Ny, int Nz);

    // Default destructor
    ~GRID_FACETS() = default; // no need for custom destructor

    // Set flag for a particular (x,y) on the facet slice
    void SetFlag(int x, int y, bool flag);

    // Get flag for a particular (x,y) on the facet slice
    bool CubeFlag(int x, int y) const;

private:
    int index(int x, int y) const;
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
std::vector<Cube> find_active_cubes(const Grid &grid, float isovalue);

std::vector<Point> load_grid_points(const Grid &grid);

bool is_point_inside_grid(const Point &p, const ScalarGrid &grid);

Point interpolate(const Point &p1, const Point &p2, float val1, float val2, float isovalue, const Grid &data_grid);
float trilinear_interpolate(const Point &p, const ScalarGrid &grid);
float trilinear_interpolate(const Point &p, const Grid &grid);

#endif