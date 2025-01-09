#include "vdc_grid.h"


//! Constructor for the ScalarGrid.
ScalarGrid::ScalarGrid(int nx, int ny, int nz, float dx, float dy, float dz, float min_x, float min_y, float min_z)
    : nx(nx), ny(ny), nz(nz), dx(dx), dy(dy), dz(dz), min_x(min_x), min_y(min_y), min_z(min_z)
{
    //! Initialize the 3D data vector with default values of 0.
    data.resize(nx, std::vector<std::vector<float>>(ny, std::vector<float>(nz, 0.0)));
}

//! Retrieve a scalar value from the grid.
float ScalarGrid::get_value(int x, int y, int z) const
{
    if (x < 0 || x >= nx || y < 0 || y >= ny || z < 0 || z >= nz)
    {
        //! Return 0 if indices are out of bounds.
        return 0;
    }
    return data[x][y][z];
}

//! Set a scalar value in the grid.
void ScalarGrid::set_value(int x, int y, int z, float value)
{
    if (x >= 0 && x < nx && y >= 0 && y < ny && z >= 0 && z < nz)
    {
        data[x][y][z] = value;
    }
}

//! Load scalar values into the grid from an external 3D vector source.
void ScalarGrid::load_from_source(const std::vector<std::vector<std::vector<float>>> &source)
{
    for (int i = 0; i < nx && i < source.size(); ++i)
    {
        for (int j = 0; j < ny && j < source[i].size(); ++j)
        {
            for (int k = 0; k < nz && k < source[i][j].size(); ++k)
            {
                data[i][j][k] = source[i][j][k];
            }
        }
    }
}

//! @brief Constructor for the GRID_FACETS structure.
GRID_FACETS::GRID_FACETS(int d, int s, const int minIdx[DIM3], const int maxIdx[DIM3])
    : orth_dir(d), side(s)
{
    //! Store `minIndex[]` and `maxIndex[]`, and compute `localSize[]`.
    for (int i = 0; i < DIM3; i++) {
        minIndex[i] = minIdx[i];
        maxIndex[i] = maxIdx[i];
        localSize[i] = (maxIndex[i] - minIndex[i] + 1);
    }

    //! Determine the 2D slice axes.
    axis_dir[0] = (orth_dir + 1) % DIM3;
    axis_dir[1] = (orth_dir + 2) % DIM3;

    //! Compute the size of the facet along the two axes.
    axis_size[0] = localSize[axis_dir[0]];
    axis_size[1] = localSize[axis_dir[1]];

    //! Allocate flags for the facet, initialized to `false`.
    cube_flag.resize(axis_size[0] * axis_size[1], false);
}

//! @brief Set the flag for a particular `(coord0, coord1)` in the facet.
void GRID_FACETS::SetFlag(int coord0, int coord1, bool flag)
{
    cube_flag[index(coord0, coord1)] = flag;
}

//! @brief Get the flag for a particular `(coord0, coord1)` in the facet.
bool GRID_FACETS::CubeFlag(int coord0, int coord1) const
{
    return cube_flag[index(coord0, coord1)];
}

//! @brief Convert a 2D coordinate `(coord0, coord1)` to a linear index.
int GRID_FACETS::index(int coord0, int coord1) const
{
    return coord1 * axis_size[0] + coord0;
}

//! @brief Print the grid's metadata and data in a human-readable format.
void Grid::print_grid()
{
    //! Print metadata (dimensions and spacing).
    std::cout << "NRRD Grid Information:\n";
    std::cout << "Dimensions: " << nx << " " << ny << " " << nz << "\n";
    std::cout << "Spacing: " << dx << " " << dy << " " << dz << "\n\n";

    //! Print the grid's data in a structured format.
    std::cout << "Data:\n";

    for (int z = 0; z < nz; ++z)
    {
        std::cout << "Slice z = " << z << ":\n";
        for (int y = 0; y < ny; ++y)
        {
            for (int x = 0; x < nx; ++x)
            {
                //! Compute the linear index for the current x, y, z coordinate.
                std::cout << std::setw(8) << data[z * nx * ny + y * nx + x] << " ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }
}

std::tuple<int, int, int> ScalarGrid::point_to_grid_index(const Point &point)
{
    int x = static_cast<int>((point.x() - min_x) / dx);
    int y = static_cast<int>((point.y() - min_y) / dy);
    int z = static_cast<int>((point.z() - min_z) / dz);
    return {x, y, z};
}

float ScalarGrid::get_scalar_value_at_point(const Point &point)
{
    auto [x, y, z] = point_to_grid_index(point);
    return get_value(x, y, z);
}


//! @brief Converts raw data of type `T` into a float vector.
/*!
 * @tparam T The type of the input data.
 * @param data_ptr Pointer to the input data.
 * @param total_size Total number of elements in the input data.
 * @return A `std::vector<float>` containing the converted data.
 */
template <typename T>
std::vector<float> convert_to_float_vector(T *data_ptr, size_t total_size)
{
    std::vector<float> data(total_size);
    for (size_t i = 0; i < total_size; ++i)
    {
        data[i] = static_cast<float>(data_ptr[i]);
    }
    return data;
}

//! @brief Loads NRRD data into a `Grid` structure.
Grid load_nrrd_data(const std::string &file_path)
{
    Nrrd *nrrd = nrrdNew();
    if (nrrdLoad(nrrd, file_path.c_str(), NULL))
    {
        char *err = biffGetDone(NRRD);
        std::cerr << "Error reading NRRD file: " << err << std::endl;
        free(err);
        nrrdNuke(nrrd);
        exit(1);
    }

    size_t total_size = nrrdElementNumber(nrrd);

    std::vector<float> data;

    if (nrrd->type == nrrdTypeFloat)
    {
        float *data_ptr = static_cast<float *>(nrrd->data);
        data = std::vector<float>(data_ptr, data_ptr + total_size);
    }
    else if (nrrd->type == nrrdTypeUChar)
    {
        unsigned char *data_ptr = static_cast<unsigned char *>(nrrd->data);
        data = convert_to_float_vector(data_ptr, total_size);
    }
    else
    {
        std::cerr << "Unsupported NRRD data type." << std::endl;
        nrrdNuke(nrrd);
        exit(1);
    }

    int nx = nrrd->axis[0].size;
    int ny = nrrd->axis[1].size;
    int nz = nrrd->axis[2].size;
    float dx = nrrd->axis[0].spacing;
    float dy = nrrd->axis[1].spacing;
    float dz = nrrd->axis[2].spacing;

    // Initialize min/max values
    float min_x = 0.0f, max_x = (nx - 1) * dx;
    float min_y = 0.0f, max_y = (ny - 1) * dy;
    float min_z = 0.0f, max_z = (nz - 1) * dz;

    // You could also check the bounds dynamically by looking at the index
    // but in this case, it's assumed the data aligns perfectly with grid structure

    std::cout << "Grid dimensions: " << nx << "x" << ny << "x" << nz << std::endl;
    std::cout << "Spacing: dx=" << dx << ", dy=" << dy << ", dz=" << dz << std::endl;
    std::cout << "X range: [" << min_x << ", " << max_x << "]" << std::endl;
    std::cout << "Y range: [" << min_y << ", " << max_y << "]" << std::endl;
    std::cout << "Z range: [" << min_z << ", " << max_z << "]" << std::endl;

    nrrdNuke(nrrd); // Properly dispose of the Nrrd structure

    // TODO: print_grid()
    return {data, nx, ny, nz, dx, dy, dz};
}

//! @brief Adjusts a point outside the bounds of the scalar grid.
Grid supersample_grid(const Grid &grid, int n)
{
    int nx2 = grid.nx * n - (n - 1);
    int ny2 = grid.ny * n - (n - 1);
    int nz2 = grid.nz * n - (n - 1);

    float dx2 = grid.dx / n;
    float dy2 = grid.dy / n;
    float dz2 = grid.dz / n;

    std::vector<float> data2(nx2 * ny2 * nz2);

    for (int z = 0; z < nz2; ++z)
    {
        for (int y = 0; y < ny2; ++y)
        {
            for (int x = 0; x < nx2; ++x)
            {
                // Convert to original grid space (ensure float division)
                float px = static_cast<float>(x) / n * grid.dx;
                float py = static_cast<float>(y) / n * grid.dy;
                float pz = static_cast<float>(z) / n * grid.dz;

                // Perform trilinear interpolation
                float interpolate_val = trilinear_interpolate(Point(px, py, pz), grid);
                data2[z * nx2 * ny2 + y * nx2 + x] = interpolate_val;
            }
        }
    }

    return {data2, nx2, ny2, nz2, dx2, dy2, dz2};
}

//! @brief Inidialize the ScalarGrid instance with a Grid instance consisting the data read from a nrrd file
void initialize_scalar_grid(ScalarGrid &grid, const Grid &nrrdGrid)
{
    // Use the dimensions from the loaded nrrdGrid
    grid.nx = nrrdGrid.nx;
    grid.ny = nrrdGrid.ny;
    grid.nz = nrrdGrid.nz;

    // Define grid dimensions
    grid.min_x = 0;
    grid.min_y = 0;
    grid.min_z = 0;

    // Define grid spacing (dx, dy, dz)
    grid.dx = nrrdGrid.dx;
    grid.dy = nrrdGrid.dy;
    grid.dz = nrrdGrid.dz;

    grid.max_x = nrrdGrid.dx * (grid.nx - 1);
    grid.max_y = nrrdGrid.dy * (grid.ny - 1);
    grid.max_z = nrrdGrid.dz * (grid.nz - 1);

    // Resizing and initializing the scalar grid data array
    grid.data.resize(grid.nx);
    for (int i = 0; i < grid.nx; ++i)
    {
        grid.data[i].resize(grid.ny);
        for (int j = 0; j < grid.ny; ++j)
        {
            grid.data[i][j].resize(grid.nz, 0.0);
        }
    }
    // Iterate through each voxel in the grid to initialize values from the nrrdGrid data
    for (int i = 0; i < grid.nx; i++)
    {
        for (int j = 0; j < grid.ny; j++)
        {
            for (int k = 0; k < grid.nz; k++)
            {
                int index = i + j * grid.nx + k * grid.nx * grid.ny;
                if (index < nrrdGrid.data.size())
                {
                    grid.data[i][j][k] = nrrdGrid.data[index]; // Assigning the value from nrrdGrid
                }
            }
        }
    }
}

//! @brief Checks if a cube is active based on scalar values.
bool is_cube_active(const Grid &grid, int x, int y, int z, float isovalue)
{
    // Offsets for each of the 8 cube vertices
    static const std::vector<std::tuple<int,int,int>> vertex_offsets = {
        {0, 0, 0},  // Vertex 0
        {1, 0, 0},  // Vertex 1
        {1, 1, 0},  // Vertex 2
        {0, 1, 0},  // Vertex 3
        {0, 0, 1},  // Vertex 4
        {1, 0, 1},  // Vertex 5
        {1, 1, 1},  // Vertex 6
        {0, 1, 1}   // Vertex 7
    };

    // Helper lambda to quickly retrieve the scalar value at a given offset
    auto get_value = [&](int dx, int dy, int dz) -> float {
        // We assume x+dx, y+dy, z+dz are valid indices here.
        int idx = (x + dx) + (y + dy) * grid.nx + (z + dz) * grid.nx * grid.ny;
        return grid.data[idx];
    };

    // Check the sign (negative vs. positive) of the first vertex
    float val0 = get_value(0, 0, 0);
    bool is_val0_negative = (val0 < isovalue);

    // Compare remaining vertices against the first vertex's sign
    for (int i = 1; i < 8; i++) {
        auto [dx, dy, dz] = vertex_offsets[i];
        float vali = get_value(dx, dy, dz);
        bool is_vali_negative = (vali < isovalue);

        // If the sign of vali differs from val0, the cube is active.
        if (is_vali_negative != is_val0_negative) {
            return true;
        }
    }

    // If all vertices have the same sign (all negative or all positive),
    // then the cube is not active.
    return false;
}

//! @brief Adjusts a point outside the bounds of the scalar grid.
Point adjust_outside_bound_points(const Point &p, const ScalarGrid &grid, const Point &v1, const Point &v2)
{
    // Convert the point to grid space
    float gx = p.x() / grid.dx;
    float gy = p.y() / grid.dy;
    float gz = p.z() / grid.dz;

    // Check if the point is outside the bounds
    if (gx < 0 || gx >= grid.nx || gy < 0 || gy >= grid.ny || gz < 0 || gz >= grid.nz)
    {
        // If the point is out of bounds, find the closest in-bound point on the segment (v1, v2)
        // Convert v1 and v2 to grid space
        float v1_gx = v1.x() / grid.dx;
        float v1_gy = v1.y() / grid.dy;
        float v1_gz = v1.z() / grid.dz;

        float v2_gx = v2.x() / grid.dx;
        float v2_gy = v2.y() / grid.dy;
        float v2_gz = v2.z() / grid.dz;

        // Get the parameter t for the projection of point p onto the segment [v1, v2]
        float t = ((gx - v1_gx) * (v2_gx - v1_gx) + (gy - v1_gy) * (v2_gy - v1_gy) + (gz - v1_gz) * (v2_gz - v1_gz)) /
                  ((v2_gx - v1_gx) * (v2_gx - v1_gx) + (v2_gy - v1_gy) * (v2_gy - v1_gy) + (v2_gz - v1_gz) * (v2_gz - v1_gz));

        // Clamp t to [0, 1] to ensure the closest point lies on the segment
        t = std::max(0.0f, std::min(t, 1.0f));

        // Compute the closest point in grid space
        float px = v1_gx + t * (v2_gx - v1_gx);
        float py = v1_gy + t * (v2_gy - v1_gy);
        float pz = v1_gz + t * (v2_gz - v1_gz);

        // Convert the grid-space point back to world-space and return
        return Point(px * grid.dx, py * grid.dy, pz * grid.dz);
    }

    // If the point is within bounds, return the original point
    return p;
}


void find_active_cubes(const Grid &grid, float isovalue, std::vector<Cube> &cubes) {
    cubes.clear();
    // (nx-1, ny-1, nz-1) cubes total
    for (int i = 0; i < grid.nx - 1; ++i)
    {
        for (int j = 0; j < grid.ny - 1; ++j)
        {
            for (int k = 0; k < grid.nz - 1; ++k)
            {
                if (is_cube_active(grid, i, j, k, isovalue))
                {
                    Point repVertex(i * grid.dx, j * grid.dy, k * grid.dz);
                    Point center((i+0.5f)*grid.dx, (j+0.5f)*grid.dy, (k+0.5f)*grid.dz);
                    cubes.push_back(Cube(repVertex, center, 1, i, j, k));
                }
            }
        }
    }
}


std::vector<Point> load_grid_points(const Grid &grid)
{
    std::vector<Point> points;
    for (int i = 0; i < grid.nx; ++i)
    {
        for (int j = 0; j < grid.ny; ++j)
        {
            for (int k = 0; k < grid.nz; ++k)
            {
                points.push_back(Point(i, j, k));
            }
        }
    }
    return points;
}

bool is_point_inside_grid(const Point &p, const ScalarGrid &grid)
{
    return (p.x() >= grid.min_x && p.x() <= grid.max_x &&
            p.y() >= grid.min_y && p.y() <= grid.max_y &&
            p.z() >= grid.min_z && p.z() <= grid.max_z);
}

Point interpolate(const Point &p1, const Point &p2, float val1, float val2, float isovalue, const Grid &data_grid)
{
    if (std::abs(val1 - val2) < 1e-6) // Avoid division by zero or near-zero differences
        return p1;
    float t = (isovalue - val1) / (val2 - val1);
    return Point(p1.x() + t * (p2.x() - p1.x()) * data_grid.dx,
                 p1.y() + t * (p2.y() - p1.y()) * data_grid.dy,
                 p1.z() + t * (p2.z() - p1.z()) * data_grid.dz);
}

float trilinear_interpolate(const Point &p, const ScalarGrid &grid)
{
    bool debug = false;
    float gx = p.x() / grid.dx;
    float gy = p.y() / grid.dy;
    float gz = p.z() / grid.dz;

    // Clamp gx, gy, gz to valid grid bounds
    gx = std::max(0.0f, std::min(gx, (float)(grid.nx - 1)));
    gy = std::max(0.0f, std::min(gy, (float)(grid.ny - 1)));
    gz = std::max(0.0f, std::min(gz, (float)(grid.nz - 1)));

    int x0 = (int)std::floor(gx);
    int x1 = std::min(x0 + 1, grid.nx - 1);
    int y0 = (int)std::floor(gy);
    int y1 = std::min(y0 + 1, grid.ny - 1);
    int z0 = (int)std::floor(gz);
    int z1 = std::min(z0 + 1, grid.nz - 1);

    // Fractional part of coordinates
    float xd = gx - x0;
    float yd = gy - y0;
    float zd = gz - z0;

    // Get the values from the grid at the 8 corner points
    float c000 = grid.get_value(x0, y0, z0);
    float c001 = grid.get_value(x0, y0, z1);
    float c010 = grid.get_value(x0, y1, z0);
    float c011 = grid.get_value(x0, y1, z1);
    float c100 = grid.get_value(x1, y0, z0);
    float c101 = grid.get_value(x1, y0, z1);
    float c110 = grid.get_value(x1, y1, z0);
    float c111 = grid.get_value(x1, y1, z1);

    // Perform trilinear interpolation
    float c00 = c000 * (1 - zd) + c001 * zd;
    float c01 = c010 * (1 - zd) + c011 * zd;
    float c10 = c100 * (1 - zd) + c101 * zd;
    float c11 = c110 * (1 - zd) + c111 * zd;

    float c0 = c00 * (1 - yd) + c01 * yd;
    float c1 = c10 * (1 - yd) + c11 * yd;

    float c = c0 * (1 - xd) + c1 * xd;

    return c;
}

float trilinear_interpolate(const Point &p, const Grid &grid)
{
    // Convert point coordinates to grid space
    float gx = p.x() / grid.dx;
    float gy = p.y() / grid.dy;
    float gz = p.z() / grid.dz;

    // Determine the indices of the eight surrounding grid points
    int x0 = static_cast<int>(std::floor(gx));
    int x1 = std::min(x0 + 1, grid.nx - 1); // Clamp x1 to grid boundary
    int y0 = static_cast<int>(std::floor(gy));
    int y1 = std::min(y0 + 1, grid.ny - 1); // Clamp y1 to grid boundary
    int z0 = static_cast<int>(std::floor(gz));
    int z1 = std::min(z0 + 1, grid.nz - 1); // Clamp z1 to grid boundary

    // Compute the differences
    float xd = gx - x0;
    float yd = gy - y0;
    float zd = gz - z0;

    // Retrieve values at the eight surrounding grid points
    float c000 = grid.data[z0 * grid.nx * grid.ny + y0 * grid.nx + x0];
    float c001 = grid.data[z1 * grid.nx * grid.ny + y0 * grid.nx + x0];
    float c010 = grid.data[z0 * grid.nx * grid.ny + y1 * grid.nx + x0];
    float c011 = grid.data[z1 * grid.nx * grid.ny + y1 * grid.nx + x0];
    float c100 = grid.data[z0 * grid.nx * grid.ny + y0 * grid.nx + x1];
    float c101 = grid.data[z1 * grid.nx * grid.ny + y0 * grid.nx + x1];
    float c110 = grid.data[z0 * grid.nx * grid.ny + y1 * grid.nx + x1];
    float c111 = grid.data[z1 * grid.nx * grid.ny + y1 * grid.nx + x1];

    // Interpolate along z-axis
    float c00 = c000 * (1 - zd) + c001 * zd;
    float c01 = c010 * (1 - zd) + c011 * zd;
    float c10 = c100 * (1 - zd) + c101 * zd;
    float c11 = c110 * (1 - zd) + c111 * zd;

    // Interpolate along y-axis
    float c0 = c00 * (1 - yd) + c01 * yd;
    float c1 = c10 * (1 - yd) + c11 * yd;

    // Interpolate along x-axis
    float c = c0 * (1 - xd) + c1 * xd;

    return c;
}

//! @brief Creates grid facets for active cubes.
std::vector<std::vector<GRID_FACETS>> create_grid_facets(const std::vector<Cube> &activeCubes) {

    int minIdx[3];
    int maxIdx[3];

    minIdx[0] = minIdx[1] = minIdx[2] = INT_MAX;
    maxIdx[0] = maxIdx[1] = maxIdx[2] = INT_MIN;

    for (auto &cube : activeCubes)
    {
        if (cube.i < minIdx[0])
            minIdx[0] = cube.i;
        if (cube.i > maxIdx[0])
            maxIdx[0] = cube.i;

        if (cube.j < minIdx[1])
            minIdx[1] = cube.j;
        if (cube.j > maxIdx[1])
            maxIdx[1] = cube.j;

        if (cube.k < minIdx[2])
            minIdx[2] = cube.k;
        if (cube.k > maxIdx[2])
            maxIdx[2] = cube.k;
    }
    std::vector<std::vector<GRID_FACETS>> grid_facets(3, std::vector<GRID_FACETS>(2,
                                                                                  GRID_FACETS(0, 0, minIdx, maxIdx)));

    // re-construct them properly with the correct (d, side):
    for (int d = 0; d < 3; d++)
    {
        for (int side = 0; side < 2; side++)
        {
            grid_facets[d][side] = GRID_FACETS(d, side, minIdx, maxIdx);
        }
    }

    // Populate them
    for (auto &cube : activeCubes)
    {
        // Global index
        int g[3] = {cube.i, cube.j, cube.k};

        for (int d = 0; d < 3; d++)
        {
            int d1 = (d + 1) % 3;
            int d2 = (d + 2) % 3;

            for (int side = 0; side < 2; side++)
            {
                GRID_FACETS &f = grid_facets[d][side];

                // Convert to local indices
                // localCoord = g - minIndex
                int coord0 = g[d1] - f.minIndex[d1];
                int coord1 = g[d2] - f.minIndex[d2];

                // Mark it
                f.SetFlag(coord0, coord1, true);
            }
        }
    }

    return grid_facets;
}