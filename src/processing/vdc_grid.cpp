#include "processing/vdc_grid.h"
#include "core/vdc_timing.h"
#include <algorithm>
#include <cmath>
#include <limits>


//! Constructor for UnifiedGrid
UnifiedGrid::UnifiedGrid(int nx, int ny, int nz, float dx, float dy, float dz, float min_x, float min_y, float min_z)
    : nx(nx), ny(ny), nz(nz),
      dx(dx), dy(dy), dz(dz),
      physical_dx(dx), physical_dy(dy), physical_dz(dz),
      min_x(min_x), min_y(min_y), min_z(min_z)
{
    update_bounds();
    flat_data.resize(nx * ny * nz, 0.0f);
    data.resize(nx, std::vector<std::vector<float>>(ny, std::vector<float>(nz, 0.0f)));
}

// Retrieve a scalar value from the grid
float UnifiedGrid::get_value(int x, int y, int z) const
{
    if (x < 0 || x >= nx || y < 0 || y >= ny || z < 0 || z >= nz)
        return 0.0f;
    return data[x][y][z];
}

// Set a scalar value in the grid
void UnifiedGrid::set_value(int x, int y, int z, float value)
{
    if (x >= 0 && x < nx && y >= 0 && y < ny && z >= 0 && z < nz)
    {
        data[x][y][z] = value;
        flat_data[z * nx * ny + y * nx + x] = value;
    }
}


// Convert point to grid index
std::tuple<int, int, int> UnifiedGrid::point_to_grid_index(const Point &point) const
{
    int x = static_cast<int>((point.x() - min_x) / dx);
    int y = static_cast<int>((point.y() - min_y) / dy);
    int z = static_cast<int>((point.z() - min_z) / dz);
    return {x, y, z};
}

// Get scalar value at a point using trilinear interpolation
float UnifiedGrid::get_scalar_value_at_point(const Point &point) const
{
    return trilinear_interpolate(point, *this);
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

// Print grid metadata and data
void UnifiedGrid::print_grid() const
{
    std::cout << "Unified Grid Information:\n";
    std::cout << "Dimensions: " << nx << "x" << ny << "x" << nz << "\n";
    std::cout << "Internal spacing (grid units): dx=" << dx << ", dy=" << dy << ", dz=" << dz << "\n";
    std::cout << "Physical spacing: dx=" << physical_dx << ", dy=" << physical_dy << ", dz=" << physical_dz << "\n";
    std::cout << "Bounds (grid units): [" << min_x << ", " << max_x << "] x [" << min_y << ", " << max_y << "] x [" << min_z << ", " << max_z << "]\n";
    const float phys_max_x = min_x + (nx - 1) * physical_dx;
    const float phys_max_y = min_y + (ny - 1) * physical_dy;
    const float phys_max_z = min_z + (nz - 1) * physical_dz;
    std::cout << "Bounds (physical): [" << min_x << ", " << phys_max_x << "] x [" << min_y << ", " << phys_max_y << "] x [" << min_z << ", " << phys_max_z << "]\n\n";
    std::cout << "Data:\n";
    for (int z = 0; z < nz; ++z)
    {
        std::cout << "Slice z = " << z << ":\n";
        for (int y = 0; y < ny; ++y)
        {
            for (int x = 0; x < nx; ++x)
                std::cout << std::setw(8) << data[x][y][z] << " ";
            std::cout << "\n";
        }
        std::cout << "\n";
    }
}

void UnifiedGrid::update_bounds()
{
    max_x = min_x + (nx - 1) * dx;
    max_y = min_y + (ny - 1) * dy;
    max_z = min_z + (nz - 1) * dz;
}

void UnifiedGrid::force_unit_spacing()
{
    dx = dy = dz = 1.0f;
    update_bounds();
}

void UnifiedGrid::zero_boundary_shell()
{
    if (nx <= 0 || ny <= 0 || nz <= 0)
        return;

    const float fill_value = flat_data.empty()
                                 ? 0.0f
                                 : *std::min_element(flat_data.begin(), flat_data.end());

    auto assign = [&](int x, int y, int z) {
        if (x < 0 || x >= nx || y < 0 || y >= ny || z < 0 || z >= nz)
            return;
        data[x][y][z] = fill_value;
        flat_data[z * nx * ny + y * nx + x] = fill_value;
    };

    for (int x = 0; x < nx; ++x)
    {
        for (int y = 0; y < ny; ++y)
        {
            assign(x, y, 0);
            assign(x, y, nz - 1);
        }
    }

    for (int x = 0; x < nx; ++x)
    {
        for (int z = 0; z < nz; ++z)
        {
            assign(x, 0, z);
            assign(x, ny - 1, z);
        }
    }

    for (int y = 0; y < ny; ++y)
    {
        for (int z = 0; z < nz; ++z)
        {
            assign(0, y, z);
            assign(nx - 1, y, z);
        }
    }
}

bool UnifiedGrid::boundary_crosses_isovalue(float isovalue) const
{
    if (nx <= 0 || ny <= 0 || nz <= 0)
        return false;

    float boundary_min = std::numeric_limits<float>::infinity();
    float boundary_max = -std::numeric_limits<float>::infinity();

    auto consider = [&](int x, int y, int z)
    {
        if (x < 0 || x >= nx || y < 0 || y >= ny || z < 0 || z >= nz)
            return;
        float val = data[x][y][z];
        boundary_min = std::min(boundary_min, val);
        boundary_max = std::max(boundary_max, val);
    };

    for (int x = 0; x < nx; ++x)
    {
        for (int y = 0; y < ny; ++y)
        {
            consider(x, y, 0);
            consider(x, y, nz - 1);
        }
    }

    for (int x = 0; x < nx; ++x)
    {
        for (int z = 1; z < nz - 1; ++z)
        {
            consider(x, 0, z);
            consider(x, ny - 1, z);
        }
    }

    for (int y = 1; y < ny - 1; ++y)
    {
        for (int z = 1; z < nz - 1; ++z)
        {
            consider(0, y, z);
            consider(nx - 1, y, z);
        }
    }

    if (!std::isfinite(boundary_min) || !std::isfinite(boundary_max))
        return false;

    return (boundary_min < isovalue) && (isovalue < boundary_max);
}


//! @brief Converts raw data of type `T` into a float vector.
/*!
 * @tparam T The type of the input data.
 * @param data_ptr Pointer to the input data.
 * @param total_size Total number of elements in the input data.
 * @return A `std::vector<float>` containing the converted data.
 */
// Convert raw data to float vector
template <typename T>
std::vector<float> convert_to_float_vector(T *data_ptr, size_t total_size)
{
    std::vector<float> data(total_size);
    for (size_t i = 0; i < total_size; ++i)
        data[i] = static_cast<float>(data_ptr[i]);
    return data;
}

// Load NRRD data
UnifiedGrid load_nrrd_data(const std::string &file_path)
{
    TimingStats& timer = TimingStats::getInstance();
    timer.startTimer("Load NRRD file", "1. Load Data and Grid Formation");

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
    int nx = nrrd->axis[0].size;
    int ny = nrrd->axis[1].size;
    int nz = nrrd->axis[2].size;
    auto sanitize_spacing = [](double spacing) -> float {
        if (!std::isfinite(spacing) || spacing <= 0.0)
            return 1.0f;
        return static_cast<float>(spacing);
    };

    float dx = sanitize_spacing(nrrd->axis[0].spacing);
    float dy = sanitize_spacing(nrrd->axis[1].spacing);
    float dz = sanitize_spacing(nrrd->axis[2].spacing);
    float min_x = 0.0f, min_y = 0.0f, min_z = 0.0f;
    timer.stopTimer("Load NRRD file");

    timer.startTimer("Grid initialization", "1. Load Data and Grid Formation");
    UnifiedGrid grid(nx, ny, nz, dx, dy, dz, min_x, min_y, min_z);

    if (nrrd->type == nrrdTypeFloat)
    {
        float *data_ptr = static_cast<float *>(nrrd->data);
        grid.flat_data = std::vector<float>(data_ptr, data_ptr + total_size);
    }
    else if (nrrd->type == nrrdTypeUChar)
    {
        unsigned char *data_ptr = static_cast<unsigned char *>(nrrd->data);
        grid.flat_data = convert_to_float_vector(data_ptr, total_size);
    }
    else
    {
        std::cerr << "Unsupported NRRD data type." << std::endl;
        nrrdNuke(nrrd);
        exit(1);
    }

    // Populate 3D data array
    for (int z = 0; z < nz; ++z)
        for (int y = 0; y < ny; ++y)
            for (int x = 0; x < nx; ++x)
                grid.data[x][y][z] = grid.flat_data[z * nx * ny + y * nx + x];

    nrrdNuke(nrrd);

    grid.force_unit_spacing();
    timer.stopTimer("Grid initialization");

    std::cout << "Grid dimensions: " << nx << "x" << ny << "x" << nz << "\n";
    std::cout << "Physical spacing: dx=" << grid.physical_dx << ", dy=" << grid.physical_dy << ", dz=" << grid.physical_dz << "\n";
    std::cout << "Internal spacing (grid units): dx=" << grid.dx << ", dy=" << grid.dy << ", dz=" << grid.dz << "\n";
    const float phys_max_x = grid.min_x + (grid.nx - 1) * grid.physical_dx;
    const float phys_max_y = grid.min_y + (grid.ny - 1) * grid.physical_dy;
    const float phys_max_z = grid.min_z + (grid.nz - 1) * grid.physical_dz;
    std::cout << "Bounds (grid units): [" << grid.min_x << ", " << grid.max_x << "] x ["
              << grid.min_y << ", " << grid.max_y << "] x [" << grid.min_z << ", " << grid.max_z << "]\n";
    std::cout << "Bounds (physical): [" << grid.min_x << ", " << phys_max_x << "] x ["
              << grid.min_y << ", " << phys_max_y << "] x [" << grid.min_z << ", " << phys_max_z << "]\n";

    return grid;
}

// Supersample grid
UnifiedGrid supersample_grid(const UnifiedGrid &grid, int n)
{
    TimingStats& timer = TimingStats::getInstance();
    timer.startTimer("Supersample", "1. Load Data and Grid Formation");

    int nx2 = grid.nx * n - (n - 1);
    int ny2 = grid.ny * n - (n - 1);
    int nz2 = grid.nz * n - (n - 1);
    float dx2 = grid.physical_dx / n;
    float dy2 = grid.physical_dy / n;
    float dz2 = grid.physical_dz / n;

    UnifiedGrid new_grid(nx2, ny2, nz2, dx2, dy2, dz2, grid.min_x, grid.min_y, grid.min_z);

    for (int z = 0; z < nz2; ++z)
    {
        for (int y = 0; y < ny2; ++y)
        {
            for (int x = 0; x < nx2; ++x)
            {
                float px = grid.min_x + (static_cast<float>(x) / n) * grid.dx;
                float py = grid.min_y + (static_cast<float>(y) / n) * grid.dy;
                float pz = grid.min_z + (static_cast<float>(z) / n) * grid.dz;
                float value = trilinear_interpolate(Point(px, py, pz), grid);
                new_grid.set_value(x, y, z, value);
            }
        }
    }

    new_grid.force_unit_spacing();
    timer.stopTimer("Supersample");

    return new_grid;
}


// Check if cube is active
bool is_cube_active(const UnifiedGrid &grid, int x, int y, int z, float isovalue)
{
    static const std::vector<std::tuple<int,int,int>> vertex_offsets = {
        {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
        {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}
    };

    auto get_value = [&](int dx, int dy, int dz) -> float {
        return grid.get_value(x + dx, y + dy, z + dz);
    };

    float val0 = get_value(0, 0, 0);
    bool is_val0_negative = (val0 < isovalue);

    for (int i = 1; i < 8; i++) {
        auto [dx, dy, dz] = vertex_offsets[i];
        float vali = get_value(dx, dy, dz);
        bool is_vali_negative = (vali < isovalue);
        if (is_vali_negative != is_val0_negative)
            return true;
    }

    return false;
}

// Adjust points outside grid bounds
Point adjust_outside_bound_points(const Point &p, const UnifiedGrid &grid, const Point &v1, const Point &v2)
{
    float gx = (p.x() - grid.min_x) / grid.dx;
    float gy = (p.y() - grid.min_y) / grid.dy;
    float gz = (p.z() - grid.min_z) / grid.dz;

    if (gx < 0 || gx >= grid.nx || gy < 0 || gy >= grid.ny || gz < 0 || gz >= grid.nz)
    {
        float v1_gx = (v1.x() - grid.min_x) / grid.dx;
        float v1_gy = (v1.y() - grid.min_y) / grid.dy;
        float v1_gz = (v1.z() - grid.min_z) / grid.dz;

        float v2_gx = (v2.x() - grid.min_x) / grid.dx;
        float v2_gy = (v2.y() - grid.min_y) / grid.dy;
        float v2_gz = (v2.z() - grid.min_z) / grid.dz;

        float t = ((gx - v1_gx) * (v2_gx - v1_gx) + (gy - v1_gy) * (v2_gy - v1_gy) + (gz - v1_gz) * (v2_gz - v1_gz)) /
                  ((v2_gx - v1_gx) * (v2_gx - v1_gx) + (v2_gy - v1_gy) * (v2_gy - v1_gy) + (v2_gz - v1_gz) * (v2_gz - v1_gz));

        t = std::max(0.0f, std::min(t, 1.0f));

        float px = v1_gx + t * (v2_gx - v1_gx);
        float py = v1_gy + t * (v2_gy - v1_gy);
        float pz = v1_gz + t * (v2_gz - v1_gz);

        return Point(px * grid.dx + grid.min_x, py * grid.dy + grid.min_y, pz * grid.dz + grid.min_z);
    }

    return p;
}


// Helper function to compute iso-crossing point in an active cube
// Note: Using cube center for stability. Edge-based centroids can cause
// degenerate Delaunay configurations when adjacent cubes share edge crossings.
Point compute_iso_crossing_point(const UnifiedGrid &grid, int i, int j, int k, float isovalue)
{
    // Return cube center - this ensures all points are distinct and well-separated,
    // avoiding degeneracies in the Delaunay triangulation
    return Point(
        (i + 0.5f) * grid.dx + grid.min_x,
        (j + 0.5f) * grid.dy + grid.min_y,
        (k + 0.5f) * grid.dz + grid.min_z
    );
}

// Find active cubes
void find_active_cubes(const UnifiedGrid &grid, float isovalue, std::vector<Cube> &cubes)
{
    cubes.clear();
    for (int i = 0; i < grid.nx - 1; ++i)
    {
        for (int j = 0; j < grid.ny - 1; ++j)
        {
            for (int k = 0; k < grid.nz - 1; ++k)
            {
                if (is_cube_active(grid, i, j, k, isovalue))
                {
                    Point repVertex(i * grid.dx + grid.min_x, j * grid.dy + grid.min_y, k * grid.dz + grid.min_z);
                    Point cubeCenter(
                        (i + 0.5f) * grid.dx + grid.min_x,
                        (j + 0.5f) * grid.dy + grid.min_y,
                        (k + 0.5f) * grid.dz + grid.min_z);
                    Cube cube(repVertex, cubeCenter, i, j, k);
                    cube.accurateIsoCrossing = compute_iso_crossing_point(grid, i, j, k, isovalue);
                    cubes.push_back(cube);
                }
            }
        }
    }
}


// Load grid points
std::vector<Point> load_grid_points(const UnifiedGrid &grid)
{
    std::vector<Point> points;
    for (int i = 0; i < grid.nx; ++i)
        for (int j = 0; j < grid.ny; ++j)
            for (int k = 0; k < grid.nz; ++k)
                points.push_back(Point(i * grid.dx + grid.min_x, j * grid.dy + grid.min_y, k * grid.dz + grid.min_z));
    return points;
}

// Check if point is inside grid
bool is_point_inside_grid(const Point &p, const UnifiedGrid &grid)
{
    return (p.x() >= grid.min_x && p.x() <= grid.max_x &&
            p.y() >= grid.min_y && p.y() <= grid.max_y &&
            p.z() >= grid.min_z && p.z() <= grid.max_z);
}

// Interpolate along an edge
Point interpolate(const Point &p1, const Point &p2, float val1, float val2, float isovalue, const UnifiedGrid &grid)
{
    if (std::abs(val1 - val2) < 1e-6)
        return p1;
    float t = (isovalue - val1) / (val2 - val1);
    return Point(p1.x() + t * (p2.x() - p1.x()),
                p1.y() + t * (p2.y() - p1.y()),
                p1.z() + t * (p2.z() - p1.z()));
}

// Trilinear interpolation
float trilinear_interpolate(const Point &p, const UnifiedGrid &grid)
{
    float gx = (p.x() - grid.min_x) / grid.dx;
    float gy = (p.y() - grid.min_y) / grid.dy;
    float gz = (p.z() - grid.min_z) / grid.dz;
    gx = std::max(0.0f, std::min(gx, (float)(grid.nx - 1)));
    gy = std::max(0.0f, std::min(gy, (float)(grid.ny - 1)));
    gz = std::max(0.0f, std::min(gz, (float)(grid.nz - 1)));

    int x0 = static_cast<int>(std::floor(gx));
    int x1 = std::min(x0 + 1, grid.nx - 1);
    int y0 = static_cast<int>(std::floor(gy));
    int y1 = std::min(y0 + 1, grid.ny - 1);
    int z0 = static_cast<int>(std::floor(gz));
    int z1 = std::min(z0 + 1, grid.nz - 1);

    float xd = gx - x0;
    float yd = gy - y0;
    float zd = gz - z0;

    float c000 = grid.get_value(x0, y0, z0);
    float c001 = grid.get_value(x0, y0, z1);
    float c010 = grid.get_value(x0, y1, z0);
    float c011 = grid.get_value(x0, y1, z1);
    float c100 = grid.get_value(x1, y0, z0);
    float c101 = grid.get_value(x1, y0, z1);
    float c110 = grid.get_value(x1, y1, z0);
    float c111 = grid.get_value(x1, y1, z1);

    float c00 = c000 * (1 - zd) + c001 * zd;
    float c01 = c010 * (1 - zd) + c011 * zd;
    float c10 = c100 * (1 - zd) + c101 * zd;
    float c11 = c110 * (1 - zd) + c111 * zd;

    float c0 = c00 * (1 - yd) + c01 * yd;
    float c1 = c10 * (1 - yd) + c11 * yd;

    return c0 * (1 - xd) + c1 * xd;
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


// Check if two cubes are adjacent in grid space
bool is_adjacent(const Cube &cubeA, const Cube &cubeB, const UnifiedGrid &grid)
{
    int di = std::abs(cubeA.i - cubeB.i);
    int dj = std::abs(cubeA.j - cubeB.j);
    int dk = std::abs(cubeA.k - cubeB.k);
    return (di <= 1 && dj <= 1 && dk <= 1) && !(di == 0 && dj == 0 && dk == 0);
}

// Calculate unique cube index
int get_cube_index(const Point &repVertex, const UnifiedGrid &grid)
{
    int i = static_cast<int>((repVertex.x() - grid.min_x) / grid.dx);
    int j = static_cast<int>((repVertex.y() - grid.min_y) / grid.dy);
    int k = static_cast<int>((repVertex.z() - grid.min_z) / grid.dz);
    return k * (grid.nx - 1) * (grid.ny - 1) + j * (grid.nx - 1) + i;
}

// Find neighbor indices
std::vector<int> find_neighbor_indices(const Point &repVertex, const UnifiedGrid &grid)
{
    std::vector<int> neighbors;
    int i = static_cast<int>((repVertex.x() - grid.min_x) / grid.dx);
    int j = static_cast<int>((repVertex.y() - grid.min_y) / grid.dy);
    int k = static_cast<int>((repVertex.z() - grid.min_z) / grid.dz);
    for (int di = -1; di <= 1; ++di)
        for (int dj = -1; dj <= 1; ++dj)
            for (int dk = -1; dk <= 1; ++dk)
                if (di != 0 || dj != 0 || dk != 0)
                {
                    int ni = i + di, nj = j + dj, nk = k + dk;
                    if (ni >= 0 && ni < grid.nx - 1 && nj >= 0 && nj < grid.ny - 1 && nk >= 0 && nk < grid.nz - 1)
                        neighbors.push_back(nk * (grid.nx - 1) * (grid.ny - 1) + nj * (grid.nx - 1) + ni);
                }
    return neighbors;
}

//! @brief Retrieves the centers of a list of cubes.
std::vector<Point> get_cube_centers(const std::vector<Cube> &cubes)
{
    std::vector<Point> cubeCenters;
    for (auto &cube : cubes)
    {
        cubeCenters.push_back(cube.cubeCenter);
    }
    return cubeCenters;
}

//! @brief Retrieves the accurate iso-crossing points of a list of cubes.
std::vector<Point> get_cube_accurate_iso_crossing_points(const std::vector<Cube> &cubes)
{
    std::vector<Point> accurateIsoCrossingPoints;
    for (auto &cube : cubes)
    {
        accurateIsoCrossingPoints.push_back(cube.accurateIsoCrossing);
    }
    return accurateIsoCrossingPoints;
}

// Separation routines have been moved to vdc_sep_isov.cpp
