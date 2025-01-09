#ifndef VDC_GRID_H
#define VDC_GRID_H

#include "vdc_cube.h"
#include "vdc_type.h"

// DIM = 3 for a 3D grid
static const int DIM3 = 3;

//! @brief A structure representing a 3D grid for storing scalar data.
/*!
 * This structure is designed for loading and processing NRRD (nearly raw raster data) format data.
 * It stores the scalar values in a 1D vector and provides metadata about the grid dimensions
 * and voxel spacing.
 */
struct Grid
{
    //! @brief 1D vector storing the scalar values of the grid.
    /*!
     * The values are stored in a flattened format for efficiency.
     */
    std::vector<float> data;

    //! @brief Number of grid cells along the x, y, and z axes.
    int nx, ny, nz;

    //! @brief Spacing of the grid cells along the x, y, and z axes.
    float dx, dy, dz;

    //! @brief Print the grid's metadata and data.
    /*!
     * This function outputs the grid's dimensions, voxel spacing, and scalar data in a human-readable format.
     */
    void print_grid();
};

//! @brief The Scalar Grid Data structure used in the algorithm.
/*!
 * This structure represents a 3D grid of scalar values. It provides methods
 * to manipulate and query the scalar data at different grid points.
 */
struct ScalarGrid
{
    //! @brief 3D vector to store scalar values.
    std::vector<std::vector<std::vector<float>>> data;

    //! @brief Dimensions of the grid along the x, y, and z axes.
    int nx, ny, nz;

    //! @brief Voxel dimensions along the x, y, and z axes.
    float dx, dy, dz;

    //! @brief Minimum and maximum coordinates of the grid.
    float min_x, min_y, min_z, max_x, max_y, max_z;

    //! @brief Constructor to initialize a scalar grid.
    /*!
     * @param nx Number of grid cells along the x-axis.
     * @param ny Number of grid cells along the y-axis.
     * @param nz Number of grid cells along the z-axis.
     * @param dx Voxel size along the x-axis.
     * @param dy Voxel size along the y-axis.
     * @param dz Voxel size along the z-axis.
     * @param min_x Minimum x-coordinate of the grid.
     * @param min_y Minimum y-coordinate of the grid.
     * @param min_z Minimum z-coordinate of the grid.
     */
    ScalarGrid(int nx, int ny, int nz, float dx, float dy, float dz, float min_x, float min_y, float min_z);

    //! @brief Get a scalar value at the specified grid index.
    /*!
     * @param x Grid index along the x-axis.
     * @param y Grid index along the y-axis.
     * @param z Grid index along the z-axis.
     * @return Scalar value at the specified index. Returns 0 if the index is out of bounds.
     */
    float get_value(int x, int y, int z) const;

    //! @brief Set a scalar value at the specified grid index.
    /*!
     * @param x Grid index along the x-axis.
     * @param y Grid index along the y-axis.
     * @param z Grid index along the z-axis.
     * @param value Scalar value to set at the specified index.
     */
    void set_value(int x, int y, int z, float value);

    //! @brief Load scalar grid data from an external source.
    /*!
     * @param source A 3D vector containing scalar values to load into the grid.
     */
    void load_from_source(const std::vector<std::vector<std::vector<float>>> &source);

    //! @brief Get the scalar value at a given point in space.
    /*!
     * @param point The 3D point at which to retrieve the scalar value.
     * @return Scalar value at the given point.
     */
    float get_scalar_value_at_point(const Point &point);

    //! @brief Convert a 3D point to its corresponding grid index.
    /*!
     * @param point The 3D point to convert.
     * @return A tuple containing the grid indices (x, y, z).
     */
    std::tuple<int, int, int> point_to_grid_index(const Point &point);
};


//! @brief A structure representing a 2D "facet" of a 3D grid, orthogonal to one axis.
/*!
 * For example, if `orth_dir = 0` (the x-axis), then this facet stores data 
 * along the y- and z-axes. The `axis_dir[0]` corresponds to the first axis, 
 * and `axis_dir[1]` corresponds to the second axis.
 */
struct GRID_FACETS 
{
    //! @brief Orthogonal direction: 0=x, 1=y, 2=z.
    int orth_dir;

    //! @brief Side of the grid facet: 0=lower (min), 1=upper (max).
    int side;

    //! @brief The two axes defining this facet (e.g., y and z if orth_dir = 0).
    int axis_dir[2];

    //! @brief Size of the facet along the two axes.
    /*!
     * - `axis_size[0]`: Size along `axis_dir[0]`.
     * - `axis_size[1]`: Size along `axis_dir[1]`.
     */
    int axis_size[2];

    //! @brief Minimum global indices for the bounding box in each dimension.
    int minIndex[DIM3];

    //! @brief Maximum global indices for the bounding box in each dimension.
    int maxIndex[DIM3];

    //! @brief Local size of the facet in each dimension.
    int localSize[DIM3];

    //! @brief Boolean flags for each `(coord0, coord1)` in the facet.
    /*!
     * The flags are stored in a linearized format for efficient memory usage.
     */
    std::vector<bool> cube_flag;

    //! @brief Constructor to initialize a grid facet.
    /*!
     * @param d Orthogonal direction (0=x, 1=y, 2=z).
     * @param s Side of the grid (0=lower, 1=upper).
     * @param minIdx Global minimum indices for the bounding box.
     * @param maxIdx Global maximum indices for the bounding box.
     */
    GRID_FACETS(int d, int s, const int minIdx[DIM3], const int maxIdx[DIM3]);

    //! @brief Default destructor (uses `std::vector<bool>` for automatic cleanup).
    ~GRID_FACETS() = default;

    //! @brief Set the flag for a particular `(coord0, coord1)` in this facet.
    /*!
     * @param coord0 Coordinate along the first axis (`axis_dir[0]`).
     * @param coord1 Coordinate along the second axis (`axis_dir[1]`).
     * @param flag Boolean value to set (true or false).
     */
    void SetFlag(int coord0, int coord1, bool flag);

    //! @brief Get the flag for a particular `(coord0, coord1)` in this facet.
    /*!
     * @param coord0 Coordinate along the first axis (`axis_dir[0]`).
     * @param coord1 Coordinate along the second axis (`axis_dir[1]`).
     * @return The boolean flag at `(coord0, coord1)`.
     */
    bool CubeFlag(int coord0, int coord1) const;

private:
    //! @brief Convert a 2D coordinate `(coord0, coord1)` to a linear index.
    /*!
     * @param coord0 Coordinate along the first axis (`axis_dir[0]`).
     * @param coord1 Coordinate along the second axis (`axis_dir[1]`).
     * @return The linearized index for accessing `cube_flag`.
     */
    int index(int coord0, int coord1) const;
};


// Functions for loading nrrd data

//! @brief Converts raw data of type `T` into a float vector.
/*!
 * @tparam T The type of the input data.
 * @param data_ptr Pointer to the input data.
 * @param total_size Total number of elements in the input data.
 * @return A `std::vector<float>` containing the converted data.
 */
template <typename T>
std::vector<float> convert_to_float_vector(T *data_ptr, size_t total_size);

//! @brief Adjusts a point that lies outside the bounds of the scalar grid.
/*!
 * @param p The input point to adjust.
 * @param grid The scalar grid defining the bounds.
 * @param v1 The start point of the segment.
 * @param v2 The end point of the segment.
 * @return A point adjusted to lie within the bounds of the grid.
 */
Point adjust_outside_bound_points(const Point &p, const ScalarGrid &grid, const Point &v1, const Point &v2);

/*Preprocessing*/

//! @brief Initializes a `ScalarGrid` from a `Grid` object.
/*!
 * @param grid The scalar grid to initialize.
 * @param nrrdGrid The input grid containing scalar data.
 */
void initialize_scalar_grid(ScalarGrid &grid, const Grid &nrrdGrid);

//! @brief Loads NRRD data into a `Grid` structure.
/*!
 * @param file_path The path to the NRRD file.
 * @return A `Grid` object containing the loaded data.
 */
Grid load_nrrd_data(const std::string &file_path);

//! @brief Supersamples a `Grid` by a factor of `n`.
/*!
 * @param grid The input grid to supersample.
 * @param n The supersampling factor.
 * @return A supersampled `Grid`.
 */
Grid supersample_grid(const Grid &grid, int n);

//! @brief Checks if a cube is active based on its scalar values.
/*!
 * @param grid The grid containing the scalar values.
 * @param x The x-coordinate of the cube.
 * @param y The y-coordinate of the cube.
 * @param z The z-coordinate of the cube.
 * @param isovalue The isovalue for activity determination.
 * @return `true` if the cube is active; `false` otherwise.
 */
bool is_cube_active(const Grid &grid, int x, int y, int z, float isovalue);

//! @brief Finds all active cubes in the grid based on an isovalue.
/*!
 * @param grid The input grid.
 * @param isovalue The isovalue for activity determination.
 * @param cubes A vector to store the active cubes.
 */
void find_active_cubes(const Grid &grid, float isovalue, std::vector<Cube> &cubes);


//! @brief Loads the grid points from a `Grid`.
/*!
 * @param grid The input grid.
 * @return A vector of points representing the grid points.
 */
std::vector<Point> load_grid_points(const Grid &grid);

//! @brief Checks if a point is inside the bounds of the grid.
/*!
 * @param p The point to check.
 * @param grid The scalar grid defining the bounds.
 * @return `true` if the point is inside the grid; `false` otherwise.
 */
bool is_point_inside_grid(const Point &p, const ScalarGrid &grid);

//! @brief Interpolates a point on a line segment based on scalar values.
/*!
 * @param p1 The first point of the segment.
 * @param p2 The second point of the segment.
 * @param val1 The scalar value at `p1`.
 * @param val2 The scalar value at `p2`.
 * @param isovalue The target scalar value for interpolation.
 * @param data_grid The grid containing the scalar data.
 * @return The interpolated point.
 */
Point interpolate(const Point &p1, const Point &p2, float val1, float val2, float isovalue, const Grid &data_grid);

//! @brief Performs trilinear interpolation at a given point in a scalar grid.
/*!
 * @param p The point for interpolation.
 * @param grid The scalar grid containing the data.
 * @return The interpolated scalar value.
 */
float trilinear_interpolate(const Point &p, const ScalarGrid &grid);

//! @brief Performs trilinear interpolation at a given point in a regular grid.
/*!
 * @param p The point for interpolation.
 * @param grid The regular grid containing the data.
 * @return The interpolated scalar value.
 */
float trilinear_interpolate(const Point &p, const Grid &grid);

//! @brief Creates grid facets for a given set of active cubes.
/*!
 * @param activeCubes A vector of active cubes.
 * @return A 3D vector of `GRID_FACETS`, grouped by direction and side.
 */
std::vector<std::vector<GRID_FACETS>> create_grid_facets(const std::vector<Cube> &activeCubes);

#endif