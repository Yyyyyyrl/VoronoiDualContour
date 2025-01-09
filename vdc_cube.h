//! @file vdc_cube.h
//! @brief Header file for the Cube structure and related functions.

#ifndef VDC_CUBE_H
#define VDC_CUBE_H

#include "vdc_type.h"

//! @brief Represents a cube in 3D space.
/*!
 * A cube is defined by its representative vertex, center, side length, and grid indices.
 */
struct Cube {
    Point repVertex; //!< Representative vertex of the cube.
    Point center;    //!< Center of the cube.
    int sidelength;  //!< Side length of the cube.
    int i, j, k;     //!< Grid indices of the cube.

    //! @brief Default constructor.
    Cube() : repVertex(0, 0, 0), center(0, 0, 0), sidelength(1), i(0), j(0), k(0) {}

    //! @brief Parameterized constructor.
    /*!
     * @param v Representative vertex.
     * @param c Center of the cube.
     * @param len Side length of the cube.
     * @param ix Grid index along the x-axis.
     * @param iy Grid index along the y-axis.
     * @param iz Grid index along the z-axis.
     */
    Cube(Point v, Point c, int len, int ix, int iy, int iz)
        : repVertex(v), center(c), sidelength(len), i(ix), j(iy), k(iz) {}
};

//! @brief Checks if two cubes are adjacent.
/*!
 * Two cubes are adjacent if their representative vertices differ by at most 1
 * along each axis.
 * 
 * @param cubeA The first cube.
 * @param cubeB The second cube.
 * @return `true` if the cubes are adjacent, `false` otherwise.
 */
bool is_adjacent(const Cube &cubeA, const Cube &cubeB);

//! @brief Calculates the unique index of a cube in a 3D grid.
/*!
 * @param repVertex Representative vertex of the cube.
 * @param nx Number of grid cells along the x-axis.
 * @param ny Number of grid cells along the y-axis.
 * @return The unique index of the cube in the grid.
 */
int get_cube_index(const Point &repVertex, int nx, int ny);

//! @brief Finds the indices of neighboring cubes in a 3D grid.
/*!
 * @param repVertex Representative vertex of the cube.
 * @param nx Number of grid cells along the x-axis.
 * @param ny Number of grid cells along the y-axis.
 * @return A vector of unique indices of neighboring cubes.
 */
std::vector<int> find_neighbor_indices(const Point3 &repVertex, int nx, int ny);

//! @brief Retrieves the centers of a list of cubes.
/*!
 * @param cubes A vector of cubes.
 * @return A vector of points representing the centers of the cubes.
 */
std::vector<Point> get_cube_centers(const std::vector<Cube> &cubes);

//! @brief Separates active cubes using a greedy approach.
/*!
 * Ensures that no two adjacent cubes are included in the resulting set.
 * 
 * @param activeCubes A vector of active cubes.
 * @param nx Number of grid cells along the x-axis.
 * @param ny Number of grid cells along the y-axis.
 * @param nz Number of grid cells along the z-axis.
 * @return A vector of separated active cubes.
 */
std::vector<Cube> separate_active_cubes_greedy(std::vector<Cube> &activeCubes, int nx, int ny, int nz);

//! @brief Separates active cubes using a graph-based approach.
/*!
 * Creates a graph where cubes are nodes and edges represent adjacency.
 * Uses graph coloring to ensure no two adjacent cubes are included in the result.
 * 
 * @param activeCubes A vector of active cubes.
 * @return A vector of separated active cubes.
 */
std::vector<Cube> separate_active_cubes_graph(std::vector<Cube> &activeCubes);

#endif // VDC_CUBE_H 