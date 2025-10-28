//! @file vdc_sep_isov.h
//! @brief Isosurface vertex separation algorithms
//! @details Provides different methods to separate active cubes to ensure
//! non-adjacent isosurface vertices for improved triangulation quality.

#ifndef VDC_SEP_ISOV_H
#define VDC_SEP_ISOV_H

#include "core/vdc_type.h"
#include "processing/vdc_grid.h"
#include <vector>
#include <unordered_map>

//! @brief Separation method I: Cube-level separation (original)
/*!
 * Filters active cubes to ensure no two selected cubes share a vertex, edge, or face.
 * Uses 26-connectivity in the cube lattice.
 *
 * @param activeCubes Input vector of active cubes (modified by sorting)
 * @param grid Grid for bounds and neighbor checks
 * @return Separated subset of cubes
 */
std::vector<Cube> separate_active_cubes_I(
    std::vector<Cube> &activeCubes,
    const UnifiedGrid &grid,
    float isovalue);

//! @brief Separation method III: Subgrid-based separation (3×3×3 refinement)
/*!
 * Divides each cube into 3×3×3 subcells and uses iso-crossing point location
 * to determine conflicts at subcell level in a global 3× refined grid.
 *
 * Uses accurate iso-crossing points (centroid of edge-isovalue intersections)
 * as Delaunay input for geometric accuracy.
 *
 * NOTE: The 3× grid mapping allows higher cube retention (70-85%) compared to
 * method I (~15-20%), but this can create dense point sets that lead to complex
 * Voronoi diagrams and occasional thin triangles in the final isosurface. The
 * thin triangle issue is caused by Voronoi topology (too many nearby points),
 * not by the choice of Delaunay input points (iso-crossings vs cube centers).
 *
 * @param activeCubes Input vector of active cubes
 * @param grid Grid for iso-crossing point computation
 * @param isovalue Isovalue for computing accurate crossing points
 * @return Separated subset of cubes with accurate iso-crossings as Delaunay points
 */
std::vector<Cube> separate_active_cubes_III(
    std::vector<Cube> &activeCubes,
    const UnifiedGrid &grid,
    float isovalue);

//! @brief Separation method III (testing variant): widened clearance in 3× grid (5×5×5 neighborhood)
/*!
 * This helper mirrors method III but rejects cubes if any neighbor within a
 * Chebyshev distance of 2 (a 5×5×5 block, 124 neighbors) in the refined grid
 * is already selected.
 */
std::vector<Cube> separate_active_cubes_III_wide(
    std::vector<Cube> &activeCubes,
    const UnifiedGrid &grid,
    float isovalue);

//! @brief Separation method III variant: exact binary fractions (1/4, 1/2, 3/4)
/*!
 * Identical to sep_isov_3 but uses exact binary fractional offsets instead of
 * thirds (1/6, 1/2, 5/6) to avoid floating-point representation errors that may
 * cause CGAL to generate degenerate Delaunay tetrahedra.
 *
 * The offset change is: 1/6→1/4, 5/6→3/4 (difference of ±1/12 ≈ 8.3% cube width).
 * All three offsets (0.25, 0.5, 0.75) have exact binary representations.
 *
 * @param activeCubes Input vector of active cubes
 * @param grid Grid for iso-crossing point computation
 * @param isovalue Isovalue for computing accurate crossing points
 * @return Separated subset of cubes with exact binary offsets
 */
std::vector<Cube> separate_active_cubes_III_exact_binary(
    std::vector<Cube> &activeCubes,
    const UnifiedGrid &grid,
    float isovalue);

//! @brief Compute accurate iso-crossing point using edge-intersection centroids
/*!
 * Finds all edge-isovalue intersections and returns their centroid.
 * Used for separation logic, not for Delaunay triangulation.
 *
 * @param grid The scalar grid
 * @param i,j,k Cube indices
 * @param isovalue Target isovalue
 * @return Centroid of edge-crossing points
 */
Point compute_iso_crossing_point_accurate(
    const UnifiedGrid &grid,
    int i, int j, int k,
    float isovalue);

//! @brief Determine subgrid index from iso-crossing point location
/*!
 * Maps the iso-crossing point position within a cube to a subcell index [0-26].
 * Uses 3×3×3 subdivision with base-3 encoding.
 *
 * @param isoCrossingPoint The iso-crossing point in world coordinates
 * @param cube The cube containing the point
 * @param grid Grid for spacing information
 * @return Subgrid index in range [0,26]
 */
unsigned char determine_subgrid_index(
    const Point &isoCrossingPoint,
    const Cube &cube,
    const UnifiedGrid &grid);

//! @brief Convert subgrid index to local subcell coordinates
/*!
 * Decodes base-3 subgrid index to (loc[0], loc[1], loc[2]) where each ∈ {0,1,2}.
 *
 * @param subgrid_index Index in range [0,26]
 * @param loc Output array of size 3
 */
void compute_subgrid_loc(int subgrid_index, int loc[3]);

//! @brief Compute linear index in 3× refined grid
/*!
 * @param i,j,k Coordinates in 3× refined grid
 * @param grid Original grid (for dimensions)
 * @return Linear index in 3× grid space
 */
int linear_cell_index3x(int i, int j, int k, const UnifiedGrid &grid);

//! @brief Check if cube conflicts with already selected cubes (sep_isov_3)
/*!
 * Checks connectivity in the 3× refined grid space within a configurable
 * Chebyshev radius (clearance) measured in subgrid units.
 *
 * @param cube Cube to check
 * @param selected_indices Set of selected subcell indices in 3× grid
 * @param grid Grid for dimension calculations
 * @param clearance Clearance radius (1 => 26-neighborhood, 2 => 5×5×5, etc.)
 * @param[out] indexA Index of this cube in 3× grid (if not nullptr)
 * @return true if conflicts with any selected cube
 */
bool does_cell_conflict_with_selected_cubes(
    const Cube &cube,
    const std::unordered_map<int, Cube> &selected_indices,
    const UnifiedGrid &grid,
    int *indexA = nullptr,
    int clearance = 1);

// ============================================================================
// DEBUG UTILITIES (temporary)
// ============================================================================

//! @brief Check all triples of iso-crossing points for near-collinearity
void check_collinear_isocrossings(const std::vector<Cube>& activeCubes,
                                  const UnifiedGrid& grid,
                                  float isovalue,
                                  double angle_threshold = 5.0);

#endif // VDC_SEP_ISOV_H
