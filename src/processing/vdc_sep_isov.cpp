//! @file vdc_sep_isov.cpp
//! @brief Implementation of isosurface vertex separation algorithms

#include "processing/vdc_sep_isov.h"
#include <algorithm>
#include <cmath>
#include <unordered_set>

// ============================================================================
// COMMON UTILITIES
// ============================================================================

//! @brief Consistent linear index for a cell (i,j,k) in original grid
static inline int linear_cell_index(int i, int j, int k, const UnifiedGrid& g) {
    const int nx = g.num_cells[0];
    const int ny = g.num_cells[1];
    const int xy_stride = (nx - 1) * (ny - 1);
    return k * xy_stride + j * (nx - 1) + i;
}

//! @brief Compute minimum distance to bounding box boundary in grid space
static inline int min_distance_to_boundary(int i, int j, int k, const UnifiedGrid& grid) {
    const int nx = grid.num_cells[0];
    const int ny = grid.num_cells[1];
    const int nz = grid.num_cells[2];
    int dist_i = std::min(i, (nx - 2) - i);  // nx-2 is the max valid cube index in i
    int dist_j = std::min(j, (ny - 2) - j);
    int dist_k = std::min(k, (nz - 2) - k);
    return std::min({dist_i, dist_j, dist_k});
}

// ============================================================================
// SEPARATION METHOD I: General Cube-Level
// ============================================================================

std::vector<Cube> separate_active_cubes_I(
    std::vector<Cube>& activeCubes,
    const UnifiedGrid& grid,
    float isovalue)
{
    // Sort cubes by distance to boundary (ascending), prioritizing boundary cubes
    std::sort(activeCubes.begin(), activeCubes.end(),
              [&grid](const Cube& a, const Cube& b) {
                  int dist_a = min_distance_to_boundary(a.indices[0], a.indices[1], a.indices[2], grid);
                  int dist_b = min_distance_to_boundary(b.indices[0], b.indices[1], b.indices[2], grid);
                  return dist_a < dist_b;
              });

    std::unordered_set<int> kept;   // key: linear cell index
    kept.reserve(activeCubes.size());
    std::vector<Cube> out;
    out.reserve(activeCubes.size());

    const float dx = grid.spacing[0];
    const float dy = grid.spacing[1];
    const float dz = grid.spacing[2];
    const float min_x = grid.min_coord[0];
    const float min_y = grid.min_coord[1];
    const float min_z = grid.min_coord[2];
    const int nx = grid.num_cells[0];
    const int ny = grid.num_cells[1];
    const int nz = grid.num_cells[2];

    for (Cube c : activeCubes) {
        const int ci = c.indices[0], cj = c.indices[1], ck = c.indices[2];

        // Reject cubes that share any vertex/edge/face with an already kept cube
        // (i.e., any of the 26-neighbors in the (i,j,k) lattice).
        bool hasAdjKept = false;
        for (int dk = -1; dk <= 1 && !hasAdjKept; ++dk) {
            for (int dj = -1; dj <= 1 && !hasAdjKept; ++dj) {
                for (int di = -1; di <= 1 && !hasAdjKept; ++di) {
                    // skip self
                    if (di == 0 && dj == 0 && dk == 0) continue;

                    const int ni = ci + di, nj = cj + dj, nk = ck + dk;
                    if (ni < 0 || nj < 0 || nk < 0) continue;
                    if (ni >= nx - 1 || nj >= ny - 1 || nk >= nz - 1) continue;

                    const int nIdx = linear_cell_index(ni, nj, nk, grid);
                    if (kept.find(nIdx) != kept.end()) hasAdjKept = true;
                }
            }
        }

        if (!hasAdjKept) {
            // Compute accurate iso-crossing for iso-vertex computation
            c.accurateIsoCrossing = compute_iso_crossing_point_accurate(grid, ci, cj, ck, isovalue);

            c.cubeCenter = Point(
                (ci + 0.5f) * dx + min_x,
                (cj + 0.5f) * dy + min_y,
                (ck + 0.5f) * dz + min_z
            );

            const int myIdx = linear_cell_index(ci, cj, ck, grid);
            kept.insert(myIdx);
            out.push_back(c);
        }
    }
    return out;
}

// ============================================================================
// ACCURATE ISO-CROSSING POINT COMPUTATION
// ============================================================================

Point compute_iso_crossing_point_accurate(
    const UnifiedGrid &grid,
    int i, int j, int k,
    float isovalue)
{
    // Cube vertex offsets
    static const int cubeVertices[8][3] = {
        {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
        {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}
    };
    static const int cubeEdges[12][2] = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0},
        {4, 5}, {5, 6}, {6, 7}, {7, 4},
        {0, 4}, {1, 5}, {2, 6}, {3, 7}
    };

    // Compute cube's base position
    const float dx = grid.spacing[0];
    const float dy = grid.spacing[1];
    const float dz = grid.spacing[2];
    const float min_x = grid.min_coord[0];
    const float min_y = grid.min_coord[1];
    const float min_z = grid.min_coord[2];

    float base_x = i * dx + min_x;
    float base_y = j * dy + min_y;
    float base_z = k * dz + min_z;

    // Get scalar values at all 8 cube vertices (exact grid values)
    std::array<Point, 8> vertices;
    std::array<float, 8> scalarValues;
    for (int v = 0; v < 8; ++v)
    {
        vertices[v] = Point(
            base_x + cubeVertices[v][0] * dx,
            base_y + cubeVertices[v][1] * dy,
            base_z + cubeVertices[v][2] * dz
        );
        // Use exact grid values at vertices (grid points)
        scalarValues[v] = grid.get_value(
            i + cubeVertices[v][0],
            j + cubeVertices[v][1],
            k + cubeVertices[v][2]
        );
    }

    // Find all edge-isovalue intersection points
    std::vector<Point> intersectionPoints;
    for (const auto &edge : cubeEdges)
    {
        int idx1 = edge[0];
        int idx2 = edge[1];
        float val1 = scalarValues[idx1];
        float val2 = scalarValues[idx2];

        // Check if edge crosses the isovalue
        if ((val1 < isovalue && val2 >= isovalue) || (val1 >= isovalue && val2 < isovalue))
        {
            // Linear interpolation along edge
            Point intersect = interpolate(vertices[idx1], vertices[idx2], val1, val2, isovalue, grid);
            intersectionPoints.push_back(intersect);
        }
    }

    // Return centroid of intersection points, or cube center as fallback
    if (!intersectionPoints.empty())
    {
        Point centroid(0, 0, 0);
        for (const auto &pt : intersectionPoints)
        {
            centroid = Point(
                centroid.x() + pt.x(),
                centroid.y() + pt.y(),
                centroid.z() + pt.z()
            );
        }
        float n = static_cast<float>(intersectionPoints.size());
        return Point(centroid.x() / n, centroid.y() / n, centroid.z() / n);
    }
    else
    {
        // Fallback to cube center
        return Point(
            base_x + 0.5f * dx,
            base_y + 0.5f * dy,
            base_z + 0.5f * dz
        );
    }
}

// ============================================================================
// SUBGRID INDEX COMPUTATION (3×3×3 SUBDIVISION)
// ============================================================================

void compute_subgrid_loc(int subgrid_index, int loc[3])
{
    // Precondition: subgrid_index ∈ [0,26]
    int index = subgrid_index;
    loc[0] = index % 3;
    index = index / 3;  // Integer division
    loc[1] = index % 3;
    index = index / 3;
    loc[2] = index;
}

unsigned char determine_subgrid_index(
    const Point &isoCrossingPoint,
    const Cube &cube,
    const UnifiedGrid &grid)
{
    // Compute relative position within cube
    const float dx = grid.spacing[0];
    const float dy = grid.spacing[1];
    const float dz = grid.spacing[2];
    const float min_x = grid.min_coord[0];
    const float min_y = grid.min_coord[1];
    const float min_z = grid.min_coord[2];
    float rel_x = (isoCrossingPoint.x() - (cube.indices[0] * dx + min_x)) / dx;
    float rel_y = (isoCrossingPoint.y() - (cube.indices[1] * dy + min_y)) / dy;
    float rel_z = (isoCrossingPoint.z() - (cube.indices[2] * dz + min_z)) / dz;

    // Clamp to [0, 1] range (in case of floating point errors)
    rel_x = std::max(0.0f, std::min(1.0f, static_cast<float>(rel_x)));
    rel_y = std::max(0.0f, std::min(1.0f, static_cast<float>(rel_y)));
    rel_z = std::max(0.0f, std::min(1.0f, static_cast<float>(rel_z)));

    // Map to subcell indices [0, 1, 2]
    int loc[3];
    loc[0] = static_cast<int>(rel_x * 3.0f);
    loc[1] = static_cast<int>(rel_y * 3.0f);
    loc[2] = static_cast<int>(rel_z * 3.0f);

    // Handle edge case: if exactly at 1.0, map to index 2
    if (loc[0] == 3) loc[0] = 2;
    if (loc[1] == 3) loc[1] = 2;
    if (loc[2] == 3) loc[2] = 2;

    // Encode to base-3 index
    unsigned char subgrid_index = loc[0] + 3 * loc[1] + 9 * loc[2];
    return subgrid_index;
}

// ============================================================================
// 3× REFINED GRID INDEX COMPUTATION
// ============================================================================

int linear_cell_index3x(int i, int j, int k, const UnifiedGrid &grid)
{
    // Linear index in 3× refined grid
    // Grid dimensions in 3× space: 3*(nx-1) × 3*(ny-1) × 3*(nz-1)
    const int nx = grid.num_cells[0];
    const int ny = grid.num_cells[1];
    const int nx1 = nx - 1;
    const int ny1 = ny - 1;
    return k * nx1 * ny1 * 9 + j * nx1 * 3 + i;
}

// ============================================================================
// CONFLICT DETECTION FOR SEP_ISOV_3
// ============================================================================

bool does_cell_conflict_with_selected_cubes(
    const Cube &cube,
    const std::unordered_map<int, Cube> &selected_indices,
    const UnifiedGrid &grid,
    int *indexA,
    int clearance)
{
    // Compute this cube's position in 3× grid
    int loc[3];
    compute_subgrid_loc(cube.isov_subgrid_index, loc);

    int grid3x_loc[3];
    grid3x_loc[0] = 3 * cube.indices[0] + loc[0];
    grid3x_loc[1] = 3 * cube.indices[1] + loc[1];
    grid3x_loc[2] = 3 * cube.indices[2] + loc[2];

    int myIndexA = linear_cell_index3x(grid3x_loc[0], grid3x_loc[1], grid3x_loc[2], grid);
    if (indexA != nullptr) {
        *indexA = myIndexA;
    }

    const int radius = std::max(1, clearance);
    const int nx = grid.num_cells[0];
    const int ny = grid.num_cells[1];
    const int nz = grid.num_cells[2];

    // Check neighbors in 3× grid within requested clearance
    for (int dk = -radius; dk <= radius; ++dk) {
        for (int dj = -radius; dj <= radius; ++dj) {
            for (int di = -radius; di <= radius; ++di) {
                // Skip self
                if (di == 0 && dj == 0 && dk == 0) continue;

                int neighbor_loc[3] = {
                    grid3x_loc[0] + di,
                    grid3x_loc[1] + dj,
                    grid3x_loc[2] + dk
                };

                // Check bounds in 3× grid
                if (neighbor_loc[0] < 0 || neighbor_loc[1] < 0 || neighbor_loc[2] < 0) continue;
                if (neighbor_loc[0] >= 3 * (nx - 1) ||
                    neighbor_loc[1] >= 3 * (ny - 1) ||
                    neighbor_loc[2] >= 3 * (nz - 1)) continue;

                int neighborIndex = linear_cell_index3x(neighbor_loc[0], neighbor_loc[1], neighbor_loc[2], grid);

                if (selected_indices.find(neighborIndex) != selected_indices.end()) {
                    return true;  // Conflict found
                }
            }
        }
    }

    return false;  // No conflict
}

// ============================================================================
// SEPARATION METHOD III: Subgrid-Based (3×3×3 Division and variants)
// ============================================================================

static std::vector<Cube> separate_active_cubes_III_with_clearance(
    std::vector<Cube> &activeCubes,
    const UnifiedGrid &grid,
    float isovalue,
    int clearance)
{
    const float dx = grid.spacing[0];
    const float dy = grid.spacing[1];
    const float dz = grid.spacing[2];
    const float min_x = grid.min_coord[0];
    const float min_y = grid.min_coord[1];
    const float min_z = grid.min_coord[2];

    // Compute accurate iso-crossing points and determine subgrid indices
    for (Cube &cube : activeCubes)
    {
        Point accurate_crossing = compute_iso_crossing_point_accurate(grid, cube.indices[0], cube.indices[1], cube.indices[2], isovalue);
        cube.accurateIsoCrossing = accurate_crossing;
        cube.isov_subgrid_index = determine_subgrid_index(accurate_crossing, cube, grid);
    }

    // Sort by distance to boundary (same as method I)
    std::sort(activeCubes.begin(), activeCubes.end(),
              [&grid](const Cube &a, const Cube &b) {
                  int dist_a = min_distance_to_boundary(a.indices[0], a.indices[1], a.indices[2], grid);
                  int dist_b = min_distance_to_boundary(b.indices[0], b.indices[1], b.indices[2], grid);
                  if (dist_a == dist_b)
                  {
                      // deterministic tie-break to keep behavior stable across platforms
                      if (a.indices[2] != b.indices[2])
                          return a.indices[2] < b.indices[2];
                      if (a.indices[1] != b.indices[1])
                          return a.indices[1] < b.indices[1];
                      return a.indices[0] < b.indices[0];
                  }
                  return dist_a < dist_b;
              });

    // Selection with 3× grid conflict detection
    std::unordered_map<int, Cube> selected_indices; // key: 3× grid index
    selected_indices.reserve(activeCubes.size());
    std::vector<Cube> out;
    out.reserve(activeCubes.size());

    for (Cube &cube : activeCubes)
    {
        int indexA;
        if (!does_cell_conflict_with_selected_cubes(cube, selected_indices, grid, &indexA, clearance))
        {
            // No conflict - select this cube

            // 1. Big cube center
/*                 cube.cubeCenter = Point(
                    (cube.i + 0.5f) * dx + min_x,
                    (cube.j + 0.5f) * dy + min_y,
                    (cube.k + 0.5f) * dz + min_z); */

            // 2. Small cube center in 3x3x3 subgrid
            //    Compute the center of the small cube containing the iso-crossing point
            int loc[3];
            compute_subgrid_loc(cube.isov_subgrid_index, loc);

            // Small cube center: big cube base + (subgrid_loc + 0.5) / 3.0 * cell_size
            cube.cubeCenter = Point(
                (cube.indices[0] + (loc[0] + 0.5f) / 3.0f) * dx + min_x,
                (cube.indices[1] + (loc[1] + 0.5f) / 3.0f) * dy + min_y,
                (cube.indices[2] + (loc[2] + 0.5f) / 3.0f) * dz + min_z);
           
            selected_indices.emplace(indexA, cube);
            out.push_back(cube);
        }
    }

    return out;
}

std::vector<Cube> separate_active_cubes_III(
    std::vector<Cube> &activeCubes,
    const UnifiedGrid &grid,
    float isovalue)
{
    return separate_active_cubes_III_with_clearance(activeCubes, grid, isovalue, 1);
}

std::vector<Cube> separate_active_cubes_III_wide(
    std::vector<Cube> &activeCubes,
    const UnifiedGrid &grid,
    float isovalue)
{
    return separate_active_cubes_III_with_clearance(activeCubes, grid, isovalue, 2);
}

std::vector<Cube> separate_active_cubes_III_exact_binary(
    std::vector<Cube> &activeCubes,
    const UnifiedGrid &grid,
    float isovalue)
{
    const float dx = grid.spacing[0];
    const float dy = grid.spacing[1];
    const float dz = grid.spacing[2];
    const float min_x = grid.min_coord[0];
    const float min_y = grid.min_coord[1];
    const float min_z = grid.min_coord[2];

    // Compute accurate iso-crossing points and determine subgrid indices
    for (Cube &cube : activeCubes)
    {
        Point accurate_crossing = compute_iso_crossing_point_accurate(grid, cube.indices[0], cube.indices[1], cube.indices[2], isovalue);
        cube.accurateIsoCrossing = accurate_crossing;
        cube.isov_subgrid_index = determine_subgrid_index(accurate_crossing, cube, grid);
    }

    // Sort by distance to boundary (same as method III)
    std::sort(activeCubes.begin(), activeCubes.end(),
              [&grid](const Cube &a, const Cube &b) {
                  int dist_a = min_distance_to_boundary(a.indices[0], a.indices[1], a.indices[2], grid);
                  int dist_b = min_distance_to_boundary(b.indices[0], b.indices[1], b.indices[2], grid);
                  if (dist_a == dist_b)
                  {
                      // Deterministic tie-break
                      if (a.indices[2] != b.indices[2]) return a.indices[2] < b.indices[2];
                      if (a.indices[1] != b.indices[1]) return a.indices[1] < b.indices[1];
                      return a.indices[0] < b.indices[0];
                  }
                  return dist_a < dist_b;
              });

    // Selection with 3× grid conflict detection
    std::unordered_map<int, Cube> selected_indices;
    selected_indices.reserve(activeCubes.size());
    std::vector<Cube> out;
    out.reserve(activeCubes.size());

    for (Cube &cube : activeCubes)
    {
        int indexA;
        if (!does_cell_conflict_with_selected_cubes(cube, selected_indices, grid, &indexA, 1))
        {
            // No conflict - select this cube
            int loc[3];
            compute_subgrid_loc(cube.isov_subgrid_index, loc);

            // KEY CHANGE: Use exact binary fractions instead of thirds
            // Map loc[d] ∈ {0, 1, 2} to {0.25, 0.5, 0.75}
            // These have exact binary representations: 1/4, 1/2, 3/4
            static const float exact_offsets[3] = {0.25f, 0.5f, 0.75f};

            cube.cubeCenter = Point(
                (cube.indices[0] + exact_offsets[loc[0]]) * dx + min_x,
                (cube.indices[1] + exact_offsets[loc[1]]) * dy + min_y,
                (cube.indices[2] + exact_offsets[loc[2]]) * dz + min_z);

            selected_indices.emplace(indexA, cube);
            out.push_back(cube);
        }
    }

    return out;
}
