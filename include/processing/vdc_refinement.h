//! @file vdc_refinement.h
//! @brief Facet-centric Delaunay refinement to improve small angles near the isosurface.

#ifndef VDC_REFINEMENT_H
#define VDC_REFINEMENT_H

#include "core/vdc_commandline.h"
#include "processing/vdc_grid.h"

//! @brief Statistics reported by the refinement pass.
struct RefinementStats
{
    int iterations_run = 0;            //!< Number of refinement iterations executed.
    std::size_t candidate_facets = 0;  //!< Facets with min angle below threshold.
    std::size_t bipolar_facets = 0;    //!< Those whose dual Voronoi edge is bipolar.
    std::size_t non_bipolar_facets = 0; //!< Candidates rejected because the dual edge is not bipolar.
    std::size_t inserted_points = 0;   //!< Points actually inserted into the Delaunay.
    std::size_t rejected_spacing = 0;  //!< Candidates discarded for being too close.
    std::size_t truncated_candidates = 0; //!< Candidates dropped due to per-iter cap.
};

//! @brief Candidate refinement point (position plus culprit angle).
struct RefinementCandidate
{
    Point position;      //!< Candidate insertion position.
    double min_angle_deg; //!< Minimum facet angle
};

//! @brief Compute the minimum internal angle (degrees) of a triangle.
double min_triangle_angle_deg(const Point &a, const Point &b, const Point &c);

//! @brief Check whether a Delaunay facet touches any dummy vertex.
bool facet_has_dummy_vertex(const Facet &facet);

//! @brief Retrieve the dual of a facet as a finite Voronoi segment if possible.
bool dual_as_segment(const Delaunay &dt, const Facet &facet, Segment3 &out_segment);

//! @brief Test whether a Voronoi segment is bipolar with respect to an isovalue.
/*!
 * Samples scalar values at both endpoints (clamped into the grid bounds) and
 * returns true if the values straddle \p isovalue.
 */
bool is_bipolar_segment(const Segment3 &seg, const UnifiedGrid &grid, float isovalue, float &out_v0, float &out_v1);

//! @brief Compute the midpoint of two points.
Point midpoint(const Point &a, const Point &b);

//! @brief Find the active cube whose center is closest to a given point.
const Cube *find_nearest_active_cube(const std::vector<Cube> &cubes, const Point &p);

//! @brief Snap a point to the center of a subcell (1×1×1, 2×2×2, or 3×3×3) within a cube.
Point snap_to_subcell_center(const Cube &cube, const UnifiedGrid &grid, const Point &target, int resolution);

//! @brief Interpolate along a segment to the given isovalue; falls back to midpoint if flat.
Point interpolate_to_isovalue(const Point &p0, const Point &p1, float v0, float v1, float isovalue);

//! @brief Check spacing from the current Delaunay vertices against a minimum squared distance.
bool is_far_from_delaunay(const Delaunay &dt, const Point &p, double min_spacing_sq);

//! @brief Check spacing from already enqueued refinement candidates.
bool is_far_from_candidates(const std::vector<RefinementCandidate> &candidates, const Point &p, double min_spacing_sq);

//! @brief Compute the spacing threshold to use for refinement (uses user value or 5% of min grid spacing).
double compute_min_spacing(const UnifiedGrid &grid, double user_value);

//! @brief Return the next available vertex index (max existing index + 1).
int next_vertex_index(const Delaunay &dt);

//! @brief Reassign sequential cell indices and clear cached dual vertex indices.
void reindex_cells(Delaunay &dt);

//! @brief Refine poorly shaped facets whose dual Voronoi edge is bipolar.
/*!
 * Scans finite Delaunay facets, checks their minimum corner angle, and if the
 * dual Voronoi edge is bipolar with respect to the target isovalue, inserts a
 * new vertex snapped to the nearest active cube (or directly on the dual
 * segment when not snapping). The pass is bounded by iteration and per-iteration
 * limits and skips any facet involving dummy vertices.
 *
 * @param dt Delaunay triangulation to refine (modified in place).
 * @param grid Scalar grid used for value interpolation and spacing metadata.
 * @param active_cubes Active cubes detected for the current isovalue.
 * @param params Global parameters controlling thresholds and budgets.
 * @param isovalue Target isovalue for bipolar testing.
 * @return Collected refinement statistics.
 */
RefinementStats refine_delaunay_small_angles(Delaunay &dt,
                                             const UnifiedGrid &grid,
                                             const std::vector<Cube> &active_cubes,
                                             const VdcParam &params,
                                             float isovalue);

#endif // VDC_REFINEMENT_H
