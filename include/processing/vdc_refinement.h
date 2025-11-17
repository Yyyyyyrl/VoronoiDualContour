//! @file vdc_refinement.h
//! @brief Surface-aware Delaunay refinement to improve small angles in the final isosurface mesh.

#ifndef VDC_REFINEMENT_H
#define VDC_REFINEMENT_H

#include "core/vdc_type.h"
#include "processing/vdc_grid.h"

#include <unordered_set>
#include <vector>

//! @brief Parameters controlling facet-centric refinement near the isosurface.
struct SurfaceRefinementParams
{
    bool enable = false;                     //!< Master switch
    double max_radius_edge_ratio = 2.0;      //!< Max allowed radius-edge ratio
    double min_dihedral_deg = 10.0;          //!< Min allowed dihedral angle (deg)
    double min_surface_angle_deg = 20.0;     //!< Target min surface angle (deg) if using continuous placement
    int insert_resolution = 3;               //!< 1=cube center, 2=2x2x2, 3=3x3x3
    bool snap_to_grid = true;                //!< Snap to cube/subcell vs continuous iso point
    double min_spacing = 0.4;                //!< Min spacing to existing vertices (world units)
    int max_iterations = 1;                  //!< Max refinement iterations
    int max_new_points_per_iter = 1000;      //!< Cap points inserted per iteration
};

//! @brief Sparse mask for active cubes to support snapped placement.
struct ActiveMask
{
    std::unordered_set<long long> activeKeys;
    int nx = 0, ny = 0, nz = 0;

    bool is_active(int i, int j, int k) const;
    bool in_bounds(int i, int j, int k) const;
};

//! @brief Build an ActiveMask from the active cube list.
ActiveMask build_active_mask_from_cubes(const std::vector<Cube> &cubes, const UnifiedGrid &grid);

//! @brief Main entry: refine Delaunay facets dual to bipolar Voronoi edges to improve small angles.
void refine_surface_mesh_small_angles(Delaunay &dt,
                                      UnifiedGrid &grid,
                                      const ActiveMask &activeMask,
                                      const CGAL::Epick::Iso_cuboid_3 &bbox,
                                      float iso,
                                      const SurfaceRefinementParams &params);

#endif // VDC_REFINEMENT_H
