//! @file vdc_refinement.cpp
//! @brief Implementation of facet-centric Delaunay refinement.

#include <algorithm>
#include <cmath>
#include <limits>
#include <CGAL/squared_distance_3.h>

#include "processing/vdc_refinement.h"
#include "core/vdc_utilities.h"

constexpr double kPi = 3.14159265358979323846;

double min_triangle_angle_deg(const Point &a, const Point &b, const Point &c)
{
    auto angle_at = [](const Point &p, const Point &q1, const Point &q2) -> double
    {
        Vector3 v1 = q1 - p;
        Vector3 v2 = q2 - p;
        const double len1 = std::sqrt(v1.squared_length());
        const double len2 = std::sqrt(v2.squared_length());
        if (len1 <= std::numeric_limits<double>::epsilon() || len2 <= std::numeric_limits<double>::epsilon())
        {
            return 0.0;
        }
        double cos_theta = (v1 * v2) / (len1 * len2);
        cos_theta = std::max(-1.0, std::min(1.0, cos_theta));
        return std::acos(cos_theta) * 180.0 / kPi;
    };

    double a0 = angle_at(a, b, c);
    double a1 = angle_at(b, a, c);
    double a2 = angle_at(c, a, b);
    return std::min({a0, a1, a2});
}

bool facet_has_dummy_vertex(const Facet &facet)
{
    const Cell_handle &cell = facet.first;
    const int facet_index = facet.second;
    for (int k = 0; k < 3; ++k)
    {
        const int cell_vertex_index = CellInfo::FacetVertexIndex(facet_index, k);
        if (cell->vertex(cell_vertex_index)->info().is_dummy)
        {
            return true;
        }
    }
    return false;
}

bool dual_as_segment(const Delaunay &dt,
                     const Facet &facet,
                     Segment3 &out_segment)
{
    Object dual_obj = dt.dual(facet);
    return CGAL::assign(out_segment, dual_obj);
}

bool is_bipolar_segment(const Segment3 &seg,
                        const UnifiedGrid &grid,
                        float isovalue,
                        float &out_v0,
                        float &out_v1)
{
    const Point &p0 = seg.source();
    const Point &p1 = seg.target();
    out_v0 = trilinear_interpolate(adjust_outside_bound_points(p0, grid, p0, p1), grid);
    out_v1 = trilinear_interpolate(adjust_outside_bound_points(p1, grid, p0, p1), grid);
    return is_bipolar(out_v0, out_v1, isovalue);
}

Point midpoint(const Point &a, const Point &b)
{
    return Point((a.x() + b.x()) * 0.5, (a.y() + b.y()) * 0.5, (a.z() + b.z()) * 0.5);
}

const Cube *find_nearest_active_cube(const std::vector<Cube> &cubes, const Point &p)
{
    const Cube *best = nullptr;
    double best_dist_sq = std::numeric_limits<double>::max();
    for (const Cube &cube : cubes)
    {
        const double dist_sq = CGAL::squared_distance(cube.cubeCenter, p);
        if (dist_sq < best_dist_sq)
        {
            best_dist_sq = dist_sq;
            best = &cube;
        }
    }
    return best;
}

Point snap_to_subcell_center(const Cube &cube,
                             const UnifiedGrid &grid,
                             const Point &target,
                             int resolution)
{
    const int res = std::max(1, std::min(3, resolution));
    const double dx = grid.spacing[0];
    const double dy = grid.spacing[1];
    const double dz = grid.spacing[2];
    const double base_x = cube.indices[0] * dx + grid.min_coord[0];
    const double base_y = cube.indices[1] * dy + grid.min_coord[1];
    const double base_z = cube.indices[2] * dz + grid.min_coord[2];

    const auto coord_center = [&](double base, double d, double value) -> double
    {
        const double local = (value - base) / d;
        const double scaled = local * res;
        int idx = static_cast<int>(std::floor(scaled));
        idx = std::max(0, std::min(res - 1, idx));
        const double cell_size = d / static_cast<double>(res);
        return base + (static_cast<double>(idx) + 0.5) * cell_size;
    };

    return Point(coord_center(base_x, dx, target.x()),
                 coord_center(base_y, dy, target.y()),
                 coord_center(base_z, dz, target.z()));
}

Point interpolate_to_isovalue(const Point &p0,
                              const Point &p1,
                              float v0,
                              float v1,
                              float isovalue)
{
    const double denom = static_cast<double>(v1) - static_cast<double>(v0);
    if (std::abs(denom) < 1e-12)
    {
        return midpoint(p0, p1);
    }
    double t = (static_cast<double>(isovalue) - static_cast<double>(v0)) / denom;
    t = std::max(0.0, std::min(1.0, t));
    return Point(p0.x() + (p1.x() - p0.x()) * t,
                 p0.y() + (p1.y() - p0.y()) * t,
                 p0.z() + (p1.z() - p0.z()) * t);
}

bool is_far_from_delaunay(const Delaunay &dt,
                          const Point &p,
                          double min_spacing_sq)
{
    if (dt.number_of_vertices() == 0)
    {
        return true;
    }
    Vertex_handle vh = dt.nearest_vertex(p);
    const double dist_sq = CGAL::squared_distance(vh->point(), p);
    return dist_sq >= min_spacing_sq;
}

bool is_far_from_candidates(const std::vector<RefinementCandidate> &candidates,
                            const Point &p,
                            double min_spacing_sq)
{
    for (const auto &cand : candidates)
    {
        if (CGAL::squared_distance(cand.position, p) < min_spacing_sq)
        {
            return false;
        }
    }
    return true;
}

double compute_min_spacing(const UnifiedGrid &grid, double user_value)
{
    if (user_value > 0.0)
    {
        return user_value;
    }
    const double min_spacing = std::min({static_cast<double>(grid.spacing[0]),
                                         static_cast<double>(grid.spacing[1]),
                                         static_cast<double>(grid.spacing[2])});
    return 0.05 * min_spacing;
}

int next_vertex_index(const Delaunay &dt)
{
    int next_index = 0;
    for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit)
    {
        next_index = std::max(next_index, vit->info().index + 1);
    }
    return next_index;
}

void reindex_cells(Delaunay &dt)
{
    int cell_index = 0;
    for (auto cit = dt.finite_cells_begin(); cit != dt.finite_cells_end(); ++cit)
    {
        cit->info().index = cell_index++;
        cit->info().dualVoronoiVertexIndex = -1;
    }
}

RefinementStats refine_delaunay_small_angles(Delaunay &dt,
                                             const UnifiedGrid &grid,
                                             const std::vector<Cube> &active_cubes,
                                             const VdcParam &params,
                                             float isovalue)
{
    RefinementStats stats;
    if (!params.refine_small_angles || active_cubes.empty())
    {
        return stats;
    }

    const double angle_threshold = params.refine_min_surface_angle_deg;
    const int max_iterations = std::max(1, params.refine_max_iterations);
    const int insert_resolution = std::max(1, std::min(3, params.refine_insert_resolution));
    const double min_spacing = compute_min_spacing(grid, params.refine_min_spacing);
    const double min_spacing_sq = min_spacing * min_spacing;

    int vertex_index = next_vertex_index(dt);

    for (int iter = 0; iter < max_iterations; ++iter)
    {
        std::vector<RefinementCandidate> candidates;
        candidates.reserve(256);

        for (auto fit = dt.finite_facets_begin(); fit != dt.finite_facets_end(); ++fit)
        {
            const Facet &facet = *fit;
            if (facet_has_dummy_vertex(facet))
            {
                continue;
            }

            const Cell_handle &cell = facet.first;
            const int facet_index = facet.second;
            const Point p0 = cell->vertex(CellInfo::FacetVertexIndex(facet_index, 0))->point();
            const Point p1 = cell->vertex(CellInfo::FacetVertexIndex(facet_index, 1))->point();
            const Point p2 = cell->vertex(CellInfo::FacetVertexIndex(facet_index, 2))->point();

            const double min_angle = min_triangle_angle_deg(p0, p1, p2);
            if (min_angle >= angle_threshold)
            {
                continue;
            }
            ++stats.candidate_facets;

            Segment3 dual_segment;
            if (!dual_as_segment(dt, facet, dual_segment))
            {
                continue;
            }

            float v0 = 0.0f, v1 = 0.0f;
            if (!is_bipolar_segment(dual_segment, grid, isovalue, v0, v1))
            {
                ++stats.non_bipolar_facets;
                continue;
            }
            ++stats.bipolar_facets;

            const Point centroid = midpoint(dual_segment.source(), dual_segment.target());
            const Cube *nearest_cube = find_nearest_active_cube(active_cubes, centroid);
            if (nearest_cube == nullptr)
            {
                continue;
            }

            Point candidate_point = params.refine_snap_to_grid
                                        ? snap_to_subcell_center(*nearest_cube, grid, centroid, insert_resolution)
                                        : interpolate_to_isovalue(dual_segment.source(), dual_segment.target(), v0, v1, isovalue);

            if (!is_far_from_delaunay(dt, candidate_point, min_spacing_sq) ||
                !is_far_from_candidates(candidates, candidate_point, min_spacing_sq))
            {
                ++stats.rejected_spacing;
                continue;
            }

            candidates.push_back({candidate_point, min_angle});
        }

        if (candidates.empty())
        {
            break;
        }

        std::sort(candidates.begin(), candidates.end(),
                  [](const RefinementCandidate &a, const RefinementCandidate &b)
                  {
                      return a.min_angle_deg < b.min_angle_deg;
                  });


        for (const auto &cand : candidates)
        {
            Vertex_handle vh = dt.insert(cand.position);
            vh->info().index = vertex_index++;
            vh->info().is_dummy = false;
            vh->info().voronoiCellIndex = -1;
            ++stats.inserted_points;
        }

        ++stats.iterations_run;
    }

    reindex_cells(dt);
    return stats;
}
