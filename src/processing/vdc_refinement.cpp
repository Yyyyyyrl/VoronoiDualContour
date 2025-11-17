#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <tuple>
#include <unordered_map>

#include <CGAL/Kernel/global_functions.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/intersections.h>

#include "processing/vdc_refinement.h"
#include "core/vdc_debug.h"
#include "core/vdc_utilities.h"
#include "processing/vdc_grid.h"
namespace
{
    using AngleMetrics = std::pair<double, double>;

    inline long long encode_key(int i, int j, int k)
    {
        return (static_cast<long long>(i) << 42) ^
               (static_cast<long long>(j) << 21) ^
               static_cast<long long>(k);
    }

    struct FacetRecord
    {
        Facet facet;
        double radius_edge_ratio = 0.0;
        double min_dihedral_deg = 180.0;
    };

    AngleMetrics evaluate_tet_quality(const Cell_handle &cell)
    {
        const Point &p0 = cell->vertex(0)->point();
        const Point &p1 = cell->vertex(1)->point();
        const Point &p2 = cell->vertex(2)->point();
        const Point &p3 = cell->vertex(3)->point();

        const double e01 = CGAL::to_double(CGAL::squared_distance(p0, p1));
        const double e02 = CGAL::to_double(CGAL::squared_distance(p0, p2));
        const double e03 = CGAL::to_double(CGAL::squared_distance(p0, p3));
        const double e12 = CGAL::to_double(CGAL::squared_distance(p1, p2));
        const double e13 = CGAL::to_double(CGAL::squared_distance(p1, p3));
        const double e23 = CGAL::to_double(CGAL::squared_distance(p2, p3));

        const double min_edge_sq = std::min({e01, e02, e03, e12, e13, e23});
        const double min_edge = std::sqrt(std::max(min_edge_sq, 1e-16));

        const Point c = CGAL::circumcenter(p0, p1, p2, p3);
        const double radius_sq = CGAL::to_double(CGAL::squared_distance(c, p0));
        const double radius = std::sqrt(std::max(radius_sq, 0.0));
        const double ratio = (min_edge > 0.0) ? (radius / min_edge) : std::numeric_limits<double>::infinity();

        auto unit_normal = [](const Point &a, const Point &b, const Point &c) -> Vector3
        {
            Vector3 n = CGAL::cross_product(b - a, c - a);
            const double len = std::sqrt(n.squared_length());
            if (len > 0.0)
            {
                n = n / len;
            }
            return n;
        };

        const Vector3 n0 = unit_normal(p1, p2, p3);
        const Vector3 n1 = unit_normal(p0, p3, p2);
        const Vector3 n2 = unit_normal(p0, p1, p3);
        const Vector3 n3 = unit_normal(p0, p2, p1);

        const std::array<Vector3, 4> normals = {n0, n1, n2, n3};
        constexpr double pi = 3.14159265358979323846;
        double min_dihedral = 180.0;
        for (int i = 0; i < 4; ++i)
        {
            for (int j = i + 1; j < 4; ++j)
            {
                const Vector3 &a = normals[i];
                const Vector3 &b = normals[j];
                const double dot = CGAL::to_double(a * b);
                const Vector3 cross = CGAL::cross_product(a, b);
                const double sin_theta = std::sqrt(CGAL::to_double(cross.squared_length()));
                const double angle = std::atan2(sin_theta, dot);
                const double deg = angle * 180.0 / pi;
                min_dihedral = std::min(min_dihedral, deg);
            }
        }

        return {ratio, min_dihedral};
    }

    bool dual_bipolar_segment(Delaunay &dt,
                              const Facet &facet,
                              UnifiedGrid &grid,
                              const CGAL::Epick::Iso_cuboid_3 &bbox,
                              float iso,
                              Segment3 &outSeg)
    {
        CGAL::Object dualObj = dt.dual(facet);

        Segment3 seg;
        Ray3 ray;
        Line3 line;

        if (CGAL::assign(seg, dualObj))
        {
            // already a finite segment
        }
        else if (CGAL::assign(ray, dualObj))
        {
            CGAL::Object inter = CGAL::intersection(bbox, ray);
            if (!CGAL::assign(seg, inter))
            {
                return false; // ray does not intersect bbox in a segment
            }
        }
        else if (CGAL::assign(line, dualObj))
        {
            CGAL::Object inter = CGAL::intersection(bbox, line);
            if (!CGAL::assign(seg, inter))
            {
                return false; // line does not intersect bbox in a segment
            }
        }
        else
        {
            return false; // unknown dual object
        }

        const Point s0 = adjust_outside_bound_points(seg.source(), grid, seg.source(), seg.target());
        const Point s1 = adjust_outside_bound_points(seg.target(), grid, seg.source(), seg.target());

        const double v0 = trilinear_interpolate(s0, grid);
        const double v1 = trilinear_interpolate(s1, grid);
        if (!is_bipolar(static_cast<float>(v0), static_cast<float>(v1), iso))
        {
            return false;
        }

        outSeg = seg;
        return true;
    }

    std::tuple<int, int, int> world_to_index(const Point &p, const UnifiedGrid &grid)
    {
        const double x = (p.x() - grid.min_coord[0]) / grid.spacing[0];
        const double y = (p.y() - grid.min_coord[1]) / grid.spacing[1];
        const double z = (p.z() - grid.min_coord[2]) / grid.spacing[2];
        return {static_cast<int>(std::floor(x)), static_cast<int>(std::floor(y)), static_cast<int>(std::floor(z))};
    }

    Point subcell_center_world(const UnifiedGrid &grid, int i, int j, int k, const Point &ref, int res)
    {
        const double hx = grid.spacing[0] / static_cast<double>(res);
        const double hy = grid.spacing[1] / static_cast<double>(res);
        const double hz = grid.spacing[2] / static_cast<double>(res);

        const double base_x = grid.min_coord[0] + (static_cast<double>(i) * grid.spacing[0]);
        const double base_y = grid.min_coord[1] + (static_cast<double>(j) * grid.spacing[1]);
        const double base_z = grid.min_coord[2] + (static_cast<double>(k) * grid.spacing[2]);

        const double lx = (ref.x() - base_x) / grid.spacing[0];
        const double ly = (ref.y() - base_y) / grid.spacing[1];
        const double lz = (ref.z() - base_z) / grid.spacing[2];

        const int sx = std::clamp(static_cast<int>(std::floor(lx * res)), 0, res - 1);
        const int sy = std::clamp(static_cast<int>(std::floor(ly * res)), 0, res - 1);
        const int sz = std::clamp(static_cast<int>(std::floor(lz * res)), 0, res - 1);

        const double cx = base_x + (sx + 0.5) * hx;
        const double cy = base_y + (sy + 0.5) * hy;
        const double cz = base_z + (sz + 0.5) * hz;
        return Point(cx, cy, cz);
    }

    struct ActiveCubeHit
    {
        Point center;
    };

    bool trace_active_cube_along_segment(const Segment3 &seg,
                                         UnifiedGrid &grid,
                                         const ActiveMask &activeMask,
                                         int res,
                                         ActiveCubeHit &hit)
    {
        const int steps = 128;
        for (int s = 0; s <= steps; ++s)
        {
            const double t = static_cast<double>(s) / static_cast<double>(steps);
            const Point p = seg.point(t);
            auto [i, j, k] = world_to_index(p, grid);
            if (!activeMask.in_bounds(i, j, k))
            {
                continue;
            }
            if (activeMask.is_active(i, j, k))
            {
                if (res <= 1)
                {
                    const double cx = grid.min_coord[0] + (static_cast<double>(i) + 0.5) * grid.spacing[0];
                    const double cy = grid.min_coord[1] + (static_cast<double>(j) + 0.5) * grid.spacing[1];
                    const double cz = grid.min_coord[2] + (static_cast<double>(k) + 0.5) * grid.spacing[2];
                    hit.center = Point(cx, cy, cz);
                }
                else
                {
                    hit.center = subcell_center_world(grid, i, j, k, p, res);
                }
                return true;
            }
        }
        return false;
    }

    Point bisection_on_segment_for_isovalue(const Segment3 &seg,
                                           UnifiedGrid &grid,
                                           float iso)
    {
        Point a = seg.source();
        Point b = seg.target();
        double va = trilinear_interpolate(a, grid) - iso;
        double vb = trilinear_interpolate(b, grid) - iso;

        for (int it = 0; it < 30; ++it)
        {
            const Point mid = CGAL::midpoint(a, b);
            const double vm = trilinear_interpolate(mid, grid) - iso;
            if (va * vm <= 0.0)
            {
                b = mid;
                vb = vm;
            }
            else
            {
                a = mid;
                va = vm;
            }
        }
        return CGAL::midpoint(a, b);
    }

    double nearest_vertex_distance(const Delaunay &dt, const Point &p)
    {
        if (dt.number_of_vertices() == 0)
        {
            return std::numeric_limits<double>::infinity();
        }
        Vertex_handle vh = dt.nearest_vertex(p);
        return std::sqrt(CGAL::to_double(CGAL::squared_distance(vh->point(), p)));
    }

    int next_vertex_index(const Delaunay &dt)
    {
        int max_idx = -1;
        for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit)
        {
            max_idx = std::max(max_idx, vit->info().index);
        }
        if (max_idx < 0)
        {
            max_idx = static_cast<int>(dt.number_of_vertices());
        }
        return max_idx + 1;
    }

    void batch_insert_points(Delaunay &dt, const std::vector<Point> &candidates)
    {
        std::unordered_map<long long, Point> seen;
        seen.reserve(candidates.size());
        std::vector<Point> dedup;
        dedup.reserve(candidates.size());

        const double cell = 1e-6;
        auto key_from_point = [cell](const Point &p) -> long long
        {
            const long long ix = static_cast<long long>(std::floor(p.x() / cell));
            const long long iy = static_cast<long long>(std::floor(p.y() / cell));
            const long long iz = static_cast<long long>(std::floor(p.z() / cell));
            return (ix << 42) ^ (iy << 21) ^ iz;
        };

        for (const Point &p : candidates)
        {
            const long long key = key_from_point(p);
            if (seen.find(key) == seen.end())
            {
                seen[key] = p;
                dedup.push_back(p);
            }
        }

        int next_idx = next_vertex_index(dt);
        for (const Point &p : dedup)
        {
            Vertex_handle vh = dt.insert(p);
            vh->info().index = next_idx++;
            vh->info().is_dummy = false;
            vh->info().voronoiCellIndex = -1;
        }
    }

    std::vector<FacetRecord> collect_bad_surface_facets(Delaunay &dt,
                                                        UnifiedGrid &grid,
                                                        const CGAL::Epick::Iso_cuboid_3 &bbox,
                                                        float iso,
                                                        const SurfaceRefinementParams &params)
    {
        std::vector<FacetRecord> bad;
        for (auto cit = dt.finite_cells_begin(); cit != dt.finite_cells_end(); ++cit)
        {
            const Cell_handle C = cit;
            auto [ratio, minDia] = evaluate_tet_quality(C);
            const bool cell_good = (ratio <= params.max_radius_edge_ratio) && (minDia >= params.min_dihedral_deg);

            for (int fi = 0; fi < 4; ++fi)
            {
                Facet f(C, fi);
                Segment3 seg;
                if (!dual_bipolar_segment(dt, f, grid, bbox, iso, seg))
                {
                    continue;
                }

                if (cell_good)
                {
                    continue;
                }

                FacetRecord rec;
                rec.facet = f;
                rec.radius_edge_ratio = ratio;
                rec.min_dihedral_deg = minDia;
                bad.push_back(std::move(rec));
            }
        }
        return bad;
    }

    bool choose_refinement_point(const Segment3 &seg,
                                 UnifiedGrid &grid,
                                 const ActiveMask &activeMask,
                                 const SurfaceRefinementParams &params,
                                 float iso,
                                 Point &outPoint)
    {
        if (params.snap_to_grid)
        {
            ActiveCubeHit hit;
            if (!trace_active_cube_along_segment(seg, grid, activeMask, params.insert_resolution, hit))
            {
                return false;
            }
            outPoint = hit.center;
            return true;
        }

        outPoint = bisection_on_segment_for_isovalue(seg, grid, iso);
        return true;
    }
} // namespace

bool ActiveMask::in_bounds(int i, int j, int k) const
{
    return (i >= 0 && j >= 0 && k >= 0 && i < nx && j < ny && k < nz);
}

bool ActiveMask::is_active(int i, int j, int k) const
{
    if (!in_bounds(i, j, k))
    {
        return false;
    }
    return activeKeys.find(encode_key(i, j, k)) != activeKeys.end();
}

ActiveMask build_active_mask_from_cubes(const std::vector<Cube> &cubes, const UnifiedGrid &grid)
{
    ActiveMask mask;
    mask.nx = grid.num_cells[0];
    mask.ny = grid.num_cells[1];
    mask.nz = grid.num_cells[2];
    mask.activeKeys.reserve(cubes.size() * 2);
    for (const Cube &c : cubes)
    {
        mask.activeKeys.insert(encode_key(c.indices[0], c.indices[1], c.indices[2]));
    }
    return mask;
}

void refine_surface_mesh_small_angles(Delaunay &dt,
                                      UnifiedGrid &grid,
                                      const ActiveMask &activeMask,
                                      const CGAL::Epick::Iso_cuboid_3 &bbox,
                                      float iso,
                                      const SurfaceRefinementParams &params)
{
    if (!params.enable)
    {
        return;
    }

    int iter = 0;
    while (iter < params.max_iterations)
    {
        std::vector<FacetRecord> bad = collect_bad_surface_facets(dt, grid, bbox, iso, params);
        if (bad.empty())
        {
            break;
        }

        std::vector<Point> inserts;
        inserts.reserve(bad.size());

        for (const FacetRecord &rec : bad)
        {
            if (static_cast<int>(inserts.size()) >= params.max_new_points_per_iter)
            {
                break;
            }

            Segment3 seg;
            if (!dual_bipolar_segment(dt, rec.facet, grid, bbox, iso, seg))
            {
                continue;
            }

            Point p_ref;
            if (!choose_refinement_point(seg, grid, activeMask, params, iso, p_ref))
            {
                continue;
            }

            if (nearest_vertex_distance(dt, p_ref) < params.min_spacing)
            {
                continue;
            }

            inserts.push_back(p_ref);
        }

        if (inserts.empty())
        {
            break;
        }

        batch_insert_points(dt, inserts);
        ++iter;
    }
}
