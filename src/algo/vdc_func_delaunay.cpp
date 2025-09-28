#include "algo/vdc_func.h"

//! @brief Adds dummy points from a facet for Voronoi diagram bounding.
std::vector<Point> add_dummy_from_facet(const GRID_FACETS &facet,
                                        const UnifiedGrid &data_grid,
                                        double supersample_multiplier)
{
    std::vector<Point> points;

    // 2D slice dimension
    int dim0 = facet.axis_size[0];
    int dim1 = facet.axis_size[1];

    // For convenience
    int d = facet.orth_dir;
    int d1 = facet.axis_dir[0];
    int d2 = facet.axis_dir[1];

    // and localSize[] = (maxIndex[i] - minIndex[i] + 1)
    // The grid spacing in each dimension
    double dx[3] = {data_grid.dx, data_grid.dy, data_grid.dz};

    // Loop over the 2D slice
    for (int coord1 = 0; coord1 < dim1; coord1++)
    {
        for (int coord0 = 0; coord0 < dim0; coord0++)
        {
            if (!facet.CubeFlag(coord0, coord1))
                continue;

            // localX[d1] = coord0, localX[d2] = coord1
            int localX[3] = {0, 0, 0};
            localX[d1] = coord0;
            localX[d2] = coord1;

            // side=0 => localX[d] = 0, side=1 => localX[d] = localSize[d]-1
            localX[d] = (facet.side == 0) ? 0 : (facet.localSize[d] - 1);

            // Convert localX -> global indices
            int g[3];
            for (int i = 0; i < 3; i++)
            {
                g[i] = localX[i] + facet.minIndex[i];
            }

            // Compute center in real-world coordinates
            double cx = (g[0] + 0.5) * dx[0];
            double cy = (g[1] + 0.5) * dx[1];
            double cz = (g[2] + 0.5) * dx[2];

            // Offset by +/- dx[d]
            // The offset multiplier is for avoid bipolar edges touching dummy voronoi cells
            const double offsetMultiplier = supersample_multiplier > 0.0 ? supersample_multiplier : 1.0;
            double offset = ((facet.side == 0) ? -4 * dx[d] : 4 * dx[d]) * offsetMultiplier;
            if (d == 0)
                cx += offset;
            else if (d == 1)
                cy += offset;
            else
                cz += offset;

            points.emplace_back(cx, cy, cz);
        }
    }

    return points;
}

//! @brief Collects points for the Delaunay triangulation.
/*!
 * Gathers original points and dummy points from grid facets for multi-isovertex mode,
 * or uses only active cube centers for single-isovertex mode.
 *
 * @param grid The grid containing data.
 * @param grid_facets The grid facets for dummy point generation.
 * @param activeCubeCenters The list of center points of active cubes.
 * @param vdc_param The VDC_PARAM instance containing user input options.
 * @param delaunay_points Output vector for all points (original + dummy).
 */
static int collect_delaunay_points(UnifiedGrid &grid,
                                   const std::vector<std::vector<GRID_FACETS>> &grid_facets,
                                   const std::vector<Point> &activeCubeCenters,
                                   VDC_PARAM &vdc_param,
                                   std::vector<Point> &delaunay_points)
{
    delaunay_points = activeCubeCenters;
    int first_dummy_index = delaunay_points.size(); // Dummies start here

    const double supersampleMultiplier = vdc_param.supersample
                                             ? static_cast<double>(vdc_param.supersample_r)
                                             : 1.0;

    if (vdc_param.multi_isov)
    {
        for (int d = 0; d < 3; ++d)
        { // Assuming 3 dimensions
            for (const auto &f : grid_facets[d])
            {
                auto pointsf = add_dummy_from_facet(f, grid, supersampleMultiplier);
                delaunay_points.insert(delaunay_points.end(), pointsf.begin(), pointsf.end());
            }
        }
    }
    return first_dummy_index;
}

//! @brief Inserts points into the Delaunay triangulation.
/*!
 * Inserts the collected points into the triangulation and writes debug output if enabled.
 *
 * @param dt The Delaunay triangulation to insert points into.
 * @param delaunay_points The points to insert.
 * @param activeCubeCenters The list of center points of active cubes.
 * @param vdc_param The VDC_PARAM instance containing user input options.
 */
static Vertex_handle insert_point_into_delaunay_triangulation(Delaunay &dt,
                                                              const Point &p,
                                                              int index,
                                                              bool is_dummy)
{
    // Insert point and retrieve handle
    Vertex_handle vh = dt.insert(p);
    // Immediately assign index and dummy status
    vh->info().index = index;
    vh->info().is_dummy = is_dummy;
    return vh;
}

//! @brief Constructs a Delaunay triangulation from a grid and grid facets.
/*!
 * Constructs a 3D Delaunay triangulation using the grid's scalar values
 * and the facets of the active cubes.
 *
 * @param dt The Delaunay triangulation instance.
 * @param grid The grid containing scalar values.
 * @param grid_facets The grid facets to use in constructing the triangulation.
 * @param vdc_param The VDC_PARAM instance holding user input options.
 * @param activeCubeCenters The list of center points of active cubes.
 */
void construct_delaunay_triangulation(Delaunay &dt, UnifiedGrid &grid, const std::vector<std::vector<GRID_FACETS>> &grid_facets, VDC_PARAM &vdc_param, std::vector<Point> &activeCubeCenters)
{
    std::clock_t start = std::clock();

    std::vector<Point> delaunay_points;
    size_t first_dummy_index = collect_delaunay_points(grid, grid_facets, activeCubeCenters, vdc_param, delaunay_points);

    std::clock_t after_collect = std::clock();
    double collect_time = static_cast<double>(after_collect - start) / CLOCKS_PER_SEC;
    std::cout << "[INFO] Time to build point list: " << collect_time << " seconds" << std::endl;

    dt.clear();

    // Batch insert all points
    dt.insert(delaunay_points.begin(), delaunay_points.end());

    std::clock_t after_insert = std::clock();
    double insert_time = static_cast<double>(after_insert - after_collect) / CLOCKS_PER_SEC;
    std::cout << "[INFO] Vertices insert time: " << insert_time << " seconds" << std::endl;

    // Create map from point to original index
    std::map<Point, size_t> point_to_index;
    for (size_t i = 0; i < delaunay_points.size(); ++i)
    {
        point_to_index[delaunay_points[i]] = i;
    }

    // Assign info to vertices
    for (auto vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit)
    {
        const Point &p = vit->point();
        auto it = point_to_index.find(p);
        if (it == point_to_index.end())
        {
            std::cerr << "[ERROR] Vertex point not found in original points!" << std::endl;
            continue;
        }
        size_t original_index = it->second;
        vit->info().index = original_index;
        vit->info().is_dummy = (original_index >= first_dummy_index);
        vit->info().voronoiCellIndex = -1; // Initialize if needed
    }

    std::clock_t after_assign = std::clock();
    double assign_time = static_cast<double>(after_assign - after_insert) / CLOCKS_PER_SEC;
    std::cout << "[INFO] Time to assign vertex info: " << assign_time << " seconds" << std::endl;
}
