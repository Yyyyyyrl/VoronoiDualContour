#include "core/vdc.h"
#include "core/vdc_timing.h"
#include "grid/vdc_sep_isov.h"
#include <algorithm>
#include <array>
#include <iomanip>
#include <iterator>
#include <limits>

namespace
{

struct SummaryStats
{
    std::size_t active_cubes = 0;
    std::size_t delaunay_vertices = 0;
    std::size_t delaunay_cells = 0;
    std::size_t voronoi_vertices = 0;
    std::size_t voronoi_edges = 0;
    std::size_t voronoi_facets = 0;
    std::size_t voronoi_cells = 0;
    std::size_t min_cell_vertices = 0;
    std::size_t max_cell_vertices = 0;
    double avg_cell_vertices = 0.0;
    std::size_t min_cell_facets = 0;
    std::size_t max_cell_facets = 0;
    double avg_cell_facets = 0.0;
    int min_cell_index = -1;
    int max_cell_index = -1;
    std::size_t min_facet_vertices = 0;
    std::size_t max_facet_vertices = 0;
    double avg_facet_vertices = 0.0;
    std::size_t min_facet_edges = 0;
    std::size_t max_facet_edges = 0;
    double avg_facet_edges = 0.0;
    int min_facet_index = -1;
    int max_facet_index = -1;
    std::array<std::size_t, 4> facet_match_counts{0, 0, 0, 0};
    std::size_t iso_vertices = 0;
    std::size_t iso_triangles = 0;
    bool multi_isov = false;
    std::size_t collapsed_vertices_removed = 0;
    std::size_t collapsed_edges_removed = 0;
    std::size_t mod_cyc_flips = 0;
    std::size_t mod_cyc_interior_flips = 0;
    std::size_t mod_cyc_boundary_flips = 0;
    std::size_t isovertex_clipped_count = 0;
    double isovertex_max_clip_distance = 0.0;
};

SummaryStats collect_summary_stats(const std::vector<Cube> &activeCubes,
                                   const Delaunay &dt,
                                   const VoronoiDiagram &vd,
                                   const IsoSurface &iso_surface,
                                   bool multi_isov,
                                   std::size_t collapsed_vertices_removed,
                                   std::size_t collapsed_edges_removed,
                                   std::size_t mod_cyc_flips,
                                   std::size_t mod_cyc_interior_flips,
                                   std::size_t mod_cyc_boundary_flips,
                                   std::size_t isovertex_clipped_count,
                                   double isovertex_max_clip_distance)
{
    SummaryStats stats;
    stats.multi_isov = multi_isov;
    stats.active_cubes = activeCubes.size();
    stats.delaunay_vertices = static_cast<std::size_t>(std::distance(dt.finite_vertices_begin(), dt.finite_vertices_end()));
    stats.delaunay_cells = static_cast<std::size_t>(std::distance(dt.finite_cells_begin(), dt.finite_cells_end()));
    stats.voronoi_vertices = vd.vertices.size();
    stats.voronoi_edges = vd.edges.size();
    stats.voronoi_facets = vd.global_facets.size();
    stats.voronoi_cells = vd.cells.size();
    stats.collapsed_vertices_removed = collapsed_vertices_removed;
    stats.collapsed_edges_removed = collapsed_edges_removed;
    stats.mod_cyc_flips = mod_cyc_flips;
    stats.mod_cyc_interior_flips = mod_cyc_interior_flips;
    stats.mod_cyc_boundary_flips = mod_cyc_boundary_flips;
    stats.isovertex_clipped_count = isovertex_clipped_count;
    stats.isovertex_max_clip_distance = isovertex_max_clip_distance;

    std::size_t total_cell_vertices = 0;
    std::size_t total_cell_facets = 0;
    std::size_t minCellVertices = std::numeric_limits<std::size_t>::max();
    std::size_t maxCellVertices = 0;
    std::size_t minCellFacets = std::numeric_limits<std::size_t>::max();
    std::size_t maxCellFacets = 0;
    int minCellIndex = -1;
    int maxCellIndex = -1;
    for (const auto &cell : vd.cells)
    {
        const std::size_t vertCount = cell.vertices_indices.size();
        const std::size_t facetCount = cell.facet_indices.size();
        total_cell_vertices += vertCount;
        total_cell_facets += facetCount;

        if (vertCount < minCellVertices)
        {
            minCellVertices = vertCount;
            minCellFacets = facetCount;
            minCellIndex = cell.cellIndex;
        }
        if (vertCount > maxCellVertices)
        {
            maxCellVertices = vertCount;
            maxCellFacets = facetCount;
            maxCellIndex = cell.cellIndex;
        }
    }
    if (!vd.cells.empty())
    {
        const double denom = static_cast<double>(vd.cells.size());
        stats.avg_cell_vertices = static_cast<double>(total_cell_vertices) / denom;
        stats.avg_cell_facets = static_cast<double>(total_cell_facets) / denom;
        stats.min_cell_vertices = minCellVertices;
        stats.max_cell_vertices = maxCellVertices;
        stats.min_cell_facets = minCellFacets;
        stats.max_cell_facets = maxCellFacets;
        stats.min_cell_index = minCellIndex;
        stats.max_cell_index = maxCellIndex;
    }

    std::size_t total_facet_vertices = 0;
    std::size_t total_facet_edges = 0;
    std::size_t minFacetVertices = std::numeric_limits<std::size_t>::max();
    std::size_t maxFacetVertices = 0;
    std::size_t minFacetEdges = std::numeric_limits<std::size_t>::max();
    std::size_t maxFacetEdges = 0;
    int minFacetIndex = -1;
    int maxFacetIndex = -1;
    for (const auto &gf : vd.global_facets)
    {
        const std::size_t vertCount = gf.vertices_indices.size();
        const std::size_t edgeCount = gf.voronoi_edge_indices.size();
        total_facet_vertices += vertCount;
        total_facet_edges += edgeCount;

        if (vertCount < minFacetVertices)
        {
            minFacetVertices = vertCount;
            minFacetEdges = edgeCount;
            minFacetIndex = gf.index;
        }
        if (vertCount > maxFacetVertices)
        {
            maxFacetVertices = vertCount;
            maxFacetEdges = edgeCount;
            maxFacetIndex = gf.index;
        }

        const auto methodIndex = static_cast<std::size_t>(gf.bipolar_match_method);
        if (methodIndex < stats.facet_match_counts.size())
            ++stats.facet_match_counts[methodIndex];
    }
    if (!vd.global_facets.empty())
    {
        const double denom = static_cast<double>(vd.global_facets.size());
        stats.avg_facet_vertices = static_cast<double>(total_facet_vertices) / denom;
        stats.avg_facet_edges = static_cast<double>(total_facet_edges) / denom;
        stats.min_facet_vertices = minFacetVertices;
        stats.max_facet_vertices = maxFacetVertices;
        stats.min_facet_edges = minFacetEdges;
        stats.max_facet_edges = maxFacetEdges;
        stats.min_facet_index = minFacetIndex;
        stats.max_facet_index = maxFacetIndex;
    }

    stats.iso_vertices = iso_surface.isosurfaceVertices.size();
    stats.iso_triangles = multi_isov ? iso_surface.isosurfaceTrianglesMulti.size()
                                     : iso_surface.isosurfaceTrianglesSingle.size();
    return stats;
}

void print_summary_report(const SummaryStats &stats)
{
    std::cout << "\n[SUMMARY] Run statistics\n";
    std::cout << "  Active cubes: " << stats.active_cubes << "\n";
    std::cout << "  Delaunay finite vertices: " << stats.delaunay_vertices
              << ", finite cells: " << stats.delaunay_cells << "\n";
    std::cout << "  Voronoi vertices: " << stats.voronoi_vertices
              << ", edges: " << stats.voronoi_edges
              << ", facets: " << stats.voronoi_facets
              << ", cells: " << stats.voronoi_cells << "\n";
    if (stats.collapsed_vertices_removed > 0 || stats.collapsed_edges_removed > 0)
    {
        std::cout << "    (collapseSmallEdges removed vertices=" << stats.collapsed_vertices_removed
                  << ", edges=" << stats.collapsed_edges_removed << ")\n";
    }
    if (stats.collapsed_vertices_removed > 0 || stats.collapsed_edges_removed > 0)
    {
        std::cout << "    (collapseSmallEdges removed vertices=" << stats.collapsed_vertices_removed
                  << ", edges=" << stats.collapsed_edges_removed << ")\n";
    }

    if (stats.max_cell_index != -1)
    {
        std::cout << "  Voronoi cell vertices (min / avg / max): "
                  << stats.min_cell_vertices << " / " << std::fixed << std::setprecision(2) << stats.avg_cell_vertices
                  << " / " << stats.max_cell_vertices << std::defaultfloat
                  << "  [indices: min=" << stats.min_cell_index
                  << ", max=" << stats.max_cell_index << "]\n";
        std::cout << "  Voronoi cell facets (min / avg / max): "
                  << stats.min_cell_facets << " / " << std::fixed << std::setprecision(2) << stats.avg_cell_facets
                  << " / " << stats.max_cell_facets << std::defaultfloat << "\n";
    }

    if (stats.max_facet_index != -1)
    {
        std::cout << "  Voronoi facet vertices (min / avg / max): "
                  << stats.min_facet_vertices << " / " << std::fixed << std::setprecision(2) << stats.avg_facet_vertices
                  << " / " << stats.max_facet_vertices << std::defaultfloat
                  << "  [indices: min=" << stats.min_facet_index
                  << ", max=" << stats.max_facet_index << "]\n";
        std::cout << "  Voronoi facet edges (min / avg / max): "
                  << stats.min_facet_edges << " / " << std::fixed << std::setprecision(2) << stats.avg_facet_edges
                  << " / " << stats.max_facet_edges << std::defaultfloat << "\n";
    }

    static const char *method_names[] = {
        "SEP_POS", "SEP_NEG", "UNCONSTRAINED_MATCH", "UNDEFINED_MATCH_TYPE"};
    bool header_printed = false;
    for (std::size_t i = 0; i < stats.facet_match_counts.size(); ++i)
    {
        if (stats.facet_match_counts[i] == 0)
            continue;
        if (!header_printed)
        {
            std::cout << "  Facet match methods:";
            header_printed = true;
        }
        std::cout << ' ' << method_names[i] << "=" << stats.facet_match_counts[i];
    }
    if (header_printed)
        std::cout << "\n";

    std::cout << "  Iso-surface vertices: " << stats.iso_vertices
              << ", triangles: " << stats.iso_triangles
              << (stats.multi_isov ? " (multi-isov)" : " (single-isov)") << "\n";
    if (stats.multi_isov)
    {
        std::cout << "  Modify-cycles flips: " << stats.mod_cyc_flips
                  << " (interior: " << stats.mod_cyc_interior_flips
                  << ", boundary: " << stats.mod_cyc_boundary_flips << ")\n";
        if (stats.isovertex_clipped_count > 0)
        {
            std::cout << "  Iso-vertex clipping: " << stats.isovertex_clipped_count << " vertices clipped"
                      << ", max distance: " << std::fixed << std::setprecision(6)
                      << stats.isovertex_max_clip_distance << std::defaultfloat << "\n";
        }
    }

    std::cout << std::defaultfloat;
}

} // namespace

int main(int argc, char *argv[])
{
    TimingStats& timer = TimingStats::getInstance();
    timer.startTimer("Total Processing");

    VoronoiDiagram vd; // Initialize an empty Voronoi diagram.
    VDC_PARAM vdc_param;
    IsoSurface iso_surface;
    Delaunay dt;
    Delaunay dt_test;

    // Parse command-line arguments to set program options and parameters.
    parse_arguments(argc, argv, vdc_param);

    // Load the NRRD data file into a grid structure.
    timer.startTimer("1. Load Data and Grid Formation", "Total Processing");
    UnifiedGrid data_grid = load_nrrd_data(vdc_param.file_path);
    timer.stopTimer("1. Load Data and Grid Formation");

    // Apply supersampling if requested.
    if (vdc_param.supersample)
    {
        timer.startTimer("1. Load Data and Grid Formation", "Total Processing");
        data_grid = supersample_grid(data_grid, vdc_param.supersample_r);
        timer.stopTimer("1. Load Data and Grid Formation");
        if (debug) {
        data_grid.print_grid();}
    }

    // Identify active cubes in the grid based on the given isovalue.
    timer.startTimer("2. Data Pre-processing", "Total Processing");
    timer.startTimer("Find active cubes", "2. Data Pre-processing");
    std::vector<Cube> activeCubes;
    find_active_cubes(data_grid, vdc_param.isovalue, activeCubes);
    timer.stopTimer("Find active cubes");

    // Separate active cubes to ensure non-adjacency if requested.
    if (vdc_param.sep_isov_1)
    {
        timer.startTimer("Separation", "2. Data Pre-processing");
        std::cout << "[INFO] Separation method I: Cube-level (26-connectivity)" << std::endl;
        std::cout << "  Original # of active cubes: " << activeCubes.size() << std::endl;
        activeCubes = separate_active_cubes_I(activeCubes, data_grid, vdc_param.isovalue);
        std::cout << "  After separation: " << activeCubes.size() << " cubes" << std::endl;
        timer.stopTimer("Separation");
    }
    else if (vdc_param.sep_isov_3)
    {
        timer.startTimer("Separation", "2. Data Pre-processing");
        std::cout << "[INFO] Separation method III: 3×3×3 subgrid-based separation" << std::endl;
        std::cout << "  Original # of active cubes: " << activeCubes.size() << std::endl;
        activeCubes = separate_active_cubes_III(activeCubes, data_grid, vdc_param.isovalue);
        std::cout << "  After separation: " << activeCubes.size() << " cubes" << std::endl;
        timer.stopTimer("Separation");
    }
    else if (vdc_param.sep_isov_3_wide)
    {
        timer.startTimer("Separation", "2. Data Pre-processing");
        std::cout << "[INFO] Separation method III-wide: 3×3×3 subgrid with 5×5×5 clearance" << std::endl;
        std::cout << "  Original # of active cubes: " << activeCubes.size() << std::endl;
        activeCubes = separate_active_cubes_III_wide(activeCubes, data_grid, vdc_param.isovalue);
        std::cout << "  After separation: " << activeCubes.size() << " cubes" << std::endl;
        timer.stopTimer("Separation");
    }
    else
    {
        // When no separation is used, compute accurate iso-crossing points for all active cubes
        // (separation methods already do this for kept cubes)
        timer.startTimer("Compute iso-crossing points", "2. Data Pre-processing");
        for (Cube &cube : activeCubes)
        {
            cube.accurateIsoCrossing = compute_iso_crossing_point_accurate(data_grid, cube.i, cube.j, cube.k, vdc_param.isovalue);
        }
        timer.stopTimer("Compute iso-crossing points");
    }

    // Create grid facets from the active cubes for further processing.
    timer.startTimer("Create grid facets", "2. Data Pre-processing");
    std::vector<std::vector<GRID_FACETS>> grid_facets = create_grid_facets(activeCubes);
    timer.stopTimer("Create grid facets");
    timer.stopTimer("2. Data Pre-processing");

    // Extract the centers of the active cubes.
    std::vector<Point> activeCubeCenters = get_cube_centers(activeCubes);
    std::vector<Point> activeCubeAccurateIsoCrossingPoints = get_cube_accurate_iso_crossing_points(activeCubes);

    std::cout << "[INFO] Number of active cube centers: " << activeCubeCenters.size() << std::endl;

    // Define the bounding box of the grid.
    Point p_min(0, 0, 0);
    Point p_max(data_grid.max_x, data_grid.max_y, data_grid.max_z);
    K::Iso_cuboid_3 bbox(p_min, p_max);

    if (debug) // Print the bounding box dimensions if debugging is enabled.
    {
        std::cout << "[DEBUG] Bounding box: ("
                  << bbox.min() << ") to ("
                  << bbox.max() << ")" << std::endl;
    }

    float cubeSideLength = data_grid.physical_dx; // Store the cube side length (equal to physical grid spacing).

    // Construct the Delaunay triangulation using the grid facets.
    if (indicator)
    {
        std::cout << "[INFO] Constructing Delaunay triangulation..." << std::endl;
    }
    timer.startTimer("3. Delaunay Triangulation Construction", "Total Processing");
    construct_delaunay_triangulation(dt, data_grid, grid_facets, vdc_param, activeCubeCenters);
    timer.stopTimer("3. Delaunay Triangulation Construction");

    //std::cout << dt << std::endl;
    // Construct the Voronoi diagram based on the Delaunay triangulation.
    if (indicator)
    {
        std::cout << "[INFO] Constructing Voronoi diagram..." << std::endl;
    }

    timer.startTimer("4. Voronoi Diagram Construction", "Total Processing");
    construct_voronoi_diagram(vd, vdc_param, data_grid, bbox, dt);
    timer.stopTimer("4. Voronoi Diagram Construction");
    // Collapse threshold: use CLI value if provided; otherwise scale to grid spacing (1% of min spacing)
    double collapse_eps = (vdc_param.collapse_eps > 0.0)
                              ? vdc_param.collapse_eps
                              : std::min({data_grid.physical_dx, data_grid.physical_dy, data_grid.physical_dz}) * 0.01;
    if (vdc_param.collapse_eps <= 0.0) {
        // Persist the resolved default so downstream stages and logs can see it.
        vdc_param.collapse_eps = collapse_eps;
    }
    std::vector<int> vertex_mapping;  // Maps old vertex indices to new after collapse
    timer.startTimer("5. Collapse Small Edges", "Total Processing");
    VoronoiDiagram vd2 = collapseSmallEdges(vd, collapse_eps, bbox, dt, vertex_mapping);

    // Re-validate and normalize facet orientations on the collapsed diagram
    timer.startTimer("Post-collapse facet validation", "5. Collapse Small Edges");
    validate_facet_orientations_and_normals(vd2);
    timer.stopTimer("Post-collapse facet validation");

    // Repopulate cell_edge_index arrays after collapse rebuilt the cellEdges
    if (vdc_param.multi_isov)
    {
        timer.startTimer("Repopulate cell edge indices", "5. Collapse Small Edges");
        populate_cell_edge_indices(vd2, dt);
        timer.stopTimer("Repopulate cell edge indices");
    }
    timer.stopTimer("5. Collapse Small Edges");

    timer.startTimer("6. Post-collapse Validation", "Total Processing");
    vd2.check(true);
    timer.stopTimer("6. Post-collapse Validation");

    //std::cout << dt <<std::endl;
    if (indicator)
    {
        std::cout << "[INFO] Constructing Iso Surface..." << std::endl;
    }

    const std::size_t collapsed_vertices_removed = (vd.vertices.size() > vd2.vertices.size())
                                                       ? (vd.vertices.size() - vd2.vertices.size())
                                                       : 0;
    const std::size_t collapsed_edges_removed = (vd.edges.size() > vd2.edges.size())
                                                    ? (vd.edges.size() - vd2.edges.size())
                                                    : 0;

    int interior_flips = 0, boundary_flips = 0;
    std::size_t clipped_count = 0;
    double max_clip_distance = 0.0;
    timer.startTimer("7. Isosurface Construction", "Total Processing");
    construct_iso_surface(dt, vd2, vdc_param, iso_surface, data_grid, activeCubeCenters, activeCubeAccurateIsoCrossingPoints, bbox, &vertex_mapping, &interior_flips, &boundary_flips, &clipped_count, &max_clip_distance);
    timer.stopTimer("7. Isosurface Construction");

    //write_voronoiDiagram(vd2, vdc_param.output_filename);

    // Handle the output mesh generation and return the appropriate status.
    timer.startTimer("8. Output Mesh", "Total Processing");
    bool retFlag;
    int retVal = handle_output_mesh(retFlag, vd2, vdc_param, iso_surface);
    timer.stopTimer("8. Output Mesh");

    if (retFlag)
        return retVal;

    if (vdc_param.summary_stats)
    {
        const int total_flips = interior_flips + boundary_flips;
        SummaryStats summary = collect_summary_stats(activeCubes, dt, vd2, iso_surface, vdc_param.multi_isov,
                                                    collapsed_vertices_removed, collapsed_edges_removed,
                                                    total_flips, interior_flips, boundary_flips,
                                                    clipped_count, max_clip_distance);
        print_summary_report(summary);
    }

    std::cout << "Finished." << std::endl;
    timer.stopTimer("Total Processing");

    // Print the comprehensive timing report
    timer.printReport();

    return EXIT_SUCCESS;
}
