#include "core/vdc.h"
#include "core/vdc_stats.h"
#include "core/vdc_timing.h"
#include "processing/vdc_sep_isov.h"

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

    if (data_grid.boundary_crosses_isovalue(vdc_param.isovalue))
    {
        data_grid.zero_boundary_shell();
        std::cout << "[INFO] Clamped boundary voxels to minimum scalar value.\n";
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

    int interior_flips = 0, boundary_flips = 0, total_flips = 0;
    std::size_t clipped_count = 0;
    double max_clip_distance = 0.0;
    timer.startTimer("7. Isosurface Construction", "Total Processing");
    construct_iso_surface(dt, vd2, vdc_param, iso_surface, data_grid, activeCubeCenters, activeCubeAccurateIsoCrossingPoints, bbox, &vertex_mapping, &interior_flips, &boundary_flips, &total_flips, &clipped_count, &max_clip_distance);
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
