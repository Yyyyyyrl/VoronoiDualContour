#include "vdc.h"

int main(int argc, char *argv[])
{
    std::clock_t initial_time = std::clock();
    VoronoiDiagram vd; // Initialize an empty Voronoi diagram.
    VDC_PARAM vdc_param;
    IsoSurface iso_surface;
    Delaunay dt;
    Delaunay dt_test;

    std::clock_t start_time = std::clock();

    // Parse command-line arguments to set program options and parameters.
    parse_arguments(argc, argv, vdc_param);

    // Load the NRRD data file into a grid structure.
    UnifiedGrid data_grid = load_nrrd_data(vdc_param.file_path);

    // Apply supersampling if requested.
    if (vdc_param.supersample)
    {
        data_grid = supersample_grid(data_grid, vdc_param.supersample_r);
        if (debug) // Print the supersampled grid if debugging is enabled.
        {
            data_grid.print_grid();
        }
    }

    // Identify active cubes in the grid based on the given isovalue.
    std::vector<Cube> activeCubes;
    find_active_cubes(data_grid, vdc_param.isovalue, activeCubes);

    // Separate active cubes to ensure non-adjacency if requested.
    if (vdc_param.sep_isov)
    {
        activeCubes = separate_active_cubes_greedy(activeCubes, data_grid);
    }

    // Create grid facets from the active cubes for further processing.
    std::vector<std::vector<GRID_FACETS>> grid_facets = create_grid_facets(activeCubes);

    // Extract the centers of the active cubes.
    std::vector<Point> activeCubeCenters = get_cube_centers(activeCubes);
    std::clock_t test1 = std::clock();
    dt_test.insert(activeCubeCenters.begin(), activeCubeCenters.end());
    std::clock_t test2 = std::clock();
    double d = static_cast<double>(test2 - test1) / CLOCKS_PER_SEC;
    std::cout << "[INFO] Delaunay test time: " << d << " seconds." << std::endl;

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

    float cubeSideLength = data_grid.dx; // Store the cube side length (equal to grid spacing, assuming regular grid so dx/dy/dz will be equal).

    std::clock_t load_data_time = std::clock(); // Get ending clock ticks
    double duration_ld = static_cast<double>(load_data_time - start_time) / CLOCKS_PER_SEC;

    std::cout << "[INFO] Loading Data processing time: " << std::to_string(duration_ld) << " seconds." << std::endl;

    // Construct the Delaunay triangulation using the grid facets.
    if (indicator)
    {
        std::cout << "[INFO] Constructing Delaunay triangulation..." << std::endl;
    }
    construct_delaunay_triangulation(dt, data_grid, grid_facets, vdc_param, activeCubeCenters);

    std::clock_t construct_dt_time = std::clock(); // Get ending clock ticks
    double duration_dt = static_cast<double>(construct_dt_time - load_data_time) / CLOCKS_PER_SEC;
    std::cout << "[INFO] Constructing Delaunay triangulation time: " << std::to_string(duration_dt) << " seconds." << std::endl;

    //std::cout << dt << std::endl;
    // Construct the Voronoi diagram based on the Delaunay triangulation.
    if (indicator)
    {
        std::cout << "[INFO] Constructing Voronoi diagram..." << std::endl;
    }

    construct_voronoi_diagram(vd, vdc_param, data_grid, bbox, dt);

    std::clock_t cons_vd_time = std::clock();
    VoronoiDiagram vd2 = collapseSmallEdges(vd, 0.01, bbox, dt);
    std::clock_t collapse_time = std::clock();
    double duration_col = static_cast<double>(collapse_time - cons_vd_time) / CLOCKS_PER_SEC;
    std::cout << "[INFO] Collapsing time: " << std::to_string(duration_col) << " seconds." << std::endl;
    //vd2.create_global_facets();
    //vd2.compute_bipolar_matches(vdc_param.isovalue);
    
    vd2.check(true);
    std::clock_t check2_time = std::clock();
    double duration_vd2check = static_cast<double>(check2_time - collapse_time) / CLOCKS_PER_SEC;
    std::cout << "[INFO] Checking vd2 time: " << std::to_string(duration_vd2check) << " seconds." << std::endl;

    if (indicator)
    {
        std::cout << "[INFO] Constructing Iso Surface..." << std::endl;
    }
    construct_iso_surface(dt, vd2, vdc_param, iso_surface, data_grid, activeCubeCenters, bbox);

    std::clock_t construct_iso_time = std::clock(); // Get ending clock ticks
    double duration_iso = static_cast<double>(construct_iso_time - check2_time) / CLOCKS_PER_SEC;
    std::cout << "[INFO] Constructing Iso Surface time: " << std::to_string(duration_iso) << " seconds." << std::endl;

    write_voronoiDiagram(vd2, vdc_param.output_filename);

    // Handle the output mesh generation and return the appropriate status.
    bool retFlag;
    int retVal = handle_output_mesh(retFlag, vd2, vdc_param, iso_surface);
    if (retFlag)
        return retVal;

    std::cout << "Finished." << std::endl;
    std::clock_t finish_time = std::clock();
    double duration = static_cast<double>(finish_time - initial_time) / CLOCKS_PER_SEC;
    std::cout << "[INFO] Total processing time: " << std::to_string(duration) << " seconds." << std::endl;

    return EXIT_SUCCESS;
}