#include "vdc.h"

int main(int argc, char *argv[])
{
    VoronoiDiagram vd; // Initialize an empty Voronoi diagram.
    VDC_PARAM vdc_param;
    IsoSurface iso_surface;
    Delaunay dt;

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

    // Construct the Delaunay triangulation using the grid facets.
    if (indicator)
    {
        std::cout << "[INFO] Constructing Delaunay triangulation..." << std::endl;
    }
    construct_delaunay_triangulation(dt, data_grid, grid_facets, vdc_param, activeCubeCenters);

    //std::cout << dt << std::endl;
    // Construct the Voronoi diagram based on the Delaunay triangulation.
    if (indicator)
    {
        std::cout << "[INFO] Constructing Voronoi diagram..." << std::endl;
    }

    construct_voronoi_diagram(vd, vdc_param, data_grid, bbox, dt);
    if (vdc_param.test_vor) {
        // If test_vor is true means in testing mode for voronoi diagram construction, no need for further move
        return EXIT_SUCCESS;
    }


    if (indicator)
    {
        std::cout << "[INFO] Constructing Iso Surface..." << std::endl;
    }
    construct_iso_surface(dt, vd, vdc_param, iso_surface, data_grid, activeCubeCenters, bbox);

    write_voronoiDiagram(vd, vdc_param.output_filename);

    // Handle the output mesh generation and return the appropriate status.
    bool retFlag;
    int retVal = handle_output_mesh(retFlag, vd, vdc_param, iso_surface);
    if (retFlag)
        return retVal;

    std::cout << "Finished." << std::endl;

    return EXIT_SUCCESS;
}