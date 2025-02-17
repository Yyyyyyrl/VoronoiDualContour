#include "vdc.h"

int main(int argc, char *argv[])
{
    VoronoiDiagram vd; // Initialize an empty Voronoi diagram.
    VDC_PARAM vdc_param;
    IsoSurface iso_surface;

    // Parse command-line arguments to set program options and parameters.
    parse_arguments(argc, argv, vdc_param);

    // Load the NRRD data file into a grid structure.
    Grid data_grid = load_nrrd_data(vdc_param.file_path);

    // Apply supersampling if requested.
    if (vdc_param.supersample)
    {
        data_grid = supersample_grid(data_grid, vdc_param.supersample_r);
        if (debug) // Print the supersampled grid if debugging is enabled.
        {
            data_grid.print_grid();
        }
    }

    // Initialize a ScalarGrid using the dimensions and spacing from the data grid.
    ScalarGrid grid(data_grid.nx, data_grid.ny, data_grid.nz, data_grid.dx, data_grid.dy, data_grid.dz, 0.0, 0.0, 0.0);
    initialize_scalar_grid(grid, data_grid); // Populate the scalar grid with data.

    // Identify active cubes in the grid based on the given isovalue.
    std::vector<Cube> activeCubes;
    find_active_cubes(data_grid, vdc_param.isovalue, activeCubes);

    // Separate active cubes to ensure non-adjacency if requested.
    if (vdc_param.sep_isov)
    {
        activeCubes = separate_active_cubes_greedy(activeCubes, data_grid.nx, data_grid.ny, data_grid.nz);
    }

    // Create grid facets from the active cubes for further processing.
    std::vector<std::vector<GRID_FACETS>> grid_facets = create_grid_facets(activeCubes);

    // Extract the centers of the active cubes.
    std::vector<Point> activeCubeCenters = get_cube_centers(activeCubes);

    // Load all grid points into a vector for later use.
    std::vector<Point> gridPoints = load_grid_points(data_grid);

    // Define the bounding box of the grid.
    Point p_min(0, 0, 0);
    Point p_max(data_grid.nx - 1, data_grid.ny - 1, data_grid.nz - 1);
    K::Iso_cuboid_3 bbox(p_min, p_max);

    if (debug) // Print the bounding box dimensions if debugging is enabled.
    {
        std::cout << "Bounding box: ("
                  << bbox.min() << ") to ("
                  << bbox.max() << ")" << std::endl;
    }

    float cubeSideLength = data_grid.dx; // Store the cube side length (equal to grid spacing, assuming regular grid so dx/dy/dz will be equal).

    // Construct the Delaunay triangulation using the grid facets.
    if (indicator)
    {
        std::cout << "Constructing Delaunay triangulation..." << std::endl;
    }
    std::map<Point, int> point_index_map; // Used in Single Iso-V Case ONLY
    construct_delaunay_triangulation(data_grid, grid_facets, vdc_param, activeCubeCenters, point_index_map);

    // Construct the Voronoi diagram based on the Delaunay triangulation.
    if (indicator)
    {
        std::cout << "Constructing Voronoi diagram..." << std::endl;
    }

    std::map<Object, std::vector<Facet>, ObjectComparator> voronoi_edge_to_delaunay_facet_map;
    std::map<Point, float> vertexValueMap;

    construct_voronoi_diagram(vd, vdc_param, voronoi_edge_to_delaunay_facet_map, grid, vertexValueMap, bbox);
    if (vdc_param.test_vor) {
        // If test_vor is true means in testing mode for voronoi diagram construction, no need for further move
        return EXIT_SUCCESS;
    }

    if (indicator)
    {
        std::cout << "Constructing Iso Surface..." << std::endl;
    }
    construct_iso_surface(vd, vdc_param, iso_surface, grid, data_grid, activeCubeCenters, voronoi_edge_to_delaunay_facet_map, vertexValueMap, bbox, point_index_map);

    // Export the Voronoi diagram to a CSV file if requested.
    if (vdc_param.out_csv)
    {
        std::cout << "Export Voronoi Diagram" << std::endl;
        export_voronoi_to_csv(vd, vdc_param.out_csv_name);
    }

    // Handle the output mesh generation and return the appropriate status.
    bool retFlag;
    int retVal = handle_output_mesh(retFlag, vd, vdc_param, iso_surface, point_index_map);
    if (retFlag)
        return retVal;

    std::cout << "Finished" << std::endl;

    return EXIT_SUCCESS;
}