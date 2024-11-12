#include "vdc.h"

int main(int argc, char *argv[])
{
    // TODO: void parse_arg()
    //  Read data points and find centers of active cubes
    parse_arguments(argc, argv);

    data_grid = load_nrrd_data(file_path);

    if (supersample)
    {
        data_grid = supersample_grid(data_grid, supersample_r);
        if (debug)
        {
            data_grid.print_grid();
        }
    }
    
    ScalarGrid grid(data_grid.nx, data_grid.ny, data_grid.nz, data_grid.dx, data_grid.dy, data_grid.dz, 0.0, 0.0, 0.0);
    // Put data from the nrrd file into the grid
    initialize_scalar_grid(grid, data_grid);
    
    
    std::vector<Cube> activeCubes = find_active_cubes(data_grid, isovalue);
    if (sep_isov)
    {
        activeCubes = separate_active_cubes_greedy(activeCubes, data_grid.nx, data_grid.ny, data_grid.nz);
    }
    activeCubeCenters = get_cube_centers(activeCubes);
    std::vector<Point> gridPoints = load_grid_points(data_grid);

    // cropAndWriteToCSV(activeCubeCenters, 48, 26, 27, 52, 30 ,30 , "../temps/fuel_cropN.csv", true);

    if (indicator)
    {
        std::cout << "Loaded data and Calculating Bounding box" << std::endl;
    }

    K::Iso_cuboid_3 bbox = CGAL::bounding_box(gridPoints.begin(), gridPoints.end());

    if (debug)
    {
        std::cout << "Bounding box: ("
                  << bbox.min() << ") to ("
                  << bbox.max() << ")" << std::endl;
    }

    float cubeSideLength = data_grid.dx;

    // Construct Delaunay Triangulation
    if (indicator)
    {
        std::cout << "Constructing Delaunay triangulation..." << std::endl;
    }

    construct_delaunay_triangulation();

    /*
    Construct Voronoi Diagram and getting the vertices, edges and cells correspondingly
    */
    std::vector<Point> voronoi_vertices;
    std::set<Point> seen_points;

    // Construct Voronoi Diagram
    if (indicator)
    {
        std::cout << "Constructing Voronoi diagram..." << std::endl;
    }
    construct_voronoi_vertices(seen_points, voronoi_vertices);

    /*
    Iterate through each facets in the DT and then calculate Voronoi Edges
    */
    std::vector<Object> voronoi_edges; // Store voronoi edges, which are dual of each facet
    std::set<std::string> seen_edges;  // Used for checking duplicate edges during iteration
    std::map<Object, std::vector<Facet>, ObjectComparator> delaunay_facet_to_voronoi_edge_map;
    construct_voronoi_edges(delaunay_facet_to_voronoi_edge_map, seen_edges, voronoi_edges);

    /*
    Compute Scalar Values at Voronoi Vertices
    */
    if (indicator)
    {
        std::cout << "Computing scalar values at Voronoi vertices..." << std::endl;
    }
    std::vector<float> voronoi_vertex_values;

    compute_voronoi_values(voronoi_vertices, grid, voronoi_vertex_values);

    /*
    Compute Isosurface Vertices
    */
    if (multi_isov)
    {
        /*
        Construct Voronoi Cells
        */

        std::vector<VoronoiCell> voronoi_cells;

        construct_voronoi_cells(voronoi_cells);

        Compute_Isosurface_Vertices_Multi(voronoi_cells);
    }
    else
    {
        Compute_Isosurface_Vertices_Single(grid);
    }

    /*
    For each bipolar edge in the Voronoi diagram, add Delaunay triangle dual to bipolar edge.
    */

    std::cout << "Checkpoint" << std::endl;

    if (multi_isov)
    {
        computeDualTrianglesMulti(voronoi_edges, bbox, delaunay_facet_to_voronoi_edge_map, grid);
    }
    else
    {
        dualTriangles = computeDualTriangles(voronoi_edges, vertexValueMap, bbox, delaunay_facet_to_voronoi_edge_map, dt, grid);
    }

    if (out_csv)
    {
        std::cout << "Export voronoi Diagram" << std::endl;
        export_voronoi_to_csv(voronoi_vertices, bipolar_voronoi_edges, out_csv_name);
    }

    bool retFlag;
    int retVal = handle_output_mesh(retFlag);
    if (retFlag)
        return retVal;

    std::cout << "Finished" << std::endl;

    return EXIT_SUCCESS;
}
