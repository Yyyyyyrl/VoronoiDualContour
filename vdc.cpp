#include "vdc.h"

int main(int argc, char *argv[])
{
    VoronoiDiagram vd;
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
    // Construct Voronoi Diagram
    if (indicator)
    {
        std::cout << "Constructing Voronoi diagram..." << std::endl;
    }
    construct_voronoi_vertices(vd);

    /*
    Iterate through each facets in the DT and then calculate Voronoi Edges
    */
    std::map<Object, std::vector<Facet>, ObjectComparator> delaunay_facet_to_voronoi_edge_map;
    construct_voronoi_edges(vd, delaunay_facet_to_voronoi_edge_map);

    /*
    Compute Scalar Values at Voronoi Vertices
    */
    if (indicator)
    {
        std::cout << "Computing scalar values at Voronoi vertices..." << std::endl;
    }
    std::vector<float> voronoi_vertex_values;

    compute_voronoi_values(vd, grid);

    /*
    Compute Isosurface Vertices
    */
    if (multi_isov)
    {
        /*
        Construct Voronoi Cells
        */

        construct_voronoi_cells(vd);

        Compute_Isosurface_Vertices_Multi(vd, isovalue);
    }
    else
    {
        Compute_Isosurface_Vertices_Single(vd, grid);
    }

    /*
    For each bipolar edge in the Voronoi diagram, add Delaunay triangle dual to bipolar edge.
    */

    std::cout << "Printing out Voronoi Diagram info: " << std::endl;
    std::ofstream log("vd_info.txt");
    log << vd;
    log.close();

    for (const auto& pair : vd.delaunayVertex_to_voronoiCell_index) {
    std::cout << "Delaunay Vertex: " << pair.first->point() << " -> Voronoi Cell Index: " << pair.second << "\n";
    }


    std::cout << "Checkpoint" << std::endl;
    if (multi_isov)
    {
        computeDualTrianglesMulti(vd, bbox, delaunay_facet_to_voronoi_edge_map, grid, isovalue);
    }
    else
    {
        dualTriangles = computeDualTriangles(vd.voronoiEdges, vertexValueMap, bbox, delaunay_facet_to_voronoi_edge_map, dt, grid);
    }

    for (const auto& triangle : isoTriangles) {
    std::cout << "IsoTriangle indices: " << std::get<0>(triangle) << ", "
              << std::get<1>(triangle) << ", " << std::get<2>(triangle) << "\n";
    }

    if (out_csv)
    {
        std::cout << "Export voronoi Diagram" << std::endl;
        export_voronoi_to_csv(vd, out_csv_name);
    }

    bool retFlag;
    int retVal = handle_output_mesh(retFlag, vd);
    if (retFlag)
        return retVal;

    std::cout << "Finished" << std::endl;

    return EXIT_SUCCESS;
}
