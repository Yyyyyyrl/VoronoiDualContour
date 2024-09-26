#include "dmr_utilities.h"
#include "debug.h"
#include "dmr_io.h"
#include <cstdlib>

/*

Body Main function

*/
// Global variables
std::string file_path;
float isovalue;
std::string output_format;
std::string output_filename;
std::string out_csv_name;
bool out_csv = false;
bool sep_isov = false;
bool supersample = false;
int supersample_r;

std::map<Point, float> vertexValueMap;
std::vector<Point> activeCubeCenters;
std::vector<Object> bipolar_voronoi_edges;


// Function to crop points based on min and max coordinates and write to CSV
void cropAndWriteToCSV(const std::vector<Point>& points, float minX, float minY, float minZ, 
                       float maxX, float maxY, float maxZ, const std::string& filename, bool save_cropped) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open the file!" << std::endl;
        return;
    }

    std::vector<Point> temp;
    
    // Write CSV header
    file << "x,y,z\n";

    // Iterate through the list and filter points within the specified range
    for (const auto& point : points) {
        std::cout << "Point: " << point << std::endl;
        float x = point.x();
        float y = point.y();
        float z = point.z();

        if (x >= minX && x <= maxX && y >= minY && y <= maxY && z >= minZ && z <= maxZ) {
            // Write the point to CSV
            file << x << "," << y << "," << z << "\n";
            temp.push_back(point);
        }
    }

    file.close();
    std::cout << "Points successfully written to " << filename << std::endl;
    if (save_cropped) {
        activeCubeCenters = temp;
    }
}

Point adjust_outside_bound_points(const Point &p, const ScalarGrid &grid, const Point &v1, const Point &v2)
{
    // Convert the point to grid space
    float gx = p.x() / grid.dx;
    float gy = p.y() / grid.dy;
    float gz = p.z() / grid.dz;

    // Check if the point is outside the bounds
    if (gx < 0 || gx >= grid.nx || gy < 0 || gy >= grid.ny || gz < 0 || gz >= grid.nz)
    {
        // If the point is out of bounds, find the closest in-bound point on the segment (v1, v2)
        // Convert v1 and v2 to grid space
        float v1_gx = v1.x() / grid.dx;
        float v1_gy = v1.y() / grid.dy;
        float v1_gz = v1.z() / grid.dz;

        float v2_gx = v2.x() / grid.dx;
        float v2_gy = v2.y() / grid.dy;
        float v2_gz = v2.z() / grid.dz;

        // Get the parameter t for the projection of point p onto the segment [v1, v2]
        float t = ((gx - v1_gx) * (v2_gx - v1_gx) + (gy - v1_gy) * (v2_gy - v1_gy) + (gz - v1_gz) * (v2_gz - v1_gz)) /
                  ((v2_gx - v1_gx) * (v2_gx - v1_gx) + (v2_gy - v1_gy) * (v2_gy - v1_gy) + (v2_gz - v1_gz) * (v2_gz - v1_gz));

        // Clamp t to [0, 1] to ensure the closest point lies on the segment
        t = std::max(0.0f, std::min(t, 1.0f));

        // Compute the closest point in grid space
        float px = v1_gx + t * (v2_gx - v1_gx);
        float py = v1_gy + t * (v2_gy - v1_gy);
        float pz = v1_gz + t * (v2_gz - v1_gz);

        // Convert the grid-space point back to world-space and return
        return Point(px * grid.dx, py * grid.dy, pz * grid.dz);
    }

    // If the point is within bounds, return the original point
    return p;
}

std::vector<DelaunayTriangle> computeDualTriangles(std::vector<CGAL::Object> &voronoi_edges, std::map<Point, float> &vertexValueMap, CGAL::Epick::Iso_cuboid_3 &bbox, std::map<Object, std::vector<Facet>, ObjectComparator> &delaunay_facet_to_voronoi_edge_map, Delaunay &dt, ScalarGrid &grid)
{
    std::vector<DelaunayTriangle> dualTriangles;
    for (const auto &edge : voronoi_edges)
    {
        Object intersectObj;
        Segment3 seg, iseg;
        Ray3 ray;
        Line3 line;
        Point3 v1, v2, ip;
        Vector3 vec1, vec2, norm;
        bool isFinite = false;

        if (CGAL::assign(seg, edge))
        {
            // If the edge is a segment
            v1 = seg.source();
            v2 = seg.target();

            // Check if it's bipolar
            // If the edge is a segment the two ends must be both in voronoi_vertices so their scalar values are pre-calculated
            if (is_bipolar(vertexValueMap[v1], vertexValueMap[v2], isovalue))
            {
                
                bipolar_voronoi_edges.push_back(edge);               // TODO: Find the Delaunay Triangle dual to the edge

                intersectObj = CGAL::intersection(bbox, Ray3(seg.source(), v2 - v1));
                CGAL::assign(iseg, intersectObj);
                Point intersection_point = iseg.target();
                CGAL::Orientation o;
                Point positive;

                if (vertexValueMap[v1] >= vertexValueMap[v2])
                {
                    positive = v1;
                }
                else
                {
                    positive = v2;
                }

                for (const auto &facet : delaunay_facet_to_voronoi_edge_map[edge])
                {
                    int iFacet = facet.second;
                    Cell_handle c = facet.first;
                    int d1, d2, d3;
                    d1 = (iFacet + 1) % 4;
                    d2 = (iFacet + 2) % 4;
                    d3 = (iFacet + 3) % 4;

                    Point p1 = c->vertex(d1)->point();
                    Point p2 = c->vertex(d2)->point();
                    Point p3 = c->vertex(d3)->point();

                    int iOrient = get_orientation(iFacet, v1, v2, vertexValueMap[v1], vertexValueMap[v2]);

                    if (dt.is_infinite(c))
                    {
                        if (iOrient < 0)
                        {
                            dualTriangles.push_back(DelaunayTriangle(p1, p2, p3));
                        }
                        else
                        {
                            dualTriangles.push_back(DelaunayTriangle(p1, p3, p2));
                        }
                    }
                    else
                    {
                        if (iOrient >= 0)
                        {
                            dualTriangles.push_back(DelaunayTriangle(p1, p2, p3));
                        }
                        else
                        {
                            dualTriangles.push_back(DelaunayTriangle(p1, p3, p2));
                        }
                    }
                }
            }
        }
        else if (CGAL::assign(ray, edge))
        {
            // If the edge is a ray
            intersectObj = CGAL::intersection(bbox, ray);

            if (CGAL::assign(iseg, intersectObj))
            {

                // assign a corresponding scalar value to the intersection point and check if the segment between the source and intersection point is bi-polar
                Point intersection_point = iseg.target();
                CGAL::Orientation o;
                Point positive;

                v1 = iseg.source();
                v2 = iseg.target();

                float iPt_value = trilinear_interpolate(adjust_outside_bound_points(intersection_point, grid, v1, v2), grid);

                if (vertexValueMap[iseg.source()] >= iPt_value)
                {
                    positive = v1;
                }
                else
                {
                    positive = v2;
                }

                // Check if it's bipolar
                if (is_bipolar(vertexValueMap[iseg.source()], iPt_value, isovalue))
                {
                    bipolar_voronoi_edges.push_back(edge);               

                    for (const auto &facet : delaunay_facet_to_voronoi_edge_map[edge])
                    {

                        Facet mirror_f = dt.mirror_facet(facet);
                        Object e = dt.dual(facet);

                        int iFacet = facet.second;
                        Cell_handle c = facet.first;
                        int d1, d2, d3;
                        d1 = (iFacet + 1) % 4;
                        d2 = (iFacet + 2) % 4;
                        d3 = (iFacet + 3) % 4;

                        Point p1 = c->vertex(d1)->point();
                        Point p2 = c->vertex(d2)->point();
                        Point p3 = c->vertex(d3)->point();

                        int iOrient = get_orientation(iFacet, v1, v2, vertexValueMap[v1], iPt_value);

                        if (dt.is_infinite(c))
                        {
                            if (iOrient < 0)
                            {
                                dualTriangles.push_back(DelaunayTriangle(p1, p2, p3));
                            }
                            else
                            {
                                dualTriangles.push_back(DelaunayTriangle(p1, p3, p2));
                            }
                        }
                        else
                        {
                            if (iOrient >= 0)
                            {
                                dualTriangles.push_back(DelaunayTriangle(p1, p2, p3));
                            }
                            else
                            {
                                dualTriangles.push_back(DelaunayTriangle(p1, p3, p2));
                            }
                        }
                    }
                }
            }
        }
        else if (CGAL::assign(line, edge))
        {
            //  If the edge is a line
            Ray3 ray1(line.point(), line.direction());
            Ray3 ray2(line.point(), -line.direction());

            intersectObj = CGAL::intersection(bbox, line);
            if (CGAL::assign(iseg, intersectObj))
            {

                Point intersection1 = iseg.source();
                Point intersection2 = iseg.target();

                float iPt1_val = trilinear_interpolate(adjust_outside_bound_points(intersection1, grid, intersection1, intersection2), grid);
                float iPt2_val = trilinear_interpolate(adjust_outside_bound_points(intersection2, grid, intersection1, intersection2), grid);

                CGAL::Orientation o;
                Point positive;

                if (iPt1_val >= iPt2_val)
                {
                    positive = intersection1;
                }
                else
                {
                    positive = intersection2;
                }

                if (is_bipolar(iPt1_val, iPt2_val, isovalue))
                {
                    bipolar_voronoi_edges.push_back(edge);
               
                    // TODO: Find the Delaunay Triangle dual to the edge

                    for (const auto &facet : delaunay_facet_to_voronoi_edge_map[edge])
                    {
                        int iFacet = facet.second;
                        Cell_handle c = facet.first;
                        int d1, d2, d3;
                        d1 = (iFacet + 1) % 4;
                        d2 = (iFacet + 2) % 4;
                        d3 = (iFacet + 3) % 4;

                        Point p1 = c->vertex(d1)->point();
                        Point p2 = c->vertex(d2)->point();
                        Point p3 = c->vertex(d3)->point();

                        int iOrient = get_orientation(iFacet, intersection1, intersection2, iPt1_val, iPt2_val);
                        if (dt.is_infinite(c))
                        {
                            if (iOrient < 0)
                            {
                                dualTriangles.push_back(DelaunayTriangle(p1, p2, p3));
                            }
                            else
                            {
                                dualTriangles.push_back(DelaunayTriangle(p1, p3, p2));
                            }
                        }
                        else
                        {
                            if (iOrient >= 0)
                            {
                                dualTriangles.push_back(DelaunayTriangle(p1, p2, p3));
                            }
                            else
                            {
                                dualTriangles.push_back(DelaunayTriangle(p1, p3, p2));
                            }
                        }
                    }
                }
            }
        }
    }
    return dualTriangles;
} // TODO: Clean up the code, and solve the orientation issue

void parse_arguments(int argc, char *argv[])
{
    // Initialize required arguments
    file_path = argv[2];
    isovalue = std::atof(argv[1]);
    output_format = argv[3];
    output_filename = argv[4];

    // Parse optional arguments
    for (int i = 5; i < argc; ++i)
    {
        std::string arg = argv[i];

        if (arg == "--out_csv" && i + 1 < argc)
        {
            out_csv = true;
            out_csv_name = argv[++i];
        }
        else if (arg == "--sep_isov")
        {
            sep_isov = true;
        }
        else if (arg == "--supersample")
        {
            supersample = true;
            supersample_r = std::atoi(argv[++i]);
        }
        else
        {
            std::cerr << "Unknown argument: " << arg << std::endl;
            exit(1);
        }
    }
}


int main(int argc, char *argv[])
{
    if (argc < 4)
    {
        std::cerr << "Usage: " << argv[0] << " <isovalue> <(nhdr/nrrd) raw data file path> <output format ( ply/off )> <output filepath> <options>\n --sep_isov : Pick a subset of non-adjacent active cubes to run \n --out_csv : Write the voronoi diagram to a csv file for visualization, follow by the path you want (--out_csv <path/to/the/output/file.csv>) \n --supersample : Supersample the input nrrd data by a factor before running the algorithm, follow by the factor(--supersample <int factor>)" << std::endl;
        return EXIT_FAILURE;
    }

    //TODO: void parse_arg()
    // Read data points and find centers of active cubes
    parse_arguments(argc, argv);
    
    Grid data_grid = load_nrrd_data(file_path);
    if (debug) {data_grid.print_grid();}
    if (supersample)
    {
        //std::cout << "Original: " << data_grid.nx << " " << data_grid.ny << " " << data_grid.nz << std::endl;
        data_grid = supersample_grid(data_grid, supersample_r);
        //std::cout << "After supersampling: " << data_grid.nx << " " << data_grid.ny << " " << data_grid.nz << std::endl;
        std::cout << "After supersample: " << std::endl;
        if (debug) {data_grid.print_grid();}
    }
    std::vector<Cube> activeCubes = find_active_cubes(data_grid, isovalue);
    if (sep_isov)
    {
        activeCubes = separate_active_cubes_greedy(activeCubes, data_grid.nx, data_grid.ny, data_grid.nz);
    }
    activeCubeCenters = get_cube_centers(activeCubes);
    std::vector<Point> gridPoints = load_grid_points(data_grid);

    ScalarGrid grid(data_grid.nx, data_grid.ny, data_grid.nz, data_grid.dx, data_grid.dy, data_grid.dz, 0.0, 0.0, 0.0);
    // Put data from the nrrd file into the grid
    initialize_scalar_grid(grid, data_grid);

    //cropAndWriteToCSV(activeCubeCenters, 27, 27, 28, 35, 29 ,33 , "../temps/fuel_crop.csv", true);

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
    Delaunay dt(activeCubeCenters.begin(), activeCubeCenters.end());
    int index = 0;

    std::map<Point, int> point_index_map;
    int i = 0;
    for (const auto &pt : activeCubeCenters)
    {
        point_index_map[pt] = i;
        i++;
    }
    

    // *** DEBUG ***
    // Print Delaunay tetrahedra.
/*         for (Delaunay::All_cells_iterator cell_it = dt.all_cells_begin();
             cell_it != dt.all_cells_end(); cell_it++)
        {

            print_cell(*cell_it);
        } */

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
    for (Delaunay::Finite_cells_iterator cit = dt.finite_cells_begin(); cit != dt.finite_cells_end(); ++cit)
    {

        Point voronoi_vertex = dt.dual(cit);

        if (seen_points.insert(voronoi_vertex).second)
        { // insert returns a pair, where .second is a bool indicating success
            voronoi_vertices.push_back(voronoi_vertex);
        }
    }

    //TODO: Use a Hash table to store the info of every voronoi vertex
    /*
    Hash Key: (p1, p2, p3, p4) points from the cit used to calculate the voronoi_vertex
    Info: Scalar Value and Whether positive/negative
    
    */
    /*
    Iterate through each facets in the DT and then calculate Voronoi Edges
    */
    std::vector<Object> voronoi_edges; // Store voronoi edges, which are dual of each facet
    std::set<std::string> seen_edges;  // Used for checking duplicate edges during iteration
    std::map<Object, std::vector<Facet>, ObjectComparator> delaunay_facet_to_voronoi_edge_map;
    for (Delaunay::Finite_facets_iterator fit = dt.finite_facets_begin(); fit != dt.finite_facets_end(); ++fit)
    {
        Facet facet = *fit;
        Segment3 seg;
        Ray3 ray;
        Line3 line;

        Delaunay::Object vEdge = dt.dual(facet);

        if (isDegenerate(vEdge))
        {
            continue;
        }
        std::string edgeRep = objectToString(vEdge);

        delaunay_facet_to_voronoi_edge_map[vEdge].push_back(facet);

        if (seen_edges.find(edgeRep) == seen_edges.end())
        {
            voronoi_edges.push_back(vEdge);
            seen_edges.insert(edgeRep);
        }
    }


    /*
    Compute Scalar Values at Voronoi Vertices
    */
    if (indicator)
    {
        std::cout << "Computing scalar values at Voronoi vertices..." << std::endl;
    }
    std::vector<float> voronoi_vertex_values;

    for (const auto &vertex : voronoi_vertices)
    {
        float value = trilinear_interpolate(vertex, grid);
        voronoi_vertex_values.push_back(value);
        vertexValueMap[vertex] = value;
    }
    std::cout << "Value at 23,30.5,35: " << vertexValueMap[Point(23,30.5,35)] << std::endl;
    std::cout << "Value at 23,30.5,33.5: " << vertexValueMap[Point(23,30.5,33.5)] << std::endl;
    std::cout << "Value at 23.5,30.5,35: " << vertexValueMap[Point(23.5,30.5,35)] << std::endl;
    std::cout << "Value at 23,31,35: " << vertexValueMap[Point(23,31,35)] << std::endl;
    std::cout << "Value at 23,30,35: " << vertexValueMap[Point(23,30,35)] << std::endl;
    std::cout << "Value at 23,30.5,43: " << vertexValueMap[Point(23,30.5,43)] << std::endl;
    std::cout << "Value at 23.5,30,35.5: " << vertexValueMap[Point(23.5,30,35.5)] << std::endl;
    std::cout << "Value at 22.5,30.5,35: " << vertexValueMap[Point(22.5,30.5,35)] << std::endl;
    std::cout << "Value at 22.5,31,34.5: " << vertexValueMap[Point(22.5,31,34.5)] << std::endl;

    /*
    For each bipolar edge in the Voronoi diagram, add Delaunay triangle dual to bipolar edge.
    */

    std::vector<DelaunayTriangle> dualTriangles = computeDualTriangles(voronoi_edges, vertexValueMap, bbox, delaunay_facet_to_voronoi_edge_map, dt, grid);

    if (out_csv)
    {
        std::cout << "Export voronoi Diagram" << std::endl;
/*         for (auto &edge : voronoi_edges) {
            std::cout << "Edge: " << edge << std::endl;   
        } */
        export_voronoi_to_csv(voronoi_vertices, voronoi_edges, out_csv_name);
    }

    /*
    Compute Isosurface Vertices
    */
    const int cubeVertices[8][3] = {
        {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, // bottom face
        {0, 0, 1},
        {1, 0, 1},
        {1, 1, 1},
        {0, 1, 1} // top face
    };

    const int cubeEdges[12][2] = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0}, // bottom face
        {4, 5},
        {5, 6},
        {6, 7},
        {7, 4}, // top face
        {0, 4},
        {1, 5},
        {2, 6},
        {3, 7} // vertical edges
    };

    std::vector<Point> isosurfaceVertices;
    if (indicator)
    {
        std::cout << "Processing each active cube to find isosurface vertices..." << std::endl;
    }

    for (const auto &center : activeCubeCenters)
    {
        std::vector<Point> intersectionPoints;
        float cubeSize = data_grid.dx;
        std::array<float, 8> scalarValues;

        /*         if (debug) {std::cout <<"Cube center: " << center.x() << " " << center.y() << " " << center.z() << std::endl;} */
        // Compute the global coordinates of the cube vertices and their scalar values
        for (int i = 0; i < 8; i++)
        {
            Point vertex(
                center.x() + (cubeVertices[i][0] - 0.5) * cubeSize,
                center.y() + (cubeVertices[i][1] - 0.5) * cubeSize,
                center.z() + (cubeVertices[i][2] - 0.5) * cubeSize);
            scalarValues[i] = grid.get_scalar_value_at_point(vertex);
        }

        // Check each edge for intersection with the isovalue
        for (const auto &edge : cubeEdges)
        {
            int idx1 = edge[0];
            int idx2 = edge[1];
            float val1 = scalarValues[idx1];
            float val2 = scalarValues[idx2];

            /*             if (debug) {std::cout << "edge: " << idx1 << " " << idx2 << " val1: " << val1 << " val2: " << val2 << std::endl;} */
            // Check if the edge crosses the isovalue
            if (((val1 >= isovalue) && (val2 < isovalue)) || ((val1 < isovalue) && (val2 >= isovalue))) // FIX: change to comparison instead of arithmatic // This is way too complicated.
            {
                Point p1(
                    center.x() + (cubeVertices[idx1][0] - 0.5) * cubeSize,
                    center.y() + (cubeVertices[idx1][1] - 0.5) * cubeSize,
                    center.z() + (cubeVertices[idx1][2] - 0.5) * cubeSize);

                Point p2(
                    center.x() + (cubeVertices[idx2][0] - 0.5) * cubeSize,
                    center.y() + (cubeVertices[idx2][1] - 0.5) * cubeSize,
                    center.z() + (cubeVertices[idx2][2] - 0.5) * cubeSize);

                Point intersect = interpolate(p1, p2, val1, val2, isovalue, data_grid);
                intersectionPoints.push_back(intersect);
            }
        }

        // Compute the centroid of the intersection points
        if (!intersectionPoints.empty())
        {
            Point centroid = compute_centroid(intersectionPoints, supersample, supersample_r);
            isosurfaceVertices.push_back(centroid);

        }
    }

    // Use locations of isosurface vertices as vertices of Delaunay triangles and write the output mesh
    if (output_format == "off")
    {
        writeOFF(output_filename, isosurfaceVertices, dualTriangles, point_index_map);
    }
    else if (output_format == "ply")
    {
        writePLY(output_filename, isosurfaceVertices, dualTriangles, point_index_map);
    }
    else
    {
        std::cerr << "Unsupported output format: " << output_format << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Finished" << std::endl;

    return EXIT_SUCCESS;
}
