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
std::map<Delaunay::Vertex_handle, std::vector<Point>> voronoi_cells;




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

    std::ofstream file("../temps/fuel-sep-sup2-record.txt");

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


                Point p1 = seg.source();
                Point p2 = seg.target();
                file << "Segment3," << p1.x() << "," << p1.y() << "," << p1.z() << ","
                 << p2.x() << "," << p2.y() << "," << p2.z() << "\n";

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
                    file << "DUAL: " << p1 << " " << p2 << " " << p3 << "\n";
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

                    Point p1 = ray.source();
                    Vector3 direction = ray.direction().vector();
                    file << "Ray3," << p1.x() << "," << p1.y() << "," << p1.z() << "," << direction.x() << "," << direction.y() << "," << direction.z() << "\n";
                    
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
                        file << "DUAL: " << p1 << " " << p2 << " " << p3 << "\n";
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

                    Point p1 = line.point(0);
                    Point p2 = line.point(1);
                    file << "Line3," << p1.x() << "," << p1.y() << "," << p1.z() << "," << p2.x() << "," << p2.y() << "," << p2.z() << "\n";
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
                        file << "DUAL: " << p1 << " " << p2 << " " << p3 << "\n";
                    }
                }
            }
        }
    }
    file.close();
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


//TODO:
/*

option: --sepPos or --sepNeg

For every Voronoi Cell Polytope P
	1. Map each bipolar edge e to a vertex v_e
	2. For each facet of P, with vertices that are either marked positive or negative (larger or smaller than isovalue)
		Let v_1, v_2, ..., v_k be vertices with v_1 being positive/negative
		If SepNeg/SepPos and there are some vertices v_i and v_j being positive/negative and every vertices in between are negative/positive ( that is, v_i+1, v_i+2, ..., v_j-1 are negative/positive):
			Take the center of edge (v_i, v_i+1) name c_i and (v_j-1, v_j) name c_j and connect the two centers to form (c_i, c_j)
	3. This process will generate cycles
	4. Add isovertex for each cycle ( calculate centroid of cycle )
	

*/

bool isPositive(double value) {
    return value >= isovalue;
}

// Helper function to get the midpoint or any point along the bipolar edge
Point getMidpointOfEdge(const Object& edge) {
    Segment3 segment;
    Ray3 ray;
    Line3 line;

    if (CGAL::assign(segment, edge)) {
        // If it's a segment, return the midpoint
        return CGAL::midpoint(segment.source(), segment.target());
    } else if (CGAL::assign(ray, edge)) {
        // If it's a ray, return a point along the ray (e.g., close to the source)
        return ray.source();
    } else if (CGAL::assign(line, edge)) {
        // If it's a line, return a point on the line (e.g., any point, say the origin)
        return line.point(0);  // Get a point on the line
    }
    throw std::runtime_error("Unknown Voronoi edge type");
}

// Map bipolar edges to a vertex
std::map<Object, Point> mapBipolarEdges(const VoronoiFacet& facet) {
    std::map<Object, Point> bipolar_edge_to_vertex;

    for (size_t i = 0; i < facet.vertices.size(); ++i) {
        Point v1 = facet.vertices[i];
        Point v2 = facet.vertices[(i + 1) % facet.vertices.size()];  // wrap around to form edges

        double scalar_v1 = facet.scalar_values[i];
        double scalar_v2 = facet.scalar_values[(i + 1) % facet.vertices.size()];

        if ((isPositive(scalar_v1) && !isPositive(scalar_v2)) || 
            (!isPositive(scalar_v1) && isPositive(scalar_v2))) {
            // This is a bipolar edge (one positive, one negative)
            // Create a CGAL::Object representing this edge
            Segment3 seg(v1, v2);
            Object edgeObject = CGAL::make_object(seg);

            // Find the midpoint or another appropriate point on the edge
            Point midpoint = CGAL::midpoint(v1, v2);

            // Map this CGAL::Object (Segment3) to the midpoint
            bipolar_edge_to_vertex[edgeObject] = midpoint;
        }
    }
    return bipolar_edge_to_vertex;
}
    
// Process facets of the Voronoi cell
void processVoronoiCell(const std::vector<VoronoiFacet>& facets) {
    std::vector<Object> new_edges;  // Store the new edges created in the process

    for (const VoronoiFacet& facet : facets) {
        // Map bipolar edges
        std::map<Object, Point> bipolar_edges = mapBipolarEdges(facet);

        // Check for sequences of alternating positive/negative vertices
        size_t n = facet.vertices.size();
        for (size_t i = 0; i < n; ++i) {
            if (isPositive(facet.scalar_values[i])) {
                size_t j = i + 1;
                while (j < n && !isPositive(facet.scalar_values[j])) {
                    ++j;
                }

                if (j < n && isPositive(facet.scalar_values[j])) {
                    // Found a sequence from v_i to v_j with all intermediate vertices negative
                    Object edge_i, edge_j;
                    Segment3 seg_i, seg_j;

                    // Extract the CGAL::Segment3 from CGAL::Object
                    if (CGAL::assign(seg_i, bipolar_edges[Object(Segment3(facet.vertices[i], facet.vertices[i + 1]))]) &&
                        CGAL::assign(seg_j, bipolar_edges[Object(Segment3(facet.vertices[j - 1], facet.vertices[j]))])) {
                        
                        Point ci = bipolar_edges[CGAL::make_object(seg_i)];
                        Point cj = bipolar_edges[CGAL::make_object(seg_j)];

                        // Create a new CGAL::Object edge between ci and cj
                        Object newEdge = CGAL::make_object(Segment3(ci, cj));
                        new_edges.push_back(newEdge);
                    }
                }
            }
        }
    }

    // Process the cycles
    // Group the new edges into cycles (closed loops)

    std::vector<std::vector<Object>> cycles = generateCycles(new_edges);

    // Step 4: Add an isovertex (centroid of each cycle)
    for (const auto& cycle : cycles) {
        Point centroid = calculateCentroid(cycle);
        // Add the isovertex to your data structure
        isosurfaceVertices.push_back(centroid);
    }
}

bool extractSegmentEndpoints(const Object& obj, Point& p1, Point& p2) {
    Segment3 seg;
    if (CGAL::assign(seg, obj)) {
        p1 = seg.source();
        p2 = seg.target();
        return true;
    }
    return false;
}

// Generate cycles from a list of edges (CGAL::Object)
std::vector<std::vector<Object>> generateCycles(const std::vector<Object>& edges) {
    // Map from each vertex to the edges connected to it
    std::map<Point, std::vector<Object>> vertex_to_edges;

    // Populate the vertex-to-edges map
    for (const auto& edge : edges) {
        Point p1, p2;
        if (extractSegmentEndpoints(edge, p1, p2)) {
            vertex_to_edges[p1].push_back(edge);
            vertex_to_edges[p2].push_back(edge);
        }
    }

    std::vector<std::vector<Object>> cycles;   // This will store the cycles
    std::set<Object> visited_edges;            // Track visited edges

    // Function to recursively follow edges and form a cycle
    auto followEdgeCycle = [&](const Object& start_edge) -> std::vector<Object> {
        std::vector<Object> cycle;
        Object current_edge = start_edge;
        Point current_vertex;
        Point next_vertex;

        // Get the first edge's endpoints
        if (!extractSegmentEndpoints(current_edge, current_vertex, next_vertex)) {
            return cycle;  // Skip if the object is not a segment
        }

        cycle.push_back(current_edge);
        visited_edges.insert(current_edge);

        // Follow the edges until a cycle is formed
        while (true) {
            // Find the next edge connected to 'next_vertex' that has not been visited
            bool found_next_edge = false;
            for (const auto& edge : vertex_to_edges[next_vertex]) {
                if (visited_edges.find(edge) == visited_edges.end()) {
                    // Found an unvisited edge connected to the next_vertex
                    visited_edges.insert(edge);
                    cycle.push_back(edge);

                    // Update the vertices for the next iteration
                    Point p1, p2;
                    extractSegmentEndpoints(edge, p1, p2);
                    if (p1 == next_vertex) {
                        current_vertex = p1;
                        next_vertex = p2;
                    } else {
                        current_vertex = p2;
                        next_vertex = p1;
                    }
                    found_next_edge = true;
                    break;
                }
            }

            // If no further connected edge is found, the cycle is complete
            if (!found_next_edge) {
                break;
            }
        }

        return cycle;
    };

    // Find cycles by following edges
    for (const auto& edge : edges) {
        if (visited_edges.find(edge) == visited_edges.end()) {
            // Start a new cycle with this unvisited edge
            std::vector<Object> cycle = followEdgeCycle(edge);
            if (!cycle.empty()) {
                cycles.push_back(cycle);  // Add the cycle to the list of cycles
            }
        }
    }

    return cycles;
}

// Helper function to calculate centroid of a cycle
Point calculateCentroid(const std::vector<Object>& cycle) {
    Point centroid(0, 0, 0);
    size_t count = 0;

    for (const auto& edgeObject : cycle) {
        Segment3 seg;
        if (CGAL::assign(seg, edgeObject)) {
            centroid = centroid + seg.source() + seg.target();
            count += 2;
        }
    }
    return centroid / count;  // Return the averaged point as the centroid
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

    //cropAndWriteToCSV(activeCubeCenters, 48, 26, 27, 52, 30 ,30 , "../temps/fuel_cropN.csv", true);

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

        
        // Iterate over the four vertices of the Delaunay cell to update their Voronoi cells
        for (int i = 0; i < 4; ++i) {
            Delaunay::Vertex_handle v = cit->vertex(i);
            
            // Store the Voronoi vertex as part of the Voronoi cell for this Delaunay vertex
            voronoi_cells[v].push_back(voronoi_vertex);
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

        // Assuming voronoi_cells is already filled with your Voronoi diagram data
    for (const auto& cell_entry : voronoi_cells) {
        const std::vector<VoronoiFacet>& facets = cell_entry.second;

        // Process each Voronoi cell
        processVoronoiCell(facets);
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
        export_voronoi_to_csv(voronoi_vertices, bipolar_voronoi_edges, out_csv_name);
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
