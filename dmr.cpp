#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>
#include <limits>
#include <algorithm>
#include <teem/nrrd.h>
#include <CGAL/bounding_box.h>
#include <CGAL/intersections.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_with_circumcenter_3.h>
#include "dmr.h"

// Define the kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

// Define the vertex base with information
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K> Vb;
typedef CGAL::Delaunay_triangulation_cell_base_with_circumcenter_3<K> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;

// Define Delaunay Triangulation
typedef CGAL::Delaunay_triangulation_3<K, Tds> Delaunay;
typedef Delaunay::Point Point;
typedef Delaunay::Cell_handle Cell_handle;
typedef Delaunay::Vertex_handle Vertex_handle;
typedef Delaunay::Cell_circulator Cell_circulator;
typedef CGAL::Object Object;
typedef Delaunay::Edge Edge;
typedef Delaunay::Facet Facet;
typedef Delaunay::Triangle Triangle;
typedef K::Segment_3 Segment3;
typedef K::Ray_3 Ray3;
typedef K::Line_3 Line3;
typedef K::Point_3 Point3;

/*Helper Functions*/
// Boolean value that controlls whether the debug print commands are executed
bool debug = true;
// Boolena value controlls whether the progress print commands are executed
bool indicator = true;

// Simplest grid used for finding active cubes and their centers
struct Grid
{
    std::vector<float> data;
    int nx, ny, nz;
    double dx, dy, dz;   // Spacings for each dimension
};


struct BoundingBox
{
    double minX, maxX;
    double minY, maxY;
    double minZ, maxZ;

    BoundingBox(double x_min, double x_max, double y_min, double y_max, double z_min, double z_max)
        : minX(x_min), maxX(x_max), minY(y_min), maxY(y_max), minZ(z_min), maxZ(z_max) {}

};

// Function to compute bounding box from active cube centers
BoundingBox compute_bounding_box(const std::vector<Point> &points)
{
    double minX = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double minY = minX;
    double maxY = maxX;
    double minZ = minX;
    double maxZ = maxX;

    for (const auto &point : points)
    {
        minX = std::min(minX, point.x())-1;
        maxX = std::max(maxX, point.x())+1;
        minY = std::min(minY, point.y())-1;
        maxY = std::max(maxY, point.y())+1;
        minZ = std::min(minZ, point.z())-1;
        maxZ = std::max(maxZ, point.z())+1;
    }
    return {minX, maxX, minY, maxY, minZ, maxZ};
}


struct DelaunayTriangle
{
    Point vertex1, vertex2, vertex3;
    DelaunayTriangle(Point v1, Point v2, Point v3) : vertex1(v1), vertex2(v2), vertex3(v3) {}
};



template <typename T>
std::vector<float> convert_to_float_vector(T* data_ptr, size_t total_size) {
    std::vector<float> data(total_size);
    for (size_t i = 0; i < total_size; ++i) {
        data[i] = static_cast<float>(data_ptr[i]);
    }
    return data;
}

Grid load_nrrd_data(const std::string& file_path) {
    Nrrd* nrrd = nrrdNew();
    if (nrrdLoad(nrrd, file_path.c_str(), NULL)) {
        char *err = biffGetDone(NRRD);
        std::cerr << "Error reading NRRD file: " << err << std::endl;
        free(err);
        nrrdNuke(nrrd);
        exit(1);
    }

    size_t total_size = nrrdElementNumber(nrrd);

    std::vector<float> data;

    if (nrrd->type == nrrdTypeFloat) {
        float* data_ptr = static_cast<float*>(nrrd->data);
        data = std::vector<float>(data_ptr, data_ptr + total_size);
    } else if (nrrd->type == nrrdTypeUChar) {
        unsigned char* data_ptr = static_cast<unsigned char*>(nrrd->data);
        data = convert_to_float_vector(data_ptr, total_size);
    } else {
        std::cerr << "Unsupported NRRD data type." << std::endl;
        nrrdNuke(nrrd);
        exit(1);
    }

    int nx = nrrd->axis[0].size;
    int ny = nrrd->axis[1].size;
    int nz = nrrd->axis[2].size;
    double dx = nrrd->axis[0].spacing;
    double dy = nrrd->axis[1].spacing;
    double dz = nrrd->axis[2].spacing;


    nrrdNuke(nrrd); // Properly dispose of the Nrrd structure

    std::ofstream out_file("output_data.txt");
    for (const auto& value : data) {
        out_file << value << std::endl;
    }
    out_file.close();


    return {data, nx, ny, nz, dx, dy, dz};
}

bool is_cube_active(const Grid &grid, int x, int y, int z, float isovalue)
{
    // Define edges and check for bipolar characteristics
    std::vector<std::pair<int, int>> edges = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0}, // Bottom face edges
        {4, 5},
        {5, 6},
        {6, 7},
        {7, 4}, // Top face edges
        {0, 4},
        {1, 5},
        {2, 6},
        {3, 7} // Vertical edges
    };

    // Edge-to-vertex mapping for cube
    std::vector<std::tuple<int, int, int>> vertex_offsets = {
        {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, // Bottom face
        {0, 0, 1},
        {1, 0, 1},
        {1, 1, 1},
        {0, 1, 1} // Top face
    };

    if (debug)
    {
        std::cout << "Checking cube at (" << x << ", " << y << ", " << z << "):\n"; // Debugging output
    }
    for (const auto &edge : edges)
    {
        auto [v1x, v1y, v1z] = vertex_offsets[edge.first];
        auto [v2x, v2y, v2z] = vertex_offsets[edge.second];

        int idx1 = (x + v1x) + (y + v1y) * grid.nx + (z + v1z) * grid.nx * grid.ny;
        int idx2 = (x + v2x) + (y + v2y) * grid.nx + (z + v2z) * grid.nx * grid.ny;

        if (idx1 < grid.data.size() && idx2 < grid.data.size())
        { // Check to ensure indices are within bounds
            float val1 = grid.data[idx1];
            float val2 = grid.data[idx2];

            if ((val1 < isovalue && val2 > isovalue) || (val1 > isovalue && val2 < isovalue))
            {
                if (debug)
                {
                    std::cout << "Active edge detected, cube is active.\n"; // Debugging output
                }
                return true;
            }
        }
    }
    if (debug)
    {
        std::cout << "No active edges, cube is inactive.\n"; // Debugging output
    }
    return false;
}

std::vector<Point> find_active_cubes(const Grid &grid, float isovalue)
{
    std::vector<Point> centers;
    for (int i = 0; i < grid.nx - 1; ++i)
    {
        for (int j = 0; j < grid.ny - 1; ++j)
        {
            for (int k = 0; k < grid.nz - 1; ++k)
            {
                if (is_cube_active(grid, i, j, k, isovalue))
                {
                    centers.push_back(Point(i + 0.5, j + 0.5, k + 0.5));
                }
            }
        }
    }
    return centers;
}

// A regular 3D scalar grid data structure
struct ScalarGrid
{
    std::vector<std::vector<std::vector<double>>> data; // 3D vector to store scalar values
    int nx, ny, nz;                                     // Dimensions of the grid
    double dx, dy, dz;                                  // Voxel dimensions
    double min_x, min_y, min_z;                         // Minimum coordinates of the grid

    ScalarGrid(int nx, int ny, int nz, double dx, double dy, double dz, double min_x, double min_y, double min_z)
        : nx(nx), ny(ny), nz(nz), dx(dx), dy(dy), dz(dz), min_x(min_x), min_y(min_y), min_z(min_z)
    {
        // Initialize the 3D grid with default values (0.0)
        data.resize(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0.0)));
    }

    // get a value from the grid; returns 0 if out of bounds
    double get_value(int x, int y, int z) const
    {
        if (x < 0 || x >= nx || y < 0 || y >= ny || z < 0 || z >= nz)
        {
            return 0; // Return a default value or handle as appropriate
        }
        return data[x][y][z];
    }

    // set values in the grid for initializing the grid
    void set_value(int x, int y, int z, double value)
    {
        if (x >= 0 && x < nx && y >= 0 && y < ny && z >= 0 && z < nz)
        {
            data[x][y][z] = value;
        }
    }

    // If there are other sources of data for the grid this method is used to load it in initialization
    void load_from_source(const std::vector<std::vector<std::vector<double>>> &source)
    {
        for (int i = 0; i < nx && i < source.size(); ++i)
        {
            for (int j = 0; j < ny && j < source[i].size(); ++j)
            {
                for (int k = 0; k < nz && k < source[i][j].size(); ++k)
                {
                    data[i][j][k] = source[i][j][k];
                }
            }
        }
    }
};

/*
Purpose: Initializes the scalar grid
*/
void initialize_scalar_grid(ScalarGrid &grid, const Grid &nrrdGrid )
{
    // Use the dimensions from the loaded nrrdGrid
    grid.nx = nrrdGrid.nx;
    grid.ny = nrrdGrid.ny;
    grid.nz = nrrdGrid.nz;


    // Define grid dimensions
    grid.min_x = 0;
    grid.min_y = 0;
    grid.min_z = 0;

    // Define grid spacing (dx, dy, dz)
    grid.dx = nrrdGrid.dx;
    grid.dy = nrrdGrid.dy;
    grid.dz = nrrdGrid.dz;
    // Resizing and initializing the scalar grid data array
    grid.data.resize(grid.nx, std::vector<std::vector<double>>(grid.ny, std::vector<double>(grid.nz, 0.0)));

    // Iterate through each voxel in the grid to initialize values from the nrrdGrid data
    for (int i = 0; i < grid.nx; i++)
    {
        for (int j = 0; j < grid.ny; j++)
        {
            for (int k = 0; k < grid.nz; k++)
            {
                int index = i + j * grid.nx + k * grid.nx * grid.ny;
                if (index < nrrdGrid.data.size()) {
                    grid.data[i][j][k] = nrrdGrid.data[index];  // Assigning the value from nrrdGrid
                }
            }
        }
    }
}

const int EDGE_CONNECTIONS[12][2] = {
    {0, 1}, {1, 2}, {2, 3}, {3, 0}, // bottom face edges
    {4, 5},
    {5, 6},
    {6, 7},
    {7, 4}, // top face edges
    {0, 4},
    {1, 5},
    {2, 6},
    {3, 7} // vertical edges connecting top and bottom faces
};

/*
Purpose: Interpolates a point on the line connecting p1 and p2 based on scalar values and an isovalue.
Inputs:
    Point p1, p2: Points defining the line segment.
    double val1, val2: Scalar values at p1 and p2, respectively.
    double isovalue: The scalar value at which to interpolate.
Process: Computes the interpolation factor and returns the interpolated point.
Output: The interpolated point as a Point.
*/
Point interpolate(const Point &p1, const Point &p2, double val1, double val2, double isovalue)
{
    if (std::abs(val1 - val2) < 1e-6) // Avoid division by zero or near-zero differences
        return p1;
    double t = (isovalue - val1) / (val2 - val1);
    return Point(p1.x() + t * (p2.x() - p1.x()),
                 p1.y() + t * (p2.y() - p1.y()),
                 p1.z() + t * (p2.z() - p1.z()));
}

/*
Purpose: Processes a cube in a 3D scalar field to find vertices of the isosurface intersecting the cube.
Inputs:
    std::array<Point, 8> cube_corners: The corners of the cube.
    std::array<double, 8> scalar_values: Scalar values at each cube corner.
    double isovalue: The scalar value defining the isosurface.
Process: Checks each edge of the cube; if it intersects the isosurface, it computes the intersection point.
Output: A vector of points where the isosurface intersects the cube.
*/
std::vector<Point> process_active_cube(const std::array<Point, 8> &cube_corners, const std::array<double, 8> &scalar_values, double isovalue)
{
    std::vector<Point> intersection_points;

    for (int i = 0; i < 12; ++i)
    {
        int idx1 = EDGE_CONNECTIONS[i][0];
        int idx2 = EDGE_CONNECTIONS[i][1];
        double val1 = scalar_values[idx1];
        double val2 = scalar_values[idx2];

        if ((val1 - isovalue) * (val2 - isovalue) < 0)
        {
            Point intersect = interpolate(cube_corners[idx1], cube_corners[idx2], val1, val2, isovalue);
            intersection_points.push_back(intersect);
        }
    }

    return intersection_points;
}
/*
Purpose: Computes the geometric centroid of a collection of points.
Inputs:
    const std::vector<Point>& points: A vector of points for which the centroid is calculated.
Process: Summates the coordinates of all points and divides by the total number of points to find the average position.
Output: Returns a Point representing the centroid of the input points.
*/
Point compute_centroid(const std::vector<Point> &points)
{
    double sumX = 0, sumY = 0, sumZ = 0;
    for (const auto &pt : points)
    {
        sumX += pt.x();
        sumY += pt.y();
        sumZ += pt.z();
    }
    int n = points.size();
    return Point(sumX / n, sumY / n, sumZ / n);
}

// The function to calculate trilinear interpolation in step 4
/*
Purpose: Performs trilinear interpolation to estimate scalar values at arbitrary points within the grid.
Inputs:
    Point p: The point at which to interpolate.
    ScalarGrid grid: The scalar grid from which values are interpolated.
Process: Calculates the relative position of the point within its containing voxel and interpolates accordingly.
Output: The interpolated scalar value at point p.*/
double trilinear_interpolate(const Point &p, const ScalarGrid &grid)
{
    double gx = (p.x() - grid.min_x) / grid.dx;
    double gy = (p.y() - grid.min_y) / grid.dy;
    double gz = (p.z() - grid.min_z) / grid.dz;

    int x0 = (int)std::floor(gx);
    int x1 = x0 + 1;
    int y0 = (int)std::floor(gy);
    int y1 = y0 + 1;
    int z0 = (int)std::floor(gz);
    int z1 = z0 + 1;

    if (x0 < 0 || x1 >= grid.nx || y0 < 0 || y1 >= grid.ny || z0 < 0 || z1 >= grid.nz)
    {
        return 0; // Handle out of bounds access
    }

    double xd = gx - x0;
    double yd = gy - y0;
    double zd = gz - z0;

    double c000 = grid.get_value(x0, y0, z0);
    double c001 = grid.get_value(x0, y0, z1);
    double c010 = grid.get_value(x0, y1, z0);
    double c011 = grid.get_value(x0, y1, z1);
    double c100 = grid.get_value(x1, y0, z0);
    double c101 = grid.get_value(x1, y0, z1);
    double c110 = grid.get_value(x1, y1, z0);
    double c111 = grid.get_value(x1, y1, z1);

    double c00 = c000 * (1 - zd) + c001 * zd;
    double c01 = c010 * (1 - zd) + c011 * zd;
    double c10 = c100 * (1 - zd) + c101 * zd;
    double c11 = c110 * (1 - zd) + c111 * zd;

    double c0 = c00 * (1 - yd) + c01 * yd;
    double c1 = c10 * (1 - yd) + c11 * yd;

    return c0 * (1 - xd) + c1 * xd;
}

// Helper function to check if an edge is bipolar
/*
Purpose: Checks if an edge between two scalar values crosses the isovalue, indicating a change in sign.
Inputs:
    double val1, val2: Scalar values at the endpoints of an edge.
    double isovalue: Isovalue to check against.
Process: Returns true if the values on either side of the edge differ in sign relative to the isovalue.
Output: Boolean indicating whether the edge is bipolar.
*/
bool is_bipolar(double val1, double val2, double isovalue = 0)
{
    return (val1 - val2) * (val2 - isovalue) <= 0;
}

/*
Purpose: Computes the eight corners of a cube given its center and side length.
Inputs:
    const Point& center: The center of the cube.
    double side_length: The length of each side of the cube.
Process: Calculates the positions of all eight corners based on the center and half of the side length.
Output: Returns an std::array of Point objects representing the cube corners.
*/
std::array<Point, 8> get_cube_corners(const Point &center, double side_length)
{
    double half_side = side_length / 2.0;
    return {{
        Point(center.x() - half_side, center.y() - half_side, center.z() - half_side), // 0
        Point(center.x() + half_side, center.y() - half_side, center.z() - half_side), // 1
        Point(center.x() + half_side, center.y() + half_side, center.z() - half_side), // 2
        Point(center.x() - half_side, center.y() + half_side, center.z() - half_side), // 3
        Point(center.x() - half_side, center.y() - half_side, center.z() + half_side), // 4
        Point(center.x() + half_side, center.y() - half_side, center.z() + half_side), // 5
        Point(center.x() + half_side, center.y() + half_side, center.z() + half_side), // 6
        Point(center.x() - half_side, center.y() + half_side, center.z() + half_side)  // 7
    }};
}

/*
Purpose: Converts a point's coordinates to grid indices in the scalar grid.
Inputs:
    const Point& point: The point whose indices are to be determined.
    const ScalarGrid& grid: The grid within which the point is located.
Process: Computes the relative position of the point within the grid and converts these positions to integer grid indices.
Output: Returns a tuple of integers (int x, int y, int z) representing the grid indices.

*/
std::tuple<int, int, int> point_to_grid_index(const Point &point, const ScalarGrid &grid)
{
    int x = static_cast<int>((point.x() - grid.min_x) / grid.dx);
    int y = static_cast<int>((point.y() - grid.min_y) / grid.dy);
    int z = static_cast<int>((point.z() - grid.min_z) / grid.dz);
    return {x, y, z};
}

/*
Purpose: Retrieves the scalar value at a specific point within the grid by interpolating values at grid points.
Inputs:
    const Point& point: The point at which the scalar value is required.
    const ScalarGrid& grid: The grid from which the value is retrieved.
Process: First converts the point's coordinates to grid indices, then retrieves the scalar value at those indices using trilinear interpolation if necessary.
Output: The scalar value at the specified point, interpolated from grid values.
*/
double get_scalar_value_at_point(const Point &point, const ScalarGrid &grid)
{
    auto [x, y, z] = point_to_grid_index(point, grid);
    return grid.get_value(x, y, z);
}


void writeOFF(const std::string &filename, const std::vector<Point> &vertices, const std::vector<DelaunayTriangle> &triangles)
{
    std::ofstream out(filename);
    if (!out)
    {
        std::cerr << "Cannot open file for writing: " << filename << std::endl;
        return;
    }

    // Write the header
    out << "OFF\n";
    out << vertices.size() << " " << triangles.size() << " 0\n";

    // Write the vertices
    for (const auto &vertex : vertices)
    {
        out << vertex.x() << " " << vertex.y() << " " << vertex.z() << "\n";
    }

    // Write the faces
    for (const auto &triangle : triangles)
    {
        out << "3 " << triangle.vertex1 << " " << triangle.vertex2 << " " << triangle.vertex3 << "\n";
    }

    out.close();
}

void writePLY(const std::string &filename, const std::vector<Point> &vertices, const std::vector<DelaunayTriangle> &triangles)
{
    std::ofstream out(filename);
    if (!out)
    {
        std::cerr << "Cannot open file for writing: " << filename << std::endl;
        return;
    }

    // Write PLY header
    out << "ply\n";
    out << "format ascii 1.0\n";
    out << "element vertex " << vertices.size() << "\n";
    out << "property float x\n";
    out << "property float y\n";
    out << "property float z\n";
    out << "element face " << triangles.size() << "\n";
    out << "property list uchar int vertex_index\n";
    out << "end_header\n";

    // Write vertices
    for (const auto &vertex : vertices)
    {
        out << vertex.x() << " " << vertex.y() << " " << vertex.z() << "\n";
    }

    // Write faces
    for (const auto &triangle : triangles)
    {
        out << "3 " << triangle.vertex1 << " " << triangle.vertex2 << " " << triangle.vertex3 << "\n";
    }

    out.close();
}

std::string objectToString(const Object &obj)
{
    std::ostringstream stream;
    Segment3 seg;
    Ray3 ray;
    Line3 line;
    if (CGAL::assign(seg, obj))
    {
        stream << "Segment: " << seg;
    }
    else if (CGAL::assign(ray, obj))
    {
        stream << "Ray: " << ray;
    }
    else if (CGAL::assign(line, obj))
    {
        stream << "Line: " << line;
    }
    return stream.str();
}

bool isDegenerate(const Object &obj)
{
    Segment3 seg;
    if (CGAL::assign(seg, obj))
    {
        return seg.source() == seg.target();
    }
    // Rays and lines cannot be degenerate in the same sense as segments.
    return false;
}


// Custom comparator for CGAL::Object
struct ObjectComparator {
    bool operator()(const Object& obj1, const Object& obj2) const {
        return objectToString(obj1) < objectToString(obj2);
    }
};
/*

Body Main function

*/

int main(int argc, char *argv[])
{
    if (argc < 4)
    {
        std::cerr << "Usage: " << argv[0] << " <(nhdr/nrrd) raw data file path> <isovalue> <output format ( ply/off )> <output filename>" << std::endl;
        return EXIT_FAILURE;
    }

    // Load the scalar grid here; for example purposes, we set some arbitrary values
    ScalarGrid grid(100, 100, 100, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0);

    // Read data points and find centers of active cubes
    std::string file_path = argv[1];
    float isovalue = std::stof(argv[2]);
    std::string output_format = argv[3];
    std::string output_filename = argv[4];


    Grid data_grid = load_nrrd_data(file_path);
    std::vector<Point> activeCubeCenters = find_active_cubes(data_grid, isovalue);
    // Put data from the nrrd file into the grid
    initialize_scalar_grid(grid, data_grid);

    if (indicator)
    {
        std::cout << "Loaded data and Calculating Bounding box" << std::endl;
    }
    K::Iso_cuboid_3 bbox = CGAL::bounding_box(activeCubeCenters.begin(), activeCubeCenters.end());
    //TODO: Use the data grid


    if (debug)
    {
        std::cout << "Bounding box: ("
                  << bbox.min() << ") to ("
                  << bbox.max() << ")" << std::endl;
    }

    float cubeSideLength = 1.0;

    // Construct Delaunay Triangulation
    if (indicator)
    {
        std::cout << "Constructing Delaunay triangulation..." << std::endl;
    }
    Delaunay dt;
    int index = 0;
    for (const Point &p : activeCubeCenters)
    {
        Delaunay::Vertex_handle v = dt.insert(p);
        if (debug)
        {
            std::cout << "Inserted vertex " << index << " at (" << p << ")" << std::endl;
        }
        v->info() = index++;
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
        if (dt.is_infinite(cit))
            continue; // Skip infinite cells

        Point voronoi_vertex = dt.dual(cit);

        if (seen_points.insert(voronoi_vertex).second)
        { // insert returns a pair, where .second is a bool indicating success
            voronoi_vertices.push_back(voronoi_vertex);
            if (debug)
            {
                std::cout << "Adding Voronoi vertex at " << dt.dual(cit) << " derived from cell with vertices: "
                          << cit->vertex(0)->point() << ", "
                          << cit->vertex(1)->point() << ", "
                          << cit->vertex(2)->point() << ", "
                          << cit->vertex(3)->point() << std::endl;
            }
        }
        else if (debug)
        {
            std::cout << "Duplicate Voronoi vertex skipped: " << voronoi_vertex << " derived from cell with vertices: "
                      << cit->vertex(0)->point() << ", "
                      << cit->vertex(1)->point() << ", "
                      << cit->vertex(2)->point() << ", "
                      << cit->vertex(3)->point() << std::endl;
        }
    }
    /*
    Iterate through each facets in the DT and then calculate Voronoi Edges
    */
    std::vector<Object> voronoi_edges; // Store voronoi edges, which are dual of each facet
    std::set<std::string> seen_edges;  // Used for checking duplicate edges during iteration
    std::map<Object, Facet, ObjectComparator> delaunay_facet_to_voronoi_edge_map;
    for (Delaunay::Finite_facets_iterator fit = dt.finite_facets_begin(); fit != dt.finite_facets_end(); ++fit)
    {
        Facet facet = *fit;
        Segment3 seg;
        Ray3 ray;
        Line3 line;

        Delaunay::Object vEdge = dt.dual(facet);

        if (isDegenerate(vEdge))
        {
            if (debug)
            {
                std::cout << "skipping degenerate voronoi edge." << std::endl;
                continue;
            }
        }
        std::string edgeRep = objectToString(vEdge);
        if (seen_edges.find(edgeRep) == seen_edges.end())
        {
            voronoi_edges.push_back(vEdge);
            delaunay_facet_to_voronoi_edge_map[vEdge] = facet;
            seen_edges.insert(edgeRep);
            std::cout << "Added Voronoi Edge: " << edgeRep << std::endl;

            // returns a line segment or ray, and it's the voronoi edge
            // For ray, make a bounding box and take intersection to determine if it's bipolar

            if (debug)
            {
                if (CGAL::assign(ray, vEdge))
                {
                    std::cout << "Voronoi Edge is a ray: " << ray << std::endl;
                }
                else if (CGAL::assign(line, vEdge))
                {
                    std::cout << "Voronoi Edge is a line: " << line << std::endl;
                }
            }
        }
    }

    /*
    Compute Scalar Values at Voronoi Vertices
    */
    if (indicator)
    {
        std::cout << "Computing scalar values at Voronoi vertices..." << std::endl;
    }
    std::map<Point, double> vertexValueMap;
    std::vector<double> voronoi_vertex_values;

    for (const auto &vertex : voronoi_vertices)
    {
        double value = trilinear_interpolate(vertex, grid);
        voronoi_vertex_values.push_back(value);
        vertexValueMap[vertex] = value;
        if (debug)
        {
            std::cout << "Interpolated scalar value at (" << vertex << "): " << value << std::endl;
        }
    }

    /*
    For each bipolar edge in the Voronoi diagram, add Delaunay triangle dual to bipolar edge.
    */
    std::vector<DelaunayTriangle> dualTriangles; // To store the dual triangles of bipolar Voronoi edges
    std::map<size_t, size_t> edgeToTriangleRelation;

    for (const auto &edge : voronoi_edges)
    {
        Object intersectObj;
        Segment3 seg, iseg;
        Ray3 ray;
        Line3 line;
        Point3 p1, p2, ip;
        bool isFinite = false;

        if (CGAL::assign(seg, edge))
        {
            // If the edge is a segment
            p1 = seg.source();
            p2 = seg.target();

            // Check if it's bipolar
            // If the edge is a segment the two ends must be both in voronoi_vertices so their scalar values are pre-calculated
            if (is_bipolar(vertexValueMap[p1], vertexValueMap[p2], isovalue)) {
                    std::cout << "Bipolar edge found: " << p1 << " to " << p2 << std::endl;
                    //TODO: Find the Delaunay Triangle dual to the edge
                    std::vector<K::Point_3> delaunay_triangle;
                    Facet facet = delaunay_facet_to_voronoi_edge_map[edge];

                    for (int i = 0; i < 4; ++i) {
                        if (i != facet.second) {
                            delaunay_triangle.push_back(facet.first->vertex(i)->point());
                        }
                    }   

                    std::cout << "Delaunay triangle vertices for Voronoi edge (" << objectToString(edge) << "):\n";
                    for (const auto& vertex : delaunay_triangle) {
                        std::cout << vertex << "\n";
                    }
            }
        }
        else if (CGAL::assign(ray, edge))
        {
            // If the edge is a ray
            intersectObj = CGAL::intersection(bbox, ray);
            if (CGAL::assign(ip, intersectObj))
            {
                std::cout << "Intersection point: " << ip << " with ray: " << ray << std::endl;
            }
            else if (CGAL::assign(iseg, intersectObj))
            {
                std::cout << "Intersection seg: " << iseg.source() << " to " << iseg.target() << " with ray: " << ray << std::endl;

                // assign a corresponding scalar value to the intersection point and check if the segment between the source and intersection point is bi-polar
                Point intersection_point = iseg.target();

                double iPt_value = trilinear_interpolate(intersection_point, grid);

                std::cout << "Checking segment (" << iseg.source() << "[value:" << vertexValueMap[iseg.source()] << "]->" << iseg.target() << "[value:" << iPt_value << "): " << std::endl;

                // Check if it's bipolar
                if (is_bipolar(vertexValueMap[iseg.source()], iPt_value, isovalue)) {
                    std::cout << "Bipolar edge found: " << p1 << " to " << p2 << std::endl;
                    //TODO: Find the Delaunay Triangle dual to the edge
                    std::vector<K::Point_3> delaunay_triangle;
                    Facet facet = delaunay_facet_to_voronoi_edge_map[edge];

                    for (int i = 0; i < 4; ++i) {
                        if (i != facet.second) {
                            delaunay_triangle.push_back(facet.first->vertex(i)->point());
                        }
                    }   

                    std::cout << "Delaunay triangle vertices for Voronoi edge (" << objectToString(edge) << "):\n";
                    for (const auto& vertex : delaunay_triangle) {
                        std::cout << vertex << "\n";
                    }
                }
            }

        }
        else if (CGAL::assign(line, edge))
        {
            // If the edge is a line
            Ray3 ray1(line.point(), line.direction());
            Ray3 ray2(line.point(), -line.direction());

            intersectObj = CGAL::intersection(bbox, line);
            if (CGAL::assign(iseg, intersectObj))
            {
                std::cout << "Intersection seg: " << iseg << std::endl;
            }
            if (CGAL::do_intersect(bbox, ray1) && CGAL::do_intersect(bbox, ray2))
            {
                isFinite = true;
                if (debug)
                {
                    std::cout << "Interset with line found : " << line << std::endl;
                }
            }
        }
    }

    /* std::cout << "Identifying and storing dual triangles for bipolar Voronoi edges..." << std::endl;

    for (const auto& edge : voronoi_edges) {
        double value1 = vertexValueMap[edge.vertex1];
        double value2 = vertexValueMap[edge.vertex2];


        // Check if the edge is bipolar
        if ((value1 > isovalue && value2 < isovalue) || (value1 < isovalue && value2 > isovalue)) {
            // Edge is bipolar, get the dual Delaunay triangle
            Cell_handle cell = edge.cell;
            int i = edge.i, j = edge.j;

            size_t vertex1 = cell->vertex((i + 1) % 4)->info(); // Ensure these wrap around correctly
            size_t vertex2 = cell->vertex((i + 2) % 4)->info();
            size_t vertex3 = cell->vertex((i + 3) % 4)->info();

            dualTriangles.push_back(DelaunayTriangle(vertex1, vertex2, vertex3));

            // Map the edge index to the correspoinding Delaunay triangle index
            edgeToTriangleRelation[edge.index] = dualTriangles.size() - 1;
            std::cout << "Bipolar edge found, triangle stored: " << vertex1 << ", " << vertex2 << ", " << vertex3 << std::endl;
        }


    }

    */
    if (indicator)
    {
        std::cout << "Dual triangles of bipolar Voronoi edges:" << std::endl;
    }
    if (debug)
    {

        for (const auto &triangle : dualTriangles)
        {
            std::cout << "Triangle vertices: " << triangle.vertex1 << ", " << triangle.vertex2 << ", " << triangle.vertex3 << std::endl;
        }
    }

    /*
        Method of accessing the Dual triangle given a Voronoi Edge

    size_t edgeIndex =  .........  ;
    if (edgeToTriangleRelation.find(edgeIndex) != edgeToTriangleRelation.end()) {
        size_t triangleIndex = edgeToTriangleRelation[edgeIndex];
        DelaunayTriangle triangle = dualTriangles[triangleIndex];
        std::cout << "Delaunay Triangle for Voronoi Edge " << edgeIndex << ": vertices " << triangle.vertex1 << ", " << triangle.vertex2 << ", " << triangle.vertex3 << std::endl;
    } else {
        std::cout << "No corresponding Delaunay Triangle found for Voronoi Edge " << edgeIndex << std::endl;
    } */

    /*
    For each active cube find the location of iso-surface vertices
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

        for (const auto& center : activeCubeCenters) {
            std::vector<Point> intersectionPoints;
            double cubeSize = 1.0;
            std::array<double, 8> scalarValues;

            // Compute the global coordinates of the cube vertices and their scalar values
            for (int i = 0; i < 8; i++) {
                Point vertex(
                    center.x() + cubeVertices[i][0] - 0.5 * cubeSize,
                    center.y() + cubeVertices[i][1] - 0.5 * cubeSize,
                    center.z() + cubeVertices[i][2] - 0.5 * cubeSize
                );
                scalarValues[i] = get_scalar_value_at_point(vertex, grid); // Assuming such a function exists
            }

            // Check each edge for intersection with the isovalue
            for (const auto& edge : cubeEdges) {
                int idx1 = edge[0];
                int idx2 = edge[1];
                double val1 = scalarValues[idx1];
                double val2 = scalarValues[idx2];

                // Check if the edge crosses the isovalue
                if ((val1 - isovalue) * (val2 - isovalue) < 0) {
                    Point p1(
                        center.x() + cubeVertices[idx1][0] - 0.5 * cubeSize,
                        center.y() + cubeVertices[idx1][1] - 0.5 * cubeSize,
                        center.z() + cubeVertices[idx1][2] - 0.5 * cubeSize
                    );
                    Point p2(
                        center.x() + cubeVertices[idx2][0] - 0.5 * cubeSize,
                        center.y() + cubeVertices[idx2][1] - 0.5 * cubeSize,
                        center.z() + cubeVertices[idx2][2] - 0.5 * cubeSize
                    );
                    Point intersect = interpolate(p1, p2, val1, val2, isovalue);
                    intersectionPoints.push_back(intersect);
                    std::cout << "Intersection at: (" << intersect << ")" << std::endl;
                }
            }

            // Compute the centroid of the intersection points
            if (!intersectionPoints.empty()) {
                Point centroid = compute_centroid(intersectionPoints);
                isosurfaceVertices.push_back(centroid);
                std::cout << "Iso surface Vertex at : (" << centroid << ")" << std::endl;
            }
        }
    
    // Use locations of isosurface vertices as vertices of Delaunay triangles and write the output mesh
    if (output_format == "off")
    {
        writeOFF(output_filename, isosurfaceVertices, dualTriangles);
    }
    else if (output_format == "ply")
    {
        writePLY(output_filename, isosurfaceVertices, dualTriangles);
    }
    else
    {
        std::cerr << "Unsupported output format: " << output_format << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}