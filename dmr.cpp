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
#include <CGAL/Vector_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_with_circumcenter_3.h>

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
typedef K::Vector_3 Vector3;

/*Helper Functions*/
// Boolean value that controlls whether the debug print commands are executed
bool debug = false;
// Boolena value controlls whether the progress print commands are executed
bool indicator = true;


// Simple cube struct
struct Cube
{
    //The vertex stands for the lowest corner vertex of the cube, for example the cube defined by (0,0,0) and (1,1,1) will have repVertex (0,0,0)
    Point repVertex;
    int sidelength;
    Cube(Point v, int len) : repVertex(v), sidelength(len) {}
};

// Simplest grid used for finding active cubes and their centers
struct Grid
{
    std::vector<float> data;
    int nx, ny, nz;
    float dx, dy, dz; // Spacings for each dimension
};

struct DelaunayTriangle
{
    Point vertex1, vertex2, vertex3;
    DelaunayTriangle(Point v1, Point v2, Point v3) : vertex1(v1), vertex2(v2), vertex3(v3) {}
};

template <typename T>
std::vector<float> convert_to_float_vector(T *data_ptr, size_t total_size)
{
    std::vector<float> data(total_size);
    for (size_t i = 0; i < total_size; ++i)
    {
        data[i] = static_cast<float>(data_ptr[i]);
    }
    return data;
}

Grid load_nrrd_data(const std::string &file_path)
{
    Nrrd *nrrd = nrrdNew();
    if (nrrdLoad(nrrd, file_path.c_str(), NULL))
    {
        char *err = biffGetDone(NRRD);
        std::cerr << "Error reading NRRD file: " << err << std::endl;
        free(err);
        nrrdNuke(nrrd);
        exit(1);
    }

    size_t total_size = nrrdElementNumber(nrrd);

    std::vector<float> data;

    if (nrrd->type == nrrdTypeFloat)
    {
        float *data_ptr = static_cast<float *>(nrrd->data);
        data = std::vector<float>(data_ptr, data_ptr + total_size);
    }
    else if (nrrd->type == nrrdTypeUChar)
    {
        unsigned char *data_ptr = static_cast<unsigned char *>(nrrd->data);
        data = convert_to_float_vector(data_ptr, total_size);
    }
    else
    {
        std::cerr << "Unsupported NRRD data type." << std::endl;
        nrrdNuke(nrrd);
        exit(1);
    }

    int nx = nrrd->axis[0].size;
    int ny = nrrd->axis[1].size;
    int nz = nrrd->axis[2].size;
    float dx = nrrd->axis[0].spacing;
    float dy = nrrd->axis[1].spacing;
    float dz = nrrd->axis[2].spacing;

    if (debug)
    {
        std::cout << "Dimension:" << nx << " " << ny << " " << nz << std::endl;
        std::cout << "Spacing: " << dx << " " << dy << " " << dz << std::endl;
    }

    nrrdNuke(nrrd); // Properly dispose of the Nrrd structure

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

std::vector<Cube> find_active_cubes(const Grid &grid, float isovalue)
{
    std::vector<Cube> activeCubes;
    for (int i = 0; i < grid.nx - 1; ++i)
    {
        for (int j = 0; j < grid.ny - 1; ++j)
        {
            for (int k = 0; k < grid.nz - 1; ++k)
            {
                if (is_cube_active(grid, i, j, k, isovalue))
                {
                    activeCubes.push_back(Cube(Point(i,j,k), 1));
                }
            }
        }
    }
    return activeCubes;
}

std::vector<Point> get_cube_centers(const std::vector<Cube> &cubes) {
    std::vector<Point> centers;
    for ( auto &cube : cubes) {
        centers.push_back(cube.repVertex);
    }
    return centers;
}

bool is_adjacent(const Cube &cubeA, const Cube &cubeB) {
    return (std::abs(cubeA.repVertex.x() - cubeB.repVertex.x()) <= 1 &&
            std::abs(cubeA.repVertex.y() - cubeB.repVertex.y()) <= 1 &&
            std::abs(cubeA.repVertex.z() - cubeB.repVertex.z()) <= 1);
}


std::vector<Cube> separate_active_cubes(std::vector<Cube> &activeCubes) {
    std::vector<Cube> separatedCubes;

    std::vector<Cube> sortedCubes = activeCubes;
    std::sort(sortedCubes.begin(), sortedCubes.end(), [](const Cube& a, const Cube& b) {
        if (a.repVertex.x() != b.repVertex.x()) return a.repVertex.x() < b.repVertex.x();
        if (a.repVertex.y() != b.repVertex.y()) return a.repVertex.y() < b.repVertex.y();
        return a.repVertex.z() < b.repVertex.z();
    });


    //Greedy select non-adjacent cubes
    for (const Cube &c : sortedCubes) {
        bool isAdj = false;
        for (const Cube &s : separatedCubes) {
            if (is_adjacent(c, s)) {
                isAdj = true;
                break;
            }
        }

        if (!isAdj) {
            separatedCubes.push_back(c);
        }
    }

    return separatedCubes;
}

std::vector<Point> load_grid_points(const Grid &grid)
{
    std::vector<Point> points;
    for (int i = 0; i < grid.nx; ++i)
    {
        for (int j = 0; j < grid.ny; ++j)
        {
            for (int k = 0; k < grid.nz; ++k)
            {
                points.push_back(Point(i, j, k));
            }
        }
    }
    return points;
}

void print_cell(Delaunay::Cell c)
{
  using namespace std;

  cerr << "Cell: [";
  for (int i = 0; i < 4; i++) {
    cerr << "(" << c.vertex(i)->point() << ")";
    if (i < 3) { cerr << ","; }
  }
  cerr << "]" << endl;
}

// A regular 3D scalar grid data structure
struct ScalarGrid
{
    std::vector<std::vector<std::vector<float>>> data; // 3D vector to store scalar values
    int nx, ny, nz;                                    // Dimensions of the grid
    float dx, dy, dz;                                  // Voxel dimensions
    float min_x, min_y, min_z;                         // Minimum coordinates of the grid

    ScalarGrid(int nx, int ny, int nz, float dx, float dy, float dz, float min_x, float min_y, float min_z)
        : nx(nx), ny(ny), nz(nz), dx(dx), dy(dy), dz(dz), min_x(min_x), min_y(min_y), min_z(min_z)
    {
        // Initialize the 3D grid with default values (0.0)
        data.resize(nx, std::vector<std::vector<float>>(ny, std::vector<float>(nz, 0.0)));
    }

    // get a value from the grid; returns 0 if out of bounds
    float get_value(int x, int y, int z) const
    {
        if (x < 0 || x >= nx || y < 0 || y >= ny || z < 0 || z >= nz)
        {
            return 0; // Return a default value
        }
        return data[x][y][z];
    }

    // set values in the grid for initializing the grid
    void set_value(int x, int y, int z, float value)
    {
        if (x >= 0 && x < nx && y >= 0 && y < ny && z >= 0 && z < nz)
        {
            data[x][y][z] = value;
        }
    }

    // If there are other sources of data for the grid this method is used to load it in initialization
    void load_from_source(const std::vector<std::vector<std::vector<float>>> &source)
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
void initialize_scalar_grid(ScalarGrid &grid, const Grid &nrrdGrid)
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
    grid.data.resize(grid.nx);
    for (int i = 0; i < grid.nx; ++i)
    {
        grid.data[i].resize(grid.ny);
        for (int j = 0; j < grid.ny; ++j)
        {
            grid.data[i][j].resize(grid.nz, 0.0);
        }
    }
    // Iterate through each voxel in the grid to initialize values from the nrrdGrid data
    for (int i = 0; i < grid.nx; i++)
    {
        for (int j = 0; j < grid.ny; j++)
        {
            for (int k = 0; k < grid.nz; k++)
            {
                int index = i + j * grid.nx + k * grid.nx * grid.ny;
                if (index < nrrdGrid.data.size())
                {
                    grid.data[i][j][k] = nrrdGrid.data[index]; // Assigning the value from nrrdGrid
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
    float val1, val2: Scalar values at p1 and p2, respectively.
    float isovalue: The scalar value at which to interpolate.
Process: Computes the interpolation factor and returns the interpolated point.
Output: The interpolated point as a Point.
*/
Point interpolate(const Point &p1, const Point &p2, float val1, float val2, float isovalue)
{
    if (std::abs(val1 - val2) < 1e-6) // Avoid division by zero or near-zero differences
        return p1;
    float t = (isovalue - val1) / (val2 - val1);
    return Point(p1.x() + t * (p2.x() - p1.x()),
                 p1.y() + t * (p2.y() - p1.y()),
                 p1.z() + t * (p2.z() - p1.z()));
}

/*
Purpose: Processes a cube in a 3D scalar field to find vertices of the isosurface intersecting the cube.
Inputs:
    std::array<Point, 8> cube_corners: The corners of the cube.
    std::array<float, 8> scalar_values: Scalar values at each cube corner.
    float isovalue: The scalar value defining the isosurface.
Process: Checks each edge of the cube; if it intersects the isosurface, it computes the intersection point.
Output: A vector of points where the isosurface intersects the cube.
*/
std::vector<Point> process_active_cube(const std::array<Point, 8> &cube_corners, const std::array<float, 8> &scalar_values, float isovalue)
{
    std::vector<Point> intersection_points;

    for (int i = 0; i < 12; ++i)
    {
        int idx1 = EDGE_CONNECTIONS[i][0];
        int idx2 = EDGE_CONNECTIONS[i][1];
        float val1 = scalar_values[idx1];
        float val2 = scalar_values[idx2];

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
    float sumX = 0, sumY = 0, sumZ = 0;
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
float trilinear_interpolate(const Point &p, const ScalarGrid &grid)
{
    bool debug = false;
    float gx = (p.x() - grid.min_x) / grid.dx;
    float gy = (p.y() - grid.min_y) / grid.dy;
    float gz = (p.z() - grid.min_z) / grid.dz;

    if (debug)
    {
        std::cout << "grid dimension: " << grid.nx << " " << grid.ny << " " << grid.nz << std::endl;
        std::cout << "(gx, gy, gz): " << gx << " " << gy << " " << gz << std::endl;
    }

    int x0 = (int)std::floor(gx);
    if (x0 == grid.nx - 1)
    {
        --x0;
    }
    int x1 = x0 + 1;
    int y0 = (int)std::floor(gy);
    if (y0 == grid.ny - 1)
    {
        --y0;
    }
    int y1 = y0 + 1;
    int z0 = (int)std::floor(gz);
    if (z0 == grid.nz - 1)
    {
        --z0;
    }
    int z1 = z0 + 1;

    if (x0 < 0 || x1 >= grid.nx || y0 < 0 || y1 >= grid.ny || z0 < 0 || z1 >= grid.nz)
    {
        return 0; // Handle out of bounds access
    }

    float xd = gx - x0;
    float yd = gy - y0;
    float zd = gz - z0;

    float c000 = grid.get_value(x0, y0, z0);
    float c001 = grid.get_value(x0, y0, z1);
    float c010 = grid.get_value(x0, y1, z0);
    float c011 = grid.get_value(x0, y1, z1);
    float c100 = grid.get_value(x1, y0, z0);
    float c101 = grid.get_value(x1, y0, z1);
    float c110 = grid.get_value(x1, y1, z0);
    float c111 = grid.get_value(x1, y1, z1);

    if (debug)
    {
        std::cout << "Point is: (" << p << ")\n Eight corners of the cube: " << c000 << " " << c001 << " " << c010 << " " << c011 << " " << c100 << " " << c101 << " " << c110 << " " << c111 << std::endl;
        std::cout << "Two corners of the cube: (" << x0 << " " << y0 << " " << z0 << ") and (" << x1 << " " << y1 << " " << z1 << ")" << std::endl;
    }

    float c00 = c000 * (1 - zd) + c001 * zd;
    float c01 = c010 * (1 - zd) + c011 * zd;
    float c10 = c100 * (1 - zd) + c101 * zd;
    float c11 = c110 * (1 - zd) + c111 * zd;

    float c0 = c00 * (1 - yd) + c01 * yd;
    float c1 = c10 * (1 - yd) + c11 * yd;

    float c = c0 * (1 - xd) + c1 * xd;

    if (debug)
    {
        std::cout << "Result: scalar value at (" << p << ") is " << c << std::endl;
    }
    return c;
}

// Helper function to check if an edge is bipolar
/*
Purpose: Checks if an edge between two scalar values crosses the isovalue, indicating a change in sign.
Inputs:
    float val1, val2: Scalar values at the endpoints of an edge.
    float isovalue: Isovalue to check against.
Process: Returns true if the values on either side of the edge differ in sign relative to the isovalue.
Output: Boolean indicating whether the edge is bipolar.
*/
bool is_bipolar(float val1, float val2, float isovalue = 0)
{
    return (val1 - isovalue) * (val2 - isovalue) < 0;
}

/*
Purpose: Computes the eight corners of a cube given its center and side length.
Inputs:
    const Point& center: The center of the cube.
    float side_length: The length of each side of the cube.
Process: Calculates the positions of all eight corners based on the center and half of the side length.
Output: Returns an std::array of Point objects representing the cube corners.
*/
std::array<Point, 8> get_cube_corners(const Point &center, float side_length)
{
    float half_side = side_length / 2.0;
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
float get_scalar_value_at_point(const Point &point, const ScalarGrid &grid)
{
    const int iprecision = std::cout.precision();
    std::cout.precision(15);
    /* if(debug) {std::cout << "Point: (" << point << ")" << std::endl;} */
    std::cout.precision(iprecision);
    auto [x, y, z] = point_to_grid_index(point, grid);
    return grid.get_value(x, y, z);
}


void writeOFF(const std::string &filename, const std::vector<Point> &vertices, const std::vector<DelaunayTriangle> &triangles, std::map<Point, int> &pointIndexMap)
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
        out << "3 " << pointIndexMap[triangle.vertex1] << " " << pointIndexMap[triangle.vertex2] << " " << pointIndexMap[triangle.vertex3] << "\n";
    }

    out.close();
}

void writePLY(const std::string &filename, const std::vector<Point> &vertices, const std::vector<DelaunayTriangle> &triangles, std::map<Point, int> &pointIndexMap)
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
        out << "3 " << pointIndexMap[triangle.vertex1] << " " << pointIndexMap[triangle.vertex2] << " " << pointIndexMap[triangle.vertex3] << "\n";
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
struct ObjectComparator
{
    bool operator()(const Object &obj1, const Object &obj2) const
    {
        return objectToString(obj1) < objectToString(obj2);
    }
};

void print_facet(Facet f)
{
    int iFacet = f.second;
    int d1, d2, d3;
    d1 = (iFacet + 1) % 4;
    d2 = (iFacet + 2) % 4;
    d3 = (iFacet + 3) % 4;
    Cell_handle c = f.first;
    std::cout << "Facet: " << c->vertex(d1)->point() << ", " << c->vertex(d2)->point() << ", " << c->vertex(d3)->point() << std::endl;
    std::cout << "ifacet: " << iFacet << std::endl;
}

// 1 for positive, -1 for negative
int get_orientation(const int iFacet, const Point v1, const Point v2, const float f1, const float f2)
{
    bool flag_v1_positive;
    if (f1 >= f2)
    {
        flag_v1_positive = true;
    }
    else
    {
        flag_v1_positive = false;
    }
    if ((iFacet % 2) == 0)
    {
        if (flag_v1_positive)
        {
            std::cout << "+ Point: v1 (" << v1 << "), Result : Positive" << std::endl;
            return 1;
        }
        else
        {
            std::cout << "+ Point: v2 (" << v2 << "), Result : Negative" << std::endl;
            return -1;
        }
    }
    else
    {
        if (flag_v1_positive)
        {
            std::cout << "+ Point: v1 (" << v1 << "), Result : Negative" << std::endl;
            return -1;
        }
        else
        {
            std::cout << "+ Point: v2 (" << v2 << "), Result : Positive" << std::endl;
            return 1;
        }
    }
}

void export_voronoi_to_csv(const std::vector<Point>& voronoi_vertices, const std::vector<Object>& voronoi_edges, const std::string& filename) {
    std::ofstream file(filename);
    
    // Export vertices
    file << "vertices\n";
    for (const auto& vertex : voronoi_vertices) {
        file << vertex.x() << "," << vertex.y() << "," << vertex.z() << "\n";
    }

    // Export edges
    file << "edges\n";
    for (const auto& edge : voronoi_edges) {
        Segment3 segment;
        Line3 line;
        Ray3 ray;

        if (CGAL::assign(segment, edge)) {
            Point p1 = segment.source();
            Point p2 = segment.target();
            file << "Segment3," << p1.x() << "," << p1.y() << "," << p1.z() << "," 
                 << p2.x() << "," << p2.y() << "," << p2.z() << "\n";
        } 
        else if (CGAL::assign(line, edge)) {
            Point p1 = line.point(0);
            Point p2 = line.point(1);
            file << "Line3," << p1.x() << "," << p1.y() << "," << p1.z() << "," 
                 << p2.x() << "," << p2.y() << "," << p2.z() << "\n";
        } 
        else if (CGAL::assign(ray, edge)) {
            Point p1 = ray.source();
            Vector3 direction = ray.direction().vector();
            file << "Ray3," << p1.x() << "," << p1.y() << "," << p1.z() << "," 
                 << direction.x() << "," << direction.y() << "," << direction.z() << "\n";
        }
    }

    file.close();
}


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
    std::string out_csv_name = argv[5];

    Grid data_grid = load_nrrd_data(file_path);
    std::vector<Cube> activeCubes = find_active_cubes(data_grid, isovalue);
    std::vector<Point> activeCubeCenters = get_cube_centers(activeCubes);
    std::vector<Point> gridPoints = load_grid_points(data_grid);
    // Put data from the nrrd file into the grid
    initialize_scalar_grid(grid, data_grid);

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

    float cubeSideLength = 1.0;

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
    for (Delaunay::All_cells_iterator cell_it = dt.all_cells_begin();
         cell_it != dt.all_cells_end(); cell_it++) {

      print_cell(*cell_it);
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
            if (debug)
            {
                std::cout << "skipping degenerate voronoi edge:" << objectToString(vEdge) << std::endl;
                continue;
            }
        }
        std::string edgeRep = objectToString(vEdge);

        if (debug)
        {
            print_facet(facet);
            CGAL::assign(ray, vEdge);

            // std::cout << "Facet: " << facet.first->vertex(0)->point() << ", " <<  facet.first->vertex(1)->point() << ", " << facet.first->vertex(2)->point() << std::endl;
            std::cout << "vEdge: " << ray << std::endl;
        }

        delaunay_facet_to_voronoi_edge_map[vEdge].push_back(facet);

        if (seen_edges.find(edgeRep) == seen_edges.end())
        {
            voronoi_edges.push_back(vEdge);
            seen_edges.insert(edgeRep);
            if (debug)
            {
                std::cout << "Added Voronoi Edge: " << edgeRep << std::endl;
            }

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

    export_voronoi_to_csv(voronoi_vertices, voronoi_edges, out_csv_name);


    /*
    Compute Scalar Values at Voronoi Vertices
    */
    if (indicator)
    {
        std::cout << "Computing scalar values at Voronoi vertices..." << std::endl;
    }
    std::map<Point, float> vertexValueMap;
    std::vector<float> voronoi_vertex_values;

    for (const auto &vertex : voronoi_vertices)
    {
        float value = trilinear_interpolate(vertex, grid);
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
                if (debug)
                {
                    std::cout << "Bipolar edge found: " << iseg.source() << " to " << iseg.target() << std::endl;
                }
                // TODO: Find the Delaunay Triangle dual to the edge

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

                    std::cout << "Voronoi Segment" << std::endl;
                    std::cout << c->vertex(0)->point() << " , " << c->vertex(1)->point() << " , " << c->vertex(2)->point() << " , " << c->vertex(3)->point() << std::endl;
                    std::cout << "Value of iFacet: " << iFacet << " Point indices: p1 (" << d1 << ") p2:(" << d2 << ") p3:(" << d3 << ")" << std::endl;

                    /*                     vec1 = p2 - p1;
                                        vec2 = p3 - p1;
                                        norm = CGAL::cross_product(vec1, vec2);
                                        float fOri = norm * (positive - p1);

                                        if (fOri >= 0)
                                        { // Above plane
                                            dualTriangles.push_back(DelaunayTriangle(p1, p3, p2));
                                            std::cout << "Triangle vertices: " << p1 << ", " << p3 << ", " << p2 << std::endl;
                                        }
                                        else
                                        {
                                            dualTriangles.push_back(DelaunayTriangle(p1, p2, p3));
                                            std::cout << "Triangle vertices: " << p1 << ", " << p2 << ", " << p3 << std::endl;
                                        } */

                    int iOrient = get_orientation(iFacet, v1, v2, vertexValueMap[v1], vertexValueMap[v2]);

                    if (dt.is_infinite(c))
                    {
                        if (iOrient < 0)
                        {
                            dualTriangles.push_back(DelaunayTriangle(p1, p2, p3));
                            std::cout << "Triangle vertices: (123)" << p1 << ", " << p2 << ", " << p3 << std::endl;
                        }
                        else
                        {
                            dualTriangles.push_back(DelaunayTriangle(p1, p3, p2));
                            std::cout << "Triangle vertices: (132)" << p1 << ", " << p3 << ", " << p2 << std::endl;
                        }
                    }
                    else
                    {
                        if (iOrient >= 0)
                        {
                            dualTriangles.push_back(DelaunayTriangle(p1, p2, p3));
                            std::cout << "Triangle vertices: (123)" << p1 << ", " << p2 << ", " << p3 << std::endl;
                        }
                        else
                        {
                            dualTriangles.push_back(DelaunayTriangle(p1, p3, p2));
                            std::cout << "Triangle vertices: (132)" << p1 << ", " << p3 << ", " << p2 << std::endl;
                        }
                    }
                }
            }
        }
        else if (CGAL::assign(ray, edge))
        {
            // If the edge is a ray
            intersectObj = CGAL::intersection(bbox, ray);
            /*             if (CGAL::assign(ip, intersectObj))
                        {
                            if (debug)
                            {
                                std::cout << "Intersection point: " << ip << " with ray: " << ray << std::endl;
                            }
                        } */
            if (CGAL::assign(iseg, intersectObj))
            {
                if (debug)
                {
                    std::cout << "Intersection seg: " << iseg.source() << " to " << iseg.target() << " with ray: " << ray << std::endl;
                }

                // assign a corresponding scalar value to the intersection point and check if the segment between the source and intersection point is bi-polar
                Point intersection_point = iseg.target();
                CGAL::Orientation o;
                Point positive;

                v1 = iseg.source();
                v2 = iseg.target();

                float iPt_value = trilinear_interpolate(intersection_point, grid);

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
                    if (debug)
                    {
                        std::cout << "Bipolar edge found: " << iseg.source() << " to " << iseg.target() << std::endl;
                    }

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

                        std::cout << "Voronoi Ray" << std::endl;
                        if (dt.is_infinite(c))
                        {
                            std::cout << "Infinite cell" << std::endl;
                        }
                        std::cout << c->vertex(0)->point() << " , " << c->vertex(1)->point() << " , " << c->vertex(2)->point() << " , " << c->vertex(3)->point() << std::endl;
                        std::cout << "Value of iFacet: " << iFacet << " Point indices: p1 (" << d1 << ") p2:(" << d2 << ") p3:(" << d3 << ")" << std::endl;

                        /*                         vec1 = p2 - p1;
                                                vec2 = p3 - p1;
                                                norm = CGAL::cross_product(vec1, vec2);
                                                float fOri = norm * (positive - p1);

                                                if (fOri >= 0)
                                                { // Above plane
                                                    dualTriangles.push_back(DelaunayTriangle(p1, p3, p2));
                                                    std::cout << "Triangle vertices: " << p1 << ", " << p3 << ", " << p2 << std::endl;
                                                }
                                                else
                                                {
                                                    dualTriangles.push_back(DelaunayTriangle(p1, p2, p3));
                                                    std::cout << "Triangle vertices: " << p1 << ", " << p2 << ", " << p3 << std::endl;
                                                } */

                        int iOrient = get_orientation(iFacet, v1, v2, vertexValueMap[v1], iPt_value);

                        if (dt.is_infinite(c))
                        {
                            if (iOrient < 0)
                            {
                                dualTriangles.push_back(DelaunayTriangle(p1, p2, p3));
                                std::cout << "Triangle vertices: (123)" << p1 << ", " << p2 << ", " << p3 << std::endl;
                            }
                            else
                            {
                                dualTriangles.push_back(DelaunayTriangle(p1, p3, p2));
                                std::cout << "Triangle vertices: (132)" << p1 << ", " << p3 << ", " << p2 << std::endl;
                            }
                        }
                        else
                        {
                            if (iOrient >= 0)
                            {
                                dualTriangles.push_back(DelaunayTriangle(p1, p2, p3));
                                std::cout << "Triangle vertices: (123)" << p1 << ", " << p2 << ", " << p3 << std::endl;
                            }
                            else
                            {
                                dualTriangles.push_back(DelaunayTriangle(p1, p3, p2));
                                std::cout << "Triangle vertices: (132)" << p1 << ", " << p3 << ", " << p2 << std::endl;
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

                if (debug)
                {
                    std::cout << "Intersection seg: " << iseg << std::endl;
                }
                Point intersection1 = iseg.source();
                Point intersection2 = iseg.target();

                float iPt1_val = trilinear_interpolate(intersection1, grid);
                float iPt2_val = trilinear_interpolate(intersection2, grid);

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
                    if (debug)
                    {
                        std::cout << "Bipolar edge found: " << iseg.source() << " to " << iseg.target() << std::endl;
                    }
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

                        std::cout << "Voronoi Line" << std::endl;
                        std::cout << c->vertex(0)->point() << " , " << c->vertex(1)->point() << " , " << c->vertex(2)->point() << " , " << c->vertex(3)->point() << std::endl;
                        std::cout << "Value of iFacet: " << iFacet << " Point indices: p1 (" << d1 << ") p2:(" << d2 << ") p3:(" << d3 << ")" << std::endl;

                        /*                         vec1 = p2 - p1;
                                                vec2 = p3 - p1;
                                                norm = CGAL::cross_product(vec1, vec2);
                                                float fOri = norm * (positive - p1);

                                                if (fOri >= 0)
                                                { // Above plane
                                                    dualTriangles.push_back(DelaunayTriangle(p1, p3, p2));
                                                    std::cout << "Triangle vertices: " << p1 << ", " << p3 << ", " << p2 << std::endl;
                                                }
                                                else
                                                {
                                                    dualTriangles.push_back(DelaunayTriangle(p1, p2, p3));
                                                    std::cout << "Triangle vertices: " << p1 << ", " << p2 << ", " << p3 << std::endl;
                                                } */

                        int iOrient = get_orientation(iFacet, intersection1, intersection2, iPt1_val, iPt2_val);
                        if (dt.is_infinite(c))
                        {
                            if (iOrient < 0)
                            {
                                dualTriangles.push_back(DelaunayTriangle(p1, p2, p3));
                                std::cout << "Triangle vertices: (123)" << p1 << ", " << p2 << ", " << p3 << std::endl;
                            }
                            else
                            {
                                dualTriangles.push_back(DelaunayTriangle(p1, p3, p2));
                                std::cout << "Triangle vertices: (132)" << p1 << ", " << p3 << ", " << p2 << std::endl;
                            }
                        }
                        else
                        {
                            if (iOrient >= 0)
                            {
                                dualTriangles.push_back(DelaunayTriangle(p1, p2, p3));
                                std::cout << "Triangle vertices: (123)" << p1 << ", " << p2 << ", " << p3 << std::endl;
                            }
                            else
                            {
                                dualTriangles.push_back(DelaunayTriangle(p1, p3, p2));
                                std::cout << "Triangle vertices: (132)" << p1 << ", " << p3 << ", " << p2 << std::endl;
                            }
                        }
                    }
                }
            }
        }
    }
    /*
        if (debug)
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
        float cubeSize = 1.0;
        std::array<float, 8> scalarValues;

        /*         if (debug) {std::cout <<"Cube center: " << center.x() << " " << center.y() << " " << center.z() << std::endl;} */
        // Compute the global coordinates of the cube vertices and their scalar values
        for (int i = 0; i < 8; i++)
        {
            Point vertex(
                center.x() + cubeVertices[i][0] - 0.5 * cubeSize,
                center.y() + cubeVertices[i][1] - 0.5 * cubeSize,
                center.z() + cubeVertices[i][2] - 0.5 * cubeSize);
            scalarValues[i] = get_scalar_value_at_point(vertex, grid);
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
            if (((val1 > isovalue) && (val2 <= isovalue)) || ((val1 >= isovalue) && (val2 < isovalue)) || ((val1 < isovalue) && (val2 >= isovalue)) || ((val1 <= isovalue) && (val2 > isovalue))) // FIX: change to comparison, no arithmatic
            {
                Point p1(
                    center.x() + cubeVertices[idx1][0] - 0.5 * cubeSize,
                    center.y() + cubeVertices[idx1][1] - 0.5 * cubeSize,
                    center.z() + cubeVertices[idx1][2] - 0.5 * cubeSize);
                Point p2(
                    center.x() + cubeVertices[idx2][0] - 0.5 * cubeSize,
                    center.y() + cubeVertices[idx2][1] - 0.5 * cubeSize,
                    center.z() + cubeVertices[idx2][2] - 0.5 * cubeSize);
                Point intersect = interpolate(p1, p2, val1, val2, isovalue);
                intersectionPoints.push_back(intersect);
                /*
                                if (debug)
                                {
                                    std::cout << "p1: (" << p1 << ")  val1: " << val1 << " idx1 : " << idx1 << std::endl;
                                    std::cout << "p2: (" << p2 << ")  val2: " << val2 << " idx2 : " << idx2 << std::endl;
                                    std::cout << "Intersection at: (" << intersect << ")" << std::endl;
                                } */
            }
        }

        // Compute the centroid of the intersection points
        if (!intersectionPoints.empty())
        {
            Point centroid = compute_centroid(intersectionPoints);
            isosurfaceVertices.push_back(centroid);
            if (debug)

            {
                std::cout << "Iso surface Vertex at : (" << centroid << ")" << std::endl;
            }
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

// TODO: Clean up the code, and solve the orientation issue