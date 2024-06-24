#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_with_circumcenter_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K> Vb;
typedef CGAL::Delaunay_triangulation_cell_base_with_circumcenter_3<K> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds> Delaunay;
typedef Delaunay::Point Point;

// A regular 3D scalar grid data structure
struct ScalarGrid {
    std::vector<std::vector<std::vector<double>>> data; // 3D vector to store scalar values
    int nx, ny, nz; // Dimensions of the grid
    double dx, dy, dz; // Voxel dimensions
    double min_x, min_y, min_z; // Minimum coordinates of the grid

    ScalarGrid(int nx, int ny, int nz, double dx, double dy, double dz, double min_x, double min_y, double min_z)
        : nx(nx), ny(ny), nz(nz), dx(dx), dy(dy), dz(dz), min_x(min_x), min_y(min_y), min_z(min_z) {
        // Initialize the 3D grid with default values (0.0)
        data.resize(nx, std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0.0)));
    }

    // get a value from the grid; returns 0 if out of bounds
    double get_value(int x, int y, int z) const {
        if (x < 0 || x >= nx || y < 0 || y >= ny || z < 0 || z >= nz) {
            return 0; // Return a default value or handle as appropriate
        }
        return data[x][y][z];
    }

    // set values in the grid for initializing the grid
    void set_value(int x, int y, int z, double value) {
        if (x >= 0 && x < nx && y >= 0 && y < ny && z >= 0 && z < nz) {
            data[x][y][z] = value;
        }
    }

    // If there are other sources of data for the grid this method is used to load it in initialization
    void load_from_source(const std::vector<std::vector<std::vector<double>>>& source) {
        for (int i = 0; i < nx && i < source.size(); ++i) {
            for (int j = 0; j < ny && j < source[i].size(); ++j) {
                for (int k = 0; k < nz && k < source[i][j].size(); ++k) {
                    data[i][j][k] = source[i][j][k];
                }
            }
        }
    }
};

/*
Purpose: Initializes the scalar grid with a simple increasing gradient.
Inputs:
    ScalarGrid& grid: A reference to the scalar grid to be initialized.
Process: Iterates over each element of the grid and sets its value to the sum of its indices.
Output: None (the grid is modified in place).
*/
void initialize_scalar_grid(ScalarGrid& grid) {
    for (int i = 0; i < grid.nx; i++) {
        for (int j = 0; j < grid.ny; j++) {
            for (int k = 0; k < grid.nz; k++) {
                grid.data[i][j][k] = i + j + k;
            }
        }
    }
}

const int EDGE_CONNECTIONS[12][2] = {
    {0, 1}, {1, 2}, {2, 3}, {3, 0}, // bottom face edges
    {4, 5}, {5, 6}, {6, 7}, {7, 4}, // top face edges
    {0, 4}, {1, 5}, {2, 6}, {3, 7}  // vertical edges connecting top and bottom faces
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
Point interpolate(const Point& p1, const Point& p2, double val1, double val2, double isovalue) {
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
std::vector<Point> process_active_cube(const std::array<Point, 8>& cube_corners, const std::array<double, 8>& scalar_values, double isovalue) {
    std::vector<Point> intersection_points;

    for (int i = 0; i < 12; ++i) {
        int idx1 = EDGE_CONNECTIONS[i][0];
        int idx2 = EDGE_CONNECTIONS[i][1];
        double val1 = scalar_values[idx1];
        double val2 = scalar_values[idx2];
        
        if ((val1 - isovalue) * (val2 - isovalue) < 0) {
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
Point compute_centroid(const std::vector<Point>& points) {
    double sumX = 0, sumY = 0, sumZ = 0;
    for (const auto& pt : points) {
        sumX += pt.x();
        sumY += pt.y();
        sumZ += pt.z();
    }
    int n = points.size();
    return Point(sumX / n, sumY / n, sumZ / n);
}


/*
Purpose: Reads points from a file, where each point is represented by its x, y, and z coordinates on separate lines.
Inputs:
    const std::string& filename: The path to the file containing the points.
Process: Opens the file, reads values line by line, constructs points, and stores them in a vector.
Output: Returns a vector of Point objects loaded from the file.
*/
std::vector<Point> read_points(const std::string& filename) {
    std::vector<Point> points;
    std::ifstream file(filename);
    double x, y, z;
    while (file >> x >> y >> z) {
        points.push_back(Point(x, y, z));
    }
    return points;
}

// The function to calculate trilinear interpolation in step 4
/*
Purpose: Performs trilinear interpolation to estimate scalar values at arbitrary points within the grid.
Inputs:
    Point p: The point at which to interpolate.
    ScalarGrid grid: The scalar grid from which values are interpolated.
Process: Calculates the relative position of the point within its containing voxel and interpolates accordingly.
Output: The interpolated scalar value at point p.*/
double trilinear_interpolate(const Point& p, const ScalarGrid& grid) {
    double gx = (p.x() - grid.min_x) / grid.dx;
    double gy = (p.y() - grid.min_y) / grid.dy;
    double gz = (p.z() - grid.min_z) / grid.dz;

    int x0 = (int)std::floor(gx);
    int x1 = x0 + 1;
    int y0 = (int)std::floor(gy);
    int y1 = y0 + 1;
    int z0 = (int)std::floor(gz);
    int z1 = z0 + 1;

    if (x0 < 0 || x1 >= grid.nx || y0 < 0 || y1 >= grid.ny || z0 < 0 || z1 >= grid.nz) {
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
bool is_bipolar(double val1, double val2, double isovalue = 0) {
    return (val1 - val2) * (val2 - isovalue) < 0;
}

/*
Purpose: Computes the eight corners of a cube given its center and side length.
Inputs:
    const Point& center: The center of the cube.
    double side_length: The length of each side of the cube.
Process: Calculates the positions of all eight corners based on the center and half of the side length.
Output: Returns an std::array of Point objects representing the cube corners.
*/
std::array<Point, 8> get_cube_corners(const Point& center, double side_length) {
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
std::tuple<int, int, int> point_to_grid_index(const Point& point, const ScalarGrid& grid) {
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
double get_scalar_value_at_point(const Point& point, const ScalarGrid& grid) {
    auto [x, y, z] = point_to_grid_index(point, grid);
    return grid.get_value(x, y, z);
}


/*
Purpose: The main execution function which reads input data, initializes structures, and executes computational tasks.

Process: Loads points, sets up the scalar grid, performs Delaunay triangulation, and identifies vertices on the isosurface.
Output: Writes the resulting data to files and potentially constructs a mesh based on the isosurface.
*/
int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <file path to input centers>" << std::endl;
        return EXIT_FAILURE;
    }

    // Load the scalar grid here; for example purposes, we set some arbitrary values
    ScalarGrid grid(100, 100, 100, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
    initialize_scalar_grid(grid);

    // Read points from file
    std::vector<Point> points = read_points(argv[1]);
    std::vector<Point> isosurface_vertices;

    // Define the cube side length;
    double cube_side_length = 1.0;
    double isovalue = 70.5;

    // Step 3: Construct Voronoi Diagram
    // Step 4: Compute scalar values at Voronoi vertices using trilinear interpolation.

    // Construct Delaunay triangulation
    Delaunay delaunay;
    int index = 0;
    for (const Point& p : points) {
        Delaunay::Vertex_handle vh = delaunay.insert(p);
        vh->info() = index++;  // Make sure this index corresponds to scalar_values
    }

    // Storage for Voronoi vertices and cells
    std::vector<Point> voronoi_vertices;
    std::vector<std::vector<Point>> voronoi_cells;
    std::vector<double> scalar_values;
    std::vector<Delaunay::Cell_handle> cells;
    
    std::ofstream vertices_file("../temps/voronoi_vertices.txt");
    std::ofstream cells_file("../temps/voronoi_cells.txt");
    std::ofstream scalars_file("../temps/scalar_values.txt");

    // Compute Voronoi vertices, cells and scalar values
    for (auto cit = delaunay.finite_cells_begin(); cit != delaunay.finite_cells_end(); ++cit) {
        Point circumcenter = cit->circumcenter();
        voronoi_vertices.push_back(circumcenter);


        vertices_file << circumcenter.x() << " " << circumcenter.y() << " " << circumcenter.z() << "\n";



        std::vector<Point> cell_vertices;
        for (int i = 0; i < 4; ++i) {  // Each Delaunay cell has four vertices
            cell_vertices.push_back(cit->vertex(i)->point());
            cells_file << "(" << cit->vertex(i)->point().x() << ", " 
                             << cit->vertex(i)->point().y() << ", " 
                             << cit->vertex(i)->point().z() << ") ";
            
        }
        voronoi_cells.push_back(cell_vertices);
        cells_file << "\n";
        cells.push_back(cit);
    }

    scalar_values.resize(points.size(), 0.0); // Resize with a default scalar value
    for (auto vit = delaunay.finite_vertices_begin(); vit != delaunay.finite_vertices_end(); ++vit) {
        Point location = vit->point();
        scalar_values[vit->info()] = trilinear_interpolate(location, grid);
        scalars_file << scalar_values[vit->info()] << "\n";
    }



    // Step 5 : Identify Bipolar Edges
    std::vector<Delaunay::Facet> delaunay_triangle_bipolar_edge;
    for (auto it = delaunay.finite_edges_begin(); it != delaunay.finite_edges_end(); ++it) {
        auto cell = it->first;
        int i = it->second, j = it->third;

        // Get the vertices of the edge directly
        auto v1 = cell->vertex(i);
        auto v2 = cell->vertex(j);

        if (v1->info() >= scalar_values.size() || v2->info() >= scalar_values.size()) {
            std::cerr << "Info index out of range" << std::endl;
            continue;
        }

        double value1 = scalar_values[v1->info()];
        double value2 = scalar_values[v2->info()];

        if (is_bipolar(value1, value2, isovalue)) {
            delaunay_triangle_bipolar_edge.push_back(Delaunay::Facet(cell, delaunay.mirror_index(cell, i)));
        }
    }

    std::vector<Point> centers = read_points(argv[1]);
    // Step 6: For each active cube, find location of isosurface vertex in active cube.
    for (const Point& center : centers) {
        auto corners = get_cube_corners(center, cube_side_length); // Example cube side length
        std::array<double, 8> values;

        for (int i = 0; i < 8; ++i) {
            values[i] = trilinear_interpolate(corners[i], grid);
        }

        std::vector<Point> intersection_points;
        for (int i = 0; i < 12; ++i) {
            int idx1 = EDGE_CONNECTIONS[i][0];
            int idx2 = EDGE_CONNECTIONS[i][1];
            if (is_bipolar(values[idx1], values[idx2], isovalue)) {
                intersection_points.push_back(interpolate(corners[idx1], corners[idx2], values[idx1], values[idx2], isovalue));
            }
        }

        if (!intersection_points.empty()) {
            Point centroid = compute_centroid(intersection_points);
            isosurface_vertices.push_back(centroid);
        }
    }


    // For checking: Save or process isosurface vertices as required
    std::ofstream output("../temps/isosurface_vertices.txt");
    for (const auto& vertex : isosurface_vertices) {
        output << vertex.x() << " " << vertex.y() << " " << vertex.z() << "\n";
    }
    output.close();


    // Step 7: Use locations of isosurface vertices as vertices of Delaunay triangles constructed in (4. / 5
    // TODO: Issues with mapping facets to isosurface vertices 

    // Output the resulting mesh in PLY format
    std::ofstream ply_file("../mesh/mesh.ply");
    if (!ply_file.is_open()) {
        std::cerr << "Failed to open file for writing." << std::endl;
        return EXIT_FAILURE;
    }


    std::ofstream mesh_file("../mesh/mesh.off");
    mesh_file << "OFF\n";
    mesh_file << delaunay.number_of_vertices() << " " << delaunay.number_of_finite_cells() << " 0\n";
    // Write PLY header
    ply_file << "ply\n";
    ply_file << "format ascii 1.0\n";
    ply_file << "element vertex " << delaunay.number_of_vertices() << "\n";
    ply_file << "property float x\n";
    ply_file << "property float y\n";
    ply_file << "property float z\n";
    ply_file << "element face " << delaunay.number_of_finite_cells() << "\n";
    ply_file << "property list uchar int vertex_index\n";
    ply_file << "end_header\n";

    // Write vertices
    std::map<Delaunay::Vertex_handle, int> vertex_indices;
    int vertex_index = 0;
    for (auto vit = delaunay.finite_vertices_begin(); vit != delaunay.finite_vertices_end(); ++vit) {
        Point p = vit->point();
        ply_file << p.x() << " " << p.y() << " " << p.z() << "\n";
        mesh_file << p.x() << " " << p.y() << " " << p.z() << "\n";
        vertex_indices[vit] = vertex_index++;
    }

    // Write faces
    // FIXME: Fix the issue in this step
    for (auto cit = delaunay.finite_cells_begin(); cit != delaunay.finite_cells_end(); ++cit) {
        ply_file << "3"; // Assuming triangular faces
        mesh_file << "3";
        for (int i = 0; i < 4; ++i) {
            if (!delaunay.is_infinite(cit->vertex(i))) {
                mesh_file << " " << vertex_indices[cit->vertex(i)];
                ply_file << " " << vertex_indices[cit->vertex(i)];
            }
        }
        ply_file << "\n";
        mesh_file << "\n";
    }

    ply_file.close();
    mesh_file.close();
    // Wrap up
    vertices_file.close();
    cells_file.close();
    scalars_file.close();
    return EXIT_SUCCESS;
}
