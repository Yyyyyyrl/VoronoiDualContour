#include "dmr_io.h"


void writeOFFSingle(const std::string &filename, const std::vector<Point> &vertices, const std::vector<DelaunayTriangle> &triangles, std::map<Point, int> &pointIndexMap)
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

void writeOFFMulti(const std::string &filename, const std::vector<Point> &isosurfaceVertices, const std::vector<IsoTriangle> isoTriangles){
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Cannot open file for writing: " << filename << std::endl;
        return;
    }

    // Header
    out << "OFF\n";
    out << isosurfaceVertices.size() << " " << isoTriangles.size() << " 0\n";

    // Write vertex coordinates
    for (const auto& vertex : isosurfaceVertices) {
        out << vertex << "\n";
    }

    // Write face indices
    for (const auto& triangle : isoTriangles) {
        out << "3 " << triangle.vertex_indices[0] << " " << triangle.vertex_indices[1] << " " << triangle.vertex_indices[2] << "\n";
    }

    out.close();
}

void writePLYSingle(const std::string &filename, const std::vector<Point> &vertices, const std::vector<DelaunayTriangle> &triangles, std::map<Point, int> &pointIndexMap)
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

void writePLYMulti(const std::string &filename, const std::vector<Point> &isosurfaceVertices, const std::vector<IsoTriangle> isoTriangles) {
    std::ofstream out(filename);
    if (!out)
    {
        std::cerr << "Cannot open file for writing: " << filename << std::endl;
        return;
    }

    // Write PLY header
    out << "ply\n";
    out << "format ascii 1.0\n";
    out << "element vertex " << isosurfaceVertices.size() << "\n";
    out << "property float x\n";
    out << "property float y\n";
    out << "property float z\n";
    out << "element face " << isoTriangles.size() << "\n";
    out << "property list uchar int vertex_index\n";
    out << "end_header\n";

    // Write vertex coordinates
    for (const auto& vertex : isosurfaceVertices) {
        out << vertex << "\n";
    }

    // Write face indices
    for (const auto& triangle : isoTriangles) {
        out << "3 " << triangle.vertex_indices[0] << " " << triangle.vertex_indices[1] << " " << triangle.vertex_indices[2] << "\n";
    }

    out.close();
}

void export_voronoi_to_csv(const std::vector<Point> &voronoi_vertices, const std::vector<Object> &voronoi_edges, const std::string &filename)
{
    std::ofstream file(filename);

    // Export vertices
    file << "vertices\n";
    for (const auto &vertex : voronoi_vertices)
    {
        //
        file << vertex.x() << "," << vertex.y() << "," << vertex.z() << "\n";
    }

    // Export edges
    file << "edges\n";
    for (const auto &edge : voronoi_edges)
    {
        Segment3 segment;
        Line3 line;
        Ray3 ray;

        if (CGAL::assign(segment, edge))
        {
            Point p1 = segment.source();
            Point p2 = segment.target();
            file << "Segment3," << p1.x() << "," << p1.y() << "," << p1.z() << ","
                 << p2.x() << "," << p2.y() << "," << p2.z() << "\n";
        }
        else if (CGAL::assign(line, edge))
        {
            Point p1 = line.point(0);
            Point p2 = line.point(1);
            file << "Line3," << p1.x() << "," << p1.y() << "," << p1.z() << ","
                 << p2.x() << "," << p2.y() << "," << p2.z() << "\n";
        }
        else if (CGAL::assign(ray, edge))
        {
            Point p1 = ray.source();
            Vector3 direction = ray.direction().vector();
            file << "Ray3," << p1.x() << "," << p1.y() << "," << p1.z() << ","
                 << direction.x() << "," << direction.y() << "," << direction.z() << "\n";
        }
    }

    file.close();
}

// Function to crop points based on min and max coordinates and write to CSV
void cropAndWriteToCSV(const std::vector<Point> &points, float minX, float minY, float minZ,
                       float maxX, float maxY, float maxZ, const std::string &filename, bool save_cropped)
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Failed to open the file!" << std::endl;
        return;
    }

    std::vector<Point> temp;

    // Write CSV header
    file << "x,y,z\n";

    // Iterate through the list and filter points within the specified range
    for (const auto &point : points)
    {
        float x = point.x();
        float y = point.y();
        float z = point.z();

        if (x >= minX && x <= maxX && y >= minY && y <= maxY && z >= minZ && z <= maxZ)
        {
            // Write the point to CSV
            file << x << "," << y << "," << z << "\n";
            temp.push_back(point);
        }
    }

    file.close();
    std::cout << "Points successfully written to " << filename << std::endl;
    if (save_cropped)
    {
        activeCubeCenters = temp;
    }
}