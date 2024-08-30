#include "io.h"



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


void export_voronoi_to_csv(const std::vector<Point> &voronoi_vertices, const std::vector<Object> &voronoi_edges, const std::string &filename)
{
    std::ofstream file(filename);

    // Export vertices
    file << "vertices\n";
    for (const auto &vertex : voronoi_vertices)
    {
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