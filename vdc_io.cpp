#include "vdc_io.h"

//! Writes a single-isovalue isosurface mesh in OFF format.
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

//! Writes a multi-isovalue isosurface mesh in OFF format.
void writeOFFMulti(const std::string &filename, const VoronoiDiagram &voronoiDiagram, const std::vector<std::tuple<int, int, int>> &isoTriangles, IsoSurface &iso_surface)
{
    std::ofstream out(filename);
    if (!out)
    {
        std::cerr << "Cannot open file for writing: " << filename << std::endl;
        return;
    }

    // Header
    out << "OFF\n";
    out << iso_surface.isosurfaceVertices.size() << " " << iso_surface.isosurfaceTrianglesMulti.size() << " 0\n";

    // Write vertex coordinates
    for (const auto &v : iso_surface.isosurfaceVertices)
    {
        out << v.x() << " " << v.y() << " " << v.z() << "\n";
    }

    // Write face indices
    for (const auto &triangle : iso_surface.isosurfaceTrianglesMulti)
    {
        int idx1 = std::get<0>(triangle);
        int idx2 = std::get<1>(triangle);
        int idx3 = std::get<2>(triangle);
        out << "3 " << idx1 << " " << idx2 << " " << idx3 << "\n";
    }

    out.close();
}

//! Writes a single-isovalue isosurface mesh in PLY format.
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

//! Writes a multi-isovalue isosurface mesh in PLY format.
void writePLYMulti(const std::string &filename, const VoronoiDiagram &voronoiDiagram, const std::vector<std::tuple<int, int, int>> &isoTriangles, IsoSurface &iso_surface)
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
    out << "element vertex " << iso_surface.isosurfaceVertices.size() << "\n";
    out << "property float x\n";
    out << "property float y\n";
    out << "property float z\n";
    out << "element face " << iso_surface.isosurfaceTrianglesMulti.size() << "\n";
    out << "property list uchar int vertex_indices\n";
    out << "end_header\n";

    // Write vertex coordinates
    for (const auto &vertex : iso_surface.isosurfaceVertices)
    {
        out << vertex.x() << " " << vertex.y() << " " << vertex.z() << "\n";
    }

    // Write face indices
    for (const auto &triangle : iso_surface.isosurfaceTrianglesMulti)
    {
        int idx1 = std::get<0>(triangle);
        int idx2 = std::get<1>(triangle);
        int idx3 = std::get<2>(triangle);
        out << "3 " << idx1 << " " << idx2 << " " << idx3 << "\n";
    }

    out.close();
}


//! Exports Voronoi diagram data to a CSV file for visualization and debugging
void export_voronoi_to_csv(const VoronoiDiagram &voronoiDiagram, const std::string &filename)
{
    std::ofstream file(filename);

    if (!file)
    {
        std::cerr << "Cannot open file for writing: " << filename << std::endl;
        return;
    }

    // Export vertices
    file << "vertices\n";
    for (const auto &vVertex : voronoiDiagram.voronoiVertices)
    {
        Point vertex = vVertex.vertex;
        file << vertex.x() << "," << vertex.y() << "," << vertex.z() << "\n";
    }

    // Export edges
    file << "edges\n";
    for (const auto &edge : voronoiDiagram.voronoiEdges)
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

