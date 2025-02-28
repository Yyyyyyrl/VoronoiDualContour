#include "test_vor.h"

void print_message()
{
    std::cout << "Usage: test_vor <input points file path ( in txt )> <output filename(without extension)>\n";
    std::cout << "- The input file should contains point coordinates in a txt file where each point takes one line\n";
    std::cout << "- This testing program will read the points from the txt file, build a delaunay triangulation and construct a voronoi diagram.\n";
    std::cout << "- The output will also be a txt file contains the voronoi diagram constructed";
}

bool readPointsFromFile(const std::string &filename, std::vector<Point> &points)
{
    std::ifstream file(filename);
    if (!file)
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return false;
    }

    std::string line;
    while (std::getline(file, line))
    {
        std::stringstream ss(line);
        double x, y, z;
        char comma1, comma2;

        if (!(ss >> x >> comma1 >> y >> comma2 >> z) || comma1 != ',' || comma2 != ',')
        {
            std::cerr << "Warning: Skipping invalid line: " << line << std::endl;
            continue;
        }

        points.emplace_back(x, y, z);
    }

    file.close();
    return true;
}

void write_triangulation(Delaunay dt, std::vector<Point> &points, std::string &input_filename)
{
    // Save points and edges to file
    std::ofstream file("triangulation.txt");
    if (!file)
    {
        std::cerr << "Error opening output file.\n";
        exit(EXIT_FAILURE);
    }

    // Write points
    file << "POINTS\n";
    for (const auto &p : points)
    {
        file << p.x() << " " << p.y() << " " << p.z() << "\n";
    }

    // Write edges
    file << "EDGES\n";
    for (auto e = dt.finite_edges_begin(); e != dt.finite_edges_end(); ++e)
    {
        auto v1 = e->first->vertex(e->second);
        auto v2 = e->first->vertex(e->third);
        file << v1->point().x() << " " << v1->point().y() << " " << v1->point().z() << " ";
        file << v2->point().x() << " " << v2->point().y() << " " << v2->point().z() << "\n";
    }

    file.close();
    std::cout << "Triangulation saved to triangulation.txt\n";
}

void write_voronoiDiagram(VoronoiDiagram &vd, std::string &output_filename) {
    std::ofstream file("voronoiDiagram.txt");
    if (!file) {
        std::cerr << "Error opening output file.\n";
        exit(EXIT_FAILURE);
    }

    //Write points
    file << "VoronoiDiagram:\n";

    // 1. Voronoi Vertices
    file << "\nVoronoiVertices:\n";
    for (size_t i = 0; i < vd.voronoiVertices.size(); ++i)
    {
        file << "Index " << i << ":\n";
        file << vd.voronoiVertices[i];
    }

    // 2. Voronoi Edges
    file << "\nVoronoiEdges:\n";
    for (const auto &edge : vd.voronoiEdges)
    {
        file << "  Edge: ";
        Segment3 segment;
        Line3 line;
        Ray3 ray;

        if (CGAL::assign(segment, edge))
        {
            file << "Segment(" << segment.source() << " - " << segment.target() << ")\n";
        }
        else if (CGAL::assign(line, edge))
        {
            file << "Line(" << line.point(0) << " - " << line.point(1) << ")\n";
        }
        else if (CGAL::assign(ray, edge))
        {
            file << "Ray(" << ray.source() << ", direction: " << ray.direction() << ")\n";
        }
        else
        {
            file << "Unknown edge type.\n";
        }
    }

    // 5. Voronoi Cells
    file << "\nVoronoiCells:\n";
    for (const auto &cell : vd.voronoiCells)
    {
        file << "\n" << cell;
    }

    // 4. Voronoi Facets
    file << "\nVoronoiFacets:\n";
    for (size_t i = 0; i < vd.voronoiFacets.size(); ++i)
    {
        file << "Index " << i << ":\n";
        file << vd.voronoiFacets[i];
    }

    file.close();
    std::cout << "voronoi diagram saved to voronoiDiagram.txt\n";
}

int main(int argc, char *argv[])
{
    Delaunay dt;
    VoronoiDiagram vd;
    std::map<Point, int> pointindexmap;

    if (argc < 2)
    {
        print_message();
        exit(EXIT_FAILURE);
    }

    std::string input_filename = argv[1];
    std::string output_filename = argv[2];

    std::vector<Point> input_points;

    if (readPointsFromFile(input_filename, input_points))
    {
        std::cout << "Successfully read " << input_points.size() << " points." << std::endl;
    }
    else
    {
        std::cerr << "Failed to read points from file." << std::endl;
    }

    dt.insert(input_points.begin(), input_points.end());

    write_triangulation(dt, input_points, input_filename);

    int index = 0;
    for (const auto &pt : input_points) {
        pointindexmap[pt] = index;
        index++;
    }

    std::map<Object, std::vector<Facet>, ObjectComparator> voronoi_edge_to_delaunay_facet_map;

    construct_voronoi_vertices(vd, dt);
    std::cout << "Checkpoint 1" << std::endl;
    construct_voronoi_edges(vd, voronoi_edge_to_delaunay_facet_map, dt);
    std::cout << "Checkpoint 2" << std::endl;
    //construct_voronoi_cells(vd,dt);
    construct_voronoi_cells_non_convex_hull(vd,dt);
    std::cout << "Checkpoint 3" << std::endl;

    write_voronoiDiagram(vd, output_filename);


}