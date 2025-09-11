#include "test/test_vor.h"
#include "core/vdc_utilities.h"  // For readPointsFromFile and write_voronoiDiagram

void print_message()
{
    std::cout << "Usage for TXT: test_vor <input points file path (in txt)> <output filename (without extension)>\n";
    std::cout << "Usage for NRRD/NHDR: test_vor <isovalue> <input nrrd/nhdr file path> <output filename (without extension)>\n";
    std::cout << "- For TXT: The input file should contain point coordinates where each point takes one line in format x,y,z\n";
    std::cout << "- For NRRD/NHDR: Provide isovalue first, then file path.\n";
    std::cout << "- This testing program will build a Delaunay triangulation and construct a Voronoi diagram.\n";
    std::cout << "- The output will be a txt file containing the Voronoi diagram if the validity check is passed, otherwise there will be an error.\n";
}

bool is_nrrd_file(const std::string& filename) {
    size_t dot_pos = filename.find_last_of('.');
    if (dot_pos != std::string::npos) {
        std::string ext = filename.substr(dot_pos);
        return (ext == ".nrrd" || ext == ".nhdr");
    }
    return false;
}

int main(int argc, char *argv[])
{
    if (argc < 3 || argc > 4) {
        print_message();
        exit(EXIT_FAILURE);
    }

    Delaunay dt;
    VoronoiDiagram vd;
    std::string output_filename;
    K::Iso_cuboid_3 bbox;
    VDC_PARAM vdc_param; // Hardcode for test
    vdc_param.multi_isov = true;
    vdc_param.convex_hull = true;

    std::vector<Point> input_points;
    UnifiedGrid data_grid;
    std::vector<Cube> activeCubes;
    std::vector<std::vector<GRID_FACETS>> grid_facets;
    std::vector<Point> activeCubeCenters;

    bool is_nrrd = false;
    float isovalue = 0.0f;
    std::string input_filename;

    if (argc == 3) {
        // TXT mode
        input_filename = argv[1];
        output_filename = argv[2];
        if (is_nrrd_file(input_filename)) {
            std::cerr << "For NRRD/NHDR, provide isovalue as first argument.\n";
            print_message();
            exit(EXIT_FAILURE);
        }
        if (!readPointsFromFile(input_filename, input_points)) {
            std::cerr << "Failed to read points from file." << std::endl;
            exit(EXIT_FAILURE);
        }
        std::cout << "Successfully read " << input_points.size() << " points." << std::endl;

        // Compute bbox from points with small epsilon padding
        if (input_points.empty()) {
            std::cerr << "No points to process." << std::endl;
            exit(EXIT_FAILURE);
        }
        double eps = 1e-6;
        double minx = input_points[0].x(), miny = input_points[0].y(), minz = input_points[0].z();
        double maxx = minx, maxy = miny, maxz = minz;
        for (const auto& p : input_points) {
            minx = std::min(minx, p.x());
            miny = std::min(miny, p.y());
            minz = std::min(minz, p.z());
            maxx = std::max(maxx, p.x());
            maxy = std::max(maxy, p.y());
            maxz = std::max(maxz, p.z());
        }
        Point p_min(minx - eps, miny - eps, minz - eps);
        Point p_max(maxx + eps, maxy + eps, maxz + eps);
        bbox = K::Iso_cuboid_3(p_min, p_max);

        // Insert points into DT (no dummies for TXT)
        for (size_t i = 0; i < input_points.size(); ++i) {
            Vertex_handle vh = dt.insert(input_points[i]);
            vh->info().index = i;
            vh->info().is_dummy = false;
        }
    } else { // argc == 4, NRRD mode
        isovalue = std::atof(argv[1]);
        input_filename = argv[2];
        output_filename = argv[3];
        if (!is_nrrd_file(input_filename)) {
            std::cerr << "Expected NRRD/NHDR file for 4 arguments.\n";
            print_message();
            exit(EXIT_FAILURE);
        }
        data_grid = load_nrrd_data(input_filename);
        vdc_param.isovalue = isovalue;
        find_active_cubes(data_grid, isovalue, activeCubes);
        grid_facets = create_grid_facets(activeCubes);
        activeCubeCenters = get_cube_centers(activeCubes);

        Point p_min(data_grid.min_x, data_grid.min_y, data_grid.min_z);
        Point p_max(data_grid.max_x, data_grid.max_y, data_grid.max_z);
        bbox = K::Iso_cuboid_3(p_min, p_max);

        construct_delaunay_triangulation(dt, data_grid, grid_facets, vdc_param, activeCubeCenters);
        is_nrrd = true;
    }

    // Construct Voronoi diagram (skip values as not needed for test)
    construct_voronoi_vertices(vd, dt);
    construct_voronoi_edges(vd, dt);
    if (vdc_param.multi_isov) {
        if (vdc_param.convex_hull) {
            construct_voronoi_cells_as_convex_hull(vd, dt);
        } else {
            construct_voronoi_cells_from_delaunay_triangulation(vd, dt);
        }
    }
    VoronoiDiagram vd2 = collapseSmallEdges(vd, 0.001, bbox, dt);
    vd2.check(true);
    vd = std::move(vd2);
    if (vdc_param.multi_isov) {
        construct_voronoi_cell_edges(vd, bbox, dt);
    }

    // Write output
    write_voronoiDiagram(vd, output_filename);

    return EXIT_SUCCESS;
}
