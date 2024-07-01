#include <iostream>
#include <vector>
#include <tuple>
#include <fstream>
#include <teem/nrrd.h>

struct Grid {
    std::vector<float> data;
    int nx, ny, nz;
};

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

    std::vector<float> data(total_size);

    float* data_ptr = static_cast<float*>(nrrd->data);

    std::copy(data_ptr, data_ptr + total_size, data.begin());

    int nx = nrrd->axis[0].size;
    int ny = nrrd->axis[1].size;
    int nz = nrrd->axis[2].size;

    nrrdNuke(nrrd); // Properly dispose of the Nrrd structure

    std::ofstream out_file("output_data.txt");
    for (const auto& value : data) {
        out_file << value << std::endl;
    }
    out_file.close();


    return {data, nx, ny, nz};
}

bool is_cube_active(const Grid& grid, int x, int y, int z, float isovalue) {
    // Define edges and check for bipolar characteristics
    std::vector<std::pair<int, int>> edges = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0}, // Bottom face edges
        {4, 5}, {5, 6}, {6, 7}, {7, 4}, // Top face edges
        {0, 4}, {1, 5}, {2, 6}, {3, 7}  // Vertical edges
    };

    // Edge-to-vertex mapping for cube
    std::vector<std::tuple<int, int, int>> vertex_offsets = {
        {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, // Bottom face
        {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}  // Top face
    };

    //std::cout << "Checking cube at (" << x << ", " << y << ", " << z << "):\n"; // Debugging output
    for (const auto& edge : edges) {
        auto [v1x, v1y, v1z] = vertex_offsets[edge.first];
        auto [v2x, v2y, v2z] = vertex_offsets[edge.second];

        int idx1 = (x + v1x) + (y + v1y) * grid.nx + (z + v1z) * grid.nx * grid.ny;
        int idx2 = (x + v2x) + (y + v2y) * grid.nx + (z + v2z) * grid.nx * grid.ny;

        if (idx1 < grid.data.size() && idx2 < grid.data.size()) { // Check to ensure indices are within bounds
            float val1 = grid.data[idx1];
            float val2 = grid.data[idx2];
            //std::cout << "Edge (" << edge.first << ", " << edge.second << ") values: " << val1 << ", " << val2 << "\n"; // Debugging output


            if ((val1 < isovalue && val2 >isovalue) || (val1 > isovalue && val2 < isovalue)) {
                //std::cout << "Active edge detected, cube is active.\n"; // Debugging output
                return true;
            }
        }

    }
    //std::cout << "No active edges, cube is inactive.\n"; // Debugging output
    return false;
}

std::vector<std::tuple<float, float, float>> find_active_cubes(const Grid& grid, float isovalue) {
    std::vector<std::tuple<float, float, float>> centers;
    for (int i = 0; i < grid.nx - 1; ++i) {
        for (int j = 0; j < grid.ny - 1; ++j) {
            for (int k = 0; k < grid.nz - 1; ++k) {
                if (is_cube_active(grid, i, j, k, isovalue)) {
                    centers.emplace_back(i + 0.5, j + 0.5, k + 0.5);
                    //std::cout << "Center added: (" << i + 0.5 << ", " << j + 0.5 << ", " << k + 0.5 << ")\n\n"; // Debugging output
                }
            }
        }
    }
    return centers;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <file_path> <isovalue>\n";
        return 1;
    }

    int count = 0;

    std::string file_path = argv[1];
    float isovalue = std::stof(argv[2]);
    
    Grid grid = load_nrrd_data(file_path);

    std::cout << "Data size: " << grid.nx << ", " << grid.ny << ", " << grid.nz << std::endl;

    std::cout << "Raw data loaded:\n";
    for (int k = 0; k < grid.nz; ++k) {
        for (int j = 0; j < grid.ny; ++j) {
            for (int i = 0; i < grid.nx; ++i) {
                //std::cout << grid.data[i + j * grid.nx + k * grid.nx * grid.ny] << " ";
            }
            //std::cout << std::endl;
        }
        //std::cout << std::endl;
    }

    auto centers = find_active_cubes(grid, isovalue);

    std::cout << "Centers of active cubes:\n";
    if (centers.empty()) {
        std::cout << "No active cubes found.\n";
    } else {
        for (auto& center : centers) {
            ++count;
            //std::cout << "(" << std::get<0>(center) << " " << std::get<1>(center) << " " << std::get<2>(center) << ")" << std::endl;
        }
        std::cout << "Total number of centers of active cubes found: " << count << std::endl;
    }

    return 0;
}
