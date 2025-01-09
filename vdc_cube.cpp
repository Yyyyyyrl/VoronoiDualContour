#include "vdc_cube.h"

//! @brief Checks if two cubes are adjacent.
bool is_adjacent(const Cube &cubeA, const Cube &cubeB)
{
    return (std::abs(cubeA.repVertex.x() - cubeB.repVertex.x()) <= 1 &&
            std::abs(cubeA.repVertex.y() - cubeB.repVertex.y()) <= 1 &&
            std::abs(cubeA.repVertex.z() - cubeB.repVertex.z()) <= 1);
}

//! @brief Calculates the unique index of a cube in a 3D grid.
int get_cube_index(const Point &repVertex, int nx, int ny) {
    return repVertex.z() * nx * ny + repVertex.y() * nx + repVertex.x();
}

//! @brief Finds the indices of neighboring cubes in a 3D grid.
std::vector<int> find_neighbor_indices(const Point3& repVertex, int nx, int ny) {
    std::vector<int> neighbors;
    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dz = -1; dz <= 1; ++dz) {
                if (dx != 0 || dy != 0 || dz != 0) {
                    Point3 neighbor(repVertex.x() + dx, repVertex.y() + dy, repVertex.z() + dz);
                    neighbors.push_back(get_cube_index(neighbor, nx, ny));
                }
            }
        }
    }
    return neighbors;
}

//! @brief Retrieves the centers of a list of cubes.
std::vector<Point> get_cube_centers(const std::vector<Cube> &cubes)
{
    std::vector<Point> centers;
    for (auto &cube : cubes)
    {
        centers.push_back(cube.center);
    }
    return centers;
}

//! @brief Separates active cubes using a greedy approach.
std::vector<Cube> separate_active_cubes_greedy(std::vector<Cube>& activeCubes, int nx, int ny, int nz) {
    std::unordered_map<int, Cube> indexToCubeMap;
    std::vector<Cube> separatedCubes;

    for (const Cube &cube : activeCubes) {
        int index = get_cube_index(cube.repVertex, nx, ny);

        bool isAdj = false;
        for (int dz = -1; dz <= 1 && !isAdj; ++dz) {
            for (int dy = -1; dy <= 1 && !isAdj; ++dy) {
                for (int dx = -1; dx <= 1 && !isAdj; ++dx) {
                    int neighborIndex = index + dz * nx * ny + dy * nx + dx;
                    if (indexToCubeMap.find(neighborIndex) != indexToCubeMap.end()) {
                        isAdj = true;
                    }
                }
            }
        }

        if (!isAdj) {
            separatedCubes.push_back(cube);
            indexToCubeMap[index] = cube;
        }
    }

    return separatedCubes;
}


//! @brief Separates active cubes using a graph-based approach.
std::vector<Cube> separate_active_cubes_graph(std::vector<Cube> &activeCubes) {

    std::vector<Cube> separatedCubes;

    int n = activeCubes.size();
    std::vector<std::vector<int>> adjList(n);


    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (is_adjacent(activeCubes[i], activeCubes[j])) {
                adjList[i].push_back(j);
                adjList[j].push_back(i);
            }
        }
    }


    std::vector<int> color(n,-1);
    std::vector<bool> available(n, true);

    color[0] = 0;

    for (int k = 1; k < n; ++k) {
        // Mark the colors of all adjacent vertices as unavailable
        for (int item : adjList[k]) {
            if (color[item] != -1) {
                available[color[item]] = false;
            }
        }
        
        // Find the first available color
        int cr;
        for (cr = 0; cr < n; ++cr) {
            if (available[cr]) {
                break;
            }
        }
        
        // Assign the found color
        color[k] = cr;
        
        // Reset the values back to true for the next iteration
        for (int adj : adjList[k]) {
            if (color[adj] != -1) {
                available[color[adj]] = true;
            }
        }
    }


    std::unordered_map<int, std::vector<Cube>> colorClasses;
    for (int i = 0; i < n; ++i) {
        colorClasses[color[i]].push_back(activeCubes[i]);
    }


    for (const auto& entry : colorClasses) {
        if (entry.second.size() > separatedCubes.size()) {
            separatedCubes = entry.second;
        }
    }


    return separatedCubes;

}
