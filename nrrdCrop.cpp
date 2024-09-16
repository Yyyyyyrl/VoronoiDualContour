#include <iostream>
#include <teem/nrrd.h>

int main(int argc, char** argv) {
    if (argc < 9) {
        std::cerr << "Usage: " << argv[0] << " input.nrrd output.nrrd minX minY minZ maxX maxY maxZ" << std::endl;
        return 1;
    }

    // Input arguments
    const char* inputFileName = argv[1];
    const char* outputFileName = argv[2];
    int minX = std::stoi(argv[3]);
    int minY = std::stoi(argv[4]);
    int minZ = std::stoi(argv[5]);
    int maxX = std::stoi(argv[6]);
    int maxY = std::stoi(argv[7]);
    int maxZ = std::stoi(argv[8]);

    // Load the NRRD file
    Nrrd* nrrd = nrrdNew();
    if (nrrdLoad(nrrd, inputFileName, NULL)) {
        std::cerr << "Error loading NRRD file: " << inputFileName << std::endl;
        return 1;
    }

    // Ensure the file is 3D
    if (nrrd->dim != 3) {
        std::cerr << "This program only supports 3D NRRD files." << std::endl;
        return 1;
    }

    // Get the size of the NRRD file's axes
    size_t nx = nrrd->axis[0].size;
    size_t ny = nrrd->axis[1].size;
    size_t nz = nrrd->axis[2].size;

    // Check bounds
    if (minX < 0 || minY < 0 || minZ < 0 || maxX >= nx || maxY >= ny || maxZ >= nz) {
        std::cerr << "Cropping region exceeds NRRD bounds!" << std::endl;
        return 1;
    }

    // Determine the data type and cast accordingly
    switch (nrrd->type) {
        case nrrdTypeFloat: {
            float* data = static_cast<float*>(nrrd->data);
            // Iterate through the NRRD data and set values outside the region to 0
            for (size_t z = 0; z < nz; ++z) {
                for (size_t y = 0; y < ny; ++y) {
                    for (size_t x = 0; x < nx; ++x) {
                        if (x < static_cast<size_t>(minX) || x > static_cast<size_t>(maxX) ||
                            y < static_cast<size_t>(minY) || y > static_cast<size_t>(maxY) ||
                            z < static_cast<size_t>(minZ) || z > static_cast<size_t>(maxZ)) {
                            data[z * nx * ny + y * nx + x] = 0.0f;  // Set value to 0 if outside the region
                        }
                    }
                }
            }
            break;
        }
        case nrrdTypeUChar: {
            unsigned char* data = static_cast<unsigned char*>(nrrd->data);
            // Iterate through the NRRD data and set values outside the region to 0
            for (size_t z = 0; z < nz; ++z) {
                for (size_t y = 0; y < ny; ++y) {
                    for (size_t x = 0; x < nx; ++x) {
                        if (x < static_cast<size_t>(minX) || x > static_cast<size_t>(maxX) ||
                            y < static_cast<size_t>(minY) || y > static_cast<size_t>(maxY) ||
                            z < static_cast<size_t>(minZ) || z > static_cast<size_t>(maxZ)) {
                            data[z * nx * ny + y * nx + x] = 0;  // Set value to 0 if outside the region
                        }
                    }
                }
            }
            break;
        }
        default:
            std::cerr << "Unsupported NRRD data type. Only float and unsigned char are supported." << std::endl;
            return 1;
    }

    // Save the modified NRRD file
    if (nrrdSave(outputFileName, nrrd, NULL)) {
        std::cerr << "Error saving NRRD file: " << outputFileName << std::endl;
        return 1;
    }

    // Clean up
    nrrdNuke(nrrd);

    std::cout << "NRRD cropping complete. Output saved to " << outputFileName << std::endl;

    return 0;
}
