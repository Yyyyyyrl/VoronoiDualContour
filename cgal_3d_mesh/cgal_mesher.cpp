/**
 * @file cgal_mesher.cpp
 * @brief CGAL 3D Mesh Generation wrapper for NRRD volumetric data
 *
 * This program provides a comparison baseline to the VDC (Voronoi-based Dual Contouring)
 * method by implementing CGAL's modern 3D mesh generation approach on the same input data.
 *
 * Input: NRRD volumetric scalar field
 * Output: OFF triangle mesh (surface only)
 *
 * Uses: CGAL's Mesh_3 package (modern, non-deprecated)
 */

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>

#include <teem/nrrd.h>

#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <string>
#include <vector>
#include <algorithm>

// ============================================================================
// CGAL Type Definitions
// ============================================================================

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Labeled_mesh_domain_3<K> Mesh_domain;
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

namespace params = CGAL::parameters;

// ============================================================================
// NRRD to CGAL::Image_3 Converter
// ============================================================================

/**
 * @brief Convert NRRD data to CGAL::Image_3 format
 *
 * CGAL::Image_3 requires data in a specific format for meshing.
 * This function handles the conversion from NRRD to CGAL's image format.
 */
class NRRDImageConverter {
public:
    /**
     * @brief Load NRRD and convert to CGAL::Image_3
     */
    static CGAL::Image_3 load(const std::string& filename) {
        Nrrd* nrrd = nrrdNew();

        if (nrrdLoad(nrrd, filename.c_str(), NULL)) {
            char* err = biffGetDone(NRRD);
            std::cerr << "Error reading NRRD file: " << err << std::endl;
            free(err);
            nrrdNuke(nrrd);
            throw std::runtime_error("Failed to load NRRD file");
        }

        // Get dimensions
        if (nrrd->dim != 3) {
            std::cerr << "Error: Expected 3D NRRD data, got " << nrrd->dim << "D" << std::endl;
            nrrdNuke(nrrd);
            throw std::runtime_error("Invalid NRRD dimensions");
        }

        unsigned int nx = nrrd->axis[0].size;
        unsigned int ny = nrrd->axis[1].size;
        unsigned int nz = nrrd->axis[2].size;

        // Get spacing (default to 1.0 if not specified)
        double dx = nrrd->axis[0].spacing > 0 ? nrrd->axis[0].spacing : 1.0;
        double dy = nrrd->axis[1].spacing > 0 ? nrrd->axis[1].spacing : 1.0;
        double dz = nrrd->axis[2].spacing > 0 ? nrrd->axis[2].spacing : 1.0;

        std::cout << "[INFO] Loaded NRRD: " << nx << " x " << ny << " x " << nz << std::endl;
        std::cout << "[INFO] Spacing: " << dx << " x " << dy << " x " << dz << std::endl;

        // Convert NRRD data to float array
        float* float_data = new float[nx * ny * nz];

        if (nrrd->type == nrrdTypeFloat) {
            float* src = static_cast<float*>(nrrd->data);
            std::copy(src, src + nx * ny * nz, float_data);
        } else if (nrrd->type == nrrdTypeDouble) {
            double* src = static_cast<double*>(nrrd->data);
            for (size_t i = 0; i < nx * ny * nz; ++i) {
                float_data[i] = static_cast<float>(src[i]);
            }
        } else if (nrrd->type == nrrdTypeUChar) {
            unsigned char* src = static_cast<unsigned char*>(nrrd->data);
            for (size_t i = 0; i < nx * ny * nz; ++i) {
                float_data[i] = static_cast<float>(src[i]);
            }
        } else if (nrrd->type == nrrdTypeShort) {
            short* src = static_cast<short*>(nrrd->data);
            for (size_t i = 0; i < nx * ny * nz; ++i) {
                float_data[i] = static_cast<float>(src[i]);
            }
        } else {
            std::cerr << "Warning: Unsupported NRRD type " << nrrd->type
                      << ", attempting conversion..." << std::endl;
            // Try to convert through float
            Nrrd* nfloat = nrrdNew();
            if (nrrdConvert(nfloat, nrrd, nrrdTypeFloat)) {
                char* err = biffGetDone(NRRD);
                std::cerr << "Error converting NRRD: " << err << std::endl;
                free(err);
                nrrdNuke(nfloat);
                nrrdNuke(nrrd);
                delete[] float_data;
                throw std::runtime_error("Failed to convert NRRD data");
            }
            float* src = static_cast<float*>(nfloat->data);
            std::copy(src, src + nx * ny * nz, float_data);
            nrrdNuke(nfloat);
        }

        // Find data range for info
        float min_val = *std::min_element(float_data, float_data + nx * ny * nz);
        float max_val = *std::max_element(float_data, float_data + nx * ny * nz);
        std::cout << "[INFO] Data range: [" << min_val << ", " << max_val << "]" << std::endl;

        nrrdNuke(nrrd);

        // Create CGAL::Image_3
        // CGAL::Image_3 expects data in a specific format
        // The data should be in the order: x varies fastest, then y, then z
        CGAL::Image_3 image;

        // Create the image with the data
        // _image is the internal pointer that Image_3 uses
        ::_image* cgal_image = ::_initImage();
        cgal_image->xdim = nx;
        cgal_image->ydim = ny;
        cgal_image->zdim = nz;
        cgal_image->vdim = 1;  // Scalar field

        cgal_image->vx = dx;
        cgal_image->vy = dy;
        cgal_image->vz = dz;

        cgal_image->tx = 0.0;
        cgal_image->ty = 0.0;
        cgal_image->tz = 0.0;

        cgal_image->wordKind = WK_FLOAT;
        cgal_image->sign = SGN_SIGNED;
        cgal_image->wdim = sizeof(float);

        cgal_image->data = float_data;
        cgal_image->imageFormat = nullptr;

        image = CGAL::Image_3(cgal_image);

        return image;
    }
};

// ============================================================================
// Helper Functions
// ============================================================================

void print_help() {
    std::cout << "CGAL Mesh Generator - Using CGAL's 3D Mesh Generation package\n\n";
    std::cout << "Usage: cgal_mesher [OPTIONS] <isovalue> <nrrd_file>\n\n";
    std::cout << "Arguments:\n";
    std::cout << "  isovalue              Isovalue for surface extraction\n";
    std::cout << "  nrrd_file             Input NRRD volumetric data file\n\n";
    std::cout << "Options:\n";
    std::cout << "  -o <filename>         Output filename (default: cgal_output.off)\n";
    std::cout << "  -facet_angle <deg>    Minimum facet angle (default: 25.0, max: 30.0)\n";
    std::cout << "  -facet_size <value>   Maximum facet size (default: auto)\n";
    std::cout << "  -facet_distance <val> Maximum distance to surface (default: auto)\n";
    std::cout << "  -cell_ratio <value>   Cell radius-edge ratio (default: 3.0)\n";
    std::cout << "  -cell_size <value>    Maximum cell size (default: auto)\n";
    std::cout << "  -surface_only         Generate surface mesh only (no volume)\n";
    std::cout << "  -h, --help            Show this help message\n\n";
}

// ============================================================================
// Main Program
// ============================================================================

int main(int argc, char* argv[]) {
    // Default parameters
    std::string output_file = "cgal_output.off";
    double facet_angle = 20.0;
    double facet_size = 0.7;      
    double facet_distance = 0.6; 
    double cell_ratio = 3.0;
    double cell_size = 4.0;      
    bool surface_only = false;
    double isovalue = 0.0;
    std::string input_file;

    // Parse command-line arguments
    int i = 1;
    while (i < argc) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            print_help();
            return 0;
        } else if (arg == "-o" && i + 1 < argc) {
            output_file = argv[++i];
        } else if (arg == "-facet_angle" && i + 1 < argc) {
            facet_angle = std::atof(argv[++i]);
            if (facet_angle > 30.0) {
                std::cerr << "Warning: Facet angle > 30° may not guarantee termination. Clamping to 30°.\n";
                facet_angle = 30.0;
            }
        } else if (arg == "-facet_size" && i + 1 < argc) {
            facet_size = std::atof(argv[++i]);
        } else if (arg == "-facet_distance" && i + 1 < argc) {
            facet_distance = std::atof(argv[++i]);
        } else if (arg == "-cell_ratio" && i + 1 < argc) {
            cell_ratio = std::atof(argv[++i]);
        } else if (arg == "-cell_size" && i + 1 < argc) {
            cell_size = std::atof(argv[++i]);
        } else if (arg == "-surface_only") {
            surface_only = true;
        } else if (arg[0] != '-') {
            // Positional arguments
            if (isovalue == 0.0 && input_file.empty()) {
                isovalue = std::atof(arg.c_str());
            } else if (input_file.empty()) {
                input_file = arg;
            } else {
                std::cerr << "Error: Too many positional arguments\n";
                print_help();
                return 1;
            }
        } else {
            std::cerr << "Error: Unknown option '" << arg << "'\n";
            print_help();
            return 1;
        }
        ++i;
    }

    if (input_file.empty()) {
        std::cerr << "Error: Missing required arguments\n";
        print_help();
        return 1;
    }

    std::cout << "=================================================================\n";
    std::cout << "CGAL 3D Mesh Generator - VDC Comparison Tool\n";
    std::cout << "=================================================================\n";
    std::cout << "[CONFIG] Input file: " << input_file << "\n";
    std::cout << "[CONFIG] Isovalue: " << isovalue << "\n";
    std::cout << "[CONFIG] Output file: " << output_file << "\n";
    std::cout << "[CONFIG] Facet angle: " << facet_angle << "°\n";
    std::cout << "[CONFIG] Surface only: " << (surface_only ? "yes" : "no") << "\n";

    auto start_time = std::chrono::high_resolution_clock::now();

    try {
        // Load NRRD data and convert to CGAL::Image_3
        std::cout << "\n[1/4] Loading NRRD data and converting to CGAL format...\n";
        CGAL::Image_3 image = NRRDImageConverter::load(input_file);

        // Create mesh domain from image
        std::cout << "[2/4] Creating mesh domain from image...\n";

        // Create domain with iso_value parameter
        // value_outside specifies the value outside the domain
        Mesh_domain domain = Mesh_domain::create_gray_image_mesh_domain(
            image,
            params::iso_value(static_cast<float>(isovalue))
                  .value_outside(0.f)
        );

        // Auto-compute sizing parameters if not specified
        // Base on image spacing (voxel size) for finer control
        double min_spacing = std::min({image.vx(), image.vy(), image.vz()});
        double max_extent = std::max({
            image.xdim() * image.vx(),
            image.ydim() * image.vy(),
            image.zdim() * image.vz()
        });


        std::cout << "[INFO] Min voxel spacing: " << min_spacing << "\n";
        std::cout << "[INFO] Max extent: " << max_extent << "\n";
        std::cout << "[CONFIG] Facet size: " << facet_size << "\n";
        std::cout << "[CONFIG] Facet distance: " << facet_distance << "\n";
        std::cout << "[CONFIG] Cell radius-edge ratio: " << cell_ratio << "\n";
        std::cout << "[CONFIG] Cell size: " << cell_size << "\n";

        // Create mesh criteria
        std::cout << "[3/4] Generating mesh...\n";

        Mesh_criteria criteria(
            params::facet_angle(facet_angle)
                  .facet_size(facet_size)
                  .facet_distance(facet_distance)
                  .cell_radius_edge_ratio(cell_ratio)
                  .cell_size(cell_size)
        );

        // Generate mesh
        auto mesh_start = std::chrono::high_resolution_clock::now();

        C3t3 c3t3;
        if (surface_only) {
            // Surface-only meshing (equivalent to deprecated make_surface_mesh)
            c3t3 = CGAL::make_mesh_3<C3t3>(
                domain,
                criteria,
                params::no_perturb().no_exude()
            );
        } else {
            // Full volume mesh generation
            c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);
        }

        auto mesh_end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> mesh_time = mesh_end - mesh_start;
        std::cout << "[INFO] Mesh generation time: " << mesh_time.count() << " seconds\n";
        std::cout << "[INFO] Number of facets in complex: " << c3t3.number_of_facets_in_complex() << "\n";
        std::cout << "[INFO] Number of cells in complex: " << c3t3.number_of_cells_in_complex() << "\n";

        // Write output (boundary surface only)
        std::cout << "[4/4] Writing output to " << output_file << "...\n";
        std::ofstream out(output_file);
        if (!out) {
            std::cerr << "Error: Cannot open output file\n";
            return 1;
        }

        // Output only the boundary surface to OFF format
        c3t3.output_boundary_to_off(out);
        out.close();

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> total_time = end_time - start_time;

        std::cout << "\n=================================================================\n";
        std::cout << "Total processing time: " << total_time.count() << " seconds\n";
        std::cout << "Output written to: " << output_file << "\n";
        std::cout << "=================================================================\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
