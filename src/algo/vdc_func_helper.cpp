//! @file vdc_func.cpp
//! @brief Implementation of functions for Voronoi Diagram and Isosurface computation.

#include "algo/vdc_func.h"

// Helper for positive mod
static int positive_mod(int val, int mod)
{
    int res = val % mod;
    return (res < 0) ? res + mod : res;
}

// Helper function to check if two directions are approximately equal
bool directions_equal(const Vector3 &d1, const Vector3 &d2, double epsilon)
{
    Vector3 n1 = d1 / std::sqrt(d1.squared_length());      // Normalize d1
    Vector3 n2 = d2 / std::sqrt(d2.squared_length());      // Normalize d2
    return (n1 - n2).squared_length() < epsilon * epsilon; // Compare squared distance of normalized vectors
}

//! @brief Handles output mesh generation.
/*!
 * Writes the final isosurface mesh to file in the specified format (OFF or PLY).
 * Supports both single and multi-isovertex modes.
 *
 * @param retFlag Output parameter indicating success/failure
 * @param vd The Voronoi diagram
 * @param vdc_param Configuration parameters
 * @param iso_surface The isosurface to output
 * @return EXIT_SUCCESS on success, EXIT_FAILURE on error
 */
int handle_output_mesh(bool &retFlag, VoronoiDiagram &vd, VDC_PARAM &vdc_param, IsoSurface &iso_surface)
{
    retFlag = true;

    std::cout << "Result file at: " << vdc_param.output_filename << std::endl;

    // Multi-isovertex mode output
    if (vdc_param.multi_isov)
    {
        if (vdc_param.output_format == "off")
        {
            writeOFFMulti(vdc_param.output_filename, vd, iso_surface.isosurfaceTrianglesMulti, iso_surface);
        }
        else if (vdc_param.output_format == "ply")
        {
            writePLYMulti(vdc_param.output_filename, vd, iso_surface.isosurfaceTrianglesMulti, iso_surface);
        }
        else
        {
            std::cerr << "Unsupported output format: " << vdc_param.output_format << std::endl;
            return EXIT_FAILURE;
        }
    }
    // Single-isovertex mode output
    else
    {
        if (vdc_param.output_format == "off")
        {
            writeOFFSingle(vdc_param.output_filename, iso_surface.isosurfaceVertices, iso_surface.isosurfaceTrianglesSingle);
        }
        else if (vdc_param.output_format == "ply")
        {
            writePLYSingle(vdc_param.output_filename, iso_surface.isosurfaceVertices, iso_surface.isosurfaceTrianglesSingle);
        }
        else
        {
            std::cerr << "Unsupported output format: " << vdc_param.output_format << std::endl;
            return EXIT_FAILURE;
        }
    }

    retFlag = false;
    return EXIT_SUCCESS;
}
