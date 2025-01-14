//! @file vdc_commandline.h
//! @brief Header file for parsing command-line arguments and providing help messages.

#ifndef VDC_COMMANDLINE_H
#define VDC_COMMANDLINE_H

#include "vdc_type.h"
#include "vdc_utilities.h"

/*!
 * @brief Holds top-level parameters (from main or command-line arguments).
 * Encapsulates what used to be global flags and file paths.
 */
struct VDC_PARAM {
    std::string file_path;
    float isovalue;
    std::string output_format;
    std::string output_filename;
    std::string out_csv_name;

    bool out_csv;
    bool sep_isov;
    bool multi_isov;
    bool supersample;
    bool add_bounding_cells;

    int supersample_r;

    // Constructor with some defaults
    VDC_PARAM()
        : file_path(""),
          isovalue(0.0f),
          output_format("off"),
          output_filename("output.off"),
          out_csv_name("voronoi.csv"),
          out_csv(false),
          sep_isov(false),
          multi_isov(false),
          supersample(false),
          add_bounding_cells(false),
          supersample_r(1)
    {}
};

//! @brief Prints the help message.
/*!
 * This function outputs the usage information and options available for the program.
 * It includes details about input/output files, processing modes, and available options.
 */
void print_help();

//! @brief Parses the command-line arguments.
/*!
 * This function processes command-line arguments to configure the program's behavior.
 * It supports various options such as output format, supersampling, and multi/single isosurface modes.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line argument strings.
 */
void parse_arguments(int argc, char *argv[], VDC_PARAM &vp);

#endif // VDC_COMMANDLINE_H
