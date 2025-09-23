//! @file vdc_commandline.h
//! @brief Header file for parsing command-line arguments and providing help messages.

#ifndef VDC_COMMANDLINE_H
#define VDC_COMMANDLINE_H

#include "core/vdc_type.h"
#include "core/vdc_utilities.h"

//! @brief Structure to hold top-level parameters parsed from command-line arguments.
/*!
 * This structure consolidates all configurable parameters for the program,
 * replacing the need for multiple global variables.
 */
struct VDC_PARAM {
    std::string file_path;         //!< Path to the input raw data file (nhdr/nrrd format).
    float isovalue;                //!< The isovalue used for isosurface extraction.
    std::string output_format;     //!< The format of the output file ("off" or "ply").
    std::string output_filename;   //!< The name of the output file.
    std::string out_csv_name;      //!< The name of the CSV file for Voronoi diagram export.
    
    bool out_csv;                  //!< Flag to enable exporting Voronoi diagram to CSV.
    bool sep_isov;                 //!< Flag to enable separation of non-adjacent active cubes.
    bool multi_isov;               //!< Flag to enable multi-isosurface mode.
    bool supersample;              //!< Flag to enable supersampling of the input data.
    bool add_bounding_cells;       //!< Flag to include bounding cells in the Voronoi diagram.
    bool convex_hull;              //!< Flag to enable convex hull computation in building voronoi cells
    bool test_vor = false;         //!< Flag for testing the Voronoi diagram construction
    bool mod_cyc = false;          //!< Guard: run modify-cycles pass (facet rematching + cycle recompute)
    bool summary_stats = false;    //!< Guard: print summary statistics at the end of the run

    int supersample_r;             //!< Factor by which the input data is supersampled.
    double collapse_eps;           //!< Absolute collapse threshold in world units (optional).

    //! @brief Constructor to initialize default parameter values.
    VDC_PARAM()
        : file_path(""),
          isovalue(0.0f),
          output_format("off"),
          output_filename(""),
          out_csv_name("voronoi.csv"),
          out_csv(false),
          sep_isov(false),
          multi_isov(false),
          supersample(false),
          add_bounding_cells(false),
          convex_hull(false),
          supersample_r(1),
          collapse_eps(-1.0),
          mod_cyc(false),
          summary_stats(false)
    {}
};

//! @brief Prints the help message to the console.
/*!
 * This function outputs usage information and available options for the program,
 * including details about input/output configurations and processing modes.
 */
void print_help();

//! @brief Parses the command-line arguments to populate program parameters.
/*!
 * This function processes command-line arguments to configure the program's behavior.
 * It sets parameters such as output format, supersampling factor, and mode selection.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line argument strings.
 * @param vp A reference to a `VDC_PARAM` object where parsed parameters are stored.
 */
void parse_arguments(int argc, char *argv[], VDC_PARAM &vp);

#endif // VDC_COMMANDLINE_H
