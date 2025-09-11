#include "core/vdc_commandline.h"

//! Prints the help message for the program.
void print_help()
{
    std::cout << "Usage: vdc [OPTIONS] <isovalue> <(nhdr/nrrd) raw data file path>\n\n";
    std::cout << "OPTIONS:\n";
    std::cout << "  -o {output_filename}        : Specify output filename (default: derived from input filename).\n";
    std::cout << "  -off                        : Generate output in .off format (default).\n";
    std::cout << "  -ply                        : Generate output in .ply format.\n";
    std::cout << "  -out_csv {output_csv_name}  : Write the Voronoi diagram to a CSV file.\n";
    std::cout << "  -sep_isov                   : Pick a subset of non-adjacent active cubes of the input data before constructing triangulation.\n";
    std::cout << "  -supersample {factor}       : Supersample the input data by the given factor.\n";
    std::cout << "  -multi_isov                 : Use multi iso-vertices mode.\n";
    std::cout << "  -single_isov                : Use single iso-vertices mode (default).\n";
    std::cout << "  -conv_H                     : Use the Convex_Hull_3 from CGAL in voronoi cell construction.\n";
    std::cout << "  -mod_cyc                    : After initial cycles, try facet rematching and recompute cycles.\n";
    std::cout << "  --help                      : Print this help message.\n";
}

//! Parses command-line arguments and configures program settings.
void parse_arguments(int argc, char *argv[], VDC_PARAM &vp)
{
    // Print help and exit if there are insufficient arguments.
    if (argc < 3)
    {
        print_help();
        exit(EXIT_FAILURE);
    }

    // Parse optional arguments (those starting with '-').
    int i = 1;
    while (i < argc && argv[i][0] == '-')
    {
        std::string arg = argv[i];

        if (arg == "-o" && i + 1 < argc)
        {
            vp.output_filename = argv[++i]; // Set custom output filename.
        }
        else if (arg == "-off")
        {
            vp.output_format = "off"; // Set output format to .off.
        }
        else if (arg == "-ply")
        {
            vp.output_format = "ply"; // Set output format to .ply.
        }
        else if (arg == "-out_csv" && i + 1 < argc)
        {
            vp.out_csv = true;                // Enable CSV output.
            vp.out_csv_name = argv[++i];      // Set CSV output filename
        }
        else if (arg == "-sep_isov")
        {
            vp.sep_isov = true; // Enable separation of non-adjacent active cubes.
        }
        else if (arg == "-supersample" && i + 1 < argc)
        {
            vp.supersample = true;                     // Enable supersampling.
            vp.supersample_r = std::atoi(argv[++i]);   // Set supersampling factor.
        }
        else if (arg == "-multi_isov")
        {
            vp.multi_isov = true; // Enable multi-isovertex mode.
        }
        else if (arg == "-single_isov")
        {
            vp.multi_isov = false; // Enable single-isovertex mode.
        }
        else if (arg == "--help")
        {
            print_help();
            exit(EXIT_SUCCESS);
        }
        else if (arg == "-bound_cells")
        {
            vp.add_bounding_cells = true; // Add bounding cells to the Voronoi diagram.
        }
        else if (arg == "-conv_H")
        {
            vp.convex_hull = true;
        }
        else if (arg == "--test_vor")
        {
            vp.test_vor = true;
        }
        else if (arg == "-mod_cyc")
        {
            vp.mod_cyc = true; // Enable modify-cycles pass (guarded)
        }
        else
        {
            // Handle unknown options.
            std::cerr << "Unknown option: " << arg << std::endl;
            print_help();
            exit(EXIT_FAILURE);
        }
        ++i;
    }

    // Parse required arguments: isovalue and file path.
    if (i + 2 > argc)
    {
        std::cerr << "Error: Missing required arguments.\n";
        print_help();
        exit(EXIT_FAILURE);
    }

    vp.isovalue = std::atof(argv[i++]); // Parse isovalue as a floating-point number.
    vp.file_path = argv[i++];           // Parse the raw data file path.

    // Generate default output filename if not specified.
    if (vp.output_filename.empty())
    {
        // Extract base name from the file path (without extension).
        std::string base_name = vp.file_path.substr(0, vp.file_path.find_last_of('.'));
        vp.output_filename = base_name;

        // Append processing details to the filename.
        if (vp.multi_isov)
        {
            vp.output_filename += "_multi-isov";
        }
        else
        {
            vp.output_filename += "_single-isov";
        }

        if (vp.supersample)
        {
            vp.output_filename += "_sup-" + std::to_string(vp.supersample_r);
        }

        if (vp.sep_isov)
        {
            vp.output_filename += "_sep-isov";
        }

        if (vp.convex_hull)
        {
            vp.output_filename += "_conv-H";
        }

        // Add file format extension.
        vp.output_filename += "." + vp.output_format;
    }
}
