#include "vdc_commandline.h"

void print_help()
{
    std::cout << "Usage: dmr [OPTIONS] <isovalue> <(nhdr/nrrd) raw data file path>\n\n";
    std::cout << "OPTIONS:\n";
    std::cout << "  -o {output_filename}        : Specify output filename (default: derived from input filename).\n";
    std::cout << "  -off                        : Generate output in .off format (default).\n";
    std::cout << "  -ply                        : Generate output in .ply format.\n";
    std::cout << "  -out_csv {output_csv_name}  : Write the Voronoi diagram to a CSV file.\n";
    std::cout << "  -sep_isov                   : Pick a subset of non-adjacent active cubes of the input data before constructing triangulation.\n";
    std::cout << "  -supersample {factor}       : Supersample the input data by the given factor.\n";
    std::cout << "  -multi_isov                 : Use multi iso-vertices mode (default).\n";
    std::cout << "  -single_isov                : Use single iso-vertices mode.\n";
    std::cout << "  --help                      : Print this help message.\n";
}

void parse_arguments(int argc, char *argv[])
{
    if (argc < 3)
    {
        print_help();
        exit(EXIT_FAILURE);
    }

    // Parse options
    int i = 1;
    while (i < argc && argv[i][0] == '-')
    {
        std::string arg = argv[i];

        if (arg == "-o" && i + 1 < argc)
        {
            output_filename = argv[++i];
        }
        else if (arg == "-off")
        {
            output_format = "off";
        }
        else if (arg == "-ply")
        {
            output_format = "ply";
        }
        else if (arg == "-out_csv" && i + 1 < argc)
        {
            out_csv = true;
            out_csv_name = argv[++i];
        }
        else if (arg == "-sep_isov")
        {
            sep_isov = true;
        }
        else if (arg == "-supersample" && i + 1 < argc)
        {
            supersample = true;
            supersample_r = std::atoi(argv[++i]);
        }
        else if (arg == "-multi_isov")
        {
            multi_isov = true;
        }
        else if (arg == "-single_isov")
        {
            multi_isov = false;
        }
        else if (arg == "--help")
        {
            print_help();
            exit(EXIT_SUCCESS);
        }
        else
        {
            std::cerr << "Unknown option: " << arg << std::endl;
            print_help();
            exit(EXIT_FAILURE);
        }
        ++i;
    }

    // Parse required arguments
    if (i + 2 > argc)
    {
        std::cerr << "Error: Missing required arguments.\n";
        print_help();
        exit(EXIT_FAILURE);
    }

    isovalue = std::atof(argv[i++]);
    file_path = argv[i++];

    // Default output filename if not specified
    if (output_filename.empty())
    {
        std::string base_name = file_path.substr(0, file_path.find_last_of('.'));
        output_filename = base_name;

        if (multi_isov) {
            output_filename += "_multi-isov";
        } else {
            output_filename += "_single-isov";
        }

        if (supersample) {
            output_filename += "_sup-" + std::to_string(supersample_r);
        }
        
        if (sep_isov) {
            output_filename += "_sep-isov";
        }

        output_filename += "." + output_format;

    }
}