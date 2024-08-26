#include "io.h"


template <typename T>
std::vector<float> convert_to_float_vector(T *data_ptr, size_t total_size)
{
    std::vector<float> data(total_size);
    for (size_t i = 0; i < total_size; ++i)
    {
        data[i] = static_cast<float>(data_ptr[i]);
    }
    return data;
}


Grid load_nrrd_data(const std::string &file_path)
{
    Nrrd *nrrd = nrrdNew();
    if (nrrdLoad(nrrd, file_path.c_str(), NULL))
    {
        char *err = biffGetDone(NRRD);
        std::cerr << "Error reading NRRD file: " << err << std::endl;
        free(err);
        nrrdNuke(nrrd);
        exit(1);
    }

    size_t total_size = nrrdElementNumber(nrrd);

    std::vector<float> data;

    if (nrrd->type == nrrdTypeFloat)
    {
        float *data_ptr = static_cast<float *>(nrrd->data);
        data = std::vector<float>(data_ptr, data_ptr + total_size);
    }
    else if (nrrd->type == nrrdTypeUChar)
    {
        unsigned char *data_ptr = static_cast<unsigned char *>(nrrd->data);
        data = convert_to_float_vector(data_ptr, total_size);
    }
    else
    {
        std::cerr << "Unsupported NRRD data type." << std::endl;
        nrrdNuke(nrrd);
        exit(1);
    }

    int nx = nrrd->axis[0].size;
    int ny = nrrd->axis[1].size;
    int nz = nrrd->axis[2].size;
    float dx = nrrd->axis[0].spacing;
    float dy = nrrd->axis[1].spacing;
    float dz = nrrd->axis[2].spacing;


    nrrdNuke(nrrd); // Properly dispose of the Nrrd structure

    return {data, nx, ny, nz, dx, dy, dz};
}


void writeOFF(const std::string &filename, const std::vector<Point> &vertices, const std::vector<DelaunayTriangle> &triangles, std::map<Point, int> &pointIndexMap)
{
    std::ofstream out(filename);
    if (!out)
    {
        std::cerr << "Cannot open file for writing: " << filename << std::endl;
        return;
    }

    // Write the header
    out << "OFF\n";
    out << vertices.size() << " " << triangles.size() << " 0\n";

    // Write the vertices
    for (const auto &vertex : vertices)
    {
        out << vertex.x() << " " << vertex.y() << " " << vertex.z() << "\n";
    }

    // Write the faces
    for (const auto &triangle : triangles)
    {
        out << "3 " << pointIndexMap[triangle.vertex1] << " " << pointIndexMap[triangle.vertex2] << " " << pointIndexMap[triangle.vertex3] << "\n";
    }

    out.close();
}

void writePLY(const std::string &filename, const std::vector<Point> &vertices, const std::vector<DelaunayTriangle> &triangles, std::map<Point, int> &pointIndexMap)
{
    std::ofstream out(filename);
    if (!out)
    {
        std::cerr << "Cannot open file for writing: " << filename << std::endl;
        return;
    }

    // Write PLY header
    out << "ply\n";
    out << "format ascii 1.0\n";
    out << "element vertex " << vertices.size() << "\n";
    out << "property float x\n";
    out << "property float y\n";
    out << "property float z\n";
    out << "element face " << triangles.size() << "\n";
    out << "property list uchar int vertex_index\n";
    out << "end_header\n";

    // Write vertices
    for (const auto &vertex : vertices)
    {
        out << vertex.x() << " " << vertex.y() << " " << vertex.z() << "\n";
    }

    // Write faces
    for (const auto &triangle : triangles)
    {
        out << "3 " << pointIndexMap[triangle.vertex1] << " " << pointIndexMap[triangle.vertex2] << " " << pointIndexMap[triangle.vertex3] << "\n";
    }

    out.close();
}