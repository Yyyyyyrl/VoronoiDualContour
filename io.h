#ifndef IO_H
#define IO_H

#include <iostream>
#include <vector>
#include <teem/nrrd.h>

#include "dt.h"
#include "utilities.h"

// Functions for loading nrrd data
template <typename T>
std::vector<float> convert_to_float_vector(T *data_ptr, size_t total_size);

Grid load_nrrd_data(const std::string &file_path);

// Processing Active Cubes

// Functions for writing output mesh

void writeOFF(const std::string &filename, const std::vector<Point> &vertices, const std::vector<DelaunayTriangle> &triangles, std::map<Point, int> &pointIndexMap);
void writePLY(const std::string &filename, const std::vector<Point> &vertices, const std::vector<DelaunayTriangle> &triangles, std::map<Point, int> &pointIndexMap);


#endif IO_H



