#ifndef DMR_IO_H
#define DMR_IO_H


#include "dmr_utilities.h"


// Processing Active Cubes

// Functions for writing output mesh

void writeOFF(const std::string &filename, const std::vector<Point> &vertices, const std::vector<DelaunayTriangle> &triangles, std::map<Point, int> &pointIndexMap);
void writePLY(const std::string &filename, const std::vector<Point> &vertices, const std::vector<DelaunayTriangle> &triangles, std::map<Point, int> &pointIndexMap);

void export_voronoi_to_csv(const std::vector<Point> &voronoi_vertices, const std::vector<Object> &voronoi_edges, const std::string &filename);
#endif