#ifndef VDC_DEBUG_H
#define VDC_DEBUG_H

#include "vdc_type.h"
#include "vdc_grid.h"

extern bool debug;
extern bool indicator;

void print_facet(Facet f);
void print_cell(Delaunay::Cell c);
void write_dummy_points(Grid &grid, std::vector<Point> dummy_points);

#endif