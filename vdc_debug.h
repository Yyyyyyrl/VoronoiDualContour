#ifndef VDC_DEBUG_H
#define VDC_DEBUG_H

#include "vdc_type.h"

extern bool debug;
extern bool indicator;

void print_facet(Facet f);
void print_cell(Delaunay::Cell c);

#endif