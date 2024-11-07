#ifndef DMR_DEBUG_H
#define DMR_DEBUG_H

#include "dmr_type.h"

extern bool debug;
extern bool indicator;

void print_facet(Facet f);
void print_cell(Delaunay::Cell c);

#endif