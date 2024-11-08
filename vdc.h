#ifndef VDC_H
#define VDC_H

#include "vdc_utilities.h"
#include "vdc_debug.h"
#include "vdc_io.h"
#include "vdc_commandline.h"
#include <cstdlib>
#include <map>

#endif // VDC_H

void Compute_Isosurface_Vertices_Multi(std::vector<VoronoiCell> &voronoi_cells);

void Compute_Isosurface_Vertices_Single(ScalarGrid &grid);
