#include "dmr_voronoi.h"

void Cycle::compute_centroid() {
    isovertex = CGAL::centroid(midpoints.begin(), midpoints.end());
}