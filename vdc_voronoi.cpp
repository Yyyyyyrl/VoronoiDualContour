#include "vdc_voronoi.h"

void Cycle::compute_centroid(const std::vector<MidpointNode>& midpoints)
{
    if (midpoint_indices.empty())
    {
        return; // No midpoints to compute centroid
    }

    // Collect the points corresponding to the midpoint indices
    std::vector<Point> points;
    for (int idx : midpoint_indices)
    {
        points.push_back(midpoints[idx].point);
    }

    // Compute the centroid using CGAL's centroid function
    isovertex = CGAL::centroid(points.begin(), points.end());
}
