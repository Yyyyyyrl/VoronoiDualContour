#include "vdc_voronoi.h"

//! @brief Computes the centroid of a cycle using the associated midpoints.
/*!
 * The centroid is calculated using CGAL's `centroid` function, which averages
 * the positions of the midpoints in the cycle. If no midpoints exist, the centroid
 * computation is skipped.
 *
 * @param midpoints Vector of all midpoints in the Voronoi diagram.
 */
void Cycle::compute_centroid(const std::vector<MidpointNode>& midpoints)
{
    if (midpoint_indices.empty())
    {
        return; // No midpoints to compute centroid.
    }

    // Collect the points corresponding to the midpoint indices.
    std::vector<Point> points;
    for (int idx : midpoint_indices)
    {
        points.push_back(midpoints[idx].point);
    }

    // Compute the centroid using CGAL's centroid function.
    isovertex = CGAL::centroid(points.begin(), points.end());
}