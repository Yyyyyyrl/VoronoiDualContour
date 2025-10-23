// Temporary debug utilities for analyzing thin triangle issues
#include "processing/vdc_sep_isov.h"
#include <cmath>
#include <iostream>

// Compute angle at vertex b in triangle (a, b, c) in degrees
static double compute_angle_at_vertex(const Point& a, const Point& b, const Point& c)
{
    // Vectors from b to a and b to c
    double ba_x = a.x() - b.x();
    double ba_y = a.y() - b.y();
    double ba_z = a.z() - b.z();

    double bc_x = c.x() - b.x();
    double bc_y = c.y() - b.y();
    double bc_z = c.z() - b.z();

    // Dot product and magnitudes
    double dot = ba_x*bc_x + ba_y*bc_y + ba_z*bc_z;
    double mag_ba = std::sqrt(ba_x*ba_x + ba_y*ba_y + ba_z*ba_z);
    double mag_bc = std::sqrt(bc_x*bc_x + bc_y*bc_y + bc_z*bc_z);

    if (mag_ba < 1e-10 || mag_bc < 1e-10) return 0.0;

    double cos_angle = dot / (mag_ba * mag_bc);
    cos_angle = std::max(-1.0, std::min(1.0, cos_angle)); // Clamp for numerical stability

    return std::acos(cos_angle) * 180.0 / M_PI;
}

// Check if three points form a thin triangle (any angle < threshold)
void check_for_thin_triangle(const Point& p1, const Point& p2, const Point& p3,
                             double angle_threshold, const char* context)
{
    double angle1 = compute_angle_at_vertex(p2, p1, p3);
    double angle2 = compute_angle_at_vertex(p1, p2, p3);
    double angle3 = compute_angle_at_vertex(p1, p3, p2);

    double min_angle = std::min({angle1, angle2, angle3});

    if (min_angle < angle_threshold) {
        std::cerr << "[THIN_TRIANGLE_WARNING] " << context << "\n";
        std::cerr << "  Point 1: (" << p1.x() << ", " << p1.y() << ", " << p1.z() << ")\n";
        std::cerr << "  Point 2: (" << p2.x() << ", " << p2.y() << ", " << p2.z() << ")\n";
        std::cerr << "  Point 3: (" << p3.x() << ", " << p3.y() << ", " << p3.z() << ")\n";
        std::cerr << "  Angles: " << angle1 << "°, " << angle2 << "°, " << angle3 << "°\n";
        std::cerr << "  Min angle: " << min_angle << "°\n\n";
    }
}

// Check all triples of iso-crossing points for near-collinearity
void check_collinear_isocrossings(const std::vector<Cube>& activeCubes,
                                  const UnifiedGrid& grid,
                                  float isovalue,
                                  double angle_threshold)
{
    std::cerr << "\n=== Checking for nearly collinear iso-crossing points ===\n";
    std::cerr << "Total active cubes: " << activeCubes.size() << "\n";
    std::cerr << "Angle threshold: " << angle_threshold << "°\n\n";

    int count = 0;

    // Only check a reasonable number of triples to avoid O(n³) explosion
    const int max_check = std::min(100, (int)activeCubes.size());

    for (int i = 0; i < max_check && i < activeCubes.size(); i++) {
        Point p1 = compute_iso_crossing_point_accurate(grid,
                                                       activeCubes[i].i,
                                                       activeCubes[i].j,
                                                       activeCubes[i].k,
                                                       isovalue);

        for (int j = i+1; j < max_check && j < activeCubes.size(); j++) {
            Point p2 = compute_iso_crossing_point_accurate(grid,
                                                           activeCubes[j].i,
                                                           activeCubes[j].j,
                                                           activeCubes[j].k,
                                                           isovalue);

            // Only check nearby cubes (within 5 grid cells)
            int di = std::abs(activeCubes[i].i - activeCubes[j].i);
            int dj = std::abs(activeCubes[i].j - activeCubes[j].j);
            int dk = std::abs(activeCubes[i].k - activeCubes[j].k);

            if (di > 5 || dj > 5 || dk > 5) continue;

            for (int k = j+1; k < max_check && k < activeCubes.size(); k++) {
                Point p3 = compute_iso_crossing_point_accurate(grid,
                                                               activeCubes[k].i,
                                                               activeCubes[k].j,
                                                               activeCubes[k].k,
                                                               isovalue);

                // Only check nearby cubes
                int dik = std::abs(activeCubes[i].i - activeCubes[k].i);
                int djk = std::abs(activeCubes[i].j - activeCubes[k].j);
                int dkk = std::abs(activeCubes[i].k - activeCubes[k].k);

                if (dik > 5 || djk > 5 || dkk > 5) continue;

                char context[256];
                snprintf(context, sizeof(context),
                        "Nearly collinear iso-crossings from cubes (%d,%d,%d), (%d,%d,%d), (%d,%d,%d)",
                        activeCubes[i].i, activeCubes[i].j, activeCubes[i].k,
                        activeCubes[j].i, activeCubes[j].j, activeCubes[j].k,
                        activeCubes[k].i, activeCubes[k].j, activeCubes[k].k);

                check_for_thin_triangle(p1, p2, p3, angle_threshold, context);
                count++;
            }
        }
    }

    std::cerr << "Checked " << count << " triples of nearby iso-crossing points\n";
    std::cerr << "=== End collinearity check ===\n\n";
}
