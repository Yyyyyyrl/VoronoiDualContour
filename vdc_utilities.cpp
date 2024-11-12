#include "vdc_utilities.h"


std::string objectToString(const Object &obj)
{
    std::ostringstream stream;
    Segment3 seg;
    Ray3 ray;
    Line3 line;
    if (CGAL::assign(seg, obj))
    {
        stream << "Segment: " << seg;
    }
    else if (CGAL::assign(ray, obj))
    {
        stream << "Ray: " << ray;
    }
    else if (CGAL::assign(line, obj))
    {
        stream << "Line: " << line;
    }
    return stream.str();
}

bool isDegenerate(const Object &obj)
{
    Segment3 seg;
    if (CGAL::assign(seg, obj))
    {
        return seg.source() == seg.target();
    }
    // Rays and lines cannot be degenerate in the same sense as segments.
    return false;
}


bool is_degenerate(Delaunay::Cell_handle cell) {
    Point p0 = cell->vertex(0)->point();
    Point p1 = cell->vertex(1)->point();
    Point p2 = cell->vertex(2)->point();
    Point p3 = cell->vertex(3)->point();

    // Compute the volume of the tetrahedron using the determinant formula
    K::FT vol = CGAL::volume(p0, p1, p2, p3);

    // If the volume is very close to zero, treat the cell as degenerate
    return CGAL::abs(vol) < 1e-6;  // Adjust the threshold if needed
}


bool is_bipolar(float val1, float val2, float isovalue)
{
    return ((val1 < isovalue) && (val2 >= isovalue)) || ((val1 >= isovalue) && (val2 < isovalue));
}


Point compute_centroid(const std::vector<Point> &points, bool supersample, int ratio)
{
    float sumX = 0, sumY = 0, sumZ = 0;
    for (const auto &pt : points)
    {
        sumX += pt.x();
        sumY += pt.y();
        sumZ += pt.z();
    }
    int n = points.size();
    float x, y, z;
    x = sumX / n;
    y = sumY / n;
    z = sumZ / n;

/*     if (supersample) {
        x = x / ratio;
        y = y / ratio;
        z = z / ratio;
    } */
    return Point(x, y, z);
}



std::array<Point, 8> get_cube_corners(const Point &center, float side_length)
{
    float half_side = side_length / 2.0;
    return {{
        Point(center.x() - half_side, center.y() - half_side, center.z() - half_side), // 0
        Point(center.x() + half_side, center.y() - half_side, center.z() - half_side), // 1
        Point(center.x() + half_side, center.y() + half_side, center.z() - half_side), // 2
        Point(center.x() - half_side, center.y() + half_side, center.z() - half_side), // 3
        Point(center.x() - half_side, center.y() - half_side, center.z() + half_side), // 4
        Point(center.x() + half_side, center.y() - half_side, center.z() + half_side), // 5
        Point(center.x() + half_side, center.y() + half_side, center.z() + half_side), // 6
        Point(center.x() - half_side, center.y() + half_side, center.z() + half_side)  // 7
    }};
}

// 1 for positive, -1 for negative
int get_orientation(const int iFacet, const Point v1, const Point v2, const float f1, const float f2)
{
    bool flag_v1_positive;
    if (f1 >= f2)
    {
        flag_v1_positive = true;
    }
    else
    {
        flag_v1_positive = false;
    }
    if ((iFacet % 2) == 0)
    {
        if (flag_v1_positive)
        {
            return 1;
        }
        else
        {
            return -1;
        }
    }
    else
    {
        if (flag_v1_positive)
        {
            return -1;
        }
        else
        {
            return 1;
        }
    }
}


// Function to check if a scalar value is positive or negative
bool isPositive(double value) {
    return value >= isovalue;
}

