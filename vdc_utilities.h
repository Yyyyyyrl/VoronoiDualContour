#ifndef VDC_UTILITIES_H
#define VDC_UTILITIES_H

#include "vdc_debug.h"
#include "vdc_type.h"
#include "vdc_cube.h"
#include "vdc_grid.h"
#include "vdc_voronoi.h"
#include "vdc_delaunay.h"
#include "vdc_globalvar.h"

std::string objectToString(const CGAL::Object &obj);

struct ObjectComparator
{
    bool operator()(const Object &obj1, const Object &obj2) const
    {
        return objectToString(obj1) < objectToString(obj2);
    }
};


bool is_bipolar(float val1, float val2, float isovalue = 0);
bool isDegenerate(const Object &obj);
bool is_degenerate(Delaunay::Cell_handle cell);
bool isPositive(double value);

/*
General Helper Functions
*/
Point compute_centroid(const std::vector<Point> &points, bool supersample, int ratio);
std::array<Point, 8> get_cube_corners(const Point &center, float side_length);
int get_orientation(const int iFacet, const Point v1, const Point v2, const float f1, const float f2);



#endif