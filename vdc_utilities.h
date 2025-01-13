//! @file vdc_utilities.h
//! @brief Utility functions for Voronoi and Delaunay computations, including centroid calculations and geometric checks.

#ifndef VDC_UTILITIES_H
#define VDC_UTILITIES_H

#include "vdc_debug.h"
#include "vdc_type.h"
#include "vdc_cube.h"
#include "vdc_grid.h"
#include "vdc_voronoi.h"
#include "vdc_delaunay.h"
#include "vdc_globalvar.h"

//! @brief Converts a CGAL::Object into a string representation.
/*!
 * This function helps in debugging by converting a geometric object (e.g., Segment, Ray, Line)
 * into a human-readable string format.
 * 
 * @param obj The CGAL::Object to convert.
 * @return A string representation of the object.
 */
std::string objectToString(const CGAL::Object &obj);

//! @brief Comparator for CGAL::Object instances.
/*!
 * Allows comparison of CGAL objects based on their string representation.
 */
struct ObjectComparator
{
    //! @brief Compare two CGAL::Object instances.
    /*!
     * @param obj1 The first object.
     * @param obj2 The second object.
     * @return `true` if obj1 < obj2, otherwise `false`.
     */
    bool operator()(const Object &obj1, const Object &obj2) const
    {
        return objectToString(obj1) < objectToString(obj2);
    }
};

//! @brief Checks if two scalar values are bipolar with respect to an isovalue.
/*!
 * Determines if one value is above and the other is below the specified isovalue.
 * 
 * @param val1 First scalar value.
 * @param val2 Second scalar value.
 * @param isovalue The isovalue used for comparison (default is 0).
 * @return `true` if the values are bipolar, otherwise `false`.
 */
bool is_bipolar(float val1, float val2, float isovalue = 0);

//! @brief Checks if a CGAL::Object is degenerate.
/*!
 * A segment is degenerate if its start and end points are the same.
 * 
 * @param obj The CGAL::Object to check.
 * @return `true` if the object is degenerate, otherwise `false`.
 */
bool isDegenerate(const Object &obj);

//! @brief Checks if a Delaunay cell is degenerate.
/*!
 * A Delaunay cell is degenerate if its volume is approximately zero.
 * 
 * @param cell The Delaunay cell handle.
 * @return `true` if the cell is degenerate, otherwise `false`.
 */
bool is_degenerate(Delaunay::Cell_handle cell);

//! @brief Checks if a value is positive or zero.
/*!
 * @param value The scalar value to check.
 * @return `true` if the value is positive or zero, otherwise `false`.
 */
bool isPositive(double value);

/* General Helper Functions */

//! @brief Computes the centroid of a set of points.
/*!
 * This function calculates the geometric centroid of a set of points.
 * 
 * @param points Vector of points to compute the centroid for.
 * @param supersample Indicates if supersampling is applied (optional).
 * @param ratio Supersampling ratio (optional).
 * @return The computed centroid as a Point.
 */
Point compute_centroid(const std::vector<Point> &points, bool supersample, int ratio);

//! @brief Computes the centroid of a set of points using CGAL.
/*!
 * @param points Vector of points to compute the centroid for.
 * @return The computed centroid as a Point.
 */
Point compute_centroid(const std::vector<Point> &points);

//! @brief Gets the corners of a cube.
/*!
 * Computes the eight corner points of a cube given its center and side length.
 * 
 * @param center The center of the cube.
 * @param side_length The side length of the cube.
 * @return An array containing the eight corner points.
 */
std::array<Point, 8> get_cube_corners(const Point &center, float side_length);

//! @brief Determines the orientation of a facet.
/*!
 * Based on the facet index, vertex positions, and scalar values, this function determines
 * the orientation of a facet (positive or negative).
 * 
 * @param iFacet The facet index.
 * @param v1 First vertex of the facet.
 * @param v2 Second vertex of the facet.
 * @param f1 Scalar value at v1.
 * @param f2 Scalar value at v2.
 * @return 1 for positive orientation, -1 for negative.
 */
int get_orientation(const int iFacet, const Point v1, const Point v2, const float f1, const float f2);

#endif // VDC_UTILITIES_H
