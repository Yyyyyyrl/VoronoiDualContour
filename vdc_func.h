#ifndef VDC_FUNC_H
#define VDC_FUNC_H

#include "vdc_type.h"
#include "vdc_utilities.h"
#include "vdc_io.h"
#include "vdc_globalvar.h"

/*Struct*/

struct PointApproxEqual
{
    bool operator()(const Point &p1, const Point &p2) const
    {
        const double EPSILON = 1e-6;
        return (CGAL::squared_distance(p1, p2) < EPSILON * EPSILON);
    }
};



/*
Functions for both single/multi isov that computes the vertices and faces of the final mesh
*/
std::vector<DelaunayTriangle> computeDualTriangles(std::vector<CGAL::Object> &voronoi_edges, std::map<Point, float> &vertexValueMap, CGAL::Epick::Iso_cuboid_3 &bbox, std::map<Object, std::vector<Facet>, ObjectComparator> &delaunay_facet_to_voronoi_edge_map, Delaunay &dt, ScalarGrid &grid);
void computeDualTrianglesMulti(
    VoronoiDiagram &voronoiDiagram,
    CGAL::Epick::Iso_cuboid_3 &bbox,
    std::map<CGAL::Object, std::vector<Facet>, ObjectComparator> &delaunay_facet_to_voronoi_edge_map,
    ScalarGrid &grid,
    float isovalue);

void Compute_Isosurface_Vertices_Multi(VoronoiDiagram &voronoiDiagram, float isovalue);
void Compute_Isosurface_Vertices_Single(VoronoiDiagram &voronoiDiagram, ScalarGrid &grid);

/*
Setting up the delaunay triangulation
*/
void construct_delaunay_triangulation(Grid &grid, const std::vector<GRID_FACETS> &gridfacets);
std::vector<Point> add_dummy_from_facet(const GRID_FACETS &facet, const Grid &grid);
/*
Functions that relates to building the voronoi diagram
*/
void construct_voronoi_vertices(VoronoiDiagram &voronoiDiagram);
void construct_voronoi_cells(VoronoiDiagram &voronoiDiagram);

void compute_voronoi_values(VoronoiDiagram &voronoiDiagram, ScalarGrid &grid);
void construct_voronoi_edges(
    VoronoiDiagram &voronoiDiagram,
    std::map<CGAL::Object, std::vector<Facet>, ObjectComparator> &delaunay_facet_to_voronoi_edge_map);

int handle_output_mesh(bool &retFlag, VoronoiDiagram &vd);

#endif