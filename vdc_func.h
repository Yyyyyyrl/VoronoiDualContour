#ifndef VDC_FUNC_H
#define VDC_FUNC_H

#include "vdc_type.h"
#include "vdc_utilities.h"
#include "vdc_io.h"


/*
Functions for both single/multi isov that computes the vertices and faces of the final mesh
*/
std::vector<DelaunayTriangle> computeDualTriangles(std::vector<CGAL::Object> &voronoi_edges, std::map<Point, float> &vertexValueMap, CGAL::Epick::Iso_cuboid_3 &bbox, std::map<Object, std::vector<Facet>, ObjectComparator> &delaunay_facet_to_voronoi_edge_map, Delaunay &dt, ScalarGrid &grid);
void computeDualTrianglesMulti(std::vector<CGAL::Object> &voronoi_edges, CGAL::Epick::Iso_cuboid_3 &bbox, std::map<Object, std::vector<Facet>, ObjectComparator> &delaunay_facet_to_voronoi_edge_map, ScalarGrid &grid);

void Compute_Isosurface_Vertices_Multi(std::vector<VoronoiCell> &voronoi_cells);
void Compute_Isosurface_Vertices_Single(ScalarGrid &grid);

/*
Setting up the delaunay triangulation
*/
void construct_delaunay_triangulation();


/*
Functions that relates to building the voronoi diagram
*/
void construct_voronoi_vertices(std::set<Point> &seen_points, std::vector<Point> &voronoi_vertices);
void construct_voronoi_cells(std::vector<VoronoiCell> &voronoi_cells);
void compute_voronoi_values(std::vector<Point> &voronoi_vertices, ScalarGrid &grid, std::vector<float> &voronoi_vertex_values);
void construct_voronoi_edges(std::map<Object, std::vector<Facet>, ObjectComparator> &delaunay_facet_to_voronoi_edge_map, std::set<std::string> &seen_edges, std::vector<CGAL::Object> &voronoi_edges);

int handle_output_mesh(bool &retFlag);

#endif