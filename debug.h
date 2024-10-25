#ifndef DEBUG_H
#define DEBUG_H

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>
#include <limits>
#include <algorithm>
#include <teem/nrrd.h>
#include <CGAL/bounding_box.h>
#include <CGAL/intersections.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_with_circumcenter_3.h>

// Define the kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

// Define the vertex base with information
typedef CGAL::Triangulation_vertex_base_with_info_3<bool, K> Vb;
typedef CGAL::Delaunay_triangulation_cell_base_3<K> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;

// Define Delaunay Triangulation
typedef CGAL::Delaunay_triangulation_3<K, Tds> Delaunay;
typedef Delaunay::Point Point;
typedef Delaunay::Cell_handle Cell_handle;
typedef Delaunay::Vertex_handle Vertex_handle;
typedef Delaunay::Cell_circulator Cell_circulator;
typedef CGAL::Object Object;
typedef Delaunay::Edge Edge;
typedef Delaunay::Facet Facet;
typedef Delaunay::Triangle Triangle;
typedef K::Segment_3 Segment3;
typedef K::Ray_3 Ray3;
typedef K::Line_3 Line3;
typedef K::Point_3 Point3;
typedef K::Vector_3 Vector3;


extern bool debug;
extern bool indicator;

void print_facet(Facet f);
void print_cell(Delaunay::Cell c);

#endif