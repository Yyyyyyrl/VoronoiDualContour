//! @file vdc_io.h
//! @brief Header file for input/output operations, including mesh and Voronoi diagram export.

#ifndef VDC_IO_H
#define VDC_IO_H

#include "core/vdc_utilities.h"
#include "voronoi/vdc_voronoi.h"

//! @brief Writes an isosurface mesh in OFF format (single-isovalue case).
/*!
 * @param filename The output file path.
 * @param vertices The list of vertices in the mesh.
 * @param triangles The list of Delaunay triangles forming the mesh.
 */
void writeOFFSingle(const std::string &filename, const std::vector<Point> &vertices,
                    const std::vector<DelaunayTriangle> &triangles);

//! @brief Writes an isosurface mesh in PLY format (single-isovalue case).
/*!
 * @param filename The output file path.
 * @param vertices The list of vertices in the mesh.
 * @param triangles The list of Delaunay triangles forming the mesh.
 */
void writePLYSingle(const std::string &filename, const std::vector<Point> &vertices,
                    const std::vector<DelaunayTriangle> &triangles);

//! @brief Writes an isosurface mesh in OFF format (multi-isovalue case).
/*!
 * @param filename The output file path.
 * @param voronoiDiagram The Voronoi diagram containing isosurface data.
 * @param isoTriangles The list of isosurface triangles.
 */
void writeOFFMulti(const std::string &filename, const VoronoiDiagram &voronoiDiagram,
                   const std::vector<std::tuple<int, int, int>> &isoTriangles, IsoSurface &iso_surface);

//! @brief Writes an isosurface mesh in PLY format (multi-isovalue case).
/*!
 * @param filename The output file path.
 * @param voronoiDiagram The Voronoi diagram containing isosurface data.
 * @param isoTriangles The list of isosurface triangles.
 */
void writePLYMulti(const std::string &filename, const VoronoiDiagram &voronoiDiagram,
                   const std::vector<std::tuple<int, int, int>> &isoTriangles, IsoSurface &iso_surface);

//! @brief Exports Voronoi diagram data to a CSV file.
/*!
 * @param voronoiDiagram The Voronoi diagram containing vertices and edges.
 * @param filename The output CSV file path.
 */
void export_voronoi_to_csv(const VoronoiDiagram &voronoiDiagram, const std::string &filename);

#endif // VDC_IO_H
