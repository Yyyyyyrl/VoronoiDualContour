//! @file vdc_io.h
//! @brief Header file for input/output operations, including mesh and Voronoi diagram export.

#ifndef VDC_IO_H
#define VDC_IO_H

#include "core/vdc_utilities.h"
#include "voronoi/vdc_voronoi.h"

//! @brief Writes an isosurface mesh in OFF format (single-isovalue case).
/*!
 * @param filename The output file path.
 * @param iso_surface The isosurface container providing vertices and triangles.
 */
void writeOFFSingle(const std::string &filename, const IsoSurface &iso_surface);

//! @brief Writes an isosurface mesh in PLY format (single-isovalue case).
/*!
 * @param filename The output file path.
 * @param iso_surface The isosurface container providing vertices and triangles.
 */
void writePLYSingle(const std::string &filename, const IsoSurface &iso_surface);

//! @brief Writes an isosurface mesh in OFF format (multi-isovalue case).
/*!
 * @param filename The output file path.
 * @param voronoiDiagram The Voronoi diagram containing isosurface data.
 * @param iso_surface The isosurface container providing vertices and triangles.
 */
void writeOFFMulti(const std::string &filename, const VoronoiDiagram &voronoiDiagram,
                   const IsoSurface &iso_surface);

//! @brief Writes an isosurface mesh in PLY format (multi-isovalue case).
/*!
 * @param filename The output file path.
 * @param voronoiDiagram The Voronoi diagram containing isosurface data.
 * @param iso_surface The isosurface container providing vertices and triangles.
 */
void writePLYMulti(const std::string &filename, const VoronoiDiagram &voronoiDiagram,
                   const IsoSurface &iso_surface);

//! @brief Exports Voronoi diagram data to a CSV file.
/*!
 * @param voronoiDiagram The Voronoi diagram containing vertices and edges.
 * @param filename The output CSV file path.
 */
void export_voronoi_to_csv(const VoronoiDiagram &voronoiDiagram, const std::string &filename);

#endif // VDC_IO_H
