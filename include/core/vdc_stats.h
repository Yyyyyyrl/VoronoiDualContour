//! @file vdc_stats.h
//! @brief Header file for collecting and reporting summary statistics.

#ifndef VDC_STATS_H
#define VDC_STATS_H

#include "core/vdc.h"
#include <array>
#include <cstddef>
#include <vector>

//! @brief Structure to hold summary statistics for the entire pipeline.
/*!
 * Collects various statistics about the Voronoi diagram, Delaunay triangulation,
 * and the extracted isosurface for reporting purposes.
 */
struct SummaryStats
{
    std::size_t active_cubes = 0;                   //!< Number of active cubes containing data
    std::size_t delaunay_vertices = 0;              //!< Number of vertices in the Delaunay triangulation
    std::size_t delaunay_cells = 0;                 //!< Number of cells in the Delaunay triangulation
    std::size_t voronoi_vertices = 0;               //!< Number of vertices in the Voronoi diagram
    std::size_t voronoi_edges = 0;                  //!< Number of edges in the Voronoi diagram
    std::size_t voronoi_facets = 0;                 //!< Number of facets in the Voronoi diagram
    std::size_t voronoi_cells = 0;                  //!< Number of cells in the Voronoi diagram
    std::size_t min_cell_vertices = 0;              //!< Minimum number of vertices in any Voronoi cell
    std::size_t max_cell_vertices = 0;              //!< Maximum number of vertices in any Voronoi cell
    double avg_cell_vertices = 0.0;                 //!< Average number of vertices per Voronoi cell
    std::size_t min_cell_facets = 0;                //!< Minimum number of facets in any Voronoi cell
    std::size_t max_cell_facets = 0;                //!< Maximum number of facets in any Voronoi cell
    double avg_cell_facets = 0.0;                   //!< Average number of facets per Voronoi cell
    int min_cell_index = -1;                        //!< Index of cell with minimum vertices
    int max_cell_index = -1;                        //!< Index of cell with maximum vertices
    std::size_t min_facet_vertices = 0;             //!< Minimum number of vertices in any Voronoi facet
    std::size_t max_facet_vertices = 0;             //!< Maximum number of vertices in any Voronoi facet
    double avg_facet_vertices = 0.0;                //!< Average number of vertices per Voronoi facet
    std::size_t min_facet_edges = 0;                //!< Minimum number of edges in any Voronoi facet
    std::size_t max_facet_edges = 0;                //!< Maximum number of edges in any Voronoi facet
    double avg_facet_edges = 0.0;                   //!< Average number of edges per Voronoi facet
    int min_facet_index = -1;                       //!< Index of facet with minimum vertices
    int max_facet_index = -1;                       //!< Index of facet with maximum vertices
    std::array<std::size_t, 4> facet_match_counts{0, 0, 0, 0}; //!< Counts of bipolar matching methods used
    std::size_t iso_vertices = 0;                   //!< Number of vertices in the extracted isosurface
    std::size_t iso_triangles = 0;                  //!< Number of triangles in the extracted isosurface
    bool multi_isov = false;                        //!< Flag indicating if multi-isovalue mode was used
    std::size_t collapsed_vertices_removed = 0;     //!< Number of vertices removed during collapse operation
    std::size_t collapsed_edges_removed = 0;        //!< Number of edges removed during collapse operation
    std::size_t mod_cyc_flips = 0;                  //!< Total number of edge flips during modify-cycles pass
    std::size_t mod_cyc_interior_flips = 0;         //!< Number of interior edge flips during modify-cycles
    std::size_t mod_cyc_boundary_flips = 0;         //!< Number of boundary edge flips during modify-cycles
    std::size_t isovertex_clipped_count = 0;        //!< Number of isosurface vertices that were clipped
    double isovertex_max_clip_distance = 0.0;       //!< Maximum clipping distance for any isosurface vertex
};

//! @brief Collects summary statistics from the pipeline data structures.
/*!
 * This function gathers various statistics from the active cubes, Delaunay triangulation,
 * Voronoi diagram, and extracted isosurface for reporting purposes.
 *
 * @param activeCubes Vector of active cubes in the grid
 * @param dt The Delaunay triangulation
 * @param vd The Voronoi diagram
 * @param iso_surface The extracted isosurface
 * @param multi_isov Flag indicating if multi-isovalue mode was used
 * @param collapsed_vertices_removed Number of vertices removed during collapse
 * @param collapsed_edges_removed Number of edges removed during collapse
 * @param mod_cyc_flips Total number of edge flips during modify-cycles
 * @param mod_cyc_interior_flips Number of interior edge flips
 * @param mod_cyc_boundary_flips Number of boundary edge flips
 * @param isovertex_clipped_count Number of clipped isosurface vertices
 * @param isovertex_max_clip_distance Maximum clipping distance
 * @return SummaryStats structure containing all collected statistics
 */
SummaryStats collect_summary_stats(const std::vector<Cube> &activeCubes,
                                   const Delaunay &dt,
                                   const VoronoiDiagram &vd,
                                   const IsoSurface &iso_surface,
                                   bool multi_isov,
                                   std::size_t collapsed_vertices_removed,
                                   std::size_t collapsed_edges_removed,
                                   std::size_t mod_cyc_flips,
                                   std::size_t mod_cyc_interior_flips,
                                   std::size_t mod_cyc_boundary_flips,
                                   std::size_t isovertex_clipped_count,
                                   double isovertex_max_clip_distance);

//! @brief Prints a formatted summary report of collected statistics.
/*!
 * @param stats The SummaryStats structure containing statistics to report
 */
void print_summary_report(const SummaryStats &stats);

#endif // VDC_STATS_H
