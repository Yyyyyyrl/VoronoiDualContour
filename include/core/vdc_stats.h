#ifndef VDC_STATS_H
#define VDC_STATS_H

#include "core/vdc.h"
#include <array>
#include <cstddef>
#include <vector>

struct SummaryStats
{
    std::size_t active_cubes = 0;
    std::size_t delaunay_vertices = 0;
    std::size_t delaunay_cells = 0;
    std::size_t voronoi_vertices = 0;
    std::size_t voronoi_edges = 0;
    std::size_t voronoi_facets = 0;
    std::size_t voronoi_cells = 0;
    std::size_t min_cell_vertices = 0;
    std::size_t max_cell_vertices = 0;
    double avg_cell_vertices = 0.0;
    std::size_t min_cell_facets = 0;
    std::size_t max_cell_facets = 0;
    double avg_cell_facets = 0.0;
    int min_cell_index = -1;
    int max_cell_index = -1;
    std::size_t min_facet_vertices = 0;
    std::size_t max_facet_vertices = 0;
    double avg_facet_vertices = 0.0;
    std::size_t min_facet_edges = 0;
    std::size_t max_facet_edges = 0;
    double avg_facet_edges = 0.0;
    int min_facet_index = -1;
    int max_facet_index = -1;
    std::array<std::size_t, 4> facet_match_counts{0, 0, 0, 0};
    std::size_t iso_vertices = 0;
    std::size_t iso_triangles = 0;
    bool multi_isov = false;
    std::size_t collapsed_vertices_removed = 0;
    std::size_t collapsed_edges_removed = 0;
    std::size_t mod_cyc_flips = 0;
    std::size_t mod_cyc_interior_flips = 0;
    std::size_t mod_cyc_boundary_flips = 0;
    std::size_t isovertex_clipped_count = 0;
    double isovertex_max_clip_distance = 0.0;
};

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

void print_summary_report(const SummaryStats &stats);

#endif // VDC_STATS_H
