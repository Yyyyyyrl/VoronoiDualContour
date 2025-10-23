#include "core/vdc_stats.h"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>

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
                                   double isovertex_max_clip_distance)
{
    SummaryStats stats;
    stats.multi_isov = multi_isov;
    stats.active_cubes = activeCubes.size();
    stats.delaunay_vertices = static_cast<std::size_t>(std::distance(dt.finite_vertices_begin(), dt.finite_vertices_end()));
    stats.delaunay_cells = static_cast<std::size_t>(std::distance(dt.finite_cells_begin(), dt.finite_cells_end()));
    stats.voronoi_vertices = vd.vertices.size();
    stats.voronoi_edges = vd.edges.size();
    stats.voronoi_facets = vd.global_facets.size();
    stats.voronoi_cells = vd.cells.size();
    stats.collapsed_vertices_removed = collapsed_vertices_removed;
    stats.collapsed_edges_removed = collapsed_edges_removed;
    stats.mod_cyc_flips = mod_cyc_flips;
    stats.mod_cyc_interior_flips = mod_cyc_interior_flips;
    stats.mod_cyc_boundary_flips = mod_cyc_boundary_flips;
    stats.isovertex_clipped_count = isovertex_clipped_count;
    stats.isovertex_max_clip_distance = isovertex_max_clip_distance;

    std::size_t total_cell_vertices = 0;
    std::size_t total_cell_facets = 0;
    std::size_t minCellVertices = std::numeric_limits<std::size_t>::max();
    std::size_t maxCellVertices = 0;
    std::size_t minCellFacets = std::numeric_limits<std::size_t>::max();
    std::size_t maxCellFacets = 0;
    int minCellIndex = -1;
    int maxCellIndex = -1;
    for (const auto &cell : vd.cells)
    {
        const std::size_t vertCount = cell.vertices_indices.size();
        const std::size_t facetCount = cell.facet_indices.size();
        total_cell_vertices += vertCount;
        total_cell_facets += facetCount;

        if (vertCount < minCellVertices)
        {
            minCellVertices = vertCount;
            minCellFacets = facetCount;
            minCellIndex = cell.cellIndex;
        }
        if (vertCount > maxCellVertices)
        {
            maxCellVertices = vertCount;
            maxCellFacets = facetCount;
            maxCellIndex = cell.cellIndex;
        }
    }
    if (!vd.cells.empty())
    {
        const double denom = static_cast<double>(vd.cells.size());
        stats.avg_cell_vertices = static_cast<double>(total_cell_vertices) / denom;
        stats.avg_cell_facets = static_cast<double>(total_cell_facets) / denom;
        stats.min_cell_vertices = minCellVertices;
        stats.max_cell_vertices = maxCellVertices;
        stats.min_cell_facets = minCellFacets;
        stats.max_cell_facets = maxCellFacets;
        stats.min_cell_index = minCellIndex;
        stats.max_cell_index = maxCellIndex;
    }

    std::size_t total_facet_vertices = 0;
    std::size_t total_facet_edges = 0;
    std::size_t minFacetVertices = std::numeric_limits<std::size_t>::max();
    std::size_t maxFacetVertices = 0;
    std::size_t minFacetEdges = std::numeric_limits<std::size_t>::max();
    std::size_t maxFacetEdges = 0;
    int minFacetIndex = -1;
    int maxFacetIndex = -1;
    for (const auto &gf : vd.global_facets)
    {
        const std::size_t vertCount = gf.vertices_indices.size();
        const std::size_t edgeCount = gf.voronoi_edge_indices.size();
        total_facet_vertices += vertCount;
        total_facet_edges += edgeCount;

        if (vertCount < minFacetVertices)
        {
            minFacetVertices = vertCount;
            minFacetEdges = edgeCount;
            minFacetIndex = gf.index;
        }
        if (vertCount > maxFacetVertices)
        {
            maxFacetVertices = vertCount;
            maxFacetEdges = edgeCount;
            maxFacetIndex = gf.index;
        }

        const auto methodIndex = static_cast<std::size_t>(gf.bipolar_match_method);
        if (methodIndex < stats.facet_match_counts.size())
            ++stats.facet_match_counts[methodIndex];
    }
    if (!vd.global_facets.empty())
    {
        const double denom = static_cast<double>(vd.global_facets.size());
        stats.avg_facet_vertices = static_cast<double>(total_facet_vertices) / denom;
        stats.avg_facet_edges = static_cast<double>(total_facet_edges) / denom;
        stats.min_facet_vertices = minFacetVertices;
        stats.max_facet_vertices = maxFacetVertices;
        stats.min_facet_edges = minFacetEdges;
        stats.max_facet_edges = maxFacetEdges;
        stats.min_facet_index = minFacetIndex;
        stats.max_facet_index = maxFacetIndex;
    }

    stats.iso_vertices = iso_surface.isosurfaceVertices.size();
    stats.iso_triangles = multi_isov ? iso_surface.isosurfaceTrianglesMulti.size()
                                     : iso_surface.isosurfaceTrianglesSingle.size();
    return stats;
}

void print_summary_report(const SummaryStats &stats)
{
    std::cout << "\n[SUMMARY] Run statistics\n";
    std::cout << "  Active cubes: " << stats.active_cubes << "\n";
    std::cout << "  Delaunay finite vertices: " << stats.delaunay_vertices
              << ", finite cells: " << stats.delaunay_cells << "\n";
    std::cout << "  Voronoi vertices: " << stats.voronoi_vertices
              << ", edges: " << stats.voronoi_edges
              << ", facets: " << stats.voronoi_facets
              << ", cells: " << stats.voronoi_cells << "\n";
    if (stats.collapsed_vertices_removed > 0 || stats.collapsed_edges_removed > 0)
    {
        std::cout << "    (collapseSmallEdges removed vertices=" << stats.collapsed_vertices_removed
                  << ", edges=" << stats.collapsed_edges_removed << ")\n";
    }
    if (stats.collapsed_vertices_removed > 0 || stats.collapsed_edges_removed > 0)
    {
        std::cout << "    (collapseSmallEdges removed vertices=" << stats.collapsed_vertices_removed
                  << ", edges=" << stats.collapsed_edges_removed << ")\n";
    }

    if (stats.max_cell_index != -1)
    {
        std::cout << "  Voronoi cell vertices (min / avg / max): "
                  << stats.min_cell_vertices << " / " << std::fixed << std::setprecision(2) << stats.avg_cell_vertices
                  << " / " << stats.max_cell_vertices << std::defaultfloat
                  << "  [indices: min=" << stats.min_cell_index
                  << ", max=" << stats.max_cell_index << "]\n";
        std::cout << "  Voronoi cell facets (min / avg / max): "
                  << stats.min_cell_facets << " / " << std::fixed << std::setprecision(2) << stats.avg_cell_facets
                  << " / " << stats.max_cell_facets << std::defaultfloat << "\n";
    }

    if (stats.max_facet_index != -1)
    {
        std::cout << "  Voronoi facet vertices (min / avg / max): "
                  << stats.min_facet_vertices << " / " << std::fixed << std::setprecision(2) << stats.avg_facet_vertices
                  << " / " << stats.max_facet_vertices << std::defaultfloat
                  << "  [indices: min=" << stats.min_facet_index
                  << ", max=" << stats.max_facet_index << "]\n";
        std::cout << "  Voronoi facet edges (min / avg / max): "
                  << stats.min_facet_edges << " / " << std::fixed << std::setprecision(2) << stats.avg_facet_edges
                  << " / " << stats.max_facet_edges << std::defaultfloat << "\n";
    }

    static const char *method_names[] = {
        "SEP_POS", "SEP_NEG", "UNCONSTRAINED_MATCH", "UNDEFINED_MATCH_TYPE"};
    bool header_printed = false;
    for (std::size_t i = 0; i < stats.facet_match_counts.size(); ++i)
    {
        if (stats.facet_match_counts[i] == 0)
            continue;
        if (!header_printed)
        {
            std::cout << "  Facet match methods:";
            header_printed = true;
        }
        std::cout << ' ' << method_names[i] << "=" << stats.facet_match_counts[i];
    }
    if (header_printed)
        std::cout << "\n";

    std::cout << "  Iso-surface vertices: " << stats.iso_vertices
              << ", triangles: " << stats.iso_triangles
              << (stats.multi_isov ? " (multi-isov)" : " (single-isov)") << "\n";
    if (stats.multi_isov)
    {
        std::cout << "  Modify-cycles flips: " << stats.mod_cyc_flips
                  << " (interior: " << stats.mod_cyc_interior_flips
                  << ", boundary: " << stats.mod_cyc_boundary_flips << ")\n";
        if (stats.isovertex_clipped_count > 0)
        {
            std::cout << "  Iso-vertex clipping: " << stats.isovertex_clipped_count << " vertices clipped"
                      << ", max distance: " << std::fixed << std::setprecision(6)
                      << stats.isovertex_max_clip_distance << std::defaultfloat << "\n";
        }
    }

    std::cout << std::defaultfloat;
}
