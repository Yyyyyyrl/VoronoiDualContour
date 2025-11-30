#include "processing/vdc_func.h"
#include "core/vdc_debug.h"
#include "core/vdc_timing.h"
#include <sstream>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

struct PairHash
{
    size_t operator()(const std::pair<int, int> &p) const noexcept
    {
        return (static_cast<size_t>(p.first) << 32) ^ static_cast<size_t>(p.second);
    }
};

using EdgeFacetMap = std::unordered_map<std::pair<int, int>, std::vector<int>, PairHash>;

//! @brief Assigns global indices to all Delaunay edges and stores them in CellInfo.
/*!
 * Each Delaunay edge is assigned a unique global index. This index is stored in
 * the edge_index[i][j] field of CellInfo for all incident cells, enabling O(1)
 * lookup of edge indices without using hash maps.
 *
 * @param dt The Delaunay triangulation.
 * @return The total number of Delaunay edges (num_del_edges).
 */
static int assign_delaunay_edge_indices(Delaunay &dt)
{
    // Finite edges iterator guarantees each edge is visited exactly once
    // (with finite endpoints), so no per-vertex filtering is needed.
    int num_del_edges = 0;

    for (auto eit = dt.finite_edges_begin(); eit != dt.finite_edges_end(); ++eit)
    {
        const Edge &ed = *eit;
        Cell_handle cell = ed.first;
        int i = ed.second;
        int j = ed.third;
        Vertex_handle v1 = cell->vertex(i);
        Vertex_handle v2 = cell->vertex(j);

        const int edge_id = num_del_edges++;

        // Store the edge id into every incident cell (skip infinite ones).
        Delaunay::Cell_circulator cc = dt.incident_cells(ed);
        Delaunay::Cell_circulator start = cc;
        do
        {
            if (!dt.is_infinite(cc))
            {
                const int ci1 = cc->index(v1);
                const int ci2 = cc->index(v2);
                cc->info().edge_index[ci1][ci2] = edge_id;
                cc->info().edge_index[ci2][ci1] = edge_id;
            }
            ++cc;
        } while (cc != start);
    }

    return num_del_edges;
}


//! @brief Constructs Voronoi vertices for the given voronoi Diagram instance.
void construct_voronoi_vertices(VoronoiDiagram &voronoiDiagram, Delaunay &dt)
{
    int vertexIndex = 0;
    voronoiDiagram.vertices.reserve(dt.number_of_finite_cells());
    for (Delaunay::Finite_cells_iterator cit = dt.finite_cells_begin(); cit != dt.finite_cells_end(); ++cit)
    {
        cit->info().dualVoronoiVertexIndex = -1; // Default to invalid until a vertex is recorded

        Point circumcenter = dt.dual(cit);

        // Check for degenerate (flat/coplanar) Delaunay tetrahedron
        if (std::isnan(circumcenter.x()) || std::isnan(circumcenter.y()) || std::isnan(circumcenter.z()))
        {
            // Print warning with Delaunay vertex details
            std::cerr << "[WARNING] Degenerate Delaunay tetrahedron detected (NaN circumcenter)\n";
            std::cerr << "  Delaunay vertices:\n";
            std::cerr << "    v0: " << cit->vertex(0)->point()
                      << " (index: " << cit->vertex(0)->info().index << ")\n";
            std::cerr << "    v1: " << cit->vertex(1)->point()
                      << " (index: " << cit->vertex(1)->info().index << ")\n";
            std::cerr << "    v2: " << cit->vertex(2)->point()
                      << " (index: " << cit->vertex(2)->info().index << ")\n";
            std::cerr << "    v3: " << cit->vertex(3)->point()
                      << " (index: " << cit->vertex(3)->info().index << ")\n";

            // Use centroid of Delaunay vertices as Voronoi vertex instead of skipping
            Point v0 = cit->vertex(0)->point();
            Point v1 = cit->vertex(1)->point();
            Point v2 = cit->vertex(2)->point();
            Point v3 = cit->vertex(3)->point();

            double cx = (v0.x() + v1.x() + v2.x() + v3.x()) / 4.0;
            double cy = (v0.y() + v1.y() + v2.y() + v3.y()) / 4.0;
            double cz = (v0.z() + v1.z() + v2.z() + v3.z()) / 4.0;

            circumcenter = Point(cx, cy, cz);
            std::cerr << "  Using centroid as Voronoi vertex: " << circumcenter << "\n";
        }
        VoronoiVertex vv(circumcenter);
        vv.index = vertexIndex;
        voronoiDiagram.vertices.push_back(vv);
        cit->info().dualVoronoiVertexIndex = vertexIndex;
        vertexIndex++;
    }
}

//! @brief Computes Voronoi Vertex values using scalar grid interpolation
void compute_voronoi_values(VoronoiDiagram &voronoiDiagram, UnifiedGrid &grid)
{
    for (size_t i = 0; i < voronoiDiagram.vertices.size(); ++i)
    {
        Point vertex = voronoiDiagram.vertices[i].coord;
        voronoiDiagram.vertices[i].value = trilinear_interpolate(vertex, grid);
    }
}

//! @brief Constructs Voronoi cells from the Delaunay triangulation.
void construct_voronoi_cells_as_convex_hull(VoronoiDiagram &voronoiDiagram, Delaunay &dt)
{
    int index = 0;
    for (auto delaunay_vertex = dt.finite_vertices_begin(); delaunay_vertex != dt.finite_vertices_end(); ++delaunay_vertex)
    {
        if (delaunay_vertex->info().is_dummy)
        {
            continue;
        }
        VoronoiCell vc(delaunay_vertex);
        vc.cellIndex = index;

        std::vector<Cell_handle> incident_cells;
        dt.finite_incident_cells(delaunay_vertex, std::back_inserter(incident_cells));

        // Collect vertex indices combinatorially
        std::set<int> unique_vertex_indices_set;
        for (Cell_handle ch : incident_cells)
        {
            if (dt.is_infinite(ch))
            {
                continue; // Skip infinite cells
            }
            // Use direct index instead of dual point + search
            int vertex_index = ch->info().dualVoronoiVertexIndex;
            unique_vertex_indices_set.insert(vertex_index);
        }

        // Copy unique indices to vector
        vc.verticesIndices.assign(unique_vertex_indices_set.begin(), unique_vertex_indices_set.end());

        // Build vertex_points and vector for lookup (allow duplicates by using first idx for matching points)
        std::vector<Point> vertex_points;
        std::vector<std::pair<Point, int>> point_index_pairs;
        for (int idx : vc.verticesIndices)
        {
            Point p = voronoiDiagram.vertices[idx].coord;
            vertex_points.push_back(p);
            // Check if point already added (approx equal), if not, add pair
            bool found = false;
            for (const auto &pair : point_index_pairs)
            {
                if (PointApproxEqual()(pair.first, p))
                {
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                point_index_pairs.emplace_back(p, idx); // Use first idx for this point
            }
        }

        // No remove duplicates: pass all to hull, it will handle

        CGAL::convex_hull_3(vertex_points.begin(), vertex_points.end(), vc.polyhedron);

        // Extract facets from polyhedron
        for (auto facet_it = vc.polyhedron.facets_begin();
             facet_it != vc.polyhedron.facets_end(); ++facet_it)
        {
            VoronoiCellFacet vf;
            auto h = facet_it->facet_begin();
            do
            {
                Point p = h->vertex()->point();
                // Linear lookup in point_index_pairs
                bool found = false;
                for (const auto &pair : point_index_pairs)
                {
                    if (PointApproxEqual()(pair.first, p))
                    {
                        vf.verticesIndices.push_back(pair.second);
                        found = true;
                        break;
                    }
                }
                if (!found)
                {
                    std::cerr << "[WARNING] Point not found during facet extraction: " << p << "\n";
                }
                ++h;
            } while (h != facet_it->facet_begin());

            // Skip degenerate facets
            if (vf.verticesIndices.size() < 3)
                continue;

            int facet_index = voronoiDiagram.cell_facets.size();
            voronoiDiagram.cell_facets.push_back(vf);
            vc.facetIndices.push_back(facet_index);
        }

        voronoiDiagram.cells.push_back(vc);
        delaunay_vertex->info().voronoiCellIndex = index;
        index++;
    }
}

//! @brief Creates a Voronoi cell for a Delaunay vertex.
/*!
 * Initializes a Voronoi cell with the given cell index and Delaunay vertex handle.
 *
 * @param delaunay_vertex The Delaunay vertex to create the cell for.
 * @param cellIndex The index to assign to the cell.
 * @return The initialized Voronoi cell.
 */
static VoronoiCell create_voronoi_cell(Vertex_handle delaunay_vertex, int cellIndex)
{
    VoronoiCell vc(delaunay_vertex);
    vc.cellIndex = cellIndex;
    return vc;
}

//! @brief Collects unique vertex indices from incident cells.
/*!
 * Retrieves the Voronoi vertex indices from cells incident to a Delaunay vertex,
 * applying the old-to-new vertex index mapping.
 *
 * @param dt The Delaunay triangulation.
 * @param delaunay_vertex The Delaunay vertex to process.
 * @param voronoiDiagram The Voronoi diagram containing vertex mappings.
 * @param vertices_indices Vector to store the collected vertex indices.
 */
static void collect_cell_vertices(
    Delaunay &dt,
    Vertex_handle delaunay_vertex,
    VoronoiDiagram &voronoiDiagram,
    std::vector<int> &vertices_indices)
{
    std::vector<Cell_handle> incidentCells;
    dt.finite_incident_cells(delaunay_vertex, std::back_inserter(incidentCells));

    vertices_indices.clear();
    vertices_indices.reserve(incidentCells.size());

    for (Cell_handle c : incidentCells)
    {
        int vertex_index = c->info().dualVoronoiVertexIndex;
        // Only include cells where the dual Voronoi vertex is defined (>= 0)
        if (vertex_index >= 0)
        {
            vertices_indices.push_back(vertex_index);
        }
    }
    std::sort(vertices_indices.begin(), vertices_indices.end());
    vertices_indices.erase(std::unique(vertices_indices.begin(), vertices_indices.end()), vertices_indices.end());
}

//! @brief Returns the index of Voronoi cell edge dual to facet and in Voronoi cell around vertex that is not vA or vB.
/*!
 * Given a Delaunay facet and two vertices vA and vB of a Delaunay edge, this function
 * finds the third vertex of the facet (the one that is neither vA nor vB) and returns
 * the corresponding dual cell edge index stored in the facet info.
 *
 * @param delFacet The Delaunay facet to query
 * @param vA First vertex of the Delaunay edge
 * @param vB Second vertex of the Delaunay edge
 * @return Index of the Voronoi cell edge, or -1 if not found
 * @precondition Facet delFacet contains vertices vA and vB
 */
static inline int get_dual_cell_edge_index(
    const Facet &delFacet,
    const Vertex_handle vA,
    const Vertex_handle vB)
{
    const int NUM_DELAUNAY_CELL_VERTICES = 4;

    const Cell_handle cc = delFacet.first;
    int facet_index = delFacet.second;

    // Iterate through the 3 vertices of the facet (k = 0, 1, 2)
    for (int k = 0; k < 3; ++k)
    {
        const int kv = (facet_index + k + 1) % NUM_DELAUNAY_CELL_VERTICES;
        const Vertex_handle delVert = cc->vertex(kv);

        // Find the vertex that is not vA or vB and not the opposite vertex
        if ((delVert != vA) && (delVert != vB) && (kv != facet_index))
        {
            return cc->info().facet_info[facet_index].dualCellEdgeIndex[k];
        }
    }

    // Should not reach here if precondition is satisfied
    if (debug)
    {
        std::cerr << "[ERROR] get_dual_cell_edge_index: Unable to find dualCellEdgeIndex for facet " << facet_index << "\n";
    }
    return -1;
}

//! @brief Builds a facet from an incident edge using cell circulators.
/*!
 * Constructs a Voronoi facet by collecting vertices around an incident edge,
 * ordering them cyclically, and assigning scalar values. The facet is built
 * into the provided reference and, if valid, appended to the diagram.
 * Mirror facets are linked inline using the vor_facet_dual_to_edge vector.
 *
 * @param dt The Delaunay triangulation.
 * @param ed The incident edge to process.
 * @param delaunay_vertex The Delaunay vertex associated with the cell.
 * @param voronoiDiagram The Voronoi diagram containing vertex and value data.
 * @param facet_indices Vector to store the facet index for the owning cell.
 * @param vor_facet_dual_to_edge Vector mapping Delaunay edge index to first Voronoi facet (-1 if none yet).
 * @param vcIdx The Voronoi cell index (for diagnostics).
 * @param outFacet Output facet to populate on success.
 * @return true if a valid facet was constructed and appended; false otherwise.
 */
static bool build_facet_from_edge(
    Delaunay &dt,
    const Edge &ed,
    Vertex_handle delaunay_vertex,
    VoronoiDiagram &voronoiDiagram,
    std::vector<int> &facet_indices,
    std::vector<int> &vor_facet_dual_to_edge,
    int vcIdx,
    VoronoiCellFacet &outFacet)
{
    Cell_handle cell_ed = ed.first;
    int ei = ed.second;
    int ej = ed.third;
    Vertex_handle v1 = cell_ed->vertex(ei);
    Vertex_handle v2 = cell_ed->vertex(ej);

    (void)delaunay_vertex;
    (void)vcIdx;

    Delaunay::Cell_circulator cc = dt.incident_cells(ed);
    Delaunay::Cell_circulator start = cc;
    std::vector<int> facetVertices;
    int finite_cell_count = 0;
    do
    {
        if (!dt.is_infinite(cc))
        {
            int vertex_index = cc->info().dualVoronoiVertexIndex;
            // Only include vertices where the Voronoi vertex is defined (>= 0)
            if (vertex_index >= 0)
            {
                facetVertices.push_back(vertex_index);
                ++finite_cell_count;
            }
        }
        ++cc;
    } while (cc != start);

    // Check for degenerate facet: need at least 3 Voronoi vertices
    if (facetVertices.size() < 3)
    {
        if (debug)
        {
            std::cout << "[DEBUG] Degenerate facet for edge with " << finite_cell_count << " finite cells\n";
        }
        outFacet.verticesIndices.clear();
        outFacet.cellEdgeIndices.clear();
        return false;
    }

    // Check for degenerate facet: need at least 3 unique Voronoi vertices
    std::unordered_set<int> unique_vertices;
    unique_vertices.reserve(facetVertices.size());
    for (int idx : facetVertices)
        unique_vertices.insert(idx);
    if (unique_vertices.size() < 3)
    {
        if (debug)
        {
            std::cout << "[DEBUG] Degenerate facet for edge with " << finite_cell_count << " finite cells (duplicate vertices)\n";
        }
        outFacet.verticesIndices.clear();
        outFacet.cellEdgeIndices.clear();
        return false;
    }

    outFacet.verticesIndices = std::move(facetVertices);

    // Single-pass iteration using stored facet_info indices
    outFacet.cellEdgeIndices.clear();
    const int n = (int)outFacet.verticesIndices.size();
    outFacet.cellEdgeIndices.reserve(n);

    // Build parallel arrays of Voronoi vertices and their corresponding Delaunay facets
    // in a single pass through the facet circulator
    std::vector<int> vorVertexIndices;
    std::vector<Facet> delaunayFacets;
    vorVertexIndices.reserve(n);
    delaunayFacets.reserve(n);

    Delaunay::Facet_circulator delFacet_circ = dt.incident_facets(ed);
    Delaunay::Facet_circulator delFacet_start = delFacet_circ;
    do
    {
        const Cell_handle cc_facet = delFacet_circ->first;
        if (!dt.is_infinite(cc_facet))
        {
            const int vor_vertex_index = cc_facet->info().dualVoronoiVertexIndex;
            // Only include facets where the dual Voronoi vertex is defined (>= 0)
            if (vor_vertex_index >= 0)
            {
                vorVertexIndices.push_back(vor_vertex_index);
                delaunayFacets.push_back(*delFacet_circ);
            }
        }
        ++delFacet_circ;
    } while (delFacet_circ != delFacet_start);

    // For each edge of the Voronoi facet, find the corresponding cell edge index
    // using the stored dualCellEdgeIndex in facet_info (accessed via get_dual_cell_edge_index)
    for (int k = 0; k < n; ++k)
    {
        const int a = outFacet.verticesIndices[k];

        // Find the Delaunay facet corresponding to Voronoi vertex 'a'
        bool found = false;
        for (size_t m = 0; m < vorVertexIndices.size(); ++m)
        {
            if (vorVertexIndices[m] == a)
            {
                const Facet &delFacet = delaunayFacets[m];
                // Get the cell edge index directly from stored facet_info
                const int dual_cell_edge_index = get_dual_cell_edge_index(delFacet, v1, v2);
                outFacet.cellEdgeIndices.push_back(dual_cell_edge_index);
                found = true;
                break;
            }
        }

        if (!found)
        {
            // should rarely happen - only for boundary/degenerate cases
            if (debug)
            {
                std::cerr << "[WARNING] build_facet_from_edge: Voronoi vertex " << a
                          << " not found in incident facets for edge ("
                          << v1->info().index << ", " << v2->info().index << ")\n";
            }
            outFacet.cellEdgeIndices.push_back(-1);
        }
    }

    // Append to diagram only after fully building facet (avoid partial copies)
    int facetIndex = (int)voronoiDiagram.cell_facets.size();
    voronoiDiagram.cell_facets.push_back(outFacet);
    facet_indices.push_back(facetIndex);

    // Get the global Delaunay edge index from CellInfo and link mirror facets inline
    int ci = cell_ed->index(v1);
    int cj = cell_ed->index(v2);
    int edge_id = cell_ed->info().edge_index[ci][cj];

    if (edge_id >= 0 && edge_id < static_cast<int>(vor_facet_dual_to_edge.size()))
    {
        if (vor_facet_dual_to_edge[edge_id] == -1)
        {
            // First facet for this edge - store its index
            vor_facet_dual_to_edge[edge_id] = facetIndex;
            voronoiDiagram.cell_facets[facetIndex].mirror_facet_index = -1;
        }
        else
        {
            // Second facet for this edge - link both as mirrors
            int mirror = vor_facet_dual_to_edge[edge_id];
            voronoiDiagram.cell_facets[facetIndex].mirror_facet_index = mirror;
            voronoiDiagram.cell_facets[mirror].mirror_facet_index = facetIndex;
        }
    }

    return true;
}

//! @brief Processes incident edges to build facets for a Voronoi cell.
/*!
 * Iterates over incident edges to construct facets and add them to the cell.
 *
 * @param dt The Delaunay triangulation.
 * @param delaunay_vertex The Delaunay vertex to process.
 * @param voronoiDiagram The Voronoi diagram to update.
 * @param vc The Voronoi cell to populate with facets.
 * @param vor_facet_dual_to_edge Vector mapping Delaunay edge index to first Voronoi facet index.
 */
static void process_incident_edges(
    Delaunay &dt,
    Vertex_handle delaunay_vertex,
    VoronoiDiagram &voronoiDiagram,
    VoronoiCell &vc,
    std::vector<int> &vor_facet_dual_to_edge)
{
    std::vector<Edge> incidentEdges;
    incidentEdges.reserve(32);
    dt.incident_edges(delaunay_vertex, std::back_inserter(incidentEdges));

    for (const Edge &ed : incidentEdges)
    {
        VoronoiCellFacet facet;
        bool ok = build_facet_from_edge(dt, ed, delaunay_vertex, voronoiDiagram, vc.facetIndices, vor_facet_dual_to_edge, vc.cellIndex, facet);
        if (!ok)
        {
            if (debug)
            {
                std::cout << "[WARNING] Facet construction failed for edge (likely boundary or degenerate)\n";
            }
            continue;
        }

        // Verify facet validity
        const size_t vertexCount = facet.verticesIndices.size();
        if (vertexCount < 3)
        {
            // Enhanced error logging for debugging degenerate facets
            std::cout << "[ERROR] Degenerate Voronoi facet detected:\n";
            std::cout << "  - Voronoi facet vertices count: " << vertexCount << "\n";

            // Extract edge endpoints
            Cell_handle cell_ed = ed.first;
            int edge_idx = ed.second;
            int i1, i2;
            // CGAL edge encoding: edge index determines which two vertices form the edge
            // Edge indices 0-5 correspond to the 6 edges of a tetrahedron
            static const int edge_vertices[6][2] = {{0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3}};
            i1 = edge_vertices[edge_idx][0];
            i2 = edge_vertices[edge_idx][1];
            Vertex_handle v1 = cell_ed->vertex(i1);
            Vertex_handle v2 = cell_ed->vertex(i2);

            std::cout << "  - Delaunay edge endpoints:\n";
            std::cout << "      v1: " << v1->point() << " (index: " << v1->info().index << ")\n";
            std::cout << "      v2: " << v2->point() << " (index: " << v2->info().index << ")\n";

            // Information about the Delaunay vertex (Voronoi cell center)
            std::cout << "  - Delaunay vertex (Voronoi cell center): " << delaunay_vertex->point()
                      << " (index: " << delaunay_vertex->info().index << ")\n";

            // Detailed incident cell information using Cell_handle << operator
            std::cout << "  - Incident Delaunay cells around Delaunay edge:\n";
            Delaunay::Cell_circulator cc_dbg = dt.incident_cells(ed);
            Delaunay::Cell_circulator start_dbg = cc_dbg;
            int cell_num = 0;
            do
            {
                std::cout << "      ===== Delaunay Cell " << cell_num << " =====\n";
                if (dt.is_infinite(cc_dbg))
                {
                    std::cout << "      INFINITE CELL\n";
                }
                else
                {
                    // Use the Cell_handle << operator for comprehensive output
                    Cell_handle ch = cc_dbg;
                    std::string cell_output;
                    std::ostringstream oss;
                    oss << ch;
                    cell_output = oss.str();

                    // Indent each line of the cell output
                    std::istringstream iss(cell_output);
                    std::string line;
                    while (std::getline(iss, line))
                    {
                        std::cout << "      " << line << "\n";
                    }
                }
                ++cc_dbg;
                ++cell_num;
            } while (cc_dbg != start_dbg);

            std::cout << "  - Voronoi cell index: " << vc.cellIndex << "\n";
            std::cout << std::flush;
            continue;
        }

        // facet index already appended by build_facet_from_edge
    }
}

//! @brief Retrieves the edge key for a pair of vertices.
static std::pair<int, int> get_edge_key(int u, int v)
{
    return {std::min(u, v), std::max(u, v)};
}

namespace
{
    //! @brief Maps each cell facet to its owning cell index (or -1).
    std::vector<int> build_facet_to_cell_map(const VoronoiDiagram &vd)
    {
        std::vector<int> facetToCell(vd.cell_facets.size(), -1);
        for (size_t cellIdx = 0; cellIdx < vd.cells.size(); ++cellIdx)
        {
            for (int fi : vd.cells[cellIdx].facetIndices)
            {
                if (fi >= 0 && fi < static_cast<int>(facetToCell.size()))
                    facetToCell[fi] = static_cast<int>(cellIdx);
            }
        }
        return facetToCell;
    }

    //! @brief Within a single cell, ensure neighboring facets disagree on the shared edge orientation.
    void propagate_facets_within_cell(size_t cellIdx, VoronoiDiagram &vd)
    {
        auto &cell = vd.cells[cellIdx];
        if (cell.facetIndices.empty())
            return;

        const size_t num_facets = cell.facetIndices.size();
        std::vector<std::map<size_t, std::pair<int, int>>> adjacency(num_facets);

        for (size_t i = 0; i < num_facets; ++i)
        {
            const int f1 = cell.facetIndices[i];
                const auto &verts1 = vd.cell_facets[f1].verticesIndices;
            std::map<std::pair<int, int>, size_t> edges1;
            for (size_t j = 0; j < verts1.size(); ++j)
            {
                const int u = verts1[j];
                const int v = verts1[(j + 1) % verts1.size()];
                edges1[get_edge_key(u, v)] = j;
            }

            for (size_t k = i + 1; k < num_facets; ++k)
            {
                const int f2 = cell.facetIndices[k];
                const auto &verts2 = vd.cell_facets[f2].verticesIndices;
                std::pair<int, int> shared = {-1, -1};
                int shared_count = 0;

                for (size_t j = 0; j < verts2.size(); ++j)
                {
                    const int u = verts2[j];
                    const int v = verts2[(j + 1) % verts2.size()];
                    auto key = get_edge_key(u, v);
                    if (edges1.count(key))
                    {
                        shared = key;
                        if (++shared_count > 1)
                            break; // only consider true neighbors sharing a single edge
                    }
                }

                if (shared_count == 1)
                {
                    adjacency[i][k] = shared;
                    adjacency[k][i] = shared;
                }
            }
        }

        std::vector<bool> visited(num_facets, false);
        std::queue<size_t> q;
        q.push(0);
        visited[0] = true;

        while (!q.empty())
        {
            const size_t curr = q.front();
            q.pop();

            const int currFacetIdx = cell.facetIndices[curr];
            auto &currVerts = vd.cell_facets[currFacetIdx].verticesIndices;

            for (const auto &entry : adjacency[curr])
            {
                const size_t next = entry.first;
                if (visited[next])
                    continue;

                visited[next] = true;
                q.push(next);

                const std::pair<int, int> shared = entry.second;
                const int nextFacetIdx = cell.facetIndices[next];
                auto &nextVerts = vd.cell_facets[nextFacetIdx].verticesIndices;

                auto has_direction = [&](const std::vector<int> &verts) {
                    for (size_t j = 0; j < verts.size(); ++j)
                    {
                        const int a = verts[j];
                        const int b = verts[(j + 1) % verts.size()];
                        if (get_edge_key(a, b) == shared)
                            return (a == shared.first && b == shared.second);
                    }
                    return false;
                };

                const bool curr_dir = has_direction(currVerts);
                const bool next_dir = has_direction(nextVerts);

                if (curr_dir == next_dir)
                    std::reverse(nextVerts.begin(), nextVerts.end());
            }
        }
    }

    //! @brief Ensure each undirected edge in the cell appears twice with opposite direction.
    bool audit_cell_edge_orientation(size_t cellIdx, const VoronoiDiagram &vd)
    {
        const auto &cell = vd.cells[cellIdx];
        if (cell.facetIndices.empty())
            return true;

        struct EdgeInfo
        {
            int count = 0;
            int orientationSum = 0;
        };

        std::map<std::pair<int, int>, EdgeInfo> usage;

        for (int fi : cell.facetIndices)
        {
            const auto &verts = vd.cell_facets[fi].verticesIndices;
            const size_t n = verts.size();
            if (n < 2)
                continue;

            for (size_t k = 0; k < n; ++k)
            {
                const int u = verts[k];
                const int v = verts[(k + 1) % n];
                if (u == v)
                    continue;

                const std::pair<int, int> key{std::min(u, v), std::max(u, v)};
                const int step = (key.first == u && key.second == v) ? +1 : -1;
                EdgeInfo &info = usage[key];
                info.count += 1;
                info.orientationSum += step;
            }
        }

        for (const auto &kv : usage)
        {
            const EdgeInfo &info = kv.second;
            if (info.count != 2 || info.orientationSum != 0)
            {
                if (debug)
                {
                    std::cerr << "[DEBUG] Cell " << cell.cellIndex
                              << " has inconsistent facet orientation along edge {"
                              << kv.first.first << "," << kv.first.second << "}: count="
                              << info.count << " orientSum=" << info.orientationSum << "\n";
                }
                return false;
            }
        }
        return true;
    }

    void flip_cell(VoronoiDiagram &vd, int cellIdx)
    {
        auto &cell = vd.cells[cellIdx];
        for (int fi : cell.facetIndices)
        {
            auto &verts = vd.cell_facets[fi].verticesIndices;
            std::reverse(verts.begin(), verts.end());
        }
    }

    //! @brief Orient the cell so its facets face away from the Delaunay site.
    void orient_cell_outward(VoronoiDiagram &vd, int cellIdx)
    {
        auto &cell = vd.cells[cellIdx];
        for (int fi : cell.facetIndices)
        {
            const auto &verts = vd.cell_facets[fi].verticesIndices;
            const size_t n = verts.size();
            if (n < 3)
                continue;

            Point centroid(0, 0, 0);
            const double invN = 1.0 / static_cast<double>(n);
            for (int idx : verts)
                centroid = centroid + (vd.vertices[idx].coord - CGAL::ORIGIN) * invN;

            Vector3 normal(0, 0, 0);
            for (size_t k = 0; k < n; ++k)
            {
                const Point &p1 = vd.vertices[verts[k]].coord;
                const Point &p2 = vd.vertices[verts[(k + 1) % n]].coord;
                normal = normal + Vector3(
                                      (p1.y() - p2.y()) * (p1.z() + p2.z()),
                                      (p1.z() - p2.z()) * (p1.x() + p2.x()),
                                      (p1.x() - p2.x()) * (p1.y() + p2.y()));
            }
            normal = normal / 2.0;

            if (normal.squared_length() <= 1e-12)
                continue;

            const Vector3 v = cell.delaunayVertex->point() - centroid;
            if (CGAL::scalar_product(normal, v) > 0)
                flip_cell(vd, cellIdx);
            return;
        }
    }

    //! @brief Spread orientations across neighboring cells via mirror facets.
    void propagate_orientation_between_cells(VoronoiDiagram &vd, const std::vector<int> &facetToCell)
    {
        std::vector<int> cellMark(vd.cells.size(), 0);
        std::queue<int> pending;

        for (size_t root = 0; root < vd.cells.size(); ++root)
        {
            if (vd.cells[root].facetIndices.empty())
                continue;
            if (cellMark[root] != 0)
                continue;

            orient_cell_outward(vd, static_cast<int>(root));
            cellMark[root] = 1;
            pending.push(static_cast<int>(root));

            while (!pending.empty())
            {
                const int curr = pending.front();
                pending.pop();

                for (int fi : vd.cells[curr].facetIndices)
                {
                    const int mirror = vd.cell_facets[fi].mirror_facet_index;
                    if (mirror < 0 || mirror >= static_cast<int>(vd.cell_facets.size()))
                        continue;

                    const int neighbor = facetToCell[mirror];
                    if (neighbor < 0 || neighbor == curr)
                        continue;

                    bool opposite = vd.haveOppositeOrientation(
                        vd.cell_facets[fi].verticesIndices,
                        vd.cell_facets[mirror].verticesIndices);

                    if (cellMark[neighbor] == 0)
                    {
                        if (!opposite)
                        {
                            flip_cell(vd, neighbor);
                            opposite = vd.haveOppositeOrientation(
                                vd.cell_facets[fi].verticesIndices,
                                vd.cell_facets[mirror].verticesIndices);
                            if (!opposite)
                            {
                                throw std::runtime_error(
                                    "Unable to orient cell " + std::to_string(neighbor) +
                                    " opposite to cell " + std::to_string(curr) + ".");
                            }
                        }
                        cellMark[neighbor] = -cellMark[curr];
                        pending.push(neighbor);
                    }
                    else if (!opposite)
                    {
                        throw std::runtime_error(
                            "Cells " + std::to_string(curr) + " and " + std::to_string(neighbor) +
                            " disagree on shared facet orientation.");
                    }
                }
            }
        }
    }
} // namespace

//! @brief Validates facet orientation using only combinatorial information with a single geometric anchor per component.
void validate_facet_orientations_and_normals(VoronoiDiagram &voronoiDiagram)
{
    const auto facetToCell = build_facet_to_cell_map(voronoiDiagram);

    for (size_t cellIdx = 0; cellIdx < voronoiDiagram.cells.size(); ++cellIdx)
    {
        const auto &cell = voronoiDiagram.cells[cellIdx];
        if (cell.facetIndices.size() < 2)
            continue;

        if (!audit_cell_edge_orientation(cellIdx, voronoiDiagram))
        {
            propagate_facets_within_cell(cellIdx, voronoiDiagram);
        }
    }

    propagate_orientation_between_cells(voronoiDiagram, facetToCell);
}

//! @brief Constructs Voronoi cells without using Convex_Hull_3 (in development).
/*!
 * Populates the Voronoi diagram with polyhedral cells derived from the Delaunay
 * triangulation by processing incident edges and cell circulators.
 *
 * Optimized to use a simple vector indexed by global Delaunay edge ID for mirror
 * facet linking, instead of an expensive hash map.
 *
 * @param voronoiDiagram The Voronoi diagram to populate with cells.
 * @param dt The Delaunay triangulation corresponding (dual) to the Voronoi diagram.
 */
void construct_voronoi_cells_from_delaunay_triangulation(VoronoiDiagram &voronoiDiagram, Delaunay &dt)
{
    TimingStats& timer = TimingStats::getInstance();

    // Assign global indices to all Delaunay edges (stored in CellInfo.edge_index)
    timer.startTimer("Assign Delaunay edge indices", "Construct Voronoi cells");
    int num_del_edges = assign_delaunay_edge_indices(dt);
    timer.stopTimer("Assign Delaunay edge indices", "Construct Voronoi cells");

    // Use a vector indexed by edge ID to track Voronoi facets for mirror linking
    // -1 means no facet yet created for this edge
    std::vector<int> vor_facet_dual_to_edge(num_del_edges, -1);

    voronoiDiagram.cells.reserve(dt.number_of_vertices());
    int cellIndex = 0;

    timer.startTimer("Build Voronoi cells", "Construct Voronoi cells");
    for (Vertex_handle v : dt.finite_vertex_handles())
    {
        if (v->info().is_dummy)
            continue;

        VoronoiCell vc = create_voronoi_cell(v, cellIndex);
        collect_cell_vertices(dt, v, voronoiDiagram, vc.verticesIndices);
        process_incident_edges(dt, v, voronoiDiagram, vc, vor_facet_dual_to_edge);

        if (vc.facetIndices.size() < 4)
        {
            std::cout << "[WARNING] Cell " << cellIndex << " has only " << vc.facetIndices.size() << " facets, skipping\n";
        }
        else
        {
            voronoiDiagram.cells.push_back(std::move(vc));
            v->info().voronoiCellIndex = cellIndex;
            cellIndex++;
        }
    }
    timer.stopTimer("Build Voronoi cells", "Construct Voronoi cells");

    // Mirror facet linking is now done inline in build_facet_from_edge()
    // No need for final loop over EdgeFacetMap
}



//! @brief Helper function to find the mirror facet index in neighboring cell.
/*!
 * Given two neighboring cells c1 and c2, finds the facet index in c2 that
 * corresponds to the shared facet (i.e., c2's facet that points back to c1).
 *
 * @param c1 The first cell.
 * @param c2 The second (neighboring) cell.
 * @return The facet index in c2, or -1 if not found.
 */
static int find_mirror_facet_index(Cell_handle c1, Cell_handle c2)
{
    for (int i = 0; i < 4; ++i)
    {
        if (c2->neighbor(i) == c1)
        {
            return i;
        }
    }
    return -1;
}

//! @brief Constructs Voronoi edges from Delaunay facets.
/*!
 * Optimized implementation that uses the dualEdgeIndex field in DelaunayFacetInfo
 * to detect if an edge has already been created, avoiding expensive map lookups.
 * Each Delaunay facet has a mirror facet in the neighboring cell. When we first
 * encounter a facet pair, we create the edge and mark both facets with the edge index.
 * When we later encounter the mirror facet, dualEdgeIndex is already set, so we skip it.
 */
void construct_voronoi_edges(VoronoiDiagram &voronoiDiagram, Delaunay &dt)
{
    voronoiDiagram.edges.clear();

    for (Delaunay::Finite_facets_iterator fit = dt.finite_facets_begin(); fit != dt.finite_facets_end(); ++fit)
    {
        Facet facet = *fit;
        Cell_handle c1 = facet.first;
        int facet1_index = facet.second;

        // Check if this facet's edge was already created (via its mirror facet)
        int existingEdgeIdx = c1->info().facet_info[facet1_index].dualEdgeIndex;
        if (existingEdgeIdx >= 0)
        {
            // Edge already exists - this facet was processed as a mirror of another
            // Check if we need to update the canonical facet selection
            VoronoiEdge &vEdge = voronoiDiagram.edges[existingEdgeIdx];
            if (!vEdge.delaunayFacets.empty())
            {
                int existingFacetIdx = vEdge.delaunayFacets[0].second;
                int newFacetIdx = facet.second;
                // Choose the facet with smallest index for canonical orientation
                if (newFacetIdx < existingFacetIdx)
                {
                    vEdge.delaunayFacets[0] = facet;
                }
            }
            continue;
        }

        // This is the first time we see this facet pair - create the edge
        CGAL::Object edgeobj = dt.dual(facet);
        Segment3 seg;
        Ray3 ray;

        Cell_handle c2 = c1->neighbor(facet1_index);
        int facet2_index = find_mirror_facet_index(c1, c2);

        if (CGAL::assign(seg, edgeobj))
        {
            int idx1 = dt.is_infinite(c1) ? -1 : c1->info().dualVoronoiVertexIndex;
            int idx2 = dt.is_infinite(c2) ? -1 : c2->info().dualVoronoiVertexIndex;

            if (idx1 != -1 && idx2 != -1 && idx1 != idx2)
            {
                // Finite segment - create new edge
                int v1 = std::min(idx1, idx2);
                int v2 = std::max(idx1, idx2);

                VoronoiEdge vEdge(edgeobj);
                vEdge.type = 0;
                vEdge.vertex1 = v1;
                vEdge.vertex2 = v2;
                int edgeIdx = voronoiDiagram.edges.size();
                vEdge.delaunayFacets.push_back(facet);
                voronoiDiagram.edges.push_back(vEdge);

                // Populate incident edges for both vertices
                voronoiDiagram.vertices[v1].incidentEdgeIndices.push_back(edgeIdx);
                voronoiDiagram.vertices[v2].incidentEdgeIndices.push_back(edgeIdx);

                // Store edge index in both Delaunay cells (this and mirror facet)
                c1->info().facet_info[facet1_index].dualEdgeIndex = edgeIdx;
                if (facet2_index >= 0)
                {
                    c2->info().facet_info[facet2_index].dualEdgeIndex = edgeIdx;
                }
            }
        }
        else if (CGAL::assign(ray, edgeobj))
        {
            int vertex1 = -1;
            if (!dt.is_infinite(c1))
            {
                vertex1 = c1->info().dualVoronoiVertexIndex;
            }
            else if (!dt.is_infinite(c2))
            {
                vertex1 = c2->info().dualVoronoiVertexIndex;
            }

            if (vertex1 != -1)
            {
                // Ray edge - create new edge
                Vector3 dir = ray.direction().vector();

                VoronoiEdge vEdge(edgeobj);
                vEdge.type = 1;
                vEdge.vertex1 = vertex1;
                vEdge.vertex2 = -1;
                vEdge.source = ray.source();
                vEdge.direction = dir;
                int edgeIdx = voronoiDiagram.edges.size();
                vEdge.delaunayFacets.push_back(facet);
                voronoiDiagram.edges.push_back(vEdge);

                // Store edge index in both Delaunay cells (this and mirror facet)
                c1->info().facet_info[facet1_index].dualEdgeIndex = edgeIdx;
                if (facet2_index >= 0)
                {
                    c2->info().facet_info[facet2_index].dualEdgeIndex = edgeIdx;
                }
            }
        }
        // Lines are not expected for finite facets
    }
}

//! @brief Stores cell edge index in all Delaunay cells for a given Voronoi cell.
/*!
 * For each Delaunay facet sharing a Voronoi edge, finds which vertex corresponds
 * to the target Voronoi cell and stores the cell edge index in the appropriate slot.
 *
 * @param dt The Delaunay triangulation.
 * @param sharedFacets Vector of Delaunay facets that share this Voronoi edge.
 * @param targetCellIdx The index of the Voronoi cell we're storing the edge for.
 * @param cell_edge_index The cell edge index to store.
 */
static void store_cell_edge_index_in_Delaunay_cell(
    Delaunay &dt,
    const std::vector<Facet> &sharedFacets,
    int targetCellIdx,
    int cell_edge_index)
{
    const int NUM_DELAUNAY_CELL_VERTICES = 4;
    for (const Facet &f : sharedFacets)
    {
        Cell_handle c = f.first;
        if (dt.is_infinite(c))
            continue;

        int facet_index = f.second;
        // Check each vertex of the facet (k = 1, 2, 3 gives the 3 vertices)
        for (int k = 1; k < NUM_DELAUNAY_CELL_VERTICES; k++)
        {
            const int kv = (facet_index + k) % NUM_DELAUNAY_CELL_VERTICES;
            Vertex_handle delaunay_vertex = c->vertex(kv);
            if (!delaunay_vertex->info().is_dummy)
            {
                int cellIdx = delaunay_vertex->info().voronoiCellIndex;
                if (cellIdx == targetCellIdx)
                {
                    c->info().facet_info[facet_index].dualCellEdgeIndex[k - 1] = cell_edge_index;
                }
            }
        }
    }
}

//! @brief Builds Voronoi cell edges for each edge in the diagram.
/*!
 * Optimized implementation that uses a single loop with a local set per edge,
 * creating VoronoiCellEdge objects and storing indices immediately when a new
 * cell is discovered, rather than using a two-pass approach with intermediate storage.
 *
 * @param voronoiDiagram The Voronoi diagram to populate with cell edges.
 * @param dt The Delaunay triangulation.
 */
static void build_cell_edges(
    VoronoiDiagram &voronoiDiagram,
    Delaunay &dt)
{
    voronoiDiagram.cellEdges.clear();

    for (size_t edgeIdx = 0; edgeIdx < voronoiDiagram.edges.size(); ++edgeIdx)
    {
        const std::vector<Facet> &sharedFacets = voronoiDiagram.edges[edgeIdx].delaunayFacets;
        std::vector<int> seenCells; // Local per-edge, small; avoids hash overhead
        seenCells.reserve(8);

        for (const Facet &f : sharedFacets)
        {
            Cell_handle c = f.first;
            if (dt.is_infinite(c))
                continue;

            int opp = f.second; // Opposite vertex index
            for (int corner = 0; corner < 4; ++corner)
            {
                if (corner == opp)
                    continue; // Skip opposite, add only facet's 3 vertices

                Vertex_handle delaunay_vertex = c->vertex(corner);
                if (!delaunay_vertex->info().is_dummy)
                {
                    int cellIdx = delaunay_vertex->info().voronoiCellIndex;

                    // Check if this cell was already processed for this edge
                    if (std::find(seenCells.begin(), seenCells.end(), cellIdx) == seenCells.end())
                    {
                        seenCells.push_back(cellIdx);

                        // Create cell edge immediately
                        VoronoiCellEdge cellEdge;
                        cellEdge.cellIndex = cellIdx;
                        cellEdge.edgeIndex = static_cast<int>(edgeIdx);
                        cellEdge.cycleIndices = {};
                        cellEdge.nextCellEdge = -1;
                        const int cell_edge_index = voronoiDiagram.cellEdges.size();
                        voronoiDiagram.cellEdges.push_back(cellEdge);

                        // Store in Delaunay cells
                        store_cell_edge_index_in_Delaunay_cell(
                            dt, sharedFacets, cellIdx, cell_edge_index);
                    }
                }
            }
        }
    }
}


//! @brief Links Voronoi cell edges in a circular ring.
/*!
 * Connects cell edges sharing the same edge index using the nextCellEdge field
 * to form a closed loop.
 *
 * @param voronoiDiagram The Voronoi diagram containing cell edges to link.
 */
static void link_cell_edges(
    VoronoiDiagram &voronoiDiagram)
{
    const size_t edgeCount = voronoiDiagram.edges.size();
    std::vector<std::vector<int>> edge_to_cellEdges(edgeCount);

    for (int ceIdx = 0; ceIdx < static_cast<int>(voronoiDiagram.cellEdges.size()); ++ceIdx)
    {
        const VoronoiCellEdge &ce = voronoiDiagram.cellEdges[ceIdx];
        if (ce.edgeIndex < 0 || ce.edgeIndex >= static_cast<int>(edgeCount))
            continue;
        edge_to_cellEdges[static_cast<size_t>(ce.edgeIndex)].push_back(ceIdx);
    }

    for (auto &cellEdgeIndices : edge_to_cellEdges)
    {
        const int N = static_cast<int>(cellEdgeIndices.size());
        if (N == 0)
            continue;
        for (int i = 0; i < N; ++i)
        {
            const int ceIdx = cellEdgeIndices[i];
            const int nextIdx = cellEdgeIndices[(i + 1) % N]; // ring (N==1 => self loop)
            voronoiDiagram.cellEdges[ceIdx].nextCellEdge = nextIdx;
        }
    }
}

//! @brief Processes edge mapping for a single Voronoi edge.
/*!
 * Updates the vertex incident edge lists for segments after intersecting with the bounding box.
 *
 * @param voronoiDiagram The Voronoi diagram to update.
 * @param edgeObj The CGAL object representing the edge.
 * @param edgeIdx The index of the edge in the diagram.
 * @param bbox The bounding box for intersection.
 */
static void process_edge_mapping(VoronoiDiagram &voronoiDiagram, VoronoiEdge &edge, int edgeIdx, CGAL::Epick::Iso_cuboid_3 &bbox)
{
    if (edge.type == 0) // Only process finite segments combinatorially
    {
        int idx1 = edge.vertex1;
        int idx2 = edge.vertex2;
        if (idx1 != -1 && idx2 != -1)
        {
            // Update incident edges for both vertices if not already present
            auto& incidentEdges1 = voronoiDiagram.vertices[idx1].incidentEdgeIndices;
            if (std::find(incidentEdges1.begin(), incidentEdges1.end(), edgeIdx) == incidentEdges1.end()) {
                incidentEdges1.push_back(edgeIdx);
            }
            auto& incidentEdges2 = voronoiDiagram.vertices[idx2].incidentEdgeIndices;
            if (std::find(incidentEdges2.begin(), incidentEdges2.end(), edgeIdx) == incidentEdges2.end()) {
                incidentEdges2.push_back(edgeIdx);
            }
        }
    }
    // Rays and lines are skipped; no mapping for infinite edges
}

//! @brief Updates edge mappings for all Voronoi edges.
/*!
 * Processes all edges to update vertex incident edge lists.
 *
 * @param voronoiDiagram The Voronoi diagram to update.
 * @param bbox The bounding box for intersection.
 */
static void update_edge_mapping(VoronoiDiagram &voronoiDiagram, CGAL::Epick::Iso_cuboid_3 &bbox)
{
    for (int edgeIdx = 0; edgeIdx < static_cast<int>(voronoiDiagram.edges.size()); ++edgeIdx)
    {
        VoronoiEdge &edge = voronoiDiagram.edges[edgeIdx];
        process_edge_mapping(voronoiDiagram, edge, edgeIdx, bbox);
    }
}

//! @brief Constructs the Voronoi cell edges in the Voronoi diagram and links them.
/*!
 * Builds cell edges for each Voronoi edge, links them in a circular ring, and
 * updates edge mappings.
 *
 * @param voronoiDiagram The Voronoi diagram to populate with edges.
 * @param bbox The bounding box used for clipping rays and lines.
 * @param dt The Delaunay triangulation.
 */
void construct_voronoi_cell_edges(
    VoronoiDiagram &voronoiDiagram,
    CGAL::Epick::Iso_cuboid_3 &bbox,
    Delaunay &dt)
{
    TimingStats& timer = TimingStats::getInstance();
    voronoiDiagram.cellEdges.clear();

    // Note: Parent timer "Construct cell edges" is started in construct_voronoi_diagram()
    timer.startTimer("Build cell edges", "Construct cell edges");
    build_cell_edges(voronoiDiagram, dt);
    timer.stopTimer("Build cell edges", "Construct cell edges");

    // No need to populate per-cell anchors; ring traversal provides lookup.

    timer.startTimer("Link cell edges", "Construct cell edges");
    link_cell_edges(voronoiDiagram);
    timer.stopTimer("Link cell edges", "Construct cell edges");

    timer.startTimer("Update edge mapping", "Construct cell edges");
    update_edge_mapping(voronoiDiagram, bbox);
    timer.stopTimer("Update edge mapping", "Construct cell edges");
}

// Resolve (cellIndex, globalEdgeIndex) → VoronoiCellEdge index
int find_cell_edge_for_cell_and_edge(const VoronoiDiagram &vd,
                                     int cellIndex,
                                     int globalEdgeIndex)
{
    if (cellIndex < 0 || cellIndex >= static_cast<int>(vd.cells.size()))
        return -1;
    if (globalEdgeIndex < 0 || globalEdgeIndex >= static_cast<int>(vd.edges.size()))
        return -1;

    // Fallback: build a per-diagram table from global edge → ring anchor once.
    {
        struct Memo {
            const VoronoiDiagram *vd_ptr = nullptr;
            std::vector<int> startByEdge; // size = vd.edges.size(), value = first ceIdx for that edge or -1
        };
        static Memo memo;

        if (memo.vd_ptr != &vd || memo.startByEdge.size() != vd.edges.size())
        {
            memo.vd_ptr = &vd;
            memo.startByEdge.assign(vd.edges.size(), -1);
            for (int i = 0; i < static_cast<int>(vd.cellEdges.size()); ++i)
            {
                const int e = vd.cellEdges[i].edgeIndex;
                if (e >= 0 && e < static_cast<int>(memo.startByEdge.size()) && memo.startByEdge[e] == -1)
                    memo.startByEdge[e] = i;
            }
        }

        int startIdx = (globalEdgeIndex >= 0 && globalEdgeIndex < (int)memo.startByEdge.size())
                            ? memo.startByEdge[globalEdgeIndex]
                            : -1;
        if (startIdx >= 0)
        {
            int ceIdx = startIdx;
            const int start = startIdx;
            for (;;)
            {
                const VoronoiCellEdge &ce = vd.cellEdges[ceIdx];
                if (ce.cellIndex == cellIndex)
                    return ceIdx;
                const int nxt = ce.nextCellEdge;
                if (nxt < 0 || nxt == start)
                    break;
                ceIdx = nxt;
            }
        }
    }

    return -1;
}


//! @brief Wrap up function of constructing voronoi diagram
void construct_voronoi_diagram(VoronoiDiagram &vd, VdcParam &vdc_param, UnifiedGrid &grid, CGAL::Epick::Iso_cuboid_3 &bbox, Delaunay &dt)
{
    TimingStats& timer = TimingStats::getInstance();

    std::cout << "[INFO] Start constructing Voronoi vertices and edges..." << std::endl;
    timer.startTimer("Construct Voronoi vertices", "4. Voronoi Diagram Construction");
    construct_voronoi_vertices(vd, dt);
    timer.stopTimer("Construct Voronoi vertices", "4. Voronoi Diagram Construction");

    timer.startTimer("Construct Voronoi edges", "4. Voronoi Diagram Construction");
    construct_voronoi_edges(vd, dt);
    timer.stopTimer("Construct Voronoi edges", "4. Voronoi Diagram Construction");

    timer.startTimer("Compute vertex values", "4. Voronoi Diagram Construction");
    compute_voronoi_values(vd, grid);
    timer.stopTimer("Compute vertex values", "4. Voronoi Diagram Construction");

    std::cout << "[INFO] Start constructing Voronoi Cells" << std::endl;
    if (vdc_param.multi_isov)
    {
        if (vdc_param.convex_hull)
        {
            timer.startTimer("Construct Voronoi cells", "4. Voronoi Diagram Construction");
            construct_voronoi_cells_as_convex_hull(vd, dt);
            timer.stopTimer("Construct Voronoi cells", "4. Voronoi Diagram Construction");
        }
        else
        {
            timer.startTimer("Construct Voronoi cells", "4. Voronoi Diagram Construction");
            construct_voronoi_cells_from_delaunay_triangulation(vd, dt);
            timer.stopTimer("Construct Voronoi cells", "4. Voronoi Diagram Construction");

            timer.startTimer("Validate facet orientations", "4. Voronoi Diagram Construction");
            validate_facet_orientations_and_normals(vd);
            timer.stopTimer("Validate facet orientations", "4. Voronoi Diagram Construction");
        }

        timer.startTimer("Construct cell edges", "4. Voronoi Diagram Construction");
        construct_voronoi_cell_edges(vd, bbox, dt);
        timer.stopTimer("Construct cell edges", "4. Voronoi Diagram Construction");

        timer.startTimer("Create global facets", "4. Voronoi Diagram Construction");
        vd.create_global_facets();
        timer.stopTimer("Create global facets", "4. Voronoi Diagram Construction");
    }

    if (!vdc_param.no_check)
    {
        timer.startTimer("VD initial check", "4. Voronoi Diagram Construction");
        vd.check(false);
        timer.stopTimer("VD initial check", "4. Voronoi Diagram Construction");
    }
}
