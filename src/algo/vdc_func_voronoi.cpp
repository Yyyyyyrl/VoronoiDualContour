#include "algo/vdc_func.h"
#include "core/vdc_debug.h"
#include "core/vdc_timing.h"


//! @brief Constructs Voronoi vertices for the given voronoi Diagram instance.
void construct_voronoi_vertices(VoronoiDiagram &voronoiDiagram, Delaunay &dt)
{
    int vertexIndex = 0;
    voronoiDiagram.vertices.reserve(dt.number_of_finite_cells());
    for (Delaunay::Finite_cells_iterator cit = dt.finite_cells_begin(); cit != dt.finite_cells_end(); ++cit)
    {
        cit->info().dualVoronoiVertexIndex = -1; // Default to invalid until a vertex is recorded
        // Skip degenerate cells to avoid invalid circumcenters (e.g., NaN coordinates)
        if (is_degenerate(cit))
        {
            if (debug)
            {
                std::cerr << "[DEBUG] Skipping degenerate cell with vertices: "
                          << cit->vertex(0)->point() << ", "
                          << cit->vertex(1)->point() << ", "
                          << cit->vertex(2)->point() << ", "
                          << cit->vertex(3)->point() << "\n";
            }
            continue;
        }

        Point circumcenter = dt.dual(cit);

        if (std::isnan(circumcenter.x()) || std::isnan(circumcenter.y()) || std::isnan(circumcenter.z()))
        {
            if (debug)
            {
                std::cerr << "[DEBUG] Skipping cell with NaN circumcenter: vertices "
                          << cit->vertex(0)->point() << ", "
                          << cit->vertex(1)->point() << ", "
                          << cit->vertex(2)->point() << ", "
                          << cit->vertex(3)->point() << "\n";
            }
            continue;
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
        vc.vertices_indices.assign(unique_vertex_indices_set.begin(), unique_vertex_indices_set.end());

        // Build vertex_points and vector for lookup (allow duplicates by using first idx for matching points)
        std::vector<Point> vertex_points;
        std::vector<std::pair<Point, int>> point_index_pairs;
        for (int idx : vc.vertices_indices)
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
                        vf.vertices_indices.push_back(pair.second);
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
            if (vf.vertices_indices.size() < 3)
                continue;

            int facet_index = voronoiDiagram.facets.size();
            voronoiDiagram.facets.push_back(vf);
            vc.facet_indices.push_back(facet_index);
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
static void collcet_cell_vertices(
    Delaunay &dt,
    Vertex_handle delaunay_vertex,
    VoronoiDiagram &voronoiDiagram,
    std::vector<int> &vertices_indices)
{
    std::vector<Cell_handle> incidentCells;
    dt.finite_incident_cells(delaunay_vertex, std::back_inserter(incidentCells));

    std::set<int> uniqueVertexIndices;
    for (Cell_handle c : incidentCells)
    {
        int vertex_index = c->info().dualVoronoiVertexIndex;
        if (vertex_index >= 0)
        {
            uniqueVertexIndices.insert(vertex_index);
        }
    }
    vertices_indices.assign(uniqueVertexIndices.begin(), uniqueVertexIndices.end());
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
 * ordering them cyclically, and assigning scalar values.
 *
 * @param dt The Delaunay triangulation.
 * @param ed The incident edge to process.
 * @param delaunay_vertex The Delaunay vertex associated with the cell.
 * @param voronoiDiagram The Voronoi diagram containing vertex and value data.
 * @param facet_indices Vector to store the facet index.
 * @return The constructed Voronoi facet, or an empty facet if invalid.
 */
static VoronoiCellFacet build_facet_from_edge(
    Delaunay &dt,
    const Edge &ed,
    Vertex_handle delaunay_vertex,
    VoronoiDiagram &voronoiDiagram,
    std::vector<int> &facet_indices,
    std::map<std::pair<int, int>, std::vector<int>> &edge_to_facets,
    int vcIdx)
{
    Cell_handle cell_ed = ed.first;
    int i = ed.second;
    int j = ed.third;
    Vertex_handle v1 = cell_ed->vertex(i);
    Vertex_handle v2 = cell_ed->vertex(j);

    Delaunay::Cell_circulator cc = dt.incident_cells(ed);
    Delaunay::Cell_circulator start = cc;
    std::vector<int> facetVertices;
    int finite_cell_count = 0;
    do
    {
        if (!dt.is_infinite(cc))
        {
            int vertex_index = cc->info().dualVoronoiVertexIndex;
            if (vertex_index >= 0)
            {
                facetVertices.push_back(vertex_index);
                ++finite_cell_count;
            }
        }
        ++cc;
    } while (cc != start);

    if (facetVertices.size() < 3)
    {
        if (debug)
        {
            std::cout << "[DEBUG] Degenerate facet for edge with " << finite_cell_count << " finite cells\n";
        }
        return VoronoiCellFacet();
    }

    std::set<int> unique_vertices(facetVertices.begin(), facetVertices.end());
    if (unique_vertices.size() < 3)
    {
        if (debug)
        {
            std::cout << "[DEBUG] Degenerate facet for edge with " << finite_cell_count << " finite cells\n";
        }
        return VoronoiCellFacet();
    }

    VoronoiCellFacet facet;
    facet.vertices_indices = std::move(facetVertices);

    int facetIndex = voronoiDiagram.facets.size();
    voronoiDiagram.facets.push_back(facet);
    facet_indices.push_back(facetIndex);
    std::pair<int, int> edge_key = std::make_pair(
        std::min(v1->info().index, v2->info().index),
        std::max(v1->info().index, v2->info().index));
    edge_to_facets[edge_key].push_back(facetIndex);

    // Map vertex loop → cell-edge loop for this cell facet
    // Optimized approach: Single-pass iteration using stored facet_info indices
    // Following the pattern from TD/build_facet_from_edge.cpp
    facet.cell_edge_indices.clear();
    const int n = (int)facet.vertices_indices.size();
    facet.cell_edge_indices.reserve(n);

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
        const Cell_handle cc = delFacet_circ->first;
        if (!dt.is_infinite(cc))
        {
            const int vor_vertex_index = cc->info().dualVoronoiVertexIndex;
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
    for (int i = 0; i < n; ++i)
    {
        const int a = facet.vertices_indices[i];

        // Find the Delaunay facet corresponding to Voronoi vertex 'a'
        // Use linear search since n is typically small (3-8 elements)
        bool found = false;
        for (size_t j = 0; j < vorVertexIndices.size(); ++j)
        {
            if (vorVertexIndices[j] == a)
            {
                const Facet &delFacet = delaunayFacets[j];
                // Get the cell edge index directly from stored facet_info
                const int dual_cell_edge_index = get_dual_cell_edge_index(delFacet, v1, v2);
                facet.cell_edge_indices.push_back(dual_cell_edge_index);
                found = true;
                break;
            }
        }

        if (!found)
        {
            // This should rarely happen - only for boundary/degenerate cases
            if (debug)
            {
                std::cerr << "[WARNING] build_facet_from_edge: Voronoi vertex " << a
                          << " not found in incident facets for edge ("
                          << v1->info().index << ", " << v2->info().index << ")\n";
            }
            facet.cell_edge_indices.push_back(-1);
        }
    }
    return facet;
}

//! @brief Processes incident edges to build facets for a Voronoi cell.
/*!
 * Iterates over incident edges to construct facets and add them to the cell.
 *
 * @param dt The Delaunay triangulation.
 * @param delaunay_vertex The Delaunay vertex to process.
 * @param voronoiDiagram The Voronoi diagram to update.
 * @param vc The Voronoi cell to populate with facets.
 */
static void process_incident_edges(
    Delaunay &dt,
    Vertex_handle delaunay_vertex,
    VoronoiDiagram &voronoiDiagram,
    VoronoiCell &vc,
    std::map<std::pair<int, int>, std::vector<int>> &edge_to_facets)
{
    std::vector<Edge> incidentEdges;
    dt.incident_edges(delaunay_vertex, std::back_inserter(incidentEdges));

    for (const Edge &ed : incidentEdges)
    {
        // Count finite incident cells
        int finite_cell_count = 0;
        Delaunay::Cell_circulator cc = dt.incident_cells(ed);
        Delaunay::Cell_circulator start = cc;
        do
        {
            if (!dt.is_infinite(cc))
            {
                finite_cell_count++;
            }
            ++cc;
        } while (cc != start);

        if (finite_cell_count < 3)
        {
            // std::cout << "[INFO] Skipping edge with " << finite_cell_count << " finite incident cells (insufficient for interior facet)\n";
            continue;
        }

        // Build facet only if edge has 3+ finite cells
        VoronoiCellFacet facet = build_facet_from_edge(dt, ed, delaunay_vertex, voronoiDiagram, vc.facet_indices, edge_to_facets, vc.cellIndex);
        if (facet.vertices_indices.empty() && debug)
        {
            std::cout << "[WARNING] Facet construction failed for edge with " << finite_cell_count << " finite cells\n";
            continue;
        }

        // Verify facet validity
        if (facet.vertices_indices.size() < 3)
        {
            // Enhanced error logging for debugging degenerate facets
            std::cout << "[ERROR] Degenerate interior facet detected:\n";
            std::cout << "  - Facet vertices count: " << facet.vertices_indices.size() << "\n";
            std::cout << "  - Finite incident cells: " << finite_cell_count << "\n";

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

            std::cout << "  - Edge endpoints:\n";
            std::cout << "      v1: " << v1->point() << " (index: " << v1->info().index << ")\n";
            std::cout << "      v2: " << v2->point() << " (index: " << v2->info().index << ")\n";

            // Information about the Delaunay vertex (Voronoi cell center)
            std::cout << "  - Delaunay vertex (Voronoi cell): " << delaunay_vertex->point()
                      << " (index: " << delaunay_vertex->info().index << ")\n";

            // Detailed incident cell information
            std::cout << "  - Incident cells around edge:\n";
            Delaunay::Cell_circulator cc_dbg = dt.incident_cells(ed);
            Delaunay::Cell_circulator start_dbg = cc_dbg;
            int cell_num = 0;
            do
            {
                std::cout << "      Cell " << cell_num << ": ";
                if (dt.is_infinite(cc_dbg))
                {
                    std::cout << "INFINITE\n";
                }
                else
                {
                    int vor_idx = cc_dbg->info().dualVoronoiVertexIndex;
                    std::cout << "dualVoronoiVertexIndex=" << vor_idx;
                    if (vor_idx >= 0)
                    {
                        std::cout << " (circumcenter: " << dt.dual(cc_dbg) << ")";
                    }
                    else
                    {
                        std::cout << " [INVALID - no Voronoi vertex assigned]";
                    }
                    std::cout << "\n";
                }
                ++cc_dbg;
                ++cell_num;
            } while (cc_dbg != start_dbg);

            std::cout << "  - Voronoi cell index: " << vc.cellIndex << "\n";
            std::cout << std::flush;
            continue;
        }

        int facetIndex = voronoiDiagram.facets.size() - 1; // Assuming facet was just added
        if (std::find(vc.facet_indices.begin(), vc.facet_indices.end(), facetIndex) == vc.facet_indices.end())
        {
            vc.facet_indices.push_back(facetIndex);
        }
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
        std::vector<int> facetToCell(vd.facets.size(), -1);
        for (size_t cellIdx = 0; cellIdx < vd.cells.size(); ++cellIdx)
        {
            for (int fi : vd.cells[cellIdx].facet_indices)
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
        if (cell.facet_indices.empty())
            return;

        const size_t num_facets = cell.facet_indices.size();
        std::vector<std::map<size_t, std::pair<int, int>>> adjacency(num_facets);

        for (size_t i = 0; i < num_facets; ++i)
        {
            const int f1 = cell.facet_indices[i];
            const auto &verts1 = vd.facets[f1].vertices_indices;
            std::map<std::pair<int, int>, size_t> edges1;
            for (size_t j = 0; j < verts1.size(); ++j)
            {
                const int u = verts1[j];
                const int v = verts1[(j + 1) % verts1.size()];
                edges1[get_edge_key(u, v)] = j;
            }

            for (size_t k = i + 1; k < num_facets; ++k)
            {
                const int f2 = cell.facet_indices[k];
                const auto &verts2 = vd.facets[f2].vertices_indices;
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

            const int currFacetIdx = cell.facet_indices[curr];
            auto &currVerts = vd.facets[currFacetIdx].vertices_indices;

            for (const auto &entry : adjacency[curr])
            {
                const size_t next = entry.first;
                if (visited[next])
                    continue;

                visited[next] = true;
                q.push(next);

                const std::pair<int, int> shared = entry.second;
                const int nextFacetIdx = cell.facet_indices[next];
                auto &nextVerts = vd.facets[nextFacetIdx].vertices_indices;

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
        if (cell.facet_indices.empty())
            return true;

        struct EdgeInfo
        {
            int count = 0;
            int orientationSum = 0;
        };

        std::map<std::pair<int, int>, EdgeInfo> usage;

        for (int fi : cell.facet_indices)
        {
            const auto &verts = vd.facets[fi].vertices_indices;
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
        for (int fi : cell.facet_indices)
        {
            auto &verts = vd.facets[fi].vertices_indices;
            std::reverse(verts.begin(), verts.end());
        }
    }

    //! @brief Orient the cell so its facets face away from the Delaunay site.
    void orient_cell_outward(VoronoiDiagram &vd, int cellIdx)
    {
        auto &cell = vd.cells[cellIdx];
        for (int fi : cell.facet_indices)
        {
            const auto &verts = vd.facets[fi].vertices_indices;
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

            const Vector3 v = cell.delaunay_vertex->point() - centroid;
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
            if (vd.cells[root].facet_indices.empty())
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

                for (int fi : vd.cells[curr].facet_indices)
                {
                    const int mirror = vd.facets[fi].mirror_facet_index;
                    if (mirror < 0 || mirror >= static_cast<int>(vd.facets.size()))
                        continue;

                    const int neighbor = facetToCell[mirror];
                    if (neighbor < 0 || neighbor == curr)
                        continue;

                    bool opposite = vd.haveOppositeOrientation(
                        vd.facets[fi].vertices_indices,
                        vd.facets[mirror].vertices_indices);

                    if (cellMark[neighbor] == 0)
                    {
                        if (!opposite)
                        {
                            flip_cell(vd, neighbor);
                            opposite = vd.haveOppositeOrientation(
                                vd.facets[fi].vertices_indices,
                                vd.facets[mirror].vertices_indices);
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
        propagate_facets_within_cell(cellIdx, voronoiDiagram);
        audit_cell_edge_orientation(cellIdx, voronoiDiagram);
    }

    propagate_orientation_between_cells(voronoiDiagram, facetToCell);
}

//! @brief Constructs Voronoi cells without using Convex_Hull_3 (in development).
/*!
 * Populates the Voronoi diagram with polyhedral cells derived from the Delaunay
 * triangulation by processing incident edges and cell circulators.
 *
 * @param voronoiDiagram The Voronoi diagram to populate with cells.
 * @param dt The Delaunay triangulation corresponding (dual) to the Voronoi diagram.
 */
void construct_voronoi_cells_from_delaunay_triangulation(VoronoiDiagram &voronoiDiagram, Delaunay &dt)
{
    std::map<std::pair<int, int>, std::vector<int>> edge_to_facets;
    int cellIndex = 0;

    for (Vertex_handle v : dt.finite_vertex_handles())
    {
        if (v->info().is_dummy)
            continue;

        VoronoiCell vc = create_voronoi_cell(v, cellIndex);
        collcet_cell_vertices(dt, v, voronoiDiagram, vc.vertices_indices);
        process_incident_edges(dt, v, voronoiDiagram, vc, edge_to_facets);

        if (vc.facet_indices.size() < 4)
        {
            std::cout << "[WARNING] Cell " << cellIndex << " has only " << vc.facet_indices.size() << " facets, skipping\n";
        }
        else
        {
            voronoiDiagram.cells.push_back(vc);
            v->info().voronoiCellIndex = cellIndex;
            cellIndex++;
        }
    }

    for (const auto &kv : edge_to_facets)
    {
        const std::vector<int> &dfacets = kv.second;
        if (dfacets.size() == 2)
        {
            int f1 = dfacets[0];
            int f2 = dfacets[1];
            voronoiDiagram.facets[f1].mirror_facet_index = f2;
            voronoiDiagram.facets[f2].mirror_facet_index = f1;

            auto &A = voronoiDiagram.facets[f1].vertices_indices;
            auto &B = voronoiDiagram.facets[f2].vertices_indices;
        }
        else if (dfacets.size() == 1)
        {
            int f = dfacets[0];
            voronoiDiagram.facets[f].mirror_facet_index = -1;
        }
    }
}



//! @brief Constructs Voronoi edges from Delaunay facets.
void construct_voronoi_edges(VoronoiDiagram &voronoiDiagram, Delaunay &dt)
{
    voronoiDiagram.edges.clear();
    std::map<std::pair<int, int>, int> segmentMap;              // Maps sorted vertex pairs to edge indices
    std::map<int, std::vector<std::pair<Vector3, int>>> rayMap; // Maps vertex index to (direction, edgeIdx) pairs
    const double EPSILON = 1e-6;

    // not filter by dummy vertices here; downstream selection
    // handles invalid facets during isosurface construction.

    for (Delaunay::Finite_facets_iterator fit = dt.finite_facets_begin(); fit != dt.finite_facets_end(); ++fit)
    {
        Facet facet = *fit;
        CGAL::Object edgeobj = dt.dual(facet);
        Segment3 seg;
        Ray3 ray;

        if (CGAL::assign(seg, edgeobj))
        {
            Cell_handle c1 = facet.first;
            Cell_handle c2 = c1->neighbor(facet.second);
            int idx1 = dt.is_infinite(c1) ? -1 : c1->info().dualVoronoiVertexIndex;
            int idx2 = dt.is_infinite(c2) ? -1 : c2->info().dualVoronoiVertexIndex;

            if (idx1 != -1 && idx2 != -1 && idx1 != idx2)
            { // Finite segment
                int v1 = std::min(idx1, idx2);
                int v2 = std::max(idx1, idx2);
                auto it = segmentMap.find({v1, v2});

                if (it == segmentMap.end())
                {
                    // New segment: create and store with first facet
                    VoronoiEdge vEdge(edgeobj);
                    vEdge.type = 0;
                    vEdge.vertex1 = v1;
                    vEdge.vertex2 = v2;
                    int edgeIdx = voronoiDiagram.edges.size();
                    vEdge.delaunayFacets.push_back(facet);
                    voronoiDiagram.edges.push_back(vEdge);
                    segmentMap[{v1, v2}] = edgeIdx;
                    voronoiDiagram.segmentVertexPairToEdgeIndex[{v1, v2}] = edgeIdx;

                    // Store the edge index in both Delaunay cells sharing this facet
                    int facet1_index = facet.second;
                    // Find the mirror facet index in c2 (the facet of c2 that points back to c1)
                    int facet2_index = -1;
                    for (int i = 0; i < 4; ++i)
                    {
                        if (c2->neighbor(i) == c1)
                        {
                            facet2_index = i;
                            break;
                        }
                    }
                    c1->info().facet_info[facet1_index].dualEdgeIndex = edgeIdx;
                    if (facet2_index >= 0)
                    {
                        c2->info().facet_info[facet2_index].dualEdgeIndex = edgeIdx;
                    }
                }
                else
                {
                    // Segment exists: replace facet if this one has better canonical properties
                    // Use the facet with smallest iFacet value for consistency
                    int edgeIdx = it->second;
                    VoronoiEdge &vEdge = voronoiDiagram.edges[edgeIdx];

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

                    // Also store edge index for this occurrence
                    int facet1_index = facet.second;
                    // Find the mirror facet index in c2
                    int facet2_index = -1;
                    for (int i = 0; i < 4; ++i)
                    {
                        if (c2->neighbor(i) == c1)
                        {
                            facet2_index = i;
                            break;
                        }
                    }
                    c1->info().facet_info[facet1_index].dualEdgeIndex = edgeIdx;
                    if (facet2_index >= 0)
                    {
                        c2->info().facet_info[facet2_index].dualEdgeIndex = edgeIdx;
                    }
                }
            }
        }
        else if (CGAL::assign(ray, edgeobj))
        {
            Cell_handle c1 = facet.first;
            Cell_handle c2 = c1->neighbor(facet.second);
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
                Vector3 dir = ray.direction().vector();
                bool found = false;
                auto it = rayMap.find(vertex1);
                if (it != rayMap.end())
                {
                    for (const auto &pair : it->second)
                    {
                        if (directions_equal(pair.first, dir, EPSILON))
                        {
                            int edgeIdx = pair.second;
                            voronoiDiagram.edges[edgeIdx].delaunayFacets.push_back(facet);

                            // Store edge index in Delaunay facets
                            int facet1_index = facet.second;
                            // Find the mirror facet index in c2
                            int facet2_index = -1;
                            for (int i = 0; i < 4; ++i)
                            {
                                if (c2->neighbor(i) == c1)
                                {
                                    facet2_index = i;
                                    break;
                                }
                            }
                            c1->info().facet_info[facet1_index].dualEdgeIndex = edgeIdx;
                            if (facet2_index >= 0)
                            {
                                c2->info().facet_info[facet2_index].dualEdgeIndex = edgeIdx;
                            }

                            found = true;
                            break;
                        }
                    }
                }
                if (!found)
                {
                    VoronoiEdge vEdge(edgeobj);
                    vEdge.type = 1;
                    vEdge.vertex1 = vertex1;
                    vEdge.vertex2 = -1;
                    vEdge.source = ray.source();
                    vEdge.direction = dir;
                    int edgeIdx = voronoiDiagram.edges.size();
                    vEdge.delaunayFacets.push_back(facet);
                    voronoiDiagram.edges.push_back(vEdge);
                    rayMap[vertex1].push_back({dir, edgeIdx});

                    // Store edge index in Delaunay facets
                    int facet1_index = facet.second;
                    // Find the mirror facet index in c2
                    int facet2_index = -1;
                    for (int i = 0; i < 4; ++i)
                    {
                        if (c2->neighbor(i) == c1)
                        {
                            facet2_index = i;
                            break;
                        }
                    }
                    c1->info().facet_info[facet1_index].dualEdgeIndex = edgeIdx;
                    if (facet2_index >= 0)
                    {
                        c2->info().facet_info[facet2_index].dualEdgeIndex = edgeIdx;
                    }
                }
            }
        }
        // Lines are not expected for finite facets
    }
}

//! @brief Builds Voronoi cell edges for each edge in the diagram.
/*!
 * Creates VoronoiCellEdge entries for cells sharing each edge, collecting cell indices
 * from associated Delaunay facets.
 *
 * @param voronoiDiagram The Voronoi diagram to populate with cell edges.
 * @param dt The Delaunay triangulation.
 */
static void build_cell_edges(
    VoronoiDiagram &voronoiDiagram,
    Delaunay &dt)
{
    const int NUM_DELAUNAY_CELL_VERTICES = 4;
    voronoiDiagram.cellEdges.clear();
    std::vector<std::unordered_set<int>> cellIndicesPerEdge(voronoiDiagram.edges.size());

#pragma omp parallel for
    for (int edgeIdx = 0; edgeIdx < voronoiDiagram.edges.size(); ++edgeIdx)
    {

        const std::vector<Facet> &sharedFacets = voronoiDiagram.edges[edgeIdx].delaunayFacets;

        std::unordered_set<int> cellIndices;

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
                    cellIndices.insert(cellIdx);
                }
            }
        }
        cellIndicesPerEdge[edgeIdx] = std::move(cellIndices);
    }

    // Collect all cell edges in a single pass and store references in Delaunay cells
    for (int edgeIdx = 0; edgeIdx < voronoiDiagram.edges.size(); ++edgeIdx)
    {
        const std::vector<Facet> &sharedFacets = voronoiDiagram.edges[edgeIdx].delaunayFacets;

        for (int cIdx : cellIndicesPerEdge[edgeIdx])
        {
            VoronoiCellEdge cellEdge;
            cellEdge.cellIndex = cIdx;
            cellEdge.edgeIndex = edgeIdx;
            cellEdge.cycleIndices = {};
            cellEdge.nextCellEdge = -1;
            const int cell_edge_index = voronoiDiagram.cellEdges.size();
            voronoiDiagram.cellEdges.push_back(cellEdge);

            // Store reference to cell edge in Delaunay cells
            // For each facet sharing this edge, find which vertex corresponds to this cell
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
                        if (cellIdx == cIdx)
                        {
                            // Found the vertex corresponding to this Voronoi cell
                            // Store the cell edge index
                            c->info().facet_info[facet_index].dualCellEdgeIndex[k - 1] = cell_edge_index;
                        }
                    }
                }
            }
        }
    }
}

//! @brief Populates cell_edge_index array in Delaunay cells.
/*!
 * For each VoronoiCellEdge, finds the corresponding Delaunay cell and facet,
 * then stores the cellEdge index in the Delaunay cell's info structure.
 * This enables O(1) lookup of cell edges during facet construction.
 *
 * @param voronoiDiagram The Voronoi diagram containing cell edges.
 * @param dt The Delaunay triangulation.
 */
void populate_cell_edge_indices(
    VoronoiDiagram &voronoiDiagram,
    Delaunay &dt)
{
    // For each cellEdge, we need to find which Delaunay cell facet it corresponds to
    for (int ceIdx = 0; ceIdx < static_cast<int>(voronoiDiagram.cellEdges.size()); ++ceIdx)
    {
        const VoronoiCellEdge &ce = voronoiDiagram.cellEdges[ceIdx];
        int voronoiCellIdx = ce.cellIndex;
        int voronoiEdgeIdx = ce.edgeIndex;

        // Get the Voronoi edge
        const VoronoiEdge &vEdge = voronoiDiagram.edges[voronoiEdgeIdx];

        // Get the two Voronoi vertices of this edge
        int vv1 = vEdge.vertex1;
        int vv2 = vEdge.vertex2;

        // Skip if either vertex is infinite
        if (vv1 < 0 || vv2 < 0) continue;

        // Voronoi vertices are dual to Delaunay cells
        // Find the Delaunay cell for the current Voronoi cell
        const VoronoiCell &vCell = voronoiDiagram.cells[voronoiCellIdx];
        Vertex_handle delaunay_vertex = vCell.delaunay_vertex;
        if (delaunay_vertex->info().is_dummy) continue;

        // Find incident cells of this Delaunay vertex
        std::vector<Cell_handle> incident_cells;
        dt.incident_cells(delaunay_vertex, std::back_inserter(incident_cells));

        // Find which incident cell has dualVoronoiVertexIndex == vv1 or vv2
        Cell_handle dual_cell1;
        Cell_handle dual_cell2;
        bool found1 = false, found2 = false;

        for (Cell_handle ch : incident_cells)
        {
            if (dt.is_infinite(ch)) continue;
            int dualIdx = ch->info().dualVoronoiVertexIndex;
            if (dualIdx == vv1) {
                dual_cell1 = ch;
                found1 = true;
            }
            if (dualIdx == vv2) {
                dual_cell2 = ch;
                found2 = true;
            }
        }

        // If we found both dual cells, find the facet between them
        if (found1 && found2)
        {
            // Check if they are neighbors
            for (int facet_idx = 0; facet_idx < 4; ++facet_idx)
            {
                Cell_handle neighbor = dual_cell1->neighbor(facet_idx);
                if (neighbor == dual_cell2)
                {
                    // Found the facet! Store the cellEdge index
                    dual_cell1->info().cell_edge_index[facet_idx] = ceIdx;
                    break;
                }
            }

            // Also store in dual_cell2 (the reverse facet)
            for (int facet_idx = 0; facet_idx < 4; ++facet_idx)
            {
                Cell_handle neighbor = dual_cell2->neighbor(facet_idx);
                if (neighbor == dual_cell1)
                {
                    dual_cell2->info().cell_edge_index[facet_idx] = ceIdx;
                    break;
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
    std::unordered_map<int, std::vector<int>> edgeIdx_to_cellEdges;
    for (int ceIdx = 0; ceIdx < (int)voronoiDiagram.cellEdges.size(); ++ceIdx)
    {
        const VoronoiCellEdge &ce = voronoiDiagram.cellEdges[ceIdx];
        edgeIdx_to_cellEdges[ce.edgeIndex].push_back(ceIdx);
    }

    for (auto &kv : edgeIdx_to_cellEdges)
    {
        auto &cellEdgeIndices = kv.second;
        int N = (int)cellEdgeIndices.size();
        for (int i = 0; i < N; i++)
        {
            int ceIdx = cellEdgeIndices[i];
            int nextIdx = cellEdgeIndices[(i + 1) % N];
            voronoiDiagram.cellEdges[ceIdx].nextCellEdge = nextIdx;
        }
    }
}

//! @brief Processes edge mapping for a single Voronoi edge.
/*!
 * Updates the segmentVertexPairToEdgeIndex map for segments, rays, and lines
 * after intersecting with the bounding box.
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
            int v1 = std::min(idx1, idx2);
            int v2 = std::max(idx1, idx2);
            voronoiDiagram.segmentVertexPairToEdgeIndex[{v1, v2}] = edgeIdx;
        }
    }
    // Rays and lines are skipped; no mapping for infinite edges
}

//! @brief Updates edge mappings for all Voronoi edges.
/*!
 * Processes all edges to update segmentVertexPairToEdgeIndex and cellEdgeLookup maps.
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

    voronoiDiagram.cellEdgeLookup.clear();
    for (int ceIdx = 0; ceIdx < static_cast<int>(voronoiDiagram.cellEdges.size()); ++ceIdx)
    {
        const VoronoiCellEdge &ce = voronoiDiagram.cellEdges[ceIdx];
        std::pair<int, int> key = std::make_pair(ce.cellIndex, ce.edgeIndex);
        voronoiDiagram.cellEdgeLookup[key] = ceIdx;
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
    timer.stopTimer("Build cell edges");

    timer.startTimer("Populate cell edge indices", "Construct cell edges");
    populate_cell_edge_indices(voronoiDiagram, dt);
    timer.stopTimer("Populate cell edge indices");

    timer.startTimer("Link cell edges", "Construct cell edges");
    link_cell_edges(voronoiDiagram);
    timer.stopTimer("Link cell edges");

    timer.startTimer("Update edge mapping", "Construct cell edges");
    update_edge_mapping(voronoiDiagram, bbox);
    timer.stopTimer("Update edge mapping");
}


//! @brief Wrap up function of constructing voronoi diagram
void construct_voronoi_diagram(VoronoiDiagram &vd, VDC_PARAM &vdc_param, UnifiedGrid &grid, CGAL::Epick::Iso_cuboid_3 &bbox, Delaunay &dt)
{
    TimingStats& timer = TimingStats::getInstance();

    std::cout << "[INFO] Start constructing Voronoi vertices and edges..." << std::endl;
    timer.startTimer("Construct Voronoi vertices", "4. Voronoi Diagram Construction");
    construct_voronoi_vertices(vd, dt);
    timer.stopTimer("Construct Voronoi vertices");

    timer.startTimer("Construct Voronoi edges", "4. Voronoi Diagram Construction");
    construct_voronoi_edges(vd, dt);
    timer.stopTimer("Construct Voronoi edges");

    timer.startTimer("Compute vertex values", "4. Voronoi Diagram Construction");
    compute_voronoi_values(vd, grid);
    timer.stopTimer("Compute vertex values");

    std::cout << "[INFO] Start constructing Voronoi Cells" << std::endl;
    if (vdc_param.multi_isov)
    {
        if (vdc_param.convex_hull)
        {
            timer.startTimer("Construct Voronoi cells", "4. Voronoi Diagram Construction");
            construct_voronoi_cells_as_convex_hull(vd, dt);
            timer.stopTimer("Construct Voronoi cells");
        }
        else
        {
            timer.startTimer("Construct Voronoi cells", "4. Voronoi Diagram Construction");
            construct_voronoi_cells_from_delaunay_triangulation(vd, dt);
            timer.stopTimer("Construct Voronoi cells");

            timer.startTimer("Validate facet orientations", "4. Voronoi Diagram Construction");
            validate_facet_orientations_and_normals(vd);
            timer.stopTimer("Validate facet orientations");
        }

        timer.startTimer("Construct cell edges", "4. Voronoi Diagram Construction");
        construct_voronoi_cell_edges(vd, bbox, dt);
        timer.stopTimer("Construct cell edges");

        timer.startTimer("Create global facets", "4. Voronoi Diagram Construction");
        vd.create_global_facets();
        timer.stopTimer("Create global facets");
    }

    timer.startTimer("VD initial check", "4. Voronoi Diagram Construction");
    vd.check(false);
    timer.stopTimer("VD initial check");
}
