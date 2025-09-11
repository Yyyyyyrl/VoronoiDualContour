#include "vdc_func.h"


//! @brief Constructs Voronoi vertices for the given voronoi Diagram instance.
void construct_voronoi_vertices(VoronoiDiagram &voronoiDiagram, Delaunay &dt)
{
    int vertexIndex = 0;
    voronoiDiagram.vertices.reserve(dt.number_of_finite_cells());
    for (Delaunay::Finite_cells_iterator cit = dt.finite_cells_begin(); cit != dt.finite_cells_end(); ++cit)
    {
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
        uniqueVertexIndices.insert(vertex_index);
    }
    vertices_indices.assign(uniqueVertexIndices.begin(), uniqueVertexIndices.end());
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
    std::vector<int> facetVertexIndices;
    int finite_cell_count = 0;
    bool involves_infinite = false;
    do
    {
        if (dt.is_infinite(cc))
        {
            involves_infinite = true;
        }
        else
        {
            int vertex_index = cc->info().dualVoronoiVertexIndex;
            facetVertexIndices.push_back(vertex_index);
            finite_cell_count++;
        }
        ++cc;
    } while (cc != start);

    std::set<int> unique_vertices(facetVertexIndices.begin(), facetVertexIndices.end());

    if (unique_vertices.size() >= 3)
    {
        // Clean duplicates
        std::vector<int> cleaned;
        for (size_t k = 0; k < facetVertexIndices.size(); ++k)
        {
            if (k == 0 || facetVertexIndices[k] != facetVertexIndices[k - 1])
            {
                cleaned.push_back(facetVertexIndices[k]);
            }
        }
        if (cleaned.size() > 1 && cleaned.front() == cleaned.back())
        {
            cleaned.pop_back();
        }
        if (cleaned.size() < 3)
        {
            std::cout << "[DEBUG] Degenerate after cleaning: " << cleaned.size() << " verts\n";
            return VoronoiCellFacet();
        }
        std::vector<int> orderedFacetVertices = std::move(cleaned);

        // Determine orientation of the facet using CGAL::orientation ( taking determinant )
        Point P1 = voronoiDiagram.vertices[orderedFacetVertices[0]].coord;
        Point P2 = voronoiDiagram.vertices[orderedFacetVertices[1]].coord;
        Point P3 = voronoiDiagram.vertices[orderedFacetVertices[2]].coord;
        Point site = delaunay_vertex->point();
        CGAL::Orientation orient = CGAL::orientation(P1, P2, P3, site);
        if (orient == CGAL::POSITIVE)
        {
            std::reverse(orderedFacetVertices.begin(), orderedFacetVertices.end());
        }
        else if (orient == CGAL::ZERO)
        {
        }

        VoronoiCellFacet facet;
        facet.vertices_indices = orderedFacetVertices;

        int facetIndex = voronoiDiagram.facets.size();
        voronoiDiagram.facets.push_back(facet);
        facet_indices.push_back(facetIndex);
        std::pair<int, int> edge_key = std::make_pair(
            std::min(v1->info().index, v2->info().index),
            std::max(v1->info().index, v2->info().index));
        edge_to_facets[edge_key].push_back(facetIndex);

        // Map vertex loop → cell-edge loop for this cell facet
        facet.cell_edge_indices.clear();
        const int n = (int)facet.vertices_indices.size();

        for (int i = 0; i < n; ++i)
        {
            const int a = facet.vertices_indices[i];
            const int b = facet.vertices_indices[(i + 1) % n];
            const int vmin = std::min(a, b), vmax = std::max(a, b);

            // 1) find global segment edge
            const auto eit = voronoiDiagram.segmentVertexPairToEdgeIndex.find({vmin, vmax});
            if (eit == voronoiDiagram.segmentVertexPairToEdgeIndex.end())
            {
                facet.cell_edge_indices.push_back(-1);
                continue;
            }
            const int globalEdge = eit->second;

            // 2) find the cell-edge for (this cell, globalEdge)
            const int thisCellIndex = vcIdx;
            const auto ceit = voronoiDiagram.cellEdgeLookup.find({thisCellIndex, globalEdge});
            facet.cell_edge_indices.push_back(ceit == voronoiDiagram.cellEdgeLookup.end() ? -1 : ceit->second);
        }
        return facet;
    }
    else
    {
        std::cout << "[DEBUG] Degenerate facet for edge with " << finite_cell_count << " finite cells\n";
        return VoronoiCellFacet();
    }
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
        if (facet.vertices_indices.empty())
        {
            std::cout << "[WARNING] Facet construction failed for edge with " << finite_cell_count << " finite cells\n";
            continue;
        }

        // Verify facet validity
        if (facet.vertices_indices.size() < 3)
        {
            std::cout << "[ERROR] Degenerate interior facet detected with " << facet.vertices_indices.size() << " vertices despite " << finite_cell_count << " finite cells\n";
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

//! @brief Validates that each Facet in the Voronoi Diagram has the correct
/// orientation and normal vector.
/*!
 *
 */
void validate_facet_orientations_and_normals(VoronoiDiagram &voronoiDiagram)
{
    // propagation per cell
    for (auto &cell : voronoiDiagram.cells)
    {
        if (cell.facet_indices.empty())
            continue;

        // Build facet adjacency: vector of maps: index in facet_indices -> {adj_index in facet_indices: shared_edge_key}
        size_t num_facets = cell.facet_indices.size();
        std::vector<std::map<size_t, std::pair<int, int>>> facet_adj(num_facets);
        for (size_t i = 0; i < num_facets; ++i)
        {
            int f1 = cell.facet_indices[i];
            const auto &verts1 = voronoiDiagram.facets[f1].vertices_indices;
            std::map<std::pair<int, int>, size_t> edge_to_pos1;
            for (size_t j = 0; j < verts1.size(); ++j)
            {
                int u = verts1[j];
                int v = verts1[(j + 1) % verts1.size()];
                edge_to_pos1[get_edge_key(u, v)] = j;
            }

            for (size_t k = i + 1; k < num_facets; ++k)
            {
                int f2 = cell.facet_indices[k];
                const auto &verts2 = voronoiDiagram.facets[f2].vertices_indices;
                std::pair<int, int> shared_edge = {-1, -1};
                int shared_count = 0;
                for (size_t j = 0; j < verts2.size(); ++j)
                {
                    int u = verts2[j];
                    int v = verts2[(j + 1) % verts2.size()];
                    auto key = get_edge_key(u, v);
                    if (edge_to_pos1.count(key))
                    {
                        shared_edge = key;
                        shared_count++;
                        if (shared_count > 1)
                            break; // Not adjacent if >1 edge shared
                    }
                }
                if (shared_count == 1)
                {
                    facet_adj[i][k] = shared_edge;
                    facet_adj[k][i] = shared_edge;
                }
            }
        }

        // BFS to propagate
        std::vector<bool> visited(num_facets, false);
        std::queue<size_t> q;
        q.push(0);
        visited[0] = true;

        while (!q.empty())
        {
            size_t curr = q.front();
            q.pop();
            int curr_f = cell.facet_indices[curr];
            auto &curr_verts = voronoiDiagram.facets[curr_f].vertices_indices;

            for (const auto &kv : facet_adj[curr])
            {
                size_t adj = kv.first;
                if (visited[adj])
                    continue;
                visited[adj] = true;
                q.push(adj);

                auto shared_edge = kv.second;
                int adj_f = cell.facet_indices[adj];
                auto &adj_verts = voronoiDiagram.facets[adj_f].vertices_indices;

                // Find direction in curr: true if u to v (min to max)
                bool curr_dir_uv = false;
                for (size_t j = 0; j < curr_verts.size(); ++j)
                {
                    int a = curr_verts[j];
                    int b = curr_verts[(j + 1) % curr_verts.size()];
                    if (get_edge_key(a, b) == shared_edge)
                    {
                        curr_dir_uv = (a == shared_edge.first && b == shared_edge.second);
                        break;
                    }
                }

                // In adj
                bool adj_dir_uv = false;
                for (size_t j = 0; j < adj_verts.size(); ++j)
                {
                    int a = adj_verts[j];
                    int b = adj_verts[(j + 1) % adj_verts.size()];
                    if (get_edge_key(a, b) == shared_edge)
                    {
                        adj_dir_uv = (a == shared_edge.first && b == shared_edge.second);
                        break;
                    }
                }

                // If same direction, reverse adj
                if (curr_dir_uv == adj_dir_uv)
                {
                    std::reverse(adj_verts.begin(), adj_verts.end());
                    // std::cout << "[INFO] Reversed intra-cell facet " << adj_f << " in cell " << cell.cellIndex << " to match opposite edge {" << shared_edge.first << "," << shared_edge.second << "} with facet " << curr_f << "\n";
                }
            }
        }

        int outward_count = 0;
        int total_non_deg = 0;
        for (size_t f = 0; f < num_facets; ++f)
        {
            int fi = cell.facet_indices[f];
            auto &V = voronoiDiagram.facets[fi].vertices_indices;
            if (V.size() < 3)
                continue;

            Point centroid(0, 0, 0);
            for (int idx : V)
            {
                centroid = centroid + (voronoiDiagram.vertices[idx].coord - CGAL::ORIGIN) / V.size();
            }

            Vector3 normal(0, 0, 0);
            size_t n = V.size();
            for (size_t k = 0; k < n; ++k)
            {
                const Point &p1 = voronoiDiagram.vertices[V[k]].coord;
                const Point &p2 = voronoiDiagram.vertices[V[(k + 1) % n]].coord;
                normal = normal + Vector3(
                                      (p1.y() - p2.y()) * (p1.z() + p2.z()),
                                      (p1.z() - p2.z()) * (p1.x() + p2.x()),
                                      (p1.x() - p2.x()) * (p1.y() + p2.y()));
            }
            normal = normal / 2.0;

            double sq_norm = normal.squared_length();
            if (sq_norm > 1e-10)
            { // Non-degenerate
                Vector3 v = cell.delaunay_vertex->point() - centroid;
                double dot = CGAL::scalar_product(normal, v);
                if (dot <= 0)
                    outward_count++;
                total_non_deg++;
            }
        }

        if (total_non_deg > 0 && outward_count < total_non_deg / 2)
        {
            // Majority inward, reverse all facets in the cell
            // std::cout << "[INFO] Reversing all facets in cell " << cell.cellIndex << " to make majority outward (outward_count: " << outward_count << " / " << total_non_deg << ")\n";
            for (size_t f = 0; f < num_facets; ++f)
            {
                int fi = cell.facet_indices[f];
                auto &V = voronoiDiagram.facets[fi].vertices_indices;
                std::reverse(V.begin(), V.end());
            }
        }
    }
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

                VoronoiEdge vEdge(edgeobj);
                vEdge.type = 0;
                vEdge.vertex1 = v1;
                vEdge.vertex2 = v2;
                int edgeIdx = voronoiDiagram.edges.size();
                vEdge.delaunayFacets.push_back(facet);
                voronoiDiagram.edges.push_back(vEdge);
                segmentMap[{v1, v2}] = edgeIdx;
                voronoiDiagram.segmentVertexPairToEdgeIndex[{v1, v2}] = edgeIdx;
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
                            voronoiDiagram.edges[pair.second].delaunayFacets.push_back(facet);
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

    // Collect all cell edges in a single pass
    for (int edgeIdx = 0; edgeIdx < voronoiDiagram.edges.size(); ++edgeIdx)
    {
        for (int cIdx : cellIndicesPerEdge[edgeIdx])
        {
            VoronoiCellEdge cellEdge;
            cellEdge.cellIndex = cIdx;
            cellEdge.edgeIndex = edgeIdx;
            cellEdge.cycleIndices = {};
            cellEdge.nextCellEdge = -1;
            voronoiDiagram.cellEdges.push_back(cellEdge);
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
    voronoiDiagram.cellEdges.clear();

    std::clock_t start = std::clock();

    build_cell_edges(voronoiDiagram, dt);
    std::clock_t check1 = std::clock();
    double duration1 = static_cast<double>(check1 - start) / CLOCKS_PER_SEC;

    std::cout << "build cell edge Execution time: " << duration1 << " seconds" << std::endl;

    link_cell_edges(voronoiDiagram);
    std::clock_t check2 = std::clock();
    double duration2 = static_cast<double>(check2 - check1) / CLOCKS_PER_SEC;

    std::cout << "link cell edge Execution time: " << duration2 << " seconds" << std::endl;

    update_edge_mapping(voronoiDiagram, bbox);

    std::clock_t check3 = std::clock();

    double duration3 = static_cast<double>(check3 - check2) / CLOCKS_PER_SEC;

    std::cout << "update edge mapping Execution time: " << duration3 << " seconds" << std::endl;
}


//! @brief Wrap up function of constructing voronoi diagram
void construct_voronoi_diagram(VoronoiDiagram &vd, VDC_PARAM &vdc_param, UnifiedGrid &grid, CGAL::Epick::Iso_cuboid_3 &bbox, Delaunay &dt)
{
    std::clock_t start = std::clock();
    std::cout << "[INFO] Start constructing Voronoi vertices and edges..." << std::endl;
    construct_voronoi_vertices(vd, dt);
    std::clock_t check0 = std::clock();
    double duration0 = static_cast<double>(check0 - start) / CLOCKS_PER_SEC;
    std::cout << "construct vertices Execution time: " << duration0 << " seconds" << std::endl;

    construct_voronoi_edges(vd, dt);
    std::clock_t check1 = std::clock();
    double duration1 = static_cast<double>(check1 - start) / CLOCKS_PER_SEC;
    std::cout << "construct edges Execution time: " << duration1 << " seconds" << std::endl;

    compute_voronoi_values(vd, grid);
    std::clock_t check2 = std::clock();
    double duration2 = static_cast<double>(check2 - check1) / CLOCKS_PER_SEC;

    if (vdc_param.multi_isov)
    {
        if (vdc_param.convex_hull)
        {
            construct_voronoi_cells_as_convex_hull(vd, dt);
        }
        else
        {
            construct_voronoi_cells_from_delaunay_triangulation(vd, dt);
            std::clock_t check3 = std::clock();
            double duration3 = static_cast<double>(check3 - check2) / CLOCKS_PER_SEC;
            std::cout << "construct cells Execution time: " << duration3 << " seconds" << std::endl;
            validate_facet_orientations_and_normals(vd);
            std::clock_t check4 = std::clock();
            double duration4 = static_cast<double>(check4 - check3) / CLOCKS_PER_SEC;
            std::cout << "validate facet orientations and normals in cells Execution time: " << duration4 << " seconds" << std::endl;
        }

        std::clock_t t = std::clock();

        construct_voronoi_cell_edges(vd, bbox, dt);
        std::clock_t check5 = std::clock();
        double duration5 = static_cast<double>(check5 - t) / CLOCKS_PER_SEC;
        std::cout << "construct cell edges Execution time: " << duration5 << " seconds" << std::endl;

        vd.create_global_facets();
    }

    std::clock_t t2 = std::clock();
    vd.check(false);
    std::clock_t check6 = std::clock();
    double duration6 = static_cast<double>(check6 - t2) / CLOCKS_PER_SEC;
    std::cout << "vd.check() Execution time: " << duration6 << " seconds" << std::endl;
}
