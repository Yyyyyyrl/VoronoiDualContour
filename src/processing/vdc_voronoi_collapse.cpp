#include "processing/vdc_voronoi.h"
#include "processing/vdc_func.h"
#include "core/vdc_timing.h"


namespace
{

    constexpr double NORMAL_SQ_EPS = 1e-4; // treat facets with smaller squared normal as degenerate

    // --------- Small utilities -------------------------------------------------

    struct DSU
    {
        std::vector<int> p, r; // parent, rank
        explicit DSU(int n = 0) : p(n), r(n, 0) { std::iota(p.begin(), p.end(), 0); }
        void reset(int n)
        {
            p.resize(n);
            r.assign(n, 0);
            std::iota(p.begin(), p.end(), 0);
        }
        int find(int x) { return p[x] == x ? x : p[x] = find(p[x]); }
        void unite(int a, int b)
        {
            a = find(a);
            b = find(b);
            if (a == b)
                return;
            if (r[a] < r[b])
                std::swap(a, b);
            p[b] = a;
            if (r[a] == r[b])
                ++r[a];
        }
    };

    inline bool isFiniteSegmentEdge(const VoronoiEdge &e)
    {
        return e.type == 0 && e.vertex1 >= 0 && e.vertex2 >= 0;
    }

    inline double sqr(double x) { return x * x; }

    inline double squaredDist(const Point &a, const Point &b)
    {
        return sqr(a.x() - b.x()) + sqr(a.y() - b.y()) + sqr(a.z() - b.z());
    }

    // Helper: undirected edge key from two vertex indices
    static inline std::pair<int, int> edge_key(int u, int v)
    {
        return {std::min(u, v), std::max(u, v)};
    }

    // After edge collapse we rebuild cell facets, but the mapping from each
    // facet boundary slot back to its VoronoiCellEdge is lost. The modify-cycles
    // module relies on this to recover per-cell cycle ids when building
    // iso-segments. Restore the association using per-edge ring traversal.
    void rebuild_cell_facet_edge_indices(VoronoiDiagram &vd)
    {
        for (auto &cell : vd.cells)
        {
            const int cellIdx = cell.cellIndex;
            for (int cfIdx : cell.facetIndices)
            {
                if (cfIdx < 0 || cfIdx >= static_cast<int>(vd.cell_facets.size()))
                    continue;
                auto &cf = vd.cell_facets[cfIdx];
                const int n = static_cast<int>(cf.verticesIndices.size());
                cf.cellEdgeIndices.clear();
                cf.cellEdgeIndices.reserve(n);
                for (int i = 0; i < n; ++i)
                {
                    const int a = cf.verticesIndices[i];
                    const int b = cf.verticesIndices[(i + 1) % n];
                    int mapped = -1;
                    if (a >= 0 && b >= 0)
                    {
                        const int globalEdge = vd.findEdgeByVertices(a, b);
                        if (globalEdge >= 0)
                        {
                            mapped = find_cell_edge_for_cell_and_edge(vd, cellIdx, globalEdge);
                        }
                    }
                    cf.cellEdgeIndices.push_back(mapped);
                }
            }
        }
    }


    // Make per-cell facet orientations consistent:
    //  1) Within each cell, ensure two facets sharing an edge traverse that edge in opposite directions.
    //  2) Then, if the majority of non-degenerate facets point inward, flip all facets in that cell.
    void fix_cell_facets_orientation_and_outwardness(VoronoiDiagram &vd)
    {
        for (auto &cell : vd.cells)
        {
            if (cell.facetIndices.empty())
                continue;

            // Build adjacency among facets that share exactly one edge
            const size_t numF = cell.facetIndices.size();
            std::vector<std::map<size_t, std::pair<int, int>>> adj(numF);

            for (size_t i = 0; i < numF; ++i)
            {
                const int f1 = cell.facetIndices[i];
                const auto &verts1 = vd.cell_facets[f1].verticesIndices;
                std::map<std::pair<int, int>, size_t> epos1;
                for (size_t j = 0; j < verts1.size(); ++j)
                {
                    epos1[edge_key(verts1[j], verts1[(j + 1) % verts1.size()])] = j;
                }
                for (size_t k = i + 1; k < numF; ++k)
                {
                    const int f2 = cell.facetIndices[k];
                    const auto &verts2 = vd.cell_facets[f2].verticesIndices;
                    std::pair<int, int> shared = {-1, -1};
                    int cnt = 0;
                    for (size_t j = 0; j < verts2.size(); ++j)
                    {
                        auto key = edge_key(verts2[j], verts2[(j + 1) % verts2.size()]);
                        if (epos1.count(key))
                        {
                            shared = key;
                            if (++cnt > 1)
                                break; // not adjacent if >1 edge shared
                        }
                    }
                    if (cnt == 1)
                    {
                        adj[i][k] = shared;
                        adj[k][i] = shared;
                    }
                }
            }

            // BFS across (possibly multiple) components to ensure opposite directions,
            // and make each component outward by flipping that component if needed.
            std::vector<bool> vis(numF, false);
            const Point site = cell.delaunayVertex->point();
            for (size_t seed = 0; seed < numF; ++seed)
            {
                if (vis[seed])
                    continue;
                // gather component
                std::queue<size_t> q;
                std::vector<size_t> comp;
                q.push(seed);
                vis[seed] = true;
                while (!q.empty())
                {
                    size_t cur = q.front();
                    q.pop();
                    comp.push_back(cur);
                    const int fcur = cell.facetIndices[cur];
                    auto &Vcur = vd.cell_facets[fcur].verticesIndices;
                    for (const auto &kv : adj[cur])
                    {
                        const size_t nb = kv.first;
                        if (!vis[nb])
                        {
                            vis[nb] = true;
                            q.push(nb);
                        }

                        const auto shared = kv.second; // undirected edge
                        const int fnb = cell.facetIndices[nb];
                        auto &Vnb = vd.cell_facets[fnb].verticesIndices;

                        // Determine traversal direction in current facet along shared edge
                        bool cur_uv = false;
                        for (size_t j = 0; j < Vcur.size(); ++j)
                        {
                            int a = Vcur[j];
                            int b = Vcur[(j + 1) % Vcur.size()];
                            if (edge_key(a, b) == shared)
                            {
                                cur_uv = (a == shared.first && b == shared.second);
                                break;
                            }
                        }

                        // Determine traversal in neighbor facet
                        bool nb_uv = false;
                        for (size_t j = 0; j < Vnb.size(); ++j)
                        {
                            int a = Vnb[j];
                            int b = Vnb[(j + 1) % Vnb.size()];
                            if (edge_key(a, b) == shared)
                            {
                                nb_uv = (a == shared.first && b == shared.second);
                                break;
                            }
                        }

                        // If the same direction, reverse neighbor to enforce opposite
                        if (cur_uv == nb_uv)
                        {
                            std::reverse(Vnb.begin(), Vnb.end());
                        }
                    }
                }

                // Choose an anchor facet in this component and flip component if inward
                bool flippedComp = false;
                for (size_t idx : comp)
                {
                    const int fi = cell.facetIndices[idx];
                    const auto &V = vd.cell_facets[fi].verticesIndices;
                    if (V.size() < 3)
                        continue;
                    // centroid
                    Point centroid(0, 0, 0);
                    for (int vid : V)
                        centroid = centroid + (vd.vertices[vid].coord - CGAL::ORIGIN) / V.size();
                    // normal
                    Vector3 normal(0, 0, 0);
                    const size_t n = V.size();
                    for (size_t k = 0; k < n; ++k)
                    {
                        const Point &p1 = vd.vertices[V[k]].coord;
                        const Point &p2 = vd.vertices[V[(k + 1) % n]].coord;
                        normal = normal + Vector3(
                                              (p1.y() - p2.y()) * (p1.z() + p2.z()),
                                              (p1.z() - p2.z()) * (p1.x() + p2.x()),
                                              (p1.x() - p2.x()) * (p1.y() + p2.y()));
                    }
                    normal = normal / 2.0;
                    if (normal.squared_length() <= NORMAL_SQ_EPS)
                        continue;
                    if (CGAL::scalar_product(normal, site - centroid) > 0)
                    {
                        if (debug && cell.cellIndex == 199723)
                        {
                            std::cerr << "[COLLAPSE DBG] Flipping component in cell 199723 (anchor facet " << fi << ") with inward normal.\n";
                        }
                        // flip all facets in this component
                        for (size_t id2 : comp)
                        {
                            const int fj = cell.facetIndices[id2];
                            auto &W = vd.cell_facets[fj].verticesIndices;
                            std::reverse(W.begin(), W.end());
                        }
                    }
                    flippedComp = true;
                    break;
                }
                (void)flippedComp; // may remain false if all facets degenerate
            }
        }
    }

    // Final safety pass: flip any individual facet whose normal is inward
    // relative to its cell. Adjacency consistency can be restored by a
    // subsequent validation pass in the caller.
    void force_outward_per_facet(VoronoiDiagram &vd)
    {
        for (auto &cell : vd.cells)
        {
            const Point site = cell.delaunayVertex->point();
            for (int fi : cell.facetIndices)
            {
                if (fi < 0 || fi >= static_cast<int>(vd.cell_facets.size()))
                    continue;
                auto &V = vd.cell_facets[fi].verticesIndices;
                if (V.size() < 3)
                    continue;

                Point centroid(0, 0, 0);
                for (int idx : V)
                    centroid = centroid + (vd.vertices[idx].coord - CGAL::ORIGIN) / V.size();

                Vector3 normal(0, 0, 0);
                const size_t n = V.size();
                for (size_t k = 0; k < n; ++k)
                {
                    const Point &p1 = vd.vertices[V[k]].coord;
                    const Point &p2 = vd.vertices[V[(k + 1) % n]].coord;
                    normal = normal + Vector3(
                                          (p1.y() - p2.y()) * (p1.z() + p2.z()),
                                          (p1.z() - p2.z()) * (p1.x() + p2.x()),
                                          (p1.x() - p2.x()) * (p1.y() + p2.y()));
                }
                normal = normal / 2.0;
                if (normal.squared_length() <= NORMAL_SQ_EPS)
                    continue;
                const double dot = CGAL::scalar_product(normal, site - centroid);
                if (dot > 0)
                {
                    std::reverse(V.begin(), V.end());
                }
            }
        }
    }

    // Keep representative vertex as the one with smallest original index.  This is
    // stable/deterministic and matches common collapse conventions.
    int chooseRepresentative(const std::vector<int> &group)
    {
        return *std::min_element(group.begin(), group.end());
    }

    // Deduplicate a sequence while keeping the first occurrence order.
    static inline std::vector<int> dedupKeepFirst(const std::vector<int> &seq)
    {
        std::vector<int> out;
        out.reserve(seq.size());
        std::unordered_set<int> seen;
        seen.reserve(seq.size() * 2 + 1);
        for (int v : seq)
            if (!seen.count(v))
            {
                seen.insert(v);
                out.push_back(v);
            }
        return out;
    }

}

// Helper for union-find
static int find(std::vector<int> &mapto, int x)
{
    if (mapto[x] != x)
    {
        mapto[x] = find(mapto, mapto[x]);
    }
    return mapto[x];
}

bool approx_equal_points(const Point &p1, const Point &p2, double eps_sq = 1e-20)
{
    return CGAL::squared_distance(p1, p2) < eps_sq;
}

Vector3 normalize_dir(const Vector3 &v)
{
    double norm = std::sqrt(v.squared_length());
    if (norm < 1e-10)
        return v;
    Vector3 nv = v / norm;
    double components[3] = {CGAL::to_double(nv.x()), CGAL::to_double(nv.y()), CGAL::to_double(nv.z())};
    for (int i = 0; i < 3; ++i)
    {
        if (std::abs(components[i]) > 1e-10)
        {
            if (components[i] < 0)
            {
                return -nv;
            }
            break;
        }
    }
    return nv;
}

bool directions_approx_equal(const Vector3 &d1, const Vector3 &d2, double eps = 1e-10)
{
    Vector3 n1 = normalize_dir(d1);
    Vector3 n2 = normalize_dir(d2);
    double dot = CGAL::scalar_product(n1, n2);
    return std::abs(dot - 1.0) < eps || std::abs(dot + 1.0) < eps; // Parallel or anti-parallel
}

bool lines_approx_equal(const Line3 &l1, const Line3 &l2, double eps_sq = 1e-20)
{
    if (!directions_approx_equal(l1.to_vector(), l2.to_vector()))
        return false;
    // Check if a point on l1 is on l2 (distance to l2 ==0)
    return CGAL::squared_distance(l1.point(), l2) < eps_sq;
}

// Standalone function to collapse small edges
//! @brief Collapses small edges in a Voronoi diagram.
/*!
 * Processes a Voronoi diagram to merge vertices connected by edges shorter than D.
 * Uses union-find data structure to track vertex merges.
 *
 * @param input_vd Input Voronoi diagram
 * @param D Distance threshold for edge collapsing
 * @param bbox Bounding box of the diagram (unused)
 * @return New Voronoi diagram with small edges collapsed
 */
void collapseSmallEdges(const VoronoiDiagram &input_vd,
                        double D,
                        const CGAL::Epick::Iso_cuboid_3 & /*bbox*/,
                        Delaunay & dt,
                        std::vector<int> &out_vertex_mapping,
                        VoronoiDiagram &vd2)
{
    TimingStats& timer = TimingStats::getInstance();

    // 0) Set up (build directly into vd2 to avoid copies)
    vd2.vertices.clear();
    vd2.edges.clear();
    vd2.cells.clear();
    vd2.cell_facets.clear();
    vd2.surface_facets.clear();
    vd2.cellEdges.clear();
    const double D2 = D * D; // squared distances

    const int nV = static_cast<int>(input_vd.vertices.size());
    const int nE = static_cast<int>(input_vd.edges.size());
    const int nC = static_cast<int>(input_vd.cells.size());
    const int nF = static_cast<int>(input_vd.cell_facets.size());
    vd2.vertices.reserve(nV);
    vd2.edges.reserve(nE);
    vd2.cells.reserve(nC);
    vd2.cell_facets.reserve(nF);
    vd2.cellEdges.reserve(input_vd.cellEdges.size());

    // 1) Decide merges: union endpoints of every segment edge shorter than D.
    timer.startTimer("Identify merges (DSU)", "5. Collapse Small Edges");
    DSU dsu(nV);
    for (int ei = 0; ei < nE; ++ei)
    {
        const VoronoiEdge &e = input_vd.edges[ei];
        if (!isFiniteSegmentEdge(e))
            continue;
        const Point &a = input_vd.vertices[e.vertex1].coord;
        const Point &b = input_vd.vertices[e.vertex2].coord;
        if (squaredDist(a, b) < D2)
        {
            dsu.unite(e.vertex1, e.vertex2); // collapse this short edge
        }
    }
    timer.stopTimer("Identify merges (DSU)", "5. Collapse Small Edges");

    // 2) Build groups and pick a representative per merged set.
    timer.startTimer("Build merge groups", "5. Collapse Small Edges");
    std::vector<int> rootToGroup(nV, -1);
    std::vector<std::vector<int>> groups;
    groups.reserve(nV);
    for (int v = 0; v < nV; ++v)
    {
        const int root = dsu.find(v);
        int idx = rootToGroup[root];
        if (idx < 0)
        {
            idx = static_cast<int>(groups.size());
            rootToGroup[root] = idx;
            groups.emplace_back();
        }
        groups[idx].push_back(v);
    }

    // 3) oldV -> newV index mapping (after ordering groups by representative)
    // Determine insertion order for new vertices by ascending representative.
    // Also select a deterministic representative element within each group.
    std::vector<std::pair<int, int>> reps; // (representative, group index)
    reps.reserve(groups.size());
    for (int gi = 0; gi < static_cast<int>(groups.size()); ++gi)
        reps.emplace_back(chooseRepresentative(groups[gi]), gi);
    std::sort(reps.begin(), reps.end(),
              [](const auto &a, const auto &b)
              { return a.first < b.first; });

    std::vector<int> oldToNewV(nV, -1);
    timer.stopTimer("Build merge groups", "5. Collapse Small Edges");

    // 4) Insert merged vertices into `vd2`.
    //    represent a merged vertex by the *representative* original
    //    vertex's coordinate/value, and union of cell membership.
    timer.startTimer("Rebuild vertices", "5. Collapse Small Edges");
    vd2.vertices.reserve(groups.size());
    for (const auto &rep : reps)
    {
        const auto &bucket = groups[rep.second];

        // Representative: smallest original index in the bucket
        const int chosen = rep.first;
        const VoronoiVertex &origVV = input_vd.vertices[chosen];

        // Merge cell memberships (if used elsewhere) and only keep unique set.
        std::vector<int> mergedCells;
        {
            std::unordered_set<int> s;
            for (int vOld : bucket)
            {
                for (int ci : input_vd.vertices[vOld].cellIndices)
                    s.insert(ci);
            }
            mergedCells.reserve(s.size());
            for (int ci : s)
                mergedCells.push_back(ci);
            std::sort(mergedCells.begin(), mergedCells.end());
        }

        const int newIdx = vd2.AddVertex(origVV.coord, origVV.value);
        // Preserve back-references to cells
        vd2.vertices[newIdx].cellIndices = std::move(mergedCells);

        // Map all old vertices in the bucket to this new index
        for (int vOld : bucket)
            oldToNewV[vOld] = newIdx;
    }
    timer.stopTimer("Rebuild vertices", "5. Collapse Small Edges");

    // 5) Rebuild edges — only keep UNcollapsed ones. Also build oldE->newE map.
    //    Preserve and MERGE edge.delaunayFacets across duplicates after collapse.
    timer.startTimer("Rebuild edges", "5. Collapse Small Edges");
    std::vector<int> oldToNewE(nE, -1);

    // Helper to append unique Facets (Delaunay::Facet is usually a pair<Cell_handle,int>)
    auto appendUniqueFacets = [](std::vector<Facet> &dst, const std::vector<Facet> &src)
    {
        for (const auto &f : src)
        {
            if (std::find(dst.begin(), dst.end(), f) == dst.end())
                dst.push_back(f);
        }
    };

    // Local cache for segment edges keyed by normalized (vmin,vmax)
    struct EdgeKeyHash
    {
        size_t operator()(const std::pair<int, int> &p) const noexcept
        {
            return (static_cast<size_t>(static_cast<unsigned int>(p.first)) << 32) ^
                   static_cast<size_t>(static_cast<unsigned int>(p.second));
        }
    };
    std::unordered_map<std::pair<int, int>, int, EdgeKeyHash> localSegMap;
    localSegMap.reserve(static_cast<size_t>(nE) * 2);

    for (int ei = 0; ei < nE; ++ei)
    {
        const VoronoiEdge &e = input_vd.edges[ei];

        if (e.type == 0)
        { // segment
            if (e.vertex1 < 0 || e.vertex2 < 0)
            {
                oldToNewE[ei] = -1;
                continue;
            }
            const int aNew = oldToNewV[e.vertex1];
            const int bNew = oldToNewV[e.vertex2];
            if (aNew == bNew || aNew < 0 || bNew < 0)
            { // collapsed or invalid
                oldToNewE[ei] = -1;
                continue;
            }

            const int vmin = std::min(aNew, bNew);
            const int vmax = std::max(aNew, bNew);
            const auto key = std::make_pair(vmin, vmax);

            // If already have this segment in the new diagram, merge facets
            auto it = localSegMap.find(key);
            if (it != localSegMap.end())
            {
                const int existing = it->second;
                oldToNewE[ei] = existing;
                appendUniqueFacets(vd2.edges[existing].delaunayFacets, e.delaunayFacets);
                continue;
            }

            // Otherwise create it and copy facets
            const Segment3 seg(vd2.vertices[aNew].coord, vd2.vertices[bNew].coord);
            const int ne = vd2.AddSegmentEdge(aNew, bNew, seg);
            oldToNewE[ei] = ne;

            // Preserve original edge's delaunayFacets
            vd2.edges[ne].delaunayFacets.clear();
            appendUniqueFacets(vd2.edges[ne].delaunayFacets, e.delaunayFacets);

            // Maintain the local lookup
            localSegMap[key] = ne;
        }
        else if (e.type == 1)
        { // ray
            Ray3 ray;
            if (!CGAL::assign(ray, e.edgeObject))
                ray = Ray3(e.source, e.direction);
            const int ne = vd2.AddRayEdge(ray);
            oldToNewE[ei] = ne;

            // Preserve facets for rays
            vd2.edges[ne].delaunayFacets = e.delaunayFacets;
        }
        else if (e.type == 2)
        { // line
            Line3 line;
            if (!CGAL::assign(line, e.edgeObject))
                line = Line3(e.source, e.direction);
            const int ne = vd2.AddLineEdge(line);
            oldToNewE[ei] = ne;

            // Preserve facets for lines
            vd2.edges[ne].delaunayFacets = e.delaunayFacets;
        }
        else
        {
            oldToNewE[ei] = -1; // unknown type
        }
    }
    timer.stopTimer("Rebuild edges", "5. Collapse Small Edges");

    // 6) Rebuild Cells & Facets with remapped vertex indices. Drop degenerate
    //    facets that end up with < 3 unique vertices after merging.
    timer.startTimer("Rebuild cells and facets", "5. Collapse Small Edges");
    std::vector<int> oldToNewCell(nC, -1);
    std::vector<int> oldToNewFacet(nF, -1);

    // clone cells in order to keep indices stable
    for (int ci = 0; ci < nC; ++ci)
    {
        const auto &oldCell = input_vd.cells[ci];
        const int nc = vd2.AddCell(oldCell.delaunayVertex);
        oldToNewCell[ci] = nc;

        // Remap the cell’s vertex list ( preserves the order of vertices and also do deduplicate )
        std::vector<int> mappedVerts;
        mappedVerts.reserve(oldCell.verticesIndices.size());
        for (int ov : oldCell.verticesIndices)
        {
            if (ov < 0)
                continue;
            int nv = oldToNewV[ov];
            if (nv >= 0)
                mappedVerts.push_back(nv);
        }
        mappedVerts = dedupKeepFirst(mappedVerts);
        vd2.cells[nc].verticesIndices = std::move(mappedVerts);

        // Copy scalar/iso bookkeeping (if any)
        vd2.cells[nc].isoVertexStartIndex = oldCell.isoVertexStartIndex;
        vd2.cells[nc].numIsoVertices = oldCell.numIsoVertices;
    }

    // Then, rebuild facets in the same order so outside code can keep indices
    for (int fi = 0; fi < nF; ++fi)
    {
        const auto &oldFacet = input_vd.cell_facets[fi];
        std::vector<int> mappedFacetVerts;
        mappedFacetVerts.reserve(oldFacet.verticesIndices.size());
        for (int ov : oldFacet.verticesIndices)
        {
            if (ov < 0)
                continue; // defensive
            const int nv = oldToNewV[ov];
            if (nv >= 0)
                mappedFacetVerts.push_back(nv);
        }
        mappedFacetVerts = dedupKeepFirst(mappedFacetVerts);

        if (mappedFacetVerts.size() < 3)
        {
            // degenerate after collapsing — drop it
            oldToNewFacet[fi] = -1;
            continue;
        }

        const int nf = vd2.AddCellFacet(mappedFacetVerts);
        oldToNewFacet[fi] = nf;

        // Carry auxiliary fields when present
        vd2.cell_facets[nf].orientation = input_vd.cell_facets[fi].orientation;
        vd2.cell_facets[nf].mirror_facet_index = -1;  // will be repaired if needed elsewhere
        vd2.cell_facets[nf].voronoi_facet_index = -1; // re-created later by create_global_facets()
        // Note: cell_edge_indices will be rebuilt by rebuild_cell_facet_edge_indices()
    }

    // Patch each cell’s facet_indices with the new facet ids, skipping dropped
    // ones; preserve order of the remaining facets.
    for (int ci = 0; ci < nC; ++ci)
    {
        const auto &oldCell = input_vd.cells[ci];
        auto &newCell = vd2.cells[oldToNewCell[ci]];
        newCell.facetIndices.clear();
        newCell.facetIndices.reserve(oldCell.facetIndices.size());
        for (int of : oldCell.facetIndices)
        {
            if (of < 0 || of >= nF)
                continue;
            const int nf = oldToNewFacet[of];
            if (nf >= 0)
                newCell.facetIndices.push_back(nf);
        }
    }
    timer.stopTimer("Rebuild cells and facets", "5. Collapse Small Edges");

    // 7) Rebuild VoronoiCellEdges by remapping & filtering
    //     collapsed edges. Also rebuild the nextCellEdge ring per edge.
    timer.startTimer("Rebuild cell edges", "5. Collapse Small Edges");
    std::vector<int> oldToNewCE;
    oldToNewCE.reserve(input_vd.cellEdges.size());
    for (size_t i = 0; i < input_vd.cellEdges.size(); ++i)
        oldToNewCE.push_back(-1);

    // First pass: insert only those cell-edges whose Voronoi edge survived
    for (size_t cei = 0; cei < input_vd.cellEdges.size(); ++cei)
    {
        const auto &oce = input_vd.cellEdges[cei];
        if (oce.edgeIndex < 0 || oce.edgeIndex >= nE)
            continue;
        const int ne = oldToNewE[oce.edgeIndex];
        if (ne < 0)
            continue; // collapsed away -> skip this cell-edge

        const int nc = (oce.cellIndex >= 0 && oce.cellIndex < nC) ? oldToNewCell[oce.cellIndex] : -1;
        if (nc < 0)
            continue;

        VoronoiCellEdge nce;
        nce.cellIndex = nc;
        nce.edgeIndex = ne;
        nce.cycleIndices.clear(); // cycles remap is optional and pipeline-specific
        nce.nextCellEdge = -1;    // will be wired in the second pass

        const int newIdx = static_cast<int>(vd2.cellEdges.size());
        vd2.cellEdges.push_back(nce);
        oldToNewCE[cei] = newIdx;
    }

    // Second pass: rebuild the per-edge ring (nextCellEdge) keeping original
    // relative ordering as much as possible.
    {
        // Group new cell-edge indices by Voronoi edge index
        std::unordered_map<int, std::vector<int>> edgeToCEs;
        edgeToCEs.reserve(vd2.cellEdges.size());
        for (int idx = 0; idx < static_cast<int>(vd2.cellEdges.size()); ++idx)
        {
            edgeToCEs[vd2.cellEdges[idx].edgeIndex].push_back(idx);
        }
        for (auto &kv : edgeToCEs)
        {
            auto &ring = kv.second;
            // Preserve insertion order
            if (ring.size() >= 1)
            {
                const int m = static_cast<int>(ring.size());
                for (int i = 0; i < m; ++i)
                {
                    const int cur = ring[i];
                    const int nxt = ring[(i + 1) % m];
                    vd2.cellEdges[cur].nextCellEdge = nxt;
                }
            }
        }
    }
    timer.stopTimer("Rebuild cell edges", "5. Collapse Small Edges");

    // 8) First ensure every facet is outward relative to its cell, then
    //    enforce edge-consistent orientations within each cell. Finally,
    //    rebuild/refresh global facets & other derived structures.
    timer.startTimer("Fix facet orientations", "5. Collapse Small Edges");
    force_outward_per_facet(vd2);
    fix_cell_facets_orientation_and_outwardness(vd2);
    rebuild_cell_facet_edge_indices(vd2);
    timer.stopTimer("Fix facet orientations", "5. Collapse Small Edges");

    timer.startTimer("Create global facets", "5. Collapse Small Edges");
    vd2.create_global_facets();
    timer.stopTimer("Create global facets", "5. Collapse Small Edges");

    // 9) Copy the vertex index mapping to output parameter
    out_vertex_mapping = oldToNewV;
}
