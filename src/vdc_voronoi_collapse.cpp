#include "vdc_voronoi.h"


namespace
{

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
VoronoiDiagram collapseSmallEdges(const VoronoiDiagram &input_vd,
                                  double D,
                                  const CGAL::Epick::Iso_cuboid_3 & /*bbox*/,
                                  Delaunay & /*dt*/)
{
    // 0) Set up
    VoronoiDiagram out;      // fresh diagram — don't mutate input_vd
    const double D2 = D * D; // squared distances

    const int nV = static_cast<int>(input_vd.vertices.size());
    const int nE = static_cast<int>(input_vd.edges.size());
    const int nC = static_cast<int>(input_vd.cells.size());
    const int nF = static_cast<int>(input_vd.facets.size());

    // 1) Decide merges: union endpoints of every segment edge shorter than D.
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

    // 2) Build groups and pick a representative per merged set.
    std::unordered_map<int, std::vector<int>> groups;
    groups.reserve(nV);
    for (int v = 0; v < nV; ++v)
        groups[dsu.find(v)].push_back(v);

    // 3) oldV -> representative oldV (root); then oldV -> newV index mapping
    std::vector<int> oldRoot(nV);
    for (int v = 0; v < nV; ++v)
        oldRoot[v] = dsu.find(v);

    // Determine insertion order for new vertices by ascending representative.
    // Also select a deterministic representative element within each group.
    std::vector<int> reps;
    reps.reserve(groups.size());
    for (auto &kv : groups)
        reps.push_back(chooseRepresentative(kv.second));
    std::sort(reps.begin(), reps.end());

    std::vector<int> oldToNewV(nV, -1);

    // 4) Insert merged vertices into `out`.
    //    represent a merged vertex by the *representative* original
    //    vertex’s coordinate/value, and union of cell membership.
    for (int rep : reps)
    {
        const auto &bucket = groups[dsu.find(rep)];

        // Representative: smallest original index in the bucket
        const int chosen = chooseRepresentative(bucket);
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

        const int newIdx = out.AddVertex(origVV.coord, origVV.value);
        // Preserve back-references to cells
        out.vertices[newIdx].cellIndices = std::move(mergedCells);

        // Map all old vertices in the bucket to this new index
        for (int vOld : bucket)
            oldToNewV[vOld] = newIdx;
    }

    // 5) Rebuild edges — only keep UNcollapsed ones. Also build oldE->newE map.
    //    Preserve and MERGE edge.delaunayFacets across duplicates after collapse.
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
    std::map<std::pair<int, int>, int> localSegMap;

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
                appendUniqueFacets(out.edges[existing].delaunayFacets, e.delaunayFacets);
                continue;
            }

            // Otherwise create it and copy facets
            const Segment3 seg(out.vertices[aNew].coord, out.vertices[bNew].coord);
            const int ne = out.AddSegmentEdge(aNew, bNew, seg);
            oldToNewE[ei] = ne;

            // Preserve original edge's delaunayFacets
            out.edges[ne].delaunayFacets.clear();
            appendUniqueFacets(out.edges[ne].delaunayFacets, e.delaunayFacets);

            // Maintain both the local and public lookups
            localSegMap[key] = ne;
            out.segmentVertexPairToEdgeIndex[key] = ne;
        }
        else if (e.type == 1)
        { // ray
            Ray3 ray;
            if (!CGAL::assign(ray, e.edgeObject))
                ray = Ray3(e.source, e.direction);
            const int ne = out.AddRayEdge(ray);
            oldToNewE[ei] = ne;

            // Preserve facets for rays
            out.edges[ne].delaunayFacets = e.delaunayFacets;
        }
        else if (e.type == 2)
        { // line
            Line3 line;
            if (!CGAL::assign(line, e.edgeObject))
                line = Line3(e.source, e.direction);
            const int ne = out.AddLineEdge(line);
            oldToNewE[ei] = ne;

            // Preserve facets for lines
            out.edges[ne].delaunayFacets = e.delaunayFacets;
        }
        else
        {
            oldToNewE[ei] = -1; // unknown type
        }
    }

    // 6) Rebuild Cells & Facets with remapped vertex indices. Drop degenerate
    //    facets that end up with < 3 unique vertices after merging.
    std::vector<int> oldToNewCell(nC, -1);
    std::vector<int> oldToNewFacet(nF, -1);

    // clone cells in order to keep indices stable
    for (int ci = 0; ci < nC; ++ci)
    {
        const auto &oldCell = input_vd.cells[ci];
        const int nc = out.AddCell(oldCell.delaunay_vertex);
        oldToNewCell[ci] = nc;

        // Remap the cell’s vertex list ( preserves the order of vertices and also do deduplicate )
        std::vector<int> mappedVerts;
        mappedVerts.reserve(oldCell.vertices_indices.size());
        for (int ov : oldCell.vertices_indices)
        {
            if (ov < 0)
                continue;
            int nv = oldToNewV[ov];
            if (nv >= 0)
                mappedVerts.push_back(nv);
        }
        mappedVerts = dedupKeepFirst(mappedVerts);
        out.cells[nc].vertices_indices = std::move(mappedVerts);

        // Copy scalar/iso bookkeeping (if any)
        out.cells[nc].isoVertexStartIndex = oldCell.isoVertexStartIndex;
        out.cells[nc].numIsoVertices = oldCell.numIsoVertices;
    }

    // Then, rebuild facets in the same order so outside code can keep indices
    for (int fi = 0; fi < nF; ++fi)
    {
        const auto &oldFacet = input_vd.facets[fi];
        std::vector<int> mappedFacetVerts;
        mappedFacetVerts.reserve(oldFacet.vertices_indices.size());
        for (int ov : oldFacet.vertices_indices)
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

        const int nf = out.AddCellFacet(mappedFacetVerts);
        oldToNewFacet[fi] = nf;

        // Carry auxiliary fields when present
        out.facets[nf].orientation = input_vd.facets[fi].orientation;
        out.facets[nf].mirror_facet_index = -1;  // will be repaired if needed elsewhere
        out.facets[nf].voronoi_facet_index = -1; // re-created later by create_global_facets()
    }

    // Patch each cell’s facet_indices with the new facet ids, skipping dropped
    // ones; preserve order of the remaining facets.
    for (int ci = 0; ci < nC; ++ci)
    {
        const auto &oldCell = input_vd.cells[ci];
        auto &newCell = out.cells[oldToNewCell[ci]];
        newCell.facet_indices.clear();
        newCell.facet_indices.reserve(oldCell.facet_indices.size());
        for (int of : oldCell.facet_indices)
        {
            if (of < 0 || of >= nF)
                continue;
            const int nf = oldToNewFacet[of];
            if (nf >= 0)
                newCell.facet_indices.push_back(nf);
        }
    }

    // 7) Rebuild VoronoiCellEdges and cellEdgeLookup by remapping & filtering
    //     collapsed edges. Also rebuild the nextCellEdge ring per edge.
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

        const int newIdx = static_cast<int>(out.cellEdges.size());
        out.cellEdges.push_back(nce);
        out.cellEdgeLookup[{nc, ne}] = newIdx;
        oldToNewCE[cei] = newIdx;
    }

    // Second pass: rebuild the per-edge ring (nextCellEdge) keeping original
    // relative ordering as much as possible.
    {
        // Group new cell-edge indices by Voronoi edge index
        std::unordered_map<int, std::vector<int>> edgeToCEs;
        edgeToCEs.reserve(out.cellEdges.size());
        for (int idx = 0; idx < static_cast<int>(out.cellEdges.size()); ++idx)
        {
            edgeToCEs[out.cellEdges[idx].edgeIndex].push_back(idx);
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
                    out.cellEdges[cur].nextCellEdge = nxt;
                }
            }
        }
    }

    // 8) Rebuild/refresh global facets & other derived structures.
    out.create_global_facets();

    // 9) Sanity checks
    out.check(false);
    out.checkAdvanced(false);

    return out;
}