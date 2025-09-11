#include "voronoi/vdc_voronoi.h"
#include "algo/vdc_func.h"
#include <iostream>
#include <fstream>
#include <iomanip>

// Build a minimal synthetic scenario with two cells sharing one global quad facet
// Vertices carry alternating values around the facet to produce 4 bipolar edges.
// We pre-seed cycleIndices to force an initial "problematic" detection, then run
// modify_cycles_pass and check that it flips the facet method and re-tags cycles.

static int add_vertex(VoronoiDiagram &vd, const Point &p, float val)
{
    int idx = (int)vd.vertices.size();
    VoronoiVertex vv(p);
    vv.index = idx;
    vv.value = val;
    vd.vertices.push_back(vv);
    return idx;
}

static int add_cell_edge(VoronoiDiagram &vd, int cellIndex, int globalEdge)
{
    VoronoiCellEdge ce;
    ce.cellIndex = cellIndex;
    ce.edgeIndex = globalEdge;
    ce.nextCellEdge = -1;
    int idx = (int)vd.cellEdges.size();
    vd.cellEdges.push_back(ce);
    vd.cellEdgeLookup[{cellIndex, globalEdge}] = idx;
    return idx;
}

static void build_quad_case(VoronoiDiagram &vd, float isovalue, int &vfi)
{
    vd.vertices.clear();
    vd.edges.clear();
    vd.cells.clear();
    vd.facets.clear();
    vd.global_facets.clear();
    vd.cellEdges.clear();
    vd.cellEdgeLookup.clear();
    vd.segmentVertexPairToEdgeIndex.clear();

    // 4 vertices in a unit square (z=0)
    int v0 = add_vertex(vd, Point(0, 0, 0), -1.0f);
    int v1 = add_vertex(vd, Point(1, 0, 0), +1.0f);
    int v2 = add_vertex(vd, Point(1, 1, 0), -1.0f);
    int v3 = add_vertex(vd, Point(0, 1, 0), +1.0f);

    // Global facet boundary as quad
    VoronoiFacet gf;
    gf.index = 0;
    gf.vertices_indices = {v0, v1, v2, v3};
    // Create actual segment edges in vd for each boundary slot
    auto add_seg = [&](int a, int b) {
        Segment3 seg(vd.vertices[a].coord, vd.vertices[b].coord);
        return vd.AddSegmentEdge(a, b, seg);
    };
    gf.voronoi_edge_indices.clear();
    gf.voronoi_edge_indices.push_back(add_seg(v0, v1));
    gf.voronoi_edge_indices.push_back(add_seg(v1, v2));
    gf.voronoi_edge_indices.push_back(add_seg(v2, v3));
    gf.voronoi_edge_indices.push_back(add_seg(v3, v0));
    gf.bipolar_match_method = BIPOLAR_MATCH_METHOD::SEP_POS; // initial
    vd.global_facets.push_back(gf);
    vfi = 0;

    // segmentVertexPairToEdgeIndex populated by AddSegmentEdge

    // Two cells that share this facet
    VoronoiCell c0(nullptr); c0.cellIndex = 0;
    VoronoiCell c1(nullptr); c1.cellIndex = 1;
    vd.cells.push_back(c0);
    vd.cells.push_back(c1);

    // Two cell facets referencing the same global facet with opposite orientation
    VoronoiCellFacet cf0; cf0.vertices_indices = {v0, v1, v2, v3}; cf0.voronoi_facet_index = vfi; cf0.orientation = 1;
    VoronoiCellFacet cf1; cf1.vertices_indices = {v0, v3, v2, v1}; cf1.voronoi_facet_index = vfi; cf1.orientation = -1;
    int cfi0 = (int)vd.facets.size(); vd.facets.push_back(cf0);
    int cfi1 = (int)vd.facets.size(); vd.facets.push_back(cf1);
    vd.cells[0].facet_indices.push_back(cfi0);
    vd.cells[1].facet_indices.push_back(cfi1);

    // Wire incident cells back to global facet
    vd.global_facets[vfi].incident_cell_indices = {0, 1};
    vd.global_facets[vfi].incident_cell_facet_indices = {cfi0, cfi1};

    // Create per-cell "cellEdges" for each boundary slot and map them into cell_facet slots
    // Slots: 0:(v0,v1), 1:(v1,v2), 2:(v2,v3), 3:(v3,v0)
    // Cell 0
    std::vector<int> ceIdx0(4);
    for (int k = 0; k < 4; ++k)
        ceIdx0[k] = add_cell_edge(vd, 0, vd.global_facets[vfi].voronoi_edge_indices[k]);
    vd.facets[cfi0].cell_edge_indices = ceIdx0;

    // Cell 1: same edges
    std::vector<int> ceIdx1(4);
    for (int k = 0; k < 4; ++k)
        ceIdx1[k] = add_cell_edge(vd, 1, vd.global_facets[vfi].voronoi_edge_indices[k]);
    vd.facets[cfi1].cell_edge_indices = ceIdx1;

    // Seed initial per-cell cycles so build_iso_segments can classify comps before we run the pass
    // Force slots 0 and 2 to map to cycle 0 in both cells (so pairs anchored at 0 and 2 collide)
    vd.cellEdges[ceIdx0[0]].cycleIndices = {0};
    vd.cellEdges[ceIdx0[2]].cycleIndices = {0};
    vd.cellEdges[ceIdx0[1]].cycleIndices = {1};
    vd.cellEdges[ceIdx0[3]].cycleIndices = {1};
    vd.cellEdges[ceIdx1[0]].cycleIndices = {0};
    vd.cellEdges[ceIdx1[2]].cycleIndices = {0};
    vd.cellEdges[ceIdx1[1]].cycleIndices = {1};
    vd.cellEdges[ceIdx1[3]].cycleIndices = {1};

    // Initial matches: pair slots (0,1) and (2,3)
    vd.global_facets[vfi].bipolar_matches = {{0,1},{2,3}};
}

static void print_iso_segments(const VoronoiDiagram &vd, int vfi)
{
    const auto &vf = vd.global_facets[vfi];
    std::cout << "iso-segments for gf#" << vfi << ":\n";
    int idx = 0;
    for (auto const &g : vf.iso_segments)
    {
        std::cout << "  seg#" << idx++ << " slotA=" << g.slotA << " slotB=" << g.slotB
                  << " comp0=" << g.comp[0] << " comp1=" << g.comp[1] << "\n";
    }
}

// Build a second case: hexagonal facet shared by cells 1 and 2, with six bipolar edges.
static void build_hex_case(VoronoiDiagram &vd, float isovalue, int &vfi_hex)
{
    // Hexagon centered near z=1, unit-ish size
    int base_vert_idx = (int)vd.vertices.size();
    auto addV = [&](double x, double y, double z, float val) {
        return add_vertex(vd, Point(x, y, z), val);
    };
    // Alternate values - + - + - +
    int h0 = addV(0.0, 0.0, 1.0, -1.0f);
    int h1 = addV(1.0, 0.0, 1.0, +1.0f);
    int h2 = addV(1.5, 0.866, 1.0, -1.0f);
    int h3 = addV(1.0, 1.732, 1.0, +1.0f);
    int h4 = addV(0.0, 1.732, 1.0, -1.0f);
    int h5 = addV(-0.5, 0.866, 1.0, +1.0f);

    VoronoiFacet gf;
    gf.index = (int)vd.global_facets.size();
    gf.vertices_indices = {h0, h1, h2, h3, h4, h5};
    // Create actual segment edges
    auto add_seg = [&](int a, int b) {
        Segment3 seg(vd.vertices[a].coord, vd.vertices[b].coord);
        return vd.AddSegmentEdge(a, b, seg);
    };
    gf.voronoi_edge_indices.clear();
    gf.voronoi_edge_indices.push_back(add_seg(h0, h1));
    gf.voronoi_edge_indices.push_back(add_seg(h1, h2));
    gf.voronoi_edge_indices.push_back(add_seg(h2, h3));
    gf.voronoi_edge_indices.push_back(add_seg(h3, h4));
    gf.voronoi_edge_indices.push_back(add_seg(h4, h5));
    gf.voronoi_edge_indices.push_back(add_seg(h5, h0));
    gf.bipolar_match_method = BIPOLAR_MATCH_METHOD::SEP_POS;
    vd.global_facets.push_back(gf);
    vfi_hex = gf.index;

    // segmentVertexPairToEdgeIndex populated by AddSegmentEdge

    // Ensure cell 1 exists from quad case; add cell 2
    int c1 = (int)vd.cells.size() > 1 ? 1 : (int)vd.cells.size()-1;
    VoronoiCell c2(nullptr); c2.cellIndex = (int)vd.cells.size();
    vd.cells.push_back(c2);

    // Two cell facets for gf: one in cell 1, one in cell 2 (opposite orientation)
    VoronoiCellFacet cf1; cf1.vertices_indices = {h0, h1, h2, h3, h4, h5}; cf1.voronoi_facet_index = vfi_hex; cf1.orientation = 1;
    VoronoiCellFacet cf2; cf2.vertices_indices = {h0, h5, h4, h3, h2, h1}; cf2.voronoi_facet_index = vfi_hex; cf2.orientation = -1;
    int cfi1 = (int)vd.facets.size(); vd.facets.push_back(cf1);
    int cfi2 = (int)vd.facets.size(); vd.facets.push_back(cf2);
    vd.cells[c1].facet_indices.push_back(cfi1);
    vd.cells[c2.cellIndex].facet_indices.push_back(cfi2);

    // Wire incident cells to global facet
    vd.global_facets[vfi_hex].incident_cell_indices = {c1, c2.cellIndex};
    vd.global_facets[vfi_hex].incident_cell_facet_indices = {cfi1, cfi2};

    // Create cellEdges and wire into cell_facet slots
    std::vector<int> ceIdx1(6), ceIdx2(6);
    for (int k = 0; k < 6; ++k)
        ceIdx1[k] = add_cell_edge(vd, c1, vd.global_facets[vfi_hex].voronoi_edge_indices[k]);
    for (int k = 0; k < 6; ++k)
        ceIdx2[k] = add_cell_edge(vd, c2.cellIndex, vd.global_facets[vfi_hex].voronoi_edge_indices[k]);
    vd.facets[cfi1].cell_edge_indices = ceIdx1;
    vd.facets[cfi2].cell_edge_indices = ceIdx2;

    // Seed cycles: even slots → cycle 0, odd slots → cycle 1 (in both cells)
    for (int k = 0; k < 6; ++k)
    {
        vd.cellEdges[ceIdx1[k]].cycleIndices = {(k % 2)};
        vd.cellEdges[ceIdx2[k]].cycleIndices = {(k % 2)};
    }

    // Initial matches: consecutive pairs (0,1), (2,3), (4,5)
    vd.global_facets[vfi_hex].bipolar_matches = {{0,1},{2,3},{4,5}};
}

// Dump a compact JSON to visualize facet boundaries and iso-segments before/after
static void write_case_json(const std::string &path,
                            const VoronoiDiagram &vd,
                            float isovalue,
                            const std::map<int, std::vector<std::pair<int,int>>> &matches_before)
{
    auto slot_midpoint = [&](const VoronoiFacet &gf, int slot) -> Point {
        const int m = (int)gf.vertices_indices.size();
        const int a = gf.vertices_indices[slot];
        const int b = gf.vertices_indices[(slot + 1) % m];
        const auto &pa = vd.vertices[a].coord;
        const auto &pb = vd.vertices[b].coord;
        const float va = vd.vertices[a].value;
        const float vb = vd.vertices[b].value;
        const double t = (isovalue - va) / (vb - va);
        return pa + (pb - pa) * t;
    };

    // Helper: compute cycles per cell given a map of facet->matches
    auto compute_cycles_for = [&](const std::map<int, std::vector<std::pair<int,int>>> &matches_map)
        -> std::vector<std::pair<int, std::vector<std::vector<Point>>>>
    {
        using EdgeKey = std::pair<int,int>;
        struct PairHash { size_t operator()(const EdgeKey &k) const noexcept { return (size_t(uint32_t(k.first))<<32) ^ uint32_t(k.second); } };
        std::vector<std::pair<int, std::vector<std::vector<Point>>>> out;

        for (const auto &cell : vd.cells)
        {
            // Build midpoint nodes per cell from all its incident global facets
            std::unordered_map<EdgeKey, int, PairHash> nodeIndex;
            std::vector<Point> nodePoints;
            std::vector<std::vector<int>> adj; // adjacency by node index

            auto ensure_node_for_slot = [&](const VoronoiFacet &gf, int slot) -> int
            {
                const int m = (int)gf.vertices_indices.size();
                int a = gf.vertices_indices[slot];
                int b = gf.vertices_indices[(slot + 1) % m];
                if (a > b) std::swap(a,b);
                EdgeKey key(a,b);
                auto it = nodeIndex.find(key);
                if (it != nodeIndex.end()) return it->second;
                // Only add if the edge is bipolar for this isovalue
                float va = vd.vertices[a].value;
                float vb = vd.vertices[b].value;
                if (!is_bipolar(va, vb, isovalue)) return -1;
                const Point &pa = vd.vertices[a].coord;
                const Point &pb = vd.vertices[b].coord;
                double t = (isovalue - va) / (vb - va);
                Point pm = pa + (pb - pa) * t;
                int idx = (int)nodePoints.size();
                nodePoints.push_back(pm);
                adj.emplace_back();
                nodeIndex.emplace(key, idx);
                return idx;
            };

            // Create nodes for all bipolar edges and connect pairs per facet matches
            for (int cfIndex : cell.facet_indices)
            {
                if (cfIndex < 0 || cfIndex >= (int)vd.facets.size()) continue;
                int vfi = vd.facets[cfIndex].voronoi_facet_index;
                if (vfi < 0 || vfi >= (int)vd.global_facets.size()) continue;
                const VoronoiFacet &gf = vd.global_facets[vfi];

                // ensure nodes for all slots that are bipolar
                const int m = (int)gf.vertices_indices.size();
                for (int s = 0; s < m; ++s) { (void)ensure_node_for_slot(gf, s); }

                // pick matches from map (fallback to current gf matches if absent)
                auto itM = matches_map.find(vfi);
                const std::vector<std::pair<int,int>> &pairs = (itM != matches_map.end()) ? itM->second : gf.bipolar_matches;
                for (const auto &pr : pairs)
                {
                    int iA = ensure_node_for_slot(gf, pr.first);
                    int iB = ensure_node_for_slot(gf, pr.second);
                    if (iA < 0 || iB < 0) continue;
                    adj[iA].push_back(iB);
                    adj[iB].push_back(iA);
                }
            }

            // Extract cycles
            const int n = (int)nodePoints.size();
            std::vector<char> used(n, 0);
            // dedup adjacency
            for (int i = 0; i < n; ++i)
            {
                std::unordered_set<int> u(adj[i].begin(), adj[i].end());
                adj[i].assign(u.begin(), u.end());
            }

            std::vector<std::vector<Point>> cyclesPts;
            for (int s = 0; s < n; ++s)
            {
                if (used[s] || adj[s].empty()) continue;
                int prev = -1, cur = s;
                std::vector<Point> cyc;
                while (true)
                {
                    used[cur] = 1;
                    cyc.push_back(nodePoints[cur]);
                    if (adj[cur].empty()) break;
                    int nxt = (adj[cur].size() == 1) ? adj[cur][0] : (adj[cur][0] == prev ? adj[cur][1] : adj[cur][0]);
                    if (nxt == s)
                    {
                        // close the loop by repeating start point for drawing convenience
                        cyc.push_back(nodePoints[s]);
                        cyclesPts.push_back(std::move(cyc));
                        break;
                    }
                    if (nxt < 0 || nxt >= n || used[nxt])
                    {
                        cyclesPts.push_back(std::move(cyc));
                        break;
                    }
                    prev = cur; cur = nxt;
                }
            }
            out.emplace_back(cell.cellIndex, std::move(cyclesPts));
        }
        return out;
    };

    // Prepare after-map from current vd state
    std::map<int, std::vector<std::pair<int,int>>> matches_after;
    for (const auto &gf : vd.global_facets)
        matches_after[gf.index] = gf.bipolar_matches;

    auto cycles_before = compute_cycles_for(matches_before);
    auto cycles_after  = compute_cycles_for(matches_after);

    std::ofstream os(path);
    os << std::fixed << std::setprecision(6);
    os << "{\n";
    os << "  \"isovalue\": " << isovalue << ",\n";
    // vertices
    os << "  \"vertices\": [\n";
    for (size_t i = 0; i < vd.vertices.size(); ++i)
    {
        const auto &v = vd.vertices[i];
        os << "    {\"id\":" << i
           << ",\"x\":" << v.coord.x()
           << ",\"y\":" << v.coord.y()
           << ",\"z\":" << v.coord.z()
           << ",\"val\":" << v.value << "}";
        if (i + 1 < vd.vertices.size()) os << ",";
        os << "\n";
    }
    os << "  ],\n";

    // facets with before/after iso-segments
    os << "  \"facets\": [\n";
    for (size_t fi = 0; fi < vd.global_facets.size(); ++fi)
    {
        const auto &gf = vd.global_facets[fi];
        os << "    {\"vfi\": " << gf.index << ", \"verts\": [";
        for (size_t k = 0; k < gf.vertices_indices.size(); ++k)
        {
            os << gf.vertices_indices[k];
            if (k + 1 < gf.vertices_indices.size()) os << ",";
        }
        os << "], ";

        // before matches
        auto itB = matches_before.find((int)fi);
        os << "\"matches_before\": [";
        if (itB != matches_before.end())
        {
            for (size_t j = 0; j < itB->second.size(); ++j)
            {
                const auto &p = itB->second[j];
                os << "[" << p.first << "," << p.second << "]";
                if (j + 1 < itB->second.size()) os << ",";
            }
        }
        os << "], ";

        // after matches (current in vd)
        os << "\"matches_after\": [";
        for (size_t j = 0; j < gf.bipolar_matches.size(); ++j)
        {
            const auto &p = gf.bipolar_matches[j];
            os << "[" << p.first << "," << p.second << "]";
            if (j + 1 < gf.bipolar_matches.size()) os << ",";
        }
        os << "], ";

        // computed segment endpoints for before/after (for convenience)
        os << "\"segments_before\": [";
        if (itB != matches_before.end())
        {
            const auto &pairs = itB->second;
            for (size_t j = 0; j < pairs.size(); ++j)
            {
                Point a = slot_midpoint(gf, pairs[j].first);
                Point b = slot_midpoint(gf, pairs[j].second);
                os << "[[" << a.x() << "," << a.y() << "," << a.z() << "],["
                   << b.x() << "," << b.y() << "," << b.z() << "]]";
                if (j + 1 < pairs.size()) os << ",";
            }
        }
        os << "], ";

        os << "\"segments_after\": [";
        for (size_t j = 0; j < gf.bipolar_matches.size(); ++j)
        {
            Point a = slot_midpoint(gf, gf.bipolar_matches[j].first);
            Point b = slot_midpoint(gf, gf.bipolar_matches[j].second);
            os << "[[" << a.x() << "," << a.y() << "," << a.z() << "],["
               << b.x() << "," << b.y() << "," << b.z() << "]]";
            if (j + 1 < gf.bipolar_matches.size()) os << ",";
        }
        os << "]}";
        if (fi + 1 < vd.global_facets.size()) os << ",";
        os << "\n";
    }
    os << "  ],\n";

    // Build simple 3D prisms for each global facet to visualize the two incident cells
    os << "  \"facet_prisms\": [\n";
    for (size_t fi = 0; fi < vd.global_facets.size(); ++fi)
    {
        const auto &gf = vd.global_facets[fi];
        // Gather facet points
        std::vector<Point> P; P.reserve(gf.vertices_indices.size());
        for (int idx : gf.vertices_indices) P.push_back(vd.vertices[idx].coord);
        // Compute normal from first triangle
        if (P.size() < 3) continue;
        Vector3 n = CGAL::cross_product(P[1] - P[0], P[2] - P[0]);
        double len = std::sqrt(n.squared_length());
        if (len < 1e-12) n = Vector3(0,0,1); else n = n / len;
        double h = 0.05; // half-thickness
        // verts for side 0 and side 1
        std::vector<Point> Vp; Vp.reserve(P.size());
        std::vector<Point> Vm; Vm.reserve(P.size());
        for (const auto &q : P)
        {
            Vp.push_back(q + n * h);
            Vm.push_back(q - n * h);
        }
        int cid0 = gf.incident_cell_indices[0];
        int cid1 = gf.incident_cell_indices[1];
        auto dump_ring = [&](const std::vector<Point> &R) {
            os << "[";
            for (size_t k = 0; k < R.size(); ++k)
            {
                os << "[" << R[k].x() << "," << R[k].y() << "," << R[k].z() << "]";
                if (k + 1 < R.size()) os << ",";
            }
            os << "]";
        };
        os << "    {\"vfi\": " << gf.index << ", \"cid0\": " << cid0 << ", \"cid1\": " << cid1 << ", \"top\": ";
        dump_ring(Vp);
        os << ", \"bottom\": ";
        dump_ring(Vm);
        os << "}";
        if (fi + 1 < vd.global_facets.size()) os << ",";
        os << "\n";
    }
    os << "  ]\n";
    // cycles per cell (before/after)
    os << ",\n  \"cell_cycles_before\": [\n";
    for (size_t i = 0; i < cycles_before.size(); ++i)
    {
        const auto &pr = cycles_before[i];
        os << "    {\"cid\": " << pr.first << ", \"cycles\": [";
        for (size_t j = 0; j < pr.second.size(); ++j)
        {
            const auto &cy = pr.second[j];
            os << "[";
            for (size_t k = 0; k < cy.size(); ++k)
            {
                os << "[" << cy[k].x() << "," << cy[k].y() << "," << cy[k].z() << "]";
                if (k + 1 < cy.size()) os << ",";
            }
            os << "]";
            if (j + 1 < pr.second.size()) os << ",";
        }
        os << "]}";
        if (i + 1 < cycles_before.size()) os << ",";
        os << "\n";
    }
    os << "  ],\n";

    os << "  \"cell_cycles_after\": [\n";
    for (size_t i = 0; i < cycles_after.size(); ++i)
    {
        const auto &pr = cycles_after[i];
        os << "    {\"cid\": " << pr.first << ", \"cycles\": [";
        for (size_t j = 0; j < pr.second.size(); ++j)
        {
            const auto &cy = pr.second[j];
            os << "[";
            for (size_t k = 0; k < cy.size(); ++k)
            {
                os << "[" << cy[k].x() << "," << cy[k].y() << "," << cy[k].z() << "]";
                if (k + 1 < cy.size()) os << ",";
            }
            os << "]";
            if (j + 1 < pr.second.size()) os << ",";
        }
        os << "]}";
        if (i + 1 < cycles_after.size()) os << ",";
        os << "\n";
    }
    os << "  ]\n";
    os << "}\n";
}

int main()
{
    VoronoiDiagram vd;
    float iso = 0.0f;
    int vfi = -1, vfi_hex = -1;

    build_quad_case(vd, iso, vfi);
    build_hex_case(vd, iso, vfi_hex);

    // Capture matches before
    std::map<int, std::vector<std::pair<int,int>>> matches_before;
    for (size_t i = 0; i < vd.global_facets.size(); ++i)
        matches_before[(int)i] = vd.global_facets[i].bipolar_matches;

    // Build iso-segments from seeded cycles and check detection on both facets
    for (int id : {vfi, vfi_hex})
    {
        build_iso_segments_for_facet(vd, id, iso);
        print_iso_segments(vd, id);
        bool problematic = facet_has_problematic_iso_segments(vd, id, nullptr);
        std::cout << "gf#" << id << " initially problematic: " << (problematic ? "yes" : "no") << "\n";
    }

    // Run the modify-cycles pass (flips method if problematic, recomputes cycles where needed)
    modify_cycles_pass(vd, iso);

    // Rebuild iso-segments after the pass and print
    for (int id : {vfi, vfi_hex})
    {
        build_iso_segments_for_facet(vd, id, iso);
        print_iso_segments(vd, id);
        std::cout << "Facet method after pass (gf#" << id << "): "
                  << matchMethodToString(vd.global_facets[id].bipolar_match_method) << "\n";
    }

    // Check that per-cell cycleIndices exist for bipolar edges (not empty)
    int num_nonempty = 0;
    for (auto const &ce : vd.cellEdges)
        if (!ce.cycleIndices.empty())
            ++num_nonempty;
    std::cout << "cellEdges with cycles after pass: " << num_nonempty << "/" << vd.cellEdges.size() << "\n";

    // Write visualization JSON (contains vertices, facet polygons, and iso-segments before/after)
    write_case_json("modcyc_case.json", vd, iso, matches_before);
    std::cout << "Wrote modcyc_case.json (for tools/plot_modcyc.py).\n";

    std::cout << "test_modcyc completed." << std::endl;
    return 0;
}
