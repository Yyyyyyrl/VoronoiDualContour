#include "algo/vdc_func.h"
#include <unordered_set>
#include <unordered_map>

// Return undirected edge key for a global facet boundary slot.
static inline std::pair<int, int> facet_slot_edge_key(const VoronoiFacet &gf, int slot)
{
    const int m = (int)gf.vertices_indices.size();
    const int a = gf.vertices_indices[slot];
    const int b = gf.vertices_indices[(slot + 1) % m];
    return (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
}

// Build midpoints for a single cell (bipolar edges only), recording global-edge ids.
// This mirrors the logic used in isosurface construction, trimmed for local recompute.
using EdgeKey = std::pair<int, int>;

struct PairHash
{
    size_t operator()(const EdgeKey &k) const noexcept
    {
        return (static_cast<size_t>(static_cast<uint32_t>(k.first)) << 32) ^ static_cast<uint32_t>(k.second);
    }
};

static void collect_midpoints_for_cell(
    const VoronoiDiagram &vd,
    int cellIndex,
    float isovalue,
    std::vector<MidpointNode> &midpoints,
    std::unordered_map<EdgeKey, int, PairHash> &edge_to_midpoint_index)
{
    const VoronoiCell &vc = vd.cells[cellIndex];
    for (int cfIndex : vc.facet_indices)
    {
        const VoronoiCellFacet &facet = vd.facets[cfIndex];
        const auto &verts = facet.vertices_indices;
        const size_t n = verts.size();
        if (n < 2)
            continue;

        for (size_t j = 0; j < n; ++j)
        {
            size_t idx1 = j;
            size_t idx2 = (j + 1) % n;

            const float val1 = vd.vertices[verts[idx1]].value;
            const float val2 = vd.vertices[verts[idx2]].value;
            if (!is_bipolar(val1, val2, isovalue))
                continue;

            const int vA = verts[idx1];
            const int vB = verts[idx2];
            const Point p1 = vd.vertices[vA].coord;
            const Point p2 = vd.vertices[vB].coord;
            const double t = (isovalue - val1) / (val2 - val1);
            const Point midpoint = p1 + (p2 - p1) * t;

            const auto edge_key = EdgeKey(std::min(vA, vB), std::max(vA, vB));
            if (edge_to_midpoint_index.find(edge_key) == edge_to_midpoint_index.end())
            {
                int globalEdgeIndex = -1;
                auto it = vd.segmentVertexPairToEdgeIndex.find(edge_key);
                if (it != vd.segmentVertexPairToEdgeIndex.end())
                    globalEdgeIndex = it->second;

                MidpointNode node;
                node.point = midpoint;
                node.connected_to.clear();
                node.facet_index = cfIndex;
                node.cycle_index = -1;
                node.global_edge_index = globalEdgeIndex;
                node.is_bipolar = true;

                midpoints.push_back(node);
                edge_to_midpoint_index[edge_key] = (int)midpoints.size() - 1;
            }
        }
    }
}

// Connect midpoints according to global-facet bipolar matches that exist in this cell.
static void connect_midpoints_via_global_matches_for_cell(
    const VoronoiDiagram &vd,
    int cellIndex,
    const std::unordered_map<EdgeKey, int, PairHash> &edge_to_midpoint_index,
    std::vector<MidpointNode> &midpoints)
{
    const VoronoiCell &vc = vd.cells[cellIndex];
    for (int cfIndex : vc.facet_indices)
    {
        const VoronoiCellFacet &cf = vd.facets[cfIndex];
        const int vfi = cf.voronoi_facet_index;
        if (vfi < 0)
            continue;
        const VoronoiFacet &gf = vd.global_facets[vfi];
        for (const auto &pr : gf.bipolar_matches)
        {
            const auto ekA = facet_slot_edge_key(gf, pr.first);
            const auto ekB = facet_slot_edge_key(gf, pr.second);
            auto itA = edge_to_midpoint_index.find(ekA);
            auto itB = edge_to_midpoint_index.find(ekB);
            if (itA == edge_to_midpoint_index.end() || itB == edge_to_midpoint_index.end())
                continue; // this cell doesn't have both midpoints
            const int iA = itA->second, iB = itB->second;
            midpoints[iA].connected_to.push_back(iB);
            midpoints[iB].connected_to.push_back(iA);
        }
    }
}

// Extract cycles from midpoint graph (unique adjacency; expect degree 2 on used nodes).
static void extract_cycles_from_midpoints(
    const std::vector<MidpointNode> &midpoints,
    std::vector<std::vector<int>> &cycles)
{
    const int n = (int)midpoints.size();
    std::vector<char> used(n, 0);

    std::vector<std::vector<int>> adj(n);
    for (int i = 0; i < n; ++i)
    {
        std::unordered_set<int> uniq(midpoints[i].connected_to.begin(), midpoints[i].connected_to.end());
        adj[i].assign(uniq.begin(), uniq.end());
    }

    for (int s = 0; s < n; ++s)
    {
        if (used[s] || adj[s].empty())
            continue;
        int prev = -1, cur = s;
        std::vector<int> cyc;
        while (true)
        {
            used[cur] = 1;
            cyc.push_back(cur);
            if (adj[cur].empty())
                break; // broken chain
            int nxt;
            if (adj[cur].size() == 1)
                nxt = adj[cur][0];
            else
                nxt = (adj[cur][0] == prev) ? adj[cur][1] : adj[cur][0];

            if (nxt == s)
            {
                cycles.push_back(cyc);
                break; // closed
            }
            if (nxt < 0 || nxt >= n || used[nxt])
            { // broken/self-intersecting
                cycles.push_back(cyc);
                break;
            }
            prev = cur;
            cur = nxt;
        }
    }
}

// Map a global-facet boundary slot to the corresponding local slot in a given cell facet.
int map_global_slot_to_cell(const VoronoiFacet &vf,
                            const VoronoiCellFacet &cf,
                            int slot_global)
{
    const int m = (int)vf.vertices_indices.size();
    if (m == 0)
        return -1;
    const int a = vf.vertices_indices[slot_global];
    const int b = vf.vertices_indices[(slot_global + 1) % m];
    const int step = (cf.orientation == 1) ? +1 : -1;
    const auto &C = cf.vertices_indices;
    for (int i = 0; i < (int)C.size(); ++i)
    {
        const int j = (i + step + (int)C.size()) % (int)C.size();
        if (C[i] == a && C[j] == b)
            return i;
    }
    return -1; // not found (shouldn't happen if facet wiring is consistent)
}

// Return the cycle id for the bipolar edge at given local slot in cell facet cf
// (single-slot lookup via cf.cell_edge_indices[local_slot]; NO hashing).
int find_cycle_for_bipolar_edge(const VoronoiDiagram &vd,
                                int cellIndex,
                                int cellFacetIndex,
                                int slot_global)
{
    if (cellFacetIndex < 0 || cellFacetIndex >= (int)vd.facets.size())
        return -1;
    const VoronoiCellFacet &cf = vd.facets[cellFacetIndex];
    const int vfi = cf.voronoi_facet_index;
    if (vfi < 0)
        return -1;
    const VoronoiFacet &vf = vd.global_facets[vfi];

    const int local_slot = map_global_slot_to_cell(vf, cf, slot_global);
    if (local_slot < 0 || local_slot >= (int)cf.cell_edge_indices.size())
        return -1;

    const int ceIdx = cf.cell_edge_indices[local_slot];
    if (ceIdx < 0 || ceIdx >= (int)vd.cellEdges.size())
        return -1;

    const VoronoiCellEdge &ce = vd.cellEdges[ceIdx];
    if (ce.cellIndex != cellIndex)
        return -1;
    if (ce.cycleIndices.empty())
        return -1;
    return ce.cycleIndices.front();
}

// Build vf.iso_segments and fill comp[2] via single-slot lookup per side.
void build_iso_segments_for_facet(VoronoiDiagram &vd,
                                  int vfi,
                                  float /*isovalue*/)
{
    if (vfi < 0 || vfi >= (int)vd.global_facets.size())
        return;
    VoronoiFacet &vf = vd.global_facets[vfi];
    vf.iso_segments.clear();

    if (vf.bipolar_matches.empty())
        return;

    for (const auto &pr : vf.bipolar_matches)
    {
        IsoSegment seg;
        seg.global_facet_index = vfi;
        seg.slotA = pr.first;
        seg.slotB = pr.second;
        seg.edgeA = (seg.slotA >= 0 && seg.slotA < (int)vf.voronoi_edge_indices.size()) ? vf.voronoi_edge_indices[seg.slotA] : -1;
        seg.edgeB = (seg.slotB >= 0 && seg.slotB < (int)vf.voronoi_edge_indices.size()) ? vf.voronoi_edge_indices[seg.slotB] : -1;

        for (int side = 0; side < 2; ++side)
        {
            const int cellIdx = vf.incident_cell_indices[side];
            const int cfIdx = vf.incident_cell_facet_indices[side];
            if (cellIdx < 0 || cfIdx < 0)
            {
                seg.comp[side] = -1;
                continue;
            }
            // Per requirement: use slotA only for component lookup
            seg.comp[side] = find_cycle_for_bipolar_edge(vd, cellIdx, cfIdx, seg.slotA);
        }

        vf.iso_segments.push_back(seg);
    }
}

// Detect whether a facet has duplicate (comp0, comp1) among its iso-segments.
bool facet_has_problematic_iso_segments(const VoronoiDiagram &vd,
                                        int vfi,
                                        std::pair<int, int> *offending_pair)
{
    if (vfi < 0 || vfi >= (int)vd.global_facets.size())
        return false;
    const auto &vf = vd.global_facets[vfi];
    const int n = (int)vf.iso_segments.size();
    for (int i = 0; i < n; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {
            const auto &a = vf.iso_segments[i];
            const auto &b = vf.iso_segments[j];
            if (a.comp[0] >= 0 && a.comp[1] >= 0 &&
                a.comp[0] == b.comp[0] && a.comp[1] == b.comp[1])
            {
                if (offending_pair)
                    *offending_pair = {i, j};
                return true;
            }
        }
    }
    return false;
}

// Flip matching method for a facet (SEP_POS <-> SEP_NEG).
void flip_bipolar_match_method(VoronoiFacet &vf)
{
    switch (vf.bipolar_match_method)
    {
    case BIPOLAR_MATCH_METHOD::SEP_POS:
        vf.bipolar_match_method = BIPOLAR_MATCH_METHOD::SEP_NEG;
        break;
    case BIPOLAR_MATCH_METHOD::SEP_NEG:
        vf.bipolar_match_method = BIPOLAR_MATCH_METHOD::SEP_POS;
        break;
    default:
        vf.bipolar_match_method = BIPOLAR_MATCH_METHOD::SEP_NEG;
        break;
    }
}

// Recompute cycles for a single cell after facet-local rematching.
void recompute_cell_cycles_for_matches_single_cell(VoronoiDiagram &vd,
                                                   int cellIndex,
                                                   float isovalue)
{
    if (cellIndex < 0 || cellIndex >= (int)vd.cells.size())
        return;

    // 1) Clear per-cell cycle tags on cellEdges and reset cell cycles
    for (auto &ce : vd.cellEdges)
    {
        if (ce.cellIndex == cellIndex)
            ce.cycleIndices.clear();
    }
    vd.cells[cellIndex].cycles.clear();

    // 2) Collect midpoints (per-cell) and connect using current global matches
    std::vector<MidpointNode> midpoints;
    std::unordered_map<EdgeKey, int, PairHash> edge_to_midpoint_index;
    collect_midpoints_for_cell(vd, cellIndex, isovalue, midpoints, edge_to_midpoint_index);
    connect_midpoints_via_global_matches_for_cell(vd, cellIndex, edge_to_midpoint_index, midpoints);

    // 3) Extract cycles
    std::vector<std::vector<int>> cycles_idx;
    extract_cycles_from_midpoints(midpoints, cycles_idx);

    // 4) Write cycles back to vd (cycle list per cell, and cycleIndices per cell-edge)
    VoronoiCell &vc = vd.cells[cellIndex];
    int cycId = 0;
    for (const auto &single_cycle : cycles_idx)
    {
        Cycle cycle;
        cycle.voronoi_cell_index = cellIndex;
        cycle.midpoint_indices = single_cycle;
        for (size_t i = 0; i < single_cycle.size(); ++i)
        {
            int a = single_cycle[i];
            int b = single_cycle[(i + 1) % single_cycle.size()];
            cycle.edges.emplace_back(a, b);
        }
        cycle.compute_centroid(midpoints);

        for (int mpIdx : single_cycle)
        {
            // Tag midpoint with its cycle id
            const int globalEdgeIdx = midpoints[mpIdx].global_edge_index;
            if (globalEdgeIdx < 0)
                continue;
            const auto it = vd.cellEdgeLookup.find({cellIndex, globalEdgeIdx});
            if (it == vd.cellEdgeLookup.end())
                continue;
            auto &cyclesVec = vd.cellEdges[it->second].cycleIndices;
            if (std::find(cyclesVec.begin(), cyclesVec.end(), cycId) == cyclesVec.end())
                cyclesVec.push_back(cycId);
        }

        vc.cycles.push_back(std::move(cycle));
        ++cycId;
    }
}

// Populate incident cell indices and corresponding cell-facet indices for each global facet.
void populate_incident_cells_for_global_facets(VoronoiDiagram &vd)
{
    // Reset arrays
    for (auto &gf : vd.global_facets)
    {
        gf.incident_cell_indices = {-1, -1};
        gf.incident_cell_facet_indices = {-1, -1};
    }

    // Build mapping by scanning cells and their facets
    for (const auto &cell : vd.cells)
    {
        for (int cfIndex : cell.facet_indices)
        {
            if (cfIndex < 0 || cfIndex >= (int)vd.facets.size())
                continue;
            const int vfi = vd.facets[cfIndex].voronoi_facet_index;
            if (vfi < 0 || vfi >= (int)vd.global_facets.size())
                continue;
            auto &gf = vd.global_facets[vfi];
            if (gf.incident_cell_indices[0] == -1)
            {
                gf.incident_cell_indices[0] = cell.cellIndex;
                gf.incident_cell_facet_indices[0] = cfIndex;
            }
            else if (gf.incident_cell_indices[1] == -1 && gf.incident_cell_indices[0] != cell.cellIndex)
            {
                gf.incident_cell_indices[1] = cell.cellIndex;
                gf.incident_cell_facet_indices[1] = cfIndex;
            }
            // else: already has two (or duplicate), ignore
        }
    }
}

// Full “modify cycles” pass over all global facets at this isovalue.
void modify_cycles_pass(VoronoiDiagram &vd, float isovalue)
{
    // Ensure global facets know their incident cells
    populate_incident_cells_for_global_facets(vd);

    vd.compute_bipolar_matches(isovalue);

    // First pass: detect problematic facets and flip locally, collecting affected cells
    std::unordered_set<int> dirty_cells;
    const int nGF = (int)vd.global_facets.size();
    for (int vfi = 0; vfi < nGF; ++vfi)
    {
        build_iso_segments_for_facet(vd, vfi, isovalue);
        if (vd.global_facets[vfi].iso_segments.empty())
            continue;

        if (facet_has_problematic_iso_segments(vd, vfi, nullptr))
        {
            // Flip method for this facet and recompute ONLY this facet’s matches
            flip_bipolar_match_method(vd.global_facets[vfi]);
            recompute_bipolar_matches_for_facet(vd, vfi, isovalue);

            // Mark both incident cells for cycle recomputation
            const auto &gf = vd.global_facets[vfi];
            for (int side = 0; side < 2; ++side)
            {
                const int cidx = gf.incident_cell_indices[side];
                if (cidx >= 0)
                    dirty_cells.insert(cidx);
            }
        }
    }

    // Second pass: recompute cycles once per affected cell
    for (int cidx : dirty_cells)
    {
        recompute_cell_cycles_for_matches_single_cell(vd, cidx, isovalue);
    }
}
