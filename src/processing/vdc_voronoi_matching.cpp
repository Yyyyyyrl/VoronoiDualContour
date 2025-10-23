#include "processing/vdc_voronoi.h"

static std::vector<int> collectFacetVoronoiEdges(const VoronoiDiagram &vd, const std::vector<int> &verts)
{
    std::vector<int> out;
    const int n = (int)verts.size();
    out.reserve(n);
    for (int i = 0; i < n; ++i)
    {
        int a = verts[i], b = verts[(i + 1) % n];
        int vmin = std::min(a, b), vmax = std::max(a, b);
        auto it = vd.segmentVertexPairToEdgeIndex.find({vmin, vmax});
        out.push_back(it == vd.segmentVertexPairToEdgeIndex.end() ? -1 : it->second);
    }
    return out;
}

void VoronoiDiagram::create_global_facets()
{
    // Group cell-facets by canonical full-vertex key (order-invariant, no collisions).
    std::map<std::vector<int>, std::vector<int>> keyToCellFacets;

    for (size_t fi = 0; fi < facets.size(); ++fi)
    {
        const auto &F = facets[fi].vertices_indices;
        auto key = getFacetHashKey(F);
        keyToCellFacets[key].push_back(static_cast<int>(fi));
    }

    // Build unique (global) facets from the groups
    for (const auto &kv : keyToCellFacets)
    {
        const auto &fvec = kv.second;
        if (fvec.size() > 2)
        {
            // With a correct key, this should never happen for a proper 3D Voronoi complex.
            // If it does, the input geometry/topology is degenerate or inconsistent.
            std::ostringstream oss;
            oss << "Facet shared by more than 2 cells. Key has " << fvec.size() << " facets; "
                << "vertices: ";
            for (int v : kv.first)
                oss << v << " ";
            throw std::runtime_error(oss.str());
        }

        VoronoiFacet vf;
        vf.index = static_cast<int>(global_facets.size());

        // Choose a primary representative facet from this group
        const int primary = fvec[0];
        vf.vertices_indices = facets[primary].vertices_indices;

        // Boundary Voronoi edges for this polygon (k -> edge(v[k], v[k+1]))
        vf.voronoi_edge_indices = collectFacetVoronoiEdges(*this, vf.vertices_indices);
        if (vf.voronoi_edge_indices.size() != vf.vertices_indices.size())
        {
            // Defensive: keep sizes consistent; mark missing with -1
            vf.voronoi_edge_indices.resize(vf.vertices_indices.size(), -1);
        }

        vf.primary_cell_facet_index = primary;
        vf.bipolar_match_method = BIPOLAR_MATCH_METHOD::SEP_POS; // default; can be changed later
        global_facets.push_back(vf);

        // Wire the primary back to this global facet
        facets[primary].voronoi_facet_index = vf.index;
        facets[primary].orientation = 1;

        // If there is exactly one mirror, it must have opposite winding
        if (fvec.size() == 2)
        {
            const int secondary = fvec[1];
            const bool opposite =
                haveOppositeOrientation(facets[primary].vertices_indices,
                                        facets[secondary].vertices_indices);
            if (!opposite)
            {
                throw std::runtime_error(
                    "Paired facets " + std::to_string(primary) + " and " +
                    std::to_string(secondary) + " do not have opposite orientations.");
            }
            facets[secondary].voronoi_facet_index = vf.index;
            facets[secondary].orientation = -1;
        }
    }
}

// --- helper: classify and pair facet bipolar edges -------------------------

// Pairs bipolar edges inside a VoronoiFacet according to the method:
//  - SEP_NEG: start from first (+,−) edge, then pair (start,next), (next2,next3), ...
//  - SEP_POS: start from first (−,+) edge, then pair likewise
//  - UNCONSTRAINED_MATCH: pair consecutive in the bipolar list
//
// Output:
//  - vf.bipolar_edge_indices: indices k (edges) around the facet boundary that are bipolar
//  - vf.bipolar_matches: pairs of indices into vf.bipolar_edge_indices (NOT global edge ids)
static void match_facet_bipolar_edges(const VoronoiDiagram &vd,
                                      VoronoiFacet &vf,
                                      float isovalue)
{
    vf.bipolar_edge_indices.clear();
    vf.bipolar_matches.clear();

    const int m = (int)vf.vertices_indices.size();
    if (m < 2)
        return;

    // Classify edges in the facet’s OWN boundary order
    // edgeType[k] = +1 for (− → +), -1 for (+ → −), 0 otherwise
    std::vector<int> edgeType(m, 0);
    for (int k = 0; k < m; ++k)
    {
        int i0 = vf.vertices_indices[k];
        int i1 = vf.vertices_indices[(k + 1) % m];

        // skip if not a finite segment edge in global map
        if (k < (int)vf.voronoi_edge_indices.size())
        {
            int ei = vf.voronoi_edge_indices[k];
            if (ei < 0 || ei >= (int)vd.edges.size() || vd.edges[ei].type != 0)
                continue;
        }

        float v0 = vd.vertices[i0].value;
        float v1 = vd.vertices[i1].value;
        if (is_bipolar(v0, v1, isovalue))
        {
            edgeType[k] = (v0 < v1) ? +1 : -1;
            vf.bipolar_edge_indices.push_back(k); // store boundary slot k
        }
    }

    // Debug: print the slot signs
    if (debug)
    {
        std::ostringstream oss;
        oss << "[ISO-MATCH] gf=" << vf.index << " slots:";
        for (int k = 0; k < m; ++k)
        {
            char c = (edgeType[k] > 0) ? '+' : (edgeType[k] < 0 ? '-' : '.');
            oss << ' ' << c;
        }
        std::cerr << oss.str() << "\n";
    }

    const int nB = (int)vf.bipolar_edge_indices.size();
    if (nB < 2)
    {
        if (debug)
        {
            std::cerr << "[ISO-MATCH] gf=" << vf.index << " nBip=" << nB << " → matches (none)\n";
        }
        return;
    }

    // ensure nB is even
    if (nB % 2)
    {
        std::cerr << "[warn] facet " << vf.index << " has odd # bipolar edges at iso =" << isovalue << " \n";
        std::cerr << " voronoi edges in facet: \n";
        for (int ei : vf.voronoi_edge_indices)
        {
            if (ei < 0 || ei >= static_cast<int>(vd.edges.size()))
            {
                std::cerr << "Invalid edge index: " << ei << "\n";
                continue;
            }

            const VoronoiEdge &edge = vd.edges[ei];
            Segment3 segment;
            Line3 line;
            Ray3 ray;

            if (CGAL::assign(segment, edge.edgeObject))
            {
                std::cerr << "Segment(" << segment.source() << " - " << segment.target() << ")\n";
            }
            else if (CGAL::assign(line, edge.edgeObject))
            {
                std::cerr << "Line(" << line.point(0) << " - " << line.point(1) << ")\n";
            }
            else if (CGAL::assign(ray, edge.edgeObject))
            {
                std::cerr << "Ray(" << ray.source() << ", direction: " << ray.direction() << ")\n";
            }
            else
            {
                std::cerr << "Unknown edge type.\n";
            }
        }
        std::cerr << " bipolar edges: ";
        for (auto be : vf.bipolar_edge_indices) {
            std::cerr << be << " ";
        }
        std::cerr << "\n";
    }

    const int want = (vf.bipolar_match_method == BIPOLAR_MATCH_METHOD::SEP_POS) ? +1 : -1;

    // Find anchor: first bipolar slot with desired sign; if missing, fall back to 0
    int startPos = -1;
    for (int i = 0; i < nB; ++i)
    {
        if (edgeType[vf.bipolar_edge_indices[i]] == want)
        {
            startPos = i;
            break;
        }
    }
    if (startPos < 0)
        startPos = 0;

    // Ensure even number by dropping one if odd (log it)
    int usable = (nB & 1) ? (nB - 1) : nB;
    if (usable < nB)
    {
        if (debug)
        {
            std::cerr << "[ISO-MATCH] gf=" << vf.index << " odd bipolar slots=" << nB << " → dropping last\n";
        }
        if (usable <= 0)
            return; // nothing to pair
    }

    std::vector<int> orderedSlots;
    orderedSlots.reserve(usable);
    for (int t = 0; t < usable; ++t)
    {
        orderedSlots.push_back(vf.bipolar_edge_indices[(startPos + t) % nB]);
    }

    // Disjoint pairing: (s0, s1), (s2, s3), ...
    for (size_t idx = 0; idx + 1 < orderedSlots.size(); idx += 2)
    {
        vf.bipolar_matches.emplace_back(orderedSlots[idx], orderedSlots[idx + 1]);
    }

    if (debug)
    {
        std::cerr << "[ISO-MATCH] gf=" << vf.index << " method=" << matchMethodToString(vf.bipolar_match_method)
                  << " nBip=" << nB << " pairs=" << vf.bipolar_matches.size() << "\n";
        if (vf.bipolar_matches.empty())
        {
            std::cerr << "[ISO-MATCH] gf=" << vf.index << " matches: (none)\n";
        }
        else
        {
            std::cerr << "[ISO-MATCH] gf=" << vf.index << " matches:";
            for (auto &pr : vf.bipolar_matches)
                std::cerr << " (" << pr.first << "," << pr.second << ")";
            std::cerr << "\n";
        }
    }
}

// Public wrapper: recompute matches for a single global facet using its current method.
void recompute_bipolar_matches_for_facet(VoronoiDiagram &vd, int vfi, float isovalue)
{
    if (vfi < 0 || vfi >= (int)vd.global_facets.size())
        return;
    match_facet_bipolar_edges(vd, vd.global_facets[vfi], isovalue);
}

void VoronoiDiagram::compute_bipolar_matches(float isovalue)
{
    for (auto &vf : global_facets)
    {
        match_facet_bipolar_edges(*this, vf, isovalue);
    }
}

std::vector<int> VoronoiDiagram::get_vertices_for_facet(int cell_facet_index) const
{
    int vfi = facets[cell_facet_index].voronoi_facet_index;
    std::vector<int> vert = global_facets[vfi].vertices_indices;
    if (facets[cell_facet_index].orientation == -1)
    {
        std::reverse(vert.begin(), vert.end());
    }
    return vert;
}

std::string matchMethodToString(BIPOLAR_MATCH_METHOD method)
{
    switch (method)
    {
    case BIPOLAR_MATCH_METHOD::SEP_POS:
        return "SEP_POS";
    case BIPOLAR_MATCH_METHOD::SEP_NEG:
        return "SEP_NEG";
    case BIPOLAR_MATCH_METHOD::UNCONSTRAINED_MATCH:
        return "UNCONSTRAINED_MATCH";
    case BIPOLAR_MATCH_METHOD::UNDEFINED_MATCH_TYPE:
        return "UNDEFINED_MATCH_TYPE";
    default:
        return "UNKNOWN";
    }
}
