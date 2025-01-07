#include "vdc_func.h"

std::vector<DelaunayTriangle> computeDualTriangles(std::vector<CGAL::Object> &voronoi_edges, std::map<Point, float> &vertexValueMap, CGAL::Epick::Iso_cuboid_3 &bbox, std::map<Object, std::vector<Facet>, ObjectComparator> &delaunay_facet_to_voronoi_edge_map, Delaunay &dt, ScalarGrid &grid)
{

    std::vector<DelaunayTriangle> dualTriangles;
    for (const auto &edge : voronoi_edges)
    {
        Object intersectObj;
        Segment3 seg, iseg;
        Ray3 ray;
        Line3 line;
        Point3 v1, v2, ip;
        Vector3 vec1, vec2, norm;
        bool isFinite = false;

        if (CGAL::assign(seg, edge))
        {
            // If the edge is a segment
            v1 = seg.source();
            v2 = seg.target();

            // Check if it's bipolar
            // If the edge is a segment the two ends must be both in voronoi_vertices so their scalar values are pre-calculated
            if (is_bipolar(vertexValueMap[v1], vertexValueMap[v2], isovalue))
            {

                bipolar_voronoi_edges.push_back(edge); // TODO: Find the Delaunay Triangle dual to the edge

                intersectObj = CGAL::intersection(bbox, Ray3(seg.source(), v2 - v1));
                CGAL::assign(iseg, intersectObj);
                Point intersection_point = iseg.target();
                CGAL::Orientation o;
                Point positive;

                Point p1 = seg.source();
                Point p2 = seg.target();

                if (vertexValueMap[v1] >= vertexValueMap[v2])
                {
                    positive = v1;
                }
                else
                {
                    positive = v2;
                }

                for (const auto &facet : delaunay_facet_to_voronoi_edge_map[edge])
                {
                    int iFacet = facet.second;
                    Cell_handle c = facet.first;
                    int d1, d2, d3;
                    d1 = (iFacet + 1) % 4;
                    d2 = (iFacet + 2) % 4;
                    d3 = (iFacet + 3) % 4;

                    Point p1 = c->vertex(d1)->point();
                    Point p2 = c->vertex(d2)->point();
                    Point p3 = c->vertex(d3)->point();

                    int iOrient = get_orientation(iFacet, v1, v2, vertexValueMap[v1], vertexValueMap[v2]);

                    if (dt.is_infinite(c))
                    {
                        if (iOrient < 0)
                        {
                            dualTriangles.push_back(DelaunayTriangle(p1, p2, p3));
                        }
                        else
                        {
                            dualTriangles.push_back(DelaunayTriangle(p1, p3, p2));
                        }
                    }
                    else
                    {
                        if (iOrient >= 0)
                        {
                            dualTriangles.push_back(DelaunayTriangle(p1, p2, p3));
                        }
                        else
                        {
                            dualTriangles.push_back(DelaunayTriangle(p1, p3, p2));
                        }
                    }
                }
            }
        }
        else if (CGAL::assign(ray, edge))
        {
            // If the edge is a ray
            intersectObj = CGAL::intersection(bbox, ray);

            if (CGAL::assign(iseg, intersectObj))
            {

                // assign a corresponding scalar value to the intersection point and check if the segment between the source and intersection point is bi-polar
                Point intersection_point = iseg.target();
                CGAL::Orientation o;
                Point positive;

                v1 = iseg.source();
                v2 = iseg.target();

                float iPt_value = trilinear_interpolate(adjust_outside_bound_points(intersection_point, grid, v1, v2), grid);

                if (vertexValueMap[iseg.source()] >= iPt_value)
                {
                    positive = v1;
                }
                else
                {
                    positive = v2;
                }

                // Check if it's bipolar
                if (is_bipolar(vertexValueMap[iseg.source()], iPt_value, isovalue))
                {

                    Point p1 = ray.source();
                    Vector3 direction = ray.direction().vector();

                    bipolar_voronoi_edges.push_back(edge);

                    for (const auto &facet : delaunay_facet_to_voronoi_edge_map[edge])
                    {

                        Facet mirror_f = dt.mirror_facet(facet);
                        Object e = dt.dual(facet);

                        int iFacet = facet.second;
                        Cell_handle c = facet.first;
                        int d1, d2, d3;
                        d1 = (iFacet + 1) % 4;
                        d2 = (iFacet + 2) % 4;
                        d3 = (iFacet + 3) % 4;

                        Point p1 = c->vertex(d1)->point();
                        Point p2 = c->vertex(d2)->point();
                        Point p3 = c->vertex(d3)->point();

                        int iOrient = get_orientation(iFacet, v1, v2, vertexValueMap[v1], iPt_value);

                        if (dt.is_infinite(c))
                        {
                            if (iOrient < 0)
                            {
                                dualTriangles.push_back(DelaunayTriangle(p1, p2, p3));
                            }
                            else
                            {
                                dualTriangles.push_back(DelaunayTriangle(p1, p3, p2));
                            }
                        }
                        else
                        {
                            if (iOrient >= 0)
                            {
                                dualTriangles.push_back(DelaunayTriangle(p1, p2, p3));
                            }
                            else
                            {
                                dualTriangles.push_back(DelaunayTriangle(p1, p3, p2));
                            }
                        }
                    }
                }
            }
        }
        else if (CGAL::assign(line, edge))
        {
            //  If the edge is a line
            Ray3 ray1(line.point(), line.direction());
            Ray3 ray2(line.point(), -line.direction());

            intersectObj = CGAL::intersection(bbox, line);
            if (CGAL::assign(iseg, intersectObj))
            {

                Point intersection1 = iseg.source();
                Point intersection2 = iseg.target();

                float iPt1_val = trilinear_interpolate(adjust_outside_bound_points(intersection1, grid, intersection1, intersection2), grid);
                float iPt2_val = trilinear_interpolate(adjust_outside_bound_points(intersection2, grid, intersection1, intersection2), grid);

                CGAL::Orientation o;
                Point positive;

                if (iPt1_val >= iPt2_val)
                {
                    positive = intersection1;
                }
                else
                {
                    positive = intersection2;
                }

                if (is_bipolar(iPt1_val, iPt2_val, isovalue))
                {

                    Point p1 = line.point(0);
                    Point p2 = line.point(1);
                    bipolar_voronoi_edges.push_back(edge);

                    // TODO: Find the Delaunay Triangle dual to the edge

                    for (const auto &facet : delaunay_facet_to_voronoi_edge_map[edge])
                    {
                        int iFacet = facet.second;
                        Cell_handle c = facet.first;
                        int d1, d2, d3;
                        d1 = (iFacet + 1) % 4;
                        d2 = (iFacet + 2) % 4;
                        d3 = (iFacet + 3) % 4;

                        Point p1 = c->vertex(d1)->point();
                        Point p2 = c->vertex(d2)->point();
                        Point p3 = c->vertex(d3)->point();

                        int iOrient = get_orientation(iFacet, intersection1, intersection2, iPt1_val, iPt2_val);
                        if (dt.is_infinite(c))
                        {
                            if (iOrient < 0)
                            {
                                dualTriangles.push_back(DelaunayTriangle(p1, p2, p3));
                            }
                            else
                            {
                                dualTriangles.push_back(DelaunayTriangle(p1, p3, p2));
                            }
                        }
                        else
                        {
                            if (iOrient >= 0)
                            {
                                dualTriangles.push_back(DelaunayTriangle(p1, p2, p3));
                            }
                            else
                            {
                                dualTriangles.push_back(DelaunayTriangle(p1, p3, p2));
                            }
                        }
                    }
                }
            }
        }
    }
    return dualTriangles;
} // TODO: Clean up the code, and solve the orientation issue

void computeDualTrianglesMulti(
    VoronoiDiagram &voronoiDiagram,
    CGAL::Epick::Iso_cuboid_3 &bbox,
    std::map<CGAL::Object, std::vector<Facet>, ObjectComparator> &delaunay_facet_to_voronoi_edge_map,
    ScalarGrid &grid,
    float isovalue)
{
    for (const auto &edge : voronoiDiagram.voronoiEdges)
    {
        Segment3 seg;
        Ray3 ray;
        Line3 line;

        if (CGAL::assign(seg, edge))
        {
            // Edge is a segment
            Point v1 = seg.source();
            Point v2 = seg.target();

            int idx_v1 = voronoiDiagram.point_to_vertex_index[v1];
            int idx_v2 = voronoiDiagram.point_to_vertex_index[v2];

            float val1 = voronoiDiagram.voronoiVertexValues[idx_v1];
            float val2 = voronoiDiagram.voronoiVertexValues[idx_v2];

            if (is_bipolar(val1, val2, isovalue))
            {
                // TODO: Rename to edge->facet map
                auto it = delaunay_facet_to_voronoi_edge_map.find(edge);
                if (it != delaunay_facet_to_voronoi_edge_map.end())
                {
                    const std::vector<Facet> &facets = it->second;
                    for (const auto &facet : facets)
                    {
                        int iFacet = facet.second;
                        Cell_handle c = facet.first;
                        int d1 = (iFacet + 1) % 4;
                        int d2 = (iFacet + 2) % 4;
                        int d3 = (iFacet + 3) % 4;

                        Vertex_handle vh1 = c->vertex(d1);
                        Vertex_handle vh2 = c->vertex(d2);
                        Vertex_handle vh3 = c->vertex(d3);

                        if (vh1->info() || vh2->info() || vh3->info())
                        {
                            continue;
                        }

                        int cellIndex1 = voronoiDiagram.delaunayVertex_to_voronoiCell_index[vh1];
                        int cellIndex2 = voronoiDiagram.delaunayVertex_to_voronoiCell_index[vh2];
                        int cellIndex3 = voronoiDiagram.delaunayVertex_to_voronoiCell_index[vh3];

                        VoronoiCell &vc1 = voronoiDiagram.voronoiCells[cellIndex1];
                        VoronoiCell &vc2 = voronoiDiagram.voronoiCells[cellIndex2];
                        VoronoiCell &vc3 = voronoiDiagram.voronoiCells[cellIndex3];

                        // For simplicity, take the first isovertex in each cell
                        int idx1 = vc1.isoVertexStartIndex;
                        int idx2 = vc2.isoVertexStartIndex;
                        int idx3 = vc3.isoVertexStartIndex;

                        if (idx1 != idx2 && idx2 != idx3 && idx1 != idx3)
                        {
                            isoTriangles.emplace_back(idx1, idx2, idx3);
                        }
                        else
                        {
                            std::cout << "Problematic triangle" << std::endl;
                            std::cout << "Vertex 1: " << idx1 << " from Cell " << cellIndex1 << std::endl;
                            std::cout << "Vertex 2: " << idx2 << " from Cell " << cellIndex2 << std::endl;
                            std::cout << "Vertex 3: " << idx3 << " from Cell " << cellIndex3 << std::endl;
                        }
                    }
                }
            }
        }
        else if (CGAL::assign(ray, edge))
        {
            // Edge is a ray
            CGAL::Object intersectObj = CGAL::intersection(bbox, ray);
            Segment3 iseg;
            if (CGAL::assign(iseg, intersectObj))
            {
                Point v1 = ray.source();
                Point v2 = iseg.target();

                int idx_v1 = voronoiDiagram.point_to_vertex_index[v1];
                float val1 = voronoiDiagram.voronoiVertexValues[idx_v1];
                float val2 = trilinear_interpolate(v2, grid);

                if (is_bipolar(val1, val2, isovalue))
                {
                    auto it = delaunay_facet_to_voronoi_edge_map.find(edge);
                    if (it != delaunay_facet_to_voronoi_edge_map.end())
                    {
                        const std::vector<Facet> &facets = it->second;
                        for (const auto &facet : facets)
                        {
                            int iFacet = facet.second;
                            Cell_handle c = facet.first;
                            int d1 = (iFacet + 1) % 4;
                            int d2 = (iFacet + 2) % 4;
                            int d3 = (iFacet + 3) % 4;

                            Vertex_handle vh1 = c->vertex(d1);
                            Vertex_handle vh2 = c->vertex(d2);
                            Vertex_handle vh3 = c->vertex(d3);

                            if (vh1->info() || vh2->info() || vh3->info())
                            {
                                continue;
                            }

                            int cellIndex1 = voronoiDiagram.delaunayVertex_to_voronoiCell_index[vh1];
                            int cellIndex2 = voronoiDiagram.delaunayVertex_to_voronoiCell_index[vh2];
                            int cellIndex3 = voronoiDiagram.delaunayVertex_to_voronoiCell_index[vh3];

                            VoronoiCell &vc1 = voronoiDiagram.voronoiCells[cellIndex1];
                            VoronoiCell &vc2 = voronoiDiagram.voronoiCells[cellIndex2];
                            VoronoiCell &vc3 = voronoiDiagram.voronoiCells[cellIndex3];

                            // For simplicity, take the first isovertex in each cell
                            int idx1 = vc1.isoVertexStartIndex;
                            int idx2 = vc2.isoVertexStartIndex;
                            int idx3 = vc3.isoVertexStartIndex;

                            if (idx1 != idx2 && idx2 != idx3 && idx1 != idx3)
                            {
                                isoTriangles.emplace_back(idx1, idx2, idx3);
                            }
                            else
                            {
                                std::cout << "Problematic triangle" << std::endl;
                                std::cout << "Vertex 1: " << idx1 << " from Cell " << cellIndex1 << std::endl;
                                std::cout << "Vertex 2: " << idx2 << " from Cell " << cellIndex2 << std::endl;
                                std::cout << "Vertex 3: " << idx3 << " from Cell " << cellIndex3 << std::endl;
                            }
                        }
                    }
                }
            }
        }
        else if (CGAL::assign(line, edge))
        {
            // Edge is a line
            CGAL::Object intersectObj = CGAL::intersection(bbox, line);
            Segment3 iseg;
            if (CGAL::assign(iseg, intersectObj))
            {
                Point v1 = iseg.source();
                Point v2 = iseg.target();

                float val1 = trilinear_interpolate(v1, grid);
                float val2 = trilinear_interpolate(v2, grid);

                if (is_bipolar(val1, val2, isovalue))
                {
                    auto it = delaunay_facet_to_voronoi_edge_map.find(edge);
                    if (it != delaunay_facet_to_voronoi_edge_map.end())
                    {
                        const std::vector<Facet> &facets = it->second;
                        for (const auto &facet : facets)
                        {
                            int iFacet = facet.second;
                            Cell_handle c = facet.first;
                            int d1 = (iFacet + 1) % 4;
                            int d2 = (iFacet + 2) % 4;
                            int d3 = (iFacet + 3) % 4;

                            Vertex_handle vh1 = c->vertex(d1);
                            Vertex_handle vh2 = c->vertex(d2);
                            Vertex_handle vh3 = c->vertex(d3);

                            if (vh1->info() || vh2->info() || vh3->info())
                            {
                                continue;
                            }

                            int cellIndex1 = voronoiDiagram.delaunayVertex_to_voronoiCell_index[vh1];
                            int cellIndex2 = voronoiDiagram.delaunayVertex_to_voronoiCell_index[vh2];
                            int cellIndex3 = voronoiDiagram.delaunayVertex_to_voronoiCell_index[vh3];

                            VoronoiCell &vc1 = voronoiDiagram.voronoiCells[cellIndex1];
                            VoronoiCell &vc2 = voronoiDiagram.voronoiCells[cellIndex2];
                            VoronoiCell &vc3 = voronoiDiagram.voronoiCells[cellIndex3];

                            // For simplicity, take the first isovertex in each cell
                            int idx1 = vc1.isoVertexStartIndex;
                            int idx2 = vc2.isoVertexStartIndex;
                            int idx3 = vc3.isoVertexStartIndex;

                            if (idx1 != idx2 && idx2 != idx3 && idx1 != idx3)
                            {
                                isoTriangles.emplace_back(idx1, idx2, idx3);
                            }
                            else
                            {
                                std::cout << "Problematic triangle" << std::endl;
                                std::cout << "Vertex 1: " << idx1 << " from Cell " << cellIndex1 << std::endl;
                                std::cout << "Vertex 2: " << idx2 << " from Cell " << cellIndex2 << std::endl;
                                std::cout << "Vertex 3: " << idx3 << " from Cell " << cellIndex3 << std::endl;
                            }
                        }
                    }
                }
            }
        }
    }
}

void Compute_Isosurface_Vertices_Single(VoronoiDiagram &voronoiDiagram, ScalarGrid &grid)
{
    const int cubeVertices[8][3] = {
        {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};

    const int cubeEdges[12][2] = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6}, {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};

    for (const auto &center : activeCubeCenters)
    {
        std::vector<Point> intersectionPoints;
        float cubeSize = grid.dx; // Assuming grid.dx is the cube size
        std::array<float, 8> scalarValues;

        // Compute scalar values at cube vertices
        for (int i = 0; i < 8; i++)
        {
            Point vertex(
                center.x() + (cubeVertices[i][0] - 0.5f) * cubeSize,
                center.y() + (cubeVertices[i][1] - 0.5f) * cubeSize,
                center.z() + (cubeVertices[i][2] - 0.5f) * cubeSize);
            scalarValues[i] = grid.get_scalar_value_at_point(vertex);
        }

        // Check each edge for intersection with the isovalue
        for (const auto &edge : cubeEdges)
        {
            int idx1 = edge[0];
            int idx2 = edge[1];
            float val1 = scalarValues[idx1];
            float val2 = scalarValues[idx2];

            if (is_bipolar(val1, val2, isovalue))
            {
                Point p1(
                    center.x() + (cubeVertices[idx1][0] - 0.5f) * cubeSize,
                    center.y() + (cubeVertices[idx1][1] - 0.5f) * cubeSize,
                    center.z() + (cubeVertices[idx1][2] - 0.5f) * cubeSize);

                Point p2(
                    center.x() + (cubeVertices[idx2][0] - 0.5f) * cubeSize,
                    center.y() + (cubeVertices[idx2][1] - 0.5f) * cubeSize,
                    center.z() + (cubeVertices[idx2][2] - 0.5f) * cubeSize);

                Point intersect = interpolate(p1, p2, val1, val2, isovalue, data_grid);
                intersectionPoints.push_back(intersect);
            }
        }

        // Compute the centroid of the intersection points
        if (!intersectionPoints.empty())
        {
            Point centroid = compute_centroid(intersectionPoints);
            voronoiDiagram.isosurfaceVertices.push_back(centroid);
        }
    }
}

void Compute_Isosurface_Vertices_Multi(VoronoiDiagram &voronoiDiagram, float isovalue)
{
    for (auto &vc : voronoiDiagram.voronoiCells)
    {
        std::vector<MidpointNode> midpoints;
        std::map<std::pair<int, int>, int> edge_to_midpoint_index;

        // First pass: Collect midpoints and build edge connectivity
        for (size_t i = 0; i < vc.facet_indices.size(); ++i)
        {
            int facet_index = vc.facet_indices[i];
            VoronoiFacet &facet = voronoiDiagram.voronoiFacets[facet_index];
            size_t num_vertices = facet.vertices_indices.size();

            // Store indices of midpoints in this facet
            std::vector<int> facet_midpoint_indices;

            for (size_t j = 0; j < num_vertices; ++j)
            {
                size_t idx1 = j;
                size_t idx2 = (j + 1) % num_vertices;

                float val1 = facet.vertex_values[idx1];
                float val2 = facet.vertex_values[idx2];

                // Check for bipolar edge
                if (is_bipolar(val1, val2, isovalue))
                {
                    int vertex_index1 = facet.vertices_indices[idx1];
                    int vertex_index2 = facet.vertices_indices[idx2];

                    Point p1 = voronoiDiagram.voronoiVertices[vertex_index1].vertex;
                    Point p2 = voronoiDiagram.voronoiVertices[vertex_index2].vertex;

                    double t = (isovalue - val1) / (val2 - val1);
                    Point midpoint = p1 + (p2 - p1) * t;

                    auto edge_key = std::make_pair(std::min(vertex_index1, vertex_index2),
                                                   std::max(vertex_index1, vertex_index2));

                    if (edge_to_midpoint_index.find(edge_key) == edge_to_midpoint_index.end())
                    {
                        MidpointNode node;
                        node.point = midpoint;
                        midpoints.push_back(node);
                        int midpoint_index = midpoints.size() - 1;
                        edge_to_midpoint_index[edge_key] = midpoint_index;
                        facet_midpoint_indices.push_back(midpoint_index);
                    }
                    else
                    {
                        int midpoint_index = edge_to_midpoint_index[edge_key];
                        facet_midpoint_indices.push_back(midpoint_index);
                    }
                }
            }

            // Connect midpoints in this facet
            size_t num_midpoints = facet_midpoint_indices.size();
            for (size_t k = 0; k + 1 < num_midpoints; k += 2)
            {
                int idx1 = facet_midpoint_indices[k];
                int idx2 = facet_midpoint_indices[k + 1];
                midpoints[idx1].connected_to.push_back(idx2);
                midpoints[idx2].connected_to.push_back(idx1);
            }
        }

        // Extract cycles from the graph of midpoints
        std::vector<std::vector<int>> cycles_indices;
        std::set<int> visited;

        for (size_t i = 0; i < midpoints.size(); ++i)
        {
            if (visited.find(i) == visited.end())
            {
                std::vector<int> cycle;
                std::stack<int> stack;
                stack.push(i);

                while (!stack.empty())
                {
                    int current = stack.top();
                    stack.pop();

                    if (visited.find(current) != visited.end())
                    {
                        continue;
                    }

                    visited.insert(current);
                    cycle.push_back(current);

                    for (int neighbor : midpoints[current].connected_to)
                    {
                        if (visited.find(neighbor) == visited.end())
                        {
                            stack.push(neighbor);
                        }
                    }
                }

                if (!cycle.empty())
                {
                    cycles_indices.push_back(cycle);
                }
            }
        }

        // For each cycle, compute the centroid and store isoVertices
        vc.isoVertexStartIndex = voronoiDiagram.isosurfaceVertices.size();
        vc.numIsoVertices = cycles_indices.size();

        for (const auto &cycle_indices : cycles_indices)
        {
            Cycle cycle;
            cycle.voronoi_cell_index = vc.cellIndex;
            cycle.midpoint_indices = cycle_indices;

            // Store edges (pairs of indices into midpoints)
            for (size_t i = 0; i < cycle_indices.size(); ++i)
            {
                int idx1 = cycle_indices[i];
                int idx2 = cycle_indices[(i + 1) % cycle_indices.size()];
                cycle.edges.emplace_back(idx1, idx2);
            }

            // Compute centroid using the midpoints
            cycle.compute_centroid(midpoints);

            vc.cycles.push_back(cycle);
            voronoiDiagram.isosurfaceVertices.push_back(cycle.isovertex);
        }

        if (debug)
        {
            std::cout << vc << std::endl;
        }
    }

    if (debug)
    {
        for (size_t i = 0; i < voronoiDiagram.isosurfaceVertices.size(); ++i)
        {
            std::cout << "  Index " << i << ": " << voronoiDiagram.isosurfaceVertices[i] << "\n";
        }
    }
}

std::vector<Point> add_dummy_from_facet(const GRID_FACETS &facet, const Grid &data_grid)
{
    std::vector<Point> points;

    // 2D slice dimension
    int dim0 = facet.axis_size[0]; 
    int dim1 = facet.axis_size[1];

    // For convenience
    int d  = facet.orth_dir;
    int d1 = facet.axis_dir[0];
    int d2 = facet.axis_dir[1];

    // We have bounding-box in facet.minIndex[], facet.maxIndex[], 
    // and localSize[] = (maxIndex[i] - minIndex[i] + 1)
    // The grid spacing in each dimension
    double dx[3] = { data_grid.dx, data_grid.dy, data_grid.dz };

    // Loop over the 2D slice
    for (int coord1 = 0; coord1 < dim1; coord1++)
    {
        for (int coord0 = 0; coord0 < dim0; coord0++)
        {
            if (!facet.CubeFlag(coord0, coord1)) 
                continue;

            // localX[d1] = coord0, localX[d2] = coord1
            int localX[3] = {0,0,0};
            localX[d1] = coord0;
            localX[d2] = coord1;

            // side=0 => localX[d] = 0, side=1 => localX[d] = localSize[d]-1
            localX[d] = (facet.side == 0) ? 0 : (facet.localSize[d] - 1);

            // Convert localX -> global indices
            int g[3];
            for (int i = 0; i < 3; i++) {
                g[i] = localX[i] + facet.minIndex[i];
            }

            // Compute center in real-world coordinates
            double cx = (g[0] + 0.5) * dx[0];
            double cy = (g[1] + 0.5) * dx[1];
            double cz = (g[2] + 0.5) * dx[2];

            // Offset by +/- dx[d]
            double offset = (facet.side == 0) ? -dx[d] : dx[d];
            if      (d == 0) cx += offset;
            else if (d == 1) cy += offset;
            else             cz += offset;

            points.emplace_back(cx, cy, cz);
        }
    }

    return points;
}

void construct_delaunay_triangulation(Grid &grid, const std::vector<std::vector<GRID_FACETS>> &grid_facets)
{
    if (multi_isov)
    {
        // Add original points
        for (const auto &p : activeCubeCenters)
        {
            delaunay_vertices.push_back({p, false});
        }

        // Add dummy points to the form of triangulation to bound the voronoi diagram
        std::vector<Point> dummy_points;

        // refers to the grid facets
        for (int d = 0; d < 3; d++)
        {
            for (const auto &f : grid_facets[d])
            {
                std::vector<Point> pointsf = add_dummy_from_facet(f, grid);
                dummy_points.insert(dummy_points.end(), pointsf.begin(), pointsf.end());
            }
        }

        /*      // brute-forcely adding dummy points on all 6 faces corresponding to grid spacing

                // Face x = xmin-dx and xmax+dx
                for ( int i = 0; i <= ny; ++i) {
                    double y = ymin + i * grid.dy;
                    for (int j = 0; j <= nz; ++j) {
                        double z = zmin + j * grid.dz;
                        dummy_points.push_back(Point( xmin - grid.dx, y, z));
                        dummy_points.push_back(Point( xmax + grid.dx, y, z));
                    }
                }

                // Face y = ymin-dy and ymax+dy
                for ( int i = 0; i <= nx; ++i) {
                    double x = xmin + i * grid.dx;
                    for (int j = 0; j <= nz; ++j) {
                        double z = zmin + j * grid.dz;
                        dummy_points.push_back(Point( x, ymin - grid.dy, z));
                        dummy_points.push_back(Point( x, ymax + grid.dy, z));
                    }
                }

                // Face z = zmin-dz and zmax+dz
                for ( int i = 0; i <= nx; ++i) {
                    double x = xmin + i * grid.dx;
                    for (int j = 0; j <= ny; ++j) {
                        double y = ymin + j * grid.dy;
                        dummy_points.push_back(Point( x, y, zmin - grid.dz));
                        dummy_points.push_back(Point( x, y, zmax + grid.dz));
                    }
                } */

        /*
         Temp method of writing dummypoints to a csv file for debug
        */
        if (true)
        {
            std::ofstream ofs("dummy_points.csv");

            ofs << grid.nx << ","
                << grid.ny << ","
                << grid.nz << ","
                << grid.dx << ","
                << grid.dy << ","
                << grid.dz << "\n";
            ofs << "x,y,z\n";

            for (const auto &p : dummy_points)
            {
                ofs << p.x() << "," << p.y() << "," << p.z() << "\n";
            }

            ofs.close();
        }

        int count = 0;
        for (const auto &dp : dummy_points)
        {
            delaunay_vertices.push_back({dp, true});
            count++;
        }

        std::cout << "Number of Dummy Points added: " << count << std::endl;

        for (const auto &lp : delaunay_vertices)
        {
            points_with_info.emplace_back(lp.point, lp.is_dummy);
        }

        dt.insert(points_with_info.begin(), points_with_info.end());
    }
    else
    {
        dt.insert(activeCubeCenters.begin(), activeCubeCenters.end());
    }

    int index = 0;

    int i = 0;
    if (multi_isov)
    {
        for (auto vh = dt.finite_vertices_begin(); vh != dt.finite_vertices_end(); ++vh)
        {
            if (vh->info() == true)
            {
                continue; // Skip dummy points
            }
            Point p = vh->point();
            point_index_map[p] = i;
            i++;
        }
    }
    else
    {
        for (const auto &pt : activeCubeCenters)
        {
            point_index_map[pt] = i;
            i++;
        }
    }
}

void construct_voronoi_vertices(VoronoiDiagram &voronoiDiagram)
{
    voronoiDiagram.voronoiVertices.clear();
    voronoiDiagram.point_to_vertex_index.clear();

    std::set<Point> seen_points;
    for (Delaunay::Finite_cells_iterator cit = dt.finite_cells_begin();
         cit != dt.finite_cells_end(); ++cit)
    {

        Point voronoi_vertex = dt.dual(cit);
        VoronoiVertex vVertex(voronoi_vertex);
        if (seen_points.insert(voronoi_vertex).second)
        {
            int vertex_index = voronoiDiagram.voronoiVertices.size();
            voronoiDiagram.voronoiVertices.push_back(vVertex);
            voronoiDiagram.point_to_vertex_index[voronoi_vertex] = vertex_index;
            voronoiDiagram.cell_to_vertex_index[cit] = vertex_index;
        }
    }
}

void compute_voronoi_values(VoronoiDiagram &voronoiDiagram, ScalarGrid &grid)
{
    voronoiDiagram.voronoiVertexValues.resize(voronoiDiagram.voronoiVertices.size());
    for (size_t i = 0; i < voronoiDiagram.voronoiVertices.size(); ++i)
    {
        Point vertex = voronoiDiagram.voronoiVertices[i].vertex;
        float value = trilinear_interpolate(vertex, grid);
        voronoiDiagram.voronoiVertexValues[i] = value;
        vertexValueMap[vertex] = value;
    }
}

void construct_voronoi_cells(VoronoiDiagram &voronoiDiagram)
{
    int index = 0;
    for (auto vh = dt.finite_vertices_begin(); vh != dt.finite_vertices_end(); ++vh)
    {
        if (vh->info())
        {
            // std::cout << "Dummy Point excluded: " << vh->point() << std::endl;
            continue;
        }
        VoronoiCell vc(vh);
        vc.cellIndex = index;

        std::vector<Cell_handle> incident_cells;
        dt.finite_incident_cells(vh, std::back_inserter(incident_cells));

        CGAL::Bbox_3 domain_bbox = delaunayBbox.bbox();

        // Collect vertex indices, ensuring uniqueness
        std::set<int> unique_vertex_indices_set;
        for (Cell_handle ch : incident_cells)
        {
            if (dt.is_infinite(ch))
            {
                continue; // Skip infinite cells
            }
            Point voronoi_vertex = dt.dual(ch);

            // Check if voronoi_vertex is within domain and exclude dummy points in the dt
            int vertex_index = voronoiDiagram.point_to_vertex_index[voronoi_vertex];
            unique_vertex_indices_set.insert(vertex_index);
        }

        // Copy unique indices to vector
        vc.vertices_indices.assign(unique_vertex_indices_set.begin(), unique_vertex_indices_set.end());

        // Build convex hull and extract facets
        std::vector<Point> vertex_points;
        for (int idx : vc.vertices_indices)
        {
            vertex_points.push_back(voronoiDiagram.voronoiVertices[idx].vertex);
        }

        // Remove duplicate points
        std::sort(vertex_points.begin(), vertex_points.end(), [](const Point &a, const Point &b)
                  { return a.x() < b.x() || (a.x() == b.x() && (a.y() < b.y() || (a.y() == b.y() && a.z() < b.z()))); });
        vertex_points.erase(std::unique(vertex_points.begin(), vertex_points.end(), PointApproxEqual()), vertex_points.end());

        CGAL::convex_hull_3(vertex_points.begin(), vertex_points.end(), vc.polyhedron);

        // Extract facets from polyhedron
        for (auto facet_it = vc.polyhedron.facets_begin();
             facet_it != vc.polyhedron.facets_end(); ++facet_it)
        {
            VoronoiFacet vf;
            auto h = facet_it->facet_begin();
            do
            {
                Point p = h->vertex()->point();
                int vertex_index = voronoiDiagram.point_to_vertex_index[p];
                vf.vertices_indices.push_back(vertex_index);
                float value = voronoiDiagram.voronoiVertexValues[vertex_index];
                vf.vertex_values.push_back(value);
                ++h;
            } while (h != facet_it->facet_begin());

            int facet_index = voronoiDiagram.voronoiFacets.size();
            voronoiDiagram.voronoiFacets.push_back(vf);
            vc.facet_indices.push_back(facet_index);
        }

        voronoiDiagram.voronoiCells.push_back(vc);
        voronoiDiagram.delaunayVertex_to_voronoiCell_index[vh] = index;
        index++;
    }
}

void construct_voronoi_edges(
    VoronoiDiagram &voronoiDiagram,
    std::map<CGAL::Object, std::vector<Facet>, ObjectComparator> &delaunay_facet_to_voronoi_edge_map)
{
    std::set<std::string> seen_edges; // Used to check for duplicate edges

    for (Delaunay::Finite_facets_iterator fit = dt.finite_facets_begin();
         fit != dt.finite_facets_end(); ++fit)
    {
        Facet facet = *fit;

        CGAL::Object vEdge = dt.dual(facet);

        if (isDegenerate(vEdge))
        {
            continue;
        }

        std::string edgeRep = objectToString(vEdge);

        delaunay_facet_to_voronoi_edge_map[vEdge].push_back(facet);

        if (seen_edges.find(edgeRep) == seen_edges.end())
        {
            voronoiDiagram.voronoiEdges.push_back(vEdge);
            seen_edges.insert(edgeRep);
        }
    }
}

int handle_output_mesh(bool &retFlag, VoronoiDiagram &vd)
{
    retFlag = true;
    // Use locations of isosurface vertices as vertices of Delaunay triangles and write the output mesh
    if (multi_isov)
    {
        if (output_format == "off")
        {
            writeOFFMulti(output_filename, vd, isoTriangles);
        }
        else if (output_format == "ply")
        {
            writePLYMulti(output_filename, vd, isoTriangles);
        }
        else
        {
            std::cerr << "Unsupported output format: " << output_format << std::endl;
            return EXIT_FAILURE;
        }
    }
    else
    {
        if (output_format == "off")
        {
            writeOFFSingle(output_filename, vd.isosurfaceVertices, dualTriangles, point_index_map);
        }
        else if (output_format == "ply")
        {
            writePLYSingle(output_filename, vd.isosurfaceVertices, dualTriangles, point_index_map);
        }
        else
        {
            std::cerr << "Unsupported output format: " << output_format << std::endl;
            return EXIT_FAILURE;
        }
    }
    retFlag = false;
    return {};
}