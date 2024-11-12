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


void computeDualTrianglesMulti(std::vector<CGAL::Object> &voronoi_edges, CGAL::Epick::Iso_cuboid_3 &bbox, std::map<Object, std::vector<Facet>, ObjectComparator> &delaunay_facet_to_voronoi_edge_map, ScalarGrid &grid)
{

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

                for (const auto& facet : delaunay_facet_to_voronoi_edge_map[edge]) {
                    // Get the Delaunay triangle vertices
                    int iFacet = facet.second;
                    Cell_handle c = facet.first;
                    int d1 = (iFacet + 1) % 4;
                    int d2 = (iFacet + 2) % 4;
                    int d3 = (iFacet + 3) % 4;

                    Point p1 = c->vertex(d1)->point();
                    Point p2 = c->vertex(d2)->point();
                    Point p3 = c->vertex(d3)->point();

                    // For each vertex, get the associated isovertices
                    const auto& iso_indices_p1 = vertex_to_isovertex_indices[p1];
                    const auto& iso_indices_p2 = vertex_to_isovertex_indices[p2];
                    const auto& iso_indices_p3 = vertex_to_isovertex_indices[p3];

                    //decide which isovertices to use for this triangle

                    // For simplicity, let's assume one isovertex per Voronoi cell
                    int idx1 = iso_indices_p1[0]; // Or select appropriately
                    int idx2 = iso_indices_p2[0];
                    int idx3 = iso_indices_p3[0];

                    // Create the isoTriangle
                    isoTriangles.emplace_back(idx1, idx2, idx3);
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

                    for (const auto& facet : delaunay_facet_to_voronoi_edge_map[edge]) {
                        // Get the Delaunay triangle vertices
                        int iFacet = facet.second;
                        Cell_handle c = facet.first;
                        int d1 = (iFacet + 1) % 4;
                        int d2 = (iFacet + 2) % 4;
                        int d3 = (iFacet + 3) % 4;

                        Point p1 = c->vertex(d1)->point();
                        Point p2 = c->vertex(d2)->point();
                        Point p3 = c->vertex(d3)->point();

                        // For each vertex, get the associated isovertices
                        const auto& iso_indices_p1 = vertex_to_isovertex_indices[p1];
                        const auto& iso_indices_p2 = vertex_to_isovertex_indices[p2];
                        const auto& iso_indices_p3 = vertex_to_isovertex_indices[p3];

                        // For simplicity, let's assume one isovertex per Voronoi cell
                        int idx1 = iso_indices_p1[0]; 
                        int idx2 = iso_indices_p2[0];
                        int idx3 = iso_indices_p3[0];

                        // Create the isoTriangle
                        isoTriangles.emplace_back(idx1, idx2, idx3);
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

                    for (const auto& facet : delaunay_facet_to_voronoi_edge_map[edge]) {
                        // Get the Delaunay triangle vertices
                        int iFacet = facet.second;
                        Cell_handle c = facet.first;
                        int d1 = (iFacet + 1) % 4;
                        int d2 = (iFacet + 2) % 4;
                        int d3 = (iFacet + 3) % 4;

                        Point p1 = c->vertex(d1)->point();
                        Point p2 = c->vertex(d2)->point();
                        Point p3 = c->vertex(d3)->point();

                        // For each vertex, get the associated isovertices
                        const auto& iso_indices_p1 = vertex_to_isovertex_indices[p1];
                        const auto& iso_indices_p2 = vertex_to_isovertex_indices[p2];
                        const auto& iso_indices_p3 = vertex_to_isovertex_indices[p3];

                        // decide which isovertices to use for this triangle

                        // For simplicity, let's assume one isovertex per Voronoi cell
                        int idx1 = iso_indices_p1[0]; // Or select appropriately
                        int idx2 = iso_indices_p2[0];
                        int idx3 = iso_indices_p3[0];

                        // Create the isoTriangle
                        isoTriangles.emplace_back(idx1, idx2, idx3);
                    }
                }
            }
        }
    }
}
void Compute_Isosurface_Vertices_Single(ScalarGrid &grid)
{
    /*
       Compute Isosurface Vertices
       */
    const int cubeVertices[8][3] = {
        {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, // bottom face
        {0, 0, 1},
        {1, 0, 1},
        {1, 1, 1},
        {0, 1, 1} // top face
    };

    const int cubeEdges[12][2] = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0}, // bottom face
        {4, 5},
        {5, 6},
        {6, 7},
        {7, 4}, // top face
        {0, 4},
        {1, 5},
        {2, 6},
        {3, 7} // vertical edges
    };
    for (const auto &center : activeCubeCenters)
    {
        std::vector<Point> intersectionPoints;
        float cubeSize = data_grid.dx;
        std::array<float, 8> scalarValues;

        // Compute the global coordinates of the cube vertices and their scalar values
        for (int i = 0; i < 8; i++)
        {
            Point vertex(
                center.x() + (cubeVertices[i][0] - 0.5) * cubeSize,
                center.y() + (cubeVertices[i][1] - 0.5) * cubeSize,
                center.z() + (cubeVertices[i][2] - 0.5) * cubeSize);
            scalarValues[i] = grid.get_scalar_value_at_point(vertex);
        }

        // Check each edge for intersection with the isovalue
        for (const auto &edge : cubeEdges)
        {
            int idx1 = edge[0];
            int idx2 = edge[1];
            float val1 = scalarValues[idx1];
            float val2 = scalarValues[idx2];

            // Check if the edge crosses the isovalue
            if (((val1 >= isovalue) && (val2 < isovalue)) || ((val1 < isovalue) && (val2 >= isovalue))) // FIX: change to comparison instead of arithmatic // This is way too complicated.
            {
                Point p1(
                    center.x() + (cubeVertices[idx1][0] - 0.5) * cubeSize,
                    center.y() + (cubeVertices[idx1][1] - 0.5) * cubeSize,
                    center.z() + (cubeVertices[idx1][2] - 0.5) * cubeSize);

                Point p2(
                    center.x() + (cubeVertices[idx2][0] - 0.5) * cubeSize,
                    center.y() + (cubeVertices[idx2][1] - 0.5) * cubeSize,
                    center.z() + (cubeVertices[idx2][2] - 0.5) * cubeSize);

                Point intersect = interpolate(p1, p2, val1, val2, isovalue, data_grid);
                intersectionPoints.push_back(intersect);
            }
        }

        // Compute the centroid of the intersection points
        if (!intersectionPoints.empty())
        {
            Point centroid = compute_centroid(intersectionPoints, supersample, supersample_r);
            isosurfaceVertices.push_back(centroid);
        }
    }
}


void Compute_Isosurface_Vertices_Multi(std::vector<VoronoiCell> &voronoi_cells)
{
    /*
    Algorithm SepNeg/SepPos
    */
    std::map<std::pair<Point, Point>, EdgeMidpoint> bipolar_edge_map;
    std::set<Point> unique_vertices;

    for (auto &vc : voronoi_cells)
    {
        std::vector<MidpointNode> midpoints;
        std::map<std::pair<Point, Point>, int> edge_to_midpoint_index;

        // First pass: Collect midpoints and build edge connectivity
        for (size_t i = 0; i < vc.facets.size(); ++i)
        {
            VoronoiFacet &facet = vc.facets[i];
            size_t num_vertices = facet.vertices.size();

            // Store indices of midpoints in this facet
            std::vector<int> facet_midpoint_indices;

            for (size_t j = 0; j < num_vertices; ++j)
            {
                // Current and next vertex indices
                size_t idx1 = j;
                size_t idx2 = (j + 1) % num_vertices;

                float val1 = facet.vertex_values[idx1];
                float val2 = facet.vertex_values[idx2];

                // Check for bipolar edge
                if (is_bipolar(val1, val2, isovalue))
                {
                    // Linear interpolation to find the exact point where the isovalue is crossed
                    Point p1 = facet.vertices[idx1];
                    Point p2 = facet.vertices[idx2];

                    double t = (isovalue - val1) / (val2 - val1);
                    Point midpoint = p1 + (p2 - p1) * t;

                    // Check if this midpoint already exists
                    auto edge_key = std::make_pair(std::min(p1, p2), std::max(p1, p2));
                    if (edge_to_midpoint_index.find(edge_key) == edge_to_midpoint_index.end())
                    {
                        // Create new midpoint node
                        MidpointNode node;
                        node.point = midpoint;
                        midpoints.push_back(node);
                        int midpoint_index = midpoints.size() - 1;
                        edge_to_midpoint_index[edge_key] = midpoint_index;
                        facet_midpoint_indices.push_back(midpoint_index);
                    }
                    else
                    {
                        // Midpoint already exists
                        int midpoint_index = edge_to_midpoint_index[edge_key];
                        facet_midpoint_indices.push_back(midpoint_index);
                    }
                }
            }

            // Connect midpoints in this facet
            /* The number of midpoints should be multiple of 2 as bipolar edges in a facet should come in pairs,
            So traverse through the vertices of the facet and connect consecutive pairs together
            */
            size_t num_midpoints = facet_midpoint_indices.size();
            for (size_t k = 0; k < num_midpoints; k += 2)
            {
                int idx1 = facet_midpoint_indices[k];
                int idx2 = facet_midpoint_indices[k + 1];
                midpoints[idx1].connected_to.push_back(idx2);
                midpoints[idx2].connected_to.push_back(idx1);
            }
        }

        // Extract cycles from the graph of midpoints
        std::vector<std::vector<int>> cycles_indices;
        std::set<int> visited; // Stores the index of the vertices that are already visited

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

                // If the cycle is closed, add it
                if (!cycle.empty())
                {
                    cycles_indices.push_back(cycle);
                }
            }
        }

        // For each cycle, compute the centroid
        for (const auto &cycle_indices : cycles_indices)
        {
            std::vector<Point> cycle_points;
            for (int idx : cycle_indices)
            {
                cycle_points.push_back(midpoints[idx].point);
            }

            // Compute centroid
            Point centroid = CGAL::centroid(cycle_points.begin(), cycle_points.end());

            // Store the cycle and its isovertex
            Cycle cycle;
            cycle.midpoints = cycle_points;
            cycle.isovertex = centroid;
            vc.cycles.push_back(cycle);

            if (unique_vertices.insert(centroid).second)
            { // insert returns a pair, where .second is a bool indicating success
                isosurfaceVertices.push_back(centroid);
            }
        }
    }

    // Populate vertex_to_isovertex_indices
    for (size_t i = 0; i < voronoi_cells.size(); ++i)
    {
        VoronoiCell &vc = voronoi_cells[i];
        Vertex_handle vh = vc.delaunay_vertex;
        Point delaunay_point = vh->point();

        std::vector<int> isovertex_indices;
        for (size_t j = 0; j < vc.cycles.size(); ++j)
        {
            Cycle &cycle = vc.cycles[j];
            // Add isovertex to isosurfaceVertices if not already added
            isosurfaceVertices.push_back(cycle.isovertex);
            isovertex_index = isosurfaceVertices.size() - 1;
            isovertex_indices.push_back(isovertex_index);
        }
        vertex_to_isovertex_indices[delaunay_point] = isovertex_indices;
    }
}

void construct_delaunay_triangulation()
{
    if (multi_isov)
    {
        K::Iso_cuboid_3 delaunayBbox = CGAL::bounding_box(activeCubeCenters.begin(), activeCubeCenters.end());

        // Six Corner Points
        double xmin = delaunayBbox.xmin();
        double xmax = delaunayBbox.xmax();
        double ymin = delaunayBbox.ymin();
        double ymax = delaunayBbox.ymax();
        double zmin = delaunayBbox.zmin();
        double zmax = delaunayBbox.zmax();

        // Bounding box side length
        double lx = xmax - xmin;
        double ly = ymax - ymin;
        double lz = zmax - zmin;

        // Add original points
        for (const auto &p : activeCubeCenters)
        {
            all_points.push_back({p, false});
        }

        // Add 24 dummy points to the point set forming triangulation
        std::vector<Point> dummy_points = {
            Point(xmin - lx, ymin, zmin),
            Point(xmin, ymin - ly, zmin),
            Point(xmin, ymin, zmin - lz),
            Point(xmin - lx, ymin, zmax),
            Point(xmin, ymin - lz, zmax),
            Point(xmin, ymin, zmax + lz),
            Point(xmin - lx, ymax, zmin),
            Point(xmin, ymax + ly, zmin),
            Point(xmin, ymax, zmin - lz),
            Point(xmin - lx, ymax, zmax),
            Point(xmin, ymax + ly, zmax),
            Point(xmin, ymax, zmax + lz),
            Point(xmax + lx, ymin, zmin),
            Point(xmax, ymin - ly, zmin),
            Point(xmax, ymin, zmin - lz),
            Point(xmax + lx, ymin, zmax),
            Point(xmax, ymin - ly, zmax),
            Point(xmax, ymin, zmax + lz),
            Point(xmax + lx, ymax, zmin),
            Point(xmax, ymax + ly, zmin),
            Point(xmax, ymax, zmin - lz),
            Point(xmax + lx, ymax, zmax),
            Point(xmax, ymax + ly, zmax),
            Point(xmax, ymax, zmax + lz)};

        for (const auto &dp : dummy_points)
        {
            all_points.push_back({dp, true});
        }

        for (const auto &lp : all_points)
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
        for (const auto &pt : points_with_info)
        {
            if (pt.second == true)
            {
                continue;
            }
            Point p = pt.first; // pt in this case is pair of <Point, bool>
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

void construct_voronoi_vertices(std::set<Point> &seen_points, std::vector<Point> &voronoi_vertices)
{
    for (Delaunay::Finite_cells_iterator cit = dt.finite_cells_begin(); cit != dt.finite_cells_end(); ++cit)
    {

        Point voronoi_vertex = dt.dual(cit);

        if (seen_points.insert(voronoi_vertex).second)
        { // insert returns a pair, where .second is a bool indicating success
            voronoi_vertices.push_back(voronoi_vertex);
        }
    }
}

void construct_voronoi_cells(std::vector<VoronoiCell> &voronoi_cells)
{
    // Iterate over all finite vertices of the Delaunay Triangulation and for each vertex, collect all finite cells incident to the vertex and append to incident_cells
    for (auto vh = dt.finite_vertices_begin(); vh != dt.finite_vertices_end(); ++vh)
    {
        if (vh->info())
        {
            continue;
        }
        VoronoiCell vc(vh);

        std::vector<Cell_handle> incident_cells;
        dt.finite_incident_cells(vh, std::back_inserter(incident_cells));

        // Collect Voronoi vertices and their scalar Values
        for (Cell_handle ch : incident_cells)
        {
            if (dt.is_infinite(ch))
            {
                continue; // Skip infinite cells
            }

            /*                 bool has_dummy_vertex = false;
                            for (int i = 0; i < 4; ++i)
                            {
                                if (ch->vertex(i)->info())
                                {
                                    has_dummy_vertex = true;
                                    break;
                                }
                            }
                            if (has_dummy_vertex)
                            {
                                // Skip this cell
                                continue;
                            } */

            Point voronoi_vertex = dt.dual(ch);
            vc.voronoi_vertices.push_back(voronoi_vertex);

            // Collect the scalar value
            float value = vertexValueMap[voronoi_vertex];
            vc.vertex_values.push_back(value);
        }

        // Build the convex hull and then extract facets from the polyhedron
        CGAL::convex_hull_3(vc.voronoi_vertices.begin(), vc.voronoi_vertices.end(), vc.polyhedron); // A function computes convex hull of a set of 3D points, stores in vc.polyhedron

        // Extract facets from polyhedrons
        for (auto facet_it = vc.polyhedron.facets_begin(); facet_it != vc.polyhedron.facets_end(); ++facet_it)
        {
            std::vector<Point> facet_vertices;
            std::vector<float> facet_values;

            // Get the halfedge around the facet
            auto h = facet_it->facet_begin(); // Returns a halfedge iterator to traverse around the facet
            do
            {
                // Get the point corresponding to the vertex of this halfedge
                Point p = h->vertex()->point();
                facet_vertices.push_back(p);

                // Retrieve the scalar value associated with this point
                float value = vertexValueMap[p];
                facet_values.push_back(value);

                ++h;
            } while (h != facet_it->facet_begin()); // Loop until we complete the cycle

            // Extract and create a Facet object and add it to the VoronoiCell
            vc.facets.emplace_back(facet_vertices, facet_values);
        }

        // Add the VoronoiCell to the collection
        voronoi_cells.push_back(vc);
    }
}

void compute_voronoi_values(std::vector<Point> &voronoi_vertices, ScalarGrid &grid, std::vector<float> &voronoi_vertex_values)
{
    for (const auto &vertex : voronoi_vertices)
    {
        float value = trilinear_interpolate(vertex, grid);
        voronoi_vertex_values.push_back(value);
        vertexValueMap[vertex] = value;
    }
}

void construct_voronoi_edges(std::map<Object, std::vector<Facet>, ObjectComparator> &delaunay_facet_to_voronoi_edge_map, std::set<std::string> &seen_edges, std::vector<CGAL::Object> &voronoi_edges)
{
    for (Delaunay::Finite_facets_iterator fit = dt.finite_facets_begin(); fit != dt.finite_facets_end(); ++fit)
    {
        Facet facet = *fit;
        Segment3 seg;
        Ray3 ray;
        Line3 line;

        Delaunay::Object vEdge = dt.dual(facet);

        if (isDegenerate(vEdge))
        {
            continue;
        }
        std::string edgeRep = objectToString(vEdge);

        delaunay_facet_to_voronoi_edge_map[vEdge].push_back(facet);

        if (seen_edges.find(edgeRep) == seen_edges.end())
        {
            voronoi_edges.push_back(vEdge);
            seen_edges.insert(edgeRep);
        }
    }
}

int handle_output_mesh(bool &retFlag)
{
    retFlag = true;
    // Use locations of isosurface vertices as vertices of Delaunay triangles and write the output mesh
    if (multi_isov)
    {
        if (output_format == "off")
        {
            writeOFFMulti(output_filename, isosurfaceVertices, isoTriangles);
        }
        else if (output_format == "ply")
        {
            writePLYMulti(output_filename, isosurfaceVertices, isoTriangles);
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
            writeOFFSingle(output_filename, isosurfaceVertices, dualTriangles, point_index_map);
        }
        else if (output_format == "ply")
        {
            writePLYSingle(output_filename, isosurfaceVertices, dualTriangles, point_index_map);
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