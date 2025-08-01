#ifndef VDC_FUNC_H
#define VDC_FUNC_H

#include "vdc_type.h"
#include "vdc_utilities.h"
#include "vdc_io.h"
#include "vdc_commandline.h"
#include "vdc_voronoi.h"

//! @brief Computes the dual triangles for the final mesh in the single isovertex case.
/*!
 * This function calculates the Delaunay triangles dual to bipolar edges in
 * the Voronoi diagram for a single isovalue.
 *
 * @param iso_surface Instance of IsoSurface to store triangles.
 * @param vd Voronoi diagram containing edges and vertices
 * @param bbox Bounding box of the computational domain.
 * @param dt Delaunay triangulation structure.
 * @param grid Scalar grid containing scalar values.
 * @param isovalue The isovalue used for computing.
 * @param pointToIndexMap Map between points and their indices in Delaunay triangulation
 */
void compute_dual_triangles(
    IsoSurface &iso_surface,
    VoronoiDiagram &vd,
    CGAL::Epick::Iso_cuboid_3 &bbox,
    Delaunay &dt,
    UnifiedGrid &grid,
    float isovalue,
    std::map<Point, int> &pointToIndexMap);

    
//! @brief Computes the dual triangles for the final mesh in the multi-isovertex case.
/*!
 * This function calculates the Delaunay triangles dual to bipolar edges in
 * the Voronoi diagram for multiple isovalues.
 *
 * @param voronoiDiagram The Voronoi diagram to compute from.
 * @param bbox Bounding box of the computational domain.
 * @param grid Scalar grid containing scalar values.
 * @param isovalue The isovalue for mesh computation.
 * @param iso_surface Instance of IsoSurface containing the isosurface vertices and faces.
 */
void compute_dual_triangles_multi(
    VoronoiDiagram &voronoiDiagram,
    CGAL::Epick::Iso_cuboid_3 &bbox,
    UnifiedGrid &grid,
    float isovalue,
    IsoSurface &iso_surface);

//! @brief Computes isosurface vertices for the multi-isovertex case.
/*!
 * Processes Voronoi cells to identify bipolar edges, form cycles, and compute
 * centroids as isosurface vertices.
 *
 * @param voronoiDiagram The Voronoi diagram to compute vertices for.
 * @param isovalue The isovalue to use for vertex computation.
 * @param iso_surface Instance of IsoSurface containing the isosurface vertices and faces.
 */
void Compute_Isosurface_Vertices_Multi(VoronoiDiagram &voronoiDiagram, float isovalue, IsoSurface &iso_surface);

//! @brief Computes isosurface vertices for the single-isovertex case.
/*!
 * Computes isosurface vertices by interpolating along edges of active cubes
 * and calculating their centroids.
 *
 * @param grid The scalar grid containing scalar values.
 * @param isovalue The isovalue to use for calculation.
 * @param iso_surface Instance of IsoSurface containing the isosurface vertices and faces.
 * @param data_grid The grid containing input data.
 * @param activeCubeCenters The list of center points of active cubes.
 */
void compute_isosurface_vertices_single(UnifiedGrid &grid, float isovalue, IsoSurface &iso_surface, std::vector<Point> &activeCubeCenters);

//! @brief Constructs a Delaunay triangulation from a grid and grid facets.
/*!
 * Constructs a 3D Delaunay triangulation using the grid's scalar values
 * and the facets of the active cubes.
 *
 * @param dt The Delaunay triangulation instance.
 * @param grid The grid containing scalar values.
 * @param grid_facets The grid facets to use in constructing the triangulation.
 * @param vdc_param The VDC_PARAM instance holding user input options.
 * @param activeCubeCenters The list of center points of active cubes.
 */
void construct_delaunay_triangulation(Delaunay &dt, UnifiedGrid &grid, const std::vector<std::vector<GRID_FACETS>> &grid_facets, VDC_PARAM &vdc_param, std::vector<Point> &activeCubeCenters);

//! @brief Adds dummy points from a facet for Voronoi diagram bounding.
/*!
 * Adds dummy points to ensure the Voronoi diagram is bounded within the computational domain.
 *
 * @param facet The grid facet to extract dummy points from.
 * @param grid The grid containing data for point computation.
 * @return A vector of points added from the facet.
 */
std::vector<Point> add_dummy_from_facet(const GRID_FACETS &facet, const UnifiedGrid &grid);

//! @brief Constructs Voronoi vertices for the given diagram.
/*!
 * Generates the vertices of the Voronoi diagram based on the Delaunay triangulation.
 *
 * @param voronoiDiagram The Voronoi diagram to construct vertices for.
 * @param dt The Delaunay triangulation corresponding (dual) to the Voronoi diagram.
 */
void construct_voronoi_vertices(VoronoiDiagram &voronoiDiagram, Delaunay &dt);

//! @brief Constructs Voronoi cells from the Delaunay triangulation.
/*!
 * Populates the Voronoi diagram with polyhedral cells derived from the Delaunay triangulation.
 *
 * @param voronoiDiagram The Voronoi diagram to populate with cells.
 * @param dt The Delaunay triangulation corresponding (dual) to the Voronoi diagram.
 */
void construct_voronoi_cells_as_convex_hull(VoronoiDiagram &voronoiDiagram, Delaunay &dt);

//! @brief Constructs Voronoi cells without using Convex_Hull_3 (in development).
/*!
 * Populates the Voronoi diagram with polyhedral cells derived from the Delaunay
 * triangulation by processing incident edges and cell circulators.
 *
 * @param voronoiDiagram The Voronoi diagram to populate with cells.
 * @param dt The Delaunay triangulation corresponding (dual) to the Voronoi diagram.
 */
void construct_voronoi_cells_from_delaunay_triangulation(VoronoiDiagram &voronoiDiagram, Delaunay &dt);

//! @brief Validates that each Facet in the Voronoi Diagram has the correct orientation and normal vector.
void validate_facet_orientations_and_normals(VoronoiDiagram &voronoiDiagram);


//! @brief Computes Voronoi vertex values using scalar grid interpolation.
/*!
 * Interpolates scalar values from the grid to each vertex of the Voronoi diagram.
 *
 * @param voronoiDiagram The Voronoi diagram to compute values for.
 * @param grid The scalar grid containing data.
 */
void compute_voronoi_values(VoronoiDiagram &voronoiDiagram, UnifiedGrid &grid);

//! @brief Constructs Voronoi edges from Delaunay facets.
/*!
 * Derives the edges of the Voronoi diagram by processing the facets of the Delaunay triangulation.
 *
 * @param voronoiDiagram The Voronoi diagram to populate with edges.
 * @param dt The Delaunay triangulation corresponding (dual) to the Voronoi diagram.
 */
void construct_voronoi_edges(
    VoronoiDiagram &voronoiDiagram,
    Delaunay &dt);

//! @brief Constructs and links Voronoi cell edges in the Voronoi diagram.
/*!
 * Builds cell edges for each Voronoi edge, links them in a circular ring, and
 * updates edge mappings.
 *
 * @param voronoiDiagram The Voronoi diagram to populate with edges.
 * @param bbox The bounding box used for clipping rays and lines.
 * @param dt The Delaunay triangulation used for edge-facet correspondence.
 */
void construct_voronoi_cell_edges(VoronoiDiagram &voronoiDiagram,
    CGAL::Epick::Iso_cuboid_3 &bbox,
    Delaunay &dt);

//! @brief Handles output mesh generation.
/*!
 * Generates the final output mesh based on the Voronoi diagram and handles file output.
 *
 * @param retFlag Reference to a flag indicating success or failure.
 * @param vd The Voronoi diagram containing mesh data.
 * @param vdc_param The VDC_PARAM instance containing user input options.
 * @param iso_surface The instance of IsoSurface containing the vertices and faces of the isosurface.
 * @return An integer representing the exit status.
 */
int handle_output_mesh(bool &retFlag, VoronoiDiagram &vd, VDC_PARAM &vdc_param, IsoSurface &iso_surface);

//! @brief Wraps up the process of building the Voronoi diagram from the Delaunay triangulation.
/*!
 * Orchestrates the construction of Voronoi vertices, edges, cells, and values.
 *
 * @param vd The Voronoi diagram to construct.
 * @param vdc_param The VDC_PARAM instance containing user input options.
 * @param grid Scalar grid containing scalar values.
 * @param bbox Bounding box of the computational domain.
 * @param dt The Delaunay triangulation.
 */
void construct_voronoi_diagram(VoronoiDiagram &vd, VDC_PARAM &vdc_param, UnifiedGrid &grid, CGAL::Epick::Iso_cuboid_3 &bbox, Delaunay &dt);

//! @brief Wraps up the process of building the isosurface from the Voronoi diagram/Delaunay triangulation.
/*!
 * Constructs isosurface vertices and triangles based on single or multi-isovertex mode.
 *
 * @param dt Delaunay triangulation used for dual computation
 * @param vd Voronoi diagram containing cell and edge data
 * @param vdc_param Processing parameters and options
 * @param iso_surface Output isosurface to store vertices and triangles
 * @param grid Scalar grid for value interpolation
 * @param data_grid Input data grid for single-isovertex mode
 * @param activeCubeCenters Active cube centers for single-isovertex mode
 * @param bbox Bounding box for clipping infinite edges
 */
void construct_iso_surface(Delaunay &dt, VoronoiDiagram &vd, VDC_PARAM &vdc_param, IsoSurface &iso_surface, UnifiedGrid &grid, std::vector<Point> &activeCubeCenters, CGAL::Epick::Iso_cuboid_3 &bbox);


// Helper function declarations (internal linkage)
//! @brief Generates a Delaunay triangle based on orientation and cell finiteness.
/*!
 * Adds a triangle to the dualTriangles vector, adjusting vertex order based on
 * the orientation and whether the cell is infinite.
 *
 * @param p1 First vertex of the triangle.
 * @param p2 Second vertex of the triangle.
 * @param p3 Third vertex of the triangle.
 * @param iOrient Orientation value determining vertex order.
 * @param isInfinite Flag indicating if the cell is infinite.
 * @param dualTriangles Vector to store the generated triangles.
 */
static void generate_triangle(const Point &p1, const Point &p2, const Point &p3, int iOrient, bool isInfinite, std::vector<DelaunayTriangle> &dualTriangles);

//! @brief Processes a segment edge for dual triangle computation.
/*!
 * Checks if the segment is bipolar and generates triangles for each associated
 * Delaunay facet.
 *
 * @param edge Voronoi edge to process
 * @param vd Voronoi diagram containing edge data
 * @param isovalue The isovalue for bipolarity check
 * @param dt Delaunay triangulation for facet lookup
 * @param dualTriangles Output vector for generated triangles
 */
static void process_segment_edge(
    VoronoiEdge &edge,
    VoronoiDiagram &vd,
    float isovalue,
    Delaunay &dt,
    std::vector<DelaunayTriangle> &dualTriangles);

//! @brief Processes a ray edge for dual triangle computation.
/*!
 * Intersects the ray with the bounding box, checks bipolarity, and generates
 * triangles for associated Delaunay facets.
 *
 * @param edge The CGAL object representing the edge.
 * @param bbox The bounding box for intersection.
 * @param grid The scalar grid for interpolation.
 * @param isovalue The isovalue for bipolarity check.
 * @param dt The Delaunay triangulation.
 * @param dualTriangles Vector to store generated triangles.
 */
static void process_ray_edge(
    VoronoiEdge &edge,
    VoronoiDiagram &vd,
    CGAL::Epick::Iso_cuboid_3 &bbox,
    UnifiedGrid &grid,
    float isovalue,
    Delaunay &dt,
    std::vector<DelaunayTriangle> &dualTriangles);

//! @brief Processes a line edge for dual triangle computation.
/*!
 * Intersects the line with the bounding box, checks bipolarity, and generates
 * triangles for associated Delaunay facets.
 *
 * @param line The line edge to process.
 * @param edge The CGAL object representing the edge.
 * @param grid The scalar grid for interpolation.
 * @param isovalue The isovalue for bipolarity check.
 * @param bbox The bounding box for intersection.
 * @param dt The Delaunay triangulation.
 * @param dualTriangles Vector to store generated triangles.
 */
static void process_line_edge(
    const Line3 &line,
    VoronoiEdge &edge,
    UnifiedGrid &grid,
    float isovalue,
    CGAL::Epick::Iso_cuboid_3 &bbox,
    Delaunay &dt,
    std::vector<DelaunayTriangle> &dualTriangles);

//! @brief Generates a triangle for the multi-isovertex case.
/*!
 * Adds a triangle to the isosurface's triangle list, adjusting vertex order based
 * on orientation, and logs problematic triangles if invalid.
 *
 * @param iso_surface The isosurface to store the triangle.
 * @param idx1 First vertex index.
 * @param idx2 Second vertex index.
 * @param idx3 Third vertex index.
 * @param iOrient Orientation value determining vertex order.
 * @param isValid Flag indicating if the triangle is valid.
 */
static void generate_triangle_multi(IsoSurface &iso_surface, int idx1, int idx2, int idx3, int iOrient, bool isValid);

//! @brief Selects isovertices for a Delaunay facet.
/*!
 * Retrieves isovertex indices for the three vertices of a facet, ensuring they
 * are valid and not dummy vertices.
 *
 * @param voronoiDiagram The Voronoi diagram containing cell data.
 * @param facet The Delaunay facet to process.
 * @param globalEdgeIndex The global edge index for isovertex selection.
 * @param idx1 Output parameter for the first vertex index.
 * @param idx2 Output parameter for the second vertex index.
 * @param idx3 Output parameter for the third vertex index.
 * @param cellIndex1 Output parameter for the first cell index.
 * @param cellIndex2 Output parameter for the second cell index.
 * @param cellIndex3 Output parameter for the third cell index.
 * @return True if the vertices are valid, false otherwise.
 */
static bool select_isovertices(const VoronoiDiagram &voronoiDiagram, const Facet &facet, int globalEdgeIndex, int &idx1, int &idx2, int &idx3, int &cellIndex1, int &cellIndex2, int &cellIndex3);

//! @brief Processes a segment edge for multi-isovertex triangle computation.
/*!
 * Checks if the segment is bipolar and generates triangles for associated
 * Delaunay facets in multi-isovertex mode.
 *
 * @param edge Voronoi edge to process
 * @param voronoiDiagram Voronoi diagram containing edge and cell data
 * @param isovalue The isovalue for bipolarity check
 * @param iso_surface The isosurface to store triangles
 */
static void process_segment_edge_multi(
    VoronoiEdge edge,
    VoronoiDiagram &voronoiDiagram,
    float isovalue,
    IsoSurface &iso_surface);

//! @brief Processes a ray edge for multi-isovertex triangle computation.
/*!
 * Intersects the ray with the bounding box, checks bipolarity, and generates
 * triangles using the first isovertex from each cell.
 *
 * @param source_pt The index of the source point of the ray edge in VoronoiDiagram
 * @param ray The ray edge to process
 * @param dualDelaunayFacets List of Delaunay facets dual to this edge
 * @param voronoiDiagram The Voronoi diagram containing edge and cell data
 * @param grid The scalar grid for interpolation
 * @param isovalue The isovalue for bipolarity check
 * @param bbox The bounding box for intersection
 * @param iso_surface The isosurface to store triangles
 */
static void process_ray_edge_multi(
    int source_pt,
    const Ray3 &ray,
    std::vector<Facet> dualDelaunayFacets,
    VoronoiDiagram &voronoiDiagram,
    UnifiedGrid &grid,
    float isovalue,
    CGAL::Epick::Iso_cuboid_3 &bbox,
    IsoSurface &iso_surface);


//! @brief Processes a line edge for multi-isovertex triangle computation.
/*!
 * Intersects the line with the bounding box, checks bipolarity, and generates
 * triangles using the first isovertex from each cell.
 *
 * @param line The line edge to process
 * @param dualDelaunayFacets List of Delaunay facets dual to this edge
 * @param voronoiDiagram The Voronoi diagram containing edge and cell data
 * @param grid The scalar grid for interpolation
 * @param isovalue The isovalue for bipolarity check
 * @param bbox The bounding box for intersection
 * @param iso_surface The isosurface to store triangles
 */
static void process_line_edge_multi(
    const Line3 &line,
    std::vector<Facet> dualDelaunayFacets,
    VoronoiDiagram &voronoiDiagram,
    UnifiedGrid &grid,
    float isovalue,
    CGAL::Epick::Iso_cuboid_3 &bbox,
    IsoSurface &iso_surface);

//! @brief Collects midpoints for bipolar edges in a Voronoi cell's facets.
/*!
 * Identifies bipolar edges, computes their midpoints, and stores them along with
 * facet associations.
 *
 * @param vc The Voronoi cell to process.
 * @param voronoiDiagram The Voronoi diagram containing facet and vertex data.
 * @param isovalue The isovalue for bipolarity check.
 * @param midpoints Vector to store computed midpoints.
 * @param edge_to_midpoint_index Map linking edge keys to midpoint indices.
 * @param facet_midpoint_indices Vector storing midpoint indices for each facet.
 */
static void collect_midpoints(VoronoiCell &vc, VoronoiDiagram &voronoiDiagram, float isovalue, std::vector<MidpointNode> &midpoints, std::map<std::pair<int, int>, int> &edge_to_midpoint_index, std::vector<std::vector<int>> &facet_midpoint_indices);

//! @brief Connects midpoints within each facet to form a graph.
/*!
 * Links pairs of midpoints in each facet to establish connectivity for cycle detection.
 *
 * @param facet_midpoint_indices Vector storing midpoint indices for each facet.
 * @param midpoints Vector of midpoints to update with connectivity.
 */
static void connect_midpoints(const std::vector<std::vector<int>> &facet_midpoint_indices, std::vector<MidpointNode> &midpoints);

//! @brief Extracts cycles from the midpoint connectivity graph.
/*!
 * Uses depth-first search to identify closed cycles in the midpoint graph.
 *
 * @param midpoints Vector of midpoints with connectivity information.
 * @param cycles Vector to store the extracted cycles as lists of midpoint indices.
 */
static void extract_cycles(const std::vector<MidpointNode> &midpoints, std::vector<std::vector<int>> &cycles);

//! @brief Computes centroids for cycles and updates the isosurface.
/*!
 * Calculates the centroid for each cycle, updates the Voronoi cell's cycles,
 * and adds the centroid to the isosurface vertices.
 *
 * @param vc The Voronoi cell to update.
 * @param voronoiDiagram The Voronoi diagram for edge lookup.
 * @param midpoints Vector of midpoints used for centroid computation.
 * @param cycles Vector of cycles as lists of midpoint indices.
 * @param iso_surface The isosurface to store vertices.
 */
static void compute_cycle_centroids(VoronoiCell &vc, VoronoiDiagram &voronoiDiagram, std::vector<MidpointNode> &midpoints, const std::vector<std::vector<int>> &cycles, IsoSurface &iso_surface);

//! @brief Builds Voronoi cell edges for each edge in the diagram.
/*!
 * Creates VoronoiCellEdge entries for cells sharing each edge, collecting cell indices
 * from associated Delaunay facets.
 *
 * @param voronoiDiagram The Voronoi diagram to populate with cell edges.
 * @param dt The Delaunay triangulation.
 */
static void build_cell_edges(VoronoiDiagram &voronoiDiagram, Delaunay &dt);

//! @brief Links Voronoi cell edges in a circular ring.
/*!
 * Connects cell edges sharing the same edge index using the nextCellEdge field
 * to form a closed loop.
 *
 * @param voronoiDiagram The Voronoi diagram containing cell edges to link.
 */
static void link_cell_edges(VoronoiDiagram &voronoiDiagram);

//! @brief Processes edge mapping for a single Voronoi edge.
/*!
 * Updates the segmentVertexPairToEdgeIndex map for segments, rays, and lines
 * after intersecting with the bounding box.
 *
 * @param voronoiDiagram The Voronoi diagram to update.
 * @param edge The VoronoiEdge representing the edge.
 * @param edgeIdx The index of the edge in the diagram.
 * @param bbox The bounding box for intersection.
 */
static void process_edge_mapping(VoronoiDiagram &voronoiDiagram, VoronoiEdge &edge, int edgeIdx, CGAL::Epick::Iso_cuboid_3 &bbox);

//! @brief Updates edge mappings for all Voronoi edges.
/*!
 * Processes all edges to update segmentVertexPairToEdgeIndex and cellEdgeLookup maps.
 *
 * @param voronoiDiagram The Voronoi diagram to update.
 * @param bbox The bounding box for intersection.
 */
static void update_edge_mapping(VoronoiDiagram &voronoiDiagram, CGAL::Epick::Iso_cuboid_3 &bbox);

//! @brief Creates a Voronoi cell for a Delaunay vertex.
/*!
 * Initializes a Voronoi cell with the given cell index and Delaunay vertex handle.
 *
 * @param delaunay_vertex The Delaunay vertex to create the cell for.
 * @param cellIndex The index to assign to the cell.
 * @return The initialized Voronoi cell.
 */
static VoronoiCell create_voronoi_cell(Vertex_handle delaunay_vertex, int cellIndex);

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
static void collect_cell_vertices(Delaunay &dt, Vertex_handle delaunay_vertex, VoronoiDiagram &voronoiDiagram, std::vector<int> &vertices_indices);

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
 * @param edge_to_facets Map tracking which facets share each edge.
 * @return The constructed Voronoi facet, or an empty facet if invalid.
 */
static VoronoiCellFacet build_facet_from_edge(Delaunay &dt, const Edge &ed, Vertex_handle delaunay_vertex, VoronoiDiagram &voronoiDiagram, std::vector<int> &facet_indices, std::map<std::pair<int, int>, std::vector<int>>& edge_to_facets);

//! @brief Processes incident edges to build facets for a Voronoi cell.
/*!
 * Iterates over incident edges to construct facets and add them to the cell.
 *
 * @param dt The Delaunay triangulation.
 * @param delaunay_vertex The Delaunay vertex to process.
 * @param voronoiDiagram The Voronoi diagram to update.
 * @param vc The Voronoi cell to populate with facets.
 * @param edge_to_facets Map tracking which facets share each edge.
 */
static void process_incident_edges(Delaunay &dt, Vertex_handle delaunay_vertex, VoronoiDiagram &voronoiDiagram, VoronoiCell &vc, std::map<std::pair<int, int>, std::vector<int>>& edge_to_facets);
//! @brief Collects points for the Delaunay triangulation.
/*!
 * Gathers original points and dummy points from grid facets for multi-isovertex mode,
 * or uses only active cube centers for single-isovertex mode.
 *
 * @param grid The grid containing data.
 * @param grid_facets The grid facets for dummy point generation.
 * @param activeCubeCenters The list of center points of active cubes.
 * @param vdc_param The VDC_PARAM instance containing user input options.
 * @param delaunay_points Output vector for all points (original + dummy).
 */
static int collect_delaunay_points(UnifiedGrid &grid,
                                   const std::vector<std::vector<GRID_FACETS>> &grid_facets,
                                   const std::vector<Point> &activeCubeCenters,
                                   VDC_PARAM &vdc_param,
                                   std::vector<Point> &delaunay_points);


//! @brief Inserts a point into Delaunay triangulation.
/*!
 * Handles point insertion with proper indexing and dummy point marking.
 *
 * @param dt Delaunay triangulation
 * @param p Point to insert
 * @param index Point index
 * @param is_dummy Whether this is a dummy point
 * @return Handle to inserted vertex
 */
static Vertex_handle insert_point_into_delaunay_triangulation(Delaunay &dt,
                                                  const Point &p,
                                                  int index,
                                                  bool is_dummy);

//! @brief Selects an isovertex from a Voronoi cell edge.
/*!
 * Retrieves the isovertex index associated with a cell edge, traversing linked edges if necessary.
 *
 * @param voronoiDiagram The Voronoi diagram containing cell edge data.
 * @param cellIndex The index of the Voronoi cell.
 * @param globalEdgeIndex The global edge index to select the isovertex for.
 * @return The isovertex index, or -1 if not found.
 */
static inline int select_isovertex_from_cell_edge(const VoronoiDiagram &voronoiDiagram, int cellIndex, int globalEdgeIndex);

//! @brief Orders facet vertices in cyclic order around a Delaunay edge.
/*!
 * Sorts vertex indices based on their angular position around the Delaunay edge
 * defined by points p0 and p1.
 *
 * @param indices The vertex indices to sort.
 * @param p0 The first point of the Delaunay edge.
 * @param p1 The second point of the Delaunay edge.
 * @param vertices The Voronoi vertices containing position data.
 */
static void order_facet_vertices(std::vector<int> &indices, const Point &p0, const Point &p1, const std::vector<VoronoiVertex> &vertices);

#endif