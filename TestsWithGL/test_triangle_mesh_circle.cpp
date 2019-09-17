#include "TestsWithGL_pcp.h"

#include "test_sim_core.h"
#include "TriangleMesh.h"

void test_triangle_mesh_circle(void)
{
	TriangleMesh mesh;
	mesh.load_mesh("..\\..\\Asset\\circle_mesh.mesh_data");
	mesh.get_bounding_box();
	mesh.find_edges_at_boundary();
	mesh.init_bg_grid();
#ifdef __DEBUG_TRIANGLE_MESH__
	//mesh.dump_mesh_info("square_mesh.txt", true);
#endif
	
	//display_triangle_mesh(mesh, true, true, true);

	Point pt1 = { 0.5, 0.2 };
	double dist, nx, ny;
	mesh.distance_to_boundary(pt1, dist, nx, ny);
#ifdef __DEBUG_TRIANGLE_MESH__
	display_triangle_mesh(mesh, true, true, true, mesh.closest_edge, &pt1);
#endif

	// point to edge distance
	//TriangleMesh::Edge e2 = { 0, 4 };
	//Point pt2 = { 0.5, 0.2 };
	//double dist2, nx2, ny2;
	//mesh.distance_to_edge(e2, pt2, dist2, nx2, ny2);
	//std::cout << "dist: " << dist2 << " normal: (" << nx2 << "," << ny2 << ")\n";
	//system("pause");
}
