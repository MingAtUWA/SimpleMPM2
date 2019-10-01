#include "TestsWithGL_pcp.h"

#include "test_sim_core.h"
#include "RigidBody.h"

void test_rigid_body_square(void)
{
	RigidBody rb;
	//rb.load_and_init_mesh("..\\..\\Asset\\square_mesh.mesh_data", 0.3);
	rb.load_and_init_mesh("..\\..\\Asset\\cap_mesh.mesh_data", 0.25);
	//display_triangle_mesh(rb.mesh, true, true, true);

	rb.x = 0.0;
	rb.y = 0.0;
	rb.theta = 3.1415 / 4.0;
	rb.init_transformation();

	Point pt1 = { 0.6, 0.2 };
	double dist, nx, ny;
	rb.distance_from_boundary(pt1, dist, nx, ny);
	Point lp1 = rb.to_local_coord(pt1);
	Point pt2 = rb.to_global_coord(lp1);
#ifdef __DEBUG_TRIANGLE_MESH__
	display_triangle_mesh(rb.mesh, true, true, true, rb.mesh.closest_edge, &lp1);
#endif
	//rb.set_params(1.0);
	//rb.add_ext_force(1.0, 0.0, 0.0, 1.0);
	//rb.init_calculation();
	//rb.predict_motion_from_ext_force(0.1);
	//std::cout << "x: " << rb.x << ", y: " << rb.y
	//		  << ", theta: " << rb.theta << "\n";

	//system("pause");
}
