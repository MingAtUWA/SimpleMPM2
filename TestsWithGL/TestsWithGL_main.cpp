#include "TestsWithGL_pcp.h"

#include <cstdlib>

#include "test_sim_core.h"
#include "test_post_processor.h"

int main(int argc, void **argv)
{
	//test_triangle_mesh_circle();
	//test_triangle_mesh_square();
	//test_rigid_body_square();
	test_get_pcls_from_mesh();

	// square
	//test_mpm_rigidbody_square();
	//test_mpm_rigid_animation_square();

	// cantilever beam
	//test_mpm_rigidbody_cantilever();
	//test_mpm_rigid_animation_can();

	// cantilever beam - frictional contact
	//test_mpm_rigidbody_cantilever_fric();
	//test_mpm_rigid_animation_can_fric();

	//system("pause");
	return 0;
}