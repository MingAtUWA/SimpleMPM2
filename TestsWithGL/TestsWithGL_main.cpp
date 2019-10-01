#include "TestsWithGL_pcp.h"

#include <cstdlib>

#include "test_sim_core.h"
#include "test_post_processor.h"

int main(int argc, void **argv)
{
	//test_triangle_mesh_circle();
	//test_triangle_mesh_square();
	//test_rigid_body_square();
	//test_get_pcls_from_mesh();

	// square
	//test_mpm_rigidbody_square();
	//test_mpm_rigid_animation_square();

	// cantilever beam
	//test_mpm_rigidbody_cantilever();
	//test_mpm_rigid_animation_can();

	// cantilever beam - frictional contact
	//test_mpm_rigidbody_cantilever_fric();
	//test_mpm_rigid_animation_can_fric();

	// block slide down slope
	//test_slide_down_frictional_slope();
	// 2
	//test_slide_down_frictional_slope2();
	//test_mpm_rigid_animation_slope_fric();

	// bar compression
	//test_mpm_rigidbody_bar_compression();
	//test_mpm_rigid_animation_bar_compression();

	// 1D consolidation
	test_mpm_chm_s_1d_consolidation();
	//test_animation_chm_s_1d_consolidation();

	//system("pause");
	return 0;
}