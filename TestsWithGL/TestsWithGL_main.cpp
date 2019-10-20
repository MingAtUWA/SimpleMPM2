#include "TestsWithGL_pcp.h"

#include <cstdlib>

#include "test_sim_core.h"
#include "test_post_processor.h"

int main(int argc, void **argv)
{
	//test_solve_functions();

	//test_triangle_mesh_circle();
	//test_triangle_mesh_square();
	//test_rigid_body_square();
	//test_get_pcls_from_mesh();

	//test_matrix_coefficient_set();
	//test_cal_stiffness_mat();

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
	
	// GIMP
	//test_init_pcl_gimp();
	// 1D consolidation
	//test_mpm_chm_s_1d_consolidation();
	//test_animation_chm_s_1d_consolidation();
	
	// Implicit 1d bar compression
	//test_mpm_me_s_up_1dbar();
	//test_animation_me_s_up_1dbar();

	test_mpm_me_s_fem_up_1dbar();

	//system("pause");
	return 0;
}