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
	
	// mixed u-p implicit FEM
	//test_fem_me_s_up_1dbar();
	//test_fem_chm_s_1d_consolidation();
	
	// Implicit ME MPM
	//test_imp_mpm_me_s_up_1dbar();
	//test_animation_me_s_up_1dbar();

	// Implicit CHM MPM
	//test_imp_mpm_chm_s_uup_1d_consolidation();
	//test_animation_chm_s_uup_1d_consolidation();

	// trianglur mesh mpm me
	//test_t2d_mpm_me_s_1d_compression();
	//test_color_animation_t2d_me_s_1d_compression();
	
	//test_t2d_mpm_me_s_t_bar_coarser();
	//test_color_animation_t2d_me_s_t_bar_coarser();

	// t bar penetrate from above ground
	//test_t2d_mpm_me_s_t_bar_above_ground();
	//test_color_animation_t2d_me_s_t_bar_above_ground();

	//test_t2d_mpm_chm_s_t_bar_above_ground();
	//test_color_animation_t2d_chm_s_t_bar_above_ground();

	// geostatic
	//test_t2d_mpm_me_s_geostatic();
	
	//test_t2d_mpm_chm_s_geostatic();
	//test_color_animation_t2d_mpm_chm_s_geostatic();

	// trianglur mesh mpm chm
	//test_t2d_mpm_square();
	//test_t2d_mpm_chm_s_1d_consolidation();
	//test_animation_t2d_chm_s_1d_consolidation();
	//test_color_animation_t2d_chm_s_1d_consolidation();

	//test_t2d_mpm_chm_s_t_bar();
	//test_animation_t2d_chm_s_t_bar();
	//test_color_animation_t2d_chm_s_t_bar();

	//test_t2d_mpm_chm_s_t_bar_coarser();
	//test_color_animation_t2d_chm_s_t_bar_coarser();

	//test_triangle_searching();

	//test_color_graph();

	// a realistic case of t-bar penetration
	test_t2d_mpm_chm_s_t_bar_real();
	test_color_animation_t2d_chm_s_t_bar_real();

	//system("pause");
	return 0;
}