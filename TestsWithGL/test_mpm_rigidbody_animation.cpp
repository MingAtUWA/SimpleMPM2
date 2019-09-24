#include "TestsWithGL_pcp.h"

#include "test_post_processor.h"

#include "GA_S2D_ME_s_RigidBody.h"
#include "GA_S2D_ME_s_RigidBody_Fric.h"

void test_mpm_rigid_animation_square(void)
{
	GA_S2D_ME_s_RigidBody gen;
	// square
	gen.generate(5.0, -2.0, 2.0, -0.5, 3.5, "mpm_rb_res_square.bin", "mpm_rb_res_square.gif");
}

void test_mpm_rigid_animation_can(void)
{
	GA_S2D_ME_s_RigidBody gen;
	// cantilever beam
	gen.generate(5.0, -1.0, 6.0, -0.5, 6.5, "mpm_rb_res_can.bin", "mpm_rb_res_can.gif");
}

void test_mpm_rigid_animation_can_fric(void)
{
	GA_S2D_ME_s_RigidBody_Fric gen;
	// cantilever beam
	gen.generate(5.0, -1.0, 6.0, -0.5, 6.5, "mpm_rb_res_can_fric.bin", "mpm_rb_res_can_fric.gif");
}

void test_mpm_rigid_animation_slope_fric(void)
{
	GA_S2D_ME_s_RigidBody_Fric gen;
	// cantilever beam
	gen.generate(5.0, -1.0, 9.5, -1.0, 7.0, "mpm_rb_res_slope_fric.bin", "mpm_rb_res_slope_fric.gif");
}