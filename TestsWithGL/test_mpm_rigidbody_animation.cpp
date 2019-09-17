#include "TestsWithGL_pcp.h"

#include "test_post_processor.h"

#include "GenerateAnimation.h"

void test_mpm_rigid_animation(void)
{
	GenerateAnimation gen;
	// square
	gen.generate(5.0, -2.0, 2.0, -0.5, 3.5, "mpm_rb_res_square.bin", "mpm_rb_res_square.gif");
	// cantilever beam
	//gen.generate(5.0, -1.0, 6.0, -0.5, 6.5, "mpm_rb_res_can.bin", "mpm_rb_res_can.gif");
}
