#include "TestsWithGL_pcp.h"

#include <fstream>

#include "Model_S2D_ME_s_RigidBody.h"
#include "Step_S2D_ME_s_RigidBody.h"
#include "TimeHistory_S2D_ME_s_RigidBody.h"
#include "ResultFile_PlainBin.h"

#include "test_sim_core.h"

void test_mpm_rigidbody_circle(void)
{
	Model_S2D_ME_s_RigidBody model;

	Step_S2D_ME_s_RigidBody step;
	step.set_model(model);
	step.set_dtime(0.01);

	ResultFile_PlainBin res_file;
	res_file.init("mpm_rb_res.bin");
	
	res_file.output(model);

	TimeHistory_S2D_ME_s_RigidBody out1;
	out1.set_res_file(res_file);
	out1.set_interval_num(10);

	step.add_time_history(out1);
	step.solve();
}
