#include "TestsWithGL_pcp.h"

#include "test_sim_core.h"

#include "Model_S2D_ME_s_RigidBody_Fric.h"
#include "Step_S2D_ME_s_RigidBody_Fric.h"
#include "ModelDataOutput_S2D_ME_s_RigidBody_Fric.h"
#include "TimeHistoryOutput_ConsoleProgressBar.h"
#include "TimeHistoryOutput_S2D_ME_s_RigidBody_Fric.h"

#include "ResultFile_PlainBin.h"
#include "ResultFile_XML.h"

#include "DisplayModel.h"

void test_mpm_rigidbody_bar_compression(void)
{
	Model_S2D_ME_s_RigidBody_Fric model;

	// background mesh
	double bgm_h = 0.5;
	size_t bgm_elem_x_num = 4, bgm_elem_y_num = 11;
	model.init_mesh(bgm_h, bgm_elem_x_num, bgm_elem_y_num, -0.5, 0.0);

	// material points
	model.init_pcl(80, 0.5 * 0.5 / 4.0, 1.0, 100.0, 0.3);
	size_t pcl_id = 0;
	for (size_t y_id = 0; y_id < 20; ++y_id)
		for (size_t x_id = 0; x_id < 4; ++x_id)
		{
			Model_S2D_ME_s_RigidBody_Fric::Particle &pcl = model.pcls[pcl_id];
			pcl.x = 0.125 + double(x_id) * 0.25;
			pcl.y = 0.125 + double(y_id) * 0.25;
			++pcl_id;
		}
	// velocity bc
	model.vx_num = 1;
	model.vxs = new VelocityBC[model.vx_num];
	model.vxs[0].node_id = 2;
	model.vxs[0].v = 0.0;
	model.vy_num = bgm_elem_x_num + 1;
	model.vys = new VelocityBC[model.vy_num];
	for (size_t n_id = 0; n_id < model.vy_num; ++n_id)
	{
		model.vys[n_id].node_id = n_id;
		model.vys[n_id].v = 0.0;
	}

	// rigid body
	RigidBody &rb = model.rigid_body;
	rb.load_and_init_mesh("..\\..\\Asset\\cap_mesh.mesh_data", 0.25);
	rb.set_params(1.0, 0.5, 5.25, 0.0);
	//rb.add_ext_force(0.0, -1.0);
	rb.set_vtheta_bc();
	rb.set_vy_bc(-0.1);

	//DisplayModel disp_model;
	//disp_model.init_win();
	//disp_model.init_rigid_body(model.rigid_body);
	//disp_model.init_bg_mesh(bgm_h, bgm_elem_x_num, bgm_elem_y_num, -0.5, 0.0);
	//disp_model.init_particles(model.pcls, model.get_pcl_num());
	//disp_model.display(-2.0, 3.0, -0.5, 6.0);

	ResultFile_PlainBin res_file_pb;
	res_file_pb.init("mpm_rb_res_bar_compression.bin");
	ResultFile_XML res_file_xml;
	res_file_xml.init("mpm_rb_res_bar_compression.xml");

	// output model
	ModelDataOutput_S2D_ME_s_RigidBody_Fric md;
	md.set_model(model);
	md.set_res_file(res_file_pb);
	md.output();
	md.set_res_file(res_file_xml);
	md.output();

	TimeHistoryOutput_S2D_ME_s_RigidBody_Fric out1;
	out1.set_res_file(res_file_pb);
	out1.set_interval_num(50);
	out1.set_output_init_state();
	TimeHistoryOutput_S2D_ME_s_RigidBody_Fric out2;
	out2.set_res_file(res_file_xml);
	out2.set_interval_num(50);
	out2.set_output_init_state();
	TimeHistoryOutput_ConsoleProgressBar cpb;

	Step_S2D_ME_s_RigidBody_Fric step;
	step.set_name("init_step");
	step.set_model(model);
	step.set_time(5.0);
	step.set_dtime(0.01);
	step.add_time_history(out1);
	step.add_time_history(out2);
	step.add_time_history(cpb);
	step.solve();
}