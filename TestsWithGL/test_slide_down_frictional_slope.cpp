#include "TestsWithGL_pcp.h"

#include "test_sim_core.h"

#include "Model_S2D_ME_s_RigidBody_Fric.h"
#include "Step_S2D_ME_s_RigidBody_Fric.h"
#include "ModelDataOutput_S2D_ME_s_RigidBody_Fric.h"
#include "TimeHistoryOutput_S2D_ME_s_RigidBody_Fric.h"
#include "TimeHistoryOutput_ConsoleProgressBar.h"

#include "ResultFile_PlainBin.h"
#include "ResultFile_XML.h"

#include "DisplayModel.h"

void test_slide_down_frictional_slope(void)
{
	Model_S2D_ME_s_RigidBody_Fric model;
	
	// rigid body
	RigidBody &rb = model.rigid_body;
	rb.load_and_init_mesh("..\\..\\Asset\\square_mesh.mesh_data", 0.3);
	rb.set_params(1.0, 0.6830127, 5.1830127, -0.52359877);
	rb.add_ext_force(0.0, -2.0);

	// background mesh
	double bgm_h = 1.0;
	size_t bgm_elem_x_num = 9, bgm_elem_y_num = 6;
	model.init_mesh(bgm_h, bgm_elem_x_num, bgm_elem_y_num);

	// material points
	model.get_pcls_from_mesh("..\\..\\Asset\\triangle_mesh.mesh_data", 1.0, 100.0, 0.0, 0.2);
	// velocity bc
	model.vx_num = model.get_node_num();
	model.vxs = new VelocityBC[model.vx_num];
	model.vy_num = model.vx_num;
	model.vys = new VelocityBC[model.vy_num];
	for (size_t n_id = 0; n_id < model.vx_num; ++n_id)
	{
		model.vxs[n_id].node_id = n_id;
		model.vxs[n_id].v = 0.0;
		model.vys[n_id].node_id = n_id;
		model.vys[n_id].v = 0.0;
	}

	//DisplayModel disp_model;
	//disp_model.init_win();
	//disp_model.init_rigid_body(model.rigid_body);
	//disp_model.init_bg_mesh(bgm_h, bgm_elem_x_num, bgm_elem_y_num);
	//disp_model.init_particles(model.pcls, model.get_pcl_num());
	//disp_model.display(-0.5, 9.5, -1.0, 7.0);

	ResultFile_PlainBin res_file_pb;
	res_file_pb.init("mpm_rb_res_slope_fric.bin");
	ResultFile_XML res_file_xml;
	res_file_xml.init("mpm_rb_res_slope_fric.xml");

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
	//out1.set_output_init_state();
	TimeHistoryOutput_S2D_ME_s_RigidBody_Fric out2;
	out2.set_res_file(res_file_xml);
	out2.set_interval_num(50);
	//out2.set_output_init_state();
	TimeHistoryOutput_ConsoleProgressBar cpb;

	Step_S2D_ME_s_RigidBody_Fric step;
	step.set_name("init_step");
	step.set_model(model);
	step.set_time(5.0);
	step.set_dtime(0.05);
	step.add_time_history(out1);
	step.add_time_history(out2);
	step.add_time_history(cpb);
	step.solve();
}
