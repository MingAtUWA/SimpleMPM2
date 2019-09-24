#include "TestsWithGL_pcp.h"

#include "test_sim_core.h"

#include "Model_S2D_ME_s_RigidBody_Fric.h"
#include "Step_S2D_ME_s_RigidBody_Fric.h"
#include "TimeHistory_S2D_ME_s_RigidBody_Fric.h"
#include "ResultFile_PlainBin.h"
#include "ResultFile_XML.h"
#include "TimeHistory_ConsoleProgressBar.h"

#include "DisplayModel.h"

#define TO_RADIAN(an) ((an) * 3.14159265359 / 180.0)

void test_slide_down_frictional_slope2(void)
{
	Model_S2D_ME_s_RigidBody_Fric model;
	
	// background mesh
	double bgm_h = 1.0;
	size_t bgm_elem_x_num = 9, bgm_elem_y_num = 6;
	model.init_mesh(bgm_h, bgm_elem_x_num, bgm_elem_y_num);

	// material points
	double slope_height = 5.0;
	double angle = TO_RADIAN(30.0);
	double pcl_vol = 0.01;
	double pcl_len = sqrt(pcl_vol);
	double dx = pcl_len * cos(angle), dy = pcl_len * sin(angle);
	size_t pcl_num_1line = size_t(ceil(slope_height / dy));
	model.init_pcl(pcl_num_1line * 2, pcl_vol, 1.0, 100.0, 0.0);
	for (size_t pcl_id = 0; pcl_id < pcl_num_1line; ++pcl_id)
	{
		Model_S2D_ME_s_RigidBody_Fric::Particle &pcl1 = model.pcls[pcl_id];
		pcl1.x = double(pcl_id) * dx;
		pcl1.y = slope_height - double(pcl_id) * dy;
		Model_S2D_ME_s_RigidBody_Fric::Particle &pcl2 = model.pcls[pcl_id + pcl_num_1line];
		pcl2.x = -pcl_len * sin(angle) + double(pcl_id) * dx;
		pcl2.y = slope_height - pcl_len * cos(angle) - double(pcl_id) * dy;
	}
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

	// rigid body
	RigidBody &rb = model.rigid_body;
	rb.load_and_init_mesh("..\\..\\Asset\\square_mesh.mesh_data", 0.3);
	rb.set_params(1.0, 0.6830127 + 0.5*dy, 5.1830127 + 0.5*dx, -0.52359877); // need to be changed
	rb.add_ext_force(0.0, -2.0);

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
	res_file_pb.output(model);
	res_file_xml.output(model);

	TimeHistory_S2D_ME_s_RigidBody_Fric out1;
	out1.set_res_file(res_file_pb);
	out1.set_interval_num(50);
	out1.set_output_init_state();
	TimeHistory_S2D_ME_s_RigidBody_Fric out2;
	out2.set_res_file(res_file_xml);
	out2.set_interval_num(50);
	out2.set_output_init_state();
	TimeHistory_ConsoleProgressBar cpb;

	Step_S2D_ME_s_RigidBody_Fric step;
	step.set_name("init_step");
	step.set_model(model);
	step.set_time(5.0);
	step.set_dtime(0.02);
	step.add_time_history(out1);
	step.add_time_history(out2);
	step.add_time_history(cpb);
	step.solve();
}
