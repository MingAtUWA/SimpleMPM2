#include "TestsWithGL_pcp.h"

#include <fstream>

#include "Model_S2D_ME_s_RigidBody_Fric.h"
#include "Step_S2D_ME_s_RigidBody_Fric.h"
#include "ModelDataOutput_S2D_ME_s_RigidBody_Fric.h"
#include "TimeHistoryOutput_S2D_ME_s_RigidBody_Fric.h"
#include "TimeHistoryOutput_ConsoleProgressBar.h"

#include "ResultFile_PlainBin.h"
#include "ResultFile_XML.h"

#include "test_sim_core.h"

void test_mpm_rigidbody_cantilever_fric(void)
{
	Model_S2D_ME_s_RigidBody_Fric model;
	model.rigid_body.load_and_init_mesh("..\\..\\Asset\\square_mesh.mesh_data", 0.3);
	model.rigid_body.set_params(1.0, 3.0, 5.0, 0.0);
	//model.rigid_body.add_ext_force(0.0, -1.0, 3.5, 5.0);
	model.rigid_body.add_ext_force(0.0, -1.0);

	model.init_mesh(1.0, 5, 5);
	model.init_pcl(16, 0.25, 1.0, 100.0, 0.0);
	size_t pcl_id = 0;
	for (size_t y_id = 0; y_id < 2; y_id++)
		for (size_t x_id = 0; x_id < 8; x_id++)
		{
			Model_S2D_ME_s_RigidBody_Fric::Particle &pcl = model.pcls[pcl_id];
			pcl.x = 0.25 + double(x_id) * 0.5;
			pcl.y = 3.25 + double(y_id) * 0.5;
			++pcl_id;
		}
	// traction
	model.ty_num = 2;
	model.tys = new TractionBC_MPM[model.ty_num];
	model.tys[0].pcl_id = 7;
	model.tys[0].t = -0.0;
	model.tys[1].pcl_id = 8;
	model.tys[1].t = -0.0;
	// vy bc
	model.vy_num = 2;
	model.vys = new VelocityBC[model.vy_num];
	model.vys[0].node_id = 18;
	model.vys[0].v = 0.0;
	model.vys[1].node_id = 24;
	model.vys[1].v = 0.0;
	// vx bc
	model.vx_num = 2;
	model.vxs = new VelocityBC[model.vx_num];
	model.vxs[0].node_id = 18;
	model.vxs[0].v = 0.0;
	model.vxs[1].node_id = 24;
	model.vxs[1].v = 0.0;

	ResultFile_PlainBin res_file_pb;
	res_file_pb.init("mpm_rb_res_can_fric.bin");
	ResultFile_XML res_file_xml;
	res_file_xml.init("mpm_rb_res_can_fric.xml");

	// output model
	ModelDataOutput_S2D_ME_s_RigidBody_Fric md;
	md.set_model(model);
	md.set_res_file(res_file_pb);
	md.output();
	md.set_res_file(res_file_xml);
	md.output();

	Step_S2D_ME_s_RigidBody_Fric step;
	step.set_name("init_step");
	step.set_model(model);
	step.set_time(10.0);
	step.set_dtime(0.01);

	TimeHistoryOutput_S2D_ME_s_RigidBody_Fric out1;
	out1.set_res_file(res_file_pb);
	out1.set_interval_num(50);
	out1.set_output_init_state();
	TimeHistoryOutput_S2D_ME_s_RigidBody_Fric out2;
	out2.set_res_file(res_file_xml);
	out2.set_interval_num(50);
	out2.set_output_init_state();
	TimeHistoryOutput_ConsoleProgressBar cpb;

	step.add_time_history(out1);
	step.add_time_history(out2);
	step.add_time_history(cpb);
	step.solve();

	//system("pause");
}
