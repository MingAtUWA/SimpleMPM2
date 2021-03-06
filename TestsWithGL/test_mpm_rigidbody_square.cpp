#include "TestsWithGL_pcp.h"

#include <fstream>

#include "Model_S2D_ME_s_RigidBody.h"
#include "Step_S2D_ME_s_RigidBody.h"
#include "ModelDataOutput_S2D_ME_s_RigidBody.h"
#include "TimeHistoryOutput_S2D_ME_s_RigidBody.h"
#include "TimeHistoryOutput_ConsoleProgressBar.h"

#include "test_sim_core.h"

void test_mpm_rigidbody_square(void)
{
	Model_S2D_ME_s_RigidBody model;

	model.rigid_body.load_and_init_mesh("..\\..\\Asset\\square_mesh.mesh_data", 0.3);
	model.rigid_body.set_params(1.0, 0.0, 2.5, 0.0);
	model.rigid_body.add_ext_force(0.0, -2.0);

	model.init_mesh(1.0, 3, 2, -1.5, 0.0);
	model.init_pcl(4, 0.25, 1.0, 100.0, 0.0);
	size_t pcl_id = 0;
	for (size_t y_id = 0; y_id < 2; y_id++)
		for (size_t x_id = 0; x_id < 2; x_id++)
		{
			Model_S2D_ME_s_RigidBody::Particle &pcl = model.pcls[pcl_id];
			pcl.x = -0.25 + double(x_id) * 0.5;
			pcl.y =  0.25 + double(y_id) * 0.5;
			++pcl_id;
		}
	// vy bc
	model.vy_num = 3;
	model.vys = new VelocityBC[model.vy_num];
	model.vys[0].node_id = 1;
	model.vys[0].v = 0.0;
	model.vys[1].node_id = 2;
	model.vys[1].v = 0.0;
	model.vys[2].node_id = 3;
	model.vys[2].v = 0.0;
	// vx bc
	model.vx_num = 1;
	model.vxs = new VelocityBC[model.vx_num];
	model.vxs[0].node_id = 2;
	model.vxs[0].v = 0.0;

	ResultFile_PlainBin res_file_pb;
	res_file_pb.init("mpm_rb_res_square.bin");
	ResultFile_XML res_file_xml;
	res_file_xml.init("mpm_rb_res_square.xml");

	// output model
	ModelDataOutput_S2D_ME_s_RigidBody md;
	md.set_model(model);
	md.set_res_file(res_file_pb);
	md.output();
	md.set_res_file(res_file_xml);
	md.output();
	
	TimeHistoryOutput_S2D_ME_s_RigidBody out1;
	out1.set_res_file(res_file_pb);
	out1.set_interval_num(50);
	TimeHistoryOutput_S2D_ME_s_RigidBody out2;
	out2.set_res_file(res_file_xml);
	out2.set_interval_num(50);
	TimeHistoryOutput_ConsoleProgressBar cpb;

	Step_S2D_ME_s_RigidBody step;
	step.set_name("init_step");
	step.set_model(model);
	step.set_time(5.0);
	step.set_dtime(0.01);
	step.add_time_history(out1);
	step.add_time_history(out2);
	//step.add_time_history(cpb);
	step.solve();

	//system("pause");
}
