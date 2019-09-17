#include "TestsWithGL_pcp.h"

#include <fstream>

#include "Model_S2D_ME_s_RigidBody.h"
#include "Step_S2D_ME_s_RigidBody.h"
#include "TimeHistory_S2D_ME_s_RigidBody.h"
#include "ResultFile_PlainBin.h"
#include "ResultFile_XML.h"

#include "test_sim_core.h"

void test_mpm_rigidbody_cantilever(void)
{
	Model_S2D_ME_s_RigidBody model;

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
			Model_S2D_ME_s_RigidBody::Particle &pcl = model.pcls[pcl_id];
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
	res_file_pb.init("mpm_rb_res_can.bin");
	ResultFile_XML res_file_xml;
	res_file_xml.init("mpm.rb_res_can.xml");

	// output model
	res_file_pb.output(model);
	res_file_xml.output(model);
	
	Step_S2D_ME_s_RigidBody step;
	step.set_name("init_step");
	step.set_model(model);
	step.set_time(10.0);
	step.set_dtime(0.01);

	TimeHistory_S2D_ME_s_RigidBody out1;
	out1.set_res_file(res_file_pb);
	out1.set_interval_num(50);
	out1.set_output_init_state();
	TimeHistory_S2D_ME_s_RigidBody out2;
	out2.set_res_file(res_file_xml);
	out2.set_interval_num(50);
	out2.set_output_init_state();

	step.add_time_history(out1);
	step.add_time_history(out2);
	step.solve();

	//system("pause");
}
