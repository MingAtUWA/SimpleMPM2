#include "TestsWithGL_pcp.h"

#include "ItemArray.hpp"

#include "test_sim_core.h"

#include "TriangleMesh.h"
#include "TriangleMeshToParticles.h"
#include "Model_T2D_CHM_s.h"

#include "Step_T2D_CHM_s_SE.h"

#include "DisplayModel_T2D.h"
#include "ModelDataOutput_T2D_CHM_s.h"
#include "TimeHistoryOutput_T2D_CHM_s_SE.h"

#include "TimeHistoryOutput_ConsoleProgressBar.h"

#include "test_post_processor.h"
#include "GA_T2D_CHM_s_hdf5.h"

namespace
{

void get_left_and_right_n_ids(Model_T2D_CHM_s &md,
	MemoryUtilities::ItemArray<size_t> &n_ids)
{
	n_ids.reset();
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Model_T2D_CHM_s::Node &n = md.nodes[n_id];
		if (n.x < 1.0e-3 || n.x > 20.0 - 1.0e-3)
			n_ids.add(&n_id);
	}
}

void get_bottom_n_ids(Model_T2D_CHM_s &md,
	MemoryUtilities::ItemArray<size_t> &n_ids)
{
	n_ids.reset();
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Model_T2D_CHM_s::Node &n = md.nodes[n_id];
		if (n.y < 1.0e-3)
			n_ids.add(&n_id);
	}
}

};


void test_t2d_mpm_chm_s_t_bar_above_ground(void)
{
	TriangleMesh tri_mesh;
	tri_mesh.load_mesh("..\\..\\Asset\\rect_half.mesh_data");
	//std::cout << "node num: " << tri_mesh.get_node_num() << "\n"
	//		  << "elem num: " << tri_mesh.get_elem_num() << "\n";

	TriangleMeshToParticles mh_2_pcl(tri_mesh);
	mh_2_pcl.replace_with_grid_points(0.0, 20.0, 0.0, 15.0, 0.2, 0.2);
	//std::cout << "pcl num: " << mh_2_pcl.get_pcl_num() << "\n";

	Model_T2D_CHM_s model;
	model.init_mesh(tri_mesh);
	tri_mesh.clear();

	model.init_bg_mesh(0.5, 0.5);

	model.init_rigid_circle(2.5, 10.0, 17.5, 0.25);
	model.set_rigid_circle_velocity(0.0, -0.25, 0.0);
	//model.set_contact_stiffness(150.0, 150.0);
	model.set_contact_stiffness(1.0e5, 1.0e5);

	// elasticity
	//model.init_pcls(mh_2_pcl, 0.4, 20.0, 10.0, 100.0, 0.0, 1000.0, 1.0e-5, 1.0);
	//model.init_pcls(mh_2_pcl, 0.5, 2700.0, 1000.0, 2.0e6, 0.3, 2.0e7, 1.0e-12, 1.0e-3);
	// mcc
	model.init_pcls(mh_2_pcl, 0.6, 2650.0, 1000.0, 5.0e6, 1.0e-10, 1.0e-3);
	ModelContainer &mc = model.model_container;
	ModifiedCamClay *cms = mc.add_ModifiedCamClay(model.pcl_num);
	double ini_stress[6] = { -40000.0, -24050.0, -24050.0, 0.0, 0.0, 0.0 };
	for (size_t p_id = 0; p_id < model.pcl_num; ++p_id)
	{
		Model_T2D_CHM_s::Particle &pcl = model.pcls[p_id];
		ModifiedCamClay &mcc = cms[p_id];
		mcc.set_param_OC(0.3, 0.044, 0.205, 23.5, 3.6677, ini_stress, 39610.0);
		pcl.set_cm(mcc);
	}
	mh_2_pcl.clear();

	MemoryUtilities::ItemArray<size_t> bc_n_ids_mem;
	bc_n_ids_mem.reserve(100);
	size_t *bc_n_ids;

	get_left_and_right_n_ids(model, bc_n_ids_mem);
	bc_n_ids = bc_n_ids_mem.get_mem();
	model.init_vsxs(bc_n_ids_mem.get_num());
	for (size_t v_id = 0; v_id < model.vsx_num; ++v_id)
	{
		VelocityBC &vbc = model.vsxs[v_id];
		vbc.node_id = bc_n_ids[v_id];
		vbc.v = 0.0;
	}

	//MemoryUtilities::ItemArray<GLfloat> pt_array;
	//pt_array.reserve(25 * 3);
	//GLfloat pt_coord;
	//for (size_t n_id = 0; n_id < model.vsx_num; ++n_id)
	//{
	//	Model_T2D_CHM_s::Node &n = model.nodes[bc_n_ids[n_id]];
	//	pt_coord = double(n.x);
	//	pt_array.add(&pt_coord);
	//	pt_coord = double(n.y);
	//	pt_array.add(&pt_coord);
	//	pt_coord = 0.0f;
	//	pt_array.add(&pt_coord);
	//}

	get_bottom_n_ids(model, bc_n_ids_mem);
	bc_n_ids = bc_n_ids_mem.get_mem();
	model.init_vsys(bc_n_ids_mem.get_num());
	for (size_t v_id = 0; v_id < model.vsy_num; ++v_id)
	{
		VelocityBC &vbc = model.vsys[v_id];
		vbc.node_id = bc_n_ids[v_id];
		vbc.v = 0.0;
	}

	//for (size_t n_id = 0; n_id < model.vsy_num; ++n_id)
	//{
	//	Model_T2D_CHM_s::Node &n = model.nodes[bc_n_ids[n_id]];
	//	pt_coord = double(n.x);
	//	pt_array.add(&pt_coord);
	//	pt_coord = double(n.y);
	//	pt_array.add(&pt_coord);
	//	pt_coord = 0.0f;
	//	pt_array.add(&pt_coord);
	//}

	//// disp only one point
	//MemoryUtilities::ItemArray<GLfloat> pt_array;
	//GLfloat pt_coord;
	//Model_T2D_CHM_s::Particle &pt = model.pcls[5];
	////Model_T2D_CHM_s::Node &pt = model.nodes[14];
	//pt_array.reserve(3);
	//pt_coord = double(pt.x);
	//pt_array.add(&pt_coord);
	//pt_coord = double(pt.y);
	//pt_array.add(&pt_coord);
	//pt_coord = 0.0f;
	//pt_array.add(&pt_coord);

	//DisplayModel_T2D disp_model;
	//disp_model.init_win();
	//disp_model.init_model(model);
	//disp_model.init_rigid_circle(model.get_rigid_circle());
	////disp_model.init_points(pt_array.get_mem(), pt_array.get_num() / 3);
	//disp_model.display(-0.5, 20.5, -0.5, 20.5);
	//return;

	//ResultFile_PlainBin res_file_pb;
	//res_file_pb.init("t2d_mpm_chm_t_bar_above_ground.bin");
	//ResultFile_XML res_file_xml;
	//res_file_xml.init("t2d_mpm_chm_t_bar_above_ground.xml");
	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_mpm_chm_t_bar_above_ground.hdf5");

	// output model
	ModelDataOutput_T2D_CHM_s md("md1");
	md.set_model(model);
	//md.set_res_file(res_file_pb);
	//md.output();
	//md.set_res_file(res_file_xml);
	//md.output();
	md.set_res_file(res_file_hdf5);
	md.output();

	//TimeHistoryOutput_T2D_CHM_s_SE out1("th1");
	//out1.set_res_file(res_file_pb);
	//out1.set_output_init_state();
	//TimeHistoryOutput_T2D_CHM_s_SE out2("th2");
	//out2.set_res_file(res_file_xml);
	//out2.set_output_init_state();
	TimeHistoryOutput_ConsoleProgressBar out3;
	TimeHistoryOutput_T2D_CHM_s_SE out4("penetration");
	out4.set_res_file(res_file_hdf5);
	out4.set_output_init_state();
	out4.set_interval_num(100);
	
	Step_T2D_CHM_s_SE step;
	step.set_model(model);
	//step.set_mass_scale(10.0, 10.0);
	//step.set_damping_ratio(0.05); // local damping
	//step.set_bv_ratio(0.0); // bulk viscosity
	step.set_time(6.0);
	step.set_dtime(3.0e-6);
	//out1.set_interval_num(100);
	//step.add_time_history(out1);
	//out2.set_interval_num(100);
	//step.add_time_history(out2);
	step.add_time_history(out3);
	step.add_time_history(out4);
	step.solve();

	//system("pause");
}

void test_color_animation_t2d_chm_s_t_bar_above_ground(void)
{
	double soil_height = 20.0;
	double soil_width = 20.0;
	double padding_height = soil_height * 0.05;
	double padding_width = soil_width * 0.05;
	// Abaqus "rainbow" spectrum scheme
	ColorGraph::Colori colors[] = {
		{ 0,   0,   255 },
		{ 0,   93,  255 },
		{ 0,   185, 255 },
		{ 0,   255, 232 },
		{ 0,   255, 139 },
		{ 0,   255, 46 },
		{ 46,  255, 0 },
		{ 139, 255, 0 },
		{ 232, 255, 0 },
		{ 255, 185, 0 },
		{ 255, 93,  0 },
		{ 255, 0,   0 }
	};
	GA_T2D_CHM_s_hdf5 gen(1000, 1000);
	gen.init_color_graph(-2500.0, 100.0, colors, sizeof(colors) / sizeof(ColorGraph::Colori));
	gen.generate(
		5.0,
		-padding_width,
		soil_width + padding_width,
		-padding_height,
		soil_height + padding_height,
		"t2d_mpm_chm_t_bar_above_ground.hdf5",
		"penetration",
		"t2d_mpm_chm_t_bar_above_ground.gif");
}