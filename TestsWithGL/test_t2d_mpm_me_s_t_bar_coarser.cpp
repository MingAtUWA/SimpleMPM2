#include "TestsWithGL_pcp.h"

#include "ItemArray.hpp"

#include "test_sim_core.h"
#include "TriangleMesh.h"
#include "TriangleMeshToParticles.h"
#include "Model_T2D_ME_s.h"
#include "Step_T2D_ME_s.h"
#include "DisplayModel_T2D.h"
#include "ModelDataOutput_T2D_ME_s.h"
#include "TimeHistoryOutput_T2D_ME_s.h"
#include "TimeHistoryOutput_ConsoleProgressBar.h"

#include "test_post_processor.h"
#include "GA_T2D_ME_s_color.h"

namespace
{

void find_bc_pcl_and_node(Model_T2D_ME_s &md)
{
	// left
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Model_T2D_ME_s::Node &n = md.nodes[n_id];
		if (n.x < 1.0e-3)
			std::cout << n.id << ", ";
	}
	std::cout << "\n";
	// right
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Model_T2D_ME_s::Node &n = md.nodes[n_id];
		if (n.x > 10.0 - 1.0e-3)
			std::cout << n.id << ", ";
	}
	std::cout << "\n";
	// bottom
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Model_T2D_ME_s::Node &n = md.nodes[n_id];
		if (n.y < 1.0e-3)
			std::cout << n.id << ", ";
	}
	std::cout << "\n";
	// top
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Model_T2D_ME_s::Node &n = md.nodes[n_id];
		if (n.y > 15.0 - 1.0e-3)
			std::cout << n.id << ", ";
	}
	std::cout << "\n";	//// corner
	//for (size_t p_id = 0; p_id < md.pcl_num; ++p_id)
	//{
	//	Model_T2D_ME_s::Particle &pcl = md.pcls[p_id];
	//	if (pcl.y > 0.9)// && pcl.x > 0.175)
	//		std::cout << pcl.id << ", ";
	//}
	//std::cout << "\n";
	//// element
	//for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	//{
	//	Model_T2D_ME_s::Element &e = md.elems[e_id];
	//	if (e.n1 == 14 || e.n2 == 14 || e.n3 == 14)
	//		std::cout << e.id << ", ";
	//}
	//std::cout << "\n";
}

void get_left_and_right_n_ids(Model_T2D_ME_s &md,
	MemoryUtilities::ItemArray<size_t> &n_ids)
{
	n_ids.reset();
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Model_T2D_ME_s::Node &n = md.nodes[n_id];
		if (n.x < 1.0e-3 || n.x > 20.0 - 1.0e-3)
			n_ids.add(&n_id);
	}
}

void get_top_and_bottom_n_ids(Model_T2D_ME_s &md,
	MemoryUtilities::ItemArray<size_t> &n_ids)
{
	n_ids.reset();
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Model_T2D_ME_s::Node &n = md.nodes[n_id];
		if (n.y < 1.0e-3 || n.y > 30.0 - 1.0e-3)
			n_ids.add(&n_id);
	}
}

};


void test_t2d_mpm_me_s_t_bar_coarser(void)
{
	TriangleMesh tri_mesh;
	tri_mesh.load_mesh("..\\..\\Asset\\rect_t_bar_coarser.mesh_data");
	//std::cout << "node num: " << tri_mesh.get_node_num() << "\n"
	//		  << "elem num: " << tri_mesh.get_elem_num() << "\n";
	
	TriangleMeshToParticles mh_2_pcl(tri_mesh);
	//mh_2_pcl.set_generator(TriangleMeshToParticles::GeneratorType::FirstOrderGaussPoint);
	//mh_2_pcl.generate_pcls();
	//mh_2_pcl.replace_with_grid_points(10.0, 20.0, 10.0, 30.0, 0.25, 0.25);
	mh_2_pcl.replace_with_grid_points(0.0, 20.0, 0.0, 30.0, 0.2, 0.2);
	//std::cout << "pcl num: " << mh_2_pcl.get_pcl_num() << "\n";

	Model_T2D_ME_s model;
	model.init_mesh(tri_mesh);
	tri_mesh.clear();
	//double area = 0.0;
	//for (size_t i = 0; i < model.elem_num; i++)
	//{
	//	Model_T2D_ME_s::Element &e = model.elems[i];
	//	area += e.area_2;
	//}
	//std::cout << area << "\n";

	model.init_rigid_circle(2.5, 10.0, 20.0, 0.25);
	model.set_rigid_circle_velocity(0.0, -1.0, 0.0);
	model.get_rigid_circle().del_pcls_in_circle(mh_2_pcl, 0.1);

	model.init_pcls(mh_2_pcl, 20.0, 100.0, 0.0);
	mh_2_pcl.clear();
	//double pcl_area = 0.0;
	//for (size_t i = 0; i < model.pcl_num; i++)
	//{
	//	Model_T2D_ME_s::Particle &pcl = model.pcls[i];
	//	pcl_area += pcl.m_s;
	//}
	//pcl_area /= (1.0 - 0.2) * 2.0;
	//std::cout << pcl_area << "\n";
	
	model.init_bg_mesh(1.0, 1.0);

	//find_bc_pcl_and_node(model);
	//system("pause");
	//return;

	MemoryUtilities::ItemArray<size_t> bc_n_ids_mem;
	bc_n_ids_mem.reserve(100);
	size_t *bc_n_ids;

	get_left_and_right_n_ids(model, bc_n_ids_mem);
	bc_n_ids = bc_n_ids_mem.get_mem();
	model.init_vxs(bc_n_ids_mem.get_num());
	for (size_t v_id = 0; v_id < model.vx_num; ++v_id)
	{
		VelocityBC &vbc = model.vxs[v_id];
		vbc.node_id = bc_n_ids[v_id];
		vbc.v = 0.0;
	}

	//MemoryUtilities::ItemArray<GLfloat> pt_array;
	//pt_array.reserve(25 * 3);
	//GLfloat pt_coord;
	//for (size_t n_id = 0; n_id < model.vx_num; ++n_id)
	//{
	//	Model_T2D_ME_s::Node &n = model.nodes[bc_n_ids[n_id]];
	//	pt_coord = double(n.x);
	//	pt_array.add(&pt_coord);
	//	pt_coord = double(n.y);
	//	pt_array.add(&pt_coord);
	//	pt_coord = 0.0f;
	//	pt_array.add(&pt_coord);
	//}

	get_top_and_bottom_n_ids(model, bc_n_ids_mem);
	bc_n_ids = bc_n_ids_mem.get_mem();
	model.init_vys(bc_n_ids_mem.get_num());
	for (size_t v_id = 0; v_id < model.vy_num; ++v_id)
	{
		VelocityBC &vbc = model.vys[v_id];
		vbc.node_id = bc_n_ids[v_id];
		vbc.v = 0.0;
	}

	//for (size_t n_id = 0; n_id < model.vy_num; ++n_id)
	//{
	//	Model_T2D_ME_s::Node &n = model.nodes[bc_n_ids[n_id]];
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
	//Model_T2D_ME_s::Particle &pt = model.pcls[5];
	////Model_T2D_ME_s::Node &pt = model.nodes[14];
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
	//disp_model.init_points(pt_array.get_mem(), pt_array.get_num() / 3);
	//disp_model.display(-0.5, 20.5, -0.5, 30.5);
	//return;

	ResultFile_PlainBin res_file_pb;
	res_file_pb.init("t2d_mpm_me_t_bar_coarser.bin");
	ResultFile_XML res_file_xml;
	res_file_xml.init("t2d_mpm_me_t_bar_coarser.xml");
	
	// output model
	ModelDataOutput_T2D_ME_s md;
	md.set_model(model);
	md.set_res_file(res_file_pb);
	md.output();
	md.set_res_file(res_file_xml);
	md.output();

	TimeHistoryOutput_T2D_ME_s out1("tbar_penetration");
	out1.set_res_file(res_file_pb);
	out1.set_output_init_state();
	TimeHistoryOutput_T2D_ME_s out2("tbar_penetration");
	out2.set_res_file(res_file_xml);
	out2.set_output_init_state();
	TimeHistoryOutput_ConsoleProgressBar out3;

	Step_T2D_ME_s step;
	step.set_model(model);
	//step.set_damping_ratio(0.05); // local damping
	//step.set_bv_ratio(0.0); // bulk viscosity
	step.set_time(4.0);
	step.set_dtime(5.0e-6);
	out1.set_interval_num(100);
	step.add_time_history(out1);
	out2.set_interval_num(100);
	step.add_time_history(out2);
	step.add_time_history(out3);
	step.solve();

	//Step_T2D_ME_s step2;
	//step2.set_prev_step(step);
	//step2.set_time(0.01);
	//step2.set_dtime(1.0e-4);
	//out1.set_interval_num(100);
	//step2.add_time_history(out1);
	//out2.set_interval_num(100);
	//step2.add_time_history(out2);
	////step2.add_time_history(out3);
	//step2.solve();

	//system("pause");
}

void test_color_animation_t2d_me_s_t_bar_coarser(void)
{
	double soil_height = 30.0;
	double soil_width  = 20.0;
	double padding_height = soil_height * 0.05;
	double padding_width  = soil_width * 0.05;
	// Abaqus "rainbow" spectrum scheme
	ColorGraph::Colori colors[] = {
		{ 0,   0,   255 },
		{ 0,   93,  255 },
		{ 0,   185, 255 },
		{ 0,   255, 232 },
		{ 0,   255, 139 },
		{ 0,   255, 46  },
		{ 46,  255, 0   },
		{ 139, 255, 0   },
		{ 232, 255, 0   },
		{ 255, 185, 0   },
		{ 255, 93,  0   },
		{ 255, 0,   0   }
	};
	GA_T2D_ME_s_color gen(1000, 1000);
	gen.init_color_graph(-200.0, 200.0, colors, sizeof(colors) / sizeof(ColorGraph::Colori));
	gen.generate(5.0, -padding_width, soil_width + padding_width,
				 -padding_height, soil_height + padding_height,
				 "t2d_mpm_me_t_bar_coarser.bin", "t2d_mpm_me_t_bar_coarser.gif");
}
