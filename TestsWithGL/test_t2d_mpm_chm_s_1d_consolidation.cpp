#include "TestsWithGL_pcp.h"

#include "test_sim_core.h"

#include "ItemArray.hpp"

#include "AdjustParticlesWithTriangleMesh.hpp"
#include "Model_T2D_CHM_s.h"
#include "Step_T2D_CHM_s_SE.h"
#include "DisplayModel_T2D.h"
#include "ModelDataOutput_T2D_CHM_s.h"
#include "TimeHistoryOutput_T2D_CHM_s_SE.h"
#include "TimeHistoryOutput_ConsoleProgressBar.h"

#include "test_post_processor.h"

#include "GA_T2D_CHM_s.h"
#include "GA_T2D_CHM_s_color.h"
#include "GA_T2D_CHM_s_hdf5.h"

void test_t2d_mpm_square(void)
{
	double n_coords[] = { 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0 };
	size_t e_ids[] = { 0, 1, 2, 0, 2, 3 };

	Model_T2D_CHM_s model;
	model.init_mesh(n_coords, 4, e_ids, 2);
	model.init_pcls(2, 0.2, 800.0, 2000.0, 1000.0, 100.0, 0.2, 5000.0, 1.0e-4, 1.0);
	model.pcls[0].x = 0.25;
	model.pcls[0].y = 0.50;
	model.pcls[1].x = 0.75;
	model.pcls[1].y = 0.50;

	model.tx_num = 2;
	model.txs = new TractionBC_MPM[model.tx_num];
	for (size_t t_id = 0; t_id < model.tx_num; ++t_id)
	{
		TractionBC_MPM &tbc = model.txs[t_id];
		tbc.pcl_id = t_id;
		tbc.t = 1.0;
	}

	//model.asx_num = 2;
	//model.asxs = new AccelerationBC[model.asx_num];
	//for (size_t a_id = 0; a_id < model.asx_num; ++a_id)
	//{
	//	AccelerationBC &abc = model.asxs[a_id];
	//	abc.node_id = ;
	//	abc.a = 0.0;
	//}
	
	Step_T2D_CHM_s_SE step;
	step.set_model(model);
	step.set_time(0.1);
	step.set_dtime(0.1);
	step.solve();

	system("pause");
}

namespace
{

void find_bc_pcl_and_node(Model_T2D_CHM_s &md)
{
	//// left
	//for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	//{
	//	Model_T2D_CHM_s::Node &n = md.nodes[n_id];
	//	if (n.x < 1.0e-3)
	//		std::cout << n.id << ", ";
	//}
	//std::cout << "\n";
	//// right
	//for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	//{
	//	Model_T2D_CHM_s::Node &n = md.nodes[n_id];
	//	if (n.x > 0.2 - 1.0e-3)
	//		std::cout << n.id << ", ";
	//}
	//std::cout << "\n";
	//// bottom
	//for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	//{
	//	Model_T2D_CHM_s::Node &n = md.nodes[n_id];
	//	if (n.y < 1.0e-3)
	//		std::cout << n.id << ", ";
	//}
	//std::cout << "\n";

	// top particles
	for (size_t p_id = 0; p_id < md.pcl_num; ++p_id)
	{
		Model_T2D_CHM_s::Particle &pcl = md.pcls[p_id];
		if (pcl.y > 1.0 - 0.02)
			std::cout << pcl.id << ", ";
	}
	std::cout << "\n";

	// bottom particles
	for (size_t p_id = 0; p_id < md.pcl_num; ++p_id)
	{
		Model_T2D_CHM_s::Particle &pcl = md.pcls[p_id];
		if (pcl.y < 0.02)// && pcl.x > 0.175)
			std::cout << pcl.id << ", ";
	}
	std::cout << "\n";

	//// element
	//for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	//{
	//	Model_T2D_CHM_s::Element &e = md.elems[e_id];
	//	if (e.n1 == 14 || e.n2 == 14 || e.n3 == 14)
	//		std::cout << e.id << ", ";
	//}
	//std::cout << "\n";
}

void get_top_pcl_ids(Model_T2D_CHM_s &md,
	MemoryUtilities::ItemArray<size_t> &pcl_ids)
{
	pcl_ids.reset();
	for (size_t p_id = 0; p_id < md.pcl_num; ++p_id)
	{
		Model_T2D_CHM_s::Particle &pcl = md.pcls[p_id];
		if (pcl.y > 0.987)
			pcl_ids.add(&p_id);
	}
}

};


void test_t2d_mpm_chm_s_1d_consolidation(void)
{
	TriangleMesh tri_mesh;
	tri_mesh.load_mesh("..\\..\\Asset\\rect.mesh_data");
	
	Model_T2D_CHM_s model;
	model.init_mesh(tri_mesh);
	model.init_bg_mesh(0.05, 0.05);

	TriangleMeshToParticles mh_2_pcl(tri_mesh);
	//mh_2_pcl.set_even_div_num(2);
	//mh_2_pcl.set_generator(TriangleMeshToParticles::GeneratorType::EvenlyDistributedPoint);
	//mh_2_pcl.generate_pcls();
	//
	mh_2_pcl.generate_grid_points(0.0, 0.2, 0.0, 1.0, 0.02, 0.02);

	AdjustParticlesWithTriangleMesh<Model_T2D_CHM_s> ap_mh(mh_2_pcl);
	ap_mh.distribute_points_area_to_mesh(model, 0.02, 0.02, 0.2, 0.2);

	model.init_pcls(mh_2_pcl, 0.4, 20.0, 10.0, 1000.0, 0.0, 40000.0, 1.0e-4, 1.0);
	
	tri_mesh.clear();
	mh_2_pcl.clear();

	size_t vx_bc_n_id[] = { 0, 3, 15, 16, 17, 18, 19, 20, 21, 22, 23,
							1, 2,  5,  6,  7,  8,  9, 10, 11, 12, 13 };
	model.init_vsxs(sizeof(vx_bc_n_id) / sizeof(vx_bc_n_id[0]));
	for (size_t v_id = 0; v_id < model.vsx_num; ++v_id)
	{
		VelocityBC &vbc = model.vsxs[v_id];
		vbc.node_id = vx_bc_n_id[v_id];
		vbc.v = 0.0;
	}
	model.init_vfxs(sizeof(vx_bc_n_id) / sizeof(vx_bc_n_id[0]));
	for (size_t v_id = 0; v_id < model.vfx_num; ++v_id)
	{
		VelocityBC &vbc = model.vfxs[v_id];
		vbc.node_id = vx_bc_n_id[v_id];
		vbc.v = 0.0;
	}

	size_t vy_bc_n_id[] = { 0, 1, 4 };
	model.init_vsys(sizeof(vy_bc_n_id) / sizeof(vy_bc_n_id[0]));
	for (size_t v_id = 0; v_id < model.vsy_num; ++v_id)
	{
		VelocityBC &vbc = model.vsys[v_id];
		vbc.node_id = vy_bc_n_id[v_id];
		vbc.v = 0.0;
	}
	model.init_vfys(sizeof(vy_bc_n_id) / sizeof(vy_bc_n_id[0]));
	for (size_t v_id = 0; v_id < model.vfy_num; ++v_id)
	{
		VelocityBC &vbc = model.vfys[v_id];
		vbc.node_id = vy_bc_n_id[v_id];
		vbc.v = 0.0;
	}
	
	MemoryUtilities::ItemArray<size_t> tbc_pcls;
	get_top_pcl_ids(model, tbc_pcls);
	size_t *tbc_pcl_id = tbc_pcls.get_mem();
	model.init_tys(tbc_pcls.get_num());
	for (size_t t_id = 0; t_id < model.ty_num; ++t_id)
	{
		TractionBC_MPM &tbc = model.tys[t_id];
		tbc.pcl_id = tbc_pcl_id[t_id];
		//tbc.t = 0.05 * -1.0;
		tbc.t = 0.02 * -10.0;
	}

	MemoryUtilities::ItemArray<GLfloat> pt_array;
	pt_array.reserve(29 * 3);
	GLfloat pt_coord;

	//for (size_t n_id = 0; n_id < sizeof(vx_bc_n_id)/sizeof(vx_bc_n_id[0]); ++n_id)
	//{
	//	Model_T2D_CHM_s::Node &n = model.nodes[vx_bc_n_id[n_id]];
	//	pt_coord = double(n.x);
	//	pt_array.add(&pt_coord);
	//	pt_coord = double(n.y);
	//	pt_array.add(&pt_coord);
	//	pt_coord = 0.0f;
	//	pt_array.add(&pt_coord);
	//}
	//for (size_t n_id = 0; n_id < sizeof(vy_bc_n_id) / sizeof(vy_bc_n_id[0]); ++n_id)
	//{
	//	Model_T2D_CHM_s::Node &n = model.nodes[vy_bc_n_id[n_id]];
	//	pt_coord = double(n.x);
	//	pt_array.add(&pt_coord);
	//	pt_coord = double(n.y);
	//	pt_array.add(&pt_coord);
	//	pt_coord = 0.0f;
	//	pt_array.add(&pt_coord);
	//}

	//for (size_t p_id = 0; p_id < model.ty_num; ++p_id)
	//{
	//	Model_T2D_CHM_s::Particle &pcl = model.pcls[tbc_pcl_id[p_id]];
	//	pt_coord = double(pcl.x);
	//	pt_array.add(&pt_coord);
	//	pt_coord = double(pcl.y);
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
	//disp_model.init_points(pt_array.get_mem(), pt_array.get_num()/3);
	//disp_model.display(-0.05, 0.25, -0.05, 1.05);
	//return;

	ResultFile_hdf5 res_file_h5;
	res_file_h5.create("t2d_mpm_1d_consolidation.h5");
	
	// output model
	ModelDataOutput_T2D_CHM_s md("md1");
	md.set_model(model);
	md.set_res_file(res_file_h5);
	md.output();

	TimeHistoryOutput_T2D_CHM_s_SE out1("consolidation");
	out1.set_interval_num(100);
	out1.set_res_file(res_file_h5);
	out1.set_output_init_state();

	TimeHistoryOutput_ConsoleProgressBar out3;

	Step_T2D_CHM_s_SE step;
	step.set_model(model);
	step.set_time(15.0);
	step.set_dtime(1.0e-5);
	step.add_time_history(out1);
	step.add_time_history(out3);
	step.solve();
}

void test_color_animation_t2d_chm_s_1d_consolidation(void)
{
	double soil_height = 1.0;
	double soil_width = 0.2;
	double padding_height = soil_height * 0.05;
	double padding_width = soil_width * 0.05;
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
	GA_T2D_CHM_s_hdf5 gen(1000, 1000); // window size
	gen.init_color_graph(
		850.0, 500.0, 50.0, 450.0,
		0.0, 10.0,
		colors, sizeof(colors) / sizeof(ColorGraph::Colori)
	);
	gen.generate(5.0,
		-padding_width, soil_width + padding_width,
		-padding_height, soil_height + padding_height,
		"t2d_mpm_1d_consolidation.h5",
		"consolidation",
		"t2d_mpm_1d_consolidation.gif"
	);

}

void test_animation_t2d_chm_s_1d_consolidation(void)
{
	double soil_height = 1.0;
	double soil_width = 0.2;
	double padding_height = soil_height * 0.05;
	double padding_width = soil_width * 0.05;
	GA_T2D_CHM_s gen;
	gen.generate(5.0, -padding_width, soil_width + padding_width,
		-padding_height, soil_height + padding_height,
		"t2d_mpm_1d_consolidation.bin",
		"t2d_mpm_1d_consolidation.gif");
}
