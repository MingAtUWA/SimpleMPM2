#include "TestsWithGL_pcp.h"

#include "ItemArray.hpp"

#include "test_sim_core.h"

#include "TriangleMesh.h"
#include "TriangleMeshToParticles.h"
#include "Model_T2D_CHM_s.h"
#include "Step_T2D_CHM_s_SE.h"

#include "ResultFile_hdf5.h"
#include "Model_T2D_CHM_s_hdf5_io_utilities.h"

#include "ModelDataOutput_T2D_CHM_s.h"

#include "TimeHistoryOutput_T2D_CHM_s_SE.h"
#include "TimeHistoryOutput_ConsoleProgressBar.h"

#include "DisplayModel_T2D.h"

#include "test_post_processor.h"
#include "GA_T2D_CHM_s_hdf5.h"


void test_t2d_chm_s_restart_from_geostatic_hdf5_mcc(void)
{
	Model_T2D_CHM_s model;
	using Model_T2D_CHM_s_hdf5_io_utilities::load_chm_s_model_from_hdf5_file;
	load_chm_s_model_from_hdf5_file(model,
		"t2d_chm_s_geostatic_hdf5_mcc.hdf5", "th3", 100);

	std::cout << "pcl num: "   << model.pcl_num
			  << " elem num: " << model.elem_num
			  << " node num: " << model.node_num << "\n";

	model.init_bg_mesh(0.05, 0.05);

	size_t vx_bc_n_id[] = { 0, 3, 15, 16, 17, 18, 19, 20, 21, 22, 23,
							1, 2, 5,  6,  7,  8,  9,  10, 11, 12, 13 };
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

	size_t tbc_pcl_id[] = { 132, 133, 168, 169 };
	model.init_tys(sizeof(tbc_pcl_id) / sizeof(tbc_pcl_id[0]));
	for (size_t t_id = 0; t_id < model.ty_num; ++t_id)
	{
		TractionBC_MPM &tbc = model.tys[t_id];
		tbc.pcl_id = tbc_pcl_id[t_id];
		tbc.t = 0.05 * -1.0;
	}

	model.init_bfys(model.pcl_num);
	for (size_t p_id = 0; p_id < model.pcl_num; ++p_id)
	{
		BodyForce &bf = model.bfys[p_id];
		bf.pcl_id = p_id;
		bf.bf = -20.0;
	}

	//MemoryUtilities::ItemArray<GLfloat> pt_array;
	//pt_array.reserve(29 * 3);
	//GLfloat pt_coord;
	//for (size_t n_id = 0; n_id < sizeof(tbc_pcl_id)/sizeof(tbc_pcl_id[0]); ++n_id)
	//{
	//	Model_T2D_CHM_s::Particle &n = model.pcls[tbc_pcl_id[n_id]];
	//	pt_coord = double(n.x);
	//	pt_array.add(&pt_coord);
	//	pt_coord = double(n.y);
	//	pt_array.add(&pt_coord);
	//	pt_coord = 0.0f;
	//	pt_array.add(&pt_coord);
	//}

	//DisplayModel_T2D disp_model;
	//disp_model.init_win();
	//disp_model.init_model(model);
	//disp_model.init_rigid_circle(model.get_rigid_circle());
	////disp_model.init_points(pt_array.get_mem(), pt_array.get_num()/3);
	//disp_model.display(-0.05, 0.25, -0.05, 1.5);
	//return;

	ResultFile_hdf5 res_file_hdf5;
	//res_file_hdf5.open("t2d_chm_s_geostatic_hdf5_mcc.hdf5", false);
	res_file_hdf5.create("t2d_chm_s_geostatic_hdf5_mcc2.hdf5");

	ModelDataOutput_T2D_CHM_s md("md1");
	md.set_model(model);
	md.set_res_file(res_file_hdf5);
	md.output();

	TimeHistoryOutput_ConsoleProgressBar out3;
	TimeHistoryOutput_T2D_CHM_s_SE out4("th4");
	out4.set_res_file(res_file_hdf5);
	out4.set_output_init_state();
	out4.set_interval_num(100);

	Step_T2D_CHM_s_SE step;
	step.set_model(model);
	//step.set_damping_ratio(0.0); // local damping
	//step.set_bv_ratio(1.0); // bulk viscosity
	step.set_time(10.0);
	step.set_dtime(1.0e-5);
	step.add_time_history(out3);
	step.add_time_history(out4);
	step.solve();

	//system("pause");
}


void test_color_animation_t2d_chm_s_restart_from_geostatic_hdf5_mcc(void)
{
	double soil_height = 1.25;
	double soil_width = 0.2;
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
	GA_T2D_CHM_s_hdf5 gen;
	gen.init_color_graph(
		500.0,
		60.0,
		40.0,
		480.0,
		-100.0,
		0.0,
		colors,
		sizeof(colors) / sizeof(ColorGraph::Colori)
	);
	gen.generate(5.0,
		-padding_width,
		soil_width + padding_width,
		-padding_height,
		soil_height + padding_height,
		"t2d_chm_s_geostatic_hdf5_mcc2.hdf5",
		"th4",
		"t2d_chm_s_geostatic_hdf5_mcc2.gif"
	);
	//system("pause");
}
