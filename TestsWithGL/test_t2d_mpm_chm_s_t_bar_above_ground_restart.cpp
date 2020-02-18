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

	void get_top_pcl_ids(Model_T2D_CHM_s &md,
		MemoryUtilities::ItemArray<size_t> &pcl_ids)
	{
		pcl_ids.reset();
		for (size_t p_id = 0; p_id < md.pcl_num; ++p_id)
		{
			Model_T2D_CHM_s::Particle &pcl = md.pcls[p_id];
			if (pcl.y > 14.85 && pcl.y < 15.0)
				pcl_ids.add(&p_id);
		}
	}
};


void test_t2d_mpm_chm_s_t_bar_above_ground_restart(void)
{
	Model_T2D_CHM_s model;

	using Model_T2D_CHM_s_hdf5_io_utilities::load_chm_s_model_from_hdf5_file;
	load_chm_s_model_from_hdf5_file(
		model,
		"t2d_mpm_chm_t_bar_real_geostatic.hdf5",
		"penetration",
		100
	);

	model.init_bg_mesh(0.5, 0.5);

	// pcl at the top
	MemoryUtilities::ItemArray<size_t> bc_pcl_ids_mem;
	bc_pcl_ids_mem.reserve(100);
	get_top_pcl_ids(model, bc_pcl_ids_mem);
	size_t *tbc_pcl_ids = bc_pcl_ids_mem.get_mem();
	model.init_tys(bc_pcl_ids_mem.get_num());
	for (size_t t_id = 0; t_id < model.ty_num; ++t_id)
	{
		TractionBC_MPM &tbc = model.tys[t_id];
		tbc.pcl_id = tbc_pcl_ids[t_id];
		tbc.t = 0.2 * -40000.0;
	}

	MemoryUtilities::ItemArray<GLfloat> pt_array;
	pt_array.reserve(25 * 3);
	GLfloat pt_coord;
	for (size_t t_id = 0; t_id < model.ty_num; ++t_id)
	{
		TractionBC_MPM &tbc = model.tys[t_id];
		Model_T2D_CHM_s::Particle &pcl = model.pcls[tbc_pcl_ids[t_id]];
		pt_coord = double(pcl.x);
		pt_array.add(&pt_coord);
		pt_coord = double(pcl.y);
		pt_array.add(&pt_coord);
		pt_coord = 0.0f;
		pt_array.add(&pt_coord);
	}

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

	DisplayModel_T2D disp_model;
	disp_model.init_win();
	disp_model.init_model(model);
	disp_model.init_rigid_circle(model.get_rigid_circle());
	disp_model.init_points(pt_array.get_mem(), pt_array.get_num() / 3);
	//disp_model.display(-0.5, 20.5, -0.5, 20.5);
	disp_model.display(5.0, 15.0, 10.0, 20.0);
	return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_mpm_chm_t_bar_above_ground.hdf5");

	// output model
	ModelDataOutput_T2D_CHM_s md("md1");
	md.set_model(model);
	md.set_res_file(res_file_hdf5);
	md.output();

	TimeHistoryOutput_T2D_CHM_s_SE out1("penetration");
	out1.set_res_file(res_file_hdf5);
	out1.set_output_init_state();
	out1.set_interval_num(100);
	TimeHistoryOutput_ConsoleProgressBar out2;

	Step_T2D_CHM_s_SE step;
	step.set_model(model);
	step.set_time(6.0);
	step.set_dtime(3.0e-6);
	step.add_time_history(out1);
	step.add_time_history(out2);
	step.solve();
}
