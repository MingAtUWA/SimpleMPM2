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
void get_left_and_right_n_ids(
	Model_T2D_CHM_s &md, double left, double right,
	MemoryUtilities::ItemArray<size_t> &n_ids)
{
	n_ids.reset();
	double left_tol = abs(left) * 1.0e-3;
	double right_tol = abs(right) * 1.0e-3;
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Model_T2D_CHM_s::Node &n = md.nodes[n_id];
		if (n.x > left - left_tol   && n.x < left + left_tol ||
			n.x > right - right_tol && n.x < right + right_tol)
			n_ids.add(&n_id);
	}
}

void get_bottom_n_ids(Model_T2D_CHM_s &md, double bottom,
	MemoryUtilities::ItemArray<size_t> &n_ids)
{
	n_ids.reset();
	double bottom_tol = abs(bottom) * 1.0e-3;
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Model_T2D_CHM_s::Node &n = md.nodes[n_id];
		if (n.y > bottom - bottom_tol && n.y < bottom + bottom_tol)
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
		if (pcl.y > -0.021 && pcl.y < 0.0)
			pcl_ids.add(&p_id);
	}
}

void get_mid_top_pcl_ids(Model_T2D_CHM_s& md,
	MemoryUtilities::ItemArray<size_t>& pcl_ids)
{
	pcl_ids.reset();
	for (size_t p_id = 0; p_id < md.pcl_num; ++p_id)
	{
		Model_T2D_CHM_s::Particle& pcl = md.pcls[p_id];
		if (pcl.y > -0.011 && pcl.y < 0.0 &&
			pcl.x > -2.5 && pcl.x < 2.5)
			pcl_ids.add(&p_id);
	}
}

void get_left_right_top_pcl_ids(Model_T2D_CHM_s& md,
	MemoryUtilities::ItemArray<size_t>& pcl_ids)
{
	pcl_ids.reset();
	for (size_t p_id = 0; p_id < md.pcl_num; ++p_id)
	{
		Model_T2D_CHM_s::Particle& pcl = md.pcls[p_id];
		if (pcl.y > -0.021 && pcl.y < 0.0 &&
			(pcl.x < -2.5 || pcl.x > 2.5))
			pcl_ids.add(&p_id);
	}
}

};


void test_t2d_mpm_chm_s_t_bar_conference_restart(void)
{
	Model_T2D_CHM_s model;

	using Model_T2D_CHM_s_hdf5_io_utilities::load_chm_s_model_from_hdf5_file;
	load_chm_s_model_from_hdf5_file(
		model,
		"t2d_mpm_chm_t_bar_conference_geo.hdf5",
		"geostatic",
		101
		);

	model.init_bg_mesh(0.05, 0.05);

	//model.set_rigid_circle_velocity(0.0, -0.5, 0.0);
	model.set_rigid_circle_velocity(0.0, -0.05, 0.0);

	MemoryUtilities::ItemArray<size_t> mid_tbc_pcl_ids_mem;
	mid_tbc_pcl_ids_mem.reserve(100);
	get_mid_top_pcl_ids(model, mid_tbc_pcl_ids_mem);
	size_t* mid_tbc_pcl_ids = mid_tbc_pcl_ids_mem.get_mem();
	size_t mid_tbc_num = mid_tbc_pcl_ids_mem.get_num();

	MemoryUtilities::ItemArray<size_t> left_right_tbc_pcl_ids_mem;
	left_right_tbc_pcl_ids_mem.reserve(100);
	get_left_right_top_pcl_ids(model, left_right_tbc_pcl_ids_mem);
	size_t* left_right_tbc_pcl_ids = left_right_tbc_pcl_ids_mem.get_mem();
	size_t left_right_tbc_num = left_right_tbc_pcl_ids_mem.get_num();

	model.init_tys(mid_tbc_num + left_right_tbc_num);
	for (size_t t_id = 0; t_id < mid_tbc_num; ++t_id)
	{
		TractionBC_MPM& tbc = model.tys[t_id];
		tbc.pcl_id = mid_tbc_pcl_ids[t_id];
		tbc.t = 0.02 * -20000.0;
	}
	for (size_t t_id = 0; t_id < left_right_tbc_num; ++t_id)
	{
		TractionBC_MPM& tbc = model.tys[mid_tbc_num + t_id];
		tbc.pcl_id = left_right_tbc_pcl_ids[t_id];
		tbc.t = 0.04 * -20000.0;
	}

	MemoryUtilities::ItemArray<GLfloat> pt_array;
	pt_array.reserve(100);
	GLfloat pt_coord;
	//for (size_t t_id = 0; t_id < mid_tbc_num; ++t_id)
	//{
	//	Model_T2D_CHM_s::Particle &pcl = model.pcls[mid_tbc_pcl_ids[t_id]];
	//	pt_coord = double(pcl.x);
	//	pt_array.add(&pt_coord);
	//	pt_coord = double(pcl.y);
	//	pt_array.add(&pt_coord);
	//	pt_coord = 0.0f;
	//	pt_array.add(&pt_coord);
	//}
	//for (size_t t_id = 0; t_id < left_right_tbc_num; ++t_id)
	//{
	//	Model_T2D_CHM_s::Particle& pcl = model.pcls[left_right_tbc_pcl_ids[t_id]];
	//	pt_coord = double(pcl.x);
	//	pt_array.add(&pt_coord);
	//	pt_coord = double(pcl.y);
	//	pt_array.add(&pt_coord);
	//	pt_coord = 0.0f;
	//	pt_array.add(&pt_coord);
	//}

	// boundary conditions
	MemoryUtilities::ItemArray<size_t> bc_n_ids_mem;
	bc_n_ids_mem.reserve(100);
	size_t *bc_n_ids;

	get_left_and_right_n_ids(model, -3.5, 3.5, bc_n_ids_mem);
	bc_n_ids = bc_n_ids_mem.get_mem();
	model.init_vsxs(bc_n_ids_mem.get_num());
	for (size_t v_id = 0; v_id < model.vsx_num; ++v_id)
	{
		VelocityBC &vbc = model.vsxs[v_id];
		vbc.node_id = bc_n_ids[v_id];
		vbc.v = 0.0;
	}
	//model.init_vfxs(bc_n_ids_mem.get_num());
	//for (size_t v_id = 0; v_id < model.vfx_num; ++v_id)
	//{
	//	VelocityBC &vbc = model.vfxs[v_id];
	//	vbc.node_id = bc_n_ids[v_id];
	//	vbc.v = 0.0;
	//}

	get_bottom_n_ids(model, -5.0, bc_n_ids_mem);
	bc_n_ids = bc_n_ids_mem.get_mem();
	model.init_vsys(bc_n_ids_mem.get_num());
	for (size_t v_id = 0; v_id < model.vsy_num; ++v_id)
	{
		VelocityBC &vbc = model.vsys[v_id];
		vbc.node_id = bc_n_ids[v_id];
		vbc.v = 0.0;
	}
	//model.init_vfys(bc_n_ids_mem.get_num());
	//for (size_t v_id = 0; v_id < model.vfy_num; ++v_id)
	//{
	//	VelocityBC &vbc = model.vfys[v_id];
	//	vbc.node_id = bc_n_ids[v_id];
	//	vbc.v = 0.0;
	//}

	//DisplayModel_T2D disp_model;
	//disp_model.init_win();
	//disp_model.init_model(model);
	//disp_model.init_rigid_circle(model.get_rigid_circle());
	//disp_model.init_points(pt_array.get_mem(), pt_array.get_num() / 3);
	// all
	//disp_model.display(-3.6, 3.6, -5.1, 1.1);
	// left
	//disp_model.display(-3.8, -2.2, -1.0, 1.0);
	// middle
	//disp_model.display(2.3, 2.7, -0.25, 0.25);
	// right
	//disp_model.display(2.2, 3.8, -1.0, 1.0);
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_mpm_chm_t_bar_conference_restart.hdf5");

	// output model
	ModelDataOutput_T2D_CHM_s md("md1");
	md.set_model(model);
	md.set_res_file(res_file_hdf5);
	md.output();

	TimeHistoryOutput_T2D_CHM_s_SE out("penetration");
	out.set_res_file(res_file_hdf5);
	out.set_output_init_state();
	TimeHistoryOutput_ConsoleProgressBar out_pb;

	// penetration step
	Step_T2D_CHM_s_SE step;
	step.set_model(model);
	step.set_time(4.0); // 1.0s
	step.set_dtime(2.0e-6); // 2.0e-7
	out.set_interval_num(1000);
	step.add_time_history(out);
	step.add_time_history(out_pb);
	step.solve();
}

void test_color_animation_t2d_chm_s_t_bar_conference_restart(void)
{
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
	GA_T2D_CHM_s_hdf5 gen(1000, 1000); // window size
	// pore pressure
	gen.init_color_graph(
		0.0,
		40000,
		colors, sizeof(colors) / sizeof(ColorGraph::Colori)
		);
	// strain
	//gen.init_color_graph(
	//	0.0,
	//	1.5,
	//	colors, sizeof(colors) / sizeof(ColorGraph::Colori)
	//);
	gen.generate(5.0,
		//-2.5, 2.5, -1.9, 1.1,
		-3.6, 3.6, -5.1, 1.1,
		"t2d_mpm_chm_t_bar_conference_restart.hdf5",
		"penetration",
		"t2d_mpm_chm_t_bar_conference_restart.gif"
		);
}