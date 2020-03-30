#include "TestsWithGL_pcp.h"

#include "ItemArray.hpp"
#include "TriangleMesh.h"
#include "TriangleMeshToParticles.h"
#include "Step_T2D_CHM_s_SE_Geostatic.h"
#include "ModelDataOutput_T2D_CHM_s.h"
#include "TimeHistoryOutput_T2D_CHM_s_SE.h"
#include "TimeHistoryOutput_ConsoleProgressBar.h"
#include "DisplayModel_T2D.h"

#include "test_sim_core.h"

#include "test_post_processor.h"

#include "GA_T2D_CHM_s_hdf5.h"

namespace
{
	void get_left_and_right_n_ids(
		Model_T2D_CHM_s &md,
		double left,
		double right,
		MemoryUtilities::ItemArray<size_t> &n_ids
		)
	{
		n_ids.reset();
		double left_tol = abs(right - left) * 1.0e-3;
		double right_tol = left_tol;
		for (size_t n_id = 0; n_id < md.node_num; ++n_id)
		{
			Model_T2D_CHM_s::Node &n = md.nodes[n_id];
			if (n.x > left - left_tol   && n.x < left + left_tol ||
				n.x > right - right_tol && n.x < right + right_tol)
				n_ids.add(&n_id);
		}
	}

	void get_bottom_n_ids(
		Model_T2D_CHM_s &md,
		double bottom,
		MemoryUtilities::ItemArray<size_t> &n_ids
		)
	{
		n_ids.reset();
		double bottom_tol = 1.0e-3;
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
			if (pcl.y > 2.0 - 0.011 && pcl.y < 2.0)
				pcl_ids.add(&p_id);
		}
	}
};

void test_t2d_chm_geostatic()
{
	TriangleMesh tri_mesh;
	tri_mesh.load_mesh("..\\..\\Asset\\rect_test_geostatic.mesh_data");

	TriangleMeshToParticles mh_2_pcl(tri_mesh);
	mh_2_pcl.generate_grid_points(0.0, 1.5, 0.4, 2.0, 0.02, 0.02);
	mh_2_pcl.generate_grid_points(0.0, 1.5, 0.0, 0.4, 0.04, 0.04);

	Model_T2D_CHM_s model;
	model.init_mesh(tri_mesh);
	tri_mesh.clear();

	model.init_pcls(mh_2_pcl, 0.5, 20.0, 10.0, 1000.0, 0.0, 40000.0, 1.0e-4, 1.0);
	double K = 0.0;
	for (size_t p_id = 0; p_id < model.pcl_num; ++p_id)
	{
		Model_T2D_CHM_s::Particle &pcl = model.pcls[p_id];
		pcl.s22 = -10.0;
		pcl.s11 = K * pcl.s22;
		pcl.s12 = 0.0;
	}

	MemoryUtilities::ItemArray<GLfloat> pt_array;
	pt_array.reserve(100);
	GLfloat pt_coord;

	MemoryUtilities::ItemArray<size_t> pcl_ids_mem;
	get_top_pcl_ids(model, pcl_ids_mem);
	size_t *pcl_ids = pcl_ids_mem.get_mem();
	model.init_tys(pcl_ids_mem.get_num());
	for (size_t t_id = 0; t_id < model.ty_num; ++t_id)
	{
		TractionBC_MPM &tbc = model.tys[t_id];
		tbc.pcl_id = pcl_ids[t_id];
		tbc.t = 0.02 * -10.0;
	}
	//for (size_t t_id = 0; t_id < model.ty_num; ++t_id)
	//{
	//	Model_T2D_CHM_s::Particle &pcl = model.pcls[model.tys[t_id].pcl_id];
	//	pt_coord = double(pcl.x);
	//	pt_array.add(&pt_coord);
	//	pt_coord = double(pcl.y);
	//	pt_array.add(&pt_coord);
	//	pt_coord = 0.0f;
	//	pt_array.add(&pt_coord);
	//}

	MemoryUtilities::ItemArray<size_t> bc_n_ids_mem;
	bc_n_ids_mem.reserve(100);
	size_t *bc_n_ids;

	get_left_and_right_n_ids(model, 0.0, 1.5, bc_n_ids_mem);
	bc_n_ids = bc_n_ids_mem.get_mem();
	model.init_vsxs(bc_n_ids_mem.get_num());
	for (size_t v_id = 0; v_id < model.vsx_num; ++v_id)
	{
		VelocityBC &vbc = model.vsxs[v_id];
		vbc.node_id = bc_n_ids[v_id];
		vbc.v = 0.0;
	}
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

	get_bottom_n_ids(model, 0.0, bc_n_ids_mem);
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

	//DisplayModel_T2D disp_model;
	//disp_model.init_win();
	//disp_model.init_model(model);
	//disp_model.init_points(pt_array.get_mem(), pt_array.get_num() / 3);
	//disp_model.display(-0.05, 1.55, -0.05, 2.05);
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_mpm_chm_geostatic.hdf5");

	ModelDataOutput_T2D_CHM_s md("md1");
	md.set_model(model);
	md.set_res_file(res_file_hdf5);
	md.output();

	TimeHistoryOutput_T2D_CHM_s_SE out("geostatic");
	out.set_res_file(res_file_hdf5);
	out.set_output_init_state();
	TimeHistoryOutput_ConsoleProgressBar out_pb;

	Step_T2D_CHM_s_SE_Geostatic step;
	step.set_model(model);
	step.set_time(1.0);
	step.set_dtime(2.0e-7);
	out.set_interval_num(1000);
	step.add_time_history(out);
	step.add_time_history(out_pb);
	step.solve();
}

void test_t2d_chm_geostatic_animation(void)
{
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
	GA_T2D_CHM_s_hdf5 gen(1000, 1000); // window size
	gen.init_color_graph(
		-10.0,
		0.0,
		colors,
		sizeof(colors) / sizeof(ColorGraph::Colori)
		);
	gen.generate(
		5.0, // time
		-0.05, 1.55, -0.05, 2.05,
		"t2d_mpm_chm_geostatic.hdf5",
		"geostatic",
		"t2d_mpm_chm_geostatic.gif"
		);
}