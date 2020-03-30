#include "TestsWithGL_pcp.h"

#include "ItemArray.hpp"

#include "test_sim_core.h"

#include "TriangleMesh.h"
#include "TriangleMeshToParticles.h"
#include "Model_T2D_CHM_s.h"

#include "Step_T2D_CHM_s_SE_Geostatic.h"

#include "DisplayModel_T2D.h"
#include "ModelDataOutput_T2D_CHM_s.h"
#include "TimeHistoryOutput_T2D_CHM_s_SE_Geostatic.h"
#include "TimeHistoryOutput_T2D_CHM_s_SE_Geostatic_fub.h"
#include "TimeHistoryOutput_ConsoleProgressBar.h"

#include "test_post_processor.h"

#include "GA_T2D_CHM_s_hdf5.h"

namespace
{

void get_left_and_right_n_ids(
	Model_T2D_CHM_s &md,
	double left,
	double right,
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

void get_bottom_n_ids(
	Model_T2D_CHM_s &md,
	double bottom,
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
};


void test_t2d_mpm_chm_s_t_bar_conference_geo(void)
{
	TriangleMesh tri_mesh;
	tri_mesh.load_mesh("..\\..\\Asset\\rect_t_bar_conference.mesh_data");
	
	Model_T2D_CHM_s model;
	model.init_mesh(tri_mesh);
	model.init_bg_mesh(0.2, 0.2);

	TriangleMeshToParticles mh_2_pcl(tri_mesh);
	mh_2_pcl.replace_with_grid_points(-3.5, 3.5, -3.5, 0.0, 0.03, 0.03);
	mh_2_pcl.replace_with_grid_points(-3.5, 3.5, -5.0, -3.5, 0.06, 0.06);

	// elastic
	//model.init_pcls(mh_2_pcl, 0.3, 2700.0, 1000.0, 1.0e5, 0.3, 5.0e6, 1.0e-12, 1.0e-3);
	//double K = 0.3 / (1.0 - 0.3);
	//for (size_t p_id = 0; p_id < model.pcl_num; ++p_id)
	//{
	//	Model_T2D_CHM_s::Particle &pcl = model.pcls[p_id];
	//	pcl.s22 = -1000.0;
	//	pcl.s11 = pcl.s22 * K;
	//	pcl.s12 = 0.0;
	//}
	// mcc
	model.init_pcls(mh_2_pcl, 0.6, 2650.0, 1000.0, 2.0e6, 1.0e-11, 1.0e-3);
	ModelContainer &mc = model.model_container;
	ModifiedCamClay *cms = mc.add_ModifiedCamClay(model.pcl_num);
	double K = 1.0 - sin(23.5 / 180.0 * 3.14159165359);
	//double ini_stress[6] = { -24267.31, -40361.43, -24267.31, 0.0, 0.0, 0.0 };
	double ini_stress[6] = { -12025.0, -20000.0, -12025.0, 0.0, 0.0, 0.0 };
	for (size_t p_id = 0; p_id < model.pcl_num; ++p_id)
	{
		Model_T2D_CHM_s::Particle &pcl = model.pcls[p_id];
		pcl.s11 = ini_stress[0];
		pcl.s22 = ini_stress[1];
		pcl.s12 = 0.0;
		ModifiedCamClay &mcc = cms[p_id];
		//mcc.set_param_OC(0.3, 0.044, 0.205, 23.5, 3.6677, ini_stress, 20030.8);
		mcc.set_param_NC(0.3, 0.044, 0.205, 23.5, 3.6677, ini_stress);
		pcl.set_cm(mcc);
	}
	std::cout << "pcl_num: " << model.pcl_num << "\n";

	tri_mesh.clear();
	mh_2_pcl.clear();

	// traction force
	MemoryUtilities::ItemArray<size_t> bc_pcl_ids_mem;
	bc_pcl_ids_mem.reserve(100);
	get_top_pcl_ids(model, bc_pcl_ids_mem);
	size_t *tbc_pcl_ids = bc_pcl_ids_mem.get_mem();
	model.init_tys(bc_pcl_ids_mem.get_num());
	for (size_t t_id = 0; t_id < model.ty_num; ++t_id)
	{
		TractionBC_MPM &tbc = model.tys[t_id];
		tbc.pcl_id = tbc_pcl_ids[t_id];
		//tbc.t = 0.03 * -40361.43;
		tbc.t = 0.03 * -20000.0;
	}
	//MemoryUtilities::ItemArray<GLfloat> pt_array;
	//pt_array.reserve(25 * 3);
	//GLfloat pt_coord;
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

	model.init_rigid_circle(0.5, 0.0, 0.5 - 0.014, 0.04);
	//model.set_rigid_circle_velocity(0.0, -0.05, 0.0);
	model.set_contact_stiffness(1.0e5, 1.0e3);
	
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

	get_bottom_n_ids(model, -5.0, bc_n_ids_mem);
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

	// body force
	//model.init_bfys(model.pcl_num);
	//for (size_t p_id = 0; p_id < model.pcl_num; ++p_id)
	//{
	//	BodyForce &bf = model.bfys[p_id];
	//	bf.pcl_id = p_id;
	//	bf.bf = -9.81;
	//}

	//DisplayModel_T2D disp_model;
	//disp_model.init_win();
	//disp_model.init_model(model);
	//disp_model.init_rigid_circle(model.get_rigid_circle());
	////disp_model.init_points(pt_array.get_mem(), pt_array.get_num() / 3);
	//disp_model.display(-3.6, 3.6, -5.1, 1.1);
	////disp_model.display(-1.2, 1.2, -1.2, 0.2);
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_mpm_chm_t_bar_conference_geo.hdf5");

	// output model
	ModelDataOutput_T2D_CHM_s md("md1");
	md.set_model(model);
	md.set_res_file(res_file_hdf5);
	md.output();

	TimeHistoryOutput_T2D_CHM_s_SE_Geostatic out("geostatic");
	out.set_res_file(res_file_hdf5);
	out.set_output_init_state();
	TimeHistoryOutput_T2D_CHM_s_SE_Geostatic_fub out_nf("nodal_force");
	out_nf.set_res_file(res_file_hdf5);
	out_nf.set_output_init_state();
	TimeHistoryOutput_ConsoleProgressBar out_pb;

	// geostatic step
	Step_T2D_CHM_s_SE_Geostatic step_gs;
	step_gs.set_model(model);
	//step_gs.set_mass_scale(10.0, 10.0);
	step_gs.set_time(3.0);
	step_gs.set_dtime(1.0e-5);
	// out
	out.set_interval_num(100);
	step_gs.add_time_history(out);
	// out_nf
	out_nf.set_interval_num(100);
	step_gs.add_time_history(out_nf);
	// out_pb
	step_gs.add_time_history(out_pb);
	step_gs.solve();
	
	//system("pause");
}

void test_color_animation_t2d_chm_s_t_bar_conference_geo(void)
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
	gen.init_color_graph(
		-25000.0,
		-15000.0,
		colors,
		sizeof(colors) / sizeof(ColorGraph::Colori)
		);
	gen.generate(
		5.0, // time
		-3.6, 3.6, -5.1, 1.1,
		"t2d_mpm_chm_t_bar_conference_geo.hdf5",
		"geostatic",
		"t2d_mpm_chm_t_bar_conference_geo.gif"
		);
	//gen.generate(
	//	5.0, // time
	//	-1.2, 1.2, -1.2, 0.2,
	//	"t2d_mpm_chm_t_bar_conference_geo.hdf5",
	//	"geostatic",
	//	"t2d_mpm_chm_t_bar_conference_geo.gif"
	//	);
}