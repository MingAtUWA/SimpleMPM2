#include "TestsWithGL_pcp.h"

#include "ItemArray.hpp"

#include "test_sim_core.h"

#include "TriangleMesh.h"
#include "TriangleMeshToParticles.h"
#include "Model_T2D_CHM_s.h"

#include "Step_T2D_CHM_s_SE.h"
#include "Step_T2D_CHM_s_SE_Geostatic.h"

#include "DisplayModel_T2D.h"
#include "ModelDataOutput_T2D_CHM_s.h"
#include "TimeHistoryOutput_T2D_CHM_s_SE.h"
#include "TimeHistoryOutput_T2D_CHM_s_SE_Geostatic.h"
#include "TimeHistoryOutput_ConsoleProgressBar.h"

#include "ResultFile_hdf5.h"

#include "Model_T2D_CHM_s_hdf5_io_utilities.h"

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

//// find small pcls in the middle
//void get_top_small_pcls(Model_T2D_CHM_s &md,
//	MemoryUtilities::ItemArray<size_t> &p_ids)
//{
//	p_ids.reset();
//	for (size_t p_id = 0; p_id < md.pcl_num; ++p_id)
//	{
//		Model_T2D_CHM_s::Particle &pcl = md.pcls[p_id];
//		if (pcl.y > -0.011 && pcl.x > -2.0 && pcl.x < 2.0)
//			p_ids.add(&p_id);
//	}
//}
//
//// find big pcls in the centre
//void get_top_big_pcls(Model_T2D_CHM_s &md,
//	MemoryUtilities::ItemArray<size_t> &p_ids)
//{
//	p_ids.reset();
//	for (size_t p_id = 0; p_id < md.pcl_num; ++p_id)
//	{
//		Model_T2D_CHM_s::Particle &pcl = md.pcls[p_id];
//		if (pcl.y > -0.041 && (pcl.x < -2.0 || pcl.x > 2.0))
//			p_ids.add(&p_id);
//	}
//}
};


void test_t2d_mpm_chm_s_t_bar_real_geostatic(void)
{
	TriangleMesh tri_mesh;
	tri_mesh.load_mesh("..\\..\\Asset\\rect_t_bar_real.mesh_data");
	//std::cout << "node num: " << tri_mesh.get_node_num() << "\n"
	//		  << "elem num: " << tri_mesh.get_elem_num() << "\n";
	
	Model_T2D_CHM_s model;
	model.init_mesh(tri_mesh);
	tri_mesh.clear();

	model.init_bg_mesh(0.1, 0.1);

	TriangleMeshToParticles mh_2_pcl(tri_mesh);
	//mh_2_pcl.set_even_div_num(2);
	mh_2_pcl.set_generator(TriangleMeshToParticles::GeneratorType::SecondOrderGaussPoint);
	mh_2_pcl.generate_pcls();
	mh_2_pcl.clear_points_in_rect(-3.0, 3.0, 0.0, 0.3);
	mh_2_pcl.replace_with_grid_points(-2.0, 2.0, -2.0, 0.0, 0.02, 0.02);
	//std::cout << "pcl num: " << mh_2_pcl.get_pcl_num() << "\n";

	model.init_pcls(mh_2_pcl, 0.706, 2700.0, 1000.0, 20.0e6, 0.3, 50.0e7, 1.0e-9, 1.0e-3);
	mh_2_pcl.clear();

	ModelContainer &mc = model.model_container;
	LinearElasticity *cms = mc.add_LinearElasticity(model.pcl_num);
	double K = 0.3 / (1.0 - 0.3);
	for (size_t p_id = 0; p_id < model.pcl_num; ++p_id)
	{
		Model_T2D_CHM_s::Particle &pcl = model.pcls[p_id];
		// init geostatic stress
		pcl.s22 = (1.0 - pcl.y) * (2700.0 - 1000.0) * (1.0 - 0.706) * -9.81;
		pcl.s11 = K * pcl.s22;
		pcl.s12 = 0.0;
		// init material model
		LinearElasticity &cm = cms[p_id];
		cm.set_param(20.0e6, 0.3);
		pcl.set_cm(cm);
	}

	model.init_rigid_circle(0.25, 0.0, 0.25, 0.025);
	//model.set_rigid_circle_velocity(0.0, -0.05, 0.0);
	model.set_contact_stiffness(200.0e6, 200.0e6);
	
	MemoryUtilities::ItemArray<size_t> bc_n_ids_mem;
	bc_n_ids_mem.reserve(100);
	size_t *bc_n_ids;

	get_left_and_right_n_ids(model, -3.0, 3.0, bc_n_ids_mem);
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

	get_bottom_n_ids(model, -3.5, bc_n_ids_mem);
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
	model.init_bfys(model.pcl_num);
	for (size_t p_id = 0; p_id < model.pcl_num; ++p_id)
	{
		BodyForce &bf = model.bfys[p_id];
		bf.pcl_id = p_id;
		bf.bf = -9.81;
	}

	//// traction force
	//// small pcls
	//MemoryUtilities::ItemArray<size_t> tbc_sp_ids_mem;
	//get_top_small_pcls(model, tbc_sp_ids_mem);
	//size_t sp_num = tbc_sp_ids_mem.get_num();
	//// body force
	//MemoryUtilities::ItemArray<size_t> tbc_bp_ids_mem;
	//get_top_big_pcls(model, tbc_bp_ids_mem);
	//size_t bp_num = tbc_bp_ids_mem.get_num();
	////
	//model.init_tys(sp_num + bp_num);
	//for (size_t t_id = 0; t_id < sp_num; ++t_id)
	//{
	//	TractionBC_MPM &tbc = model.tys[t_id];
	//	tbc.pcl_id = tbc_sp_ids_mem[t_id];
	//	tbc.t = 0.02 * -1000.0;
	//}
	//for (size_t t_id = 0; t_id < bp_num; ++t_id)
	//{
	//	TractionBC_MPM &tbc = model.tys[sp_num + t_id];
	//	tbc.pcl_id = tbc_bp_ids_mem[t_id];
	//	tbc.t = 0.1 * -1000.0;
	//}

	//for (size_t t_id = 0; t_id < sp_num; ++t_id)
	//{
	//	Model_T2D_CHM_s::Particle &pcl = model.pcls[tbc_sp_ids_mem[t_id]];
	//	pt_coord = double(pcl.x);
	//	pt_array.add(&pt_coord);
	//	pt_coord = double(pcl.y);
	//	pt_array.add(&pt_coord);
	//	pt_coord = 0.0f;
	//	pt_array.add(&pt_coord);
	//}
	//for (size_t t_id = 0; t_id < bp_num; ++t_id)
	//{
	//	Model_T2D_CHM_s::Particle &pcl = model.pcls[tbc_bp_ids_mem[t_id]];
	//	pt_coord = double(pcl.x);
	//	pt_array.add(&pt_coord);
	//	pt_coord = double(pcl.y);
	//	pt_array.add(&pt_coord);
	//	pt_coord = 0.0f;
	//	pt_array.add(&pt_coord);
	//}
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
	//disp_model.init_points(pt_array.get_mem(), pt_array.get_num() / 3);
	//disp_model.display(-3.2, 3.2, -3.7, 0.5);
	////disp_model.display(-1.2, 1.2, -1.2, 0.2);
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_mpm_chm_t_bar_real_geostatic.hdf5");

	// output model
	ModelDataOutput_T2D_CHM_s md("md1");
	md.set_model(model);
	md.set_res_file(res_file_hdf5);
	md.output();

	TimeHistoryOutput_T2D_CHM_s_SE_Geostatic out("geostatic");
	out.set_res_file(res_file_hdf5);
	out.set_output_init_state();
	TimeHistoryOutput_ConsoleProgressBar out_pb;

	// geostatic step
	Step_T2D_CHM_s_SE_Geostatic step_gs;
	step_gs.set_model(model);
	step_gs.set_mass_scale(10.0, 10.0);
	step_gs.set_time(1.0);
	step_gs.set_dtime(1.0e-6);
	// out
	out.set_interval_num(100);
	step_gs.add_time_history(out);
	// out_pb
	step_gs.add_time_history(out_pb);
	step_gs.solve();
	
	//system("pause");
}

void test_color_animation_t2d_chm_s_t_bar_real_geostatic(void)
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
		-5.0e3,
		5.0e3,
		colors,
		sizeof(colors) / sizeof(ColorGraph::Colori));
	gen.generate(5.0, -3.2, 3.2, -3.7, 0.5,
		"t2d_mpm_chm_t_bar_real_geostatic.hdf5",
		"geostatic",
		"t2d_mpm_chm_t_bar_real_geostatic.gif");
}