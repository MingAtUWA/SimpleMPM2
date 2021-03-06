#include "TestsWithGL_pcp.h"

#include "ItemArray.hpp"

#include "test_sim_core.h"

#include "TriangleMesh.h"
#include "TriangleMeshToParticles.h"
#include "AdjustParticlesWithTriangleMesh.hpp"
#include "AdjustParticleSizeWithMesh.hpp"
#include "Model_T2D_CHM_s.h"
#include "Step_T2D_CHM_s_SE_Geostatic.h"

#include "ResultFile_PlainBin.h"
#include "ResultFile_XML.h"
#include "ResultFile_hdf5.h"
#include "Model_T2D_CHM_s_hdf5_io_utilities.h"

#include "ModelDataOutput_T2D_CHM_s.h"
#include "TimeHistoryOutput_T2D_CHM_s_SE_Geostatic.h"
#include "TimeHistoryOutput_T2D_CHM_s_SE_Geostatic_fub.h"
#include "TimeHistoryOutput_ConsoleProgressBar.h"

#include "DisplayModel_T2D.h"

#include "test_post_processor.h"
#include "GA_T2D_CHM_s_hdf5.h"

namespace
{
void get_top_pcl_ids(Model_T2D_CHM_s &md,
	MemoryUtilities::ItemArray<size_t> &pcl_ids)
{
	pcl_ids.reset();
	for (size_t p_id = 0; p_id < md.pcl_num; ++p_id)
	{
		Model_T2D_CHM_s::Particle &pcl = md.pcls[p_id];
		if (pcl.y > 0.985)
			pcl_ids.add(&p_id);
	}
}
};

void test_t2d_chm_s_geostatic_hdf5(void)
{
	TriangleMesh tri_mesh;
	tri_mesh.load_mesh("..\\..\\Asset\\rect.mesh_data");

	Model_T2D_CHM_s model;
	model.init_mesh(tri_mesh);
	model.init_bg_mesh(0.05, 0.05);

	TriangleMeshToParticles mh_2_pcl(tri_mesh);
	//
	//mh_2_pcl.set_even_div_num(2);
	//mh_2_pcl.set_generator(TriangleMeshToParticles::GeneratorType::EvenlyDistributedPoint);
	//mh_2_pcl.generate_pcls();
	//
	//mh_2_pcl.set_generator(TriangleMeshToParticles::GeneratorType::SecondOrderGaussPoint);
	//mh_2_pcl.generate_pcls();
	//
	mh_2_pcl.generate_grid_points(0.0, 0.2, 0.0, 1.0, 0.025, 0.025);
	
	//AdjustParticlesWithTriangleMesh<Model_T2D_CHM_s> ap_mh(mh_2_pcl);
	//ap_mh.distribute_points_area_to_mesh(model, 0.025, 0.025, 0.2, 0.2);

	AdjustParticleSizeWithMesh<Model_T2D_CHM_s> ap_mh(mh_2_pcl);
	//ap_mh.adjust_particles(model, 4, 0.025, 0.025);
	ap_mh.adjust_particles2(model);

	// elasticity
	model.init_pcls(mh_2_pcl, 0.5, 20.0, 10.0, 1000.0, 0.0, 40000.0, 1.0e-4, 1.0);
	for (size_t p_id = 0; p_id < model.pcl_num; ++p_id)
	{
		Model_T2D_CHM_s::Particle &pcl = model.pcls[p_id];
		pcl.s22 = -10;
	}
	// mcc
	//model.init_pcls(mh_2_pcl, );

	tri_mesh.clear();
	mh_2_pcl.clear();

	// rigid objects
	model.init_rigid_circle(0.1, 0.1, 1.15, 0.01);
	model.set_rigid_circle_velocity(0.0, 0.0, 0.0);
	model.set_contact_stiffness(100.0, 100.0);

	MemoryUtilities::ItemArray<size_t> tbc_pcls;
	get_top_pcl_ids(model, tbc_pcls);
	size_t* tbc_pcl_id = tbc_pcls.get_mem();
	model.init_tys(tbc_pcls.get_num());
	for (size_t t_id = 0; t_id < model.ty_num; ++t_id)
	{
		TractionBC_MPM& tbc = model.tys[t_id];
		tbc.pcl_id = tbc_pcl_id[t_id];
		tbc.t = 0.025 * -10.0;
		//tbc.t = 0.05 * -10.0;
	}

	size_t vx_bc_n_id[] = { 0, 3, 15, 16, 17, 18, 19, 20, 21, 22, 23,
							1, 2, 5,  6,  7,  8,  9,  10, 11, 12, 13 };
	model.init_vsxs(sizeof(vx_bc_n_id) / sizeof(vx_bc_n_id[0]));
	for (size_t v_id = 0; v_id < model.vsx_num; ++v_id)
	{
		VelocityBC &vbc = model.vsxs[v_id];
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

	MemoryUtilities::ItemArray<GLfloat> pt_array;
	pt_array.reserve(25 * 3);
	GLfloat pt_coord;
	for (size_t t_id = 0; t_id < model.ty_num; ++t_id)
	{
		Model_T2D_CHM_s::Particle &pcl = model.pcls[model.tys[t_id].pcl_id];
		pt_coord = double(pcl.x);
		pt_array.add(&pt_coord);
		pt_coord = double(pcl.y);
		pt_array.add(&pt_coord);
		pt_coord = 0.0f;
		pt_array.add(&pt_coord);
	}
	
	//model.init_bfys(model.pcl_num);
	//for (size_t p_id = 0; p_id < model.pcl_num; ++p_id)
	//{
	//	BodyForce &bf = model.bfys[p_id];
	//	bf.pcl_id = p_id;
	//	bf.bf = -20.0;
	//}
	
	//DisplayModel_T2D disp_model;
	//disp_model.init_win();
	//disp_model.init_model(model);
	////disp_model.init_rigid_circle(model.get_rigid_circle());
	////disp_model.init_points(pt_array.get_mem(), pt_array.get_num() / 3);
	//disp_model.display(-0.05, 0.25,-0.05, 1.05);
	//return;

	ResultFile_XML res_file_xml;
	res_file_xml.init("t2d_chm_s_geostatic_hdf5.xml");
	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_chm_s_geostatic_hdf5.hdf5");

	// output model
	ModelDataOutput_T2D_CHM_s md("md1");
	md.set_model(model);
	//md.set_res_file(res_file_pb);
	//md.output();
	md.set_res_file(res_file_xml);
	md.output();
	md.set_res_file(res_file_hdf5);
	md.output();

	//TimeHistoryOutput_T2D_CHM_s_SE_Geostatic out1("th1");
	//out1.set_res_file(res_file_pb);
	//out1.set_output_init_state();
	TimeHistoryOutput_T2D_CHM_s_SE_Geostatic out2("th2");
	out2.set_res_file(res_file_xml);
	out2.set_output_init_state();
	TimeHistoryOutput_T2D_CHM_s_SE_Geostatic out3("th3");
	out3.set_res_file(res_file_hdf5);
	out3.set_output_init_state();
	TimeHistoryOutput_T2D_CHM_s_SE_Geostatic_fub out_nf("nodal_force");
	out_nf.set_res_file(res_file_hdf5);
	out_nf.set_output_init_state();
	TimeHistoryOutput_ConsoleProgressBar out4;

	Step_T2D_CHM_s_SE_Geostatic step;
	step.set_model(model);
	step.set_time(1.0);
	step.set_dtime(1.0e-4);
	//out1.set_interval_num(100);
	//step.add_time_history(out1);
	out2.set_interval_num(100);
	step.add_time_history(out2);
	out3.set_interval_num(100);
	step.add_time_history(out3);
	out_nf.set_interval_num(100);
	step.add_time_history(out_nf);
	step.add_time_history(out4);
	step.solve();
}

void test_color_animation_t2d_chm_s_geostatic_hdf5(void)
{
	double soil_height = 1.0;
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
		-11.0,
		-9.0,
		colors,
		sizeof(colors) / sizeof(ColorGraph::Colori)
		);
	gen.generate(5.0,
		-padding_width,
		soil_width + padding_width,
		-padding_height,
		soil_height + padding_height,
		"t2d_chm_s_geostatic_hdf5.hdf5",
		"th3",
		"t2d_chm_s_geostatic_hdf5.gif"
		);
}
