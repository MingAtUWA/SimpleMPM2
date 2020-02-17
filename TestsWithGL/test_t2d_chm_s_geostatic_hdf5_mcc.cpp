#include "TestsWithGL_pcp.h"

#include <cmath>
#include "ItemArray.hpp"

#include "test_sim_core.h"

#include "TriangleMesh.h"
#include "TriangleMeshToParticles.h"
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

void test_t2d_chm_s_geostatic_hdf5_mcc(void)
{
	TriangleMesh tri_mesh;
	tri_mesh.load_mesh("..\\..\\Asset\\rect.mesh_data");

	TriangleMeshToParticles mh_2_pcl(tri_mesh);
	mh_2_pcl.set_even_div_num(2);
	mh_2_pcl.set_generator(TriangleMeshToParticles::GeneratorType::EvenlyDistributedPoint);
	mh_2_pcl.generate_pcls();

	Model_T2D_CHM_s model;
	model.init_mesh(tri_mesh);
	tri_mesh.clear();

	// elasticity
	//model.init_pcls(mh_2_pcl, 0.5, 20.0, 10.0, 1000.0, 0.0, 40000.0, 1.0e-4, 1.0);
	//double K = 0.0;
	//for (size_t p_id = 0; p_id < model.pcl_num; ++p_id)
	//{
	//	Model_T2D_CHM_s::Particle &pcl = model.pcls[p_id];
	//	pcl.s22 = -(1.0 - pcl.y) * 100.0;
	//	pcl.s11 = K * pcl.s22;
	//	pcl.s12 = 0.0;
	//}
	// mcc
	model.init_pcls(mh_2_pcl, 0.5, 20.0, 10.0, 1.0e4, 1.0e-4, 1.0);
	ModifiedCamClay *cms = model.model_container.add_ModifiedCamClay(model.pcl_num);
	double K = 1.0 - sin(25.0 / 180.0 * 3.14159265359);
	double ini_stress[6];
	for (size_t p_id = 0; p_id < model.pcl_num; ++p_id)
	{
		Model_T2D_CHM_s::Particle &pcl = model.pcls[p_id];
		pcl.s22 = -(1.0 - pcl.y) * 100.0;
		pcl.s11 = K * pcl.s22;
		pcl.s12 = 0.0;
		ini_stress[0] = pcl.s11;
		ini_stress[1] = pcl.s22;
		ini_stress[2] = pcl.s11;
		ini_stress[3] = 0.0;
		ini_stress[4] = 0.0;
		ini_stress[5] = 0.0;
		ModifiedCamClay &cm = cms[p_id];
		cm.set_param_NC(0.25, 0.04, 0.15, 25, 1.5, ini_stress);
		pcl.cm = &cm;
		cm.ext_data = &pcl;
	}
	mh_2_pcl.clear();

	model.init_bg_mesh(0.05, 0.05);

	// rigid objects
	model.init_rigid_circle(0.1, 0.1, 1.15, 0.01);
	model.set_rigid_circle_velocity(0.0, 0.0, 0.0);
	model.set_contact_stiffness(100.0, 100.0);

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

	//size_t tbc_pcl_id[] = { 132, 133, 168, 169 };
	//model.init_tys(sizeof(tbc_pcl_id) / sizeof(tbc_pcl_id[0]));
	//for (size_t t_id = 0; t_id < model.ty_num; ++t_id)
	//{
	//	TractionBC_MPM &tbc = model.tys[t_id];
	//	tbc.pcl_id = tbc_pcl_id[t_id];
	//	tbc.t = 0.05 * -100.0;
	//}

	model.init_bfys(model.pcl_num);
	for (size_t p_id = 0; p_id < model.pcl_num; ++p_id)
	{
		BodyForce &bf = model.bfys[p_id];
		bf.pcl_id = p_id;
		bf.bf = -20.0;
	}
	
	ResultFile_XML res_file_xml;
	res_file_xml.init("t2d_chm_s_geostatic_hdf5_mcc.xml");
	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_chm_s_geostatic_hdf5_mcc.hdf5");

	// output model
	ModelDataOutput_T2D_CHM_s md("md1");
	md.set_model(model);
	md.set_res_file(res_file_xml);
	md.output();
	md.set_res_file(res_file_hdf5);
	md.output();

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
	step.set_dtime(1.e-5);
	out2.set_interval_num(100);
	step.add_time_history(out2);
	out3.set_interval_num(100);
	step.add_time_history(out3);
	out_nf.set_interval_num(100);
	step.add_time_history(out_nf);
	step.add_time_history(out4);
	step.solve();

	//system("pause");
}


void test_color_animation_t2d_chm_s_geostatic_hdf5_mcc(void)
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
		{ 0,   255, 46  },
		{ 46,  255, 0   },
		{ 139, 255, 0   },
		{ 232, 255, 0   },
		{ 255, 185, 0   },
		{ 255, 93,  0   },
		{ 255, 0,   0   }
	};
	GA_T2D_CHM_s_hdf5 gen;
	gen.init_color_graph(
		500.0,
		60.0,
		40.0,
		480.0,
		-100,
		0.0,
		colors, 
		sizeof(colors) / sizeof(ColorGraph::Colori)
		);
	gen.generate(5.0, 
		-padding_width, 
		soil_width + padding_width,
		-padding_height,
		soil_height + padding_height,
		"t2d_chm_s_geostatic_hdf5_mcc.hdf5",
		"th3",
		"t2d_chm_s_geostatic_hdf5_mcc.gif"
		);
}
