#include "TestsWithGL_pcp.h"

#include "ItemArray.hpp"

#include "Model_S2D_ME_s.h"
#include "Step_S2D_ME_s_Geostatic.h"
#include "ResultFile_hdf5.h"
#include "ModelDataOutput_S2D_ME_s.h"
#include "TimeHistoryOutput_S2D_ME_s_Geostatic.h"
#include "TimeHistoryOutput_ConsoleProgressBar.h"

#include "DisplayModel_S2D.h"

#include "test_sim_core.h"

#include "GA_S2D_ME_s_hdf5.h"
#include "test_post_processor.h"

void test_s2d_mpm_me_s_geostatic(void)
{
	Model_S2D_ME_s model;
	
	size_t elem_x_num = 3;
	size_t elem_y_num = 10;
	model.init_mesh(0.0, 0.0, 1.5, 5.0, elem_x_num, elem_y_num);

	// pcl doesn't overlaps the boundary
	//size_t pcl_x_num = elem_x_num * 2;
	//size_t pcl_y_num = elem_y_num * 2;
	// pcl overlaps the boundary
	size_t pcl_x_num = 7;
	size_t pcl_y_num = 21;
	model.init_pcls(0.0, 0.0, 1.5, 5.0, pcl_x_num, pcl_y_num, 10.0);

	model.init_mat_model_le(100.0, 0.0);

	model.init_tys(pcl_x_num);
	for (size_t t_id = 0; t_id < pcl_x_num; ++t_id)
	{
		TractionBC_MPM &tbc = model.get_tys()[t_id];
		tbc.pcl_id = (pcl_y_num - 1) * pcl_x_num + t_id;
		tbc.t = 1.5 / double(pcl_x_num) * -10.0;
	}
	for (size_t p_id = 0; p_id < model.get_pcl_num(); ++p_id)
	{
		Model_S2D_ME_s::Particle &pcl = model.get_pcls()[p_id];
		pcl.s22 = -10.0;
	}

	size_t node_x_num = elem_x_num + 1;
	size_t node_y_num = elem_y_num + 1;
	model.init_vxs(node_y_num * 2);
	for (size_t v_id = 0; v_id < node_y_num; ++v_id)
	{
		VelocityBC &vbc1 = model.get_vxs()[v_id];
		vbc1.node_id = v_id * node_x_num;
		vbc1.v = 0.0;
		VelocityBC &vbc2 = model.get_vxs()[v_id + node_y_num];
		vbc2.node_id = v_id * node_x_num + node_x_num - 1;
		vbc2.v = 0.0;
	}
	model.init_vys(node_x_num);
	for (size_t v_id = 0; v_id < node_x_num; ++v_id)
	{
		VelocityBC &vbc = model.get_vys()[v_id];
		vbc.node_id = v_id;
		vbc.v = 0.0;
	}

	MemoryUtilities::ItemArray<GLfloat> pt_array;
	pt_array.reserve(100);
	GLfloat pt_coord;
	//for (size_t v_id = 0; v_id < model.get_vx_num(); ++v_id)
	//{
	//	VelocityBC &vbc = model.get_vxs()[v_id];
	//	Model_S2D_ME_s::Node &n = model.get_nodes()[vbc.node_id];
	//	pt_coord = GLfloat(model.get_node_x(n));
	//	pt_array.add(&pt_coord);
	//	pt_coord = GLfloat(model.get_node_y(n));
	//	pt_array.add(&pt_coord);
	//	pt_coord = 0.0f;
	//	pt_array.add(&pt_coord);
	//}
	//for (size_t v_id = 0; v_id < model.get_vy_num(); ++v_id)
	//{
	//	VelocityBC &vbc = model.get_vys()[v_id];
	//	Model_S2D_ME_s::Node &n = model.get_nodes()[vbc.node_id];
	//	pt_coord = GLfloat(model.get_node_x(n));
	//	pt_array.add(&pt_coord);
	//	pt_coord = GLfloat(model.get_node_y(n));
	//	pt_array.add(&pt_coord);
	//	pt_coord = 0.0f;
	//	pt_array.add(&pt_coord);
	//}

	//for (size_t t_id = 0; t_id < model.get_ty_num(); ++t_id)
	//{
	//	TractionBC_MPM &tbc = model.get_tys()[t_id];
	//	Model_S2D_ME_s::Particle &pcl = model.get_pcls()[tbc.pcl_id];
	//	pt_coord = double(pcl.x);
	//	pt_array.add(&pt_coord);
	//	pt_coord = double(pcl.y);
	//	pt_array.add(&pt_coord);
	//	pt_coord = 0.0f;
	//	pt_array.add(&pt_coord);
	//}

	DisplayModel_S2D disp_model;
	disp_model.init_win();
	disp_model.init_model(model);
	disp_model.init_points(pt_array.get_mem(), pt_array.get_num()/3);
	disp_model.display(-0.1, 1.6, -0.1, 5.1);
	return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("s2d_me_s_geostatic.h5");

	// output model
	ModelDataOutput_S2D_ME_s md("md_out");
	md.set_model(model);
	md.set_res_file(res_file_hdf5);
	md.output();

	TimeHistoryOutput_S2D_ME_s_Geostatic out1("geostatic");
	out1.set_res_file(res_file_hdf5);
	out1.set_interval_num(100);
	out1.set_output_init_state();

	TimeHistoryOutput_ConsoleProgressBar out3;

	Step_S2D_ME_s_Geostatic step;
	step.use_volume_enhancement();
	step.set_model(model);
	step.set_time(5.0);
	step.set_dtime(1.0e-4);
	step.add_time_history(out1);
	step.add_time_history(out3);
	step.solve();
}

void display_s2d_mpm_me_s_geostatic()
{
	double soil_height = 5.0;
	double soil_width = 1.5;
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
	GA_S2D_ME_s_hdf5 gen(500, 1000);
	gen.init_color_graph(
		-12.0,
		-8.0,
		colors,
		sizeof(colors) / sizeof(ColorGraph::Colori)
		);
	gen.generate(5.0, 
		-padding_width,
		soil_width + padding_width,
		-padding_height,
		soil_height + padding_height,
		"s2d_me_s_geostatic.h5",
		"geostatic",
		"s2d_me_s_geostatic.gif"
		);
}
