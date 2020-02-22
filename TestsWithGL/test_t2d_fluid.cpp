#include "TestsWithGL_pcp.h"

#include "Model_T2D_fluid.h"
#include "Step_T2D_fluid.h"
#include "ModelDataOutput_T2D_fluid.h"
#include "TimeHistoryOutput_T2D_fluid.h"
#include "TimeHistoryOutput_ConsoleProgressBar.h"
#include "test_sim_core.h"

#include "GA_T2D_fluid_hdf5.h"
#include "test_post_processor.h"

void test_t2d_fluid(void)
{
	TriangleMesh tri_mesh;
	tri_mesh.load_mesh("..\\..\\Asset\\square.mesh_data");

	Model_T2D_fluid model;
	model.init_mesh(tri_mesh); // calculation mesh
	model.init_bg_mesh(0.05, 0.05); // background mesh

	TriangleMeshToParticles mh_2_pcl(tri_mesh);
	mh_2_pcl.set_even_div_num(2);
	mh_2_pcl.set_generator(TriangleMeshToParticles::GeneratorType::EvenlyDistributedPoint);
	mh_2_pcl.generate_pcls();
	model.init_pcls(mh_2_pcl, 20.0, 1.0, 1.0, 100.0);
	
	tri_mesh.clear();
	mh_2_pcl.clear();

	size_t vx_bc_n_id[] = { 0, 3, 15, 16, 17, 18, 19, 20, 21, 22, 23,
							1, 2,  5,  6,  7,  8,  9, 10, 11, 12, 13 };
	model.init_vxs(sizeof(vx_bc_n_id) / sizeof(vx_bc_n_id[0]));
	for (size_t v_id = 0; v_id < model.vx_num; ++v_id)
	{
		VelocityBC &vbc = model.vxs[v_id];
		vbc.node_id = vx_bc_n_id[v_id];
		vbc.v = 0.0;
	}

	size_t vy_bc_n_id[] = { 0, 1, 4 };
	model.init_vys(sizeof(vy_bc_n_id) / sizeof(vy_bc_n_id[0]));
	for (size_t v_id = 0; v_id < model.vy_num; ++v_id)
	{
		VelocityBC &vbc = model.vys[v_id];
		vbc.node_id = vy_bc_n_id[v_id];
		vbc.v = 0.0;
	}

	//MemoryUtilities::ItemArray<GLfloat> pt_array;
	//pt_array.reserve(29 * 3);
	//GLfloat pt_coord;
	//for (size_t n_id = 0; n_id < sizeof(vx_bc_n_id)/sizeof(vx_bc_n_id[0]); ++n_id)
	//{
	//	Model_T2D_ME_s::Node &n = model.nodes[vx_bc_n_id[n_id]];
	//	pt_coord = double(n.x);
	//	pt_array.add(&pt_coord);
	//	pt_coord = double(n.y);
	//	pt_array.add(&pt_coord);
	//	pt_coord = 0.0f;
	//	pt_array.add(&pt_coord);
	//}
	//for (size_t n_id = 0; n_id < sizeof(vy_bc_n_id) / sizeof(vy_bc_n_id[0]); ++n_id)
	//{
	//	Model_T2D_ME_s::Node &n = model.nodes[vy_bc_n_id[n_id]];
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
	////disp_model.init_points(pt_array.get_mem(), pt_array.get_num()/3);
	//disp_model.display(-0.05, 0.25, -0.05, 1.05);
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_fluid.hdf5");

	// output model
	ModelDataOutput_T2D_fluid md("fluid_md");
	md.set_model(model);
	md.set_res_file(res_file_hdf5);
	md.output();

	TimeHistoryOutput_T2D_fluid out1("result");
	out1.set_res_file(res_file_hdf5);
	out1.set_output_init_state();
	out1.set_interval_num(100);
	TimeHistoryOutput_ConsoleProgressBar out2;

	Step_T2D_fluid step;
	step.set_model(model);
	step.set_time(1.0);
	step.set_dtime(1.0e-5);
	step.add_time_history(out1);
	step.add_time_history(out2);
	step.solve();

	//system("pause");
}

void test_postprocess_t2d_fluid(void)
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
	GA_T2D_fluid_hdf5 gen;
	gen.init_color_graph(
		-1.0,
		1.0,
		colors,
		sizeof(colors) / sizeof(ColorGraph::Colori)
		);
	gen.generate(
		5.0,
		-padding_width,
		soil_width + padding_width,
		-padding_height,
		soil_height + padding_height,
		"t2d_fluid.hdf5",
		"result",
		"t2d_fluid.gif"
		);
}
