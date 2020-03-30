#include "TestsWithGL_pcp.h"

#include "Model_T2D_fluid.h"
#include "Step_T2D_fluid.h"
#include "ModelDataOutput_T2D_fluid.h"
#include "TimeHistoryOutput_T2D_fluid.h"
#include "TimeHistoryOutput_ConsoleProgressBar.h"
#include "DisplayModel_T2D.h"
#include "test_sim_core.h"

#include "GA_T2D_fluid_hdf5.h"
#include "test_post_processor.h"

namespace
{
void get_left_and_right_n_ids(
	Model_T2D_fluid &md,
	double left,
	double right,
	MemoryUtilities::ItemArray<size_t> &n_ids)
{
	n_ids.reset();
	double left_tol = abs(right - left) * 1.0e-3;
	double right_tol = left_tol;
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Model_T2D_fluid::Node &n = md.nodes[n_id];
		if (n.x > left - left_tol   && n.x < left + left_tol ||
			n.x > right - right_tol && n.x < right + right_tol)
			n_ids.add(&n_id);
	}
}

void get_top_and_bottom_n_ids(
	Model_T2D_fluid &md,
	double top,
	double bottom,
	MemoryUtilities::ItemArray<size_t> &n_ids)
{
	n_ids.reset();
	double top_tol = abs(top - bottom) * 1.0e-3;
	double bottom_tol = top_tol;
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Model_T2D_fluid::Node &n = md.nodes[n_id];
		if (n.y > top    - top_tol    && n.y < top    + top_tol ||
			n.y > bottom - bottom_tol && n.y < bottom + bottom_tol)
			n_ids.add(&n_id);
	}
}
}

void test_t2d_fluid(void)
{
	TriangleMesh tri_mesh;
	tri_mesh.load_mesh("..\\..\\Asset\\rect_fluid.mesh_data");

	Model_T2D_fluid model;
	model.init_mesh(tri_mesh); // calculation mesh
	model.init_bg_mesh(0.1, 0.1); // background mesh
	std::cout << model.node_num << " " << model.elem_num << "\n";

	TriangleMeshToParticles mh_2_pcl(tri_mesh);
	mh_2_pcl.replace_with_grid_points(0.0, 1.0, 0.0, 2.0, 0.075, 0.075);
	model.init_pcls(mh_2_pcl, 10.0, 1.0e-3, 0.0, 100.0);
	
	tri_mesh.clear();
	mh_2_pcl.clear();

	model.init_bfys(model.pcl_num);
	for (size_t pcl_id = 0; pcl_id < model.pcl_num; ++pcl_id)
	{
		BodyForce &bf = model.bfys[pcl_id];
		bf.bf = -1.0;
		bf.pcl_id = pcl_id;
	}

	MemoryUtilities::ItemArray<size_t> bc_n_ids_mem;
	bc_n_ids_mem.reserve(100);
	size_t *bc_n_ids;
	//MemoryUtilities::ItemArray<GLfloat> pt_array;
	//pt_array.reserve(29 * 3);
	//GLfloat pt_coord;
	
	get_left_and_right_n_ids(model, 0.0, 4.0, bc_n_ids_mem);
	bc_n_ids = bc_n_ids_mem.get_mem();
	model.init_vxs(bc_n_ids_mem.get_num());
	for (size_t v_id = 0; v_id < model.vx_num; ++v_id)
	{
		VelocityBC &vbc = model.vxs[v_id];
		vbc.node_id = bc_n_ids[v_id];
		vbc.v = 0.0;
	}
	//for (size_t n_id = 0; n_id < model.vx_num; ++n_id)
	//{
	//	Model_T2D_fluid::Node &n = model.nodes[bc_n_ids[n_id]];
	//	pt_coord = double(n.x);
	//	pt_array.add(&pt_coord);
	//	pt_coord = double(n.y);
	//	pt_array.add(&pt_coord);
	//	pt_coord = 0.0f;
	//	pt_array.add(&pt_coord);
	//}
	
	get_top_and_bottom_n_ids(model, 0.0, 4.0, bc_n_ids_mem);
	bc_n_ids = bc_n_ids_mem.get_mem();
	model.init_vys(bc_n_ids_mem.get_num());
	for (size_t v_id = 0; v_id < model.vy_num; ++v_id)
	{
		VelocityBC &vbc = model.vys[v_id];
		vbc.node_id = bc_n_ids[v_id];
		vbc.v = 0.0;
	}
	//for (size_t n_id = 0; n_id < model.vy_num; ++n_id)
	//{
	//	Model_T2D_fluid::Node &n = model.nodes[bc_n_ids[n_id]];
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
	//disp_model.display(-0.5, 4.5, -0.5, 4.5);
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_fluid.hdf5");

	// output model
	ModelDataOutput_T2D_fluid md("fluid_md");
	md.set_model(model);
	md.set_res_file(res_file_hdf5);
	md.output();

	TimeHistoryOutput_T2D_fluid out1("flow");
	out1.set_res_file(res_file_hdf5);
	out1.set_output_init_state();
	out1.set_interval_num(100);
	TimeHistoryOutput_ConsoleProgressBar out2;

	Step_T2D_fluid step;
	step.set_model(model);
	step.set_time(5.0);
	step.set_dtime(5.0e-5);
	step.add_time_history(out1);
	step.add_time_history(out2);
	step.solve();

	//system("pause");
}

void test_postprocess_t2d_fluid(void)
{
	double soil_height = 4.0;
	double soil_width = 4.0;
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
	GA_T2D_fluid_hdf5 gen;
	gen.init_color_graph(
		-10.0, 1.0,
		colors, sizeof(colors) / sizeof(ColorGraph::Colori)
		);
	gen.generate(
		5.0,
		-padding_width,  soil_width  + padding_width,
		-padding_height, soil_height + padding_height,
		"t2d_fluid.hdf5",
		"flow",
		"t2d_fluid.gif"
		);
}
