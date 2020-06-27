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

void test_t2d_mpm_chm_s_t_bar_conference_restart2()
{
	Model_T2D_CHM_s model;

	using Model_T2D_CHM_s_hdf5_io_utilities::load_chm_s_model_from_hdf5_file;
	load_chm_s_model_from_hdf5_file(
		model,
		"t2d_mpm_chm_t_bar_conference_restart1.h5",
		"penetration",
		501
		);

	//MemoryUtilities::ItemArray<GLfloat> pt_array;
	//pt_array.reserve(100);
	//GLfloat pt_coord;
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

	//DisplayModel_T2D disp_model;
	//disp_model.init_win();
	//disp_model.init_model(model);
	//disp_model.init_rigid_circle(model.get_rigid_circle());
	//disp_model.init_points(pt_array.get_mem(), pt_array.get_num() / 3);
	// all
	//disp_model.display(-3.6, 3.6, -5.1, 1.1);

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_mpm_chm_t_bar_conference_restart2.h5");

	ModelDataOutput_T2D_CHM_s md("md1");
	md.set_model(model);
	md.set_res_file(res_file_hdf5);
	md.output();

	TimeHistoryOutput_T2D_CHM_s_SE out("penetration");
	out.set_res_file(res_file_hdf5);
	out.set_output_init_state();
	TimeHistoryOutput_ConsoleProgressBar out_pb;

	Step_T2D_CHM_s_SE step;
	step.set_model(model);
	step.set_time(5.0); // 1.0s
	step.set_dtime(2.0e-6); // 2.0e-7
	out.set_interval_num(500);
	step.add_time_history(out);
	step.add_time_history(out_pb);
	step.solve();
}

void test_color_animation_t2d_chm_s_t_bar_conference_restart2()
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
		"t2d_mpm_chm_t_bar_conference_restart2.h5",
		"penetration",
		"t2d_mpm_chm_t_bar_conference_restart2.gif"
		);
}