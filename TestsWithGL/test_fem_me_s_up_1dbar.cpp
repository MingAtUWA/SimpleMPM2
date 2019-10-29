#include "TestsWithGL_pcp.h"

#include "Model_S2D_ME_s_FEM_up.h"
#include "Step_S2D_ME_s_FEM_up.h"
#include "ModelDataOutput_S2D_ME_s_FEM_up.h"
#include "TimeHistoryOutput_S2D_ME_s_FEM_up.h"
#include "TimeHistoryOutput_ConsoleProgressBar.h"

#include "test_sim_core.h"
#include "test_post_processor.h"

static double bgm_h = 0.25;
static size_t len = 1;

// implicit 1d bar compression
void test_fem_me_s_up_1dbar(void)
{
	Model_S2D_ME_s_FEM_up model;
	
	// background mesh
	model.init_mesh(bgm_h, 1, len);
	model.init_mat_param(1.0, 100.0, 0.0, 100.0/3.0);
	
	// velocity bc
	model.uy_num = 2;
	model.uys = new DisplacementBC[model.uy_num];
	model.uys[0].node_id = 0;
	model.uys[0].u = 0.0;
	model.uys[1].node_id = 1;
	model.uys[1].u = 0.0;
	model.ux_num = (len + 1) * 2;
	model.uxs = new DisplacementBC[model.ux_num];
	for (size_t n_id = 0; n_id < len + 1; ++n_id)
	{
		model.uxs[n_id].node_id = 2 * n_id;
		model.uxs[n_id].u = 0.0;
		model.uxs[len + 1 + n_id].node_id = 2 * n_id + 1;
		model.uxs[len + 1 + n_id].u = 0.0;
	}
	// traction bc
	model.ty_num = 1;
	model.tys = new TractionBC_2DFEM[model.ty_num];
	TractionBC_2DFEM &ty = model.tys[0];
	ty.elem_id = len - 1;
	ty.xi0 = -1.0;
	ty.xi1 = 1.0;
	ty.eta0 = 1.0;
	ty.eta1 = 1.0;
	ty.t0 = -1.0;
	ty.t1 = -1.0;

	//ResultFile_PlainBin res_file_pb;
	//res_file_pb.init("fem_me_up_res_1dbar.bin");
	ResultFile_XML res_file_xml;
	res_file_xml.init("fem_me_up_res_1dbar.xml");

	// output model
	ModelDataOutput_S2D_ME_s_FEM_up md;
	md.set_model(model);
	//md.set_res_file(res_file_pb);
	//md.output();
	md.set_res_file(res_file_xml);
	md.output();

	//TimeHistoryOutput_S2D_ME_s_FEM_up out1;
	//out1.set_res_file(res_file_pb);
	//out1.set_interval_num(50);
	//out1.set_output_init_state();
	TimeHistoryOutput_S2D_ME_s_FEM_up out2;
	out2.set_res_file(res_file_xml);
	out2.set_interval_num(50);
	out2.set_output_init_state();
	TimeHistoryOutput_ConsoleProgressBar cpb;
	
	Step_S2D_ME_s_FEM_up step;
	step.set_name("init_step");
	step.set_model(model);
	step.set_time(0.01);
	step.set_dtime(0.01);
	//step.add_time_history(out1);
	step.add_time_history(out2);
	//step.add_time_history(cpb);
	step.solve();

	system("pause");
}

//void test_animation_me_s_fem_up_1dbar(void)
//{
//	double height = double(len) * bgm_h;
//	double width = bgm_h;
//	double padding_height = height * 0.5;
//	double padding_width = width * 0.5;
//	GA_S2D_ME_s_FEM_up gen;
//	gen.generate(5.0, -padding_width, width + padding_width,
//		-padding_height, height + padding_height,
//		"mpm_me_up_res_1dbar.bin",
//		"mpm_me_up_res_1dbar.gif");
//}
