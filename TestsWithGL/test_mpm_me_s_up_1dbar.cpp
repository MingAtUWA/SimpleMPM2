#include "TestsWithGL_pcp.h"

#include "Model_S2D_ME_s_up.h"
#include "Step_S2D_ME_s_up.h"

#include "ModelDataOutput_S2D_ME_s_up.h"
#include "TimeHistoryOutput_S2D_ME_s_up.h"
#include "TimeHistoryOutput_ConsoleProgressBar.h"

#include "GA_S2D_ME_s_up.h"

#include "test_sim_core.h"
#include "test_post_processor.h"

using namespace Model_S2D_ME_s_up_Internal;

static double bgm_h = 0.25;
static size_t width = 1;
static size_t len = 2;

// implicit 1d bar compression
void test_mpm_me_s_up_1dbar(void)
{
	Model_S2D_ME_s_up model;
	
	// background mesh
	model.init_mesh(bgm_h, width, len);

	// material points
	size_t pcl_per_elem_len = 2;
	size_t pcl_num = (width * pcl_per_elem_len) * (len * pcl_per_elem_len);
	double pcl_area = bgm_h * bgm_h / double(pcl_per_elem_len * pcl_per_elem_len);
	model.init_pcl(pcl_num, pcl_area, 1.0, 100.0, 0.0, 100.0 / 3.0);
	size_t pcl_id = 0;
	for (size_t y_id = 0; y_id < len * pcl_per_elem_len; ++y_id)
		for (size_t x_id = 0; x_id < width * pcl_per_elem_len; ++x_id)
		{
			Model_S2D_ME_s_up::Particle &pcl = model.pcls[pcl_id++];
			pcl.x = (0.25 + double(x_id) * 0.5) * bgm_h;
			pcl.y = (0.25 + double(y_id) * 0.5) * bgm_h;
		}
	//model.pcls[0].x = 0.0   * N1(-0.5773502692, -0.5773502692) + bgm_h * N2(-0.5773502692, -0.5773502692)
	//				+ bgm_h * N3(-0.5773502692, -0.5773502692) + 0.0   * N4(-0.5773502692, -0.5773502692);
	//model.pcls[0].y = 0.0   * N1(-0.5773502692, -0.5773502692) + 0.0   * N2(-0.5773502692, -0.5773502692)
	//				+ bgm_h * N3(-0.5773502692, -0.5773502692) + bgm_h * N4(-0.5773502692, -0.5773502692);
	//model.pcls[1].x = 0.0   * N1( 0.5773502692, -0.5773502692) + bgm_h * N2( 0.5773502692, -0.5773502692)
	//				+ bgm_h * N3( 0.5773502692, -0.5773502692) + 0.0   * N4( 0.5773502692, -0.5773502692);
	//model.pcls[1].y = 0.0   * N1( 0.5773502692, -0.5773502692) + 0.0   * N2( 0.5773502692, -0.5773502692)
	//				+ bgm_h * N3( 0.5773502692, -0.5773502692) + bgm_h * N4( 0.5773502692, -0.5773502692);
	//model.pcls[2].x = 0.0   * N1( 0.5773502692,  0.5773502692) + bgm_h * N2( 0.5773502692,  0.5773502692)
	//				+ bgm_h * N3( 0.5773502692,  0.5773502692) + 0.0   * N4( 0.5773502692,  0.5773502692);
	//model.pcls[2].y = 0.0   * N1( 0.5773502692,  0.5773502692) + 0.0   * N2( 0.5773502692,  0.5773502692)
	//				+ bgm_h * N3( 0.5773502692,  0.5773502692) + bgm_h * N4( 0.5773502692,  0.5773502692);
	//model.pcls[3].x = 0.0   * N1(-0.5773502692,  0.5773502692) + bgm_h * N2(-0.5773502692,  0.5773502692)
	//				+ bgm_h * N3(-0.5773502692,  0.5773502692) + 0.0   * N4(-0.5773502692,  0.5773502692);
	//model.pcls[3].y = 0.0   * N1(-0.5773502692,  0.5773502692) + 0.0   * N2(-0.5773502692,  0.5773502692)
	//				+ bgm_h * N3(-0.5773502692,  0.5773502692) + bgm_h * N4(-0.5773502692,  0.5773502692);
	
	// velocity bc
	model.uy_num = width + 1;
	model.uys = new DisplacementBC[model.uy_num];
	for (size_t n_id = 0; n_id < model.uy_num; ++n_id)
	{
		model.uys[n_id].node_id = n_id;
		model.uys[n_id].u = 0.0;
	}
	model.ux_num = (len + 1) * 2;
	model.uxs = new DisplacementBC[model.ux_num];
	for (size_t n_id = 0; n_id < len + 1; ++n_id)
	{
		model.uxs[n_id].node_id = (width + 1) * n_id;
		model.uxs[n_id].u = 0.0;
		model.uxs[len + 1 + n_id].node_id = (width + 1) * (n_id + 1) - 1;
		model.uxs[len + 1 + n_id].u = 0.0;
	}
	// traction bc
	model.ty_num = width * pcl_per_elem_len;
	model.tys = new TractionBC_MPM[model.ty_num];
	for (size_t t_id = 0; t_id < model.ty_num; ++t_id)
	{
		model.tys[t_id].pcl_id = (len * pcl_per_elem_len - 1) * width * pcl_per_elem_len + t_id;
		model.tys[t_id].t = bgm_h / double(pcl_per_elem_len) * -10.0;
	}

	ResultFile_PlainBin res_file_pb;
	res_file_pb.init("mpm_me_up_res_1dbar.bin");
	ResultFile_XML res_file_xml;
	res_file_xml.init("mpm_me_up_res_1dbar.xml");

	// output model
	ModelDataOutput_S2D_ME_s_up md;
	md.set_model(model);
	md.set_res_file(res_file_pb);
	md.output();
	md.set_res_file(res_file_xml);
	md.output();

	TimeHistoryOutput_S2D_ME_s_up out1;
	out1.set_res_file(res_file_pb);
	out1.set_interval_num(78);
	out1.set_output_init_state();
	TimeHistoryOutput_S2D_ME_s_up out2;
	out2.set_res_file(res_file_xml);
	out2.set_interval_num(78);
	out2.set_output_init_state();
	TimeHistoryOutput_ConsoleProgressBar cpb;
	
	Step_S2D_ME_s_up step;
	step.set_name("init_step");
	step.set_model(model);
	step.set_time(0.78);
	step.set_dtime(0.01);
	step.add_time_history(out1);
	step.add_time_history(out2);
	//step.add_time_history(cpb);
	step.solve();

	//system("pause");
}

void test_animation_me_s_up_1dbar(void)
{
	double height = double(len) * bgm_h;
	double width = bgm_h;
	double padding_height = height * 0.5;
	double padding_width = width * 0.5;
	GA_S2D_ME_s_up gen;
	gen.generate(5.0, -padding_width, width + padding_width,
		-padding_height, height + padding_height,
		"mpm_me_up_res_1dbar.bin",
		"mpm_me_up_res_1dbar.gif");
}
