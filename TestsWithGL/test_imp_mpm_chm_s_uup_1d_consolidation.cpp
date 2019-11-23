#include "TestsWithGL_pcp.h"

#include "Model_S2D_CHM_s_uUp.h"
#include "Step_S2D_CHM_s_uUp.h"

#include "ModelDataOutput_S2D_CHM_s_uUp.h"
#include "TimeHistoryOutput_S2D_CHM_s_uUp.h"
#include "TimeHistoryOutput_ConsoleProgressBar.h"

#include "GA_S2D_CHM_s_uUp.h"

#include "test_sim_core.h"
#include "test_post_processor.h"

using namespace Model_S2D_CHM_s_uUp_Internal;

static size_t width = 1;
static size_t len = 5;
static size_t pcl_per_elem_len = 2;
static double bgm_h = 1.0 / double(len);

// implicit 1d consolidation
void test_imp_mpm_chm_s_uup_1d_consolidation(void)
{
	Model_S2D_CHM_s_uUp model;
	
	// background mesh
	model.init_mesh(bgm_h, width, len);

	// material points
	size_t pcl_num = (width * pcl_per_elem_len) * (len * pcl_per_elem_len);
	double pcl_area = bgm_h * bgm_h / double(pcl_per_elem_len * pcl_per_elem_len);
	model.init_pcl(pcl_num, 0.2, pcl_area*(1.0-0.2)*20.0, 20.0,
				   10.0, 1000.0, 0.2, 50000.0, 1.0e-4, 1.0);
	size_t pcl_id = 0;
	for (size_t y_id = 0; y_id < len * pcl_per_elem_len; ++y_id)
		for (size_t x_id = 0; x_id < width * pcl_per_elem_len; ++x_id)
		{
			Model_S2D_CHM_s_uUp::Particle &pcl = model.pcls[pcl_id++];
			pcl.x = (0.25 + double(x_id) * 0.5) * bgm_h;
			pcl.y = (0.25 + double(y_id) * 0.5) * bgm_h;
		}
	//model.pcls[3].x = 0.0   * N1(-0.5773502692, -0.5773502692) + bgm_h * N2(-0.5773502692, -0.5773502692)
	//				+ bgm_h * N3(-0.5773502692, -0.5773502692) + 0.0   * N4(-0.5773502692, -0.5773502692);
	//model.pcls[3].y = 0.0   * N1(-0.5773502692, -0.5773502692) + 0.0   * N2(-0.5773502692, -0.5773502692)
	//				+ bgm_h * N3(-0.5773502692, -0.5773502692) + bgm_h * N4(-0.5773502692, -0.5773502692);
	//model.pcls[2].x = 0.0   * N1( 0.5773502692, -0.5773502692) + bgm_h * N2( 0.5773502692, -0.5773502692)
	//				+ bgm_h * N3( 0.5773502692, -0.5773502692) + 0.0   * N4( 0.5773502692, -0.5773502692);
	//model.pcls[2].y = 0.0   * N1( 0.5773502692, -0.5773502692) + 0.0   * N2( 0.5773502692, -0.5773502692)
	//				+ bgm_h * N3( 0.5773502692, -0.5773502692) + bgm_h * N4( 0.5773502692, -0.5773502692);
	//model.pcls[1].x = 0.0   * N1( 0.5773502692,  0.5773502692) + bgm_h * N2( 0.5773502692,  0.5773502692)
	//				+ bgm_h * N3( 0.5773502692,  0.5773502692) + 0.0   * N4( 0.5773502692,  0.5773502692);
	//model.pcls[1].y = 0.0   * N1( 0.5773502692,  0.5773502692) + 0.0   * N2( 0.5773502692,  0.5773502692)
	//				+ bgm_h * N3( 0.5773502692,  0.5773502692) + bgm_h * N4( 0.5773502692,  0.5773502692);
	//model.pcls[0].x = 0.0   * N1(-0.5773502692,  0.5773502692) + bgm_h * N2(-0.5773502692,  0.5773502692)
	//				+ bgm_h * N3(-0.5773502692,  0.5773502692) + 0.0   * N4(-0.5773502692,  0.5773502692);
	//model.pcls[0].y = 0.0   * N1(-0.5773502692,  0.5773502692) + 0.0   * N2(-0.5773502692,  0.5773502692)
	//				+ bgm_h * N3(-0.5773502692,  0.5773502692) + bgm_h * N4(-0.5773502692,  0.5773502692);
	//model.pcls[3].vy_f = 0.1;
	//model.pcls[2].vy_f = 0.1;
	//model.pcls[1].vy_f = 0.1;
	//model.pcls[0].vy_f = 0.1;
	//model.pcls[3].s11 = 0.5;
	//model.pcls[2].s11 = 0.5;
	//model.pcls[1].s11 = 0.5;
	//model.pcls[0].s11 = 0.5;
	//model.pcls[3].s22 = 0.5;
	//model.pcls[2].s22 = 0.5;
	//model.pcls[1].s22 = 0.5;
	//model.pcls[0].s22 = 0.5;
	//model.pcls[3].s12 = 0.5;
	//model.pcls[2].s12 = 0.5;
	//model.pcls[1].s12 = 0.5;
	//model.pcls[0].s12 = 0.5;

	// displacement bc
	// solid phase
	model.usx_num = (len + 1) * 2;
	model.usxs = new DisplacementBC[model.usx_num];
	for (size_t n_id = 0; n_id < len + 1; ++n_id)
	{
		model.usxs[n_id].node_id = (width + 1) * n_id;
		model.usxs[n_id].u = 0.0;
		model.usxs[len + 1 + n_id].node_id = (width + 1) * (n_id + 1) - 1;
		model.usxs[len + 1 + n_id].u = 0.0;
	}
	model.usy_num = width + 1;
	model.usys = new DisplacementBC[model.usy_num];
	for (size_t n_id = 0; n_id < model.usy_num; ++n_id)
	{
		model.usys[n_id].node_id = n_id;
		model.usys[n_id].u = 0.0;
	}
	// fluid phase
	model.ufx_num = (len + 1) * 2;
	model.ufxs = new DisplacementBC[model.ufx_num];
	for (size_t n_id = 0; n_id < len + 1; ++n_id)
	{
		model.ufxs[n_id].node_id = (width + 1) * n_id;
		model.ufxs[n_id].u = 0.0;
		model.ufxs[len + 1 + n_id].node_id = (width + 1) * (n_id + 1) - 1;
		model.ufxs[len + 1 + n_id].u = 0.0;
	}
	model.ufy_num = width + 1;
	model.ufys = new DisplacementBC[model.ufy_num];
	for (size_t n_id = 0; n_id < model.ufy_num; ++n_id)
	{
		model.ufys[n_id].node_id = n_id;
		model.ufys[n_id].u = 0.0;
	}
	// pore pressure
	model.pbc_num = width + 1;
	model.pbcs = new PressureBC[model.pbc_num];
	for (size_t n_id = 0; n_id < model.pbc_num; ++n_id)
	{
		model.pbcs[n_id].node_id = (width + 1) * len + n_id;
		model.pbcs[n_id].p = 0.0;
	}

	// traction bc
	model.ty_num = width * pcl_per_elem_len;
	model.tys = new TractionBC_MPM[model.ty_num];
	for (size_t t_id = 0; t_id < model.ty_num; ++t_id)
	{
		//model.tys[t_id].pcl_id = t_id;
		model.tys[t_id].pcl_id = (len * pcl_per_elem_len - 1) * width * pcl_per_elem_len + t_id;
		model.tys[t_id].t = bgm_h / double(pcl_per_elem_len) * -400.0;
	}

	ResultFile_PlainBin res_file_pb;
	res_file_pb.init("mpm_chm_uup_res_1d_consoldation.bin");
	ResultFile_XML res_file_xml;
	res_file_xml.init("mpm_chm_uup_res_1d_consolidation.xml");

	// output model
	ModelDataOutput_S2D_CHM_s_uUp md;
	md.set_model(model);
	md.set_res_file(res_file_pb);
	md.output();
	md.set_res_file(res_file_xml);
	md.output();

	TimeHistoryOutput_S2D_CHM_s_uUp out1;
	out1.set_res_file(res_file_pb);
	out1.set_interval_num(50);
	out1.set_output_init_state();
	TimeHistoryOutput_S2D_CHM_s_uUp out2;
	out2.set_res_file(res_file_xml);
	out2.set_interval_num(50);
	out2.set_output_init_state();
	TimeHistoryOutput_ConsoleProgressBar cpb;
	
	Step_S2D_CHM_s_uUp step;
	step.set_name("init_step");
	step.set_model(model);
	step.set_time(15.0);
	step.set_dtime(0.01);
	step.add_time_history(out1);
	step.add_time_history(out2);
	//step.add_time_history(cpb);
	step.solve();

	system("pause");
}

void test_animation_chm_s_uup_1d_consolidation(void)
{
	double height = double(len) * bgm_h;
	double width = bgm_h;
	double padding_height = height * 0.5;
	double padding_width = width * 0.5;
	GA_S2D_CHM_s_uUp gen;
	gen.generate(5.0, -padding_width, width + padding_width,
		-padding_height, height + padding_height,
		"mpm_chm_uup_res_1d_consoldation.bin",
		"mpm_chm_uup_res_1d_consoldation.gif");
}
