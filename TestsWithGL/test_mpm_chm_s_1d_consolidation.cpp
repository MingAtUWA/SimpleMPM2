#include "TestsWithGL_pcp.h"

#include "test_sim_core.h"

#include "Model_S2D_CHM_s.h"
#include "Step_S2D_CHM_s.h"
#include "ModelDataOutput_S2D_CHM_s.h"
#include "TimeHistoryOutput_S2D_CHM_s.h"
#include "TimeHistoryOutput_ConsoleProgressBar.h"

void test_mpm_chm_s_1d_consolidation(void)
{
	Model_S2D_CHM_s model;
	size_t elem_x_num = 2;
	size_t elem_y_num = 10;
	double elem_len = 1.0 / (elem_y_num - 1);
	model.init_mesh(elem_len, elem_x_num, elem_y_num);

	size_t pcl_per_elem_len = 2;
	double m_s = elem_len * elem_len / (pcl_per_elem_len * pcl_per_elem_len) * 2650.0;
	model.init_pcl(elem_x_num * elem_y_num * pcl_per_elem_len * pcl_per_elem_len,
				   0.3, m_s, 2650.0, 1000.0, 1.0e3, 0.25, 50.0e3, 1.0e-4, 1.0);
	
	model.ty_num = model.elem_x_num * 2;
	model.tys = new TractionBC_MPM[model.ty_num];
	for (size_t i = 0; i < model.ty_num; i++)
	{
		model.tys[i].pcl_id = (model.elem_y_num * 2 - 1) * model.elem_x_num * 2 + i;
		model.tys[i].t = -10.0 * 0.5 * elem_len;
	}
	
	model.vsx_num = model.node_y_num * 2;
	model.vsxs = new VelocityBC[model.vsx_num];
	for (size_t i = 0; i < model.node_y_num; ++i)
	{
		model.vsxs[i].node_id = i * model.node_x_num;
		model.vsxs[i].v = 0.0;
		model.vsxs[i + model.node_y_num].node_id = (i + 1) * model.node_x_num - 1;
		model.vsxs[i + model.node_y_num].v = 0.0;
	}
	model.vsy_num = model.node_x_num;
	model.vsys = new VelocityBC[model.vsy_num];
	for (size_t i = 0; i < model.node_x_num; ++i)
	{
		model.vsys[i].node_id = i;
		model.vsys[i].v = 0.0;
	}

	model.vfx_num = model.node_y_num * 2;
	model.vfxs = new VelocityBC[model.vfx_num];
	for (size_t i = 0; i < model.node_y_num; ++i)
	{
		model.vfxs[i].node_id = i * model.node_x_num;
		model.vfxs[i].v = 0.0;
		model.vfxs[i + model.node_y_num].node_id = (i + 1) * model.node_x_num - 1;
		model.vfxs[i + model.node_y_num].v = 0.0;
	}
	model.vfy_num = model.node_x_num * 2;
	model.vfys = new VelocityBC[model.vfy_num];
	for (size_t i = 0; i < model.node_x_num; ++i)
	{
		model.vfys[i].node_id = i;
		model.vfys[i].v = 0.0;
		model.vfys[i + model.node_x_num].node_id = (model.node_y_num - 1) * model.node_x_num + i;
		model.vfys[i + model.node_x_num].v = 0.0;
	}

	ResultFile_PlainBin res_file_pb;
	res_file_pb.init("mpm_rb_res_bar_compression.bin");
	ResultFile_XML res_file_xml;
	res_file_xml.init("mpm_rb_res_bar_compression.xml");

	// output model
	ModelDataOutput_S2D_CHM_s md;
	md.set_model(model);
	md.set_res_file(res_file_pb);
	md.output();
	md.set_res_file(res_file_xml);
	md.output();

	TimeHistoryOutput_S2D_CHM_s out1;
	out1.set_res_file(res_file_pb);
	out1.set_interval_num(10);
	out1.set_output_init_state();
	TimeHistoryOutput_S2D_CHM_s out2;
	out2.set_res_file(res_file_xml);
	out2.set_interval_num(50);
	out2.set_output_init_state();
	TimeHistoryOutput_ConsoleProgressBar cpb;

	Step_S2D_CHM_s step1;
	step1.set_name("geostatic");
	step1.set_model(model);
	step1.set_time(5.0);
	step1.set_dtime(1.0e-5);
	step1.add_time_history(out1);
	step1.add_time_history(out2);
	step1.add_time_history(cpb);
	step1.solve();

	Step_S2D_CHM_s step2;
	step2.set_name("consolidation");
	step2.set_model(model);
	step2.set_time(30.0);
	step2.set_dtime(1.0e-5);
	step2.add_time_history(out1);
	step2.add_time_history(out2);
	step2.add_time_history(cpb);
	// free drainage bcs
	model.vfy_num = model.node_x_num;
	step2.solve();
}
