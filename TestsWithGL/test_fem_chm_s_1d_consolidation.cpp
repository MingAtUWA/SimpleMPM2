#include "TestsWithGL_pcp.h"

#include "Model_S2D_CHM_s_FEM_uUp.h"
#include "Step_S2D_CHM_s_FEM_uUp.h"
#include "ModelDataOutput_S2D_CHM_s_FEM_uUp.h"
#include "TimeHistoryOutput_S2D_CHM_s_FEM_uUp.h"
#include "TimeHistoryOutput_ConsoleProgressBar.h"

#include "test_sim_core.h"

//#include "test_post_processor.h"
//#include "GA_S2D_CHM_s.h"

static size_t width = 1;
static size_t len = 10;
static double bgm_h = 1.0 / double(len);

void test_fem_chm_s_1d_consolidation(void)
{
	Model_S2D_CHM_s_FEM_uUp model;
	// background mesh
	model.init_mesh(bgm_h, width, len);
	model.init_mat_param(0.3, 2650.0, 1000.0, 1000.0, 0.25, 50000.0, 1.0e-4, 1.0);
	// traction bc
	model.ty_num = width;
	model.tys = new TractionBC_2DFEM[model.ty_num];
	for (size_t i = 0; i < model.ty_num; i++)
	{
		TractionBC_2DFEM &ty = model.tys[i];
		ty.elem_id = len - 1;
		ty.xi0 = -1.0;
		ty.xi1 = 1.0;
		ty.eta0 = 1.0;
		ty.eta1 = 1.0;
		ty.t0 = -1.0;
		ty.t1 = -1.0;
	}

	// velocity bc
	model.usy_num = width + 1;
	model.usys = new DisplacementBC[model.usy_num];
	for (size_t bc_id = 0; bc_id < model.usy_num; ++bc_id)
	{
		DisplacementBC &dbc = model.usys[bc_id];
		dbc.node_id = bc_id;
		dbc.u = 0.0;
	}
	model.usx_num = (len + 1) * 2;
	model.usxs = new DisplacementBC[model.usx_num];
	model.ufx_num = model.usx_num;
	model.ufxs = new DisplacementBC[model.ufx_num];
	for (size_t bc_id = 0; bc_id < len + 1; ++bc_id)
	{
		DisplacementBC &dbc1 = model.usxs[bc_id];
		dbc1.node_id = (width + 1) * bc_id;
		dbc1.u = 0.0;
		DisplacementBC &dbc2 = model.usxs[len + 1 + bc_id];
		dbc2.node_id = (width + 1) * (bc_id + 1) - 1;
		dbc2.u = 0.0;
		DisplacementBC &dbc3 = model.ufxs[bc_id];
		dbc3.node_id = (width + 1) * bc_id;
		dbc3.u = 0.0;
		DisplacementBC &dbc4 = model.ufxs[len + 1 + bc_id];
		dbc4.node_id = (width + 1) * (bc_id + 1) - 1;
		dbc4.u = 0.0;
	}
	model.ufy_num = (width + 1) * 2;
	model.ufys = new DisplacementBC[model.ufy_num];
	for (size_t bc_id = 0; bc_id < width + 1; ++bc_id)
	{
		DisplacementBC &dbc1 = model.ufys[bc_id];
		dbc1.node_id = bc_id;
		dbc1.u = 0.0;
		DisplacementBC &dbc2 = model.ufys[width + 1 + bc_id];
		dbc2.node_id = (width + 1) * len + bc_id;
		dbc2.u = 0.0;
	}

	//ResultFile_PlainBin res_file_pb;
	//res_file_pb.init("mpm_me_up_res_1dbar.bin");
	ResultFile_XML res_file_xml;
	res_file_xml.init("fem_chm_uup_res_1d_conso.xml");

	// output model
	ModelDataOutput_S2D_CHM_s_FEM_uUp md;
	md.set_model(model);
	//md.set_res_file(res_file_pb);
	//md.output();
	md.set_res_file(res_file_xml);
	md.output();

	//TimeHistoryOutput_S2D_CHM_s_FEM_uUp out1;
	//out1.set_res_file(res_file_pb);
	//out1.set_interval_num(50);
	//out1.set_output_init_state();
	TimeHistoryOutput_S2D_CHM_s_FEM_uUp out2;
	out2.set_res_file(res_file_xml);
	out2.set_interval_num(150);
	out2.set_output_init_state();
	//TimeHistoryOutput_ConsoleProgressBar cpb;

	Step_S2D_CHM_s_FEM_uUp step;
	step.set_name("consolidation");
	step.set_model(model);
	step.set_time(15.0);
	step.set_dtime(0.05);
	//step.add_time_history(out1);
	step.add_time_history(out2);
	//step.add_time_history(cpb);

	// freely drainage bc
	model.ufy_num = width + 1;
	model.pbc_num = width + 1;
	model.pbcs = new PressureBC[model.pbc_num];
	for (size_t bc_id = 0; bc_id < model.pbc_num; ++bc_id)
	{
		PressureBC &pbc = model.pbcs[bc_id];
		pbc.node_id = (width + 1) * len + bc_id;
		pbc.p = 0.0;
	}

	step.solve();

	system("pause");
}
