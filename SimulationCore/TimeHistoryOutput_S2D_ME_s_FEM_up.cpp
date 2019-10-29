#include "SimulationCore_pcp.h"

#include "Model_S2D_ME_s_FEM_up.h"
#include "Step_S2D_ME_s_FEM_up.h"

#include "ResultFile_PlainBin_DataStruct.h"

#include "TimeHistoryOutput_S2D_ME_s_FEM_up.h"

int TimeHistoryOutput_S2D_ME_s_FEM_up::init_per_step(void)
{
	if (is_plain_bin)
	{
		Model_S2D_ME_s_FEM_up &md = *static_cast<Model_S2D_ME_s_FEM_up *>(model);
		node_data = new double[md.node_num];
		gp_data = new double[md.elem_num];
	}
	return 0;
}

void TimeHistoryOutput_S2D_ME_s_FEM_up::finalize_per_step(void)
{
	if (node_data)
	{
		delete[] node_data;
		node_data = nullptr;
	}
	if (gp_data)
	{
		delete[] gp_data;
		gp_data = nullptr;
	}
}

int time_history_output_func_s2d_me_s_fem_up_to_plain_bin_res_file(TimeHistoryOutput &_self)
{
	TimeHistoryOutput_S2D_ME_s_FEM_up &th = static_cast<TimeHistoryOutput_S2D_ME_s_FEM_up &>(_self);
	Model_S2D_ME_s_FEM_up &model = static_cast<Model_S2D_ME_s_FEM_up &>(th.get_model());
	Step_S2D_ME_s_FEM_up &step = static_cast<Step_S2D_ME_s_FEM_up &>(th.get_step());
	ResultFile_PlainBin &rf = static_cast<ResultFile_PlainBin &>(*th.res_file);
	std::fstream &file = rf.get_file();

	typedef ResultFile_PlainBin_DataStruct::TimeHistoryHeader TimeHistoryHeader;
	typedef ResultFile_PlainBin_DataStruct::MeshObjectHeader_2D4R MeshObjectHeader;
	
	// time history header
	TimeHistoryHeader thh;
	thh.init();
	thh.substep_num = step.get_substep_num();
	thh.total_substep_num = step.get_total_substep_num();
	thh.current_time = step.get_current_time();
	thh.total_time = step.get_total_time();
	file.write(reinterpret_cast<char *>(&thh), sizeof(thh));

	// output node data
	MeshObjectHeader mh;
	mh.node_num = model.node_num;
	mh.node_fld_num = 2;
	mh.elem_num = model.elem_num;
	mh.elem_fld_num = 6;
	file.write(reinterpret_cast<char *>(&mh), sizeof(mh));
	// node data
	// ux
	for (size_t n_id = 0; n_id < model.node_num; ++n_id)
		th.node_data[n_id] = model.nodes[n_id].ux;
	file.write(reinterpret_cast<char *>(&th.node_data), sizeof(double) * model.node_num);
	// uy
	for (size_t n_id = 0; n_id < model.node_num; ++n_id)
		th.node_data[n_id] = model.nodes[n_id].uy;
	file.write(reinterpret_cast<char *>(&th.node_data), sizeof(double) * model.node_num);
	// gauss point data
	//// ux
	//for (size_t gp_id = 0; gp_id < model.elem_num; ++gp_id)
	//	th.gp_data[gp_id] = model.elems[gp_id].gp.ux;
	//file.write(reinterpret_cast<char *>(&th.gp_data), sizeof(double) * model.elem_num);
	//// uy
	//for (size_t gp_id = 0; gp_id < model.elem_num; ++gp_id)
	//	th.gp_data[gp_id] = model.elems[gp_id].gp.uy;
	//file.write(reinterpret_cast<char *>(&th.gp_data), sizeof(double) * model.elem_num);
	//// p
	//for (size_t gp_id = 0; gp_id < model.elem_num; ++gp_id)
	//	th.gp_data[gp_id] = model.elems[gp_id].gp.p;
	//file.write(reinterpret_cast<char *>(&th.gp_data), sizeof(double) * model.elem_num);
	//// s11
	//for (size_t gp_id = 0; gp_id < model.elem_num; ++gp_id)
	//	th.gp_data[gp_id] = model.elems[gp_id].gp.s11;
	//file.write(reinterpret_cast<char *>(&th.gp_data), sizeof(double) * model.elem_num);
	//// s22
	//for (size_t gp_id = 0; gp_id < model.elem_num; ++gp_id)
	//	th.gp_data[gp_id] = model.elems[gp_id].gp.s22;
	//file.write(reinterpret_cast<char *>(&th.gp_data), sizeof(double) * model.elem_num);
	//// s12
	//for (size_t gp_id = 0; gp_id < model.elem_num; ++gp_id)
	//	th.gp_data[gp_id] = model.elems[gp_id].gp.s12;
	//file.write(reinterpret_cast<char *>(&th.gp_data), sizeof(double) * model.elem_num);

	return 0;
}

int time_history_output_func_s2d_me_s_fem_up_to_xml_res_file(TimeHistoryOutput &_self)
{
	TimeHistoryOutput_S2D_ME_s_FEM_up &th = static_cast<TimeHistoryOutput_S2D_ME_s_FEM_up &>(_self);
	ResultFile_XML &rf = static_cast<ResultFile_XML &>(*th.res_file);
	std::fstream &file = rf.get_file();
	
	char str_buffer[512];
#define str_buffer_len (sizeof(str_buffer) / sizeof(str_buffer[0]))

	// time history
	Step_S2D_ME_s_FEM_up &step = static_cast<Step_S2D_ME_s_FEM_up &>(th.get_step());
	const char *time_history_info = ""
		"<TimeHistory>\n"
		"    <substep_num> %zu </substep_num>\n"
		"    <total_substep_num> %zu </total_substep_num>\n"
		"    <current_time> %16.10e </current_time>\n"
		"    <total_time> %16.10e </total_time>\n";
	snprintf(str_buffer, str_buffer_len, time_history_info,
		step.get_substep_num(), step.get_total_substep_num(),
		step.get_current_time(), step.get_total_time());
	file.write(str_buffer, strlen(str_buffer));

	// output mesh data (nodes and gauss points)
	Model_S2D_ME_s_FEM_up &model = static_cast<Model_S2D_ME_s_FEM_up &>(th.get_model());
	file << "    <MeshObject>\n";
	// node data
	const char *node_info = "        <node_num> %zu </node_num>\n";
	snprintf(str_buffer, str_buffer_len, node_info, model.node_num);
	file.write(str_buffer, strlen(str_buffer));
	file << "        <node_data>\n"
			"        <!-- ux, uy, p -->\n";
	for (size_t n_id = 0; n_id < model.node_num; ++n_id)
	{
		Model_S2D_ME_s_FEM_up::Node &n = model.nodes[n_id];
		file << "        ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", n.ux);
		file << str_buffer << ", ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", n.uy);
		file << str_buffer << ", ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", n.p);
		file << str_buffer << "\n";
	}
	file << "        </node_data>\n";
	// gauss point data
	const char *gauss_point_info = "        <gauss_point_num> %zu </gauss_point_num>\n";
	snprintf(str_buffer, str_buffer_len, gauss_point_info, model.elem_num * 4);
	file.write(str_buffer, strlen(str_buffer));
	file << "        <gauss_point_data>\n"
			"        <!-- ux, uy, p, s11, s22, s12 -->\n";
	for (size_t e_id = 0; e_id < model.elem_num; ++e_id)
	{
		Model_S2D_ME_s_FEM_up::Element &e = model.elems[e_id];
		for (size_t gp_id = 0; gp_id < 4; ++gp_id)
		{
			Model_S2D_ME_s_FEM_up::GaussPoint &gp = e.gps[gp_id];
			file << "        ";
			snprintf(str_buffer, str_buffer_len, "%16.10e", gp.ux);
			file << str_buffer << ", ";
			snprintf(str_buffer, str_buffer_len, "%16.10e", gp.uy);
			file << str_buffer << ", ";
			snprintf(str_buffer, str_buffer_len, "%16.10e", gp.p);
			file << str_buffer << ", ";
			snprintf(str_buffer, str_buffer_len, "%16.10e", gp.s11);
			file << str_buffer << ", ";
			snprintf(str_buffer, str_buffer_len, "%16.10e", gp.s22);
			file << str_buffer << ", ";
			snprintf(str_buffer, str_buffer_len, "%16.10e", gp.s12);
			file << str_buffer << "\n";
		}
	}
	file << "        </gauss_point_data>\n";
	// ending
	file << "    </MeshObject>\n";

	// ending
	file << "</TimeHistory>\n";

	return 0;
}
