#include "SimulationCore_pcp.h"

#include "Model_S2D_CHM_s_FEM_uUp.h"
#include "Step_S2D_CHM_s_FEM_uUp.h"

#include "ResultFile_PlainBin_DataStruct.h"

#include "TimeHistoryOutput_S2D_CHM_s_FEM_uUp.h"

int time_history_output_func_s2d_chm_s_fem_uup_to_plain_bin_res_file(TimeHistoryOutput &_self)
{
	//TimeHistoryOutput_S2D_CHM_s_FEM_uUp &th = static_cast<TimeHistoryOutput_S2D_CHM_s_FEM_uUp &>(_self);
	//Model_S2D_CHM_s_FEM_uUp &model = static_cast<Model_S2D_CHM_s_FEM_uUp &>(th.get_model());
	//Step_S2D_CHM_s_FEM_uUp &step = static_cast<Step_S2D_CHM_s_FEM_uUp &>(th.get_step());
	//ResultFile_PlainBin &rf = static_cast<ResultFile_PlainBin &>(*th.res_file);
	//std::fstream &file = rf.get_file();

	//typedef ResultFile_PlainBin_DataStruct::TimeHistoryHeader TimeHistoryHeader;
	//typedef ResultFile_PlainBin_DataStruct::MPObjectHeader MPObjectHeader;

	//// time history header
	//TimeHistoryHeader thh;
	//thh.init();
	//thh.substep_num = step.get_substep_num();
	//thh.total_substep_num = step.get_total_substep_num();
	//thh.current_time = step.get_current_time();
	//thh.total_time = step.get_total_time();
	//file.write(reinterpret_cast<char *>(&thh), sizeof(thh));

	//// output particles data
	//MPObjectHeader mph;
	//mph.init();
	//mph.pcl_num = model.pcl_num;
	//mph.fld_num = 4;
	//file.write(reinterpret_cast<char *>(&mph), sizeof(mph));
	//double *pcl_data = new double[mph.pcl_num];
	//size_t data_len = sizeof(double) * mph.pcl_num;
	//// x
	//for (size_t pcl_id = 0; pcl_id < mph.pcl_num; ++pcl_id)
	//	pcl_data[pcl_id] = model.pcls[pcl_id].x;
	//file.write(reinterpret_cast<char *>(pcl_data), data_len);
	//// y
	//for (size_t pcl_id = 0; pcl_id < mph.pcl_num; ++pcl_id)
	//	pcl_data[pcl_id] = model.pcls[pcl_id].y;
	//file.write(reinterpret_cast<char *>(pcl_data), data_len);
	//// vol
	//for (size_t pcl_id = 0; pcl_id < mph.pcl_num; ++pcl_id)
	//{
	//	Model_S2D_CHM_s_FEM_uUp::Particle &pcl = model.pcls[pcl_id];
	//	pcl_data[pcl_id] = pcl.m_s / ((1.0 - pcl.n) * pcl.density_s);
	//}
	//file.write(reinterpret_cast<char *>(pcl_data), data_len);
	//// p
	//for (size_t pcl_id = 0; pcl_id < mph.pcl_num; ++pcl_id)
	//{
	//	Model_S2D_CHM_s_FEM_uUp::Particle &pcl = model.pcls[pcl_id];
	//	pcl_data[pcl_id] = pcl.p;
	//}
	//file.write(reinterpret_cast<char *>(pcl_data), data_len);
	//delete[] pcl_data;

	return 0;
}

int time_history_output_func_s2d_chm_s_fem_uup_to_xml_res_file(TimeHistoryOutput &_self)
{
	TimeHistoryOutput_S2D_CHM_s_FEM_uUp &th = static_cast<TimeHistoryOutput_S2D_CHM_s_FEM_uUp &>(_self);
	ResultFile_XML &rf = static_cast<ResultFile_XML &>(*th.res_file);
	std::fstream &file = rf.get_file();

	char str_buffer[512];
#define str_buffer_len (sizeof(str_buffer) / sizeof(str_buffer[0]))

	// time history
	Step_S2D_CHM_s_FEM_uUp &step = static_cast<Step_S2D_CHM_s_FEM_uUp &>(th.get_step());
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
	Model_S2D_CHM_s_FEM_uUp &model = static_cast<Model_S2D_CHM_s_FEM_uUp &>(th.get_model());
	file << "    <MeshObject>\n";
	// node data
	const char *node_info = "        <node_num> %zu </node_num>\n";
	snprintf(str_buffer, str_buffer_len, node_info, model.node_num);
	file.write(str_buffer, strlen(str_buffer));
	file << "        <node_data>\n"
			"        <!-- ux_s, uy_s, ux_f, uy_f, p -->\n";
	for (size_t n_id = 0; n_id < model.node_num; ++n_id)
	{
		Model_S2D_CHM_s_FEM_uUp::Node &n = model.nodes[n_id];
		file << "        ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", n.ux_s);
		file << str_buffer << ", ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", n.uy_s);
		file << str_buffer << ", ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", n.ux_f);
		file << str_buffer << ", ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", n.uy_f);
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
			"        <!-- ux_s, uy_s, ux_f, uy_f, p -->\n";
	for (size_t e_id = 0; e_id < model.elem_num; ++e_id)
	{
		Model_S2D_CHM_s_FEM_uUp::Element &e = model.elems[e_id];
		for (size_t gp_id = 0; gp_id < 4; ++gp_id)
		{
			Model_S2D_CHM_s_FEM_uUp::GaussPoint &gp = e.gps[gp_id];
			file << "        ";
			snprintf(str_buffer, str_buffer_len, "%16.10e", gp.ux_s);
			file << str_buffer << ", ";
			snprintf(str_buffer, str_buffer_len, "%16.10e", gp.uy_s);
			file << str_buffer << ", ";
			snprintf(str_buffer, str_buffer_len, "%16.10e", gp.ux_f);
			file << str_buffer << ", ";
			snprintf(str_buffer, str_buffer_len, "%16.10e", gp.uy_f);
			file << str_buffer << ", ";
			snprintf(str_buffer, str_buffer_len, "%16.10e", gp.p);
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
