#include "SimulationCore_pcp.h"

#include "Model_T2D_CHM_s.h"
#include "Step_T2D_CHM_s.h"

#include "ResultFile_PlainBin_DataStruct.h"

#include "TimeHistoryOutput_T2D_CHM_s.h"

int time_history_output_func_t2d_chm_s_to_plain_bin_res_file(TimeHistoryOutput &_self)
{
	TimeHistoryOutput_T2D_CHM_s &th = static_cast<TimeHistoryOutput_T2D_CHM_s &>(_self);
	Model_T2D_CHM_s &model = static_cast<Model_T2D_CHM_s &>(th.get_model());
	Step_T2D_CHM_s &step = static_cast<Step_T2D_CHM_s &>(th.get_step());
	ResultFile_PlainBin &rf = static_cast<ResultFile_PlainBin &>(*th.res_file);
	std::fstream &file = rf.get_file();

	typedef ResultFile_PlainBin_DataStruct::TimeHistoryHeader TimeHistoryHeader;
	typedef ResultFile_PlainBin_DataStruct::MPObjectHeader MPObjectHeader;
	
	// time history header
	TimeHistoryHeader thh;
	thh.init();
	thh.substep_num = step.get_substep_num();
	thh.total_substep_num = step.get_total_substep_num();
	thh.current_time = step.get_current_time();
	thh.total_time = step.get_total_time();
	file.write(reinterpret_cast<char *>(&thh), sizeof(thh));

	// output particles data
	MPObjectHeader mph;
	mph.init();
	mph.pcl_num = model.pcl_num;
	mph.fld_num = 4;
	file.write(reinterpret_cast<char *>(&mph), sizeof(mph));
	double *pcl_data = new double[mph.pcl_num];
	size_t data_len = sizeof(double) * mph.pcl_num;
	// x
	for (size_t pcl_id = 0; pcl_id < mph.pcl_num; ++pcl_id)
		pcl_data[pcl_id] = model.pcls[pcl_id].x;
	file.write(reinterpret_cast<char *>(pcl_data), data_len);
	// y
	for (size_t pcl_id = 0; pcl_id < mph.pcl_num; ++pcl_id)
		pcl_data[pcl_id] = model.pcls[pcl_id].y;
	file.write(reinterpret_cast<char *>(pcl_data), data_len);
	// vol
	for (size_t pcl_id = 0; pcl_id < mph.pcl_num; ++pcl_id)
	{
		Model_T2D_CHM_s::Particle &pcl = model.pcls[pcl_id];
		pcl_data[pcl_id] = pcl.m_s / ((1.0 - pcl.n) * pcl.density_s);
	}
	file.write(reinterpret_cast<char *>(pcl_data), data_len);
	// p
	for (size_t pcl_id = 0; pcl_id < mph.pcl_num; ++pcl_id)
		pcl_data[pcl_id] = model.pcls[pcl_id].p;
	file.write(reinterpret_cast<char *>(pcl_data), data_len);
	delete[] pcl_data;

	return 0;
}

int time_history_output_func_t2d_chm_s_to_xml_res_file(TimeHistoryOutput &_self)
{
	TimeHistoryOutput_T2D_CHM_s &th = static_cast<TimeHistoryOutput_T2D_CHM_s &>(_self);
	ResultFile_XML &rf = static_cast<ResultFile_XML &>(*th.res_file);
	std::fstream &file = rf.get_file();
	
	char str_buffer[512];
#define str_buffer_len (sizeof(str_buffer) / sizeof(str_buffer[0]))

	// time history
	Step_T2D_CHM_s &step = static_cast<Step_T2D_CHM_s &>(th.get_step());
	const char *time_history_info = ""
		"<TimeHistory>\n"
		"    <substep_num> %zu </substep_num>\n"
		"    <total_substep_num> %zu </total_substep_num>\n"
		"    <current_time> %16.10e </current_time>\n"
		"    <total_time> %16.10e </total_time>\n";
	snprintf(str_buffer, str_buffer_len, time_history_info,
		step.get_substep_num(),  step.get_total_substep_num(),
		step.get_current_time(), step.get_total_time());
	file.write(str_buffer, strlen(str_buffer));

	// output material points data
	Model_T2D_CHM_s &model = static_cast<Model_T2D_CHM_s &>(th.get_model());
	const char *material_point_info = ""
		"    <MaterialPointObject>\n"
		"        <pcl_num> %zu </pcl_num>\n";
	snprintf(str_buffer, str_buffer_len, material_point_info, model.pcl_num);
	file.write(str_buffer, strlen(str_buffer));
	// field data: x, y, vol
	file << "        <field_data>\n"
			"        <!-- x, y, vol, p, n, s11, s22, s12 -->\n";
	for (size_t pcl_id = 0; pcl_id < model.pcl_num; ++pcl_id)
	{
		Model_T2D_CHM_s::Particle &pcl = model.pcls[pcl_id];
		file << "        ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", pcl.x);
		file << str_buffer << ", ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", pcl.y);
		file << str_buffer << ", ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", pcl.m_s / ((1.0 - pcl.n) * pcl.density_s));
		file << str_buffer << ", ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", pcl.p);
		file << str_buffer << ", ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", pcl.n);
		file << str_buffer << ", ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", pcl.s11);
		file << str_buffer << ", ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", pcl.s22);
		file << str_buffer << ", ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", pcl.s12);
		file << str_buffer << "\n";
	}
	file << "        </field_data>\n"
			"    </MaterialPointObject>\n";

	// ending
	file << "</TimeHistory>\n";

	return 0;
}
