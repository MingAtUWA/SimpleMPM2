#include "SimulationCore_pcp.h"

#include "Model_T2D_CHM_s.h"
#include "Step_T2D_CHM_s_SE_Geostatic.h"

#include "ResultFile_PlainBin_DataStruct.h"
#include "Model_T2D_CHM_s_hdf5_io_utilities.h"

#include "TimeHistoryOutput_T2D_CHM_s_SE_Geostatic.h"

int TimeHistoryOutput_T2D_CHM_s_SE_Geostatic::init(void)
{
	switch (res_file->get_type())
	{
	case ResultFileType::PlainBin:
		break;
	case ResultFileType::XML:
		break;
	case ResultFileType::Hdf5:
	{
		ResultFile_hdf5 &rf = *static_cast<ResultFile_hdf5 *>(res_file);
		hid_t th_grp_id = rf.get_time_history_grp_id();
		th_id = rf.create_group(th_grp_id, name.c_str());
	}
	break;
	default:
		break;
	}
	return 0;
}

int TimeHistoryOutput_T2D_CHM_s_SE_Geostatic::close(void)
{
	switch (res_file->get_type())
	{
	case ResultFileType::PlainBin:
		break;
	case ResultFileType::XML:
		break;
	case ResultFileType::Hdf5:
		{
		ResultFile_hdf5 &rf = *static_cast<ResultFile_hdf5 *>(res_file);
		rf.write_attribute(th_id, "output_num", output_id);
		if (th_id > 0)
		{
			rf.close_group(th_id);
			th_id = -1;
		}
		}
		break;
	default:
		break;
	}
	return 0;
}

int time_history_output_func_t2d_chm_s_SE_geostatic_to_plain_bin_res_file(TimeHistoryOutput &_self)
{
	TimeHistoryOutput_T2D_CHM_s_SE_Geostatic &th
		= static_cast<TimeHistoryOutput_T2D_CHM_s_SE_Geostatic &>(_self);
	Model_T2D_CHM_s &model = static_cast<Model_T2D_CHM_s &>(th.get_model());
	Step_T2D_CHM_s_SE_Geostatic &step
		= static_cast<Step_T2D_CHM_s_SE_Geostatic &>(th.get_step());
	ResultFile_PlainBin &rf = static_cast<ResultFile_PlainBin &>(*th.res_file);
	std::fstream &file = rf.get_file();
	double *pcl_data;
	size_t data_len;

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
	pcl_data = new double[mph.pcl_num];
	data_len = sizeof(double) * mph.pcl_num;
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
	// s22
	for (size_t pcl_id = 0; pcl_id < mph.pcl_num; ++pcl_id)
		pcl_data[pcl_id] = model.pcls[pcl_id].s22;
	file.write(reinterpret_cast<char *>(pcl_data), data_len);
	delete[] pcl_data;

	return 0;
}

int time_history_output_func_t2d_chm_s_SE_geostatic_to_xml_res_file(TimeHistoryOutput &_self)
{
	TimeHistoryOutput_T2D_CHM_s_SE_Geostatic &th
		= static_cast<TimeHistoryOutput_T2D_CHM_s_SE_Geostatic &>(_self);
	Step_T2D_CHM_s_SE_Geostatic &step
		= static_cast<Step_T2D_CHM_s_SE_Geostatic &>(th.get_step());
	ResultFile_XML &rf = static_cast<ResultFile_XML &>(*th.res_file);
	std::fstream &file = rf.get_file();
	
	char str_buffer[512];
#define str_buffer_len (sizeof(str_buffer) / sizeof(str_buffer[0]))

	// time history
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


int time_history_output_func_t2d_chm_s_SE_geostatic_to_hdf5_res_file(TimeHistoryOutput &_self)
{
	TimeHistoryOutput_T2D_CHM_s_SE_Geostatic &th
		= static_cast<TimeHistoryOutput_T2D_CHM_s_SE_Geostatic &>(_self);
	Step_T2D_CHM_s_SE_Geostatic &step
		= static_cast<Step_T2D_CHM_s_SE_Geostatic &>(th.get_step());
	Model_T2D_CHM_s &md = static_cast<Model_T2D_CHM_s &>(step.get_model());
	ResultFile_hdf5 &rf = static_cast<ResultFile_hdf5 &>(*th.res_file);

	char frame_name[30];
	snprintf(frame_name, 30, "frame_%zu", th.output_id);
	hid_t frame_grp_id = rf.create_group(th.th_id, frame_name);
	rf.write_attribute(frame_grp_id, "current_time", step.get_current_time());
	rf.write_attribute(frame_grp_id, "total_time", step.get_total_time());
	rf.write_attribute(frame_grp_id, "substep_num", step.get_substep_num());
	rf.write_attribute(frame_grp_id, "total_substep_num", step.get_total_substep_num());
	// output particle data
	ouput_pcl_data_to_hdf5_file(md, rf, frame_grp_id);
	// output consititutive model
	output_model_container_to_hdf5_file(md.model_container, rf, frame_grp_id);
	// output rigid body
	output_rigid_ciricle_to_hdf5_file(md.get_rigid_circle(), rf, frame_grp_id);
	rf.close_group(frame_grp_id);

	++th.output_id;
	return 0;
}
