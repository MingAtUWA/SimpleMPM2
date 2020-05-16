#include "SimulationCore_pcp.h"

#include "Model_S2D_ME_s.h"
#include "Step_S2D_ME_s_Geostatic.h"
#include "Model_S2D_ME_s_hdf5_io_utilities.h"

#include "TimeHistoryOutput_S2D_ME_s_Geostatic.h"

int TimeHistoryOutput_S2D_ME_s_Geostatic::init(void)
{
	if (is_init)
		return 0;
	is_init = true;

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

int TimeHistoryOutput_S2D_ME_s_Geostatic::close(void)
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

	is_init = false;
	return 0;
}

int time_history_output_func_s2d_me_s_geostatic_to_hdf5_res_file(TimeHistoryOutput &_self)
{
	TimeHistoryOutput_S2D_ME_s_Geostatic &th
		= static_cast<TimeHistoryOutput_S2D_ME_s_Geostatic &>(_self);
	Step_S2D_ME_s_Geostatic &step = static_cast<Step_S2D_ME_s_Geostatic &>(th.get_step());
	Model_S2D_ME_s &md = static_cast<Model_S2D_ME_s &>(step.get_model());
	ResultFile_hdf5 &rf = static_cast<ResultFile_hdf5 &>(*th.res_file);

	char frame_name[30];
	snprintf(frame_name, 30, "frame_%zu", th.output_id);
	hid_t frame_grp_id = rf.create_group(th.th_id, frame_name);
	rf.write_attribute(frame_grp_id, "current_time", step.get_current_time());
	rf.write_attribute(frame_grp_id, "total_time", step.get_total_time());
	rf.write_attribute(frame_grp_id, "substep_num", step.get_substep_num());
	rf.write_attribute(frame_grp_id, "total_substep_num", step.get_total_substep_num());
	// output particle data
	using Model_S2D_ME_s_hdf5_io_utilities::output_pcl_data_to_hdf5_file;
	output_pcl_data_to_hdf5_file(md, rf, frame_grp_id);

	++th.output_id;
	return 0;
}

int time_history_output_func_s2d_me_s_geostatic_to_xml_res_file(TimeHistoryOutput &_self)
{
	TimeHistoryOutput_S2D_ME_s_Geostatic &th
		= static_cast<TimeHistoryOutput_S2D_ME_s_Geostatic &>(_self);
	Step_S2D_ME_s_Geostatic &step = static_cast<Step_S2D_ME_s_Geostatic &>(th.get_step());
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
		step.get_substep_num(), step.get_total_substep_num(),
		step.get_current_time(), step.get_total_time());
	file.write(str_buffer, strlen(str_buffer));

	// output material points data
	Model_S2D_ME_s &model = static_cast<Model_S2D_ME_s &>(th.get_model());
	size_t pcl_num = model.get_pcl_num();
	Model_S2D_ME_s::Particle *pcls = model.get_pcls();
	const char *material_point_info = ""
		"    <MaterialPointObject>\n"
		"        <pcl_num> %zu </pcl_num>\n";
	snprintf(str_buffer, str_buffer_len, material_point_info, pcl_num);
	file.write(str_buffer, strlen(str_buffer));
	// field data: x, y, vol
	file << "        <field_data>\n"
		"        <!-- x, y, vol, s11, s22, s12 -->\n";
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		Model_S2D_ME_s::Particle &pcl = pcls[pcl_id];
		file << "        ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", pcl.x);
		file << str_buffer << ", ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", pcl.y);
		file << str_buffer << ", ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", pcl.m / pcl.density);
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
