#include "SimulationCore_pcp.h"

#include "Model_T2D_fluid.h"
#include "Step_T2D_fluid.h"

#include "Model_T2D_fluid_hdf5_io_utilities.h"

#include "TimeHistoryOutput_T2D_fluid.h"

int TimeHistoryOutput_T2D_fluid::init(void)
{
	if (is_init) return 0;
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

int TimeHistoryOutput_T2D_fluid::close(void)
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

int time_history_output_func_t2d_fluid_to_xml_res_file(TimeHistoryOutput &_self)
{
	TimeHistoryOutput_T2D_fluid &th = static_cast<TimeHistoryOutput_T2D_fluid &>(_self);
	ResultFile_XML &rf = static_cast<ResultFile_XML &>(*th.res_file);
	std::fstream &file = rf.get_file();
	
	char str_buffer[512];
#define str_buffer_len (sizeof(str_buffer) / sizeof(str_buffer[0]))

	// time history
	Step_T2D_fluid &step = static_cast<Step_T2D_fluid &>(th.get_step());
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
	Model_T2D_fluid &model = static_cast<Model_T2D_fluid &>(th.get_model());
	const char *material_point_info = ""
		"    <MaterialPointObject>\n"
		"        <pcl_num> %zu </pcl_num>\n";
	snprintf(str_buffer, str_buffer_len, material_point_info, model.pcl_num);
	file.write(str_buffer, strlen(str_buffer));
	// field data: x, y, vol
	file << "        <field_data>\n"
			"        <!-- x, y, vol, t11, t22, t12, p -->\n";
	for (size_t pcl_id = 0; pcl_id < model.pcl_num; ++pcl_id)
	{
		Model_T2D_fluid::Particle &pcl = model.pcls[pcl_id];
		file << "        ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", pcl.x);
		file << str_buffer << ", ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", pcl.y);
		file << str_buffer << ", ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", pcl.m / pcl.density);
		file << str_buffer << ", ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", pcl.t11);
		file << str_buffer << ", ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", pcl.t22);
		file << str_buffer << ", ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", pcl.t12);
		file << str_buffer << ", ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", pcl.p);
		file << str_buffer << "\n";
	}
	file << "        </field_data>\n"
			"    </MaterialPointObject>\n";

	// ending
	file << "</TimeHistory>\n";

	return 0;
}

int time_history_output_func_t2d_fluid_to_hdf5_res_file(TimeHistoryOutput &_self)
{
	TimeHistoryOutput_T2D_fluid &th
		= static_cast<TimeHistoryOutput_T2D_fluid &>(_self);
	Step_T2D_fluid &step = static_cast<Step_T2D_fluid &>(th.get_step());
	Model_T2D_fluid &md = static_cast<Model_T2D_fluid &>(step.get_model());
	ResultFile_hdf5 &rf = static_cast<ResultFile_hdf5 &>(*th.res_file);

	char frame_name[30];
	snprintf(frame_name, 30, "frame_%zu", th.output_id);
	hid_t frame_grp_id = rf.create_group(th.th_id, frame_name);
	rf.write_attribute(frame_grp_id, "current_time", step.get_current_time());
	rf.write_attribute(frame_grp_id, "total_time", step.get_total_time());
	rf.write_attribute(frame_grp_id, "substep_num", step.get_substep_num());
	rf.write_attribute(frame_grp_id, "total_substep_num", step.get_total_substep_num());
	// output particle data
	using Model_T2D_fluid_hdf5_io_utilities::output_pcl_data_to_hdf5_file;
	output_pcl_data_to_hdf5_file(md, rf, frame_grp_id);

	++th.output_id;
	return 0;

}
