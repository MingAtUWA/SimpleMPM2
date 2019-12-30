#include "SimulationCore_pcp.h"

#include "Model_T2D_ME_s.h"
#include "Step_T2D_ME_s.h"

#include "ResultFile_PlainBin_DataStruct.h"

#include "TimeHistoryOutput_T2D_ME_s.h"

int time_history_output_func_t2d_me_s_to_plain_bin_res_file(TimeHistoryOutput &_self)
{
	TimeHistoryOutput_T2D_ME_s &th = static_cast<TimeHistoryOutput_T2D_ME_s &>(_self);
	Model_T2D_ME_s &model = static_cast<Model_T2D_ME_s &>(th.get_model());
	Step_T2D_ME_s &step = static_cast<Step_T2D_ME_s &>(th.get_step());
	ResultFile_PlainBin &rf = static_cast<ResultFile_PlainBin &>(*th.res_file);
	std::fstream &file = rf.get_file();
	double *pcl_data;
	size_t data_len;

	typedef ResultFile_PlainBin_DataStruct::TimeHistoryHeader TimeHistoryHeader;
	typedef ResultFile_PlainBin_DataStruct::DispConRigidCircleMotionHeader DispConRigidCircleMotionHeader;
	typedef ResultFile_PlainBin_DataStruct::MPObjectHeader MPObjectHeader;
	
	// time history header
	TimeHistoryHeader thh;
	thh.init();
	thh.substep_num = step.get_substep_num();
	thh.total_substep_num = step.get_total_substep_num();
	thh.current_time = step.get_current_time();
	thh.total_time = step.get_total_time();
	file.write(reinterpret_cast<char *>(&thh), sizeof(thh));

	// rigid circle data
	Model_T2D_ME_s &md = static_cast<Model_T2D_ME_s &>(step.get_model());
	DispConRigidCircle &rc = md.get_rigid_circle();
	size_t rc_pcl_num = rc.get_pcl_num();
	if (rc_pcl_num != 0)
	{
		DispConRigidCircle::State state = rc.get_state();
		DispConRigidCircle::Particle *rb_pcls = state.pcls;
		DispConRigidCircleMotionHeader rcmh;
		rcmh.init();
		rcmh.x = state.cen_x;
		rcmh.y = state.cen_y;
		rcmh.theta = state.theta;
		rcmh.vx = state.vx;
		rcmh.vy = state.vy;
		rcmh.w = state.w;
		rcmh.rfx = state.rfx;
		rcmh.rfy = state.rfy;
		rcmh.rm = state.rm;
		file.write(reinterpret_cast<char *>(&rcmh), sizeof(rcmh));
		// rigid body pcl data
		pcl_data = new double[rc_pcl_num];
		data_len = sizeof(double) * rc_pcl_num;
		// x
		for (size_t pcl_id = 0; pcl_id < rc_pcl_num; ++pcl_id)
			pcl_data[pcl_id] = rb_pcls[pcl_id].x;
		file.write(reinterpret_cast<char *>(pcl_data), data_len);
		// y
		for (size_t pcl_id = 0; pcl_id < rc_pcl_num; ++pcl_id)
			pcl_data[pcl_id] = rb_pcls[pcl_id].y;
		file.write(reinterpret_cast<char *>(pcl_data), data_len);
		// vol
		for (size_t pcl_id = 0; pcl_id < rc_pcl_num; ++pcl_id)
			pcl_data[pcl_id] = rb_pcls[pcl_id].vol;
		file.write(reinterpret_cast<char *>(pcl_data), data_len);
		delete[] pcl_data;
	}

	// output particles data
	MPObjectHeader mph;
	mph.init();
	mph.pcl_num = model.pcl_num;
	mph.fld_num = 6;
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
		Model_T2D_ME_s::Particle &pcl = model.pcls[pcl_id];
		pcl_data[pcl_id] = pcl.m / pcl.density;
	}
	file.write(reinterpret_cast<char *>(pcl_data), data_len);
	// s11
	for (size_t pcl_id = 0; pcl_id < mph.pcl_num; ++pcl_id)
		pcl_data[pcl_id] = model.pcls[pcl_id].s11;
	file.write(reinterpret_cast<char *>(pcl_data), data_len);
	// s22
	for (size_t pcl_id = 0; pcl_id < mph.pcl_num; ++pcl_id)
		pcl_data[pcl_id] = model.pcls[pcl_id].s22;
	file.write(reinterpret_cast<char *>(pcl_data), data_len);
	// s12
	for (size_t pcl_id = 0; pcl_id < mph.pcl_num; ++pcl_id)
		pcl_data[pcl_id] = model.pcls[pcl_id].s12;
	file.write(reinterpret_cast<char *>(pcl_data), data_len);
	delete[] pcl_data;

	return 0;
}

int time_history_output_func_t2d_me_s_to_xml_res_file(TimeHistoryOutput &_self)
{
	TimeHistoryOutput_T2D_ME_s &th = static_cast<TimeHistoryOutput_T2D_ME_s &>(_self);
	ResultFile_XML &rf = static_cast<ResultFile_XML &>(*th.res_file);
	std::fstream &file = rf.get_file();
	
	char str_buffer[512];
#define str_buffer_len (sizeof(str_buffer) / sizeof(str_buffer[0]))

	// time history
	Step_T2D_ME_s &step = static_cast<Step_T2D_ME_s &>(th.get_step());
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

	// output rigid body data
	Model_T2D_ME_s &md = static_cast<Model_T2D_ME_s &>(step.get_model());
	DispConRigidCircle &rb = md.get_rigid_circle();
	DispConRigidCircle::State state = rb.get_state();
	const char *rigid_body_info = ""
			"    <RigidBody>\n"
			"        <x> %16.10e </x>\n"
			"        <y> %16.10e </y>\n"
			"        <theta> %16.10e </theta>\n"
			"        <vx> %16.10e </vx>\n"
			"        <vy> %16.10e </vy>\n"
			"        <w> %16.10e </w>\n"
			"        <rfx> %16.10e </rfx>\n"
			"        <rfy> %16.10e </rfy>\n"
			"        <rm> %16.10e </rm>\n"
			"        <pcl_num> %zu </pcl_num>\n"
			"        <pcl_data>\n"
			"        <!-- x, y, vol, vx, vy -->\n";
	snprintf(str_buffer, str_buffer_len, rigid_body_info,
			 state.cen_x, state.cen_y, state.theta,
			 state.vx, state.vy, state.w, 
			 state.rfx, state.rfy, state.rm, state.pcl_num);
	file.write(str_buffer, strlen(str_buffer));
	for (size_t pcl_id = 0; pcl_id < state.pcl_num; ++pcl_id)
	{
		DispConRigidCircle::Particle &pcl = state.pcls[pcl_id];
		file << "        ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", pcl.x);
		file << str_buffer << ", ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", pcl.y);
		file << str_buffer << ", ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", pcl.vol);
		file << str_buffer << ", ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", pcl.vx);
		file << str_buffer << ", ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", pcl.vy);
		file << str_buffer << "\n";
	}
	file << "        </pcl_data>\n"
			"    </RigidBody>\n";

	// output material points data
	Model_T2D_ME_s &model = static_cast<Model_T2D_ME_s &>(th.get_model());
	const char *material_point_info = ""
		"    <MaterialPointObject>\n"
		"        <pcl_num> %zu </pcl_num>\n";
	snprintf(str_buffer, str_buffer_len, material_point_info, model.pcl_num);
	file.write(str_buffer, strlen(str_buffer));
	// field data: x, y, vol
	file << "        <field_data>\n"
			"        <!-- x, y, vol, s11, s22, s12 -->\n";
	for (size_t pcl_id = 0; pcl_id < model.pcl_num; ++pcl_id)
	{
		Model_T2D_ME_s::Particle &pcl = model.pcls[pcl_id];
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
