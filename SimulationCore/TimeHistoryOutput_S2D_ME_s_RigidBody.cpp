#include "SimulationCore_pcp.h"

#include "Model_S2D_ME_s_RigidBody.h"
#include "Step_S2D_ME_s_RigidBody.h"

#include "ResultFile_PlainBin_DataStruct.h"

#include "TimeHistoryOutput_S2D_ME_s_RigidBody.h"

int time_history_output_func_s2d_me_s_rigid_body_to_plain_bin_res_file(TimeHistoryOutput &_self)
{
	TimeHistoryOutput_S2D_ME_s_RigidBody &th
		= static_cast<TimeHistoryOutput_S2D_ME_s_RigidBody &>(_self);
	ResultFile_PlainBin &rf = static_cast<ResultFile_PlainBin &>(*th.res_file);
	std::fstream &file = rf.get_file();

	typedef ResultFile_PlainBin_DataStruct::TimeHistoryHeader TimeHistoryHeader;
	typedef ResultFile_PlainBin_DataStruct::RigidBodyMotionHeader RigidBodyMotionHeader;
	typedef ResultFile_PlainBin_DataStruct::MPObjectHeader MPObjectHeader;

	TimeHistoryHeader thh;
	Step_S2D_ME_s_RigidBody &step = static_cast<Step_S2D_ME_s_RigidBody &>(th.get_step());
	thh.substep_num = step.get_substep_num();
	thh.total_substep_num = step.get_total_substep_num();
	thh.current_time = step.get_current_time();
	thh.total_time = step.get_total_time();
	file.write(reinterpret_cast<char *>(&thh), sizeof(thh));

	RigidBodyMotionHeader rbmh;
	Model_S2D_ME_s_RigidBody &model = static_cast<Model_S2D_ME_s_RigidBody &>(th.get_model());
	RigidBody &rb = model.rigid_body;
	rbmh.x = rb.x;
	rbmh.y = rb.y;
	rbmh.theta = rb.theta;
	rbmh.vx = rb.vx;
	rbmh.vy = rb.vy;
	rbmh.v_theta = rb.vtheta;
	rbmh.fx_con = rb.Fx_con;
	rbmh.fy_con = rb.Fy_con;
	rbmh.m_con = rb.M_con;
	file.write(reinterpret_cast<char *>(&rbmh), sizeof(rbmh));

	// output particles data
	MPObjectHeader mph;
	mph.pcl_num = model.pcl_num;
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
		pcl_data[pcl_id] = model.pcls[pcl_id].m / model.pcls[pcl_id].density;
	file.write(reinterpret_cast<char *>(pcl_data), data_len);
	delete[] pcl_data;

	return 0;
}

int time_history_output_func_s2d_me_s_rigid_body_to_xml_res_file(TimeHistoryOutput &_self)
{
	TimeHistoryOutput_S2D_ME_s_RigidBody &th
		= static_cast<TimeHistoryOutput_S2D_ME_s_RigidBody &>(_self);
	ResultFile_XML &rf = static_cast<ResultFile_XML &>(*th.res_file);
	std::fstream &file = rf.get_file();
	
	char str_buffer[512];

	// time history
	Step_S2D_ME_s_RigidBody &step = static_cast<Step_S2D_ME_s_RigidBody &>(th.get_step());
	const char *time_history_info = ""
		"<TimeHistory>\n"
		"    <substep_num> %zu </substep_num>\n"
		"    <total_substep_num> %zu </total_substep_num>\n"
		"    <current_time> %16.10e </current_time>\n"
		"    <total_time> %16.10e </total_time>\n";
	snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]), time_history_info,
		step.get_substep_num(), step.get_total_substep_num(),
		step.get_current_time(), step.get_total_time());
	file.write(str_buffer, strlen(str_buffer));

	// rigid object
	Model_S2D_ME_s_RigidBody &model = static_cast<Model_S2D_ME_s_RigidBody &>(th.get_model());
	RigidBody &rb = model.rigid_body;
	const char *rigid_object_info = ""
		"    <RigidObject>\n"
		"        <x> %16.10e </x>\n"
		"        <y> %16.10e </y>\n"
		"        <theta> %16.10e </theta>\n"
		"        <vx> %16.10e </vx>\n"
		"        <vy> %16.10e </vy>\n"
		"        <vtheta> %16.10e </vtheta>\n"
		"        <Fx_contact> %16.10e </Fx_contact>\n"
		"        <Fy_contact> %16.10e </Fy_contact>\n"
		"        <M_contact> %16.10e </M_contact>\n"
		"    </RigidObject>\n";
	// resistence caused by contact
	snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]), rigid_object_info,
		rb.x, rb.y, rb.theta, rb.vx, rb.vy, rb.vtheta, rb.Fx_con, rb.Fy_con, rb.M_con);
	file.write(str_buffer, strlen(str_buffer));

	// output material points data
	const char *material_point_info = ""
		"    <MaterialPointObject>\n"
		"        <pcl_num> %zu </pcl_num>\n";
	snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]), material_point_info, model.pcl_num);
	file.write(str_buffer, strlen(str_buffer));
	// field data: x, y, vol
	file << "        <field_data>\n"
		"        <!-- x, y, vol -->\n";
	for (size_t pcl_id = 0; pcl_id < model.pcl_num; ++pcl_id)
	{
		file << "            " << model.pcls[pcl_id].x << ", " << model.pcls[pcl_id].y << ", "
			<< model.pcls[pcl_id].m / model.pcls[pcl_id].density << "\n";
	}
	file << "        </field_data>\n";
	// ending
	file << "    </MaterialPointObject>\n";

	// ending
	file << "</TimeHistory>\n";

	return 0;
}
