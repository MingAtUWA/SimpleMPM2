#ifndef __TIME_HISTORY_OUTPUT_S2D_ME_S_RIGIDBODY_H__
#define __TIME_HISTORY_OUTPUT_S2D_ME_S_RIGIDBODY_H__

#include "ResultFile_PlainBin.h"
#include "ResultFile_XML.h"

#include "TimeHistoryOutput.h"

int time_history_output_func_s2d_me_s_rigid_body_to_plain_bin_res_file(TimeHistoryOutput &_self);
int time_history_output_func_s2d_me_s_rigid_body_to_xml_res_file(TimeHistoryOutput &_self);
/* ===========================================================
Class TimeHistoryOutput_S2D_ME_s_RigidBody
=========================================================== */
class TimeHistoryOutput_S2D_ME_s_RigidBody : public TimeHistoryOutput
{
public:
	TimeHistoryOutput_S2D_ME_s_RigidBody(const char *_name = "TimeHistory") : 
		TimeHistoryOutput(_name, "S2D_ME_RigidBody") {}
	~TimeHistoryOutput_S2D_ME_s_RigidBody() {}

	friend int time_history_output_func_s2d_me_s_rigid_body_to_plain_bin_res_file(TimeHistoryOutput &_self);
	inline void set_res_file(ResultFile_PlainBin &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &time_history_output_func_s2d_me_s_rigid_body_to_plain_bin_res_file;
	}

	friend int time_history_output_func_s2d_me_s_rigid_body_to_xml_res_file(TimeHistoryOutput &_self);
	inline void set_res_file(ResultFile_XML &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &time_history_output_func_s2d_me_s_rigid_body_to_xml_res_file;
	}
};

#endif