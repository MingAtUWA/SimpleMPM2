#ifndef __Time_History_Output_T2D_CHM_s_H__
#define __Time_History_Output_T2D_CHM_s_H__

#include "ResultFile_PlainBin.h"
#include "ResultFile_XML.h"

#include "TimeHistoryOutput.h"

int time_history_output_func_t2d_chm_s_to_plain_bin_res_file(TimeHistoryOutput &_self);
int time_history_output_func_t2d_chm_s_to_xml_res_file(TimeHistoryOutput &_self);
/* ===========================================================
Class TimeHistoryOutput_T2D_CHM_s
=========================================================== */
class TimeHistoryOutput_T2D_CHM_s : public TimeHistoryOutput
{
public:
	TimeHistoryOutput_T2D_CHM_s(const char *_name = "TimeHistory") :
		TimeHistoryOutput(_name, "T2D_CHM_s") {}
	~TimeHistoryOutput_T2D_CHM_s() {}

	friend int time_history_output_func_t2d_chm_s_to_plain_bin_res_file(TimeHistoryOutput &_self);
	inline void set_res_file(ResultFile_PlainBin &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &time_history_output_func_t2d_chm_s_to_plain_bin_res_file;
	}

	friend int time_history_output_func_t2d_chm_s_to_xml_res_file(TimeHistoryOutput &_self);
	inline void set_res_file(ResultFile_XML &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &time_history_output_func_t2d_chm_s_to_xml_res_file;
	}
};

#endif