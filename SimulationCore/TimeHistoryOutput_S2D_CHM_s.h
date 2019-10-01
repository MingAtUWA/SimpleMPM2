#ifndef __TIME_HISTORY_OUTPUT_S2D_CHM_S_H__
#define __TIME_HISTORY_OUTPUT_S2D_CHM_S_H__

#include "ResultFile_PlainBin.h"
#include "ResultFile_XML.h"

#include "TimeHistoryOutput.h"

int time_history_output_func_s2d_chm_s_to_plain_bin_res_file(TimeHistoryOutput &_self);
int time_history_output_func_s2d_chm_s_to_xml_res_file(TimeHistoryOutput &_self);
/* ===========================================================
Class TimeHistoryOutput_S2D_CHM_s
=========================================================== */
class TimeHistoryOutput_S2D_CHM_s : public TimeHistoryOutput
{
public:
	TimeHistoryOutput_S2D_CHM_s() : TimeHistoryOutput("S2D_CHM_s") {}
	~TimeHistoryOutput_S2D_CHM_s() {}

	friend int time_history_output_func_s2d_chm_s_to_plain_bin_res_file(TimeHistoryOutput &_self);
	inline void set_res_file(ResultFile_PlainBin &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &time_history_output_func_s2d_chm_s_to_plain_bin_res_file;
	}

	friend int time_history_output_func_s2d_chm_s_to_xml_res_file(TimeHistoryOutput &_self);
	inline void set_res_file(ResultFile_XML &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &time_history_output_func_s2d_chm_s_to_xml_res_file;
	}
};

#endif