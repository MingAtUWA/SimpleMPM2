#ifndef __Time_History_Output_T2D_CHM_s_SE_Geostatic_H__
#define __Time_History_Output_T2D_CHM_s_SE_Geostatic_H__

#include "ResultFile_PlainBin.h"
#include "ResultFile_XML.h"

#include "TimeHistoryOutput.h"

int time_history_output_func_t2d_chm_s_SE_geostatic_to_plain_bin_res_file(TimeHistoryOutput &_self);
int time_history_output_func_t2d_chm_s_SE_geostatic_to_xml_res_file(TimeHistoryOutput &_self);
/* ===========================================================
Class TimeHistoryOutput_T2D_CHM_s_SE_Geostatic
=========================================================== */
class TimeHistoryOutput_T2D_CHM_s_SE_Geostatic : public TimeHistoryOutput
{
public:
	TimeHistoryOutput_T2D_CHM_s_SE_Geostatic() : TimeHistoryOutput("T2D_CHM_s_SE_Geostatic") {}
	~TimeHistoryOutput_T2D_CHM_s_SE_Geostatic() {}

	friend int time_history_output_func_t2d_chm_s_SE_geostatic_to_plain_bin_res_file(TimeHistoryOutput &_self);
	inline void set_res_file(ResultFile_PlainBin &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &time_history_output_func_t2d_chm_s_SE_geostatic_to_plain_bin_res_file;
	}

	friend int time_history_output_func_t2d_chm_s_SE_geostatic_to_xml_res_file(TimeHistoryOutput &_self);
	inline void set_res_file(ResultFile_XML &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &time_history_output_func_t2d_chm_s_SE_geostatic_to_xml_res_file;
	}
};

#endif