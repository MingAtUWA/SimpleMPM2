#ifndef __Time_History_Output_S2D_CHM_s_FEM_uUp_H__
#define __Time_History_Output_S2D_CHM_s_FEM_uUp_H__

#include "ResultFile_PlainBin.h"
#include "ResultFile_XML.h"

#include "TimeHistoryOutput.h"

int time_history_output_func_s2d_chm_s_fem_uup_to_plain_bin_res_file(TimeHistoryOutput &_self);
int time_history_output_func_s2d_chm_s_fem_uup_to_xml_res_file(TimeHistoryOutput &_self);
/* ===========================================================
Class TimeHistoryOutput_S2D_CHM_s_FEM_uUp
=========================================================== */
class TimeHistoryOutput_S2D_CHM_s_FEM_uUp : public TimeHistoryOutput
{
public:
	TimeHistoryOutput_S2D_CHM_s_FEM_uUp() : TimeHistoryOutput("S2D_CHM_s") {}
	~TimeHistoryOutput_S2D_CHM_s_FEM_uUp() {}

	friend int time_history_output_func_s2d_chm_s_fem_uup_to_plain_bin_res_file(TimeHistoryOutput &_self);
	inline void set_res_file(ResultFile_PlainBin &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &time_history_output_func_s2d_chm_s_fem_uup_to_plain_bin_res_file;
	}

	friend int time_history_output_func_s2d_chm_s_fem_uup_to_xml_res_file(TimeHistoryOutput &_self);
	inline void set_res_file(ResultFile_XML &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &time_history_output_func_s2d_chm_s_fem_uup_to_xml_res_file;
	}
};

#endif