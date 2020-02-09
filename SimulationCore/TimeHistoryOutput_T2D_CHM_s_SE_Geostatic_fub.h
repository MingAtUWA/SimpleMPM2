#ifndef __Time_History_Output_T2D_CHM_s_SE_Geostatic_fub_H__
#define __Time_History_Output_T2D_CHM_s_SE_Geostatic_fub_H__

#include "ResultFile_PlainBin.h"
#include "ResultFile_XML.h"
#include "ResultFile_hdf5.h"

#include "TimeHistoryOutput.h"

int time_history_output_func_t2d_chm_s_SE_geostatic_fub_to_plain_bin_res_file(TimeHistoryOutput &_self);
int time_history_output_func_t2d_chm_s_SE_geostatic_fub_to_xml_res_file(TimeHistoryOutput &_self);
int time_history_output_func_t2d_chm_s_SE_geostatic_fub_to_hdf5_res_file(TimeHistoryOutput &_self);

/* ===========================================================
Class TimeHistoryOutput_T2D_CHM_s_SE_Geostatic_fub
	1. only output to xml and hdf5 file
	2. output unbalance nodal force to check geostatic condition
=========================================================== */
class TimeHistoryOutput_T2D_CHM_s_SE_Geostatic_fub : public TimeHistoryOutput
{
protected:
	size_t output_id;
	bool is_init;
	int init(void);
	int close(void);

public:
	TimeHistoryOutput_T2D_CHM_s_SE_Geostatic_fub(const char *_name) :
		TimeHistoryOutput(_name, "T2D_CHM_s_SE_Geostatic_fub"),
		output_id(0), is_init(false),
		th_id(-1) {}
	~TimeHistoryOutput_T2D_CHM_s_SE_Geostatic_fub() { close(); }

	int init_per_step(void) { return init(); }

	friend int time_history_output_func_t2d_chm_s_SE_geostatic_fub_to_plain_bin_res_file(TimeHistoryOutput &_self);
	inline void set_res_file(ResultFile_PlainBin &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &time_history_output_func_t2d_chm_s_SE_geostatic_fub_to_plain_bin_res_file;
	}

	friend int time_history_output_func_t2d_chm_s_SE_geostatic_fub_to_xml_res_file(TimeHistoryOutput &_self);
	inline void set_res_file(ResultFile_XML &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &time_history_output_func_t2d_chm_s_SE_geostatic_fub_to_xml_res_file;
	}

protected:
	hid_t th_id;
public:
	friend int time_history_output_func_t2d_chm_s_SE_geostatic_fub_to_hdf5_res_file(TimeHistoryOutput &_self);
	inline void set_res_file(ResultFile_hdf5 &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &time_history_output_func_t2d_chm_s_SE_geostatic_fub_to_hdf5_res_file;
	}
};

#endif