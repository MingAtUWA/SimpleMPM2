#ifndef __Time_History_Output_T2D_ME_s_Geostatic_H__
#define __Time_History_Output_T2D_ME_s_Geostatic_H__

#include "ResultFile_PlainBin.h"
#include "ResultFile_XML.h"
#include "ResultFile_hdf5.h"

#include "TimeHistoryOutput.h"

int time_history_output_func_t2d_me_s_geostatic_to_plain_bin_res_file(TimeHistoryOutput &_self);
int time_history_output_func_t2d_me_s_geostatic_to_xml_res_file(TimeHistoryOutput &_self);
int time_history_output_func_t2d_me_s_geostatic_to_hdf5_res_file(TimeHistoryOutput &_self);

/* ===========================================================
Class TimeHistoryOutput_T2D_ME_s_Geostatic
=========================================================== */
class TimeHistoryOutput_T2D_ME_s_Geostatic : public TimeHistoryOutput
{
protected:
	size_t output_id;
	bool is_init;
	int init(void);
	int close(void);

public:
	TimeHistoryOutput_T2D_ME_s_Geostatic(const char *_name) :
		TimeHistoryOutput(_name, "T2D_CHM_s_SE_Geostatic"),
		output_id(0), is_init(false),
		th_id(-1) {}
	~TimeHistoryOutput_T2D_ME_s_Geostatic() { close(); }

	int init_per_step(void) { return init(); }

	friend int time_history_output_func_t2d_me_s_geostatic_to_plain_bin_res_file(TimeHistoryOutput &_self);
	inline void set_res_file(ResultFile_PlainBin &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &time_history_output_func_t2d_me_s_geostatic_to_plain_bin_res_file;
	}

	friend int time_history_output_func_t2d_me_s_geostatic_to_xml_res_file(TimeHistoryOutput &_self);
	inline void set_res_file(ResultFile_XML &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &time_history_output_func_t2d_me_s_geostatic_to_xml_res_file;
	}

protected:
	hid_t th_id;
public:
	friend int time_history_output_func_t2d_me_s_geostatic_to_hdf5_res_file(TimeHistoryOutput &_self);
	inline void set_res_file(ResultFile_hdf5 &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &time_history_output_func_t2d_me_s_geostatic_to_hdf5_res_file;
	}
};

#endif