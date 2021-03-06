#ifndef __Time_History_Output_S2D_ME_s_Geostatic_h__
#define __Time_History_Output_S2D_ME_s_Geostatic_h__

#include "ResultFile_hdf5.h"
#include "ResultFile_XML.h"

#include "TimeHistoryOutput.h"

int time_history_output_func_s2d_me_s_geostatic_to_hdf5_res_file(TimeHistoryOutput &_self);
int time_history_output_func_s2d_me_s_geostatic_to_xml_res_file(TimeHistoryOutput &_self);

/* ===========================================================
Class TimeHistoryOutput_S2D_ME_s_Geostatic
=========================================================== */
class TimeHistoryOutput_S2D_ME_s_Geostatic : public TimeHistoryOutput
{
protected:
	size_t output_id;
	bool is_init;
	int init(void);
	int close(void);

public:
	TimeHistoryOutput_S2D_ME_s_Geostatic(const char *_name) :
		TimeHistoryOutput(_name, "S2D_ME_s_Geostatic"),
		output_id(0), is_init(false), th_id(-1) {}
	~TimeHistoryOutput_S2D_ME_s_Geostatic() { close(); }

	int init_per_step(void) { return init(); }

	friend int time_history_output_func_s2d_me_s_geostatic_to_hdf5_res_file(TimeHistoryOutput &_self);
	inline void set_res_file(ResultFile_hdf5 &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &time_history_output_func_s2d_me_s_geostatic_to_hdf5_res_file;
	}

protected:
	hid_t th_id;
public:
	friend int time_history_output_func_s2d_me_s_geostatic_to_xml_res_file(TimeHistoryOutput &_self);
	inline void set_res_file(ResultFile_XML &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &time_history_output_func_s2d_me_s_geostatic_to_xml_res_file;
	}
};

#endif