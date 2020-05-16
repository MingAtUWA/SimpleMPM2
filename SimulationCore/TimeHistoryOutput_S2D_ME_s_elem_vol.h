#ifndef __Time_History_Output_S2D_ME_s_elem_vol_h__
#define __Time_History_Output_S2D_ME_s_elem_vol_h__

#include "ResultFile_hdf5.h"
#include "ResultFile_XML.h"

#include "TimeHistoryOutput.h"

int time_history_output_func_s2d_me_s_elem_vol_to_hdf5_res_file(TimeHistoryOutput &_self);
int time_history_output_func_s2d_me_s_elem_vol_to_xml_res_file(TimeHistoryOutput &_self);

/* ===========================================================
Class TimeHistoryOutput_S2D_ME_s_elem_vol
=========================================================== */
class TimeHistoryOutput_S2D_ME_s_elem_vol : public TimeHistoryOutput
{
protected:
	size_t output_id;
	bool is_init;
	int init(void);
	int close(void);

public:
	TimeHistoryOutput_S2D_ME_s_elem_vol(const char *_name) :
		TimeHistoryOutput(_name, "S2D_ME_s"),
		output_id(0), is_init(false), th_id(-1) {}
	~TimeHistoryOutput_S2D_ME_s_elem_vol() { close(); }

	int init_per_step(void) { return init(); }
	inline void set_output_init_state(bool _need = true) noexcept {}

	friend int time_history_output_func_s2d_me_s_elem_vol_to_hdf5_res_file(TimeHistoryOutput &_self);
	inline void set_res_file(ResultFile_hdf5 &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &time_history_output_func_s2d_me_s_elem_vol_to_hdf5_res_file;
	}

protected:
	hid_t th_id;
public:
	friend int time_history_output_func_s2d_me_s_elem_vol_to_xml_res_file(TimeHistoryOutput &_self);
	inline void set_res_file(ResultFile_XML &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &time_history_output_func_s2d_me_s_elem_vol_to_xml_res_file;
	}
};

#endif