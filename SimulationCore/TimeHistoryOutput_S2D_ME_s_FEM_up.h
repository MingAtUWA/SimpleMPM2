#ifndef __Time_History_Output_S2D_ME_s_FEM_up_H__
#define __Time_History_Output_S2D_ME_s_FEM_up_H__

#include "ResultFile_PlainBin.h"
#include "ResultFile_XML.h"

#include "TimeHistoryOutput.h"

int time_history_output_func_s2d_me_s_fem_up_to_plain_bin_res_file(TimeHistoryOutput &_self);
int time_history_output_func_s2d_me_s_fem_up_to_xml_res_file(TimeHistoryOutput &_self);
/* ===========================================================
Class TimeHistoryOutput_S2D_ME_s_FEM_up
=========================================================== */
class TimeHistoryOutput_S2D_ME_s_FEM_up : public TimeHistoryOutput
{
protected:
	double *node_data;
	double *gp_data;
	bool is_plain_bin;

public:
	TimeHistoryOutput_S2D_ME_s_FEM_up(const char *_name = "TimeHistory") :
		TimeHistoryOutput(_name, "S2D_ME_s_up"),
		node_data(nullptr), gp_data(nullptr),
		is_plain_bin(false) {}
	~TimeHistoryOutput_S2D_ME_s_FEM_up() {}

	// Initialize before each steps
	int init_per_step(void) override;
	// Finalize after each steps
	void finalize_per_step(void) override;

	friend int time_history_output_func_s2d_me_s_fem_up_to_plain_bin_res_file(TimeHistoryOutput &_self);
	inline void set_res_file(ResultFile_PlainBin &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &time_history_output_func_s2d_me_s_fem_up_to_plain_bin_res_file;
		is_plain_bin = true;
	}

	friend int time_history_output_func_s2d_me_s_fem_up_to_xml_res_file(TimeHistoryOutput &_self);
	inline void set_res_file(ResultFile_XML &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &time_history_output_func_s2d_me_s_fem_up_to_xml_res_file;
	}
};

#endif