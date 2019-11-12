#ifndef __Model_Data_Output_S2D_ME_s_up_H__
#define __Model_Data_Output_S2D_ME_s_up_H__

#include <string>

#include "ResultFile_PlainBin.h"
#include "ResultFile_XML.h"

#include "ModelDataOutput.h"

int model_data_output_func_s2d_me_s_up_to_plain_bin_res_file(ModelDataOutput &_self);
int model_data_output_func_s2d_me_s_up_to_xml_res_file(ModelDataOutput &_self);
/*=============================================================
Class ModelDataOutput_S2D_ME_s_up
==============================================================*/
class ModelDataOutput_S2D_ME_s_up : public ModelDataOutput
{
public:
	ModelDataOutput_S2D_ME_s_up() : ModelDataOutput("ModelDataOutput_S2D_ME_s_up") {}
	~ModelDataOutput_S2D_ME_s_up() {}

	friend int model_data_output_func_s2d_me_s_up_to_plain_bin_res_file(ModelDataOutput &_self);
	inline void set_res_file(ResultFile_PlainBin &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &model_data_output_func_s2d_me_s_up_to_plain_bin_res_file;
	}

	friend int model_data_output_func_s2d_me_s_up_to_xml_res_file(ModelDataOutput &_self);
	inline void set_res_file(ResultFile_XML &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &model_data_output_func_s2d_me_s_up_to_xml_res_file;
	}
};

#endif