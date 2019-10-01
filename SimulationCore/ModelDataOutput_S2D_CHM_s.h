#ifndef __MODEL_DATA_OUTPUT_S2D_CHM_S_H__
#define __MODEL_DATA_OUTPUT_S2D_CHM_S_H__

#include <string>

#include "ResultFile_PlainBin.h"
#include "ResultFile_XML.h"

#include "ModelDataOutput.h"

int model_data_output_func_s2d_chm_s_to_plain_bin_res_file(ModelDataOutput &_self);
int model_data_output_func_s2d_chm_s_to_xml_res_file(ModelDataOutput &_self);
/*=============================================================
Class ModelDataOutput_S2D_CHM_s
==============================================================*/
class ModelDataOutput_S2D_CHM_s : public ModelDataOutput
{
public:
	ModelDataOutput_S2D_CHM_s() : ModelDataOutput("ModelDataOutput_S2D_CHM_s") {}
	~ModelDataOutput_S2D_CHM_s() {}

	friend int model_data_output_func_s2d_chm_s_to_plain_bin_res_file(ModelDataOutput &_self);
	inline void set_res_file(ResultFile_PlainBin &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &model_data_output_func_s2d_chm_s_to_plain_bin_res_file;
	}

	friend int model_data_output_func_s2d_chm_s_to_xml_res_file(ModelDataOutput &_self);
	inline void set_res_file(ResultFile_XML &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &model_data_output_func_s2d_chm_s_to_xml_res_file;
	}
};

#endif