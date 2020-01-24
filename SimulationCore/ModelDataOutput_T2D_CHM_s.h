#ifndef __Model_Data_Output_T2D_CHM_s_H__
#define __Model_Data_Output_T2D_CHM_s_H__

#include <string>

#include "ResultFile_PlainBin.h"
#include "ResultFile_XML.h"
#include "ResultFile_hdf5.h"

#include "ModelDataOutput.h"

int model_data_output_func_t2d_chm_s_to_plain_bin_res_file(ModelDataOutput &_self);
int model_data_output_func_t2d_chm_s_to_xml_res_file(ModelDataOutput &_self);
int model_data_output_func_t2d_chm_s_to_hdf5_res_file(ModelDataOutput &_self);
/*=============================================================
Class ModelDataOutput_T2D_CHM_s
==============================================================*/
class ModelDataOutput_T2D_CHM_s : public ModelDataOutput
{
public:
	ModelDataOutput_T2D_CHM_s(const char *_name) :
		ModelDataOutput(_name, "ModelDataOutput_T2D_CHM_s") {}
	~ModelDataOutput_T2D_CHM_s() {}

	friend int model_data_output_func_t2d_chm_s_to_plain_bin_res_file(ModelDataOutput &_self);
	inline void set_res_file(ResultFile_PlainBin &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &model_data_output_func_t2d_chm_s_to_plain_bin_res_file;
	}

	friend int model_data_output_func_t2d_chm_s_to_xml_res_file(ModelDataOutput &_self);
	inline void set_res_file(ResultFile_XML &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &model_data_output_func_t2d_chm_s_to_xml_res_file;
	}

	friend int model_data_output_func_t2d_chm_s_to_hdf5_res_file(ModelDataOutput &_self);
	inline void set_res_file(ResultFile_hdf5 &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &model_data_output_func_t2d_chm_s_to_hdf5_res_file;
	}
};

#endif