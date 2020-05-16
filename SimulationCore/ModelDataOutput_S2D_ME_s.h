#ifndef __Model_Data_Output_S2D_ME_s_h__
#define __Model_Data_Output_S2D_ME_s_h__

#include <string>

#include "ResultFile_hdf5.h"
#include "ResultFile_XML.h"

#include "ModelDataOutput.h"

int model_data_output_func_s2d_me_s_geostatic_to_hdf5_res_file(ModelDataOutput &_self);
int model_data_output_func_s2d_me_s_geostatic_to_xml_res_file(ModelDataOutput &_self);

/*=============================================================
Class ModelDataOutput_S2D_ME_s
==============================================================*/
class ModelDataOutput_S2D_ME_s : public ModelDataOutput
{
public:
	ModelDataOutput_S2D_ME_s(const char *_name) :
		ModelDataOutput(_name, "ModelDataOutput_S2D_ME_s") {}
	~ModelDataOutput_S2D_ME_s() {}

	friend int model_data_output_func_s2d_me_s_geostatic_to_hdf5_res_file(ModelDataOutput &_self);
	inline void set_res_file(ResultFile_hdf5 &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &model_data_output_func_s2d_me_s_geostatic_to_hdf5_res_file;
	}

	friend int model_data_output_func_s2d_me_s_geostatic_to_xml_res_file(ModelDataOutput &_self);
	inline void set_res_file(ResultFile_XML &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &model_data_output_func_s2d_me_s_geostatic_to_xml_res_file;
	}
};

#endif