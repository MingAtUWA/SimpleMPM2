#ifndef __MODEL_DATA_OUTPUT_S2D_ME_S_RIGID_BODY_FRIC_H__
#define __MODEL_DATA_OUTPUT_S2D_ME_S_RIGID_BODY_FRIC_H__

#include <string>

#include "ResultFile_PlainBin.h"
#include "ResultFile_XML.h"

#include "ModelDataOutput.h"

int model_data_output_func_s2d_me_s_rigid_body_fric_to_plain_bin_res_file(ModelDataOutput &_self);
int model_data_output_func_s2d_me_s_rigid_body_fric_to_xml_res_file(ModelDataOutput &_self);
/*=============================================================
Class ModelDataOutput_S2D_ME_s_RigidBody_Fric
==============================================================*/
class ModelDataOutput_S2D_ME_s_RigidBody_Fric : public ModelDataOutput
{
public:
	ModelDataOutput_S2D_ME_s_RigidBody_Fric() :
		ModelDataOutput("ModelDataOutput_S2D_ME_s_RigidBody") {}
	~ModelDataOutput_S2D_ME_s_RigidBody_Fric() {}

	friend int model_data_output_func_s2d_me_s_rigid_body_fric_to_plain_bin_res_file(ModelDataOutput &_self);
	inline void set_res_file(ResultFile_PlainBin &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &model_data_output_func_s2d_me_s_rigid_body_fric_to_plain_bin_res_file;
	}

	friend int model_data_output_func_s2d_me_s_rigid_body_fric_to_xml_res_file(ModelDataOutput &_self);
	inline void set_res_file(ResultFile_XML &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &model_data_output_func_s2d_me_s_rigid_body_fric_to_xml_res_file;
	}
};

#endif