#include "SimulationCore_pcp.h"

#include "Model_S2D_ME_s_RigidBody.h"
#include "Step_S2D_ME_s_RigidBody.h"

#include "ResultFile_PlainBin_DataStruct.h"

#include "ModelDataOutput_S2D_ME_s_RigidBody.h"

int model_data_output_func_s2d_me_s_rigid_body_to_plain_bin_res_file(ModelDataOutput &_self)
{
	ModelDataOutput_S2D_ME_s_RigidBody &md 
		= static_cast<ModelDataOutput_S2D_ME_s_RigidBody &>(_self);
	Model_S2D_ME_s_RigidBody &model = static_cast<Model_S2D_ME_s_RigidBody &>(*md.model);
	ResultFile_PlainBin &rf = static_cast<ResultFile_PlainBin &>(*md.res_file);
	std::fstream &file = rf.get_file();

	typedef ResultFile_PlainBin_DataStruct::MeshHeader MeshHeader;
	typedef ResultFile_PlainBin_DataStruct::MPObjectHeader MPObjectHeader;

	// mesh
	MeshHeader mh;
	mh.h = model.h;
	mh.x0 = model.x0;
	mh.xn = model.xn;
	mh.y0 = model.y0;
	mh.yn = model.yn;
	mh.elem_x_num = model.elem_x_num;
	mh.elem_y_num = model.elem_y_num;
	file.write(reinterpret_cast<char *>(&mh), sizeof(mh));

	// material point object
	MPObjectHeader mph;
	mph.pcl_num = model.pcl_num;
	file.write(reinterpret_cast<char *>(&mph), sizeof(mph));

	return 0;
}

int model_data_output_func_s2d_me_s_rigid_body_to_xml_res_file(ModelDataOutput &_self)
{
	ModelDataOutput_S2D_ME_s_RigidBody &md
		= static_cast<ModelDataOutput_S2D_ME_s_RigidBody &>(_self);
	Model_S2D_ME_s_RigidBody &model = static_cast<Model_S2D_ME_s_RigidBody &>(*md.model);
	ResultFile_XML &rf = static_cast<ResultFile_XML &>(*md.res_file);
	std::fstream &file = rf.get_file();

	char str_buffer[512];

	// mesh
	const char *mesh_info = ""
		"<BackGroundMesh type = \"S2D\">\n"
		"    <h> %16.10e </h>\n"
		"    <x0> %16.10e </x0>\n"
		"    <xn> %16.10e </xn>\n"
		"    <y0> %16.10e </y0>\n"
		"    <yn> %16.10e </yn>\n"
		"    <elem_x_num> %zu </elem_x_num>\n"
		"    <elem_y_num> %zu </elem_y_num>\n"
		"</BackGroundMesh>\n";
	snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]), mesh_info, model.h,
		model.x0, model.xn, model.y0, model.yn,
		model.elem_x_num, model.elem_y_num);
	file.write(str_buffer, strlen(str_buffer));
	
	// material point object
	const char *mp_obj_info = ""
		"<MaterialPointObject type = \"ME_2D\">\n"
		"    <pcl_num> %zu </pcl_num>\n"
		"</MaterialPointObject>\n";
	snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]), mp_obj_info, model.pcl_num);
	file.write(str_buffer, strlen(str_buffer));

	return 0;
}
