#include "SimulationCore_pcp.h"

#include "Model_S2D_ME_s.h"
#include "Model_S2D_ME_s_hdf5_io_utilities.h"

#include "ModelDataOutput_S2D_ME_s.h"

int model_data_output_func_s2d_me_s_geostatic_to_hdf5_res_file(ModelDataOutput &_self)
{
	ModelDataOutput_S2D_ME_s &md = static_cast<ModelDataOutput_S2D_ME_s &>(_self);
	Model_S2D_ME_s &model = static_cast<Model_S2D_ME_s &>(*md.model);
	ResultFile_hdf5 &rf = static_cast<ResultFile_hdf5 &>(*md.res_file);

	using Model_S2D_ME_s_hdf5_io_utilities::output_model_data_to_hdf5_file;
	output_model_data_to_hdf5_file(model, rf);

	return 0;
}

int model_data_output_func_s2d_me_s_geostatic_to_xml_res_file(ModelDataOutput &_self)
{
	ModelDataOutput_S2D_ME_s &md
		= static_cast<ModelDataOutput_S2D_ME_s &>(_self);
	Model_S2D_ME_s &model = static_cast<Model_S2D_ME_s &>(*md.model);
	ResultFile_XML &rf = static_cast<ResultFile_XML &>(*md.res_file);
	std::fstream &file = rf.get_file();

	char str_buffer[512];

	// model data
	const char *model_data_info = ""
		"<ModelData>\n"
		"    <current_time> %16.10e </current_time>\n";
		"    <total_time> %16.10e </total_time>\n";
	snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]), model_data_info,
			 md.current_time, md.total_time);
	file.write(str_buffer, strlen(str_buffer));

	// mesh
	const char *mesh_info = ""
		"    <BackGroundMesh type = \"S2D\">\n"
		"        <x0> %16.10e </x0>\n"
		"        <y0> %16.10e </y0>\n"
		"        <elem_x_num> %zu </elem_x_num>\n"
		"        <elem_y_num> %zu </elem_y_num>\n"
		"        <hx> %16.10e </hx>\n"
		"        <hy> %16.10e </hy>\n"
		"    </BackGroundMesh>\n";
	snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]), mesh_info,
			 model.get_x0(), model.get_y0(), 
			 model.get_elem_x_num(), model.get_elem_y_num(),
			 model.get_hx(), model.get_hy());
	file.write(str_buffer, strlen(str_buffer));

	// material point object
	const char *mp_obj_info = ""
		"    <MaterialPointObject type = \"ME_2D\">\n"
		"        <pcl_num> %zu </pcl_num>\n"
		"    </MaterialPointObject>\n";
	snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]), mp_obj_info, model.get_pcl_num());
	file.write(str_buffer, strlen(str_buffer));

	// ending
	file << "</ModelData>\n";

	return 0;
}
