#include "SimulationCore_pcp.h"

#include "Model_T2D_CHM_s.h"
#include "Step_T2D_CHM_s.h"

#include "ResultFile_PlainBin_DataStruct.h"
#include "Model_T2D_CHM_s_hdf5_io_utilities.h"

#include "ModelDataOutput_T2D_CHM_s.h"

int model_data_output_func_t2d_chm_s_to_plain_bin_res_file(ModelDataOutput &_self)
{
	ModelDataOutput_T2D_CHM_s &md = static_cast<ModelDataOutput_T2D_CHM_s &>(_self);
	Model_T2D_CHM_s &model = static_cast<Model_T2D_CHM_s &>(*md.model);
	ResultFile_PlainBin &rf = static_cast<ResultFile_PlainBin &>(*md.res_file);
	std::fstream &file = rf.get_file();
	
	typedef ResultFile_PlainBin_DataStruct::ModelDataHeader ModelDataHeader;
	typedef ResultFile_PlainBin_DataStruct::BackgroundMeshHeader_T2D BackgroundMeshHeader;
	typedef ResultFile_PlainBin_DataStruct::DispConRigidCircleHeader DispConRigidCircleHeader;
	typedef ResultFile_PlainBin_DataStruct::MPObjectHeader MPObjectHeader;
	
	// model data
	ModelDataHeader mdh;
	mdh.init();
	mdh.current_time = md.current_time;
	mdh.total_time = md.total_time;
	file.write(reinterpret_cast<char *>(&mdh), sizeof(mdh));

	// bgmesh
	BackgroundMeshHeader mh;
	mh.init();
	mh.node_num = model.node_num;
	mh.elem_num = model.elem_num;
	file.write(reinterpret_cast<char *>(&mh), sizeof(mh));
	// nodal coordinates
	double *n_coords = new double[model.node_num * 2];
	for (size_t n_id = 0; n_id < model.node_num; ++n_id)
	{
		Model_T2D_CHM_s::Node &n = model.nodes[n_id];
		n_coords[n_id*2]     = n.x;
		n_coords[n_id*2 + 1] = n.y;
	}
	file.write(reinterpret_cast<char *>(n_coords), model.node_num * 2 * sizeof(double));
	delete[] n_coords;
	// element node indices
	unsigned long long *e_n_ids = new unsigned long long[model.elem_num * 3];
	for (size_t e_id = 0; e_id < model.elem_num; ++e_id)
	{
		Model_T2D_CHM_s::Element &e = model.elems[e_id];
		e_n_ids[3*e_id]     = e.n1;
		e_n_ids[3*e_id + 1] = e.n2;
		e_n_ids[3*e_id + 2] = e.n3;
	}
	file.write(reinterpret_cast<char *>(e_n_ids), model.elem_num * 3 * sizeof(unsigned long long));
	delete[] e_n_ids;

	// rigid circle object
	DispConRigidCircle &rc = model.get_rigid_circle();
	DispConRigidCircleHeader rch;
	rch.init();
	rch.pcl_num = rc.get_pcl_num();
	file.write(reinterpret_cast<char *>(&rch), sizeof(rch));

	// material point object
	MPObjectHeader mph;
	mph.init();
	mph.pcl_num = model.pcl_num;
	mph.fld_num = 0;
	file.write(reinterpret_cast<char *>(&mph), sizeof(mph));

	return 0;
}

int model_data_output_func_t2d_chm_s_to_xml_res_file(ModelDataOutput &_self)
{
	ModelDataOutput_T2D_CHM_s &md = static_cast<ModelDataOutput_T2D_CHM_s &>(_self);
	Model_T2D_CHM_s &model = static_cast<Model_T2D_CHM_s &>(*md.model);
	ResultFile_PlainBin &rf = static_cast<ResultFile_PlainBin &>(*md.res_file);
	std::fstream &file = rf.get_file();

	char str_buffer[512];

	// model data
	const char *model_data_info = ""
		"<ModelData>\n"
		"    <current_time> %16.10e </current_time>\n"
		"    <total_time> %16.10e </total_time>\n";
	snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]), model_data_info,
			 md.current_time, md.total_time);
	file.write(str_buffer, strlen(str_buffer));

	// mesh
	const char *mesh_info1 = ""
		"    <BackGroundMesh type = \"T2D\">\n"
		"        <Nodes>\n"
		"            <num> %zu </num>\n"
		"            <coordinates>\n";
	snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]), mesh_info1, model.node_num);
	file.write(str_buffer, strlen(str_buffer));
	const char *mesh_coords = "            %le, %le\n";
	for (size_t n_id = 0; n_id < model.node_num; ++n_id)
	{
		Model_T2D_CHM_s::Node &n = model.nodes[n_id];
		snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]), mesh_coords, n.x, n.y);
		file.write(str_buffer, strlen(str_buffer));
	}
	const char *mesh_info2 = ""
		"            </coordinates>\n"
		"        </Nodes>\n"
		"        <Elements>\n"
		"            <num> %zu </num>\n"
		"            <node_indices>\n";
	snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]), mesh_info2, model.elem_num);
	file.write(str_buffer, strlen(str_buffer));
	const char *mesh_indices = "            %zu, %zu, %zu\n";
	for (size_t e_id = 0; e_id < model.elem_num; ++e_id)
	{
		Model_T2D_CHM_s::Element &e = model.elems[e_id];
		snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]), mesh_indices, e.n1, e.n2, e.n3);
		file.write(str_buffer, strlen(str_buffer));
	}
	const char *mesh_info3 = ""
		"            </node_indices>\n"
		"        </Elements>\n"
		"    </BackGroundMesh>\n";
	file.write(mesh_info3, strlen(mesh_info3));

	// material point object
	const char *mp_obj_info = ""
		"    <MaterialPointObject type = \"ME_2D\">\n"
		"        <pcl_num> %zu </pcl_num>\n"
		"    </MaterialPointObject>\n";
	snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]), mp_obj_info, model.pcl_num);
	file.write(str_buffer, strlen(str_buffer));

	// ending
	file << "</ModelData>\n";
	
	return 0;
}


int model_data_output_func_t2d_chm_s_to_hdf5_res_file(ModelDataOutput &_self)
{
	ModelDataOutput_T2D_CHM_s &md = static_cast<ModelDataOutput_T2D_CHM_s &>(_self);
	Model_T2D_CHM_s &model = static_cast<Model_T2D_CHM_s &>(*md.model);
	ResultFile_hdf5 &rf = static_cast<ResultFile_hdf5 &>(*md.res_file);
	
	using Model_T2D_CHM_s_hdf5_io_utilities::output_model_data_to_hdf5_file;
	output_model_data_to_hdf5_file(model, rf);
	return 0;
}