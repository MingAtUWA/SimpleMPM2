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
	typedef ResultFile_PlainBin_DataStruct::RigidBodyHeader RigidBodyHeader;
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
	
	// rigid body
	TriangleMesh &rb_mesh = model.rigid_body.mesh;
	RigidBodyHeader rbh;
	rbh.node_num = rb_mesh.get_node_num();
	rbh.elem_num = rb_mesh.get_elem_num();
	rbh.x_mc = rb_mesh.get_x_mc();
	rbh.y_mc = rb_mesh.get_y_mc();
	file.write(reinterpret_cast<char *>(&rbh), sizeof(rbh));
	// node coordinates
	TriangleMesh::Node *rb_mesh_nodes = rb_mesh.get_nodes();
	double *rb_node_coords = new double[rbh.node_num * 2];
	for (size_t n_id = 0; n_id < rbh.node_num; ++n_id)
	{
		rb_node_coords[n_id * 2] = rb_mesh_nodes[n_id].x;
		rb_node_coords[n_id * 2 + 1] = rb_mesh_nodes[n_id].y;
	}
	file.write(reinterpret_cast<char *>(rb_node_coords), rbh.node_num * 2 * sizeof(double));
	delete[] rb_node_coords;
	// element topology
	TriangleMesh::Element *rb_mesh_elems = rb_mesh.get_elems();
	unsigned long long *rb_elem_indices = new unsigned long long[rbh.elem_num * 3];
	for (size_t e_id = 0; e_id < rbh.elem_num; ++e_id)
	{
		rb_elem_indices[e_id * 3] = rb_mesh_elems[e_id].n1;
		rb_elem_indices[e_id * 3 + 1] = rb_mesh_elems[e_id].n2;
		rb_elem_indices[e_id * 3 + 2] = rb_mesh_elems[e_id].n3;
	}
	file.write(reinterpret_cast<char *>(rb_elem_indices), rbh.elem_num * 3 * sizeof(unsigned long long));
	delete[] rb_elem_indices;

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
	
	// rigid body
	TriangleMesh &rb_mesh = model.rigid_body.mesh;
	size_t node_num = rb_mesh.get_node_num();
	size_t elem_num = rb_mesh.get_elem_num();
	const char *rigid_obj_info = ""
		"<RigidObject type = \"TriangleMesh\">\n"
		"    <node_num> %zu </node_num>\n"
		"    <elem_num> %zu </elem_num>\n"
		"    <x_mc> %16.10e </x_mc>\n"
		"    <y_mc> %16.10e </y_mc>\n";
	snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]), rigid_obj_info,
		node_num, elem_num, rb_mesh.get_x_mc(), rb_mesh.get_y_mc());
	file.write(str_buffer, strlen(str_buffer));
	// node coordinates
	file << "    <nodes>\n"
		"    <!-- index, x, y -->\n";
	const TriangleMesh::Node *rb_mesh_nodes = rb_mesh.get_nodes();
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]),
			"        %zu, %16.10e, %16.10e\n", n_id,
			rb_mesh_nodes[n_id].x, rb_mesh_nodes[n_id].y);
		file.write(str_buffer, strlen(str_buffer));
	}
	file << "    </nodes>\n";
	// element topology
	file << "    <elements>\n"
		"    <!-- index, node1, node2, node3 -->\n";
	const TriangleMesh::Element *rb_mesh_elems = rb_mesh.get_elems();
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		file << "        " << e_id
			<< ", " << rb_mesh_elems[e_id].n1
			<< ", " << rb_mesh_elems[e_id].n2
			<< ", " << rb_mesh_elems[e_id].n3 << "\n";
	}
	file << "    </elements>\n";
	// ending
	file << "</RigidObject>\n";

	// material point object
	const char *mp_obj_info = ""
		"<MaterialPointObject type = \"ME_2D\">\n"
		"    <pcl_num> %zu </pcl_num>\n"
		"</MaterialPointObject>\n";
	snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]), mp_obj_info, model.pcl_num);
	file.write(str_buffer, strlen(str_buffer));

	return 0;
}
