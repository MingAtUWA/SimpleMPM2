#include "SimulationCore_pcp.h"

#include "ResultFile_hdf5_DataStruct.h"

#include "Model_T2D_fluid_hdf5_io_utilities.h"

namespace Model_T2D_fluid_hdf5_io_utilities
{
using namespace ResultFile_hdf5_DataStruct;

int output_model_data_to_hdf5_file(Model_T2D_fluid &md, ResultFile_hdf5 &rf)
{
	hid_t md_id = rf.get_model_data_grp_id();
	hid_t bg_mesh_id = rf.create_group(md_id, "BackgroundMesh");

	// bg mesh attributes
	const char *bg_mesh_type = "T2D";
	rf.write_attribute(bg_mesh_id, "type", strlen(bg_mesh_type), bg_mesh_type);
	rf.write_attribute(bg_mesh_id, "node_num", md.node_num);
	rf.write_attribute(bg_mesh_id, "element_num", md.elem_num);

	// node coordinates
	NodeData *nodes_data = new NodeData[md.node_num];
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		NodeData &node_data = nodes_data[n_id];
		Model_T2D_fluid::Node n = md.nodes[n_id];
		node_data.id = n.id;
		node_data.x = n.x;
		node_data.y = n.y;
	}
	hid_t nd_dt_id = get_nd_dt_id();
	rf.write_dataset(
		bg_mesh_id,
		"NodeCoordinate",
		md.node_num,
		nodes_data,
		nd_dt_id
		);
	H5Tclose(nd_dt_id);
	delete[] nodes_data;

	// element indices
	ElemData *elems_data = new ElemData[md.elem_num];
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		ElemData &elem_data = elems_data[e_id];
		Model_T2D_fluid::Element &e = md.elems[e_id];
		elem_data.id = e.id;
		elem_data.n1 = e.n1;
		elem_data.n2 = e.n2;
		elem_data.n3 = e.n3;
	}
	hid_t ed_dt_id = get_ed_dt_id();
	rf.write_dataset(
		bg_mesh_id,
		"ElementTopology",
		md.elem_num,
		elems_data,
		ed_dt_id
		);
	H5Tclose(ed_dt_id);
	delete[] elems_data;

	rf.close_group(bg_mesh_id);
	return 0;
}

int load_model_data_from_hdf5_file(Model_T2D_fluid &md, ResultFile_hdf5 &rf)
{
	// model data output
	size_t elem_num, node_num;
	hid_t md_out_id = rf.get_model_data_grp_id();

	TriangleMesh tri_mesh;
	hid_t bg_mesh_id = rf.open_group(md_out_id, "BackgroundMesh");
	rf.read_attribute(bg_mesh_id, "node_num", node_num);
	rf.read_attribute(bg_mesh_id, "element_num", elem_num);

	// nodes
	NodeData *nodes_data = new NodeData[node_num];
	hid_t nd_dt_id = get_nd_dt_id();
	rf.read_dataset(
		bg_mesh_id,
		"NodeCoordinate",
		node_num,
		nodes_data,
		nd_dt_id
		);
	H5Tclose(nd_dt_id);
	tri_mesh.alloc_nodes(node_num);
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		NodeData &node_data = nodes_data[n_id];
		TriangleMesh::Node &n = tri_mesh.get_nodes()[n_id];
		n.id = node_data.id;
		n.x = node_data.x;
		n.y = node_data.y;
	}
	delete[] nodes_data;

	// elements
	ElemData *elems_data = new ElemData[elem_num];
	hid_t ed_dt_id = get_ed_dt_id();
	rf.read_dataset(
		bg_mesh_id,
		"ElementTopology",
		elem_num,
		elems_data,
		ed_dt_id
		);
	H5Tclose(ed_dt_id);
	tri_mesh.alloc_elements(elem_num);
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		ElemData &elem_data = elems_data[e_id];
		TriangleMesh::Element &elem = tri_mesh.get_elems()[e_id];
		elem.id = elem_data.id;
		elem.n1 = elem_data.n1;
		elem.n2 = elem_data.n2;
		elem.n3 = elem_data.n3;
	}
	delete[] elems_data;

	md.init_mesh(tri_mesh);
	rf.close_group(bg_mesh_id);
	return 0;
}

int output_pcl_data_to_hdf5_file(Model_T2D_fluid &md, ResultFile_hdf5 &rf, hid_t frame_id)
{
	ParticleData *pcl_data = new ParticleData[md.pcl_num];
	for (size_t p_id = 0; p_id < md.pcl_num; ++p_id)
	{
		Model_T2D_fluid::Particle &pcl = md.pcls[p_id];
		ParticleData &pd = pcl_data[p_id];
		pd.id = pcl.id;
		pd.from_pcl(pcl);
	}
	hid_t pcl_dt_id = get_pcl_dt_id();
	int res = rf.write_dataset(
		frame_id,
		"ParticleData",
		md.pcl_num,
		pcl_data,
		pcl_dt_id
	);
	H5Tclose(pcl_dt_id);
	delete[] pcl_data;
	// particle num
	hid_t pcl_data_id = rf.open_dataset(frame_id, "ParticleData");
	rf.write_attribute(pcl_data_id, "pcl_num", md.pcl_num);
	rf.close_dataset(pcl_data_id);
	return res;
}

int load_pcl_data_from_hdf5_file(Model_T2D_fluid &md, ResultFile_hdf5 &rf, hid_t frame_id)
{
	// particle num
	size_t pcl_num;
	hid_t pcl_data_id = rf.open_dataset(frame_id, "ParticleData");
	rf.read_attribute(pcl_data_id, "pcl_num", pcl_num);
	rf.close_dataset(pcl_data_id);

	ParticleData *pcls_data = new ParticleData[pcl_num];
	hid_t pcl_dt_id = get_pcl_dt_id();
	int res = rf.read_dataset(
		frame_id,
		"ParticleData",
		pcl_num,
		pcls_data,
		pcl_dt_id
	);
	H5Tclose(pcl_dt_id);
	if (res) return res;
	md.alloc_pcls(pcl_num);
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		ParticleData &pcl_data = pcls_data[p_id];
		Model_T2D_fluid::Particle &pcl = md.pcls[p_id];
		pcl.id = pcl_data.id;
		pcl_data.to_pcl(pcl);
	}
	delete[] pcls_data;
	return 0;
}


//int load_fluid_model_from_hdf5_file(Model_T2D_fluid &md,
//	const char *hdf5_name, const char *th_name, size_t frame_id)
//{
//	ResultFile_hdf5 rf;
//	rf.open(hdf5_name);
//	hid_t file_id = rf.get_file_id();
//	if (file_id < 0) return -1;
//
//	// model data
//	load_model_data_from_hdf5_file(md, rf);
//
//	// time history output
//	hid_t th_grp_id = rf.get_time_history_grp_id();
//	hid_t th_id = rf.open_group(th_grp_id, th_name);
//	char th_frame_name[30];
//	snprintf(th_frame_name, 30, "frame_%zu", frame_id);
//	hid_t th_frame_id = rf.open_group(th_id, th_frame_name);
//	// particle data
//	load_pcl_data_from_hdf5_file(md, rf, th_frame_id);
//	// constitutive model data
//	load_model_container_from_hdf5_file(md, rf, th_frame_id);
//	// rigid body
//	load_rigid_ciricle_to_hdf5_file(md.get_rigid_circle(), rf, th_frame_id);
//	rf.close_group(th_frame_id);
//	rf.close_group(th_id);
//
//	return 0;
//}

};
