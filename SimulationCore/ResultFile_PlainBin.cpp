#include "SimulationCore_pcp.h"

#include "ResultFileData.h"

#include "Step_S2D_ME_s_RigidBody.h"
#include "TimeHistory_S2D_ME_s_RigidBody.h"

#include "Step_S2D_ME_s_RigidBody_Fric.h"
#include "TimeHistory_S2D_ME_s_RigidBody_Fric.h"

#include "ResultFile_PlainBin.h"

ResultFile_PlainBin::ResultFile_PlainBin() {}

ResultFile_PlainBin::~ResultFile_PlainBin() { finalize(); }

int ResultFile_PlainBin::init(const char *file_name)
{
	file.open(file_name, std::ios::binary | std::ios::out);
	if (file.is_open())
		return 0;
	return -1;
}

void ResultFile_PlainBin::finalize(void)
{
	file.close();
}

// MPM - Rigid body contact problem
// Output model data
int ResultFile_PlainBin::output(Model_S2D_ME_s_RigidBody &model)
{
	typedef ResultFileData::ME_s_RigidBody::MeshHeader MeshHeader;
	typedef ResultFileData::ME_s_RigidBody::RigidBodyHeader RigidBodyHeader;
	typedef ResultFileData::ME_s_RigidBody::MPObjectHeader MPObjectHeader;
	
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
	RigidBodyHeader rbh;
	TriangleMesh &rb_mesh = model.rigid_body.mesh;
	rbh.node_num = rb_mesh.get_node_num();
	rbh.elem_num = rb_mesh.get_elem_num();
	rbh.x_mc = rb_mesh.get_x_mc();
	rbh.y_mc = rb_mesh.get_y_mc();
	file.write(reinterpret_cast<char *>(&rbh), sizeof(rbh));
	// node coordinates
	double *node_coords = new double[rbh.node_num * 2];
	const TriangleMesh::Node *rb_mesh_nodes = rb_mesh.get_nodes();
	for (size_t n_id = 0; n_id < rbh.node_num; ++n_id)
	{
		node_coords[n_id * 2] = rb_mesh_nodes[n_id].x;
		node_coords[n_id * 2 + 1] = rb_mesh_nodes[n_id].y;
	}
	file.write(reinterpret_cast<char *>(node_coords), rbh.node_num * 2 * sizeof(double));
	delete[] node_coords;
	// element topology
	unsigned long long *elem_indices = new unsigned long long[rbh.elem_num * 3];
	const TriangleMesh::Element *rb_mesh_elems = rb_mesh.get_elems();
	for (size_t e_id = 0; e_id < rbh.elem_num; ++e_id)
	{
		elem_indices[e_id * 3] = rb_mesh_elems[e_id].n1;
		elem_indices[e_id * 3 + 1] = rb_mesh_elems[e_id].n2;
		elem_indices[e_id * 3 + 2] = rb_mesh_elems[e_id].n3;
	}
	file.write(reinterpret_cast<char *>(elem_indices), rbh.elem_num * 3 * sizeof(unsigned long long));
	delete[] elem_indices;

	// material point object
	MPObjectHeader mph;
	mph.pcl_num = model.pcl_num;
	file.write(reinterpret_cast<char *>(&mph), sizeof(mph));

	return 0;
}

// Output time history
int rf_out_func_imp_th_MPM_RigidBody_PlainBin(TimeHistory &_th, ResultFile &_rf)
{
	TimeHistory_S2D_ME_s_RigidBody &th = static_cast<TimeHistory_S2D_ME_s_RigidBody &>(_th);
	ResultFile_PlainBin &rf = static_cast<ResultFile_PlainBin &>(_rf);
	std::fstream &file = rf.file;

	typedef ResultFileData::TimeHistoryHeader TimeHistoryHeader;
	typedef ResultFileData::ME_s_RigidBody::RigidBodyMotionHeader RigidBodyMotionHeader;
	typedef ResultFileData::ME_s_RigidBody::MPObjectHeader MPObjectHeader;
	
	TimeHistoryHeader thh;
	Step_S2D_ME_s_RigidBody &step = static_cast<Step_S2D_ME_s_RigidBody &>(th.get_step());
	thh.substep_num = step.get_substep_num();
	thh.total_substep_num = step.get_total_substep_num();
	thh.current_time = step.get_current_time();
	thh.total_time = step.get_total_time();
	file.write(reinterpret_cast<char *>(&thh), sizeof(thh));

	RigidBodyMotionHeader rbmh;
	Model_S2D_ME_s_RigidBody &model = static_cast<Model_S2D_ME_s_RigidBody &>(th.get_model());
	RigidBody &rb = model.rigid_body;
	rbmh.x = rb.x;
	rbmh.y = rb.y;
	rbmh.theta = rb.theta;
	rbmh.vx = rb.vx;
	rbmh.vy = rb.vy;
	rbmh.v_theta = rb.vtheta;
	rbmh.fx_con = rb.Fx_con;
	rbmh.fy_con = rb.Fy_con;
	rbmh.m_con = rb.M_con;
	file.write(reinterpret_cast<char *>(&rbmh), sizeof(rbmh));

	// output particles data
	MPObjectHeader mph;
	mph.pcl_num = model.pcl_num;
	file.write(reinterpret_cast<char *>(&mph), sizeof(mph));
	double *pcl_data = new double[mph.pcl_num];
	size_t data_len = sizeof(double) * mph.pcl_num;
	// x
	for (size_t pcl_id = 0; pcl_id < mph.pcl_num; ++pcl_id)
		pcl_data[pcl_id] = model.pcls[pcl_id].x;
	file.write(reinterpret_cast<char *>(pcl_data), data_len);
	// y
	for (size_t pcl_id = 0; pcl_id < mph.pcl_num; ++pcl_id)
		pcl_data[pcl_id] = model.pcls[pcl_id].y;
	file.write(reinterpret_cast<char *>(pcl_data), data_len);
	// vol
	for (size_t pcl_id = 0; pcl_id < mph.pcl_num; ++pcl_id)
		pcl_data[pcl_id] = model.pcls[pcl_id].m / model.pcls[pcl_id].density;
	file.write(reinterpret_cast<char *>(pcl_data), data_len);
	delete[] pcl_data;

	return 0;
}

const OutputFunc ResultFile_PlainBin::out_func_th_MPM_RigidBody_PlainBin = &rf_out_func_imp_th_MPM_RigidBody_PlainBin;

// MPM - Rigid body frictional contact problem
// Output model data
int ResultFile_PlainBin::output(Model_S2D_ME_s_RigidBody_Fric &model)
{
	typedef ResultFileData::ME_s_RigidBody::MeshHeader MeshHeader;
	typedef ResultFileData::ME_s_RigidBody::RigidBodyHeader RigidBodyHeader;
	typedef ResultFileData::ME_s_RigidBody::MPObjectHeader MPObjectHeader;

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
	RigidBodyHeader rbh;
	TriangleMesh &rb_mesh = model.rigid_body.mesh;
	rbh.node_num = rb_mesh.get_node_num();
	rbh.elem_num = rb_mesh.get_elem_num();
	rbh.x_mc = rb_mesh.get_x_mc();
	rbh.y_mc = rb_mesh.get_y_mc();
	file.write(reinterpret_cast<char *>(&rbh), sizeof(rbh));
	// node coordinates
	double *node_coords = new double[rbh.node_num * 2];
	const TriangleMesh::Node *rb_mesh_nodes = rb_mesh.get_nodes();
	for (size_t n_id = 0; n_id < rbh.node_num; ++n_id)
	{
		node_coords[n_id * 2] = rb_mesh_nodes[n_id].x;
		node_coords[n_id * 2 + 1] = rb_mesh_nodes[n_id].y;
	}
	file.write(reinterpret_cast<char *>(node_coords), rbh.node_num * 2 * sizeof(double));
	delete[] node_coords;
	// element topology
	unsigned long long *elem_indices = new unsigned long long[rbh.elem_num * 3];
	const TriangleMesh::Element *rb_mesh_elems = rb_mesh.get_elems();
	for (size_t e_id = 0; e_id < rbh.elem_num; ++e_id)
	{
		elem_indices[e_id * 3] = rb_mesh_elems[e_id].n1;
		elem_indices[e_id * 3 + 1] = rb_mesh_elems[e_id].n2;
		elem_indices[e_id * 3 + 2] = rb_mesh_elems[e_id].n3;
	}
	file.write(reinterpret_cast<char *>(elem_indices), rbh.elem_num * 3 * sizeof(unsigned long long));
	delete[] elem_indices;

	// material point object
	MPObjectHeader mph;
	mph.pcl_num = model.pcl_num;
	file.write(reinterpret_cast<char *>(&mph), sizeof(mph));

	return 0;
}

// Output time history
int rf_out_func_imp_th_MPM_RigidBody_Fric_PlainBin(TimeHistory &_th, ResultFile &_rf)
{
	TimeHistory_S2D_ME_s_RigidBody_Fric &th = static_cast<TimeHistory_S2D_ME_s_RigidBody_Fric &>(_th);
	ResultFile_PlainBin &rf = static_cast<ResultFile_PlainBin &>(_rf);
	std::fstream &file = rf.file;

	typedef ResultFileData::TimeHistoryHeader TimeHistoryHeader;
	typedef ResultFileData::ME_s_RigidBody::RigidBodyMotionHeader RigidBodyMotionHeader;
	typedef ResultFileData::ME_s_RigidBody::MPObjectHeader MPObjectHeader;

	TimeHistoryHeader thh;
	Step_S2D_ME_s_RigidBody_Fric &step = static_cast<Step_S2D_ME_s_RigidBody_Fric &>(th.get_step());
	thh.substep_num = step.get_substep_num();
	thh.total_substep_num = step.get_total_substep_num();
	thh.current_time = step.get_current_time();
	thh.total_time = step.get_total_time();
	file.write(reinterpret_cast<char *>(&thh), sizeof(thh));

	RigidBodyMotionHeader rbmh;
	Model_S2D_ME_s_RigidBody_Fric &model = static_cast<Model_S2D_ME_s_RigidBody_Fric &>(th.get_model());
	RigidBody &rb = model.rigid_body;
	rbmh.x = rb.x;
	rbmh.y = rb.y;
	rbmh.theta = rb.theta;
	rbmh.vx = rb.vx;
	rbmh.vy = rb.vy;
	rbmh.v_theta = rb.vtheta;
	rbmh.fx_con = rb.Fx_con;
	rbmh.fy_con = rb.Fy_con;
	rbmh.m_con = rb.M_con;
	file.write(reinterpret_cast<char *>(&rbmh), sizeof(rbmh));

	// output particles data
	MPObjectHeader mph;
	mph.pcl_num = model.pcl_num;
	file.write(reinterpret_cast<char *>(&mph), sizeof(mph));
	double *pcl_data = new double[mph.pcl_num];
	size_t data_len = sizeof(double) * mph.pcl_num;
	// x
	for (size_t pcl_id = 0; pcl_id < mph.pcl_num; ++pcl_id)
		pcl_data[pcl_id] = model.pcls[pcl_id].x;
	file.write(reinterpret_cast<char *>(pcl_data), data_len);
	// y
	for (size_t pcl_id = 0; pcl_id < mph.pcl_num; ++pcl_id)
		pcl_data[pcl_id] = model.pcls[pcl_id].y;
	file.write(reinterpret_cast<char *>(pcl_data), data_len);
	// vol
	for (size_t pcl_id = 0; pcl_id < mph.pcl_num; ++pcl_id)
		pcl_data[pcl_id] = model.pcls[pcl_id].m / model.pcls[pcl_id].density;
	file.write(reinterpret_cast<char *>(pcl_data), data_len);
	delete[] pcl_data;

	return 0;
}

const OutputFunc ResultFile_PlainBin::out_func_th_MPM_RigidBody_Fric_PlainBin = &rf_out_func_imp_th_MPM_RigidBody_Fric_PlainBin;
