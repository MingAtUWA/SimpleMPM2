#include "SimulationCore_pcp.h"

#include "Model_T2D_CHM_s_hdf5_io_utilities.h"

namespace
{
struct NodeData
{
	unsigned long long id;
	double x;
	double y;
};
struct ElemData
{
	unsigned long long id;
	unsigned long long n1;
	unsigned long long n2;
	unsigned long long n3;
};
};

int output_model_data_to_hdf5_file(Model_T2D_CHM_s &md, ResultFile_hdf5 &rf)
{
	hid_t md_id = rf.get_model_data_grp_id();
	rf.write_attribute(md_id, "Kf", md.Kf);
	rf.write_attribute(md_id, "k", md.k);
	rf.write_attribute(md_id, "miu", md.miu);

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
		Model_T2D_CHM_s::Node n = md.nodes[n_id];
		node_data.id = n.id;
		node_data.x = n.x;
		node_data.y = n.y;
	}
	rf.write_dataset(bg_mesh_id, "NodeCoordinate", md.node_num, 
					 size_t(sizeof(NodeData)/sizeof(double)), 
					 (double *)nodes_data);
	delete[] nodes_data;
	// element indices
	ElemData *elems_data = new ElemData[md.elem_num];
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		ElemData &elem_data = elems_data[e_id];
		Model_T2D_CHM_s::Element &e = md.elems[e_id];
		elem_data.id = e.id;
		elem_data.n1 = e.n1;
		elem_data.n2 = e.n2;
		elem_data.n3 = e.n3;
	}
	rf.write_dataset(bg_mesh_id, "ElementTopology", md.elem_num, 
					 size_t(sizeof(ElemData)/sizeof(unsigned long long)),
					 (unsigned long long *)elems_data);
	delete[] elems_data;
	rf.close_group(bg_mesh_id);

	return 0;
}


int load_model_data_from_hdf5_file(Model_T2D_CHM_s &md, ResultFile_hdf5 &rf)
{
	// model data output
	size_t elem_num, node_num;
	hid_t md_out_id = rf.get_model_data_grp_id();
	
	rf.read_attribute(md_out_id, "Kf", md.Kf);
	rf.read_attribute(md_out_id, "k", md.k);
	rf.read_attribute(md_out_id, "miu", md.miu);

	TriangleMesh tri_mesh;
	hid_t bg_mesh_id = rf.open_group(md_out_id, "BackgroundMesh");
	rf.read_attribute(bg_mesh_id, "node_num", node_num);
	rf.read_attribute(bg_mesh_id, "element_num", elem_num);
	
	// nodes
	NodeData *node_data = new NodeData[node_num];
	rf.read_dataset(bg_mesh_id, "NodeCoordinate", node_num,
					size_t(sizeof(NodeData) / sizeof(double)),
					(double *)node_data);
	tri_mesh.alloc_nodes(node_num);
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		TriangleMesh::Node &n = tri_mesh.get_nodes()[n_id];
		n.id = node_data[n_id].id;
		n.x = node_data[n_id].x;
		n.y = node_data[n_id].y;
	}
	delete[] node_data;

	// elements
	ElemData *elem_data = new ElemData[elem_num];
	rf.read_dataset(bg_mesh_id, "ElementTopology", elem_num,
		size_t(sizeof(ElemData) / sizeof(unsigned long long)),
		(unsigned long long *)elem_data);
	tri_mesh.alloc_elements(elem_num);
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		TriangleMesh::Element &elem = tri_mesh.get_elems()[e_id];
		elem.id = elem_data[e_id].id;
		elem.n1 = elem_data[e_id].n1;
		elem.n2 = elem_data[e_id].n2;
		elem.n3 = elem_data[e_id].n3;
	}
	delete[] elem_data;
	md.init_mesh(tri_mesh);
	// close
	rf.close_group(bg_mesh_id);

	return 0;
}


namespace
{
// Linear Elasticity
struct LinearElasticityStateData
{
	unsigned long long pcl_id;
	double E; // Young's modulus
	double niu; // possion ratio
	void from_cm(LinearElasticity &cm)
	{
		E = cm.E;
		niu = cm.niu;
	}
	void to_cm(LinearElasticity &cm)
	{
		cm.set_param(E, niu);
	}
};

hid_t get_le_hdf5_dt_id(void)
{
	hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(LinearElasticityStateData));
	H5Tinsert(res, "pcl_id", HOFFSET(LinearElasticityStateData, pcl_id), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "E", HOFFSET(LinearElasticityStateData, E), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "niu", HOFFSET(LinearElasticityStateData, niu), H5T_NATIVE_DOUBLE);
	return res;
}

// ModifiedCamClay
struct ModifiedCamClayStateData
{
	unsigned long long pcl_id;
	double niu; // possion ratio
	double kappa; // logrithmic recompression modulus
	double lambda; // logrithmic compression modulus
	double fric_angle; // friction angle
	double e; // void ratio
	double pc; // pre-consolidation stress
	double N; // normal consolidation line
	double Gamma; // critical state line
	double M; // critical state line q = Mp 
	double s11, s22, s33, s12, s23, s31; // stress state
	void from_cm(ModifiedCamClay &cm)
	{
		niu = cm.niu;
		kappa = cm.kappa;
		lambda = cm.lambda;
		fric_angle = cm.fric_angle;
		e = cm.e;
		pc = cm.pc;
		N = cm.N;
		Gamma = cm.Gamma;
		M = cm.M;
		s11 = cm.s11;
		s22 = cm.s22;
		s33 = cm.s33;
		s12 = cm.s12;
		s23 = cm.s23;
		s31 = cm.s31;
	}
	void to_cm(ModifiedCamClay &cm)
	{
		cm.niu = niu;
		cm.kappa = kappa;
		cm.lambda = lambda;
		cm.fric_angle = fric_angle;
		cm.e = e;
		cm.pc = pc;
		cm.N = N;
		cm.Gamma = Gamma;
		cm.M = M;
		cm.M2 = M * M;
		cm.s11 = s11;
		cm.s22 = s22;
		cm.s33 = s33;
		cm.s12 = s12;
		cm.s23 = s23;
		cm.s31 = s31;
	}
};

hid_t get_mcc_hdf5_dt_id(void)
{
	hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(ModifiedCamClayStateData));
	H5Tinsert(res, "pcl_id", HOFFSET(ModifiedCamClayStateData, pcl_id), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "niu", HOFFSET(ModifiedCamClayStateData, niu), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "kappa", HOFFSET(ModifiedCamClayStateData, kappa), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "lambda", HOFFSET(ModifiedCamClayStateData, lambda), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "fric_angle", HOFFSET(ModifiedCamClayStateData, fric_angle), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "e", HOFFSET(ModifiedCamClayStateData, e), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "pc", HOFFSET(ModifiedCamClayStateData, pc), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "N", HOFFSET(ModifiedCamClayStateData, N), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "Gamma", HOFFSET(ModifiedCamClayStateData, Gamma), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "M", HOFFSET(ModifiedCamClayStateData, M), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "s11", HOFFSET(ModifiedCamClayStateData, s11), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "s22", HOFFSET(ModifiedCamClayStateData, s22), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "s33", HOFFSET(ModifiedCamClayStateData, s33), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "s12", HOFFSET(ModifiedCamClayStateData, s12), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "s23", HOFFSET(ModifiedCamClayStateData, s23), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "s31", HOFFSET(ModifiedCamClayStateData, s31), H5T_NATIVE_DOUBLE);
	return res;
}

template <typename Particle>
inline size_t get_pcl_id(void *ext_data) { return static_cast<Particle *>(ext_data)->id; }

};

int output_model_container_to_hdf5_file(ModelContainer &mc, ResultFile_hdf5 &rf, hid_t frame_id)
{
	size_t cm_id, cm_num;
	hid_t cm_grp_id = rf.create_group(frame_id, "ConstitutiveModel");

	// linear elasticity
	cm_num = mc.get_num_LinearElasticity();
	if (cm_num)
	{
		LinearElasticityStateData *cm_data = new LinearElasticityStateData[cm_num];
		cm_id = 0;
		for (LinearElasticity *iter = mc.first_LinearElasticity();
			 iter; iter = mc.next_LinearElasticity(iter))
		{
			cm_data[cm_id].pcl_id = get_pcl_id<Model_T2D_CHM_s::Particle>(iter->ext_data);
			cm_data[cm_id].from_cm(*iter);
			++cm_id;
		}
		hid_t le_dt_id = get_le_hdf5_dt_id();
		rf.write_dataset(cm_grp_id, "LinearElasticity", cm_num,
						 cm_data, le_dt_id);
		H5Tclose(le_dt_id);
		delete[] cm_data;
		hid_t cm_dset_id = rf.open_dataset(cm_grp_id, "LinearElasticity");
		rf.write_attribute(cm_dset_id, "cm_num", cm_num);
		rf.close_dataset(cm_dset_id);
	}

	cm_num = mc.get_num_ModifiedCamClay();
	if (cm_num)
	{
		ModifiedCamClayStateData *cm_data = new ModifiedCamClayStateData[cm_num];
		cm_id = 0;
		for (ModifiedCamClay *iter = mc.first_ModifiedCamClay();
			 iter; iter = mc.next_ModifiedCamClay(iter))
		{
			cm_data[cm_id].pcl_id
				= get_pcl_id<Model_T2D_CHM_s::Particle>(iter->ext_data);
			cm_data[cm_id].from_cm(*iter);
			++cm_id;
		}
		hid_t mcc_dt_id = get_mcc_hdf5_dt_id();
		rf.write_dataset(cm_grp_id, "ModifiedCamClay", cm_num,
						 cm_data, mcc_dt_id);
		H5Tclose(mcc_dt_id);
		delete[] cm_data;
		hid_t cm_dset_id = rf.open_dataset(cm_grp_id, "ModifiedCamClay");
		rf.write_attribute(cm_dset_id, "cm_num", cm_num);
		rf.close_dataset(cm_dset_id);
	}

	rf.close_group(cm_grp_id);
	return 0;
}

int load_model_container_from_hdf5_file(Model_T2D_CHM_s &md, ResultFile_hdf5 &rf, hid_t frame_id)
{
	hid_t cm_dset_id;
	size_t cm_num;
	hid_t cm_grp_id = rf.open_group(frame_id, "ConstitutiveModel");
	ModelContainer &mc = md.model_container;

	// linear elasticity
	if (rf.has_dataset(cm_grp_id, "LinearElasticity"))
	{
		// cm_num
		cm_dset_id = rf.open_dataset(cm_grp_id, "LinearElasticity");
		rf.read_attribute(cm_dset_id, "cm_num", cm_num);
		rf.close_dataset(cm_dset_id);
		// get data
		LinearElasticityStateData *cm_data = new LinearElasticityStateData[cm_num];
		hid_t le_dt_id = get_le_hdf5_dt_id();
		rf.read_dataset(cm_grp_id, "LinearElasticity", cm_num,
						cm_data, le_dt_id);
		H5Tclose(le_dt_id);
		LinearElasticity *cms = mc.add_LinearElasticity(cm_num);
		for (size_t cm_id = 0; cm_id < cm_num; ++cm_id)
		{
			LinearElasticityStateData &cmd = cm_data[cm_id];
			LinearElasticity &cm = cms[cm_id];
			//std::cout << cmd.pcl_id << " " << md.pcls[cmd.pcl_id].id << "\n";
			cmd.to_cm(cm);
			cm.ext_data = &md.pcls[cmd.pcl_id];
			md.pcls[cmd.pcl_id].cm = &cm;
		}
		delete[] cm_data;
	}

	// modified cam clay
	if (rf.has_dataset(cm_grp_id, "ModifiedCamClay"))
	{
		// cm_num
		cm_dset_id = rf.open_dataset(cm_grp_id, "ModifiedCamClay");
		rf.read_attribute(cm_dset_id, "cm_num", cm_num);
		rf.close_dataset(cm_dset_id);
		// get data
		ModifiedCamClayStateData *cm_data = new ModifiedCamClayStateData[cm_num];
		hid_t mcc_dt_id = get_mcc_hdf5_dt_id();
		rf.read_dataset(cm_grp_id, "ModifiedCamClay", cm_num,
						cm_data, mcc_dt_id);
		H5Tclose(mcc_dt_id);
		ModifiedCamClay *cms = mc.add_ModifiedCamClay(cm_num);
		for (size_t cm_id = 0; cm_id < cm_num; ++cm_id)
		{
			ModifiedCamClayStateData &cmd = cm_data[cm_id];
			ModifiedCamClay &cm = cms[cm_id];
			cmd.to_cm(cm);
			cm.ext_data = &md.pcls[cmd.pcl_id];
			md.pcls[cmd.pcl_id].cm = &cm;
		}
		delete[] cm_data;
	}

	rf.close_group(cm_grp_id);
	return 0;
}


namespace
{
struct RigidBodyParticleData
{
	double xr;
	double yr;
	double vol;
};
};
int output_rigid_ciricle_to_hdf5_file(DispConRigidCircle &rc, ResultFile_hdf5 &rf, hid_t frame_id)
{
	size_t pcl_num = rc.get_pcl_num();
	if (pcl_num == 0) return 0;

	DispConRigidCircle::State state = rc.get_state();
	DispConRigidCircle::Particle *rb_pcls = state.pcls;

	hid_t rb_id = rf.create_group(frame_id, "RigidBody");

	// rigid body data
	rf.write_attribute(rb_id, "radius", state.r);
	rf.write_attribute(rb_id, "cen_x", state.cen_x);
	rf.write_attribute(rb_id, "cen_y", state.cen_y);
	rf.write_attribute(rb_id, "theta", state.theta);
	rf.write_attribute(rb_id, "vx", state.vx);
	rf.write_attribute(rb_id, "vy", state.vy);
	rf.write_attribute(rb_id, "w", state.w);
	rf.write_attribute(rb_id, "rfx", state.rfx);
	rf.write_attribute(rb_id, "rfy", state.rfy);
	rf.write_attribute(rb_id, "rm", state.rm);
	rf.write_attribute(rb_id, "pcl_num", pcl_num);

	// particle data
	RigidBodyParticleData *rb_pcls_data = new RigidBodyParticleData[pcl_num];
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		RigidBodyParticleData &rb_pcl_data = rb_pcls_data[p_id];
		DispConRigidCircle::Particle &rb_pcl = rb_pcls[p_id];
		rb_pcl_data.xr = rb_pcl.xr;
		rb_pcl_data.yr = rb_pcl.yr;
		rb_pcl_data.vol = rb_pcl.vol;
	}
	rf.write_dataset(rb_id, "ParticleData", pcl_num, 
					 size_t(sizeof(RigidBodyParticleData) / sizeof(double)),
					 (double *)rb_pcls_data);

	rf.close_group(rb_id);
	return 0;
}

int load_rigid_ciricle_to_hdf5_file(DispConRigidCircle &rc, ResultFile_hdf5 &rf, hid_t frame_id)
{
	if (!rf.has_group(frame_id, "RigidBody")) return 0;

	hid_t rb_id = rf.open_group(frame_id, "RigidBody");

	// rigid body data
	DispConRigidCircle::State &rc_state 
		= const_cast<DispConRigidCircle::State &>(rc.get_state());
	rf.read_attribute(rb_id, "radius", rc_state.r);
	rc_state.r2 = rc_state.r * rc_state.r;
	double cen_x, cen_y, theta;
	size_t pcl_num;
	rf.read_attribute(rb_id, "cen_x", cen_x);
	rc_state.cen_x = cen_x;
	rf.read_attribute(rb_id, "cen_y", cen_y);
	rc_state.cen_y = cen_y;
	rf.read_attribute(rb_id, "theta", theta);
	rc_state.theta = theta;
	rf.read_attribute(rb_id, "vx", rc_state.vx);
	rf.read_attribute(rb_id, "vy", rc_state.vy);
	rf.read_attribute(rb_id, "w", rc_state.w);
	rf.read_attribute(rb_id, "rfx", rc_state.rfx);
	rf.read_attribute(rb_id, "rfy", rc_state.rfy);
	rf.read_attribute(rb_id, "rm", rc_state.rm);
	rf.read_attribute(rb_id, "pcl_num", pcl_num);

	// particle data
	struct RigidBodyParticleData *rb_pcls_data
		= new RigidBodyParticleData[pcl_num];
	rf.read_dataset(rb_id, "ParticleData", pcl_num,
		size_t(sizeof(RigidBodyParticleData) / sizeof(double)),
		(double *)rb_pcls_data);
	DispConRigidCircle::Particle *rb_pcls = rc.alloc_pcls(pcl_num);
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		RigidBodyParticleData &rb_pcl_data = rb_pcls_data[p_id];
		DispConRigidCircle::Particle &rb_pcl = rb_pcls[p_id];
		rb_pcl.xr = rb_pcl_data.xr;
		rb_pcl.yr = rb_pcl_data.yr;
		rb_pcl.vol = rb_pcl_data.vol;
		rb_pcl.x = cen_x + rb_pcl.xr * cos(theta) + rb_pcl.yr * -sin(theta);
		rb_pcl.y = cen_y + rb_pcl.xr * sin(theta) + rb_pcl.yr *  cos(theta);
	}
	delete[] rb_pcls_data;

	rf.close_group(rb_id);
	return 0;
}

namespace
{
struct ParticleData
{
	unsigned long long id;
	double m_s;
	double density_s;
	double density_f;
	double n;
	double x;
	double y;
	double vx_s;
	double vy_s;
	double vx_f;
	double vy_f;
	double s11;
	double s22;
	double s12;
	double p;
	double e11;
	double e22;
	double e12;
	void from_pcl(Model_T2D_CHM_s::Particle &pcl)
	{
		m_s = pcl.m_s;
		density_s = pcl.density_s;
		density_f = pcl.density_f;
		n = pcl.n;
		x = pcl.x;
		y = pcl.y;
		vx_s = pcl.vx_s;
		vy_s = pcl.vy_s;
		vx_f = pcl.vx_f;
		vy_f = pcl.vy_f;
		s11 = pcl.s11;
		s22 = pcl.s22;
		s12 = pcl.s12;
		p = pcl.p;
		e11 = pcl.e11;
		e22 = pcl.e22;
		e12 = pcl.e12;
	}
	void to_pcl(Model_T2D_CHM_s::Particle &pcl)
	{
		pcl.m_s = m_s;
		pcl.density_s = density_s;
		pcl.density_f = density_f;
		pcl.n = n;
		pcl.x = x;
		pcl.y = y;
		pcl.vx_s = vx_s;
		pcl.vy_s = vy_s;
		pcl.vx_f = vx_f;
		pcl.vy_f = vy_f;
		pcl.s11 = s11;
		pcl.s22 = s22;
		pcl.s12 = s12;
		pcl.p = p;
		pcl.e11 = e11;
		pcl.e22 = e22;
		pcl.e12 = e12;
	}
};
};

int ouput_pcl_data_to_hdf5_file(Model_T2D_CHM_s &md, ResultFile_hdf5 &rf, hid_t frame_id)
{
	ParticleData *pcl_data = new ParticleData[md.pcl_num];
	for (size_t p_id = 0; p_id < md.pcl_num; ++p_id)
	{
		Model_T2D_CHM_s::Particle &pcl = md.pcls[p_id];
		ParticleData &pd = pcl_data[p_id];
		pd.id = pcl.id;
		pd.from_pcl(pcl);
	}
	int res = rf.write_dataset(frame_id, "ParticleData", md.pcl_num,
							  size_t(sizeof(ParticleData)/sizeof(double)),
							  (double *)pcl_data);
	delete[] pcl_data;
	// particle num
	hid_t pcl_data_id = rf.open_dataset(frame_id, "ParticleData");
	rf.write_attribute(pcl_data_id, "pcl_num", md.pcl_num);
	rf.close_dataset(pcl_data_id);
	return res;
}

int load_pcl_data_from_hdf5_file(Model_T2D_CHM_s &md, ResultFile_hdf5 &rf, hid_t frame_id)
{
	// particle num
	size_t pcl_num;
	hid_t pcl_data_id = rf.open_dataset(frame_id, "ParticleData");
	rf.read_attribute(pcl_data_id, "pcl_num", pcl_num);
	rf.close_dataset(pcl_data_id);

	ParticleData *pcls_data = new ParticleData[pcl_num];
	int res = rf.read_dataset(frame_id, "ParticleData", pcl_num,
							  size_t(sizeof(ParticleData) / sizeof(double)),
							  (double *)pcls_data);
	if (res) return res;
	md.alloc_pcls(pcl_num);
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		ParticleData &pcl_data = pcls_data[p_id];
		Model_T2D_CHM_s::Particle &pcl = md.pcls[p_id];
		pcl.id = pcl_data.id;
		pcl_data.to_pcl(pcl);
	}
	delete[] pcls_data;
	return 0;
}


int load_chm_s_model_from_hdf5_file(Model_T2D_CHM_s &md,
	const char *hdf5_name, const char *th_name, size_t frame_id)
{
	ResultFile_hdf5 rf;
	rf.open(hdf5_name);
	hid_t file_id = rf.get_file_id();
	if (file_id < 0) return -1;
	
	// model data
	load_model_data_from_hdf5_file(md, rf);

	// time history output
	hid_t th_grp_id = rf.get_time_history_grp_id();
	hid_t th_id = rf.open_group(th_grp_id, th_name);
	char th_frame_name[30];
	snprintf(th_frame_name, 30, "frame_%zu", frame_id);
	hid_t th_frame_id = rf.open_group(th_id, th_frame_name);
	// particle data
	load_pcl_data_from_hdf5_file(md, rf, th_frame_id);
	// constitutive model data
	load_model_container_from_hdf5_file(md, rf, th_frame_id);
	// rigid body
	load_rigid_ciricle_to_hdf5_file(md.get_rigid_circle(), rf, th_frame_id);
	rf.close_group(th_frame_id);
	rf.close_group(th_id);

	return 0;
}
