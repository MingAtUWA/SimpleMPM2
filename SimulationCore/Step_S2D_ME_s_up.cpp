#include "SimulationCore_pcp.h"

#include <cstdio>
#include <Eigen/Sparse>

#define KEEP_NEWMARK_BETA_COEFFICIENT
#include "Step_S2D_ME_s_up.h"

using namespace Model_S2D_ME_s_up_Internal;

namespace
{
	typedef Model_S2D_ME_s_up::Particle Particle_mpm;
	typedef Model_S2D_ME_s_up::GaussPoint GaussPoint_mpm;
	typedef Model_S2D_ME_s_up::Element Element_mpm;
	typedef Model_S2D_ME_s_up::Node Node_mpm;
	typedef Model_S2D_ME_s_up::ShapeFuncValue ShapeFuncValue;
	typedef Model_S2D_ME_s_up::DOF DOF;

	void print_sparse_mat(Eigen::SparseMatrix<double> &mat,
		std::fstream &out_file, const char *mat_name = nullptr)
	{
		if (mat_name)
			out_file << mat_name << "\n";
		size_t row_num = mat.rows();
		size_t col_num = mat.cols();
		double value;
		for (size_t i = 0; i < row_num; ++i)
		{
			for (size_t j = 0; j < col_num; ++j)
			{
				value = mat.coeff(i, j);
				out_file << value << ", ";
			}
			out_file << "\n";
		}
	}

	void print_vec(Eigen::VectorXd &vec, size_t num,
				   std::fstream &out_file,
				   const char *vec_name = nullptr)
	{
		if (vec_name)
			out_file << vec_name << "\n";
		for (size_t i = 0; i < num; ++i)
			out_file << vec[i] << ", ";
		out_file << "\n";
	}
};

Step_S2D_ME_s_up::Step_S2D_ME_s_up() :
	//Step(&solve_substep_S2D_ME_s_up),
	Step(&solve_substep_S2D_ME_s_up_pure),
	model(nullptr) {}

Step_S2D_ME_s_up::~Step_S2D_ME_s_up() {}

int Step_S2D_ME_s_up::init_calculation(void)
{
	if (is_first_step)
	{
		for (size_t pcl_id = 0; pcl_id < model->pcl_num; ++pcl_id)
		{
			Particle_mpm &pcl = model->pcls[pcl_id];
			pcl.pe = model->elems;
		}
	}

	for (size_t pcl_id = 0; pcl_id < model->pcl_num; ++pcl_id)
	{
		Particle_mpm &pcl = model->pcls[pcl_id];
		pcl.ux = 0.0;
		pcl.uy = 0.0;
		pcl.x_ori = pcl.x;
		pcl.y_ori = pcl.y;
	}

	// for debug
	out_file.open("debug_mat_out_ME_s_up.csv", std::ios::binary | std::ios::out);
	
	return 0;
}

int Step_S2D_ME_s_up::finalize_calculation(void)
{
	// for debug
	out_file.close();

	return 0;
}

int solve_substep_S2D_ME_s_up_pure(void *_self)
{
	Step_S2D_ME_s_up &self = *(Step_S2D_ME_s_up *)(_self);
	Model_S2D_ME_s_up &model = *self.model;
	double dt = self.dtime;
	double dt2 = dt * dt;

	// init nodes
	size_t node_num = model.node_num;
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node_mpm &n = model.nodes[n_id];
		n.g_id = node_num;
	}
	// init elements
	for (size_t e_id = 0; e_id < model.elem_num; ++e_id)
	{
		Element_mpm &elem = model.elems[e_id];
		elem.pcls = nullptr;
	}

	// init particles
	for (size_t pcl_id = 0; pcl_id < model.pcl_num; ++pcl_id)
	{
		Particle_mpm &pcl = model.pcls[pcl_id];
		if (pcl.pe)
		{
			if (!model.init_pcl_cal_var(pcl))
				continue;
			pcl.vol = pcl.m / pcl.density;
			pcl.pe->add_pcl(pcl);
		}
	}
	
	size_t cal_node_num = 0;
	for (size_t e_id = 0; e_id < model.elem_num; ++e_id)
	{
		Element_mpm &e = model.elems[e_id];
		if (e.pcls)
		{
			Node_mpm &n1 = model.nodes[e.n1_id];
			if (n1.g_id == node_num)
			{
				n1.g_id = cal_node_num;
				++cal_node_num;
			}
			Node_mpm &n2 = model.nodes[e.n2_id];
			if (n2.g_id == node_num)
			{
				n2.g_id = cal_node_num;
				++cal_node_num;
			}
			Node_mpm &n3 = model.nodes[e.n3_id];
			if (n3.g_id == node_num)
			{
				n3.g_id = cal_node_num;
				++cal_node_num;
			}
			Node_mpm &n4 = model.nodes[e.n4_id];
			if (n4.g_id == node_num)
			{
				n4.g_id = cal_node_num;
				++cal_node_num;
			}
		}
	}

	// form g_id_map
	size_t *g_id_map = self.node_g_id_map_mem.resize(cal_node_num);
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node_mpm &n = model.nodes[n_id];
		if (n.g_id != node_num)
			g_id_map[n.g_id] = n_id;
	}

	size_t dof_num = cal_node_num * 3; // ux, uy, p for each node
	// list of non-zeros coefficients
	MatrixCoefficientSet<> &g_kmat_coefs = self.g_kmat_coefs; 
	g_kmat_coefs.init(dof_num);
	Eigen::VectorXd g_fvec(dof_num);
	g_fvec.setZero();

	size_t node_g_id, l2g_id_map[12];
	double e_kmat[12][12], e_fvec[12];
	for (size_t e_id = 0; e_id < model.elem_num; ++e_id)
	{
		Element_mpm &e = model.elems[e_id];
		if (e.pcls)
		{
			// form elemental stiffness matrix and force vector
			self.form_elem_stiffness_mat_and_force_vec_pure(e, e_kmat, e_fvec);

			// form map from global id to local id
			node_g_id = model.nodes[e.n1_id].g_id;
			l2g_id_map[0] = model.n_id_to_dof_id(node_g_id, DOF::ux);
			l2g_id_map[4] = model.n_id_to_dof_id(node_g_id, DOF::uy);
			l2g_id_map[8] = model.n_id_to_dof_id(node_g_id, DOF::p);
			node_g_id = model.nodes[e.n2_id].g_id;
			l2g_id_map[1] = model.n_id_to_dof_id(node_g_id, DOF::ux);
			l2g_id_map[5] = model.n_id_to_dof_id(node_g_id, DOF::uy);
			l2g_id_map[9] = model.n_id_to_dof_id(node_g_id, DOF::p);
			node_g_id = model.nodes[e.n3_id].g_id;
			l2g_id_map[2] = model.n_id_to_dof_id(node_g_id, DOF::ux);
			l2g_id_map[6] = model.n_id_to_dof_id(node_g_id, DOF::uy);
			l2g_id_map[10] = model.n_id_to_dof_id(node_g_id, DOF::p);
			node_g_id = model.nodes[e.n4_id].g_id;
			l2g_id_map[3] = model.n_id_to_dof_id(node_g_id, DOF::ux);
			l2g_id_map[7] = model.n_id_to_dof_id(node_g_id, DOF::uy);
			l2g_id_map[11] = model.n_id_to_dof_id(node_g_id, DOF::p);
			//print_vec(l2g_id_map);
			// add to global matrix and vector
			size_t g_id1, g_id2;
			for (size_t l_id1 = 0; l_id1 < 12; ++l_id1)
			{
				g_id1 = l2g_id_map[l_id1];
				for (size_t l_id2 = 0; l_id2 < 12; ++l_id2)
				{
					g_id2 = l2g_id_map[l_id2];
					g_kmat_coefs.add_coefficient(g_id1, g_id2, e_kmat[l_id1][l_id2]);
				}
				g_fvec[g_id1] += e_fvec[l_id1];
			}
		}
	}

	// apply external force
	// traction
	size_t dof_g_id;
	for (size_t t_id = 0; t_id < model.tx_num; ++t_id)
	{
		TractionBC_MPM &tx = model.txs[t_id];
		Particle_mpm &pcl = model.pcls[tx.pcl_id];
		if (pcl.pe)
		{
			Element_mpm &elem = *pcl.pe;
			// node 1
			node_g_id = model.nodes[elem.n1_id].g_id;
			dof_g_id = model.n_id_to_dof_id(node_g_id, DOF::ux);
			g_fvec[dof_g_id] += pcl.N1 * tx.t;
			// node 2
			node_g_id = model.nodes[elem.n2_id].g_id;
			dof_g_id = model.n_id_to_dof_id(node_g_id, DOF::ux);
			g_fvec[dof_g_id] += pcl.N2 * tx.t;
			// node 3
			node_g_id = model.nodes[elem.n3_id].g_id;
			dof_g_id = model.n_id_to_dof_id(node_g_id, DOF::ux);
			g_fvec[dof_g_id] += pcl.N3 * tx.t;
			// node 4
			node_g_id = model.nodes[elem.n4_id].g_id;
			dof_g_id = model.n_id_to_dof_id(node_g_id, DOF::ux);
			g_fvec[dof_g_id] += pcl.N4 * tx.t;
		}
	}
	for (size_t t_id = 0; t_id < model.ty_num; ++t_id)
	{
		TractionBC_MPM &ty = model.tys[t_id];
		Particle_mpm &pcl = model.pcls[ty.pcl_id];
		if (pcl.pe)
		{
			Element_mpm &elem = *pcl.pe;
			// node 1
			node_g_id = model.nodes[elem.n1_id].g_id;
			dof_g_id = model.n_id_to_dof_id(node_g_id, DOF::uy);
			g_fvec[dof_g_id] += pcl.N1 * ty.t;
			// node 2
			node_g_id = model.nodes[elem.n2_id].g_id;
			dof_g_id = model.n_id_to_dof_id(node_g_id, DOF::uy);
			g_fvec[dof_g_id] += pcl.N2 * ty.t;
			// node 3
			node_g_id = model.nodes[elem.n3_id].g_id;
			dof_g_id = model.n_id_to_dof_id(node_g_id, DOF::uy);
			g_fvec[dof_g_id] += pcl.N3 * ty.t;
			// node 4
			node_g_id = model.nodes[elem.n4_id].g_id;
			dof_g_id = model.n_id_to_dof_id(node_g_id, DOF::uy);
			g_fvec[dof_g_id] += pcl.N4 * ty.t;
		}
	}
	// body force to be finished...

	// apply displacement boundary condition
	double dig_term;
	double *kmat_col = self.kmat_col_mem.resize(dof_num);
	for (size_t bc_id = 0; bc_id < model.ux_num; ++bc_id)
	{
		DisplacementBC &dbc = model.uxs[bc_id];
		node_g_id = model.nodes[dbc.node_id].g_id;
		if (node_g_id != node_num)
		{
			dof_g_id = model.n_id_to_dof_id(node_g_id, DOF::ux);
			dig_term = g_kmat_coefs.del_col_and_row(dof_g_id, kmat_col);
			for (size_t row_id = 0; row_id < dof_num; ++row_id)
				g_fvec[row_id] -= kmat_col[row_id] * dbc.u;
			g_fvec[dof_g_id] = dig_term * dbc.u;
		}
	}
	for (size_t bc_id = 0; bc_id < model.uy_num; ++bc_id)
	{
		DisplacementBC &dbc = model.uys[bc_id];
		node_g_id = model.nodes[dbc.node_id].g_id;
		if (node_g_id != node_num)
		{
			dof_g_id = model.n_id_to_dof_id(node_g_id, DOF::uy);
			dig_term = g_kmat_coefs.del_col_and_row(dof_g_id, kmat_col);
			for (size_t row_id = 0; row_id < dof_num; ++row_id)
				g_fvec[row_id] -= kmat_col[row_id] * dbc.u;
			g_fvec[dof_g_id] = dig_term * dbc.u;
		}
	}
	for (size_t bc_id = 0; bc_id < model.pbc_num; ++bc_id)
	{
		PressureBC &pbc = model.pbcs[bc_id];
		node_g_id = model.nodes[pbc.node_id].g_id;
		if (node_g_id != node_num)
		{
			dof_g_id = model.n_id_to_dof_id(node_g_id, DOF::p);
			dig_term = g_kmat_coefs.del_col_and_row(dof_g_id, kmat_col);
			for (size_t row_id = 0; row_id < dof_num; ++row_id)
				g_fvec[row_id] -= kmat_col[row_id] * pbc.p;
			g_fvec[dof_g_id] = dig_term * pbc.p;
		}
	}

	Eigen::SparseMatrix<double> g_kmat(dof_num, dof_num);
	g_kmat.setFromTriplets(g_kmat_coefs.begin(), g_kmat_coefs.end());
	//print_sparse_mat(g_kmat, self.out_file, nullptr);
	Eigen::SparseLU<Eigen::SparseMatrix<double> > solver(g_kmat);
	Eigen::VectorXd g_du_vec = solver.solve(g_fvec);
	char time_str[30];
	sprintf(time_str, "%11.5lf", self.get_current_time());
	print_vec(g_du_vec, dof_num, self.out_file, time_str);

	// Update particle variables
	double dux1, duy1, dp1;
	double dux2, duy2, dp2;
	double dux3, duy3, dp3;
	double dux4, duy4, dp4;
	double dux, duy, dvx, dvy, dax, day;
	double de11, de22, de12;
	double ds11, ds22, ds12;
	double de_vol;
	for (size_t e_id = 0; e_id < model.elem_num; ++e_id)
	{
		Element_mpm &e = model.elems[e_id];
		if (e.pcls)
		{
			node_g_id = model.nodes[e.n1_id].g_id;
			dux1 = g_du_vec[model.n_id_to_dof_id(node_g_id, DOF::ux)];
			duy1 = g_du_vec[model.n_id_to_dof_id(node_g_id, DOF::uy)];
			dp1  = g_du_vec[model.n_id_to_dof_id(node_g_id, DOF::p)];
			node_g_id = model.nodes[e.n2_id].g_id;
			dux2 = g_du_vec[model.n_id_to_dof_id(node_g_id, DOF::ux)];
			duy2 = g_du_vec[model.n_id_to_dof_id(node_g_id, DOF::uy)];
			dp2  = g_du_vec[model.n_id_to_dof_id(node_g_id, DOF::p)];
			node_g_id = model.nodes[e.n3_id].g_id;
			dux3 = g_du_vec[model.n_id_to_dof_id(node_g_id, DOF::ux)];
			duy3 = g_du_vec[model.n_id_to_dof_id(node_g_id, DOF::uy)];
			dp3  = g_du_vec[model.n_id_to_dof_id(node_g_id, DOF::p)];
			node_g_id = model.nodes[e.n4_id].g_id;
			dux4 = g_du_vec[model.n_id_to_dof_id(node_g_id, DOF::ux)];
			duy4 = g_du_vec[model.n_id_to_dof_id(node_g_id, DOF::uy)];
			dp4  = g_du_vec[model.n_id_to_dof_id(node_g_id, DOF::p)];

			// update particle variables
			for (Particle_mpm *pcl_iter = e.pcls; pcl_iter; pcl_iter = pcl_iter->next)
			{
				Particle_mpm &pcl = *pcl_iter;

				// displacement increment
				dux = dux1 * pcl.N1 + dux2 * pcl.N2 + dux3 * pcl.N3 + dux4 * pcl.N4;
				duy = duy1 * pcl.N1 + duy2 * pcl.N2 + duy3 * pcl.N3 + duy4 * pcl.N4;
				// velocity increment
				dvx = gamma / (beta * dt) * dux - gamma / beta * pcl.vx
					+ dt * (1.0 - gamma / (2.0 * beta)) * pcl.ax;
				dvy = gamma / (beta * dt) * duy - gamma / beta * pcl.vy
					+ dt * (1.0 - gamma / (2.0 * beta)) * pcl.ay;
				// acceleration increment
				dax = dux / (beta * dt2) - pcl.vx / (beta * dt) - pcl.ax / (2.0 * beta);
				day = duy / (beta * dt2) - pcl.vy / (beta * dt) - pcl.ay / (2.0 * beta);
				// displacement
				pcl.ux += dux;
				pcl.uy += duy;
				// position
				pcl.x = pcl.x_ori + pcl.ux;
				pcl.y = pcl.y_ori + pcl.uy;
				// velocity
				pcl.vx += dvx;
				pcl.vy += dvy;
				// acceleration
				pcl.ax += dax;
				pcl.ay += day;

				// strain increment
				de11 = dux1 * pcl.dN1_dx + dux2 * pcl.dN2_dx
					 + dux3 * pcl.dN3_dx + dux4 * pcl.dN4_dx;
				de22 = duy1 * pcl.dN1_dy + duy2 * pcl.dN2_dy
					 + duy3 * pcl.dN3_dy + duy4 * pcl.dN4_dy;
				de12 = (dux1 * pcl.dN1_dy + dux2 * pcl.dN2_dy
					 + dux3 * pcl.dN3_dy + dux4 * pcl.dN4_dy
					 + duy1 * pcl.dN1_dx + duy2 * pcl.dN2_dx
					 + duy3 * pcl.dN3_dx + duy4 * pcl.dN4_dx) * 0.5;
				pcl.e11 += de11;
				pcl.e22 += de22;
				pcl.e12 += de12;

				// update stress
				double E_tmp = model.E / (3.0 * (1.0 + model.niu));
				ds11 = 2.0 * E_tmp * de11 - E_tmp * de22;
				ds22 = -E_tmp * de11 + 2.0 * E_tmp * de22;
				ds12 = model.E / (2.0 * (1.0 + model.niu)) * 2.0 * de12;
				pcl.s11 += ds11;
				pcl.s22 += ds22;
				pcl.s12 += ds12;

				pcl.p += dp1 * pcl.N1 + dp2 * pcl.N2 + dp3 * pcl.N3 + dp4 * pcl.N4;

				// density
				de_vol = de11 + de22;
				pcl.density /= 1.0 + de_vol;
			}
		}
	}

	return 0;
}

int solve_substep_S2D_ME_s_up(void *_self)
{
	Step_S2D_ME_s_up &self = *(Step_S2D_ME_s_up *)(_self);
	Model_S2D_ME_s_up &model = *self.model;
	double dt = self.dtime;
	double dt2 = dt * dt;
	
	// init nodes
	size_t node_num = model.node_num;
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node_mpm &n = model.nodes[n_id];
		n.g_id = node_num;
		n.vol = 0.0;
		n.density = 0.0;
		n.s11 = 0.0;
		n.s22 = 0.0;
		n.s12 = 0.0;
		n.p = 0.0;
	}
	// init elements
	for (size_t e_id = 0; e_id < model.elem_num; ++e_id)
	{
		Element_mpm &elem = model.elems[e_id];
		elem.pcls = nullptr;
		elem.vf = 0.0;
	}
	// init particles
	for (size_t pcl_id = 0; pcl_id < model.pcl_num; ++pcl_id)
	{
		Particle_mpm &pcl = model.pcls[pcl_id];
		if (pcl.pe)
		{
			if (!model.init_pcl_cal_var(pcl))
				continue;
			pcl.vol = pcl.m / pcl.density;
			// element
			Element_mpm &elem = *pcl.pe;
			elem.add_pcl(pcl);
			elem.vf += pcl.vol;
			// map variables from pcls to nodes
			Node_mpm &n1 = model.nodes[elem.n1_id];
			double vol_N1 = pcl.vol * pcl.N1;
			n1.vol += vol_N1;
			n1.density += pcl.density * vol_N1;
			n1.s11 += pcl.s11 * vol_N1;
			n1.s22 += pcl.s22 * vol_N1;
			n1.s12 += pcl.s12 * vol_N1;
			n1.p += pcl.p * vol_N1;
			Node_mpm &n2 = model.nodes[elem.n2_id];
			double vol_N2 = pcl.vol * pcl.N2;
			n2.vol += vol_N2;
			n2.density += pcl.density * vol_N2;
			n2.s11 += pcl.s11 * vol_N2;
			n2.s22 += pcl.s22 * vol_N2;
			n2.s12 += pcl.s12 * vol_N2;
			n2.p += pcl.p * vol_N2;
			Node_mpm &n3 = model.nodes[elem.n3_id];
			double vol_N3 = pcl.vol * pcl.N3;
			n3.vol += vol_N3;
			n3.density += pcl.density * vol_N3;
			n3.s11 += pcl.s11 * vol_N3;
			n3.s22 += pcl.s22 * vol_N3;
			n3.s12 += pcl.s12 * vol_N3;
			n3.p += pcl.p * vol_N3;
			Node_mpm &n4 = model.nodes[elem.n4_id];
			double vol_N4 = pcl.vol * pcl.N4;
			n4.vol += vol_N4;
			n4.density += pcl.density * vol_N4;
			n4.s11 += pcl.s11 * vol_N4;
			n4.s22 += pcl.s22 * vol_N4;
			n4.s12 += pcl.s12 * vol_N4;
			n4.p += pcl.p * vol_N4;
		}
	}
	
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node_mpm &n = model.nodes[n_id];
		if (n.vol != 0.0)
		{
			n.density /= n.vol;
			n.s11 /= n.vol;
			n.s22 /= n.vol;
			n.s12 /= n.vol;
			n.p /= n.vol;
		}
	}

	size_t cal_node_num = 0;
	ShapeFuncValue &sf1 = model.sf1;
	ShapeFuncValue &sf2 = model.sf2;
	ShapeFuncValue &sf3 = model.sf3;
	ShapeFuncValue &sf4 = model.sf4;
	for (size_t e_id = 0; e_id < model.elem_num; ++e_id)
	{
		Element_mpm &e = model.elems[e_id];
		if (e.pcls)
		{
			e.vf /= model.elem_vol;
			Node_mpm &n1 = model.nodes[e.n1_id];
			if (n1.g_id == node_num)
			{
				n1.g_id = cal_node_num;
				++cal_node_num;
			}
			Node_mpm &n2 = model.nodes[e.n2_id];
			if (n2.g_id == node_num)
			{
				n2.g_id = cal_node_num;
				++cal_node_num;
			}
			Node_mpm &n3 = model.nodes[e.n3_id];
			if (n3.g_id == node_num)
			{
				n3.g_id = cal_node_num;
				++cal_node_num;
			}
			Node_mpm &n4 = model.nodes[e.n4_id];
			if (n4.g_id == node_num)
			{
				n4.g_id = cal_node_num;
				++cal_node_num;
			}
			// gauss point
			GaussPoint_mpm &gp1 = e.gp1;
			gp1.density = sf1.N1 * n1.density + sf2.N1 * n2.density
						+ sf3.N1 * n3.density + sf4.N1 * n4.density;
			gp1.s11 = sf1.N1 * n1.s11 + sf2.N1 * n2.s11
					+ sf3.N1 * n3.s11 + sf4.N1 * n4.s11;
			gp1.s22 = sf1.N1 * n1.s22 + sf2.N1 * n2.s22
					+ sf3.N1 * n3.s22 + sf4.N1 * n4.s22;
			gp1.s12 = sf1.N1 * n1.s12 + sf2.N1 * n2.s12
					+ sf3.N1 * n3.s12 + sf4.N1 * n4.s12;
			gp1.p = sf1.N1 * n1.p + sf2.N1 * n2.p + sf3.N1 * n3.p + sf4.N1 * n4.p;
			GaussPoint_mpm &gp2 = e.gp2;
			gp2.density = sf1.N2 * n1.density + sf2.N2 * n2.density
						+ sf3.N2 * n3.density + sf4.N2 * n4.density;
			gp2.s11 = sf1.N2 * n1.s11 + sf2.N2 * n2.s11
					+ sf3.N2 * n3.s11 + sf4.N2 * n4.s11;
			gp2.s22 = sf1.N2 * n1.s22 + sf2.N2 * n2.s22
					+ sf3.N2 * n3.s22 + sf4.N2 * n4.s22;
			gp2.s12 = sf1.N2 * n1.s12 + sf2.N2 * n2.s12
					+ sf3.N2 * n3.s12 + sf4.N2 * n4.s12;
			gp2.p = sf1.N2 * n1.p + sf2.N2 * n2.p + sf3.N2 * n3.p + sf4.N2 * n4.p;
			GaussPoint_mpm &gp3 = e.gp3;
			gp3.density = sf1.N3 * n1.density + sf2.N3 * n2.density
						+ sf3.N3 * n3.density + sf4.N3 * n4.density;
			gp3.s11 = sf1.N3 * n1.s11 + sf2.N3 * n2.s11
					+ sf3.N3 * n3.s11 + sf4.N3 * n4.s11;
			gp3.s22 = sf1.N3 * n1.s22 + sf2.N3 * n2.s22
					+ sf3.N3 * n3.s22 + sf4.N3 * n4.s22;
			gp3.s12 = sf1.N3 * n1.s12 + sf2.N3 * n2.s12
					+ sf3.N3 * n3.s12 + sf4.N3 * n4.s12;
			gp3.p = sf1.N3 * n1.p + sf2.N3 * n2.p + sf3.N3 * n3.p + sf4.N3 * n4.p;
			GaussPoint_mpm &gp4 = e.gp4;
			gp4.density = sf1.N4 * n1.density + sf2.N4 * n2.density
						+ sf3.N4 * n3.density + sf4.N4 * n4.density;
			gp4.s11 = sf1.N4 * n1.s11 + sf2.N4 * n2.s11
					+ sf3.N4 * n3.s11 + sf4.N4 * n4.s11;
			gp4.s22 = sf1.N4 * n1.s22 + sf2.N4 * n2.s22
					+ sf3.N4 * n3.s22 + sf4.N4 * n4.s22;
			gp4.s12 = sf1.N4 * n1.s12 + sf2.N4 * n2.s12
					+ sf3.N4 * n3.s12 + sf4.N4 * n4.s12;
			gp4.p = sf1.N4 * n1.p + sf2.N4 * n2.p + sf3.N4 * n3.p + sf4.N4 * n4.p;
		}
	}
	
	// form g_id_map
	size_t *g_id_map = self.node_g_id_map_mem.resize(cal_node_num);
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node_mpm &n = model.nodes[n_id];
		if (n.g_id != node_num)
			g_id_map[n.g_id] = n_id;
	}
	
	size_t dof_num = cal_node_num * 3; // ux, uy, p for each node
	// list of non-zeros coefficients
	MatrixCoefficientSet<> &g_kmat_coefs = self.g_kmat_coefs;
	g_kmat_coefs.init(dof_num);
	Eigen::VectorXd g_fvec(dof_num);
	g_fvec.setZero();

	size_t node_g_id;
	size_t n_ids[4], l2g_id_map[12];
	double e_kmat[12][12], e_fvec[12];
	for (size_t e_id = 0; e_id < model.elem_num; ++e_id)
	{
		Element_mpm &e = model.elems[e_id];
		if (e.pcls)
		{
			// form elemental stiffness matrix and force vector
			self.form_elem_stiffness_mat_and_force_vec(e, e_kmat, e_fvec);

			// form map from global id to local id
			node_g_id = model.nodes[e.n1_id].g_id;
			l2g_id_map[0] = model.n_id_to_dof_id(node_g_id, DOF::ux);
			l2g_id_map[4] = model.n_id_to_dof_id(node_g_id, DOF::uy);
			l2g_id_map[8] = model.n_id_to_dof_id(node_g_id, DOF::p);
			node_g_id = model.nodes[e.n2_id].g_id;
			l2g_id_map[1] = model.n_id_to_dof_id(node_g_id, DOF::ux);
			l2g_id_map[5] = model.n_id_to_dof_id(node_g_id, DOF::uy);
			l2g_id_map[9] = model.n_id_to_dof_id(node_g_id, DOF::p);
			node_g_id = model.nodes[e.n3_id].g_id;
			l2g_id_map[2] = model.n_id_to_dof_id(node_g_id, DOF::ux);
			l2g_id_map[6] = model.n_id_to_dof_id(node_g_id, DOF::uy);
			l2g_id_map[10] = model.n_id_to_dof_id(node_g_id, DOF::p);
			node_g_id = model.nodes[e.n4_id].g_id;
			l2g_id_map[3] = model.n_id_to_dof_id(node_g_id, DOF::ux);
			l2g_id_map[7] = model.n_id_to_dof_id(node_g_id, DOF::uy);
			l2g_id_map[11] = model.n_id_to_dof_id(node_g_id, DOF::p);
			//print_vec(l2g_id_map);
			// add to global matrix and vector
			size_t g_id1, g_id2;
			for (size_t l_id1 = 0; l_id1 < 12; ++l_id1)
			{
				g_id1 = l2g_id_map[l_id1];
				for (size_t l_id2 = 0; l_id2 < 12; ++l_id2)
				{
					g_id2 = l2g_id_map[l_id2];
					g_kmat_coefs.add_coefficient(g_id1, g_id2, e_kmat[l_id1][l_id2]);
				}
				g_fvec[g_id1] += e_fvec[l_id1];
			}
		}
	}
	//std::cout << g_fvec << "\n";

	// apply external force
	// traction
	size_t dof_g_id;
	for (size_t t_id = 0; t_id < model.tx_num; ++t_id)
	{
		TractionBC_MPM &tx = model.txs[t_id];
		Particle_mpm &pcl = model.pcls[tx.pcl_id];
		if (pcl.pe)
		{
			Element_mpm &elem = *pcl.pe;
			// node 1
			node_g_id = model.nodes[elem.n1_id].g_id;
			dof_g_id = model.n_id_to_dof_id(node_g_id, DOF::ux);
			g_fvec[dof_g_id] += pcl.N1 * tx.t;
			// node 2
			node_g_id = model.nodes[elem.n2_id].g_id;
			dof_g_id = model.n_id_to_dof_id(node_g_id, DOF::ux);
			g_fvec[dof_g_id] += pcl.N2 * tx.t;
			// node 3
			node_g_id = model.nodes[elem.n3_id].g_id;
			dof_g_id = model.n_id_to_dof_id(node_g_id, DOF::ux);
			g_fvec[dof_g_id] += pcl.N3 * tx.t;
			// node 4
			node_g_id = model.nodes[elem.n4_id].g_id;
			dof_g_id = model.n_id_to_dof_id(node_g_id, DOF::ux);
			g_fvec[dof_g_id] += pcl.N4 * tx.t;
		}
	}
	for (size_t t_id = 0; t_id < model.ty_num; ++t_id)
	{
		TractionBC_MPM &ty = model.tys[t_id];
		Particle_mpm &pcl = model.pcls[ty.pcl_id];
		if (pcl.pe)
		{
			Element_mpm &elem = *pcl.pe;
			// node 1
			node_g_id = model.nodes[elem.n1_id].g_id;
			dof_g_id = model.n_id_to_dof_id(node_g_id, DOF::uy);
			g_fvec[dof_g_id] += pcl.N1 * ty.t;
			// node 2
			node_g_id = model.nodes[elem.n2_id].g_id;
			dof_g_id = model.n_id_to_dof_id(node_g_id, DOF::uy);
			g_fvec[dof_g_id] += pcl.N2 * ty.t;
			// node 3
			node_g_id = model.nodes[elem.n3_id].g_id;
			dof_g_id = model.n_id_to_dof_id(node_g_id, DOF::uy);
			g_fvec[dof_g_id] += pcl.N3 * ty.t;
			// node 4
			node_g_id = model.nodes[elem.n4_id].g_id;
			dof_g_id = model.n_id_to_dof_id(node_g_id, DOF::uy);
			g_fvec[dof_g_id] += pcl.N4 * ty.t;
		}
	}
	// body force to be finished...

	// apply displacement boundary condition
	double dig_term;
	double *kmat_col = self.kmat_col_mem.resize(dof_num);
	for (size_t bc_id = 0; bc_id < model.ux_num; ++bc_id)
	{
		DisplacementBC &dbc = model.uxs[bc_id];
		node_g_id = model.nodes[dbc.node_id].g_id;
		if (node_g_id != node_num)
		{
			dof_g_id = model.n_id_to_dof_id(node_g_id, DOF::ux);
			dig_term = g_kmat_coefs.del_col_and_row(dof_g_id, kmat_col);
			for (size_t row_id = 0; row_id < dof_num; ++row_id)
				g_fvec[row_id] -= kmat_col[row_id] * dbc.u;
			g_fvec[dof_g_id] = dig_term * dbc.u;
		}
	}
	for (size_t bc_id = 0; bc_id < model.uy_num; ++bc_id)
	{
		DisplacementBC &dbc = model.uys[bc_id];
		node_g_id = model.nodes[dbc.node_id].g_id;
		if (node_g_id != node_num)
		{
			dof_g_id = model.n_id_to_dof_id(node_g_id, DOF::uy);
			dig_term = g_kmat_coefs.del_col_and_row(dof_g_id, kmat_col);
			for (size_t row_id = 0; row_id < dof_num; ++row_id)
				g_fvec[row_id] -= kmat_col[row_id] * dbc.u;
			g_fvec[dof_g_id] = dig_term * dbc.u;
		}
	}
	for (size_t bc_id = 0; bc_id < model.pbc_num; ++bc_id)
	{
		PressureBC &pbc = model.pbcs[bc_id];
		node_g_id = model.nodes[pbc.node_id].g_id;
		if (node_g_id != node_num)
		{
			dof_g_id = model.n_id_to_dof_id(node_g_id, DOF::p);
			dig_term = g_kmat_coefs.del_col_and_row(dof_g_id, kmat_col);
			for (size_t row_id = 0; row_id < dof_num; ++row_id)
				g_fvec[row_id] -= kmat_col[row_id] * pbc.p;
			g_fvec[dof_g_id] = dig_term * pbc.p;
		}
	}

	Eigen::SparseMatrix<double> g_kmat(dof_num, dof_num);
	g_kmat.setFromTriplets(g_kmat_coefs.begin(), g_kmat_coefs.end());
	//print_sparse_mat(g_kmat, self.out_file, nullptr);
	Eigen::SparseLU<Eigen::SparseMatrix<double> > solver(g_kmat);
	Eigen::VectorXd g_du_vec = solver.solve(g_fvec);
	//std::cout << g_fvec << "\n";
	//std::cout << g_du_vec << "\n";

	// Update particle variables
	double dux1, duy1, dp1;
	double dux2, duy2, dp2;
	double dux3, duy3, dp3;
	double dux4, duy4, dp4;
	double dux, duy, dvx, dvy, dax, day;
	double de11, de22, de12;
	double ds11, ds22, ds12;
	double de_vol;
	for (size_t e_id = 0; e_id < model.elem_num; ++e_id)
	{
		Element_mpm &e = model.elems[e_id];
		if (e.pcls)
		{
			node_g_id = model.nodes[e.n1_id].g_id;
			dux1 = g_du_vec[model.n_id_to_dof_id(node_g_id, DOF::ux)];
			duy1 = g_du_vec[model.n_id_to_dof_id(node_g_id, DOF::uy)];
			dp1  = g_du_vec[model.n_id_to_dof_id(node_g_id, DOF::p)];
			node_g_id = model.nodes[e.n2_id].g_id;
			dux2 = g_du_vec[model.n_id_to_dof_id(node_g_id, DOF::ux)];
			duy2 = g_du_vec[model.n_id_to_dof_id(node_g_id, DOF::uy)];
			dp2  = g_du_vec[model.n_id_to_dof_id(node_g_id, DOF::p)];
			node_g_id = model.nodes[e.n3_id].g_id;
			dux3 = g_du_vec[model.n_id_to_dof_id(node_g_id, DOF::ux)];
			duy3 = g_du_vec[model.n_id_to_dof_id(node_g_id, DOF::uy)];
			dp3  = g_du_vec[model.n_id_to_dof_id(node_g_id, DOF::p)];
			node_g_id = model.nodes[e.n4_id].g_id;
			dux4 = g_du_vec[model.n_id_to_dof_id(node_g_id, DOF::ux)];
			duy4 = g_du_vec[model.n_id_to_dof_id(node_g_id, DOF::uy)];
			dp4  = g_du_vec[model.n_id_to_dof_id(node_g_id, DOF::p)];

			// update particle variables
			for (Particle_mpm *pcl_iter = e.pcls; pcl_iter; pcl_iter = pcl_iter->next)
			{
				Particle_mpm &pcl = *pcl_iter;

				// displacement increment
				dux = dux1 * pcl.N1 + dux2 * pcl.N2 + dux3 * pcl.N3 + dux4 * pcl.N4;
				duy = duy1 * pcl.N1 + duy2 * pcl.N2 + duy3 * pcl.N3 + duy4 * pcl.N4;
				// velocity increment
				dvx = gamma / (beta * dt) * dux - gamma / beta * pcl.vx
					+ dt * (1.0 - gamma / (2.0 * beta)) * pcl.ax;
				dvy = gamma / (beta * dt) * duy	- gamma / beta * pcl.vy
					+ dt * (1.0 - gamma / (2.0 * beta)) * pcl.ay;
				// acceleration increment
				dax = dux / (beta * dt2) - pcl.vx / (beta * dt) - pcl.ax / (2.0 * beta);
				day = duy / (beta * dt2) - pcl.vy / (beta * dt) - pcl.ay / (2.0 * beta);
				// displacement
				pcl.ux += dux;
				pcl.uy += duy;
				// position
				pcl.x = pcl.x_ori + pcl.ux;
				pcl.y = pcl.y_ori + pcl.uy;
				// velocity
				pcl.vx += dvx;
				pcl.vy += dvy;
				// acceleration
				pcl.ax += dax;
				pcl.ay += day;

				// strain increment
				de11 = dux1 * pcl.dN1_dx + dux2 * pcl.dN2_dx
					 + dux3 * pcl.dN3_dx + dux4 * pcl.dN4_dx;
				de22 = duy1 * pcl.dN1_dy + duy2 * pcl.dN2_dy
					 + duy3 * pcl.dN3_dy + duy4 * pcl.dN4_dy;
				de12 = (dux1 * pcl.dN1_dy + dux2 * pcl.dN2_dy
					  + dux3 * pcl.dN3_dy + dux4 * pcl.dN4_dy
					  + duy1 * pcl.dN1_dx + duy2 * pcl.dN2_dx
					  + duy3 * pcl.dN3_dx + duy4 * pcl.dN4_dx) * 0.5;
				pcl.e11 += de11;
				pcl.e22 += de22;
				pcl.e12 += de12;

				// update stress
				double E_tmp = model.E / (3.0 * (1.0 + model.niu));
				ds11 =  2.0 * E_tmp * de11 - E_tmp * de22;
				ds22 = -E_tmp * de11 + 2.0 * E_tmp * de22;
				ds12 = model.E / (2.0 * (1.0 + model.niu)) * 2.0 * de12;
				pcl.s11 += ds11;
				pcl.s22 += ds22;
				pcl.s12 += ds12;

				pcl.p += dp1 * pcl.N1 + dp2 * pcl.N2 + dp3 * pcl.N3 + dp4 * pcl.N4;

				// density
				de_vol = de11 + de22;
				pcl.density /= (1.0 + de_vol);
			}
		}
	}

	return 0;
}
