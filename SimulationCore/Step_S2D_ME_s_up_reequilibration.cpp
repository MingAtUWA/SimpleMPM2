#include "SimulationCore_pcp.h"

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

	void print_vec(double *vec, size_t num,
				   std::fstream &out_file,
				   const char *vec_name = nullptr)
	{
		if (vec_name)
			out_file << vec_name << "\n";
		for (size_t i = 0; i < num; ++i)
			out_file << vec[i] << ", ";
		out_file << "\n";
	}

	void print_vec(size_t *vec, size_t num)
	{
		std::cout << "\n";
		for (size_t i = 0; i < num; ++i)
			std::cout << vec[i] << ", ";
		std::cout << "\n";
	}
};

int substep_requilibration_S2D_ME_s_up(void *_self, Eigen::VectorXd& g_fvec)
{
	Step_S2D_ME_s_up &self = *(Step_S2D_ME_s_up *)(_self);
	Model_S2D_ME_s_up &model = *self.model;
	size_t node_g_id, l2g_id_map[12];
	double e_kmat[12][12], e_fvec[12];

	size_t dof_num = self.cal_node_num * 3; // ux, uy, p for each node
	MatrixCoefficientSet<> &g_kmat_coefs = self.g_kmat_coefs;
	g_kmat_coefs.init(dof_num);
	g_fvec.setZero();
	
	for (size_t e_id = 0; e_id < model.elem_num; ++e_id)
	{
		Element_mpm &e = model.elems[e_id];
		if (e.pcls)
		{
			// form stiffness matrix and force vector
			self.reequilibration_elem_stiffness_mat_and_force_vec(e, e_kmat, e_fvec);

			// form map from global id to local id
			node_g_id = model.nodes[e.n1_id].g_id;
			l2g_id_map[0] = self.n_id_to_dof_id(node_g_id, DOF::ux);
			l2g_id_map[4] = self.n_id_to_dof_id(node_g_id, DOF::uy);
			l2g_id_map[8] = self.n_id_to_dof_id(node_g_id, DOF::p);
			node_g_id = model.nodes[e.n2_id].g_id;
			l2g_id_map[1] = self.n_id_to_dof_id(node_g_id, DOF::ux);
			l2g_id_map[5] = self.n_id_to_dof_id(node_g_id, DOF::uy);
			l2g_id_map[9] = self.n_id_to_dof_id(node_g_id, DOF::p);
			node_g_id = model.nodes[e.n3_id].g_id;
			l2g_id_map[2] = self.n_id_to_dof_id(node_g_id, DOF::ux);
			l2g_id_map[6] = self.n_id_to_dof_id(node_g_id, DOF::uy);
			l2g_id_map[10] = self.n_id_to_dof_id(node_g_id, DOF::p);
			node_g_id = model.nodes[e.n4_id].g_id;
			l2g_id_map[3] = self.n_id_to_dof_id(node_g_id, DOF::ux);
			l2g_id_map[7] = self.n_id_to_dof_id(node_g_id, DOF::uy);
			l2g_id_map[11] = self.n_id_to_dof_id(node_g_id, DOF::p);
			// add to global matrix and vector
			for (size_t l_id1 = 0; l_id1 < 12; ++l_id1)
			{
				size_t g_id1 = l2g_id_map[l_id1];
				for (size_t l_id2 = 0; l_id2 < 12; ++l_id2)
				{
					size_t g_id2 = l2g_id_map[l_id2];
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
			Element_mpm &e = *pcl.pe;
			// node 1
			node_g_id = model.nodes[e.n1_id].g_id;
			dof_g_id = self.n_id_to_dof_id(node_g_id, DOF::ux);
			g_fvec[dof_g_id] += pcl.N1 * tx.t;
			// node 2
			node_g_id = model.nodes[e.n2_id].g_id;
			dof_g_id = self.n_id_to_dof_id(node_g_id, DOF::ux);
			g_fvec[dof_g_id] += pcl.N2 * tx.t;
			// node 3
			node_g_id = model.nodes[e.n3_id].g_id;
			dof_g_id = self.n_id_to_dof_id(node_g_id, DOF::ux);
			g_fvec[dof_g_id] += pcl.N3 * tx.t;
			// node 4
			node_g_id = model.nodes[e.n4_id].g_id;
			dof_g_id = self.n_id_to_dof_id(node_g_id, DOF::ux);
			g_fvec[dof_g_id] += pcl.N4 * tx.t;
		}
	}
	for (size_t t_id = 0; t_id < model.ty_num; ++t_id)
	{
		TractionBC_MPM &ty = model.tys[t_id];
		Particle_mpm &pcl = model.pcls[ty.pcl_id];
		if (pcl.pe)
		{
			Element_mpm &e = *pcl.pe;
			// node 1
			node_g_id = model.nodes[e.n1_id].g_id;
			dof_g_id = self.n_id_to_dof_id(node_g_id, DOF::uy);
			g_fvec[dof_g_id] += pcl.N1 * ty.t;
			// node 2
			node_g_id = model.nodes[e.n2_id].g_id;
			dof_g_id = self.n_id_to_dof_id(node_g_id, DOF::uy);
			g_fvec[dof_g_id] += pcl.N2 * ty.t;
			// node 3
			node_g_id = model.nodes[e.n3_id].g_id;
			dof_g_id = self.n_id_to_dof_id(node_g_id, DOF::uy);
			g_fvec[dof_g_id] += pcl.N3 * ty.t;
			// node 4
			node_g_id = model.nodes[e.n4_id].g_id;
			dof_g_id = self.n_id_to_dof_id(node_g_id, DOF::uy);
			g_fvec[dof_g_id] += pcl.N4 * ty.t;
		}
	}
	// body force to be finished...

	// apply displacement boundary condition
	size_t node_num = model.node_num;
	double dig_term;
	double *kmat_col = self.kmat_col_mem.resize(dof_num);
	for (size_t bc_id = 0; bc_id < model.ux_num; ++bc_id)
	{
		DisplacementBC &dbc = model.uxs[bc_id];
		node_g_id = model.nodes[dbc.node_id].g_id;
		if (node_g_id != node_num)
		{
			dof_g_id = self.n_id_to_dof_id(node_g_id, DOF::ux);
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
			dof_g_id = self.n_id_to_dof_id(node_g_id, DOF::uy);
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
			dof_g_id = self.n_id_to_dof_id(node_g_id, DOF::p);
			dig_term = g_kmat_coefs.del_col_and_row(dof_g_id, kmat_col);
			for (size_t row_id = 0; row_id < dof_num; ++row_id)
				g_fvec[row_id] -= kmat_col[row_id] * pbc.p;
			g_fvec[dof_g_id] = dig_term * pbc.p;
		}
	}

	//if (self.substep_num == 11)
	//	int efef = 0;

	Eigen::SparseMatrix<double> g_kmat(dof_num, dof_num);
	g_kmat.setFromTriplets(g_kmat_coefs.begin(), g_kmat_coefs.end());

	//// output for debug
	//char time_str[30];
	//if (self.substep_num >= 10)
	//{
	//	sprintf(time_str, "kmat %11.5lf s", self.get_current_time());
	//	print_sparse_mat(g_kmat, self.out_file, time_str);
	//	sprintf(time_str, "nf %11.5lf s", self.get_current_time());
	//	print_vec(g_fvec, dof_num, self.out_file, time_str);
	//}

	Eigen::SparseLU<Eigen::SparseMatrix<double> > solver(g_kmat);
	Eigen::VectorXd g_du_vec = solver.solve(g_fvec);
	//sprintf(time_str, "du %11.5lf s", self.get_current_time());
	//print_vec(g_du_vec, dof_num, self.out_file, time_str);

	// reapply disp bc
	for (size_t bc_id = 0; bc_id < model.ux_num; ++bc_id)
	{
		DisplacementBC &dbc = model.uxs[bc_id];
		node_g_id = model.nodes[dbc.node_id].g_id;
		if (node_g_id != node_num)
		{
			dof_g_id = self.n_id_to_dof_id(node_g_id, DOF::ux);
			g_du_vec[dof_g_id] = dbc.u;
		}
	}
	for (size_t bc_id = 0; bc_id < model.uy_num; ++bc_id)
	{
		DisplacementBC &dbc = model.uys[bc_id];
		node_g_id = model.nodes[dbc.node_id].g_id;
		if (node_g_id != node_num)
		{
			dof_g_id = self.n_id_to_dof_id(node_g_id, DOF::uy);
			g_du_vec[dof_g_id] = dbc.u;
		}
	}

	// Update particle variables
	double dux1, duy1, dp1;
	double dux2, duy2, dp2;
	double dux3, duy3, dp3;
	double dux4, duy4, dp4;
	double de11, de22, de12;
	double ds11, ds22, ds12, dp;
	for (size_t e_id = 0; e_id < model.elem_num; ++e_id)
	{
		Element_mpm &e = model.elems[e_id];
		if (e.pcls)
		{
			node_g_id = model.nodes[e.n1_id].g_id;
			dux1 = g_du_vec[self.n_id_to_dof_id(node_g_id, DOF::ux)];
			duy1 = g_du_vec[self.n_id_to_dof_id(node_g_id, DOF::uy)];
			dp1  = g_du_vec[self.n_id_to_dof_id(node_g_id, DOF::p)];
			node_g_id = model.nodes[e.n2_id].g_id;
			dux2 = g_du_vec[self.n_id_to_dof_id(node_g_id, DOF::ux)];
			duy2 = g_du_vec[self.n_id_to_dof_id(node_g_id, DOF::uy)];
			dp2  = g_du_vec[self.n_id_to_dof_id(node_g_id, DOF::p)];
			node_g_id = model.nodes[e.n3_id].g_id;
			dux3 = g_du_vec[self.n_id_to_dof_id(node_g_id, DOF::ux)];
			duy3 = g_du_vec[self.n_id_to_dof_id(node_g_id, DOF::uy)];
			dp3  = g_du_vec[self.n_id_to_dof_id(node_g_id, DOF::p)];
			node_g_id = model.nodes[e.n4_id].g_id;
			dux4 = g_du_vec[self.n_id_to_dof_id(node_g_id, DOF::ux)];
			duy4 = g_du_vec[self.n_id_to_dof_id(node_g_id, DOF::uy)];
			dp4  = g_du_vec[self.n_id_to_dof_id(node_g_id, DOF::p)];

			// update particle variables
			for (Particle_mpm *pcl_iter = e.pcls; pcl_iter; pcl_iter = pcl_iter->next)
			{
				Particle_mpm &pcl = *pcl_iter;

				// strain increment
				de11 = dux1 * pcl.dN1_dx + dux2 * pcl.dN2_dx
					 + dux3 * pcl.dN3_dx + dux4 * pcl.dN4_dx;
				de22 = duy1 * pcl.dN1_dy + duy2 * pcl.dN2_dy
					 + duy3 * pcl.dN3_dy + duy4 * pcl.dN4_dy;
				de12 = (dux1 * pcl.dN1_dy + dux2 * pcl.dN2_dy
					  + dux3 * pcl.dN3_dy + dux4 * pcl.dN4_dy
					  + duy1 * pcl.dN1_dx + duy2 * pcl.dN2_dx
					  + duy3 * pcl.dN3_dx + duy4 * pcl.dN4_dx) * 0.5;

				// update stress
				double G = model.E / (2.0 * (1.0 + model.niu));
				ds11 =  4.0 / 3.0 * G * de11 - 2.0 / 3.0 * G * de22;
				ds22 = -2.0 / 3.0 * G * de11 + 4.0 / 3.0 * G * de22;
				ds12 = G * 2.0 * de12;
				pcl.s11 += ds11;
				pcl.s22 += ds22;
				pcl.s12 += ds12;
				pcl.p += dp1 * pcl.N1 + dp2 * pcl.N2 + dp3 * pcl.N3 + dp4 * pcl.N4;
			}
		}
	}

	return 0;
}
