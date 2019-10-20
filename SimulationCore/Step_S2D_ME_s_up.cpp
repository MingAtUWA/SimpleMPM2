#include "SimulationCore_pcp.h"

#include <cstdio>
#include <Eigen/Sparse>

#include "Step_S2D_ME_s_up.h"

// constant of Newmark-beta method
#define beta 0.25
#define gamma 0.5

using namespace Model_S2D_ME_s_up_Internal;

namespace
{
	typedef Model_S2D_ME_s_up::Particle Particle_mpm;
	typedef Model_S2D_ME_s_up::Element Element_mpm;
	typedef Model_S2D_ME_s_up::Node Node_mpm;

	void print_mat(double mat[12][12])
	{
		std::cout << "mat\n";
		for (size_t i = 0; i < 12; i++)
		{
			for (size_t j = 0; j < 12; j++)
				printf("%+8.2e ", mat[i][j]);
			std::cout << "\n";
		}
	}

	void print_vec(double *vec, size_t num)
	{
		std::cout << "vec\n";
		for (size_t i = 0; i < num; i++)
			printf("%+8.2e ", vec[i]);
		std::cout << "\n";
	}

	void print_vec(size_t vec[12])
	{
		std::cout << "vec\n";
		for (size_t i = 0; i < 12; i++)
			std::cout << vec[i] << " ";
		std::cout << "\n";
	}
};

Step_S2D_ME_s_up::Step_S2D_ME_s_up() :
	Step(&solve_substep_S2D_ME_s_up),
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

	return 0;
}

int Step_S2D_ME_s_up::finalize_calculation(void) { return 0; }

int solve_substep_S2D_ME_s_up(void *_self)
{
	typedef Step_S2D_ME_s_up::DOF DOF;
	Step_S2D_ME_s_up &self = *(Step_S2D_ME_s_up *)(_self);
	Model_S2D_ME_s_up &model = *self.model;
	double dt = self.dtime;
	double dt2 = dt * dt;
	
	// init nodes
	size_t node_num = model.node_num;
	for (size_t n_id = 0; n_id < node_num; ++n_id)
		model.nodes[n_id].g_id = node_num;
	// init elements
	for (size_t e_id = 0; e_id < model.elem_num; ++e_id)
		model.elems[e_id].pcls = nullptr;
	// init particles
	size_t cal_node_num = 0;
	for (size_t pcl_id = 0; pcl_id < model.pcl_num; ++pcl_id)
	{
		Particle_mpm &pcl = model.pcls[pcl_id];
		if (pcl.pe)
		{
			if (!model.init_pcl_cal_var(pcl))
				continue;
			pcl.vol = pcl.m / pcl.density;
			pcl.pe->add_pcl(pcl);
			Node_mpm &n1 = model.nodes[pcl.n1_id];
			if (n1.g_id == node_num)
			{
				n1.g_id = cal_node_num;
				++cal_node_num;
			}
			Node_mpm &n2 = model.nodes[pcl.n2_id];
			if (n2.g_id == node_num)
			{
				n2.g_id = cal_node_num;
				++cal_node_num;
			}
			Node_mpm &n3 = model.nodes[pcl.n3_id];
			if (n3.g_id == node_num)
			{
				n3.g_id = cal_node_num;
				++cal_node_num;
			}
			Node_mpm &n4 = model.nodes[pcl.n4_id];
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
	Eigen::SparseMatrix<double> g_kmat(dof_num, dof_num);
	Eigen::VectorXd g_fvec = Eigen::VectorXd::Zero(dof_num);

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
			//print_mat(e_kmat);
			//print_vec(e_fvec);

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
			//print_vec(l2g_id_map);
			// add to global matrix and vector
			size_t g_id1, g_id2, g_id;
			for (size_t l_id1 = 0; l_id1 < 12; ++l_id1)
			{
				g_id1 = l2g_id_map[l_id1];
				for (size_t l_id2 = 0; l_id2 < 12; ++l_id2)
				{
					g_id2 = l2g_id_map[l_id2];
					g_kmat_coefs.add_coefficient(g_id1, g_id2, e_kmat[l_id1][l_id2]);
				}
			}
			for (size_t l_id = 0; l_id < 12; ++l_id)
			{
				g_id = l2g_id_map[l_id];
				g_fvec[g_id] += e_fvec[l_id];
			}
		}
	}
	//std::cout << "f:\n" << g_fvec << "\n";

	// apply external force body force and traction
	// traction force
	size_t dof_g_id;
	double tmp1, tmp2;
	double len, xi_g1, xi_g2, eta_g1, eta_g2, t_g1, t_g2;
	for (size_t t_id = 0; t_id < model.tx_num; ++t_id)
	{
		TractionBC_MPM &tx = model.txs[t_id];
		Particle_mpm &pcl = model.pcls[tx.pcl_id];
		Node_mpm &n1 = model.nodes[pcl.n1_id];
		dof_g_id = self.n_id_to_dof_id(n1.g_id, DOF::ux);
		g_fvec[dof_g_id] += pcl.N1 * tx.t;
		Node_mpm &n2 = model.nodes[pcl.n2_id];
		dof_g_id = self.n_id_to_dof_id(n2.g_id, DOF::ux);
		g_fvec[dof_g_id] += pcl.N2 * tx.t;
		Node_mpm &n3 = model.nodes[pcl.n3_id];
		dof_g_id = self.n_id_to_dof_id(n3.g_id, DOF::ux);
		g_fvec[dof_g_id] += pcl.N3 * tx.t;
		Node_mpm &n4 = model.nodes[pcl.n4_id];
		dof_g_id = self.n_id_to_dof_id(n4.g_id, DOF::ux);
		g_fvec[dof_g_id] += pcl.N4 * tx.t;
	}
	for (size_t t_id = 0; t_id < model.ty_num; ++t_id)
	{
		TractionBC_MPM &ty = model.tys[t_id];
		Particle_mpm &pcl = model.pcls[ty.pcl_id];
		Node_mpm &n1 = model.nodes[pcl.n1_id];
		dof_g_id = self.n_id_to_dof_id(n1.g_id, DOF::uy);
		g_fvec[dof_g_id] += pcl.N1 * ty.t;
		Node_mpm &n2 = model.nodes[pcl.n2_id];
		dof_g_id = self.n_id_to_dof_id(n2.g_id, DOF::uy);
		g_fvec[dof_g_id] += pcl.N2 * ty.t;
		Node_mpm &n3 = model.nodes[pcl.n3_id];
		dof_g_id = self.n_id_to_dof_id(n3.g_id, DOF::uy);
		g_fvec[dof_g_id] += pcl.N3 * ty.t;
		Node_mpm &n4 = model.nodes[pcl.n4_id];
		dof_g_id = self.n_id_to_dof_id(n4.g_id, DOF::uy);
		g_fvec[dof_g_id] += pcl.N4 * ty.t;
	}
	// body force to be finished...
	//std::cout << "f:\n" << g_fvec << "\n";

	// apply displacement boundary condition
	double dig_term;
	double *kmat_col = self.kmat_col_mem.resize(dof_num);
	for (size_t bc_id = 0; bc_id < model.ux_num; ++bc_id)
	{
		DisplacementBC &dbc = model.uxs[bc_id];
		Node_mpm &n = model.nodes[dbc.node_id];
		if (n.g_id != node_num) // this node need calculation
		{
			dof_g_id = self.n_id_to_dof_id(n.g_id, DOF::ux);
			//dig_term = g_kmat_coefs.del_row(dof_g_id);
			//g_fvec[dof_g_id] = dig_term * dbc.u;
			dig_term = g_kmat_coefs.del_col_and_row(dof_g_id, kmat_col);
			for (size_t row_id = 0; row_id < dof_num; ++row_id)
				g_fvec[row_id] -= kmat_col[row_id] * dbc.u;
			g_fvec[dof_g_id] = dig_term * dbc.u;
		}
	}
	for (size_t bc_id = 0; bc_id < model.uy_num; ++bc_id)
	{
		DisplacementBC &dbc = model.uys[bc_id];
		Node_mpm &n = model.nodes[dbc.node_id];
		if (n.g_id != node_num)
		{
			dof_g_id = self.n_id_to_dof_id(n.g_id, DOF::uy);
			dig_term = g_kmat_coefs.del_col_and_row(dof_g_id, kmat_col);
			for (size_t row_id = 0; row_id < dof_num; ++row_id)
				g_fvec[row_id] -= kmat_col[row_id] * dbc.u;
			g_fvec[dof_g_id] = dig_term * dbc.u;
		}
	}
	for (size_t bc_id = 0; bc_id < model.pbc_num; ++bc_id)
	{
		PressureBC &pbc = model.pbcs[bc_id];
		Node_mpm &n = model.nodes[pbc.node_id];
		if (n.g_id != node_num)
		{
			dof_g_id = self.n_id_to_dof_id(n.g_id, DOF::p);
			dig_term = g_kmat_coefs.del_col_and_row(dof_g_id, kmat_col);
			for (size_t row_id = 0; row_id < dof_num; ++row_id)
				g_fvec[row_id] -= kmat_col[row_id] * pbc.p;
			g_fvec[dof_g_id] = dig_term * pbc.p;
		}
	}

	g_kmat.setFromTriplets(g_kmat_coefs.begin(), g_kmat_coefs.end());
	// performs a Cholesky factorization
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > chol(g_kmat);
	// solve with the factorization 
	Eigen::VectorXd g_du_vec = chol.solve(g_fvec);

	// Update particle variables
	double dux1, duy1, dp1;
	double dux2, duy2, dp2;
	double dux3, duy3, dp3;
	double dux4, duy4, dp4;
	Particle_mpm *pcl_iter;
	double dux, duy, dvx, dvy, dax, day;
	double de11, de22, de12, dw12;
	double ds11, ds22, ds12;
	double de_vol;
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
			pcl_iter = e.pcls;
			while (pcl_iter)
			{
				// displacement increment
				dux = dux1 * pcl_iter->N1 + dux2 * pcl_iter->N2
					+ dux3 * pcl_iter->N3 + dux4 * pcl_iter->N4;
				duy = duy1 * pcl_iter->N1 + duy2 * pcl_iter->N2
					+ duy3 * pcl_iter->N3 + duy4 * pcl_iter->N4;
				// velocity increment
				dvx = gamma / (beta * dt) * dux
					- gamma / beta * pcl_iter->vx
					+ dt * (1.0 - gamma / (2.0 * beta)) * pcl_iter->ax;
				dvy = gamma / (beta * dt) * duy
					- gamma / beta * pcl_iter->vy
					+ dt * (1.0 - gamma / (2.0 * beta)) * pcl_iter->ay;
				dax = 1.0 / (beta * dt2) * dux
					- 1.0 / (beta * dt) * pcl_iter->vx
					- 1.0 / (2.0 * beta) * pcl_iter->ax;
				day = 1.0 / (beta * dt2) * duy
					- 1.0 / (beta * dt) * pcl_iter->vy
					- 1.0 / (2.0 * beta) * pcl_iter->ay;

				// displacement
				pcl_iter->ux += dux;
				pcl_iter->uy += duy;
				// position
				pcl_iter->x = pcl_iter->x_ori + pcl_iter->ux;
				pcl_iter->y = pcl_iter->y_ori + pcl_iter->uy;
				// velocity
				pcl_iter->vx += dvx;
				pcl_iter->vy += dvy;
				// acceleration
				pcl_iter->ax += dax;
				pcl_iter->ay += day;

				// strain increment
				de11 = dux1 * pcl_iter->dN1_dx + dux2 * pcl_iter->dN2_dx
					 + dux3 * pcl_iter->dN3_dx + dux4 * pcl_iter->dN4_dx;
				de22 = duy1 * pcl_iter->dN1_dy + duy2 * pcl_iter->dN2_dy
					 + duy3 * pcl_iter->dN3_dy + duy4 * pcl_iter->dN4_dy;
				de12 = (dux1 * pcl_iter->dN1_dy + dux2 * pcl_iter->dN2_dy
					  + dux3 * pcl_iter->dN3_dy + dux4 * pcl_iter->dN4_dy
					  + duy1 * pcl_iter->dN1_dx + duy2 * pcl_iter->dN2_dx
					  + duy3 * pcl_iter->dN3_dx + duy4 * pcl_iter->dN4_dx) * 0.5;
				dw12 = (dux1 * pcl_iter->dN1_dy + dux2 * pcl_iter->dN2_dy
					  + dux3 * pcl_iter->dN3_dy + dux4 * pcl_iter->dN4_dy
					  - duy1 * pcl_iter->dN1_dx - duy2 * pcl_iter->dN2_dx
					  - duy3 * pcl_iter->dN3_dx - duy4 * pcl_iter->dN4_dx) * 0.5;
				// update strain (also assume that strain increment is Jaumann rate)
				//de11 +=  dw12 * e12 * 2.0;
				//de22 += -dw12 * e12 * 2.0;
				//de12 +=  dw12 * (e22 - e11);
				pcl_iter->e11 += de11;
				pcl_iter->e22 += de22;
				pcl_iter->e12 += de12;

				// update stress
				double E_tmp = pcl_iter->E / (3.0 * (1.0 + pcl_iter->niu));
				ds11 = 2.0 * E_tmp * de11 - E_tmp * de22;
				ds22 = -E_tmp * de11 + 2.0 * E_tmp * de22;
				ds12 = pcl_iter->E / (1.0 + pcl_iter->niu) * de12;
				/* ------------------------------------------------------------------
				Rotate as Jaumann rate:
				tensor_rate = tensor_Jaumann_rate + tensor * dW_T + dW * tensor
				------------------------------------------------------------------- */
				//ds11 +=  dw12 * pcl_var.s12 * 2.0;
				//ds22 += -dw12 * pcl_var.s12 * 2.0;
				//ds12 +=  dw12 * (pcl_var.s22 - pcl_var.s11);
				pcl_iter->s11 += ds11;
				pcl_iter->s22 += ds22;
				pcl_iter->s12 += ds12;

				pcl_iter->p += dp1 * pcl_iter->N1 + dp2 * pcl_iter->N2
							 + dp3 * pcl_iter->N3 + dp4 * pcl_iter->N4;

				// density
				de_vol = de11 + de22;
				pcl_iter->density /= (1.0 + de_vol);

				// next particle
				pcl_iter = pcl_iter->next;
			}
		}
	}

	return 0;
}

namespace
{
	void cal_stiffness_mat(double k_mat[12][12], double E[3][3], double dN_dx[3][8], double pcl_vol)
	{
		double mid_mat[3][8];
		for (size_t i = 0; i < 3; ++i)
			for (size_t j = 0; j < 8; ++j)
				mid_mat[i][j] = E[i][0] * dN_dx[0][j]
				+ E[i][1] * dN_dx[1][j]
				+ E[i][2] * dN_dx[2][j];
		for (size_t i = 0; i < 8; ++i)
			for (size_t j = 0; j < 8; ++j)
				k_mat[i][j] += (dN_dx[0][i] * mid_mat[0][j]
					+ dN_dx[1][i] * mid_mat[1][j]
					+ dN_dx[2][i] * mid_mat[2][j]) * pcl_vol;
	}
};

void Step_S2D_ME_s_up::form_elem_stiffness_mat_and_force_vec(
	Model_S2D_ME_s_up::Element &e, double kmat[12][12], double fvec[12])
{
	double pcl_vol;
	Particle_mpm *ppcl;
	double mat_term;
	union
	{
		struct
		{
			double m_mat[12][12];
			double k_mat[12][12];
			double f_int[8];
			double f_a[8];
			double f_v[8];
		};
		char tmp_mem[1];
	};

	memset(tmp_mem, 0, sizeof(m_mat) + sizeof(k_mat) + sizeof(f_int) + sizeof(f_a) + sizeof(f_v));
	for (ppcl = e.pcls; ppcl; ppcl = ppcl->next)
	{
		pcl_vol = ppcl->vol;

		// mass matrix
		mat_term = ppcl->N1 * ppcl->m * ppcl->N1;
		m_mat[0][0] += mat_term;
		m_mat[4][4] += mat_term;
		mat_term = ppcl->N1 * ppcl->m * ppcl->N2;
		m_mat[0][1] += mat_term;
		m_mat[1][0] += mat_term;
		m_mat[4][5] += mat_term;
		m_mat[5][4] += mat_term;
		mat_term = ppcl->N1 * ppcl->m * ppcl->N3;
		m_mat[0][2] += mat_term;
		m_mat[2][0] += mat_term;
		m_mat[4][6] += mat_term;
		m_mat[6][4] += mat_term;
		mat_term = ppcl->N1 * ppcl->m * ppcl->N4;
		m_mat[0][3] += mat_term;
		m_mat[3][0] += mat_term;
		m_mat[4][7] += mat_term;
		m_mat[7][4] += mat_term;
		mat_term = ppcl->N2 * ppcl->m * ppcl->N2;
		m_mat[1][1] += mat_term;
		m_mat[5][5] += mat_term;
		mat_term = ppcl->N2 * ppcl->m * ppcl->N3;
		m_mat[1][2] += mat_term;
		m_mat[2][1] += mat_term;
		m_mat[5][6] += mat_term;
		m_mat[6][5] += mat_term;
		mat_term = ppcl->N2 * ppcl->m * ppcl->N4;
		m_mat[1][3] += mat_term;
		m_mat[3][1] += mat_term;
		m_mat[5][7] += mat_term;
		m_mat[7][5] += mat_term;
		mat_term = ppcl->N3 * ppcl->m * ppcl->N3;
		m_mat[2][2] += mat_term;
		m_mat[6][6] += mat_term;
		mat_term = ppcl->N3 * ppcl->m * ppcl->N4;
		m_mat[2][3] += mat_term;
		m_mat[3][2] += mat_term;
		m_mat[6][7] += mat_term;
		m_mat[7][6] += mat_term;
		mat_term = ppcl->N4* ppcl->m * ppcl->N4;
		m_mat[3][3] += mat_term;
		m_mat[7][7] += mat_term;

		// stiffness matrix
		// B * D * B
		double E_mat[3][3], dN_dx_mat[3][8];
		double G = ppcl->E / (2.0 * (1.0 + ppcl->niu));
		E_mat[0][0] = 4.0 / 3.0 * G;
		E_mat[0][1] = -2.0 / 3.0 * G;
		E_mat[0][2] = 0.0;
		E_mat[1][0] = -2.0 / 3.0 * G;
		E_mat[1][1] = 4.0 / 3.0 * G;
		E_mat[1][2] = 0.0;
		E_mat[2][0] = 0.0;
		E_mat[2][1] = 0.0;
		E_mat[2][2] = G;
		dN_dx_mat[0][0] = ppcl->dN1_dx;
		dN_dx_mat[0][1] = 0.0;
		dN_dx_mat[0][2] = ppcl->dN2_dx;
		dN_dx_mat[0][3] = 0.0;
		dN_dx_mat[0][4] = ppcl->dN3_dx;
		dN_dx_mat[0][5] = 0.0;
		dN_dx_mat[0][6] = ppcl->dN4_dx;
		dN_dx_mat[0][7] = 0.0;
		dN_dx_mat[1][0] = 0.0;
		dN_dx_mat[1][1] = ppcl->dN1_dy;
		dN_dx_mat[1][2] = 0.0;
		dN_dx_mat[1][3] = ppcl->dN2_dy;
		dN_dx_mat[1][4] = 0.0;
		dN_dx_mat[1][5] = ppcl->dN3_dy;
		dN_dx_mat[1][6] = 0.0;
		dN_dx_mat[1][7] = ppcl->dN4_dy;
		dN_dx_mat[2][0] = ppcl->dN1_dy;
		dN_dx_mat[2][1] = ppcl->dN1_dx;
		dN_dx_mat[2][2] = ppcl->dN2_dy;
		dN_dx_mat[2][3] = ppcl->dN2_dx;
		dN_dx_mat[2][4] = ppcl->dN3_dy;
		dN_dx_mat[2][5] = ppcl->dN3_dx;
		dN_dx_mat[2][6] = ppcl->dN4_dy;
		dN_dx_mat[2][7] = ppcl->dN4_dx;
		cal_stiffness_mat(k_mat, E_mat, dN_dx_mat, pcl_vol);
		// dNi_dx * Nj
		mat_term = ppcl->dN1_dx * ppcl->N1 * pcl_vol;
		k_mat[0][8] += mat_term;
		mat_term = ppcl->dN1_dx * ppcl->N2 * pcl_vol;
		k_mat[0][9] += mat_term;
		k_mat[1][8] += mat_term;
		mat_term = ppcl->dN1_dx * ppcl->N3 * pcl_vol;
		k_mat[0][10] += mat_term;
		k_mat[2][8] += mat_term;
		mat_term = ppcl->dN1_dx * ppcl->N4 * pcl_vol;
		k_mat[0][11] += mat_term;
		k_mat[3][8] += mat_term;
		mat_term = ppcl->dN2_dx * ppcl->N2 * pcl_vol;
		k_mat[1][9] += mat_term;
		mat_term = ppcl->dN2_dx * ppcl->N3 * pcl_vol;
		k_mat[1][10] += mat_term;
		k_mat[2][9] += mat_term;
		mat_term = ppcl->dN2_dx * ppcl->N4 * pcl_vol;
		k_mat[1][11] += mat_term;
		k_mat[3][9] += mat_term;
		mat_term = ppcl->dN3_dx * ppcl->N3 * pcl_vol;
		k_mat[2][10] += mat_term;
		mat_term = ppcl->dN3_dx * ppcl->N4 * pcl_vol;
		k_mat[2][11] += mat_term;
		k_mat[3][10] += mat_term;
		mat_term = ppcl->dN4_dx * ppcl->N4 * pcl_vol;
		k_mat[3][11] += mat_term;
		// dNi_dy * Nj
		mat_term = ppcl->dN1_dy * ppcl->N1 * pcl_vol;
		k_mat[4][8] += mat_term;
		mat_term = ppcl->dN1_dy * ppcl->N2 * pcl_vol;
		k_mat[4][9] += mat_term;
		k_mat[5][8] += mat_term;
		mat_term = ppcl->dN1_dy * ppcl->N3 * pcl_vol;
		k_mat[4][10] += mat_term;
		k_mat[6][8] += mat_term;
		mat_term = ppcl->dN1_dy * ppcl->N4 * pcl_vol;
		k_mat[4][11] += mat_term;
		k_mat[7][8] += mat_term;
		mat_term = ppcl->dN2_dy * ppcl->N2 * pcl_vol;
		k_mat[5][9] += mat_term;
		mat_term = ppcl->dN2_dy * ppcl->N3 * pcl_vol;
		k_mat[5][10] += mat_term;
		k_mat[6][9] += mat_term;
		mat_term = ppcl->dN2_dy * ppcl->N4 * pcl_vol;
		k_mat[5][11] += mat_term;
		k_mat[7][9] += mat_term;
		mat_term = ppcl->dN3_dy * ppcl->N3 * pcl_vol;
		k_mat[6][10] += mat_term;
		mat_term = ppcl->dN3_dy * ppcl->N4 * pcl_vol;
		k_mat[6][11] += mat_term;
		k_mat[7][10] += mat_term;
		mat_term = ppcl->dN4_dy * ppcl->N4 * pcl_vol;
		k_mat[7][11] += mat_term;
		// Ni * K * dNj_dx
		mat_term = ppcl->N1 * ppcl->K * ppcl->dN1_dx * pcl_vol;
		k_mat[8][0] += mat_term;
		mat_term = ppcl->N1 * ppcl->K * ppcl->dN2_dx * pcl_vol;
		k_mat[8][1] += mat_term;
		k_mat[9][0] += mat_term;
		mat_term = ppcl->N1 * ppcl->K * ppcl->dN3_dx * pcl_vol;
		k_mat[8][2] += mat_term;
		k_mat[10][0] += mat_term;
		mat_term = ppcl->N1 * ppcl->K * ppcl->dN4_dx * pcl_vol;
		k_mat[8][3] += mat_term;
		k_mat[11][0] += mat_term;
		mat_term = ppcl->N2 * ppcl->K * ppcl->dN2_dx * pcl_vol;
		k_mat[9][1] += mat_term;
		mat_term = ppcl->N2 * ppcl->K * ppcl->dN3_dx * pcl_vol;
		k_mat[9][2] += mat_term;
		k_mat[10][1] += mat_term;
		mat_term = ppcl->N2 * ppcl->K * ppcl->dN4_dx * pcl_vol;
		k_mat[9][3] += mat_term;
		k_mat[11][1] += mat_term;
		mat_term = ppcl->N3 * ppcl->K * ppcl->dN3_dx * pcl_vol;
		k_mat[10][2] += mat_term;
		mat_term = ppcl->N3 * ppcl->K * ppcl->dN4_dx * pcl_vol;
		k_mat[10][3] += mat_term;
		k_mat[11][2] += mat_term;
		mat_term = ppcl->N4 * ppcl->K * ppcl->dN4_dx * pcl_vol;
		k_mat[11][3] += mat_term;
		// Ni * K * dNj_dy
		mat_term = ppcl->N1 * ppcl->K * ppcl->dN1_dy * pcl_vol;
		k_mat[8][4] += mat_term;
		mat_term = ppcl->N1 * ppcl->K * ppcl->dN2_dy * pcl_vol;
		k_mat[8][5] += mat_term;
		k_mat[9][4] += mat_term;
		mat_term = ppcl->N1 * ppcl->K * ppcl->dN3_dy * pcl_vol;
		k_mat[8][6] += mat_term;
		k_mat[10][4] += mat_term;
		mat_term = ppcl->N1 * ppcl->K * ppcl->dN4_dy * pcl_vol;
		k_mat[8][7] += mat_term;
		k_mat[11][4] += mat_term;
		mat_term = ppcl->N2 * ppcl->K * ppcl->dN2_dy * pcl_vol;
		k_mat[9][5] += mat_term;
		mat_term = ppcl->N2 * ppcl->K * ppcl->dN3_dy * pcl_vol;
		k_mat[9][6] += mat_term;
		k_mat[10][5] += mat_term;
		mat_term = ppcl->N2 * ppcl->K * ppcl->dN4_dy * pcl_vol;
		k_mat[9][7] += mat_term;
		k_mat[11][5] += mat_term;
		mat_term = ppcl->N3 * ppcl->K * ppcl->dN3_dy * pcl_vol;
		k_mat[10][6] += mat_term;
		mat_term = ppcl->N3 * ppcl->K * ppcl->dN4_dy * pcl_vol;
		k_mat[10][7] += mat_term;
		k_mat[11][6] += mat_term;
		mat_term = ppcl->N4 * ppcl->K * ppcl->dN4_dy * pcl_vol;
		k_mat[11][7] += mat_term;
		// Ni * Nj
		mat_term = -ppcl->N1 * ppcl->N1 * pcl_vol;
		k_mat[8][8] += mat_term;
		mat_term = -ppcl->N1 * ppcl->N2 * pcl_vol;
		k_mat[8][9] += mat_term;
		k_mat[9][8] += mat_term;
		mat_term = -ppcl->N1 * ppcl->N3 * pcl_vol;
		k_mat[8][10] += mat_term;
		k_mat[10][8] += mat_term;
		mat_term = -ppcl->N1 * ppcl->N4 * pcl_vol;
		k_mat[8][11] += mat_term;
		k_mat[11][8] += mat_term;
		mat_term = -ppcl->N2 * ppcl->N2 * pcl_vol;
		k_mat[9][9] += mat_term;
		mat_term = -ppcl->N2 * ppcl->N3 * pcl_vol;
		k_mat[9][10] += mat_term;
		k_mat[10][9] += mat_term;
		mat_term = -ppcl->N2 * ppcl->N4 * pcl_vol;
		k_mat[9][11] += mat_term;
		k_mat[11][9] += mat_term;
		mat_term = -ppcl->N3 * ppcl->N3 * pcl_vol;
		k_mat[10][10] += mat_term;
		mat_term = -ppcl->N3 * ppcl->N4 * pcl_vol;
		k_mat[10][11] += mat_term;
		k_mat[11][10] += mat_term;
		mat_term = -ppcl->N4 * ppcl->N4 * pcl_vol;
		k_mat[11][11] += mat_term;

		// force vector
		f_int[0] += (ppcl->dN1_dx * (ppcl->s11 + ppcl->p) + ppcl->dN1_dy * ppcl->s12) * pcl_vol;
		f_int[1] += (ppcl->dN2_dx * (ppcl->s11 + ppcl->p) + ppcl->dN2_dy * ppcl->s12) * pcl_vol;
		f_int[2] += (ppcl->dN3_dx * (ppcl->s11 + ppcl->p) + ppcl->dN3_dy * ppcl->s12) * pcl_vol;
		f_int[3] += (ppcl->dN4_dx * (ppcl->s11 + ppcl->p) + ppcl->dN4_dy * ppcl->s12) * pcl_vol;
		f_int[4] += (ppcl->dN1_dx * ppcl->s12 + ppcl->dN1_dy * (ppcl->s22 + ppcl->p)) * pcl_vol;
		f_int[5] += (ppcl->dN2_dx * ppcl->s12 + ppcl->dN2_dy * (ppcl->s22 + ppcl->p)) * pcl_vol;
		f_int[6] += (ppcl->dN3_dx * ppcl->s12 + ppcl->dN3_dy * (ppcl->s22 + ppcl->p)) * pcl_vol;
		f_int[7] += (ppcl->dN4_dx * ppcl->s12 + ppcl->dN4_dy * (ppcl->s22 + ppcl->p)) * pcl_vol;

		f_a[0] += ppcl->N1 * ppcl->m * ppcl->ax;
		f_a[1] += ppcl->N2 * ppcl->m * ppcl->ax;
		f_a[2] += ppcl->N3 * ppcl->m * ppcl->ax;
		f_a[3] += ppcl->N4 * ppcl->m * ppcl->ax;
		f_a[4] += ppcl->N1 * ppcl->m * ppcl->ay;
		f_a[5] += ppcl->N2 * ppcl->m * ppcl->ay;
		f_a[6] += ppcl->N3 * ppcl->m * ppcl->ay;
		f_a[7] += ppcl->N4 * ppcl->m * ppcl->ay;

		f_v[0] += ppcl->N1 * ppcl->m * ppcl->vx;
		f_v[1] += ppcl->N2 * ppcl->m * ppcl->vx;
		f_v[2] += ppcl->N3 * ppcl->m * ppcl->vx;
		f_v[3] += ppcl->N4 * ppcl->m * ppcl->vx;
		f_v[4] += ppcl->N1 * ppcl->m * ppcl->vy;
		f_v[5] += ppcl->N2 * ppcl->m * ppcl->vy;
		f_v[6] += ppcl->N3 * ppcl->m * ppcl->vy;
		f_v[7] += ppcl->N4 * ppcl->m * ppcl->vy;
	}

	for (size_t i = 0; i < 12; ++i)
		for (size_t j = 0; j < 12; ++j)
			kmat[i][j] = m_mat[i][j] / (beta * dtime * dtime) + k_mat[i][j];

	for (size_t i = 0; i < 8; ++i)
		fvec[i] = -f_int[i] + (1.0 / (2.0*beta) - 1.0) * f_a[i] + 1.0 / (beta*dtime) * f_v[i];
	fvec[8] = 0.0;
	fvec[9] = 0.0;
	fvec[10] = 0.0;
	fvec[11] = 0.0;
}
