#include "SimulationCore_pcp.h"

#include <cmath>
#include <vector>

#include <Eigen/Sparse>
#include "Step_S2D_CHM_s_FEM_uUp.h"

// constant of Newmark-beta method
#define beta 0.5
#define gamma 0.25

using namespace Model_S2D_CHM_s_FEM_uUp_Internal;

namespace
{
	typedef Model_S2D_CHM_s_FEM_uUp::Element Element_fem;
	typedef Model_S2D_CHM_s_FEM_uUp::Node Node_fem;
};

void cal_stiffness_mat(double kmat[20][20], double E[3][3], double dN_dx[3][8])
{
	double mid_mat[3][8];
	for (size_t i = 0; i < 3; ++i)
		for (size_t j = 0; j < 8; ++j)
			mid_mat[i][j] = E[i][0] * dN_dx[0][j]
						  + E[i][1] * dN_dx[1][j]
						  + E[i][2] * dN_dx[2][j];
	for (size_t i = 0; i < 8; ++i)
		for (size_t j = 0; j < 8; ++j)
			kmat[i][j] = dN_dx[0][i] * mid_mat[0][j]
					   + dN_dx[1][i] * mid_mat[1][j]
					   + dN_dx[2][i] * mid_mat[2][j];
}

Step_S2D_CHM_s_FEM_uUp::Step_S2D_CHM_s_FEM_uUp() :
	Step(&solve_substep_S2D_CHM_s_FEM_uUp),
	model(nullptr) {}

Step_S2D_CHM_s_FEM_uUp::~Step_S2D_CHM_s_FEM_uUp() {}

int Step_S2D_CHM_s_FEM_uUp::init_calculation(void)
{
	if (is_first_step)
	{

	}

	return 0;
}

int Step_S2D_CHM_s_FEM_uUp::finalize_calculation(void) { return 0; }

int solve_substep_S2D_CHM_s_FEM_uUp(void *_self)
{
	typedef Step_S2D_CHM_s_FEM_uUp::DOF DOF;
	Step_S2D_CHM_s_FEM_uUp &self = *(Step_S2D_CHM_s_FEM_uUp *)(_self);
	Model_S2D_CHM_s_FEM_uUp &model = *self.model;

	double dt = self.dtime;
	double dt2 = dt * dt;

	// init nodes
	for (size_t n_id = 0; n_id < model.node_num; ++n_id)
	{
		Node_fem &n = model.nodes[n_id];
		// solid phase
		n.ax_s = 0.0;
		n.ay_s = 0.0;
		n.vx_s = 0.0;
		n.vy_s = 0.0;
		n.ux_s = 0.0;
		n.uy_s = 0.0;
		// fluid phase
		n.ax_f = 0.0;
		n.ay_f = 0.0;
		n.vx_f = 0.0;
		n.vy_f = 0.0;
		n.ux_f = 0.0;
		n.uy_f = 0.0;
	}

	std::vector<Eigen::Triplet<double> > g_kmat_coefs; // list of non-zeros coefficients
	size_t dof_num = model.node_num * 5;
	Eigen::SparseMatrix<double> g_kmat(dof_num, dof_num);
	Eigen::VectorXd g_fvec(dof_num);
	size_t n_ids[4], l2g_id_map[20];
	double e_kmat[20][20], e_fvec[20];
	for (size_t e_id = 0; e_id < model.elem_num; ++e_id)
	{
		Element_fem &e = model.elems[e_id];
		
		// form elemental stiffness matrix and force vector
		model.cal_shape_func_value(e.sf_var, 0.0, 0.0);
		self.form_elem_stiffness_mat(e, e_kmat);
		self.form_elem_force_vec(e, dt, e_fvec);

		// form map from global id to local id
		n_ids[0] = e.index_y * model.node_x_num + e.index_x;
		n_ids[1] = n_ids[0] + 1;
		n_ids[2] = n_ids[1] + model.node_x_num;
		n_ids[3] = n_ids[2] - 1;
		size_t k = 0;
		for (size_t dof_id = 0; dof_id < 5; ++dof_id)
			for (size_t l_id = 0; l_id < 4; ++l_id)
				l2g_id_map[k++] = self.n_id_to_dof_id(n_ids[l_id], DOF(dof_id));
		// add to global matrix and vector
		size_t g_id1, g_id2, g_id;
		for (size_t l_id1 = 0; l_id1 < 20; ++l_id1)
		{
			g_id1 = l2g_id_map[l_id1];
			for (size_t l_id2 = 0; l_id2 < 20; ++l_id2)
			{
				g_id2 = l2g_id_map[l_id2];
				Eigen::Triplet<double> coefficent(g_id1, g_id2, e_kmat[l_id1][l_id2]);
				g_kmat_coefs.push_back(coefficent);
			}
		}
		for (size_t l_id = 0; l_id < 20; ++l_id)
		{
			g_id = l2g_id_map[l_id];
			g_fvec[g_id] += e_fvec[l_id];
		}
	}

	// apply external force body force and traction
	// traction force
#define p_g1 -0.57735
#define p_g2  0.57735
#define w_g1 1.0
#define w_g2 1.0
	for (size_t t_id = 0; t_id < model.tx_num; ++t_id)
	{
		TractionBC_2DFEM &tx = model.txs[t_id];
		Element_fem &e = model.elems[tx.elem_id];
		double xi_g1, xi_g2, eta_g1, eta_g2, t_g1, t_g2;
		double tmp1, tmp2;
		tmp1 = tx.xi1 - tx.xi0;
		tmp2 = tx.eta1 - tx.eta0;
		double len = model.h * 0.25 * sqrt(tmp1 * tmp1 + tmp2 * tmp2);
		tmp1 = 0.5 * (tx.xi0 + tx.xi1);
		tmp2 = 0.5 * (tx.xi1 - tx.xi0);
		xi_g1 = tmp1 + p_g1 * tmp2;
		xi_g2 = tmp1 + p_g2 * tmp2;
		tmp1 = 0.5 * (tx.eta0 + tx.eta1);
		tmp2 = 0.5 * (tx.eta1 - tx.eta0);
		eta_g1 = tmp1 + p_g1 * tmp2;
		eta_g2 = tmp1 + p_g2 * tmp2;
		tmp1 = 0.5 * (tx.t0 + tx.t1);
		tmp2 = 0.5 * (tx.t1 - tx.t0);
		t_g1 = tmp1 + p_g1 * tmp2;
		t_g2 = tmp1 + p_g2 * tmp2;
		size_t n1_id, n2_id, n3_id, n4_id;
		model.get_node_id(e, n1_id, n2_id, n3_id, n4_id);
		g_fvec[self.n_id_to_dof_id(n1_id, DOF::usx)] += (N1(xi_g1, eta_g1) * t_g1 * w_g1 + N1(xi_g2, eta_g2) * t_g2 * w_g2) * len;
		g_fvec[self.n_id_to_dof_id(n2_id, DOF::usx)] += (N2(xi_g1, eta_g1) * t_g1 * w_g1 + N2(xi_g2, eta_g2) * t_g2 * w_g2) * len;
		g_fvec[self.n_id_to_dof_id(n3_id, DOF::usx)] += (N3(xi_g1, eta_g1) * t_g1 * w_g1 + N3(xi_g2, eta_g2) * t_g2 * w_g2) * len;
		g_fvec[self.n_id_to_dof_id(n4_id, DOF::usx)] += (N4(xi_g1, eta_g1) * t_g1 * w_g1 + N4(xi_g2, eta_g2) * t_g2 * w_g2) * len;
	}
	for (size_t t_id = 0; t_id < model.ty_num; ++t_id)
	{
		TractionBC_2DFEM &ty = model.tys[t_id];
		Element_fem &e = model.elems[ty.elem_id];
		double xi_g1, xi_g2, eta_g1, eta_g2, t_g1, t_g2;
		double tmp1, tmp2;
		tmp1 = ty.xi1 - ty.xi0;
		tmp2 = ty.eta1 - ty.eta0;
		double len = model.h * 0.25 * sqrt(tmp1 * tmp1 + tmp2 * tmp2);
		tmp1 = 0.5 * (ty.xi0 + ty.xi1);
		tmp2 = 0.5 * (ty.xi1 - ty.xi0);
		xi_g1 = tmp1 + p_g1 * tmp2;
		xi_g2 = tmp1 + p_g2 * tmp2;
		tmp1 = 0.5 * (ty.eta0 + ty.eta1);
		tmp2 = 0.5 * (ty.eta1 - ty.eta0);
		eta_g1 = tmp1 + p_g1 * tmp2;
		eta_g2 = tmp1 + p_g2 * tmp2;
		tmp1 = 0.5 * (ty.t0 + ty.t1);
		tmp2 = 0.5 * (ty.t1 - ty.t0);
		t_g1 = tmp1 + p_g1 * tmp2;
		t_g2 = tmp1 + p_g2 * tmp2;
		size_t n1_id, n2_id, n3_id, n4_id;
		model.get_node_id(e, n1_id, n2_id, n3_id, n4_id);
		g_fvec[self.n_id_to_dof_id(n1_id, DOF::usy)] += (N1(xi_g1, eta_g1) * t_g1 * w_g1 + N1(xi_g2, eta_g2) * t_g2 * w_g2) * len;
		g_fvec[self.n_id_to_dof_id(n2_id, DOF::usy)] += (N2(xi_g1, eta_g1) * t_g1 * w_g1 + N2(xi_g2, eta_g2) * t_g2 * w_g2) * len;
		g_fvec[self.n_id_to_dof_id(n3_id, DOF::usy)] += (N3(xi_g1, eta_g1) * t_g1 * w_g1 + N3(xi_g2, eta_g2) * t_g2 * w_g2) * len;
		g_fvec[self.n_id_to_dof_id(n4_id, DOF::usy)] += (N4(xi_g1, eta_g1) * t_g1 * w_g1 + N4(xi_g2, eta_g2) * t_g2 * w_g2) * len;
	}
	// body force ...  to be finished
	
	// apply displacement boundary condition
	size_t dof_g_id;
	for (size_t bc_id = 0; bc_id < model.usx_num; ++bc_id)
	{
		dof_g_id = self.n_id_to_dof_id(model.usxs[bc_id].node_id, DOF::usx);
	}
	for (size_t bc_id = 0; bc_id < model.usy_num; ++bc_id)
	{
		dof_g_id = self.n_id_to_dof_id(model.usys[bc_id].node_id, DOF::usx);
	}
	for (size_t bc_id = 0; bc_id < model.ufx_num; ++bc_id)
	{
		dof_g_id = self.n_id_to_dof_id(model.ufxs[bc_id].node_id, DOF::usx);
	}
	for (size_t bc_id = 0; bc_id < model.ufy_num; ++bc_id)
	{
		dof_g_id = self.n_id_to_dof_id(model.ufys[bc_id].node_id, DOF::usx);
	}
	for (size_t bc_id = 0; bc_id < model.pbc_num; ++bc_id)
	{
		dof_g_id = self.n_id_to_dof_id(model.pbcs[bc_id].node_id, DOF::usx);
	}

	// solve
	g_kmat.setFromTriplets(g_kmat_coefs.begin(), g_kmat_coefs.end());
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > chol(g_kmat);  // performs a Cholesky factorization of A
	Eigen::VectorXd g_du_vec = chol.solve(g_fvec);         // use the factorization to solve for the given right hand side

	// update nodal variables
	size_t dof_id;
	double da, dv, du;
	for (size_t n_id = 0; n_id < model.node_num; ++n_id)
	{
		Node_fem &n = model.nodes[n_id];
		// sx
		dof_id = self.n_id_to_dof_id(n_id, DOF::usx);
		du = g_du_vec[dof_id];
		n.ux_s += du;
		da = 1.0 / (gamma * dt2) * du - 1.0 / (gamma * dt) * n.vx_s - 1.0 / (2.0 * gamma) * n.ax_s;
		n.ax_s += da;
		dv = gamma / (beta * dt) * du - gamma / beta * n.vx_s + dt * (1.0 - gamma / (2.0 * beta)) * n.ax_s;
		n.vx_s += dv;
		// sy
		dof_id = self.n_id_to_dof_id(n_id, DOF::usy);
		du = g_du_vec[dof_id];
		n.uy_s += du;
		da = 1.0 / (gamma * dt2) * du - 1.0 / (gamma * dt) * n.vy_s - 1.0 / (2.0 * gamma) * n.ay_s;
		n.ay_s += da;
		dv = gamma / (beta * dt) * du - gamma / beta * n.vy_s + dt * (1.0 - gamma / (2.0 * beta)) * n.ay_s;
		n.vy_s += dv;
		// fx
		dof_id = self.n_id_to_dof_id(n_id, DOF::ufx);
		du = g_du_vec[dof_id];
		n.ux_f += du;
		da = 1.0 / (gamma * dt2) * du - 1.0 / (gamma * dt) * n.vx_f - 1.0 / (2.0 * gamma) * n.ax_f;
		n.ax_f += da;
		dv = gamma / (beta * dt) * du - gamma / beta * n.vx_f + dt * (1.0 - gamma / (2.0 * beta)) * n.ax_f;
		n.vx_f += dv;
		// fy
		dof_id = self.n_id_to_dof_id(n_id, DOF::ufy);
		du = g_du_vec[dof_id];
		n.uy_f += du;
		da = 1.0 / (gamma * dt2) * du - 1.0 / (gamma * dt) * n.vy_f - 1.0 / (2.0 * gamma) * n.ay_s;
		n.ay_f += da;
		dv = gamma / (beta * dt) * du - gamma / beta * n.vy_f + dt * (1.0 - gamma / (2.0 * beta)) * n.ay_s;
		n.vy_f += dv;
		// p
		dof_id = self.n_id_to_dof_id(n_id, DOF::p);
		n.p += g_du_vec[dof_id];
	}

	// Update element variables
	double dux_s1, duy_s1, dux_f1, duy_f1;
	double dux_s2, duy_s2, dux_f2, duy_f2;
	double dux_s3, duy_s3, dux_f3, duy_f3;
	double dux_s4, duy_s4, dux_f4, duy_f4;
	size_t n_id1, n_id2, n_id3, n_id4;
	double de11, de22, de12, dw12;
	double ds11, ds22, ds12;
	double de_vol_s, de_vol_f;
	for (size_t e_id = 0; e_id < model.elem_num; ++e_id)
	{
		Element_fem &e = model.elems[e_id];
		n_id1 = model.node_x_num * e.index_y + e.index_x;
		n_id2 = n_id1 + 1;
		n_id3 = n_id2 + model.node_x_num;
		n_id4 = n_id3 - 1;
		dux_s1 = g_du_vec[self.n_id_to_dof_id(n_id1, DOF::usx)];
		duy_s1 = g_du_vec[self.n_id_to_dof_id(n_id1, DOF::usy)];
		dux_f1 = g_du_vec[self.n_id_to_dof_id(n_id1, DOF::ufx)];
		duy_f1 = g_du_vec[self.n_id_to_dof_id(n_id1, DOF::ufy)];
		dux_s2 = g_du_vec[self.n_id_to_dof_id(n_id2, DOF::usx)];
		duy_s2 = g_du_vec[self.n_id_to_dof_id(n_id2, DOF::usy)];
		dux_f2 = g_du_vec[self.n_id_to_dof_id(n_id2, DOF::ufx)];
		duy_f2 = g_du_vec[self.n_id_to_dof_id(n_id2, DOF::ufy)];
		dux_s3 = g_du_vec[self.n_id_to_dof_id(n_id3, DOF::usx)];
		duy_s3 = g_du_vec[self.n_id_to_dof_id(n_id3, DOF::usy)];
		dux_f3 = g_du_vec[self.n_id_to_dof_id(n_id3, DOF::ufx)];
		duy_f3 = g_du_vec[self.n_id_to_dof_id(n_id3, DOF::ufy)];
		dux_s4 = g_du_vec[self.n_id_to_dof_id(n_id4, DOF::usx)];
		duy_s4 = g_du_vec[self.n_id_to_dof_id(n_id4, DOF::usy)];
		dux_f4 = g_du_vec[self.n_id_to_dof_id(n_id4, DOF::ufx)];
		duy_f4 = g_du_vec[self.n_id_to_dof_id(n_id4, DOF::ufy)];

		// strain increment
		de11 = dux_s1 * e.dN1_dx + dux_s2 * e.dN2_dx
			 + dux_s3 * e.dN3_dx + dux_s4 * e.dN4_dx;
		de22 = duy_s1 * e.dN1_dy + duy_s2 * e.dN2_dy
			 + duy_s3 * e.dN3_dy + duy_s4 * e.dN4_dy;
		de12 = (dux_s1 * e.dN1_dy + dux_s2 * e.dN2_dy
			  + dux_s3 * e.dN3_dy + dux_s4 * e.dN4_dy
			  + duy_s1 * e.dN1_dx + duy_s2 * e.dN2_dx
			  + duy_s3 * e.dN3_dx + duy_s4 * e.dN4_dx) * 0.5;
		dw12 = (dux_s1 * e.dN1_dy + dux_s2 * e.dN2_dy
			  + dux_s3 * e.dN3_dy + dux_s4 * e.dN4_dy
			  - duy_s1 * e.dN1_dx + duy_s2 * e.dN2_dx
			  - duy_s3 * e.dN3_dx + duy_s4 * e.dN4_dx) * 0.5;
		// update strain (also assume that strain increment is Jaumann rate)
		//de11 +=  dw12 * e12 * 2.0;
		//de22 += -dw12 * e12 * 2.0;
		//de12 +=  dw12 * (e22 - e11);
		e.e11 += de11;
		e.e22 += de22;
		e.e12 += de12;

		// update stress
		double E_tmp = e.E / (1.0 + e.niu) / (1.0 - 2.0 * e.niu);
		ds11 = E_tmp * ((1.0 - e.niu) * de11 + e.niu * de22);
		ds22 = E_tmp * (e.niu * de11 + (1.0 - e.niu) * de22);
		ds12 = 2.0 * e.E / (2.0 * (1.0 + e.niu)) * de12;
		/* ------------------------------------------------------------------
		Rotate as Jaumann rate:
		tensor_rate = tensor_Jaumann_rate + tensor * dW_T + dW * tensor
		------------------------------------------------------------------- */
		//ds11 +=  dw12 * pcl_var.s12 * 2.0;
		//ds22 += -dw12 * pcl_var.s12 * 2.0;
		//ds12 +=  dw12 * (pcl_var.s22 - pcl_var.s11);
		e.s11 += ds11;
		e.s22 += ds22;
		e.s12 += ds12;

		// volumetric strain of solid phase
		de_vol_s = de11 + de22;
		// "volumetric strain" of fluid phase
		de_vol_f = -(1.0 - e.n) / e.n * de_vol_s
			- (dux_f1 * e.dN1_dx + dux_f2 * e.dN2_dx + dux_f3 * e.dN3_dx + dux_f4 * e.dN4_dx)
			- (duy_f1 * e.dN1_dy + duy_f2 * e.dN2_dy + duy_f3 * e.dN3_dy + duy_f4 * e.dN4_dy);
		// porosity
		e.n = (de_vol_s + e.n) / (1.0 + de_vol_s);

		// fluid density
		e.density_f += e.density_f * de_vol_f;
		// pore pressure (can use EOS and get p from density_f instead)
		e.p += e.Kf * de_vol_f; // ???? 
	}


	return 0;
}

void Step_S2D_CHM_s_FEM_uUp::form_elem_stiffness_mat(Model_S2D_CHM_s_FEM_uUp::Element &e, double kmat[20][20])
{
	double elem_vol = model->elem_vol;
	double mat_term;
	memset(kmat, 0, sizeof(double) * 20 * 20);
	
	// mass matrix
	// solid
	double e_m_s = (1.0 - e.n) * e.density_s * elem_vol;
	mat_term = e.N1 * e.N1 * e_m_s;
	kmat[0][0] = mat_term;
	kmat[4][4] = mat_term;
	mat_term = e.N1 * e.N2 * e_m_s;
	kmat[0][1] = mat_term;
	kmat[1][0] = mat_term;
	kmat[4][5] = mat_term;
	kmat[5][4] = mat_term;
	mat_term = e.N1 * e.N3 * e_m_s;
	kmat[0][2] = mat_term;
	kmat[2][0] = mat_term;
	kmat[4][6] = mat_term;
	kmat[6][4] = mat_term;
	mat_term = e.N1 * e.N4 * e_m_s;
	kmat[0][3] = mat_term;
	kmat[3][0] = mat_term;
	kmat[4][7] = mat_term;
	kmat[7][4] = mat_term;
	mat_term = e.N2 * e.N2 * e_m_s;
	kmat[1][1] = mat_term;
	kmat[5][5] = mat_term;
	mat_term = e.N2 * e.N3 * e_m_s;
	kmat[1][2] = mat_term;
	kmat[2][1] = mat_term;
	kmat[5][6] = mat_term;
	kmat[6][5] = mat_term;
	mat_term = e.N2 * e.N4 * e_m_s;
	kmat[1][3] = mat_term;
	kmat[3][1] = mat_term;
	kmat[5][7] = mat_term;
	kmat[7][5] = mat_term;
	mat_term = e.N3 * e.N3 * e_m_s;
	kmat[2][2] = mat_term;
	kmat[6][6] = mat_term;
	mat_term = e.N3 * e.N4 * e_m_s;
	kmat[2][3] = mat_term;
	kmat[3][2] = mat_term;
	kmat[6][7] = mat_term;
	kmat[7][6] = mat_term;
	mat_term = e.N4 * e.N4 * e_m_s;
	kmat[3][3] = mat_term;
	kmat[7][7] = mat_term;
	// fluid
	double e_m_f = e.n * e.density_f * elem_vol;
	mat_term = e.N1 * e.N1 * e_m_f;
	kmat[8][8] = mat_term;
	kmat[12][12] = mat_term;
	mat_term = e.N1 * e.N2 * e_m_f;
	kmat[8][9] = mat_term;
	kmat[9][8] = mat_term;
	kmat[12][13] = mat_term;
	kmat[13][12] = mat_term;
	mat_term = e.N1 * e.N3 * e_m_f;
	kmat[8][10] = mat_term;
	kmat[10][8] = mat_term;
	kmat[12][14] = mat_term;
	kmat[14][12] = mat_term;
	mat_term = e.N1 * e.N4 * e_m_f;
	kmat[8][11] = mat_term;
	kmat[11][8] = mat_term;
	kmat[12][15] = mat_term;
	kmat[15][12] = mat_term;
	mat_term = e.N2 * e.N2 * e_m_f;
	kmat[9][9] = mat_term;
	kmat[13][13] = mat_term;
	mat_term = e.N2 * e.N3 * e_m_f;
	kmat[9][10] = mat_term;
	kmat[10][9] = mat_term;
	kmat[13][14] = mat_term;
	kmat[14][13] = mat_term;
	mat_term = e.N2 * e.N4 * e_m_f;
	kmat[9][11] = mat_term;
	kmat[11][9] = mat_term;
	kmat[13][15] = mat_term;
	kmat[15][13] = mat_term;
	mat_term = e.N3 * e.N3 * e_m_f;
	kmat[10][10] = mat_term;
	kmat[14][14] = mat_term;
	mat_term = e.N3 * e.N4 * e_m_f;
	kmat[10][11] = mat_term;
	kmat[11][10] = mat_term;
	kmat[14][15] = mat_term;
	kmat[15][14] = mat_term;
	mat_term = e.N4 * e.N4 * e_m_f;
	kmat[11][11] = mat_term;
	kmat[15][15] = mat_term;

	// damping matrix (seepage)
	double seepage_force = e.n * e.n * e.miu / e.k * elem_vol;
#define FILL_SEEP_MAT(x_id, y_id) \
	kmat[x_id][y_id] += mat_term; \
	kmat[x_id+4][y_id+4] += mat_term; \
	kmat[x_id+8][y_id+8] += mat_term; \
	kmat[x_id+12][y_id+12] += mat_term; \
	kmat[x_id+8][y_id] -= mat_term; \
	kmat[x_id+12][y_id+4] -= mat_term; \
	kmat[x_id][y_id+8] -= mat_term; \
	kmat[x_id+4][y_id+12] -= mat_term
	mat_term = e.N1 * e.N1 * seepage_force;
	FILL_SEEP_MAT(0, 0);
	mat_term = e.N1 * e.N2 * seepage_force;
	FILL_SEEP_MAT(0, 1);
	mat_term = e.N1 * e.N3 * seepage_force;
	FILL_SEEP_MAT(0, 2);
	mat_term = e.N1 * e.N4 * seepage_force;
	FILL_SEEP_MAT(0, 3);
	mat_term = e.N2 * e.N1 * seepage_force;
	FILL_SEEP_MAT(1, 0);
	mat_term = e.N2 * e.N2 * seepage_force;
	FILL_SEEP_MAT(1, 1);
	mat_term = e.N2 * e.N3 * seepage_force;
	FILL_SEEP_MAT(1, 2);
	mat_term = e.N2 * e.N4 * seepage_force;
	FILL_SEEP_MAT(1, 3);
	mat_term = e.N3 * e.N1 * seepage_force;
	FILL_SEEP_MAT(2, 0);
	mat_term = e.N3 * e.N2 * seepage_force;
	FILL_SEEP_MAT(2, 1);
	mat_term = e.N3 * e.N3 * seepage_force;
	FILL_SEEP_MAT(2, 2);
	mat_term = e.N3 * e.N4 * seepage_force;
	FILL_SEEP_MAT(2, 3);
	mat_term = e.N4 * e.N1 * seepage_force;
	FILL_SEEP_MAT(3, 0);
	mat_term = e.N4 * e.N2 * seepage_force;
	FILL_SEEP_MAT(3, 1);
	mat_term = e.N4 * e.N3 * seepage_force;
	FILL_SEEP_MAT(3, 2);
	mat_term = e.N4 * e.N4 * seepage_force;
	FILL_SEEP_MAT(3, 3);

	// stiffness matrix
	double tmp;
	// K
	double E_mat_lambda = e.miu * e.E / ((1.0 + e.miu) * (1.0 - 2.0 * e.miu));
	double E_mat_miu = e.E / (2.0 * (1.0 + e.miu));
	Eigen::Matrix<double, 3, 3> E_mat;
	E_mat(0, 0) = E_mat_lambda + 2.0 * E_mat_miu;
	E_mat(0, 1) = E_mat_lambda;
	E_mat(0, 2) = 0.0;
	E_mat(1, 0) = E_mat_lambda;
	E_mat(1, 1) = E_mat_lambda + 2.0 * E_mat_miu;
	E_mat(1, 2) = 0.0;
	E_mat(2, 0) = 0.0;
	E_mat(2, 1) = 0.0;
	E_mat(2, 2) = E_mat_miu;
	Eigen::Matrix<double, 3, 8> dN_dx;
	dN_dx(0, 0) = e.dN1_dx;
	dN_dx(0, 1) = 0.0;
	dN_dx(0, 2) = e.dN2_dx;
	dN_dx(0, 3) = 0.0;
	dN_dx(0, 4) = e.dN3_dx;
	dN_dx(0, 5) = 0.0;
	dN_dx(0, 6) = e.dN4_dx;
	dN_dx(0, 7) = 0.0;
	dN_dx(1, 0) = 0.0;
	dN_dx(1, 1) = e.dN1_dy;
	dN_dx(1, 2) = 0.0;
	dN_dx(1, 3) = e.dN2_dy;
	dN_dx(1, 4) = 0.0;
	dN_dx(1, 5) = e.dN3_dy;
	dN_dx(1, 6) = 0.0;
	dN_dx(1, 7) = e.dN4_dy;
	dN_dx(2, 0) = e.dN1_dy;
	dN_dx(2, 1) = e.dN1_dx;
	dN_dx(2, 2) = e.dN2_dy;
	dN_dx(2, 3) = e.dN2_dx;
	dN_dx(2, 4) = e.dN3_dy;
	dN_dx(2, 5) = e.dN3_dx;
	dN_dx(2, 6) = e.dN4_dy;
	dN_dx(2, 7) = e.dN4_dx;
	Eigen::Matrix<double, 8, 8> K_mat = dN_dx.transpose() * E_mat * dN_dx;
	for (size_t i = 0; i < 8; ++i)
		for (size_t j = 0; j < 8; ++j)
			kmat[i][j] += K_mat(i, j) * elem_vol;
	// Gsx
	tmp = -(1.0 - e.n) * elem_vol;
	mat_term = e.dN1_dx * e.N1 * tmp;
	kmat[0][16] += mat_term;
	mat_term = e.dN1_dx * e.N2 * tmp;
	kmat[0][17] += mat_term;
	mat_term = e.dN1_dx * e.N3 * tmp;
	kmat[0][18] += mat_term;
	mat_term = e.dN1_dx * e.N4 * tmp;
	kmat[0][19] += mat_term;
	mat_term = e.dN2_dx * e.N1 * tmp;
	kmat[1][16] += mat_term;
	mat_term = e.dN2_dx * e.N2 * tmp;
	kmat[1][17] += mat_term;
	mat_term = e.dN2_dx * e.N3 * tmp;
	kmat[1][18] += mat_term;
	mat_term = e.dN2_dx * e.N4 * tmp;
	kmat[1][19] += mat_term;
	mat_term = e.dN3_dx * e.N1 * tmp;
	kmat[2][16] += mat_term;
	mat_term = e.dN3_dx * e.N2 * tmp;
	kmat[2][17] += mat_term;
	mat_term = e.dN3_dx * e.N3 * tmp;
	kmat[2][18] += mat_term;
	mat_term = e.dN3_dx * e.N4 * tmp;
	kmat[2][19] += mat_term;
	mat_term = e.dN4_dx * e.N1 * tmp;
	kmat[3][16] += mat_term;
	mat_term = e.dN4_dx * e.N2 * tmp;
	kmat[3][17] += mat_term;
	mat_term = e.dN4_dx * e.N3 * tmp;
	kmat[3][18] += mat_term;
	mat_term = e.dN4_dx * e.N4 * tmp;
	kmat[3][19] += mat_term;
	// Gsy
	mat_term = e.dN1_dy * e.N1 * tmp;
	kmat[4][16] += mat_term;
	mat_term = e.dN1_dy * e.N2 * tmp;
	kmat[4][17] += mat_term;
	mat_term = e.dN1_dy * e.N3 * tmp;
	kmat[4][18] += mat_term;
	mat_term = e.dN1_dy * e.N4 * tmp;
	kmat[4][19] += mat_term;
	mat_term = e.dN2_dy * e.N1 * tmp;
	kmat[5][16] += mat_term;
	mat_term = e.dN2_dy * e.N2 * tmp;
	kmat[5][17] += mat_term;
	mat_term = e.dN2_dy * e.N3 * tmp;
	kmat[5][18] += mat_term;
	mat_term = e.dN2_dy * e.N4 * tmp;
	kmat[5][19] += mat_term;
	mat_term = e.dN3_dy * e.N1 * tmp;
	kmat[6][16] += mat_term;
	mat_term = e.dN3_dy * e.N2 * tmp;
	kmat[6][17] += mat_term;
	mat_term = e.dN3_dy * e.N3 * tmp;
	kmat[6][18] += mat_term;
	mat_term = e.dN3_dy * e.N4 * tmp;
	kmat[6][19] += mat_term;
	mat_term = e.dN4_dy * e.N1 * tmp;
	kmat[7][16] += mat_term;
	mat_term = e.dN4_dy * e.N2 * tmp;
	kmat[7][17] += mat_term;
	mat_term = e.dN4_dy * e.N3 * tmp;
	kmat[7][18] += mat_term;
	mat_term = e.dN4_dy * e.N4 * tmp;
	kmat[7][19] += mat_term;
	// Gfx
	tmp = -e.n * elem_vol;
	mat_term = e.dN1_dx * e.N1 * tmp;
	kmat[8][16] += mat_term;
	mat_term = e.dN1_dx * e.N2 * tmp;
	kmat[8][17] += mat_term;
	mat_term = e.dN1_dx * e.N3 * tmp;
	kmat[8][18] += mat_term;
	mat_term = e.dN1_dx * e.N4 * tmp;
	kmat[8][19] += mat_term;
	mat_term = e.dN2_dx * e.N1 * tmp;
	kmat[9][16] += mat_term;
	mat_term = e.dN2_dx * e.N2 * tmp;
	kmat[9][17] += mat_term;
	mat_term = e.dN2_dx * e.N3 * tmp;
	kmat[9][18] += mat_term;
	mat_term = e.dN2_dx * e.N4 * tmp;
	kmat[9][19] += mat_term;
	mat_term = e.dN3_dx * e.N1 * tmp;
	kmat[10][16] += mat_term;
	mat_term = e.dN3_dx * e.N2 * tmp;
	kmat[10][17] += mat_term;
	mat_term = e.dN3_dx * e.N3 * tmp;
	kmat[10][18] += mat_term;
	mat_term = e.dN3_dx * e.N4 * tmp;
	kmat[10][19] += mat_term;
	mat_term = e.dN4_dx * e.N1 * tmp;
	kmat[11][16] += mat_term;
	mat_term = e.dN4_dx * e.N2 * tmp;
	kmat[11][17] += mat_term;
	mat_term = e.dN4_dx * e.N3 * tmp;
	kmat[11][18] += mat_term;
	mat_term = e.dN4_dx * e.N4 * tmp;
	kmat[11][19] += mat_term;
	// Gfy
	mat_term = e.dN1_dy * e.N1 * tmp;
	kmat[12][16] += mat_term;
	mat_term = e.dN1_dy * e.N2 * tmp;
	kmat[12][17] += mat_term;
	mat_term = e.dN1_dy * e.N3 * tmp;
	kmat[12][18] += mat_term;
	mat_term = e.dN1_dy * e.N4 * tmp;
	kmat[12][19] += mat_term;
	mat_term = e.dN2_dy * e.N1 * tmp;
	kmat[13][16] += mat_term;
	mat_term = e.dN2_dy * e.N2 * tmp;
	kmat[13][17] += mat_term;
	mat_term = e.dN2_dy * e.N3 * tmp;
	kmat[13][18] += mat_term;
	mat_term = e.dN2_dy * e.N4 * tmp;
	kmat[13][19] += mat_term;
	mat_term = e.dN3_dy * e.N1 * tmp;
	kmat[14][16] += mat_term;
	mat_term = e.dN3_dy * e.N2 * tmp;
	kmat[14][17] += mat_term;
	mat_term = e.dN3_dy * e.N3 * tmp;
	kmat[14][18] += mat_term;
	mat_term = e.dN3_dy * e.N4 * tmp;
	kmat[14][19] += mat_term;
	mat_term = e.dN4_dy * e.N1 * tmp;
	kmat[15][16] += mat_term;
	mat_term = e.dN4_dy * e.N2 * tmp;
	kmat[15][17] += mat_term;
	mat_term = e.dN4_dy * e.N3 * tmp;
	kmat[15][18] += mat_term;
	mat_term = e.dN4_dy * e.N4 * tmp;
	kmat[15][19] += mat_term;
	// Hsx
	tmp = e.Kf * (1.0 - e.n) / e.n * elem_vol;
	mat_term = e.N1 * e.dN1_dx * tmp;
	kmat[16][0] += mat_term;
	mat_term = e.N1 * e.dN2_dx * tmp;
	kmat[16][1] += mat_term;
	mat_term = e.N1 * e.dN3_dx * tmp;
	kmat[16][2] += mat_term;
	mat_term = e.N1 * e.dN4_dx * tmp;
	kmat[16][3] += mat_term;
	mat_term = e.N2 * e.dN1_dx * tmp;
	kmat[17][0] += mat_term;
	mat_term = e.N2 * e.dN2_dx * tmp;
	kmat[17][1] += mat_term;
	mat_term = e.N2 * e.dN3_dx * tmp;
	kmat[17][2] += mat_term;
	mat_term = e.N2 * e.dN4_dx * tmp;
	kmat[17][3] += mat_term;
	mat_term = e.N3 * e.dN1_dx * tmp;
	kmat[18][0] += mat_term;
	mat_term = e.N3 * e.dN2_dx * tmp;
	kmat[18][1] += mat_term;
	mat_term = e.N3 * e.dN3_dx * tmp;
	kmat[18][2] += mat_term;
	mat_term = e.N3 * e.dN4_dx * tmp;
	kmat[18][3] += mat_term;
	mat_term = e.N4 * e.dN1_dx * tmp;
	kmat[19][0] += mat_term;
	mat_term = e.N4 * e.dN2_dx * tmp;
	kmat[19][1] += mat_term;
	mat_term = e.N4 * e.dN3_dx * tmp;
	kmat[19][2] += mat_term;
	mat_term = e.N4 * e.dN4_dx * tmp;
	kmat[19][3] += mat_term;
	// Hsy
	mat_term = e.N1 * e.dN1_dy * tmp;
	kmat[16][4] += mat_term;
	mat_term = e.N1 * e.dN2_dy * tmp;
	kmat[16][5] += mat_term;
	mat_term = e.N1 * e.dN3_dy * tmp;
	kmat[16][6] += mat_term;
	mat_term = e.N1 * e.dN4_dy * tmp;
	kmat[16][7] += mat_term;
	mat_term = e.N2 * e.dN1_dy * tmp;
	kmat[17][4] += mat_term;
	mat_term = e.N2 * e.dN2_dy * tmp;
	kmat[17][5] += mat_term;
	mat_term = e.N2 * e.dN3_dy * tmp;
	kmat[17][6] += mat_term;
	mat_term = e.N2 * e.dN4_dy * tmp;
	kmat[17][7] += mat_term;
	mat_term = e.N3 * e.dN1_dy * tmp;
	kmat[18][4] += mat_term;
	mat_term = e.N3 * e.dN2_dy * tmp;
	kmat[18][5] += mat_term;
	mat_term = e.N3 * e.dN3_dy * tmp;
	kmat[18][6] += mat_term;
	mat_term = e.N3 * e.dN4_dy * tmp;
	kmat[18][7] += mat_term;
	mat_term = e.N4 * e.dN1_dy * tmp;
	kmat[19][4] += mat_term;
	mat_term = e.N4 * e.dN2_dy * tmp;
	kmat[19][5] += mat_term;
	mat_term = e.N4 * e.dN3_dy * tmp;
	kmat[19][6] += mat_term;
	mat_term = e.N4 * e.dN4_dy * tmp;
	kmat[19][7] += mat_term;
	// Hfx
	tmp = e.Kf * elem_vol;
	mat_term = e.N1 * e.dN1_dx * tmp;
	kmat[16][8] += mat_term;
	mat_term = e.N1 * e.dN2_dx * tmp;
	kmat[16][9] += mat_term;
	mat_term = e.N1 * e.dN3_dx * tmp;
	kmat[16][10] += mat_term;
	mat_term = e.N1 * e.dN4_dx * tmp;
	kmat[16][11] += mat_term;
	mat_term = e.N2 * e.dN1_dx * tmp;
	kmat[17][8] += mat_term;
	mat_term = e.N2 * e.dN2_dx * tmp;
	kmat[17][9] += mat_term;
	mat_term = e.N2 * e.dN3_dx * tmp;
	kmat[17][10] += mat_term;
	mat_term = e.N2 * e.dN4_dx * tmp;
	kmat[17][11] += mat_term;
	mat_term = e.N3 * e.dN1_dx * tmp;
	kmat[18][8] += mat_term;
	mat_term = e.N3 * e.dN2_dx * tmp;
	kmat[18][9] += mat_term;
	mat_term = e.N3 * e.dN3_dx * tmp;
	kmat[18][10] += mat_term;
	mat_term = e.N3 * e.dN4_dx * tmp;
	kmat[18][11] += mat_term;
	mat_term = e.N4 * e.dN1_dx * tmp;
	kmat[19][8] += mat_term;
	mat_term = e.N4 * e.dN2_dx * tmp;
	kmat[19][9] += mat_term;
	mat_term = e.N4 * e.dN3_dx * tmp;
	kmat[19][10] += mat_term;
	mat_term = e.N4 * e.dN4_dx * tmp;
	kmat[19][11] += mat_term;
	// Hfy
	mat_term = e.N1 * e.dN1_dy * tmp;
	kmat[16][12] += mat_term;
	mat_term = e.N1 * e.dN2_dy * tmp;
	kmat[16][13] += mat_term;
	mat_term = e.N1 * e.dN3_dy * tmp;
	kmat[16][14] += mat_term;
	mat_term = e.N1 * e.dN4_dy * tmp;
	kmat[16][15] += mat_term;
	mat_term = e.N2 * e.dN1_dy * tmp;
	kmat[17][12] += mat_term;
	mat_term = e.N2 * e.dN2_dy * tmp;
	kmat[17][13] += mat_term;
	mat_term = e.N2 * e.dN3_dy * tmp;
	kmat[17][14] += mat_term;
	mat_term = e.N2 * e.dN4_dy * tmp;
	kmat[17][15] += mat_term;
	mat_term = e.N3 * e.dN1_dy * tmp;
	kmat[18][12] += mat_term;
	mat_term = e.N3 * e.dN2_dy * tmp;
	kmat[18][13] += mat_term;
	mat_term = e.N3 * e.dN3_dy * tmp;
	kmat[18][14] += mat_term;
	mat_term = e.N3 * e.dN4_dy * tmp;
	kmat[18][15] += mat_term;
	mat_term = e.N4 * e.dN1_dy * tmp;
	kmat[19][12] += mat_term;
	mat_term = e.N4 * e.dN2_dy * tmp;
	kmat[19][13] += mat_term;
	mat_term = e.N4 * e.dN3_dy * tmp;
	kmat[19][14] += mat_term;
	mat_term = e.N4 * e.dN4_dy * tmp;
	kmat[19][15] += mat_term;
	// P
	mat_term = e.N1 * e.N1 * elem_vol;
	kmat[16][16] += mat_term;
	mat_term = e.N1 * e.N2 * elem_vol;
	kmat[16][17] += mat_term;
	mat_term = e.N1 * e.N3 * elem_vol;
	kmat[16][18] += mat_term;
	mat_term = e.N1 * e.N4 * elem_vol;
	kmat[16][19] += mat_term;
	mat_term = e.N2 * e.N1 * elem_vol;
	kmat[17][16] += mat_term;
	mat_term = e.N2 * e.N2 * elem_vol;
	kmat[17][17] += mat_term;
	mat_term = e.N2 * e.N3 * elem_vol;
	kmat[17][18] += mat_term;
	mat_term = e.N2 * e.N4 * elem_vol;
	kmat[17][19] += mat_term;
	mat_term = e.N3 * e.N1 * elem_vol;
	kmat[18][16] += mat_term;
	mat_term = e.N3 * e.N2 * elem_vol;
	kmat[18][17] += mat_term;
	mat_term = e.N3 * e.N3 * elem_vol;
	kmat[18][18] += mat_term;
	mat_term = e.N3 * e.N4 * elem_vol;
	kmat[18][19] += mat_term;
	mat_term = e.N4 * e.N1 * elem_vol;
	kmat[19][16] += mat_term;
	mat_term = e.N4 * e.N2 * elem_vol;
	kmat[19][17] += mat_term;
	mat_term = e.N4 * e.N3 * elem_vol;
	kmat[19][18] += mat_term;
	mat_term = e.N4 * e.N4 * elem_vol;
	kmat[19][19] += mat_term;
}

void Step_S2D_CHM_s_FEM_uUp::form_elem_force_vec(
	Model_S2D_CHM_s_FEM_uUp::Element &e, double dt, double fvec[20])
{
	double elem_vol = model->elem_vol;
	Node_fem &n1 = model->nodes[e.index_y * model->node_x_num + e.index_x];
	Node_fem &n2 = *(&n1 + 1);
	Node_fem &n3 = *(&n2 + model->node_x_num);
	Node_fem &n4 = *(&n3 - 1);

	// f_int
	double elem_p = n1.p * e.N1 + n2.p * e.N2 + n3.p * e.N3 + n4.p * e.N4;
	fvec[0] = -(e.dN1_dx * (e.s11 - (1.0 - e.n) * elem_p) + e.dN1_dy * e.s12) * elem_vol;
	fvec[1] = -(e.dN2_dx * (e.s11 - (1.0 - e.n) * elem_p) + e.dN2_dy * e.s12) * elem_vol;
	fvec[2] = -(e.dN3_dx * (e.s11 - (1.0 - e.n) * elem_p) + e.dN3_dy * e.s12) * elem_vol;
	fvec[3] = -(e.dN4_dx * (e.s11 - (1.0 - e.n) * elem_p) + e.dN4_dy * e.s12) * elem_vol;
	fvec[4] = -(e.dN1_dx * e.s12 + e.dN1_dy * (e.s22 - (1.0 - e.n) * elem_p)) * elem_vol;
	fvec[5] = -(e.dN2_dx * e.s12 + e.dN2_dy * (e.s22 - (1.0 - e.n) * elem_p)) * elem_vol;
	fvec[6] = -(e.dN3_dx * e.s12 + e.dN3_dy * (e.s22 - (1.0 - e.n) * elem_p)) * elem_vol;
	fvec[7] = -(e.dN4_dx * e.s12 + e.dN4_dy * (e.s22 - (1.0 - e.n) * elem_p)) * elem_vol;
	fvec[8] = -e.dN1_dx * e.n * -elem_p * elem_vol;
	fvec[9] = -e.dN2_dx * e.n * -elem_p * elem_vol;
	fvec[10] = -e.dN3_dx * e.n * -elem_p * elem_vol;
	fvec[11] = -e.dN4_dx * e.n * -elem_p * elem_vol;
	fvec[12] = -e.dN1_dy * e.n * -elem_p * elem_vol;
	fvec[13] = -e.dN2_dy * e.n * -elem_p * elem_vol;
	fvec[14] = -e.dN3_dy * e.n * -elem_p * elem_vol;
	fvec[15] = -e.dN4_dy * e.n * -elem_p * elem_vol;
	fvec[16] = 0.0;
	fvec[17] = 0.0;
	fvec[18] = 0.0;
	fvec[19] = 0.0;

	double tmp;
	// a
	double elem_a;
	// solid
	tmp = ((1.0 / (2.0 * beta) - 1.0) * (1.0 - e.n) * e.density_s - dt * (1.0 - gamma / (2.0 * beta)) * e.n * e.n * e.miu / e.k) * elem_vol;
	elem_a = n1.ax_s * e.N1 + n2.ax_s * e.N2 + n3.ax_s * e.N3 + n4.ax_s * e.N4;
	fvec[0] += e.N1 * tmp * elem_a;
	fvec[1] += e.N2 * tmp * elem_a;
	fvec[2] += e.N3 * tmp * elem_a;
	fvec[3] += e.N4 * tmp * elem_a;
	elem_a = n1.ay_s * e.N1 + n2.ay_s * e.N2 + n3.ay_s * e.N3 + n4.ay_s * e.N4;
	fvec[4] += e.N1 * tmp * elem_a;
	fvec[5] += e.N2 * tmp * elem_a;
	fvec[6] += e.N3 * tmp * elem_a;
	fvec[7] += e.N4 * tmp * elem_a;
	// fluid
	tmp = ((1.0 / (2.0 * beta) - 1.0) * e.n * e.density_f - dt * (1.0 - gamma / (2.0 * beta)) * e.n * e.n * e.miu / e.k) * elem_vol;
	elem_a = n1.ax_f * e.N1 + n2.ax_f * e.N2 + n3.ax_f * e.N3 + n4.ax_f * e.N4;
	fvec[8] += e.N1 * tmp * elem_a;
	fvec[9] += e.N2 * tmp * elem_a;
	fvec[10] += e.N3 * tmp * elem_a;
	fvec[11] += e.N4 * tmp * elem_a;
	elem_a = n1.ay_f * e.N1 + n2.ay_f * e.N2 + n3.ay_f * e.N3 + n4.ay_f * e.N4;
	fvec[12] += e.N1 * tmp * elem_a;
	fvec[13] += e.N2 * tmp * elem_a;
	fvec[14] += e.N3 * tmp * elem_a;
	fvec[15] += e.N4 * tmp * elem_a;

	// v
	double elem_v;
	// solid
	tmp = (1.0 / (beta * dt) * (1.0 - e.n) * e.density_s + (gamma / beta - 1.0) * e.n * e.n * e.miu / e.k) * elem_vol;
	elem_v = n1.vx_s * e.N1 + n2.vx_s * e.N2 + n3.vx_s * e.N3 + n4.vx_s * e.N4;
	fvec[0] += e.N1 * tmp * elem_v;
	fvec[1] += e.N2 * tmp * elem_v;
	fvec[2] += e.N3 * tmp * elem_v;
	fvec[3] += e.N4 * tmp * elem_v;
	elem_v = n1.vy_s * e.N1 + n2.vy_s * e.N2 + n3.vy_s * e.N3 + n4.vy_s * e.N4;
	fvec[4] += e.N1 * tmp * elem_v;
	fvec[5] += e.N2 * tmp * elem_v;
	fvec[6] += e.N3 * tmp * elem_v;
	fvec[7] += e.N4 * tmp * elem_v;
	// fluid
	tmp = (1.0 / (beta * dt) * e.n * e.density_f + (gamma / beta - 1.0) * e.n * e.n * e.miu / e.k) * elem_vol;
	elem_v = n1.vx_f * e.N1 + n2.vx_f * e.N2 + n3.vx_f * e.N3 + n4.vx_f * e.N4;
	fvec[8] += e.N1 * tmp * elem_v;
	fvec[9] += e.N2 * tmp * elem_v;
	fvec[10] += e.N3 * tmp * elem_v;
	fvec[11] += e.N4 * tmp * elem_v;
	elem_v = n1.vy_f * e.N1 + n2.vy_f * e.N2 + n3.vy_f * e.N3 + n4.vy_f * e.N4;
	fvec[12] += e.N1 * tmp * elem_v;
	fvec[13] += e.N2 * tmp * elem_v;
	fvec[14] += e.N3 * tmp * elem_v;
	fvec[15] += e.N4 * tmp * elem_v;
}