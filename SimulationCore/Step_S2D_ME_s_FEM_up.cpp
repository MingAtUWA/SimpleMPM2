#include "SimulationCore_pcp.h"

#include <cstdio>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include "Step_S2D_ME_s_FEM_up.h"

// Newmark-beta method parameters
#define beta 0.25
#define gamma 0.5

namespace
{
	typedef Model_S2D_ME_s_FEM_up::Node Node_fem;
	typedef Model_S2D_ME_s_FEM_up::Element Element_fem;
	typedef Model_S2D_ME_s_FEM_up::GaussPoint GaussPoint_fem;
	typedef Model_S2D_ME_s_FEM_up::ShapeFuncValue ShapeFuncValue;
	typedef Model_S2D_ME_s_FEM_up::DOF DOF;

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

	void update_gauss_point(GaussPoint_fem &gp,
			ShapeFuncValue &sf, double dt,
			double dux1, double dux2, double dux3, double dux4,
			double duy1, double duy2, double duy3, double duy4,
			double dp1, double dp2, double dp3, double dp4)
	{
		double dux, duy, dvx, dvy, dax, day;
		double de11, de22, de12;
		double ds11, ds22, ds12;
		double de_vol;
		double dt2 = dt * dt;
		// displacement increment
		dux = dux1 * sf.N1 + dux2 * sf.N2 + dux3 * sf.N3 + dux4 * sf.N4;
		duy = duy1 * sf.N1 + duy2 * sf.N2 + duy3 * sf.N3 + duy4 * sf.N4;
		// velocity increment
		dvx = gamma / (beta * dt) * dux - gamma / beta * gp.vx
			+ dt * (1.0 - gamma / (2.0 * beta)) * gp.ax;
		dvy = gamma / (beta * dt) * duy - gamma / beta * gp.vy
			+ dt * (1.0 - gamma / (2.0 * beta)) * gp.ay;
		dax = 1.0 / (beta * dt2) * dux - 1.0 / (beta * dt) * gp.vx
			- 1.0 / (2.0 * beta) * gp.ax;
		day = 1.0 / (beta * dt2) * duy - 1.0 / (beta * dt) * gp.vy
			- 1.0 / (2.0 * beta) * gp.ay;
		// displacement
		gp.ux += dux;
		gp.uy += duy;
		// velocity
		gp.vx += dvx;
		gp.vy += dvy;
		// acceleration
		gp.ax += dax;
		gp.ay += day;
		// strain increment
		de11 = dux1 * sf.dN1_dx + dux2 * sf.dN2_dx + dux3 * sf.dN3_dx + dux4 * sf.dN4_dx;
		de22 = duy1 * sf.dN1_dy + duy2 * sf.dN2_dy + duy3 * sf.dN3_dy + duy4 * sf.dN4_dy;
		de12 = (dux1 * sf.dN1_dy + dux2 * sf.dN2_dy + dux3 * sf.dN3_dy + dux4 * sf.dN4_dy
			  + duy1 * sf.dN1_dx + duy2 * sf.dN2_dx + duy3 * sf.dN3_dx + duy4 * sf.dN4_dx) * 0.5;
		gp.e11 += de11;
		gp.e22 += de22;
		gp.e12 += de12;
		// update stress
		double G = gp.E / (2.0 * (1.0 + gp.niu));
		ds11 =  4.0 / 3.0 * G * de11 - 2.0 / 3.0 * G * de22;
		ds22 = -2.0 / 3.0 * G * de11 + 4.0 / 3.0 * G * de22;
		ds12 = G * 2.0 * de12;
		gp.s11 += ds11;
		gp.s22 += ds22;
		gp.s12 += ds12;
		// pore pressure
		gp.p += dp1 * sf.N1 + dp2 * sf.N2 + dp3 * sf.N3 + dp4 * sf.N4;
		// density
		de_vol = de11 + de22;
		gp.density /= (1.0 + de_vol);
	}
};

Step_S2D_ME_s_FEM_up::Step_S2D_ME_s_FEM_up() :
	Step(&solve_substep_S2D_ME_s_FEM_up),
	model(nullptr), kmat_col(nullptr) {}

Step_S2D_ME_s_FEM_up::~Step_S2D_ME_s_FEM_up() {}

int Step_S2D_ME_s_FEM_up::init_calculation(void)
{
	if (is_first_step) {}

	kmat_col = new double[model->dof_num];

	is_first_substep = true;

	return 0;
}

int Step_S2D_ME_s_FEM_up::finalize_calculation(void)
{
	if (kmat_col)
	{
		delete[] kmat_col;
		kmat_col = nullptr;
	}

	return 0;
}

int solve_substep_S2D_ME_s_FEM_up(void *_self)
{
	Step_S2D_ME_s_FEM_up &self = *(Step_S2D_ME_s_FEM_up *)(_self);
	Model_S2D_ME_s_FEM_up &model = *self.model;
	double dt = self.dtime;
	double dt2 = dt * dt;
	
	// list of non-zeros coefficients
	MatrixCoefficientSet<> &g_kmat_coefs = self.g_kmat_coefs;
	g_kmat_coefs.init(model.dof_num);
	Eigen::SparseMatrix<double> g_kmat(model.dof_num, model.dof_num);
	Eigen::VectorXd g_fvec = Eigen::VectorXd::Zero(model.dof_num);

	size_t l2g_dof_id_map[12];
	double e_kmat[12][12], e_fvec[12];
	for (size_t e_id = 0; e_id < model.elem_num; ++e_id)
	{
		Element_fem &e = model.elems[e_id];
		// form elemental stiffness matrix and force vector
		self.form_elem_stiffness_mat_and_force_vec(e, e_kmat, e_fvec);
		//print_mat(e_kmat);
		//print_vec(e_fvec);

		// map from local dof id to global dof id
		// ux
		l2g_dof_id_map[0] = model.n_id_to_dof_id(e.n1_id, DOF::ux);
		l2g_dof_id_map[1] = model.n_id_to_dof_id(e.n2_id, DOF::ux);
		l2g_dof_id_map[2] = model.n_id_to_dof_id(e.n3_id, DOF::ux);
		l2g_dof_id_map[3] = model.n_id_to_dof_id(e.n4_id, DOF::ux);
		// uy
		l2g_dof_id_map[4] = model.n_id_to_dof_id(e.n1_id, DOF::uy);
		l2g_dof_id_map[5] = model.n_id_to_dof_id(e.n2_id, DOF::uy);
		l2g_dof_id_map[6] = model.n_id_to_dof_id(e.n3_id, DOF::uy);
		l2g_dof_id_map[7] = model.n_id_to_dof_id(e.n4_id, DOF::uy);
		// p
		l2g_dof_id_map[8] = model.n_id_to_dof_id(e.n1_id, DOF::p);
		l2g_dof_id_map[9] = model.n_id_to_dof_id(e.n2_id, DOF::p);
		l2g_dof_id_map[10] = model.n_id_to_dof_id(e.n3_id, DOF::p);
		l2g_dof_id_map[11] = model.n_id_to_dof_id(e.n4_id, DOF::p);
		//print_vec(l2g_dof_id_map);
		// add to global matrix and vector
		size_t g_dof_id1, g_dof_id2;
		for (size_t l_id1 = 0; l_id1 < 12; ++l_id1)
		{
			g_dof_id1 = l2g_dof_id_map[l_id1];
			for (size_t l_id2 = 0; l_id2 < 12; ++l_id2)
			{
				g_dof_id2 = l2g_dof_id_map[l_id2];
				g_kmat_coefs.add_coefficient(g_dof_id1, g_dof_id2, e_kmat[l_id1][l_id2]);
			}
		}
		size_t g_dof_id;
		for (size_t l_id = 0; l_id < 12; ++l_id)
		{
			g_dof_id = l2g_dof_id_map[l_id];
			g_fvec[g_dof_id] += e_fvec[l_id];
		}
	}

	// apply external force
	// traction force
	size_t dof_g_id;
	double nfs[4];
	if (self.is_first_substep)
	{
		self.is_first_substep = false;
		for (size_t t_id = 0; t_id < model.tx_num; ++t_id)
		{
			TractionBC_2DFEM &tx = model.txs[t_id];
			model.cal_traction_bc(tx, nfs);
			Element_fem &e = model.elems[tx.elem_id];
			dof_g_id = model.n_id_to_dof_id(e.n1_id, DOF::ux);
			g_fvec[dof_g_id] += nfs[0];
			dof_g_id = model.n_id_to_dof_id(e.n2_id, DOF::ux);
			g_fvec[dof_g_id] += nfs[1];
			dof_g_id = model.n_id_to_dof_id(e.n3_id, DOF::ux);
			g_fvec[dof_g_id] += nfs[2];
			dof_g_id = model.n_id_to_dof_id(e.n4_id, DOF::ux);
			g_fvec[dof_g_id] += nfs[3];
		}
		for (size_t t_id = 0; t_id < model.ty_num; ++t_id)
		{
			TractionBC_2DFEM &ty = model.tys[t_id];
			model.cal_traction_bc(ty, nfs);
			Element_fem &e = model.elems[ty.elem_id];
			dof_g_id = model.n_id_to_dof_id(e.n1_id, DOF::uy);
			g_fvec[dof_g_id] += nfs[0];
			dof_g_id = model.n_id_to_dof_id(e.n2_id, DOF::uy);
			g_fvec[dof_g_id] += nfs[1];
			dof_g_id = model.n_id_to_dof_id(e.n3_id, DOF::uy);
			g_fvec[dof_g_id] += nfs[2];
			dof_g_id = model.n_id_to_dof_id(e.n4_id, DOF::uy);
			g_fvec[dof_g_id] += nfs[3];
		}
		// body force to be finished...
	}
	//std::cout << "f:\n" << g_fvec << "\n";

	// apply displacement boundary condition
	double dig_term;
	for (size_t bc_id = 0; bc_id < model.ux_num; ++bc_id)
	{
		DisplacementBC &dbc = model.uxs[bc_id];
		dof_g_id = model.n_id_to_dof_id(dbc.node_id, DOF::ux);
		dig_term = g_kmat_coefs.del_col_and_row(dof_g_id, self.kmat_col);
		for (size_t row_id = 0; row_id < model.dof_num; ++row_id)
			g_fvec[row_id] -= self.kmat_col[row_id] * dbc.u;
		g_fvec[dof_g_id] = dig_term * dbc.u;
	}
	for (size_t bc_id = 0; bc_id < model.uy_num; ++bc_id)
	{
		DisplacementBC &dbc = model.uys[bc_id];
		dof_g_id = model.n_id_to_dof_id(dbc.node_id, DOF::uy);
		dig_term = g_kmat_coefs.del_col_and_row(dof_g_id, self.kmat_col);
		for (size_t row_id = 0; row_id < model.dof_num; ++row_id)
			g_fvec[row_id] -= self.kmat_col[row_id] * dbc.u;
		g_fvec[dof_g_id] = dig_term * dbc.u;
	}
	//for (size_t bc_id = 0; bc_id < model.pbc_num; ++bc_id)
	//{
	//	PressureBC &pbc = model.pbcs[bc_id];
	//	dof_g_id = model.n_id_to_dof_id(pbc.node_id, DOF::p);
	//	dig_term = g_kmat_coefs.del_col_and_row(dof_g_id, self.kmat_col);
	//	for (size_t row_id = 0; row_id < model.dof_num; ++row_id)
	//		g_fvec[row_id] -= self.kmat_col[row_id] * pbc.p;
	//	g_fvec[dof_g_id] = dig_term * pbc.p;
	//}

	g_kmat.setFromTriplets(g_kmat_coefs.begin(), g_kmat_coefs.end());
	Eigen::SparseLU<Eigen::SparseMatrix<double> > solver(g_kmat);
	Eigen::VectorXd g_du_vec = solver.solve(g_fvec);
	//std::cout << g_kmat << "\n";
	//std::cout << g_fvec << "\n";
	//std::cout << g_du_vec << "\n";

	// reapply disp bc
	for (size_t bc_id = 0; bc_id < model.ux_num; ++bc_id)
	{
		DisplacementBC &dbc = model.uxs[bc_id];
		dof_g_id = model.n_id_to_dof_id(dbc.node_id, DOF::ux);
		g_du_vec[dof_g_id] = dbc.u;
	}
	for (size_t bc_id = 0; bc_id < model.uy_num; ++bc_id)
	{
		DisplacementBC &dbc = model.uys[bc_id];
		dof_g_id = model.n_id_to_dof_id(dbc.node_id, DOF::uy);
		g_du_vec[dof_g_id] = dbc.u;
	}
	//for (size_t bc_id = 0; bc_id < model.pbc_num; ++bc_id)
	//{
	//	PressureBC &pbc = model.pbcs[bc_id];
	//	dof_g_id = model.n_id_to_dof_id(pbc.node_id, DOF::p);
	//	g_du_vec[dof_g_id] = pbc.p;
	//}

	// Update node variables
	for (size_t n_id = 0; n_id < model.node_num; ++n_id)
	{
		Node_fem &n = model.nodes[n_id];
		n.ux += g_du_vec[model.n_id_to_dof_id(n_id, DOF::ux)];
		n.uy += g_du_vec[model.n_id_to_dof_id(n_id, DOF::uy)];
		n.p  += g_du_vec[model.n_id_to_dof_id(n_id, DOF::p)];
	}

	// Update guass point variables
	double dux1, duy1, dp1;
	double dux2, duy2, dp2;
	double dux3, duy3, dp3;
	double dux4, duy4, dp4;
	ShapeFuncValue &sf1 = model.gp1_sf;
	ShapeFuncValue &sf2 = model.gp2_sf;
	ShapeFuncValue &sf3 = model.gp3_sf;
	ShapeFuncValue &sf4 = model.gp4_sf;
	for (size_t e_id = 0; e_id < model.elem_num; ++e_id)
	{
		Element_fem &e = model.elems[e_id];
		// dux
		dux1 = g_du_vec[model.n_id_to_dof_id(e.n1_id, DOF::ux)];
		dux2 = g_du_vec[model.n_id_to_dof_id(e.n2_id, DOF::ux)];
		dux3 = g_du_vec[model.n_id_to_dof_id(e.n3_id, DOF::ux)];
		dux4 = g_du_vec[model.n_id_to_dof_id(e.n4_id, DOF::ux)];
		// duy
		duy1 = g_du_vec[model.n_id_to_dof_id(e.n1_id, DOF::uy)];
		duy2 = g_du_vec[model.n_id_to_dof_id(e.n2_id, DOF::uy)];
		duy3 = g_du_vec[model.n_id_to_dof_id(e.n3_id, DOF::uy)];
		duy4 = g_du_vec[model.n_id_to_dof_id(e.n4_id, DOF::uy)];
		// dp
		dp1 = g_du_vec[model.n_id_to_dof_id(e.n1_id, DOF::p)];
		dp2 = g_du_vec[model.n_id_to_dof_id(e.n2_id, DOF::p)];
		dp3 = g_du_vec[model.n_id_to_dof_id(e.n3_id, DOF::p)];
		dp4 = g_du_vec[model.n_id_to_dof_id(e.n4_id, DOF::p)];

		update_gauss_point(e.gp1, sf1, self.dtime,
			dux1, dux2, dux3, dux4, duy1, duy2, duy3, duy4,
			dp1, dp2, dp3, dp4);
		update_gauss_point(e.gp2, sf2, self.dtime,
			dux1, dux2, dux3, dux4, duy1, duy2, duy3, duy4,
			dp1, dp2, dp3, dp4);
		update_gauss_point(e.gp3, sf3, self.dtime,
			dux1, dux2, dux3, dux4, duy1, duy2, duy3, duy4,
			dp1, dp2, dp3, dp4);
		update_gauss_point(e.gp4, sf4, self.dtime,
			dux1, dux2, dux3, dux4, duy1, duy2, duy3, duy4,
			dp1, dp2, dp3, dp4);
	}

	return 0;
}

namespace
{
	void cal_stiffness_mat(double k_mat[12][12], double E[3][3], double dN_dx[3][8], double weight)
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
							  + dN_dx[2][i] * mid_mat[2][j]) * weight;
	}

	void form_E_matrix(double E_mat[3][3], GaussPoint_fem &gp)
	{
		double G = gp.E / (2.0 * (1.0 + gp.niu));
		E_mat[0][0] =  4.0 / 3.0 * G;
		E_mat[0][1] = -2.0 / 3.0 * G;
		E_mat[0][2] = 0.0;
		E_mat[1][0] = -2.0 / 3.0 * G;
		E_mat[1][1] =  4.0 / 3.0 * G;
		E_mat[1][2] = 0.0;
		E_mat[2][0] = 0.0;
		E_mat[2][1] = 0.0;
		E_mat[2][2] = G;
	}
};

void Step_S2D_ME_s_FEM_up::form_elem_stiffness_mat_and_force_vec(
	Model_S2D_ME_s_FEM_up::Element &e, double kmat[12][12], double fvec[12])
{
	ShapeFuncValue &sf1 = model->gp1_sf;
	ShapeFuncValue &sf2 = model->gp2_sf;
	ShapeFuncValue &sf3 = model->gp3_sf;
	ShapeFuncValue &sf4 = model->gp4_sf;
	GaussPoint_fem &gp1 = e.gp1;
	GaussPoint_fem &gp2 = e.gp2;
	GaussPoint_fem &gp3 = e.gp3;
	GaussPoint_fem &gp4 = e.gp4;
	double det_dx_dxi = model->h * model->h * 0.25;
	double mat_term;
	union
	{
		struct
		{
			double m_mat[12][12];
			double k_mat[12][12];
			double f_a[8];
			double f_v[8];
		};
		char tmp_mem[1];
	};

	memset(tmp_mem, 0, sizeof(m_mat) + sizeof(k_mat) + sizeof(f_a) + sizeof(f_v));
	// mass matrix
	mat_term = (sf1.N1 * gp1.density * sf1.N1 + sf2.N1 * gp2.density * sf2.N1
			  + sf3.N1 * gp3.density * sf3.N1 + sf4.N1 * gp4.density * sf4.N1) * det_dx_dxi;
	m_mat[0][0] = mat_term;
	m_mat[4][4] = mat_term;
	mat_term = (sf1.N1 * gp1.density * sf1.N2 + sf2.N1 * gp2.density * sf2.N2
			  + sf3.N1 * gp3.density * sf3.N2 + sf4.N1 * gp4.density * sf4.N2) * det_dx_dxi;
	m_mat[0][1] = mat_term;
	m_mat[1][0] = mat_term;
	m_mat[4][5] = mat_term;
	m_mat[5][4] = mat_term;
	mat_term = (sf1.N1 * gp1.density * sf1.N3 + sf2.N1 * gp2.density * sf2.N3
			  + sf3.N1 * gp3.density * sf3.N3 + sf4.N1 * gp4.density * sf4.N3) * det_dx_dxi;
	m_mat[0][2] = mat_term;
	m_mat[2][0] = mat_term;
	m_mat[4][6] = mat_term;
	m_mat[6][4] = mat_term;
	mat_term = (sf1.N1 * gp1.density * sf1.N4 + sf2.N1 * gp2.density * sf2.N4
			  + sf3.N1 * gp3.density * sf3.N4 + sf4.N1 * gp4.density * sf4.N4) * det_dx_dxi;
	m_mat[0][3] = mat_term;
	m_mat[3][0] = mat_term;
	m_mat[4][7] = mat_term;
	m_mat[7][4] = mat_term;
	mat_term = (sf1.N2 * gp1.density * sf1.N2 + sf2.N2 * gp2.density * sf2.N2
			  + sf3.N2 * gp3.density * sf3.N2 + sf4.N2 * gp4.density * sf4.N2) * det_dx_dxi;
	m_mat[1][1] = mat_term;
	m_mat[5][5] = mat_term;
	mat_term = (sf1.N2 * gp1.density * sf1.N3 + sf2.N2 * gp2.density * sf2.N3
			  + sf3.N2 * gp3.density * sf3.N3 + sf4.N2 * gp4.density * sf4.N3) * det_dx_dxi;
	m_mat[1][2] = mat_term;
	m_mat[2][1] = mat_term;
	m_mat[5][6] = mat_term;
	m_mat[6][5] = mat_term;
	mat_term = (sf1.N2 * gp1.density * sf1.N4 + sf2.N2 * gp2.density * sf2.N4
			  + sf3.N2 * gp3.density * sf3.N4 + sf4.N2 * gp4.density * sf4.N4) * det_dx_dxi;
	m_mat[1][3] = mat_term;
	m_mat[3][1] = mat_term;
	m_mat[5][7] = mat_term;
	m_mat[7][5] = mat_term;
	mat_term = (sf1.N3 * gp1.density * sf1.N3 + sf2.N3 * gp2.density * sf2.N3
			  + sf3.N3 * gp3.density * sf3.N3 + sf4.N3 * gp4.density * sf4.N3) * det_dx_dxi;
	m_mat[2][2] = mat_term;
	m_mat[6][6] = mat_term;
	mat_term = (sf1.N3 * gp1.density * sf1.N4 + sf2.N3 * gp2.density * sf2.N4
			  + sf3.N3 * gp3.density * sf3.N4 + sf4.N3 * gp4.density * sf4.N4) * det_dx_dxi;
	m_mat[2][3] = mat_term;
	m_mat[3][2] = mat_term;
	m_mat[6][7] = mat_term;
	m_mat[7][6] = mat_term;
	mat_term = (sf1.N4 * gp1.density * sf1.N4 + sf2.N4 * gp2.density * sf2.N4
			  + sf3.N4 * gp3.density * sf3.N4 + sf4.N4 * gp4.density * sf4.N4) * det_dx_dxi;
	m_mat[3][3] = mat_term;
	m_mat[7][7] = mat_term;

	//print_mat(m_mat);

	// stiffness matrix
	// B * D * B
	double E_mat[3][3];
	// gauss point 1
	form_E_matrix(E_mat, gp1);
	cal_stiffness_mat(k_mat, E_mat, model->dN_dx_mat1, det_dx_dxi);
	// gauss point 2
	form_E_matrix(E_mat, gp2);
	cal_stiffness_mat(k_mat, E_mat, model->dN_dx_mat2, det_dx_dxi);
	// gauss point 3
	form_E_matrix(E_mat, gp3);
	cal_stiffness_mat(k_mat, E_mat, model->dN_dx_mat3, det_dx_dxi);
	// gauss point 4
	form_E_matrix(E_mat, gp4);
	cal_stiffness_mat(k_mat, E_mat, model->dN_dx_mat4, det_dx_dxi);
	
	// dNi_dx * Nj
	mat_term = (sf1.dN1_dx * sf1.N1 + sf2.dN1_dx * sf2.N1
			  + sf3.dN1_dx * sf3.N1 + sf4.dN1_dx * sf4.N1) * det_dx_dxi;
	k_mat[0][8] = mat_term;
	k_mat[8][0] = mat_term;
	mat_term = (sf1.dN1_dx * sf1.N2 + sf2.dN1_dx * sf2.N2
			  + sf3.dN1_dx * sf3.N2 + sf4.dN1_dx * sf4.N2) * det_dx_dxi;
	k_mat[0][9] = mat_term;
	k_mat[9][0] = mat_term;
	mat_term = (sf1.dN1_dx * sf1.N3 + sf2.dN1_dx * sf2.N3
			  + sf3.dN1_dx * sf3.N3 + sf4.dN1_dx * sf4.N3) * det_dx_dxi;
	k_mat[0][10] = mat_term;
	k_mat[10][0] = mat_term;
	mat_term = (sf1.dN1_dx * sf1.N4 + sf2.dN1_dx * sf2.N4
			  + sf3.dN1_dx * sf3.N4 + sf4.dN1_dx * sf4.N4) * det_dx_dxi;
	k_mat[0][11] = mat_term;
	k_mat[11][0] = mat_term;
	mat_term = (sf1.dN2_dx * sf1.N1 + sf2.dN2_dx * sf2.N1
			  + sf3.dN2_dx * sf3.N1 + sf4.dN2_dx * sf4.N1) * det_dx_dxi;
	k_mat[1][8] = mat_term;
	k_mat[8][1] = mat_term;
	mat_term = (sf1.dN2_dx * sf1.N2 + sf2.dN2_dx * sf2.N2
			  + sf3.dN2_dx * sf3.N2 + sf4.dN2_dx * sf4.N2) * det_dx_dxi;
	k_mat[1][9] = mat_term;
	k_mat[9][1] = mat_term;
	mat_term = (sf1.dN2_dx * sf1.N3 + sf2.dN2_dx * sf2.N3
			  + sf3.dN2_dx * sf3.N3 + sf4.dN2_dx * sf4.N3) * det_dx_dxi;
	k_mat[1][10] = mat_term;
	k_mat[10][1] = mat_term;
	mat_term = (sf1.dN2_dx * sf1.N4 + sf2.dN2_dx * sf2.N4
			  + sf3.dN2_dx * sf3.N4 + sf4.dN2_dx * sf4.N4) * det_dx_dxi;
	k_mat[1][11] = mat_term;
	k_mat[11][1] = mat_term;
	mat_term = (sf1.dN3_dx * sf1.N1 + sf2.dN3_dx * sf2.N1
			  + sf3.dN3_dx * sf3.N1 + sf4.dN3_dx * sf4.N1) * det_dx_dxi;
	k_mat[2][8] = mat_term;
	k_mat[8][2] = mat_term;
	mat_term = (sf1.dN3_dx * sf1.N2 + sf2.dN3_dx * sf2.N2
			  + sf3.dN3_dx * sf3.N2 + sf4.dN3_dx * sf4.N2) * det_dx_dxi;
	k_mat[2][9] = mat_term;
	k_mat[9][2] = mat_term;
	mat_term = (sf1.dN3_dx * sf1.N3 + sf2.dN3_dx * sf2.N3
			  + sf3.dN3_dx * sf3.N3 + sf4.dN3_dx * sf4.N3) * det_dx_dxi;
	k_mat[2][10] = mat_term;
	k_mat[10][2] = mat_term;
	mat_term = (sf1.dN3_dx * sf1.N4 + sf2.dN3_dx * sf2.N4
			  + sf3.dN3_dx * sf3.N4 + sf4.dN3_dx * sf4.N4) * det_dx_dxi;
	k_mat[2][11] = mat_term;
	k_mat[11][2] = mat_term;
	mat_term = (sf1.dN4_dx * sf1.N1 + sf2.dN4_dx * sf2.N1
			  + sf3.dN4_dx * sf3.N1 + sf4.dN4_dx * sf4.N1) * det_dx_dxi;
	k_mat[3][8] = mat_term;
	k_mat[8][3] = mat_term;
	mat_term = (sf1.dN4_dx * sf1.N2 + sf2.dN4_dx * sf2.N2
			  + sf3.dN4_dx * sf3.N2 + sf4.dN4_dx * sf4.N2) * det_dx_dxi;
	k_mat[3][9] = mat_term;
	k_mat[9][3] = mat_term;
	mat_term = (sf1.dN4_dx * sf1.N3 + sf2.dN4_dx * sf2.N3
			  + sf3.dN4_dx * sf3.N3 + sf4.dN4_dx * sf4.N3) * det_dx_dxi;
	k_mat[3][10] = mat_term;
	k_mat[10][3] = mat_term;
	mat_term = (sf1.dN4_dx * sf1.N4 + sf2.dN4_dx * sf2.N4
			  + sf3.dN4_dx * sf3.N4 + sf4.dN4_dx * sf4.N4) * det_dx_dxi;
	k_mat[3][11] = mat_term;
	k_mat[11][3] = mat_term;
	// dNi_dy * Nj
	mat_term = (sf1.dN1_dy * sf1.N1 + sf2.dN1_dy * sf2.N1
			  + sf3.dN1_dy * sf3.N1 + sf4.dN1_dy * sf4.N1) * det_dx_dxi;
	k_mat[4][8] = mat_term;
	k_mat[8][4] = mat_term;
	mat_term = (sf1.dN1_dy * sf1.N2 + sf2.dN1_dy * sf2.N2
			  + sf3.dN1_dy * sf3.N2 + sf4.dN1_dy * sf4.N2) * det_dx_dxi;
	k_mat[4][9] = mat_term;
	k_mat[9][4] = mat_term;
	mat_term = (sf1.dN1_dy * sf1.N3 + sf2.dN1_dy * sf2.N3
			  + sf3.dN1_dy * sf3.N3 + sf4.dN1_dy * sf4.N3) * det_dx_dxi;
	k_mat[4][10] = mat_term;
	k_mat[10][4] = mat_term;
	mat_term = (sf1.dN1_dy * sf1.N4 + sf2.dN1_dy * sf2.N4
			  + sf3.dN1_dy * sf3.N4 + sf4.dN1_dy * sf4.N4) * det_dx_dxi;
	k_mat[4][11] = mat_term;
	k_mat[11][4] = mat_term;
	mat_term = (sf1.dN2_dy * sf1.N1 + sf2.dN2_dy * sf2.N1
			  + sf3.dN2_dy * sf3.N1 + sf4.dN2_dy * sf4.N1) * det_dx_dxi;
	k_mat[5][8] = mat_term;
	k_mat[8][5] = mat_term;
	mat_term = (sf1.dN2_dy * sf1.N2 + sf2.dN2_dy * sf2.N2
			  + sf3.dN2_dy * sf3.N2 + sf4.dN2_dy * sf4.N2) * det_dx_dxi;
	k_mat[5][9] = mat_term;
	k_mat[9][5] = mat_term;
	mat_term = (sf1.dN2_dy * sf1.N3 + sf2.dN2_dy * sf2.N3
			  + sf3.dN2_dy * sf3.N3 + sf4.dN2_dy * sf4.N3) * det_dx_dxi;
	k_mat[5][10] = mat_term;
	k_mat[10][5] = mat_term;
	mat_term = (sf1.dN2_dy * sf1.N4 + sf2.dN2_dy * sf2.N4
			  + sf3.dN2_dy * sf3.N4 + sf4.dN2_dy * sf4.N4) * det_dx_dxi;
	k_mat[5][11] = mat_term;
	k_mat[11][5] = mat_term;
	mat_term = (sf1.dN3_dy * sf1.N1 + sf2.dN3_dy * sf2.N1
			  + sf3.dN3_dy * sf3.N1 + sf4.dN3_dy * sf4.N1) * det_dx_dxi;
	k_mat[6][8] = mat_term;
	k_mat[8][6] = mat_term;
	mat_term = (sf1.dN3_dy * sf1.N2 + sf2.dN3_dy * sf2.N2
			  + sf3.dN3_dy * sf3.N2 + sf4.dN3_dy * sf4.N2) * det_dx_dxi;
	k_mat[6][9] = mat_term;
	k_mat[9][6] = mat_term;
	mat_term = (sf1.dN3_dy * sf1.N3 + sf2.dN3_dy * sf2.N3
			  + sf3.dN3_dy * sf3.N3 + sf4.dN3_dy * sf4.N3) * det_dx_dxi;
	k_mat[6][10] = mat_term;
	k_mat[10][6] = mat_term;
	mat_term = (sf1.dN3_dy * sf1.N4 + sf2.dN3_dy * sf2.N4
			  + sf3.dN3_dy * sf3.N4 + sf4.dN3_dy * sf4.N4) * det_dx_dxi;
	k_mat[6][11] = mat_term;
	k_mat[11][6] = mat_term;
	mat_term = (sf1.dN4_dy * sf1.N1 + sf2.dN4_dy * sf2.N1
			  + sf3.dN4_dy * sf3.N1 + sf4.dN4_dy * sf4.N1) * det_dx_dxi;
	k_mat[7][8] = mat_term;
	k_mat[8][7] = mat_term;
	mat_term = (sf1.dN4_dy * sf1.N2 + sf2.dN4_dy * sf2.N2
			  + sf3.dN4_dy * sf3.N2 + sf4.dN4_dy * sf4.N2) * det_dx_dxi;
	k_mat[7][9] = mat_term;
	k_mat[9][7] = mat_term;
	mat_term = (sf1.dN4_dy * sf1.N3 + sf2.dN4_dy * sf2.N3
			  + sf3.dN4_dy * sf3.N3 + sf4.dN4_dy * sf4.N3) * det_dx_dxi;
	k_mat[7][10] = mat_term;
	k_mat[10][7] = mat_term;
	mat_term = (sf1.dN4_dy * sf1.N4 + sf2.dN4_dy * sf2.N4
			  + sf3.dN4_dy * sf3.N4 + sf4.dN4_dy * sf4.N4) * det_dx_dxi;
	k_mat[7][11] = mat_term;
	k_mat[11][7] = mat_term;
	// Ni * Nj / K
	mat_term = (sf1.N1 * sf1.N1 / gp1.K + sf2.N1 * sf2.N1 / gp2.K
			  + sf3.N1 * sf3.N1 / gp3.K + sf4.N1 * sf4.N1 / gp4.K) * -det_dx_dxi;
	k_mat[8][8] = mat_term;
	mat_term = (sf1.N1 * sf1.N2 / gp1.K + sf2.N1 * sf2.N2 / gp2.K
			  + sf3.N1 * sf3.N2 / gp3.K + sf4.N1 * sf4.N2 / gp4.K) * -det_dx_dxi;
	k_mat[8][9] = mat_term;
	k_mat[9][8] = mat_term;
	mat_term = (sf1.N1 * sf1.N3 / gp1.K + sf2.N1 * sf2.N3 / gp2.K
			  + sf3.N1 * sf3.N3 / gp3.K + sf4.N1 * sf4.N3 / gp4.K) * -det_dx_dxi;
	k_mat[8][10] = mat_term;
	k_mat[10][8] = mat_term;
	mat_term = (sf1.N1 * sf1.N4 / gp1.K + sf2.N1 * sf2.N4 / gp2.K
			  + sf3.N1 * sf3.N4 / gp3.K + sf4.N1 * sf4.N4 / gp4.K) * -det_dx_dxi;
	k_mat[8][11] = mat_term;
	k_mat[11][8] = mat_term;
	mat_term = (sf1.N2 * sf1.N2 / gp1.K + sf2.N2 * sf2.N2 / gp2.K
			  + sf3.N2 * sf3.N2 / gp3.K + sf4.N2 * sf4.N2 / gp4.K) * -det_dx_dxi;
	k_mat[9][9] = mat_term;
	mat_term = (sf1.N2 * sf1.N3 / gp1.K + sf2.N2 * sf2.N3 / gp2.K
			  + sf3.N2 * sf3.N3 / gp3.K + sf4.N2 * sf4.N3 / gp4.K) * -det_dx_dxi;
	k_mat[9][10] = mat_term;
	k_mat[10][9] = mat_term;
	mat_term = (sf1.N2 * sf1.N4 / gp1.K + sf2.N2 * sf2.N4 / gp2.K
			  + sf3.N2 * sf3.N4 / gp3.K + sf4.N2 * sf4.N4 / gp4.K) * -det_dx_dxi;
	k_mat[9][11] = mat_term;
	k_mat[11][9] = mat_term;
	mat_term = (sf1.N3 * sf1.N3 / gp1.K + sf2.N3 * sf2.N3 / gp2.K
			  + sf3.N3 * sf3.N3 / gp3.K + sf4.N3 * sf4.N3 / gp4.K) * -det_dx_dxi;
	k_mat[10][10] = mat_term;
	mat_term = (sf1.N3 * sf1.N4 / gp1.K + sf2.N3 * sf2.N4 / gp2.K
			  + sf3.N3 * sf3.N4 / gp3.K + sf4.N3 * sf4.N4 / gp4.K) * -det_dx_dxi;
	k_mat[10][11] = mat_term;
	k_mat[11][10] = mat_term;
	mat_term = (sf1.N4 * sf1.N4 / gp1.K + sf2.N4 * sf2.N4 / gp2.K
			  + sf3.N4 * sf3.N4 / gp3.K + sf4.N4 * sf4.N4 / gp4.K) * -det_dx_dxi;
	k_mat[11][11] = mat_term;

	//print_mat(k_mat);

	// force vector
	f_a[0] = (sf1.N1 * gp1.density * gp1.ax + sf2.N1 * gp2.density * gp2.ax
			+ sf3.N1 * gp3.density * gp3.ax + sf4.N1 * gp4.density * gp4.ax) * det_dx_dxi;
	f_a[1] = (sf1.N2 * gp1.density * gp1.ax + sf2.N2 * gp2.density * gp2.ax
			+ sf3.N2 * gp3.density * gp3.ax + sf4.N2 * gp4.density * gp4.ax) * det_dx_dxi;
	f_a[2] = (sf1.N3 * gp1.density * gp1.ax + sf2.N3 * gp2.density * gp2.ax
			+ sf3.N3 * gp3.density * gp3.ax + sf4.N3 * gp4.density * gp4.ax) * det_dx_dxi;
	f_a[3] = (sf1.N4 * gp1.density * gp1.ax + sf2.N4 * gp2.density * gp2.ax
			+ sf3.N4 * gp3.density * gp3.ax + sf4.N4 * gp4.density * gp4.ax) * det_dx_dxi;
	f_a[4] = (sf1.N1 * gp1.density * gp1.ay + sf2.N1 * gp2.density * gp2.ay
			+ sf3.N1 * gp3.density * gp3.ay + sf4.N1 * gp4.density * gp4.ay) * det_dx_dxi;
	f_a[5] = (sf1.N2 * gp1.density * gp1.ay + sf2.N2 * gp2.density * gp2.ay
			+ sf3.N2 * gp3.density * gp3.ay + sf4.N2 * gp4.density * gp4.ay) * det_dx_dxi;
	f_a[6] = (sf1.N3 * gp1.density * gp1.ay + sf2.N3 * gp2.density * gp2.ay
			+ sf3.N3 * gp3.density * gp3.ay + sf4.N3 * gp4.density * gp4.ay) * det_dx_dxi;
	f_a[7] = (sf1.N4 * gp1.density * gp1.ay + sf2.N4 * gp2.density * gp2.ay
			+ sf3.N4 * gp3.density * gp3.ay + sf4.N4 * gp4.density * gp4.ay) * det_dx_dxi;

	f_v[0] = (sf1.N1 * gp1.density * gp1.vx + sf2.N1 * gp2.density * gp2.vx
			+ sf3.N1 * gp3.density * gp3.vx + sf4.N1 * gp4.density * gp4.vx) * det_dx_dxi;
	f_v[1] = (sf1.N2 * gp1.density * gp1.vx + sf2.N2 * gp2.density * gp2.vx
			+ sf3.N2 * gp3.density * gp3.vx + sf4.N2 * gp4.density * gp4.vx) * det_dx_dxi;
	f_v[2] = (sf1.N3 * gp1.density * gp1.vx + sf2.N3 * gp2.density * gp2.vx
			+ sf3.N3 * gp3.density * gp3.vx + sf4.N3 * gp4.density * gp4.vx) * det_dx_dxi;
	f_v[3] = (sf1.N4 * gp1.density * gp1.vx + sf2.N4 * gp2.density * gp2.vx
			+ sf3.N4 * gp3.density * gp3.vx + sf4.N4 * gp4.density * gp4.vx) * det_dx_dxi;
	f_v[4] = (sf1.N1 * gp1.density * gp1.vy + sf2.N1 * gp2.density * gp2.vy
			+ sf3.N1 * gp3.density * gp3.vy + sf4.N1 * gp4.density * gp4.vy) * det_dx_dxi;
	f_v[5] = (sf1.N2 * gp1.density * gp1.vy + sf2.N2 * gp2.density * gp2.vy
			+ sf3.N2 * gp3.density * gp3.vy + sf4.N2 * gp4.density * gp4.vy) * det_dx_dxi;
	f_v[6] = (sf1.N3 * gp1.density * gp1.vy + sf2.N3 * gp2.density * gp2.vy
			+ sf3.N3 * gp3.density * gp3.vy + sf4.N3 * gp4.density * gp4.vy) * det_dx_dxi;
	f_v[7] = (sf1.N4 * gp1.density * gp1.vy + sf2.N4 * gp2.density * gp2.vy
			+ sf3.N4 * gp3.density * gp3.vy + sf4.N4 * gp4.density * gp4.vy) * det_dx_dxi;

	for (size_t i = 0; i < 12; ++i)
		for (size_t j = 0; j < 12; ++j)
			kmat[i][j] = m_mat[i][j] / (beta * dtime * dtime) + k_mat[i][j];

	for (size_t i = 0; i < 8; ++i)
		fvec[i] = f_a[i] / (2.0 * beta) + f_v[i] / (beta * dtime);
	fvec[8] = 0.0;
	fvec[9] = 0.0;
	fvec[10] = 0.0;
	fvec[11] = 0.0;
}
