#include "SimulationCore_pcp.h"

#include <Eigen/SparseLU>

#define Keep_Newmark_Coefficients
#include "Step_S2D_CHM_s_FEM_uUp.h"

namespace
{
	typedef Model_S2D_CHM_s_FEM_uUp::ShapeFuncValue ShapeFuncValue;
	typedef Model_S2D_CHM_s_FEM_uUp::GaussPoint GaussPoint_fem;
	typedef Model_S2D_CHM_s_FEM_uUp::Element Element_fem;
	typedef Model_S2D_CHM_s_FEM_uUp::Node Node_fem;
	typedef Model_S2D_CHM_s_FEM_uUp::DOF DOF;
	
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
};

Step_S2D_CHM_s_FEM_uUp::Step_S2D_CHM_s_FEM_uUp() :
	Step(&solve_substep_S2D_CHM_s_FEM_uUp),
	model(nullptr), kmat_col(nullptr) {}

Step_S2D_CHM_s_FEM_uUp::~Step_S2D_CHM_s_FEM_uUp() {}

int Step_S2D_CHM_s_FEM_uUp::init_calculation(void)
{
	if (is_first_step) {}

	kmat_col = new double[model->dof_num];

	// for debug
	out_file.open("debug_mat_out_CHM_s_FEM_uUp.csv", std::ios::binary | std::ios::out);
	
	return 0;
}

int Step_S2D_CHM_s_FEM_uUp::finalize_calculation(void)
{
	if (kmat_col)
	{
		delete[] kmat_col;
		kmat_col = nullptr;
	}
	
	// for debug
	out_file.close();

	return 0;
}

int solve_substep_S2D_CHM_s_FEM_uUp(void *_self)
{
	// for debug
	//static std::fstream out_f2("out_f2.txt", std::ios::binary | std::ios::out);

	Step_S2D_CHM_s_FEM_uUp &self = *(Step_S2D_CHM_s_FEM_uUp *)(_self);
	Model_S2D_CHM_s_FEM_uUp &model = *self.model;
	double dt = self.dtime;
	double dt2 = dt * dt;

	// list of non-zeros coefficients
	MatrixCoefficientSet<> &g_kmat_coefs = self.g_kmat_coefs;
	g_kmat_coefs.init(model.dof_num);
	Eigen::VectorXd g_fvec(model.dof_num);
	g_fvec.setZero();

	size_t l2g_dof_id_map[20];
	double e_kmat[20][20], e_fvec[20];
	for (size_t e_id = 0; e_id < model.elem_num; ++e_id)
	{
		Element_fem &e = model.elems[e_id];
		// form elemental stiffness matrix and force vector
		self.form_elem_stiffness_mat_and_force_vec(e, e_kmat, e_fvec);
		
		// map from local dof id to global dof id
		// solid phase ux
		l2g_dof_id_map[0] = model.n_id_to_dof_id(e.n1_id, DOF::usx);
		l2g_dof_id_map[1] = model.n_id_to_dof_id(e.n2_id, DOF::usx);
		l2g_dof_id_map[2] = model.n_id_to_dof_id(e.n3_id, DOF::usx);
		l2g_dof_id_map[3] = model.n_id_to_dof_id(e.n4_id, DOF::usx);
		// solid phase uy
		l2g_dof_id_map[4] = model.n_id_to_dof_id(e.n1_id, DOF::usy);
		l2g_dof_id_map[5] = model.n_id_to_dof_id(e.n2_id, DOF::usy);
		l2g_dof_id_map[6] = model.n_id_to_dof_id(e.n3_id, DOF::usy);
		l2g_dof_id_map[7] = model.n_id_to_dof_id(e.n4_id, DOF::usy);
		// fluid phase ux
		l2g_dof_id_map[8]  = model.n_id_to_dof_id(e.n1_id, DOF::ufx);
		l2g_dof_id_map[9]  = model.n_id_to_dof_id(e.n2_id, DOF::ufx);
		l2g_dof_id_map[10] = model.n_id_to_dof_id(e.n3_id, DOF::ufx);
		l2g_dof_id_map[11] = model.n_id_to_dof_id(e.n4_id, DOF::ufx);
		// fluid phase uy
		l2g_dof_id_map[12] = model.n_id_to_dof_id(e.n1_id, DOF::ufy);
		l2g_dof_id_map[13] = model.n_id_to_dof_id(e.n2_id, DOF::ufy);
		l2g_dof_id_map[14] = model.n_id_to_dof_id(e.n3_id, DOF::ufy);
		l2g_dof_id_map[15] = model.n_id_to_dof_id(e.n4_id, DOF::ufy);
		// p
		l2g_dof_id_map[16] = model.n_id_to_dof_id(e.n1_id, DOF::p);
		l2g_dof_id_map[17] = model.n_id_to_dof_id(e.n2_id, DOF::p);
		l2g_dof_id_map[18] = model.n_id_to_dof_id(e.n3_id, DOF::p);
		l2g_dof_id_map[19] = model.n_id_to_dof_id(e.n4_id, DOF::p);
		//print_vec(l2g_dof_id_map);
		// add to global matrix and vector
		size_t g_dof_id1, g_dof_id2;
		for (size_t l_id1 = 0; l_id1 < 20; ++l_id1)
		{
			// add to global matrix
			g_dof_id1 = l2g_dof_id_map[l_id1];
			for (size_t l_id2 = 0; l_id2 < 20; ++l_id2)
			{
				g_dof_id2 = l2g_dof_id_map[l_id2];
				g_kmat_coefs.add_coefficient(g_dof_id1, g_dof_id2, e_kmat[l_id1][l_id2]);
			}
			// add to global force vector
			g_fvec[g_dof_id1] += e_fvec[l_id1];
		}
	}

	// apply external force body force and traction
	size_t dof_id;
	// traction
	double tfs[4];
	for (size_t t_id = 0; t_id < model.tx_num; ++t_id)
	{
		TractionBC_2DFEM &tx = model.txs[t_id];
		model.cal_traction_bc(tx, tfs);
		Element_fem &e = model.elems[tx.elem_id];
		dof_id = model.n_id_to_dof_id(e.n1_id, DOF::usx);
		g_fvec[dof_id] += tfs[0];
		dof_id = model.n_id_to_dof_id(e.n2_id, DOF::usx);
		g_fvec[dof_id] += tfs[1];
		dof_id = model.n_id_to_dof_id(e.n3_id, DOF::usx);
		g_fvec[dof_id] += tfs[2];
		dof_id = model.n_id_to_dof_id(e.n4_id, DOF::usx);
		g_fvec[dof_id] += tfs[3];
	}
	for (size_t t_id = 0; t_id < model.ty_num; ++t_id)
	{
		TractionBC_2DFEM &ty = model.tys[t_id];
		model.cal_traction_bc(ty, tfs);
		Element_fem &e = model.elems[ty.elem_id];
		dof_id = model.n_id_to_dof_id(e.n1_id, DOF::usy);
		g_fvec[dof_id] += tfs[0];
		dof_id = model.n_id_to_dof_id(e.n2_id, DOF::usy);
		g_fvec[dof_id] += tfs[1];
		dof_id = model.n_id_to_dof_id(e.n3_id, DOF::usy);
		g_fvec[dof_id] += tfs[2];
		dof_id = model.n_id_to_dof_id(e.n4_id, DOF::usy);
		g_fvec[dof_id] += tfs[3];
	}
	// body force ...  to be finished
	
	// apply displacement boundary condition
	double dig_term;
	for (size_t bc_id = 0; bc_id < model.usx_num; ++bc_id)
	{
		DisplacementBC &dbc = model.usxs[bc_id];
		dof_id = model.n_id_to_dof_id(dbc.node_id, DOF::usx);
		dig_term = g_kmat_coefs.del_col_and_row(dof_id, self.kmat_col);
		for (size_t row_id = 0; row_id < model.dof_num; ++row_id)
			g_fvec[row_id] -= self.kmat_col[row_id] * dbc.u;
		g_fvec[dof_id] = dig_term * dbc.u;
	}
	for (size_t bc_id = 0; bc_id < model.usy_num; ++bc_id)
	{
		DisplacementBC &dbc = model.usys[bc_id];
		dof_id = model.n_id_to_dof_id(dbc.node_id, DOF::usy);
		dig_term = g_kmat_coefs.del_col_and_row(dof_id, self.kmat_col);
		for (size_t row_id = 0; row_id < model.dof_num; ++row_id)
			g_fvec[row_id] -= self.kmat_col[row_id] * dbc.u;
		g_fvec[dof_id] = dig_term * dbc.u;
	}
	for (size_t bc_id = 0; bc_id < model.ufx_num; ++bc_id)
	{
		DisplacementBC &dbc = model.ufxs[bc_id];
		dof_id = model.n_id_to_dof_id(dbc.node_id, DOF::ufx);
		dig_term = g_kmat_coefs.del_col_and_row(dof_id, self.kmat_col);
		for (size_t row_id = 0; row_id < model.dof_num; ++row_id)
			g_fvec[row_id] -= self.kmat_col[row_id] * dbc.u;
		g_fvec[dof_id] = dig_term * dbc.u;
	}
	for (size_t bc_id = 0; bc_id < model.ufy_num; ++bc_id)
	{
		DisplacementBC &dbc = model.ufys[bc_id];
		dof_id = model.n_id_to_dof_id(dbc.node_id, DOF::ufy);
		dig_term = g_kmat_coefs.del_col_and_row(dof_id, self.kmat_col);
		for (size_t row_id = 0; row_id < model.dof_num; ++row_id)
			g_fvec[row_id] -= self.kmat_col[row_id] * dbc.u;
		g_fvec[dof_id] = dig_term * dbc.u;
	}
	for (size_t bc_id = 0; bc_id < model.pbc_num; ++bc_id)
	{
		PressureBC &pbc = model.pbcs[bc_id];
		dof_id = model.n_id_to_dof_id(pbc.node_id, DOF::p);
		dig_term = g_kmat_coefs.del_col_and_row(dof_id, self.kmat_col);
		for (size_t row_id = 0; row_id < model.dof_num; ++row_id)
			g_fvec[row_id] -= self.kmat_col[row_id] * pbc.p;
		g_fvec[dof_id] = dig_term * pbc.p;
	}
	
	// solve
	//g_kmat_coefs.print_with_iter();
	Eigen::SparseMatrix<double> g_kmat(model.dof_num, model.dof_num);
	g_kmat.setFromTriplets(g_kmat_coefs.begin(), g_kmat_coefs.end());
	//print_sparse_mat(g_kmat, self.out_file,nullptr);
	Eigen::SparseLU<Eigen::SparseMatrix<double> > solver(g_kmat);
	Eigen::VectorXd g_du_vec = solver.solve(g_fvec);
	//std::cout << g_fvec << "\n";
	//std::cout << g_du_vec << "\n";

	// reapply disp bc
	for (size_t bc_id = 0; bc_id < model.usx_num; ++bc_id)
	{
		DisplacementBC &dbc = model.usxs[bc_id];
		dof_id = model.n_id_to_dof_id(dbc.node_id, DOF::usx);
		g_du_vec[dof_id] = dbc.u;
	}
	for (size_t bc_id = 0; bc_id < model.usy_num; ++bc_id)
	{
		DisplacementBC &dbc = model.usys[bc_id];
		dof_id = model.n_id_to_dof_id(dbc.node_id, DOF::usy);
		g_du_vec[dof_id] = dbc.u;
	}
	for (size_t bc_id = 0; bc_id < model.ufx_num; ++bc_id)
	{
		DisplacementBC &dbc = model.ufxs[bc_id];
		dof_id = model.n_id_to_dof_id(dbc.node_id, DOF::ufx);
		g_du_vec[dof_id] = dbc.u;
	}
	for (size_t bc_id = 0; bc_id < model.ufy_num; ++bc_id)
	{
		DisplacementBC &dbc = model.ufys[bc_id];
		dof_id = model.n_id_to_dof_id(dbc.node_id, DOF::ufy);
		g_du_vec[dof_id] = dbc.u;
	}
	for (size_t bc_id = 0; bc_id < model.pbc_num; ++bc_id)
	{
		PressureBC &pbc = model.pbcs[bc_id];
		dof_id = model.n_id_to_dof_id(pbc.node_id, DOF::p);
		g_du_vec[dof_id] = pbc.p;
	}
	
	// update nodal variables
	double da, dv, du;
	for (size_t n_id = 0; n_id < model.node_num; ++n_id)
	{
		Node_fem &n = model.nodes[n_id];
		// sx
		dof_id = model.n_id_to_dof_id(n_id, DOF::usx);
		du = g_du_vec[dof_id];
		n.ux_s += du;
		dv = gamma / (beta * dt) * du - gamma / beta * n.vx_s + dt * (1.0 - gamma / (2.0 * beta)) * n.ax_s;
		n.vx_s += dv;
		da = du / (gamma * dt2) - n.vx_s / (gamma * dt) - n.ax_s / (2.0 * gamma);
		n.ax_s += da;
		// sy
		dof_id = model.n_id_to_dof_id(n_id, DOF::usy);
		du = g_du_vec[dof_id];
		n.uy_s += du;
		dv = gamma / (beta * dt) * du - gamma / beta * n.vy_s + dt * (1.0 - gamma / (2.0 * beta)) * n.ay_s;
		n.vy_s += dv;
		da = du / (gamma * dt2) -  n.vy_s / (gamma * dt) - n.ay_s / (2.0 * gamma);
		n.ay_s += da;
		// fx
		dof_id = model.n_id_to_dof_id(n_id, DOF::ufx);
		du = g_du_vec[dof_id];
		n.ux_f += du;
		dv = gamma / (beta * dt) * du - gamma / beta * n.vx_f + dt * (1.0 - gamma / (2.0 * beta)) * n.ax_f;
		n.vx_f += dv;
		da = 1.0 / (gamma * dt2) * du - 1.0 / (gamma * dt) * n.vx_f - 1.0 / (2.0 * gamma) * n.ax_f;
		n.ax_f += da;
		// fy
		dof_id = model.n_id_to_dof_id(n_id, DOF::ufy);
		du = g_du_vec[dof_id];
		n.uy_f += du;
		dv = gamma / (beta * dt) * du - gamma / beta * n.vy_f + dt * (1.0 - gamma / (2.0 * beta)) * n.ay_s;
		n.vy_f += dv;
		da = 1.0 / (gamma * dt2) * du - 1.0 / (gamma * dt) * n.vy_f - 1.0 / (2.0 * gamma) * n.ay_s;
		n.ay_f += da;
		// p
		dof_id = model.n_id_to_dof_id(n_id, DOF::p);
		n.p += g_du_vec[dof_id];
	}

	// Update element variables
	ShapeFuncValue &sf1 = model.gp1_sf;
	ShapeFuncValue &sf2 = model.gp2_sf;
	ShapeFuncValue &sf3 = model.gp3_sf;
	ShapeFuncValue &sf4 = model.gp4_sf;
	double dux_s1, duy_s1, dux_f1, duy_f1, dp1;
	double dux_s2, duy_s2, dux_f2, duy_f2, dp2;
	double dux_s3, duy_s3, dux_f3, duy_f3, dp3;
	double dux_s4, duy_s4, dux_f4, duy_f4, dp4;
	for (size_t e_id = 0; e_id < model.elem_num; ++e_id)
	{
		Element_fem &e = model.elems[e_id];
		dux_s1 = g_du_vec[model.n_id_to_dof_id(e.n1_id, DOF::usx)];
		duy_s1 = g_du_vec[model.n_id_to_dof_id(e.n1_id, DOF::usy)];
		dux_f1 = g_du_vec[model.n_id_to_dof_id(e.n1_id, DOF::ufx)];
		duy_f1 = g_du_vec[model.n_id_to_dof_id(e.n1_id, DOF::ufy)];
		dp1    = g_du_vec[model.n_id_to_dof_id(e.n1_id, DOF::p)];
		dux_s2 = g_du_vec[model.n_id_to_dof_id(e.n2_id, DOF::usx)];
		duy_s2 = g_du_vec[model.n_id_to_dof_id(e.n2_id, DOF::usy)];
		dux_f2 = g_du_vec[model.n_id_to_dof_id(e.n2_id, DOF::ufx)];
		duy_f2 = g_du_vec[model.n_id_to_dof_id(e.n2_id, DOF::ufy)];
		dp2    = g_du_vec[model.n_id_to_dof_id(e.n2_id, DOF::p)];
		dux_s3 = g_du_vec[model.n_id_to_dof_id(e.n3_id, DOF::usx)];
		duy_s3 = g_du_vec[model.n_id_to_dof_id(e.n3_id, DOF::usy)];
		dux_f3 = g_du_vec[model.n_id_to_dof_id(e.n3_id, DOF::ufx)];
		duy_f3 = g_du_vec[model.n_id_to_dof_id(e.n3_id, DOF::ufy)];
		dp3    = g_du_vec[model.n_id_to_dof_id(e.n3_id, DOF::p)];
		dux_s4 = g_du_vec[model.n_id_to_dof_id(e.n4_id, DOF::usx)];
		duy_s4 = g_du_vec[model.n_id_to_dof_id(e.n4_id, DOF::usy)];
		dux_f4 = g_du_vec[model.n_id_to_dof_id(e.n4_id, DOF::ufx)];
		duy_f4 = g_du_vec[model.n_id_to_dof_id(e.n4_id, DOF::ufy)];
		dp4    = g_du_vec[model.n_id_to_dof_id(e.n4_id, DOF::p)];

		self.update_gauss_point(e.gp1, sf1, dt, dt2,
			dux_s1, dux_s2, dux_s3, dux_s4,
			duy_s1, duy_s2, duy_s3, duy_s4,
			dux_f1, dux_f2, dux_f3, dux_f4,
			duy_f1, duy_f2, duy_f3, duy_f4,
			dp1, dp2, dp3, dp4);

		self.update_gauss_point(e.gp2, sf2, dt, dt2,
			dux_s1, dux_s2, dux_s3, dux_s4,
			duy_s1, duy_s2, duy_s3, duy_s4,
			dux_f1, dux_f2, dux_f3, dux_f4,
			duy_f1, duy_f2, duy_f3, duy_f4,
			dp1, dp2, dp3, dp4);

		self.update_gauss_point(e.gp3, sf3, dt, dt2,
			dux_s1, dux_s2, dux_s3, dux_s4,
			duy_s1, duy_s2, duy_s3, duy_s4,
			dux_f1, dux_f2, dux_f3, dux_f4,
			duy_f1, duy_f2, duy_f3, duy_f4,
			dp1, dp2, dp3, dp4);

		self.update_gauss_point(e.gp4, sf4, dt, dt2,
			dux_s1, dux_s2, dux_s3, dux_s4,
			duy_s1, duy_s2, duy_s3, duy_s4,
			dux_f1, dux_f2, dux_f3, dux_f4,
			duy_f1, duy_f2, duy_f3, duy_f4,
			dp1, dp2, dp3, dp4);
	}
	
	return 0;
}

void Step_S2D_CHM_s_FEM_uUp::update_gauss_point(
	Model_S2D_CHM_s_FEM_uUp::GaussPoint &gp,
	Model_S2D_CHM_s_FEM_uUp::ShapeFuncValue &sf,
	double dt, double dt2,
	double dux_s1, double dux_s2, double dux_s3, double dux_s4,
	double duy_s1, double duy_s2, double duy_s3, double duy_s4,
	double dux_f1, double dux_f2, double dux_f3, double dux_f4,
	double duy_f1, double duy_f2, double duy_f3, double duy_f4,
	double dp1, double dp2, double dp3, double dp4)
{
	double dux_s, duy_s, dux_f, duy_f, dp;
	double dvx_s, dvy_s, dvx_f, dvy_f;
	double dax_s, day_s, dax_f, day_f;
	double de11, de22, de12, ds11, ds22, ds12;
	double de_vol_s, de_vol_f;

	// displacement
	dux_s = dux_s1 * sf.N1 + dux_s2 * sf.N2 + dux_s3 * sf.N3 + dux_s4 * sf.N4;
	duy_s = duy_s1 * sf.N1 + duy_s2 * sf.N2 + duy_s3 * sf.N3 + duy_s4 * sf.N4;
	dux_f = dux_f1 * sf.N1 + dux_f2 * sf.N2 + dux_f3 * sf.N3 + dux_f4 * sf.N4;
	duy_f = duy_f1 * sf.N1 + duy_f2 * sf.N2 + duy_f3 * sf.N3 + duy_f4 * sf.N4;
	// velocity
	dvx_s = gamma / (beta * dt) * dux_s - gamma / beta * gp.vx_s
		  + dt * (1.0 - gamma / (2.0 * beta)) * gp.ax_s;
	dvy_s = gamma / (beta * dt) * duy_s - gamma / beta * gp.vy_s
		  + dt * (1.0 - gamma / (2.0 * beta)) * gp.ay_s;
	dvx_f = gamma / (beta * dt) * dux_f - gamma / beta * gp.vx_f
		  + dt * (1.0 - gamma / (2.0 * beta)) * gp.ax_f;
	dvy_f = gamma / (beta * dt) * duy_f - gamma / beta * gp.vy_f
		  + dt * (1.0 - gamma / (2.0 * beta)) * gp.ay_f;
	// acceleration
	dax_s = dux_s / (beta * dt2) - gp.vx_s / (beta * dt) - gp.ax_s / (2.0 * beta);
	day_s = duy_s / (beta * dt2) - gp.vy_s / (beta * dt) - gp.ay_s / (2.0 * beta);
	dax_f = dux_f / (beta * dt2) - gp.vx_f / (beta * dt) - gp.ax_f / (2.0 * beta);
	day_f = duy_f / (beta * dt2) - gp.vy_f / (beta * dt) - gp.ay_f / (2.0 * beta);
	gp.ux_s += dux_s;
	gp.uy_s += duy_s;
	gp.ux_f += dux_f;
	gp.uy_f += duy_f;
	gp.vx_s += dvx_s;
	gp.vy_s += dvy_s;
	gp.vx_f += dvx_f;
	gp.vy_f += dvy_f;
	gp.ax_s += dax_s;
	gp.ay_s += day_s;
	gp.ax_f += dax_f;
	gp.ay_f += day_f;

	// pore pressure
	gp.p += dp1 * sf.N1 + dp2 * sf.N2 + dp3 * sf.N3 + dp4 * sf.N4;

	// strain increment
	de11 = dux_s1 * sf.dN1_dx + dux_s2 * sf.dN2_dx
		 + dux_s3 * sf.dN3_dx + dux_s4 * sf.dN4_dx;
	de22 = duy_s1 * sf.dN1_dy + duy_s2 * sf.dN2_dy
		 + duy_s3 * sf.dN3_dy + duy_s4 * sf.dN4_dy;
	de12 = (dux_s1 * sf.dN1_dy + dux_s2 * sf.dN2_dy
		  + dux_s3 * sf.dN3_dy + dux_s4 * sf.dN4_dy
		  + duy_s1 * sf.dN1_dx + duy_s2 * sf.dN2_dx
		  + duy_s3 * sf.dN3_dx + duy_s4 * sf.dN4_dx) * 0.5;
	gp.e11 += de11;
	gp.e22 += de22;
	gp.e12 += de12;

	// update stress
	double E_tmp = gp.E / (1.0 + gp.niu) / (1.0 - 2.0 * gp.niu);
	ds11 = E_tmp * ((1.0 - gp.niu) * de11 + gp.niu * de22);
	ds22 = E_tmp * (gp.niu * de11 + (1.0 - gp.niu) * de22);
	ds12 = 2.0 * gp.E / (2.0 * (1.0 + gp.niu)) * de12;
	gp.s11 += ds11;
	gp.s22 += ds22;
	gp.s12 += ds12;
	
	// volumetric strain of solid phase
	de_vol_s = de11 + de22;
	// "volumetric strain" of fluid phase
	de_vol_f = -(1.0 - gp.n) / gp.n * de_vol_s
		- (dux_f1 * sf.dN1_dx + dux_f2 * sf.dN2_dx + dux_f3 * sf.dN3_dx + dux_f4 * sf.dN4_dx)
		- (duy_f1 * sf.dN1_dy + duy_f2 * sf.dN2_dy + duy_f3 * sf.dN3_dy + duy_f4 * sf.dN4_dy);
	
	// porosity
	//gp.n = (de_vol_s + gp.n) / (1.0 + de_vol_s);
	// fluid density
	//gp.density_f += gp.density_f * de_vol_f;
}
