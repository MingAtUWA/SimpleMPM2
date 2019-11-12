#include "SimulationCore_pcp.h"

#define KEEP_NEWMARK_BETA_COEFFICIENT
#include "Step_S2D_ME_s_up.h"

namespace
{
	typedef Model_S2D_ME_s_up::Node Node_mpm;
	typedef Model_S2D_ME_s_up::Element Element_mpm;
	typedef Model_S2D_ME_s_up::GaussPoint GaussPoint_mpm;
	typedef Model_S2D_ME_s_up::Particle Particle_mpm;
	typedef Model_S2D_ME_s_up::ShapeFuncValue ShapeFuncValue;

	void print_mat(double mat[12][12],
		std::fstream &out_file,
		const char *mat_name = nullptr)
	{
		if (mat_name)
			out_file << mat_name << "\n";
		for (size_t i = 0; i < 12; i++)
		{
			for (size_t j = 0; j < 12; j++)
				out_file << mat[i][j] << ", ";
			out_file << "\n";
		}
	}

	void add_stiffness_mat(double k_mat[12][12], 
		double E[3][3], double dN_dx[3][8], double weight)
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

	void form_E_matrix(double E_mat[3][3], double E, double niu)
	{
		double G = E / (2.0 * (1.0 + niu));
		E_mat[0][0] = 4.0 / 3.0 * G;
		E_mat[0][1] = -2.0 / 3.0 * G;
		E_mat[0][2] = 0.0;
		E_mat[1][0] = -2.0 / 3.0 * G;
		E_mat[1][1] = 4.0 / 3.0 * G;
		E_mat[1][2] = 0.0;
		E_mat[2][0] = 0.0;
		E_mat[2][1] = 0.0;
		E_mat[2][2] = G;
	}
};

void Step_S2D_ME_s_up::reequilibration_elem_stiffness_mat_and_force_vec(
	Model_S2D_ME_s_up::Element &e, double kmat[12][12], double fvec[12])
{
	double mat_term;
	double E_mat[3][3], dN_dx_mat[3][8];
	
	memset(kmat, 0, sizeof(double) * 12 * 12);
	memset(fvec, 0, sizeof(double) * 12);
	for (Model_S2D_ME_s_up::Particle *pcl_iter = e.pcls;
		 pcl_iter; pcl_iter = pcl_iter->next)
	{
		Model_S2D_ME_s_up::Particle &pcl = *pcl_iter;

		// stiffness matrix
		// B * D * B
		form_E_matrix(E_mat, model->E, model->niu);
		dN_dx_mat[0][0] = pcl.dN1_dx;
		dN_dx_mat[0][1] = pcl.dN2_dx;
		dN_dx_mat[0][2] = pcl.dN3_dx;
		dN_dx_mat[0][3] = pcl.dN4_dx;
		dN_dx_mat[0][4] = 0.0;
		dN_dx_mat[0][5] = 0.0;
		dN_dx_mat[0][6] = 0.0;
		dN_dx_mat[0][7] = 0.0;
		dN_dx_mat[1][0] = 0.0;
		dN_dx_mat[1][1] = 0.0;
		dN_dx_mat[1][2] = 0.0;
		dN_dx_mat[1][3] = 0.0;
		dN_dx_mat[1][4] = pcl.dN1_dy;
		dN_dx_mat[1][5] = pcl.dN2_dy;
		dN_dx_mat[1][6] = pcl.dN3_dy;
		dN_dx_mat[1][7] = pcl.dN4_dy;
		dN_dx_mat[2][0] = pcl.dN1_dy;
		dN_dx_mat[2][1] = pcl.dN2_dy;
		dN_dx_mat[2][2] = pcl.dN3_dy;
		dN_dx_mat[2][3] = pcl.dN4_dy;
		dN_dx_mat[2][4] = pcl.dN1_dx;
		dN_dx_mat[2][5] = pcl.dN2_dx;
		dN_dx_mat[2][6] = pcl.dN3_dx;
		dN_dx_mat[2][7] = pcl.dN4_dx;
		add_stiffness_mat(kmat, E_mat, dN_dx_mat, pcl.vol);

		// dNi_dx * Nj
		mat_term = pcl.dN1_dx * pcl.N1 * pcl.vol;
		kmat[0][8] += mat_term;
		kmat[8][0] += mat_term;
		mat_term = pcl.dN1_dx * pcl.N2 * pcl.vol;
		kmat[0][9] += mat_term;
		kmat[9][0] += mat_term;
		mat_term = pcl.dN1_dx * pcl.N3 * pcl.vol;
		kmat[0][10] += mat_term;
		kmat[10][0] += mat_term;
		mat_term = pcl.dN1_dx * pcl.N4 * pcl.vol;
		kmat[0][11] += mat_term;
		kmat[11][0] += mat_term;
		mat_term = pcl.dN2_dx * pcl.N1 * pcl.vol;
		kmat[1][8] += mat_term;
		kmat[8][1] += mat_term;
		mat_term = pcl.dN2_dx * pcl.N2 * pcl.vol;
		kmat[1][9] += mat_term;
		kmat[9][1] += mat_term;
		mat_term = pcl.dN2_dx * pcl.N3 * pcl.vol;
		kmat[1][10] += mat_term;
		kmat[10][1] += mat_term;
		mat_term = pcl.dN2_dx * pcl.N4 * pcl.vol;
		kmat[1][11] += mat_term;
		kmat[11][1] += mat_term;
		mat_term = pcl.dN3_dx * pcl.N1 * pcl.vol;
		kmat[2][8] += mat_term;
		kmat[8][2] += mat_term;
		mat_term = pcl.dN3_dx * pcl.N2 * pcl.vol;
		kmat[2][9] += mat_term;
		kmat[9][2] += mat_term;
		mat_term = pcl.dN3_dx * pcl.N3 * pcl.vol;
		kmat[2][10] += mat_term;
		kmat[10][2] += mat_term;
		mat_term = pcl.dN3_dx * pcl.N4 * pcl.vol;
		kmat[2][11] += mat_term;
		kmat[11][2] += mat_term;
		mat_term = pcl.dN4_dx * pcl.N1 * pcl.vol;
		kmat[3][8] += mat_term;
		kmat[8][3] += mat_term;
		mat_term = pcl.dN4_dx * pcl.N2 * pcl.vol;
		kmat[3][9] += mat_term;
		kmat[9][3] += mat_term;
		mat_term = pcl.dN4_dx * pcl.N3 * pcl.vol;
		kmat[3][10] += mat_term;
		kmat[10][3] += mat_term;
		mat_term = pcl.dN4_dx * pcl.N4 * pcl.vol;
		kmat[3][11] += mat_term;
		kmat[11][3] += mat_term;

		// dNi_dy * Nj
		mat_term = pcl.dN1_dy * pcl.N1 * pcl.vol;
		kmat[4][8] += mat_term;
		kmat[8][4] += mat_term;
		mat_term = pcl.dN1_dy * pcl.N2 * pcl.vol;
		kmat[4][9] += mat_term;
		kmat[9][4] += mat_term;
		mat_term = pcl.dN1_dy * pcl.N3 * pcl.vol;
		kmat[4][10] += mat_term;
		kmat[10][4] += mat_term;
		mat_term = pcl.dN1_dy * pcl.N4 * pcl.vol;
		kmat[4][11] += mat_term;
		kmat[11][4] += mat_term;
		mat_term = pcl.dN2_dy * pcl.N1 * pcl.vol;
		kmat[5][8] += mat_term;
		kmat[8][5] += mat_term;
		mat_term = pcl.dN2_dy * pcl.N2 * pcl.vol;
		kmat[5][9] += mat_term;
		kmat[9][5] += mat_term;
		mat_term = pcl.dN2_dy * pcl.N3 * pcl.vol;
		kmat[5][10] += mat_term;
		kmat[10][5] += mat_term;
		mat_term = pcl.dN2_dy * pcl.N4 * pcl.vol;
		kmat[5][11] += mat_term;
		kmat[11][5] += mat_term;
		mat_term = pcl.dN3_dy * pcl.N1 * pcl.vol;
		kmat[6][8] += mat_term;
		kmat[8][6] += mat_term;
		mat_term = pcl.dN3_dy * pcl.N2 * pcl.vol;
		kmat[6][9] += mat_term;
		kmat[9][6] += mat_term;
		mat_term = pcl.dN3_dy * pcl.N3 * pcl.vol;
		kmat[6][10] += mat_term;
		kmat[10][6] += mat_term;
		mat_term = pcl.dN3_dy * pcl.N4 * pcl.vol;
		kmat[6][11] += mat_term;
		kmat[11][6] += mat_term;
		mat_term = pcl.dN4_dy * pcl.N1 * pcl.vol;
		kmat[7][8] += mat_term;
		kmat[8][7] += mat_term;
		mat_term = pcl.dN4_dy * pcl.N2 * pcl.vol;
		kmat[7][9] += mat_term;
		kmat[9][7] += mat_term;
		mat_term = pcl.dN4_dy * pcl.N3 * pcl.vol;
		kmat[7][10] += mat_term;
		kmat[10][7] += mat_term;
		mat_term = pcl.dN4_dy * pcl.N4 * pcl.vol;
		kmat[7][11] += mat_term;
		kmat[11][7] += mat_term;

		// Ni * Nj / K
		double K = model->K;
		mat_term = -pcl.N1 * pcl.N1 / K * pcl.vol;
		kmat[8][8] += mat_term;
		mat_term = -pcl.N1 * pcl.N2 / K * pcl.vol;
		kmat[8][9] += mat_term;
		kmat[9][8] += mat_term;
		mat_term = -pcl.N1 * pcl.N3 / K * pcl.vol;
		kmat[8][10] += mat_term;
		kmat[10][8] += mat_term;
		mat_term = -pcl.N1 * pcl.N4 / K * pcl.vol;
		kmat[8][11] += mat_term;
		kmat[11][8] += mat_term;
		mat_term = -pcl.N2 * pcl.N2 / K * pcl.vol;
		kmat[9][9] += mat_term;
		mat_term = -pcl.N2 * pcl.N3 / K * pcl.vol;
		kmat[9][10] += mat_term;
		kmat[10][9] += mat_term;
		mat_term = -pcl.N2 * pcl.N4 / K * pcl.vol;
		kmat[9][11] += mat_term;
		kmat[11][9] += mat_term;
		mat_term = -pcl.N3 * pcl.N3 / K * pcl.vol;
		kmat[10][10] += mat_term;
		mat_term = -pcl.N3 * pcl.N4 / K * pcl.vol;
		kmat[10][11] += mat_term;
		kmat[11][10] += mat_term;
		mat_term = -pcl.N4 * pcl.N4 / K * pcl.vol;
		kmat[11][11] += mat_term;

		// force vector
		// internal force
		fvec[0] -= (pcl.dN1_dx * (pcl.s11 + pcl.p) + pcl.dN1_dy * pcl.s12) * pcl.vol + pcl.N1 * pcl.m * pcl.ax;
		fvec[1] -= (pcl.dN2_dx * (pcl.s11 + pcl.p) + pcl.dN2_dy * pcl.s12) * pcl.vol + pcl.N2 * pcl.m * pcl.ax;
		fvec[2] -= (pcl.dN3_dx * (pcl.s11 + pcl.p) + pcl.dN3_dy * pcl.s12) * pcl.vol + pcl.N3 * pcl.m * pcl.ax;
		fvec[3] -= (pcl.dN4_dx * (pcl.s11 + pcl.p) + pcl.dN4_dy * pcl.s12) * pcl.vol + pcl.N4 * pcl.m * pcl.ax;
		fvec[4] -= (pcl.dN1_dx * pcl.s12 + pcl.dN1_dy * (pcl.s22 + pcl.p)) * pcl.vol + pcl.N1 * pcl.m * pcl.ay;
		fvec[5] -= (pcl.dN2_dx * pcl.s12 + pcl.dN2_dy * (pcl.s22 + pcl.p)) * pcl.vol + pcl.N2 * pcl.m * pcl.ay;
		fvec[6] -= (pcl.dN3_dx * pcl.s12 + pcl.dN3_dy * (pcl.s22 + pcl.p)) * pcl.vol + pcl.N3 * pcl.m * pcl.ay;
		fvec[7] -= (pcl.dN4_dx * pcl.s12 + pcl.dN4_dy * (pcl.s22 + pcl.p)) * pcl.vol + pcl.N4 * pcl.m * pcl.ay;
	}
}

void Step_S2D_ME_s_up::form_elem_stiffness_mat_and_force_vec(
	Model_S2D_ME_s_up::Element &e, double kmat[12][12], double fvec[12])
{
	Model_S2D_ME_s_up &md = *static_cast<Model_S2D_ME_s_up *>(model);
	double mat_term;
	double E_mat[3][3], dN_dx_mat[3][8];
	double coef_m = 1.0 / (beta * dtime * dtime);
	double coef_a = 1.0 / (2.0 * beta) - 1.0;
	double coef_v = 1.0 / (beta * dtime);

	memset(kmat, 0, sizeof(double) * 12 * 12);
	memset(fvec, 0, sizeof(double) * 12);
	for (Model_S2D_ME_s_up::Particle *pcl_iter = e.pcls;
		 pcl_iter; pcl_iter = pcl_iter->next)
	{
		Model_S2D_ME_s_up::Particle &pcl = *pcl_iter;

		// mass matrix
		mat_term = pcl.N1 * pcl.density * pcl.N1 * pcl.vol * coef_m;
		kmat[0][0] += mat_term;
		kmat[4][4] += mat_term;
		mat_term = pcl.N1 * pcl.density * pcl.N2 * pcl.vol * coef_m;
		kmat[0][1] += mat_term;
		kmat[1][0] += mat_term;
		kmat[4][5] += mat_term;
		kmat[5][4] += mat_term;
		mat_term = pcl.N1 * pcl.density * pcl.N3 * pcl.vol * coef_m;
		kmat[0][2] += mat_term;
		kmat[2][0] += mat_term;
		kmat[4][6] += mat_term;
		kmat[6][4] += mat_term;
		mat_term = pcl.N1 * pcl.density * pcl.N4 * pcl.vol * coef_m;
		kmat[0][3] += mat_term;
		kmat[3][0] += mat_term;
		kmat[4][7] += mat_term;
		kmat[7][4] += mat_term;
		mat_term = pcl.N2 * pcl.density * pcl.N2 * pcl.vol * coef_m;
		kmat[1][1] += mat_term;
		kmat[5][5] += mat_term;
		mat_term = pcl.N2 * pcl.density * pcl.N3 * pcl.vol * coef_m;
		kmat[1][2] += mat_term;
		kmat[2][1] += mat_term;
		kmat[5][6] += mat_term;
		kmat[6][5] += mat_term;
		mat_term = pcl.N2 * pcl.density * pcl.N4 * pcl.vol * coef_m;
		kmat[1][3] += mat_term;
		kmat[3][1] += mat_term;
		kmat[5][7] += mat_term;
		kmat[7][5] += mat_term;
		mat_term = pcl.N3 * pcl.density * pcl.N3 * pcl.vol * coef_m;
		kmat[2][2] += mat_term;
		kmat[6][6] += mat_term;
		mat_term = pcl.N3 * pcl.density * pcl.N4 * pcl.vol * coef_m;
		kmat[2][3] += mat_term;
		kmat[3][2] += mat_term;
		kmat[6][7] += mat_term;
		kmat[7][6] += mat_term;
		mat_term = pcl.N4 * pcl.density * pcl.N4 * pcl.vol * coef_m;
		kmat[3][3] += mat_term;
		kmat[7][7] += mat_term;

		// stiffness matrix
		// B * D * B
		form_E_matrix(E_mat, model->E, model->niu);
		dN_dx_mat[0][0] = pcl.dN1_dx;
		dN_dx_mat[0][1] = pcl.dN2_dx;
		dN_dx_mat[0][2] = pcl.dN3_dx;
		dN_dx_mat[0][3] = pcl.dN4_dx;
		dN_dx_mat[0][4] = 0.0;
		dN_dx_mat[0][5] = 0.0;
		dN_dx_mat[0][6] = 0.0;
		dN_dx_mat[0][7] = 0.0;
		dN_dx_mat[1][0] = 0.0;
		dN_dx_mat[1][1] = 0.0;
		dN_dx_mat[1][2] = 0.0;
		dN_dx_mat[1][3] = 0.0;
		dN_dx_mat[1][4] = pcl.dN1_dy;
		dN_dx_mat[1][5] = pcl.dN2_dy;
		dN_dx_mat[1][6] = pcl.dN3_dy;
		dN_dx_mat[1][7] = pcl.dN4_dy;
		dN_dx_mat[2][0] = pcl.dN1_dy;
		dN_dx_mat[2][1] = pcl.dN2_dy;
		dN_dx_mat[2][2] = pcl.dN3_dy;
		dN_dx_mat[2][3] = pcl.dN4_dy;
		dN_dx_mat[2][4] = pcl.dN1_dx;
		dN_dx_mat[2][5] = pcl.dN2_dx;
		dN_dx_mat[2][6] = pcl.dN3_dx;
		dN_dx_mat[2][7] = pcl.dN4_dx;
		add_stiffness_mat(kmat, E_mat, dN_dx_mat, pcl.vol);

		// dNi_dx * Nj
		mat_term = pcl.dN1_dx * pcl.N1 * pcl.vol;
		kmat[0][8] += mat_term;
		kmat[8][0] += mat_term;
		mat_term = pcl.dN1_dx * pcl.N2 * pcl.vol;
		kmat[0][9] += mat_term;
		kmat[9][0] += mat_term;
		mat_term = pcl.dN1_dx * pcl.N3 * pcl.vol;
		kmat[0][10] += mat_term;
		kmat[10][0] += mat_term;
		mat_term = pcl.dN1_dx * pcl.N4 * pcl.vol;
		kmat[0][11] += mat_term;
		kmat[11][0] += mat_term;
		mat_term = pcl.dN2_dx * pcl.N1 * pcl.vol;
		kmat[1][8] += mat_term;
		kmat[8][1] += mat_term;
		mat_term = pcl.dN2_dx * pcl.N2 * pcl.vol;
		kmat[1][9] += mat_term;
		kmat[9][1] += mat_term;
		mat_term = pcl.dN2_dx * pcl.N3 * pcl.vol;
		kmat[1][10] += mat_term;
		kmat[10][1] += mat_term;
		mat_term = pcl.dN2_dx * pcl.N4 * pcl.vol;
		kmat[1][11] += mat_term;
		kmat[11][1] += mat_term;
		mat_term = pcl.dN3_dx * pcl.N1 * pcl.vol;
		kmat[2][8] += mat_term;
		kmat[8][2] += mat_term;
		mat_term = pcl.dN3_dx * pcl.N2 * pcl.vol;
		kmat[2][9] += mat_term;
		kmat[9][2] += mat_term;
		mat_term = pcl.dN3_dx * pcl.N3 * pcl.vol;
		kmat[2][10] += mat_term;
		kmat[10][2] += mat_term;
		mat_term = pcl.dN3_dx * pcl.N4 * pcl.vol;
		kmat[2][11] += mat_term;
		kmat[11][2] += mat_term;
		mat_term = pcl.dN4_dx * pcl.N1 * pcl.vol;
		kmat[3][8] += mat_term;
		kmat[8][3] += mat_term;
		mat_term = pcl.dN4_dx * pcl.N2 * pcl.vol;
		kmat[3][9] += mat_term;
		kmat[9][3] += mat_term;
		mat_term = pcl.dN4_dx * pcl.N3 * pcl.vol;
		kmat[3][10] += mat_term;
		kmat[10][3] += mat_term;
		mat_term = pcl.dN4_dx * pcl.N4 * pcl.vol;
		kmat[3][11] += mat_term;
		kmat[11][3] += mat_term;

		// dNi_dy * Nj
		mat_term = pcl.dN1_dy * pcl.N1 * pcl.vol;
		kmat[4][8] += mat_term;
		kmat[8][4] += mat_term;
		mat_term = pcl.dN1_dy * pcl.N2 * pcl.vol;
		kmat[4][9] += mat_term;
		kmat[9][4] += mat_term;
		mat_term = pcl.dN1_dy * pcl.N3 * pcl.vol;
		kmat[4][10] += mat_term;
		kmat[10][4] += mat_term;
		mat_term = pcl.dN1_dy * pcl.N4 * pcl.vol;
		kmat[4][11] += mat_term;
		kmat[11][4] += mat_term;
		mat_term = pcl.dN2_dy * pcl.N1 * pcl.vol;
		kmat[5][8] += mat_term;
		kmat[8][5] += mat_term;
		mat_term = pcl.dN2_dy * pcl.N2 * pcl.vol;
		kmat[5][9] += mat_term;
		kmat[9][5] += mat_term;
		mat_term = pcl.dN2_dy * pcl.N3 * pcl.vol;
		kmat[5][10] += mat_term;
		kmat[10][5] += mat_term;
		mat_term = pcl.dN2_dy * pcl.N4 * pcl.vol;
		kmat[5][11] += mat_term;
		kmat[11][5] += mat_term;
		mat_term = pcl.dN3_dy * pcl.N1 * pcl.vol;
		kmat[6][8] += mat_term;
		kmat[8][6] += mat_term;
		mat_term = pcl.dN3_dy * pcl.N2 * pcl.vol;
		kmat[6][9] += mat_term;
		kmat[9][6] += mat_term;
		mat_term = pcl.dN3_dy * pcl.N3 * pcl.vol;
		kmat[6][10] += mat_term;
		kmat[10][6] += mat_term;
		mat_term = pcl.dN3_dy * pcl.N4 * pcl.vol;
		kmat[6][11] += mat_term;
		kmat[11][6] += mat_term;
		mat_term = pcl.dN4_dy * pcl.N1 * pcl.vol;
		kmat[7][8] += mat_term;
		kmat[8][7] += mat_term;
		mat_term = pcl.dN4_dy * pcl.N2 * pcl.vol;
		kmat[7][9] += mat_term;
		kmat[9][7] += mat_term;
		mat_term = pcl.dN4_dy * pcl.N3 * pcl.vol;
		kmat[7][10] += mat_term;
		kmat[10][7] += mat_term;
		mat_term = pcl.dN4_dy * pcl.N4 * pcl.vol;
		kmat[7][11] += mat_term;
		kmat[11][7] += mat_term;

		// Ni * Nj / K
		double K = model->K;
		mat_term = -pcl.N1 * pcl.N1 / K * pcl.vol;
		kmat[8][8] += mat_term;
		mat_term = -pcl.N1 * pcl.N2 / K * pcl.vol;
		kmat[8][9] += mat_term;
		kmat[9][8] += mat_term;
		mat_term = -pcl.N1 * pcl.N3 / K * pcl.vol;
		kmat[8][10] += mat_term;
		kmat[10][8] += mat_term;
		mat_term = -pcl.N1 * pcl.N4 / K * pcl.vol;
		kmat[8][11] += mat_term;
		kmat[11][8] += mat_term;
		mat_term = -pcl.N2 * pcl.N2 / K * pcl.vol;
		kmat[9][9] += mat_term;
		mat_term = -pcl.N2 * pcl.N3 / K * pcl.vol;
		kmat[9][10] += mat_term;
		kmat[10][9] += mat_term;
		mat_term = -pcl.N2 * pcl.N4 / K * pcl.vol;
		kmat[9][11] += mat_term;
		kmat[11][9] += mat_term;
		mat_term = -pcl.N3 * pcl.N3 / K * pcl.vol;
		kmat[10][10] += mat_term;
		mat_term = -pcl.N3 * pcl.N4 / K * pcl.vol;
		kmat[10][11] += mat_term;
		kmat[11][10] += mat_term;
		mat_term = -pcl.N4 * pcl.N4 / K * pcl.vol;
		kmat[11][11] += mat_term;

		// force vector
		// internal force
		fvec[0] += -(pcl.dN1_dx * (pcl.s11 + pcl.p) + pcl.dN1_dy * pcl.s12) * pcl.vol;
		fvec[1] += -(pcl.dN2_dx * (pcl.s11 + pcl.p) + pcl.dN2_dy * pcl.s12) * pcl.vol;
		fvec[2] += -(pcl.dN3_dx * (pcl.s11 + pcl.p) + pcl.dN3_dy * pcl.s12) * pcl.vol;
		fvec[3] += -(pcl.dN4_dx * (pcl.s11 + pcl.p) + pcl.dN4_dy * pcl.s12) * pcl.vol;
		fvec[4] += -(pcl.dN1_dx * pcl.s12 + pcl.dN1_dy * (pcl.s22 + pcl.p)) * pcl.vol;
		fvec[5] += -(pcl.dN2_dx * pcl.s12 + pcl.dN2_dy * (pcl.s22 + pcl.p)) * pcl.vol;
		fvec[6] += -(pcl.dN3_dx * pcl.s12 + pcl.dN3_dy * (pcl.s22 + pcl.p)) * pcl.vol;
		fvec[7] += -(pcl.dN4_dx * pcl.s12 + pcl.dN4_dy * (pcl.s22 + pcl.p)) * pcl.vol;
		
		// a vector
		Element_mpm &e = *pcl.pe;
		Node_mpm &n1 = md.nodes[e.n1_id];
		Node_mpm &n2 = md.nodes[e.n2_id];
		Node_mpm &n3 = md.nodes[e.n3_id];
		Node_mpm &n4 = md.nodes[e.n4_id];
		double pcl_ax = n1.ax * pcl.N1 + n2.ax * pcl.N2 + n3.ax * pcl.N3 + n4.ax * pcl.N4;
		fvec[0] += pcl.N1 * pcl.m * pcl_ax * coef_a;
		fvec[1] += pcl.N2 * pcl.m * pcl_ax * coef_a;
		fvec[2] += pcl.N3 * pcl.m * pcl_ax * coef_a;
		fvec[3] += pcl.N4 * pcl.m * pcl_ax * coef_a;
		double pcl_ay = n1.ay * pcl.N1 + n2.ay * pcl.N2 + n3.ay * pcl.N3 + n4.ay * pcl.N4;
		fvec[4] += pcl.N1 * pcl.m * pcl_ay * coef_a;
		fvec[5] += pcl.N2 * pcl.m * pcl_ay * coef_a;
		fvec[6] += pcl.N3 * pcl.m * pcl_ay * coef_a;
		fvec[7] += pcl.N4 * pcl.m * pcl_ay * coef_a;

		// v vector
		double pcl_vx = n1.vx * pcl.N1 + n2.vx * pcl.N2 + n3.vx * pcl.N3 + n4.vx * pcl.N4;
		fvec[0] += pcl.N1 * pcl.m * pcl_vx * coef_v;
		fvec[1] += pcl.N2 * pcl.m * pcl_vx * coef_v;
		fvec[2] += pcl.N3 * pcl.m * pcl_vx * coef_v;
		fvec[3] += pcl.N4 * pcl.m * pcl_vx * coef_v;
		double pcl_vy = n1.vy * pcl.N1 + n2.vy * pcl.N2 + n3.vy * pcl.N3 + n4.vy * pcl.N4;
		fvec[4] += pcl.N1 * pcl.m * pcl_vy * coef_v;
		fvec[5] += pcl.N2 * pcl.m * pcl_vy * coef_v;
		fvec[6] += pcl.N3 * pcl.m * pcl_vy * coef_v;
		fvec[7] += pcl.N4 * pcl.m * pcl_vy * coef_v;
	}

	//print_mat(kmat, out_file);
}
