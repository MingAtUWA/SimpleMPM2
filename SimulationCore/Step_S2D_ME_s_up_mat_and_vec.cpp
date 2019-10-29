#include "SimulationCore_pcp.h"

#define KEEP_NEWMARK_BETA_COEFFICIENT
#include "Step_S2D_ME_s_up.h"

namespace
{
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

	void print_vec(double *vec, size_t num)
	{
		std::cout << "vec\n";
		for (size_t i = 0; i < num; i++)
			printf("%+8.2e ", vec[i]);
		std::cout << "\n";
	}

	void cal_stiffness_mat(double k_mat[12][12], 
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


void Step_S2D_ME_s_up::form_elem_stiffness_mat_and_force_vec_pure(
	Model_S2D_ME_s_up::Element &e, double kmat[12][12], double fvec[12])
{
	double mat_term;
	double coef_m = 1.0 / (beta * dtime * dtime);
	double coef_a = 1.0 / (2.0 * beta) - 1.0;
	double coef_v = 1.0 / (beta * dtime);
	memset(kmat, 0, sizeof(double) * 12 * 12);
	memset(fvec, 0, sizeof(double) * 12);
	double E_mat[3][3], dN_dx_mat[3][8];

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
		cal_stiffness_mat(kmat, E_mat, dN_dx_mat, pcl.vol);

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
		fvec[0] += pcl.N1 * pcl.m * pcl.ax * coef_a;
		fvec[1] += pcl.N2 * pcl.m * pcl.ax * coef_a;
		fvec[2] += pcl.N3 * pcl.m * pcl.ax * coef_a;
		fvec[3] += pcl.N4 * pcl.m * pcl.ax * coef_a;
		fvec[4] += pcl.N1 * pcl.m * pcl.ay * coef_a;
		fvec[5] += pcl.N2 * pcl.m * pcl.ay * coef_a;
		fvec[6] += pcl.N3 * pcl.m * pcl.ay * coef_a;
		fvec[7] += pcl.N4 * pcl.m * pcl.ay * coef_a;

		// v vector
		fvec[0] += pcl.N1 * pcl.m * pcl.vx * coef_v;
		fvec[1] += pcl.N2 * pcl.m * pcl.vx * coef_v;
		fvec[2] += pcl.N3 * pcl.m * pcl.vx * coef_v;
		fvec[3] += pcl.N4 * pcl.m * pcl.vx * coef_v;
		fvec[4] += pcl.N1 * pcl.m * pcl.vy * coef_v;
		fvec[5] += pcl.N2 * pcl.m * pcl.vy * coef_v;
		fvec[6] += pcl.N3 * pcl.m * pcl.vy * coef_v;
		fvec[7] += pcl.N4 * pcl.m * pcl.vy * coef_v;
	}

	//print_mat(kmat, out_file);
}


void Step_S2D_ME_s_up::form_elem_stiffness_mat_and_force_vec(
	Model_S2D_ME_s_up::Element &e, double kmat[12][12], double fvec[12])
{
	ShapeFuncValue &sf1 = model->sf1;
	ShapeFuncValue &sf2 = model->sf2;
	ShapeFuncValue &sf3 = model->sf3;
	ShapeFuncValue &sf4 = model->sf4;
	GaussPoint_mpm &gp1 = e.gp1;
	GaussPoint_mpm &gp2 = e.gp2;
	GaussPoint_mpm &gp3 = e.gp3;
	GaussPoint_mpm &gp4 = e.gp4;
	double gp_w = e.vf * model->gp_w;
	double coef, mat_term;

	memset(kmat, 0, sizeof(double) * 12 * 12);
	// mass matrix
	coef = 1.0 / (beta * dtime * dtime);
	mat_term = (sf1.N1 * gp1.density * sf1.N1 + sf2.N1 * gp2.density * sf2.N1
			  + sf3.N1 * gp3.density * sf3.N1 + sf4.N1 * gp4.density * sf4.N1) * gp_w * coef;
	kmat[0][0] = mat_term;
	kmat[4][4] = mat_term;
	mat_term = (sf1.N1 * gp1.density * sf1.N2 + sf2.N1 * gp2.density * sf2.N2
			  + sf3.N1 * gp3.density * sf3.N2 + sf4.N1 * gp4.density * sf4.N2) * gp_w * coef;
	kmat[0][1] = mat_term;
	kmat[1][0] = mat_term;
	kmat[4][5] = mat_term;
	kmat[5][4] = mat_term;
	mat_term = (sf1.N1 * gp1.density * sf1.N3 + sf2.N1 * gp2.density * sf2.N3
			  + sf3.N1 * gp3.density * sf3.N3 + sf4.N1 * gp4.density * sf4.N3) * gp_w * coef;
	kmat[0][2] = mat_term;
	kmat[2][0] = mat_term;
	kmat[4][6] = mat_term;
	kmat[6][4] = mat_term;
	mat_term = (sf1.N1 * gp1.density * sf1.N4 + sf2.N1 * gp2.density * sf2.N4
			  + sf3.N1 * gp3.density * sf3.N4 + sf4.N1 * gp4.density * sf4.N4) * gp_w * coef;
	kmat[0][3] = mat_term;
	kmat[3][0] = mat_term;
	kmat[4][7] = mat_term;
	kmat[7][4] = mat_term;
	mat_term = (sf1.N2 * gp1.density * sf1.N2 + sf2.N2 * gp2.density * sf2.N2
			  + sf3.N2 * gp3.density * sf3.N2 + sf4.N2 * gp4.density * sf4.N2) * gp_w * coef;
	kmat[1][1] = mat_term;
	kmat[5][5] = mat_term;
	mat_term = (sf1.N2 * gp1.density * sf1.N3 + sf2.N2 * gp2.density * sf2.N3
			  + sf3.N2 * gp3.density * sf3.N3 + sf4.N2 * gp4.density * sf4.N3) * gp_w * coef;
	kmat[1][2] = mat_term;
	kmat[2][1] = mat_term;
	kmat[5][6] = mat_term;
	kmat[6][5] = mat_term;
	mat_term = (sf1.N2 * gp1.density * sf1.N4 + sf2.N2 * gp2.density * sf2.N4
			  + sf3.N2 * gp3.density * sf3.N4 + sf4.N2 * gp4.density * sf4.N4) * gp_w * coef;
	kmat[1][3] = mat_term;
	kmat[3][1] = mat_term;
	kmat[5][7] = mat_term;
	kmat[7][5] = mat_term;
	mat_term = (sf1.N3 * gp1.density * sf1.N3 + sf2.N3 * gp2.density * sf2.N3
			  + sf3.N3 * gp3.density * sf3.N3 + sf4.N3 * gp4.density * sf4.N3) * gp_w * coef;
	kmat[2][2] = mat_term;
	kmat[6][6] = mat_term;
	mat_term = (sf1.N3 * gp1.density * sf1.N4 + sf2.N3 * gp2.density * sf2.N4
			  + sf3.N3 * gp3.density * sf3.N4 + sf4.N3 * gp4.density * sf4.N4) * gp_w * coef;
	kmat[2][3] = mat_term;
	kmat[3][2] = mat_term;
	kmat[6][7] = mat_term;
	kmat[7][6] = mat_term;
	mat_term = (sf1.N4 * gp1.density * sf1.N4 + sf2.N4 * gp2.density * sf2.N4
			  + sf3.N4 * gp3.density * sf3.N4 + sf4.N4 * gp4.density * sf4.N4) * gp_w * coef;
	kmat[3][3] = mat_term;
	kmat[7][7] = mat_term;

	// stiffness matrix
	// B * D * B
	double E_mat[3][3];
	form_E_matrix(E_mat, model->E, model->niu);
	// gauss point 1
	cal_stiffness_mat(kmat, E_mat, model->dN_dx_mat1, gp_w);
	// gauss point 2
	cal_stiffness_mat(kmat, E_mat, model->dN_dx_mat2, gp_w);
	// gauss point 3
	cal_stiffness_mat(kmat, E_mat, model->dN_dx_mat3, gp_w);
	// gauss point 4
	cal_stiffness_mat(kmat, E_mat, model->dN_dx_mat4, gp_w);

	// dNi_dx * Nj
	mat_term = (sf1.dN1_dx * sf1.N1 + sf2.dN1_dx * sf2.N1
			  + sf3.dN1_dx * sf3.N1 + sf4.dN1_dx * sf4.N1) * gp_w;
	kmat[0][8] += mat_term;
	kmat[8][0] += mat_term;
	mat_term = (sf1.dN1_dx * sf1.N2 + sf2.dN1_dx * sf2.N2
			  + sf3.dN1_dx * sf3.N2 + sf4.dN1_dx * sf4.N2) * gp_w;
	kmat[0][9] += mat_term;
	kmat[9][0] += mat_term;
	mat_term = (sf1.dN1_dx * sf1.N3 + sf2.dN1_dx * sf2.N3
			  + sf3.dN1_dx * sf3.N3 + sf4.dN1_dx * sf4.N3) * gp_w;
	kmat[0][10] += mat_term;
	kmat[10][0] += mat_term;
	mat_term = (sf1.dN1_dx * sf1.N4 + sf2.dN1_dx * sf2.N4
			  + sf3.dN1_dx * sf3.N4 + sf4.dN1_dx * sf4.N4) * gp_w;
	kmat[0][11] += mat_term;
	kmat[11][0] += mat_term;
	mat_term = (sf1.dN2_dx * sf1.N1 + sf2.dN2_dx * sf2.N1
			  + sf3.dN2_dx * sf3.N1 + sf4.dN2_dx * sf4.N1) * gp_w;
	kmat[1][8] += mat_term;
	kmat[8][1] += mat_term;
	mat_term = (sf1.dN2_dx * sf1.N2 + sf2.dN2_dx * sf2.N2
			  + sf3.dN2_dx * sf3.N2 + sf4.dN2_dx * sf4.N2) * gp_w;
	kmat[1][9] += mat_term;
	kmat[9][1] += mat_term;
	mat_term = (sf1.dN2_dx * sf1.N3 + sf2.dN2_dx * sf2.N3
			  + sf3.dN2_dx * sf3.N3 + sf4.dN2_dx * sf4.N3) * gp_w;
	kmat[1][10] += mat_term;
	kmat[10][1] += mat_term;
	mat_term = (sf1.dN2_dx * sf1.N4 + sf2.dN2_dx * sf2.N4
			  + sf3.dN2_dx * sf3.N4 + sf4.dN2_dx * sf4.N4) * gp_w;
	kmat[1][11] += mat_term;
	kmat[11][1] += mat_term;
	mat_term = (sf1.dN3_dx * sf1.N1 + sf2.dN3_dx * sf2.N1
			  + sf3.dN3_dx * sf3.N1 + sf4.dN3_dx * sf4.N1) * gp_w;
	kmat[2][8] += mat_term;
	kmat[8][2] += mat_term;
	mat_term = (sf1.dN3_dx * sf1.N2 + sf2.dN3_dx * sf2.N2
			  + sf3.dN3_dx * sf3.N2 + sf4.dN3_dx * sf4.N2) * gp_w;
	kmat[2][9] += mat_term;
	kmat[9][2] += mat_term;
	mat_term = (sf1.dN3_dx * sf1.N3 + sf2.dN3_dx * sf2.N3
			  + sf3.dN3_dx * sf3.N3 + sf4.dN3_dx * sf4.N3) * gp_w;
	kmat[2][10] += mat_term;
	kmat[10][2] += mat_term;
	mat_term = (sf1.dN3_dx * sf1.N4 + sf2.dN3_dx * sf2.N4
			  + sf3.dN3_dx * sf3.N4 + sf4.dN3_dx * sf4.N4) * gp_w;
	kmat[2][11] += mat_term;
	kmat[11][2] += mat_term;
	mat_term = (sf1.dN4_dx * sf1.N1 + sf2.dN4_dx * sf2.N1
			  + sf3.dN4_dx * sf3.N1 + sf4.dN4_dx * sf4.N1) * gp_w;
	kmat[3][8] += mat_term;
	kmat[8][3] += mat_term;
	mat_term = (sf1.dN4_dx * sf1.N2 + sf2.dN4_dx * sf2.N2
			  + sf3.dN4_dx * sf3.N2 + sf4.dN4_dx * sf4.N2) * gp_w;
	kmat[3][9] += mat_term;
	kmat[9][3] += mat_term;
	mat_term = (sf1.dN4_dx * sf1.N3 + sf2.dN4_dx * sf2.N3
			  + sf3.dN4_dx * sf3.N3 + sf4.dN4_dx * sf4.N3) * gp_w;
	kmat[3][10] += mat_term;
	kmat[10][3] += mat_term;
	mat_term = (sf1.dN4_dx * sf1.N4 + sf2.dN4_dx * sf2.N4
			  + sf3.dN4_dx * sf3.N4 + sf4.dN4_dx * sf4.N4) * gp_w;
	kmat[3][11] += mat_term;
	kmat[11][3] += mat_term;
	// dNi_dy * Nj
	mat_term = (sf1.dN1_dy * sf1.N1 + sf2.dN1_dy * sf2.N1
			  + sf3.dN1_dy * sf3.N1 + sf4.dN1_dy * sf4.N1) * gp_w;
	kmat[4][8] += mat_term;
	kmat[8][4] += mat_term;
	mat_term = (sf1.dN1_dy * sf1.N2 + sf2.dN1_dy * sf2.N2
			  + sf3.dN1_dy * sf3.N2 + sf4.dN1_dy * sf4.N2) * gp_w;
	kmat[4][9] += mat_term;
	kmat[9][4] += mat_term;
	mat_term = (sf1.dN1_dy * sf1.N3 + sf2.dN1_dy * sf2.N3
			  + sf3.dN1_dy * sf3.N3 + sf4.dN1_dy * sf4.N3) * gp_w;
	kmat[4][10] += mat_term;
	kmat[10][4] += mat_term;
	mat_term = (sf1.dN1_dy * sf1.N4 + sf2.dN1_dy * sf2.N4
			  + sf3.dN1_dy * sf3.N4 + sf4.dN1_dy * sf4.N4) * gp_w;
	kmat[4][11] += mat_term;
	kmat[11][4] += mat_term;
	mat_term = (sf1.dN2_dy * sf1.N1 + sf2.dN2_dy * sf2.N1
			  + sf3.dN2_dy * sf3.N1 + sf4.dN2_dy * sf4.N1) * gp_w;
	kmat[5][8] += mat_term;
	kmat[8][5] += mat_term;
	mat_term = (sf1.dN2_dy * sf1.N2 + sf2.dN2_dy * sf2.N2
			  + sf3.dN2_dy * sf3.N2 + sf4.dN2_dy * sf4.N2) * gp_w;
	kmat[5][9] += mat_term;
	kmat[9][5] += mat_term;
	mat_term = (sf1.dN2_dy * sf1.N3 + sf2.dN2_dy * sf2.N3
			  + sf3.dN2_dy * sf3.N3 + sf4.dN2_dy * sf4.N3) * gp_w;
	kmat[5][10] += mat_term;
	kmat[10][5] += mat_term;
	mat_term = (sf1.dN2_dy * sf1.N4 + sf2.dN2_dy * sf2.N4
			  + sf3.dN2_dy * sf3.N4 + sf4.dN2_dy * sf4.N4) * gp_w;
	kmat[5][11] += mat_term;
	kmat[11][5] += mat_term;
	mat_term = (sf1.dN3_dy * sf1.N1 + sf2.dN3_dy * sf2.N1
			  + sf3.dN3_dy * sf3.N1 + sf4.dN3_dy * sf4.N1) * gp_w;
	kmat[6][8] += mat_term;
	kmat[8][6] += mat_term;
	mat_term = (sf1.dN3_dy * sf1.N2 + sf2.dN3_dy * sf2.N2
			  + sf3.dN3_dy * sf3.N2 + sf4.dN3_dy * sf4.N2) * gp_w;
	kmat[6][9] += mat_term;
	kmat[9][6] += mat_term;
	mat_term = (sf1.dN3_dy * sf1.N3 + sf2.dN3_dy * sf2.N3
			  + sf3.dN3_dy * sf3.N3 + sf4.dN3_dy * sf4.N3) * gp_w;
	kmat[6][10] += mat_term;
	kmat[10][6] += mat_term;
	mat_term = (sf1.dN3_dy * sf1.N4 + sf2.dN3_dy * sf2.N4
			  + sf3.dN3_dy * sf3.N4 + sf4.dN3_dy * sf4.N4) * gp_w;
	kmat[6][11] += mat_term;
	kmat[11][6] += mat_term;
	mat_term = (sf1.dN4_dy * sf1.N1 + sf2.dN4_dy * sf2.N1
			  + sf3.dN4_dy * sf3.N1 + sf4.dN4_dy * sf4.N1) * gp_w;
	kmat[7][8] += mat_term;
	kmat[8][7] += mat_term;
	mat_term = (sf1.dN4_dy * sf1.N2 + sf2.dN4_dy * sf2.N2
			  + sf3.dN4_dy * sf3.N2 + sf4.dN4_dy * sf4.N2) * gp_w;
	kmat[7][9] += mat_term;
	kmat[9][7] += mat_term;
	mat_term = (sf1.dN4_dy * sf1.N3 + sf2.dN4_dy * sf2.N3
			  + sf3.dN4_dy * sf3.N3 + sf4.dN4_dy * sf4.N3) * gp_w;
	kmat[7][10] += mat_term;
	kmat[10][7] += mat_term;
	mat_term = (sf1.dN4_dy * sf1.N4 + sf2.dN4_dy * sf2.N4
			  + sf3.dN4_dy * sf3.N4 + sf4.dN4_dy * sf4.N4) * gp_w;
	kmat[7][11] += mat_term;
	kmat[11][7] += mat_term;
	// Ni * Nj / K
	double K = model->K;
	mat_term = (sf1.N1 * sf1.N1 / K + sf2.N1 * sf2.N1 / K
			  + sf3.N1 * sf3.N1 / K + sf4.N1 * sf4.N1 / K) * -gp_w;
	kmat[8][8] += mat_term;
	mat_term = (sf1.N1 * sf1.N2 / K + sf2.N1 * sf2.N2 / K
			  + sf3.N1 * sf3.N2 / K + sf4.N1 * sf4.N2 / K) * -gp_w;
	kmat[8][9] += mat_term;
	kmat[9][8] += mat_term;
	mat_term = (sf1.N1 * sf1.N3 / K + sf2.N1 * sf2.N3 / K
			  + sf3.N1 * sf3.N3 / K + sf4.N1 * sf4.N3 / K) * -gp_w;
	kmat[8][10] += mat_term;
	kmat[10][8] += mat_term;
	mat_term = (sf1.N1 * sf1.N4 / K + sf2.N1 * sf2.N4 / K
			  + sf3.N1 * sf3.N4 / K + sf4.N1 * sf4.N4 / K) * -gp_w;
	kmat[8][11] += mat_term;
	kmat[11][8] += mat_term;
	mat_term = (sf1.N2 * sf1.N2 / K + sf2.N2 * sf2.N2 / K
			  + sf3.N2 * sf3.N2 / K + sf4.N2 * sf4.N2 / K) * -gp_w;
	kmat[9][9] += mat_term;
	mat_term = (sf1.N2 * sf1.N3 / K + sf2.N2 * sf2.N3 / K
			  + sf3.N2 * sf3.N3 / K + sf4.N2 * sf4.N3 / K) * -gp_w;
	kmat[9][10] += mat_term;
	kmat[10][9] += mat_term;
	mat_term = (sf1.N2 * sf1.N4 / K + sf2.N2 * sf2.N4 / K
			  + sf3.N2 * sf3.N4 / K + sf4.N2 * sf4.N4 / K) * -gp_w;
	kmat[9][11] += mat_term;
	kmat[11][9] += mat_term;
	mat_term = (sf1.N3 * sf1.N3 / K + sf2.N3 * sf2.N3 / K
			  + sf3.N3 * sf3.N3 / K + sf4.N3 * sf4.N3 / K) * -gp_w;
	kmat[10][10] += mat_term;
	mat_term = (sf1.N3 * sf1.N4 / K + sf2.N3 * sf2.N4 / K
			  + sf3.N3 * sf3.N4 / K + sf4.N3 * sf4.N4 / K) * -gp_w;
	kmat[10][11] += mat_term;
	kmat[11][10] += mat_term;
	mat_term = (sf1.N4 * sf1.N4 / K + sf2.N4 * sf2.N4 / K
			  + sf3.N4 * sf3.N4 / K + sf4.N4 * sf4.N4 / K) * -gp_w;
	kmat[11][11] += mat_term;
	
	// force vector
	fvec[0] = (sf1.dN1_dx * (gp1.s11 + gp1.p) + sf1.dN1_dy * gp1.s12
			 + sf2.dN1_dx * (gp2.s11 + gp2.p) + sf2.dN1_dy * gp2.s12
			 + sf3.dN1_dx * (gp3.s11 + gp3.p) + sf3.dN1_dy * gp3.s12
			 + sf4.dN1_dx * (gp4.s11 + gp4.p) + sf4.dN1_dy * gp4.s12) * -gp_w;
	fvec[1] = (sf1.dN2_dx * (gp1.s11 + gp1.p) + sf1.dN2_dy * gp1.s12
			 + sf2.dN2_dx * (gp2.s11 + gp2.p) + sf2.dN2_dy * gp2.s12
			 + sf3.dN2_dx * (gp3.s11 + gp3.p) + sf3.dN2_dy * gp3.s12
			 + sf4.dN2_dx * (gp4.s11 + gp4.p) + sf4.dN2_dy * gp4.s12) * -gp_w;
	fvec[2] = (sf1.dN3_dx * (gp1.s11 + gp1.p) + sf1.dN3_dy * gp1.s12
			 + sf2.dN3_dx * (gp2.s11 + gp2.p) + sf2.dN3_dy * gp2.s12
			 + sf3.dN3_dx * (gp3.s11 + gp3.p) + sf3.dN3_dy * gp3.s12
			 + sf4.dN3_dx * (gp4.s11 + gp4.p) + sf4.dN3_dy * gp4.s12) * -gp_w;
	fvec[3] = (sf1.dN4_dx * (gp1.s11 + gp1.p) + sf1.dN4_dy * gp1.s12
			 + sf2.dN4_dx * (gp2.s11 + gp2.p) + sf2.dN4_dy * gp2.s12
			 + sf3.dN4_dx * (gp3.s11 + gp3.p) + sf3.dN4_dy * gp3.s12
			 + sf4.dN4_dx * (gp4.s11 + gp4.p) + sf4.dN4_dy * gp4.s12) * -gp_w;
	fvec[4] = (sf1.dN1_dx * gp1.s12 + sf1.dN1_dy * (gp1.s22 + gp1.p)
			 + sf2.dN1_dx * gp2.s12 + sf2.dN1_dy * (gp2.s22 + gp2.p)
			 + sf3.dN1_dx * gp3.s12 + sf3.dN1_dy * (gp3.s22 + gp3.p)
			 + sf4.dN1_dx * gp4.s12 + sf4.dN1_dy * (gp4.s22 + gp4.p)) * -gp_w;
	fvec[5] = (sf1.dN2_dx * gp1.s12 + sf1.dN2_dy * (gp1.s22 + gp1.p)
			 + sf2.dN2_dx * gp2.s12 + sf2.dN2_dy * (gp2.s22 + gp2.p)
			 + sf3.dN2_dx * gp3.s12 + sf3.dN2_dy * (gp3.s22 + gp3.p)
			 + sf4.dN2_dx * gp4.s12 + sf4.dN2_dy * (gp4.s22 + gp4.p)) * -gp_w;
	fvec[6] = (sf1.dN3_dx * gp1.s12 + sf1.dN3_dy * (gp1.s22 + gp1.p)
			 + sf2.dN3_dx * gp2.s12 + sf2.dN3_dy * (gp2.s22 + gp2.p)
			 + sf3.dN3_dx * gp3.s12 + sf3.dN3_dy * (gp3.s22 + gp3.p)
			 + sf4.dN3_dx * gp4.s12 + sf4.dN3_dy * (gp4.s22 + gp4.p)) * -gp_w;
	fvec[7] = (sf1.dN4_dx * gp1.s12 + sf1.dN4_dy * (gp1.s22 + gp1.p)
			 + sf2.dN4_dx * gp2.s12 + sf2.dN4_dy * (gp2.s22 + gp2.p)
			 + sf3.dN4_dx * gp3.s12 + sf3.dN4_dy * (gp3.s22 + gp3.p)
			 + sf4.dN4_dx * gp4.s12 + sf4.dN4_dy * (gp4.s22 + gp4.p)) * -gp_w;
	fvec[8] = 0.0;
	fvec[9] = 0.0;
	fvec[10] = 0.0;
	fvec[11] = 0.0;

	double coef_a = 1.0 / (2.0 * beta) - 1.0;
	double coef_v = 1.0 / (beta * dtime);
	for (Particle_mpm *pcl_iter = e.pcls; pcl_iter; pcl_iter = pcl_iter->next)
	{
		Particle_mpm &pcl = *pcl_iter;

		fvec[0] += pcl.N1 * pcl.m * pcl.ax * coef_a;
		fvec[1] += pcl.N2 * pcl.m * pcl.ax * coef_a;
		fvec[2] += pcl.N3 * pcl.m * pcl.ax * coef_a;
		fvec[3] += pcl.N4 * pcl.m * pcl.ax * coef_a;
		fvec[4] += pcl.N1 * pcl.m * pcl.ay * coef_a;
		fvec[5] += pcl.N2 * pcl.m * pcl.ay * coef_a;
		fvec[6] += pcl.N3 * pcl.m * pcl.ay * coef_a;
		fvec[7] += pcl.N4 * pcl.m * pcl.ay * coef_a;

		fvec[0] += pcl.N1 * pcl.m * pcl.vx * coef_v;
		fvec[1] += pcl.N2 * pcl.m * pcl.vx * coef_v;
		fvec[2] += pcl.N3 * pcl.m * pcl.vx * coef_v;
		fvec[3] += pcl.N4 * pcl.m * pcl.vx * coef_v;
		fvec[4] += pcl.N1 * pcl.m * pcl.vy * coef_v;
		fvec[5] += pcl.N2 * pcl.m * pcl.vy * coef_v;
		fvec[6] += pcl.N3 * pcl.m * pcl.vy * coef_v;
		fvec[7] += pcl.N4 * pcl.m * pcl.vy * coef_v;
	}
	//print_vec(fvec, 12);
}
