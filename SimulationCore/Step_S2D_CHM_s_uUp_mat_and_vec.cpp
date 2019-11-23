#include "SimulationCore_pcp.h"

#define Keep_Newmark_Coefficients
#include "Step_S2D_CHM_s_uUp.h"

namespace
{

void form_E_matrix(double E_mat[3][3], double E, double niu)
{
	double E_tmp = E / ((1.0 + niu) * (1.0 - 2.0 * niu));
	E_mat[0][0] = E_tmp * (1.0 - niu);
	E_mat[0][1] = E_tmp * niu;
	E_mat[0][2] = 0.0;
	E_mat[1][0] = E_tmp * niu;
	E_mat[1][1] = E_tmp * (1.0 - niu);
	E_mat[1][2] = 0.0;
	E_mat[2][0] = 0.0;
	E_mat[2][1] = 0.0;
	E_mat[2][2] = E / (2.0 * (1.0 + niu));
}

};

void Step_S2D_CHM_s_uUp::add_internal_force_vec(double fvec[20], Particle &pcl)
{
	Element &e = *pcl.pe;
	Node &n1 = model->nodes[e.n1_id];
	Node &n2 = model->nodes[e.n2_id];
	Node &n3 = model->nodes[e.n3_id];
	Node &n4 = model->nodes[e.n4_id];
	double pcl_p = n1.p * pcl.N1 + n2.p * pcl.N2 + n3.p * pcl.N3 + n4.p * pcl.N4;
	// solid x
	fvec[0] -= (pcl.dN1_dx * (pcl.s11 - (1.0 - pcl.n) * pcl_p) + pcl.dN1_dy * pcl.s12) * pcl.vol;
	fvec[1] -= (pcl.dN2_dx * (pcl.s11 - (1.0 - pcl.n) * pcl_p) + pcl.dN2_dy * pcl.s12) * pcl.vol;
	fvec[2] -= (pcl.dN3_dx * (pcl.s11 - (1.0 - pcl.n) * pcl_p) + pcl.dN3_dy * pcl.s12) * pcl.vol;
	fvec[3] -= (pcl.dN4_dx * (pcl.s11 - (1.0 - pcl.n) * pcl_p) + pcl.dN4_dy * pcl.s12) * pcl.vol;
	// solid y
	fvec[4] -= (pcl.dN1_dx * pcl.s12 + pcl.dN1_dy * (pcl.s22 - (1.0 - pcl.n) * pcl_p)) * pcl.vol;
	fvec[5] -= (pcl.dN2_dx * pcl.s12 + pcl.dN2_dy * (pcl.s22 - (1.0 - pcl.n) * pcl_p)) * pcl.vol;
	fvec[6] -= (pcl.dN3_dx * pcl.s12 + pcl.dN3_dy * (pcl.s22 - (1.0 - pcl.n) * pcl_p)) * pcl.vol;
	fvec[7] -= (pcl.dN4_dx * pcl.s12 + pcl.dN4_dy * (pcl.s22 - (1.0 - pcl.n) * pcl_p)) * pcl.vol;
	// fluid x
	fvec[8]  -= pcl.dN1_dx * pcl.n * -pcl_p * pcl.vol;
	fvec[9]  -= pcl.dN2_dx * pcl.n * -pcl_p * pcl.vol;
	fvec[10] -= pcl.dN3_dx * pcl.n * -pcl_p * pcl.vol;
	fvec[11] -= pcl.dN4_dx * pcl.n * -pcl_p * pcl.vol;
	// fluid y
	fvec[12] -= pcl.dN1_dy * pcl.n * -pcl_p * pcl.vol;
	fvec[13] -= pcl.dN2_dy * pcl.n * -pcl_p * pcl.vol;
	fvec[14] -= pcl.dN3_dy * pcl.n * -pcl_p * pcl.vol;
	fvec[15] -= pcl.dN4_dy * pcl.n * -pcl_p * pcl.vol;
}


void Step_S2D_CHM_s_uUp::form_elem_stiffness_mat_and_force_vec(
	Element &e, double kmat[20][20], double fvec[20])
{
	// kmat
	double coef_m = 1.0 / (beta * dtime * dtime);
	double coef_c = gamma / (beta * dtime);
	double mat_term;
	double E_mat[3][3], dN_dx[3][8], mid_mat[3][8];
	memset(kmat, 0, sizeof(double) * 20 * 20);
	// fvec
	double coef_ma = 1.0 / (2.0 * beta) - 1.0;
	double coef_ca = -dtime * (1.0 - gamma / (2.0 * beta));
	double coef_mv = 1.0 / (beta * dtime);
	double coef_cv = gamma / beta - 1.0;
	Node &n1 = model->nodes[e.n1_id];
	Node &n2 = model->nodes[e.n2_id];
	Node &n3 = model->nodes[e.n3_id];
	Node &n4 = model->nodes[e.n4_id];
	double vec_term1, vec_term2;
	double pcl_m_f, n2_miu_div_k_vol;
	double pcl_ax_s, pcl_ay_s, pcl_ax_f, pcl_ay_f;
	double pcl_max_s, pcl_may_s, pcl_max_f, pcl_may_f;
	double pcl_n2_miu_div_k_vol_ax, pcl_n2_miu_div_k_vol_ay;
	double pcl_vx_s, pcl_vy_s, pcl_vx_f, pcl_vy_f;
	double pcl_mvx_s, pcl_mvy_s, pcl_mvx_f, pcl_mvy_f;
	double pcl_n2_miu_div_k_vol_vx, pcl_n2_miu_div_k_vol_vy;
	memset(fvec, 0, sizeof(double) * 20);

	for (Particle *pcl_iter = e.pcls; pcl_iter; pcl_iter = pcl_iter->next)
	{
		Particle &pcl = *pcl_iter;

		// mass matrix
		// solid
		mat_term = pcl.N1 * pcl.m_s * pcl.N1 * coef_m;
		kmat[0][0] += mat_term;
		kmat[4][4] += mat_term;
		mat_term = pcl.N1 * pcl.m_s * pcl.N2 * coef_m;
		kmat[0][1] += mat_term;
		kmat[1][0] += mat_term;
		kmat[4][5] += mat_term;
		kmat[5][4] += mat_term;
		mat_term = pcl.N1 * pcl.m_s * pcl.N3 * coef_m;
		kmat[0][2] += mat_term;
		kmat[2][0] += mat_term;
		kmat[4][6] += mat_term;
		kmat[6][4] += mat_term;
		mat_term = pcl.N1 * pcl.m_s * pcl.N4 * coef_m;
		kmat[0][3] += mat_term;
		kmat[3][0] += mat_term;
		kmat[4][7] += mat_term;
		kmat[7][4] += mat_term;
		mat_term = pcl.N2 * pcl.m_s * pcl.N2 * coef_m;
		kmat[1][1] += mat_term;
		kmat[5][5] += mat_term;
		mat_term = pcl.N2 * pcl.m_s * pcl.N3 * coef_m;
		kmat[1][2] += mat_term;
		kmat[2][1] += mat_term;
		kmat[5][6] += mat_term;
		kmat[6][5] += mat_term;
		mat_term = pcl.N2 * pcl.m_s * pcl.N4 * coef_m;
		kmat[1][3] += mat_term;
		kmat[3][1] += mat_term;
		kmat[5][7] += mat_term;
		kmat[7][5] += mat_term;
		mat_term = pcl.N3 * pcl.m_s * pcl.N3 * coef_m;
		kmat[2][2] += mat_term;
		kmat[6][6] += mat_term;
		mat_term = pcl.N3 * pcl.m_s * pcl.N4 * coef_m;
		kmat[2][3] += mat_term;
		kmat[3][2] += mat_term;
		kmat[6][7] += mat_term;
		kmat[7][6] += mat_term;
		mat_term = pcl.N4 * pcl.m_s * pcl.N4 * coef_m;
		kmat[3][3] += mat_term;
		kmat[7][7] += mat_term;
		// fluid
		pcl_m_f = pcl.n * pcl.density_f * pcl.vol;
		mat_term = pcl.N1 * pcl_m_f * pcl.N1 * coef_m;
		kmat[8][8]   += mat_term;
		kmat[12][12] += mat_term;
		mat_term = pcl.N1 * pcl_m_f * pcl.N2 * coef_m;
		kmat[8][9]   += mat_term;
		kmat[9][8]   += mat_term;
		kmat[12][13] += mat_term;
		kmat[13][12] += mat_term;
		mat_term = pcl.N1 * pcl_m_f * pcl.N3 * coef_m;
		kmat[8][10]  += mat_term;
		kmat[10][8]  += mat_term;
		kmat[12][14] += mat_term;
		kmat[14][12] += mat_term;
		mat_term = pcl.N1 * pcl_m_f * pcl.N4 * coef_m;
		kmat[8][11]  += mat_term;
		kmat[11][8]  += mat_term;
		kmat[12][15] += mat_term;
		kmat[15][12] += mat_term;
		mat_term = pcl.N2 * pcl_m_f * pcl.N2 * coef_m;
		kmat[9][9]   += mat_term;
		kmat[13][13] += mat_term;
		mat_term = pcl.N2 * pcl_m_f * pcl.N3 * coef_m;
		kmat[9][10]  += mat_term;
		kmat[10][9]  += mat_term;
		kmat[13][14] += mat_term;
		kmat[14][13] += mat_term;
		mat_term = pcl.N2 * pcl_m_f * pcl.N4 * coef_m;
		kmat[9][11]  += mat_term;
		kmat[11][9]  += mat_term;
		kmat[13][15] += mat_term;
		kmat[15][13] += mat_term;
		mat_term = pcl.N3 * pcl_m_f * pcl.N3 * coef_m;
		kmat[10][10] += mat_term;
		kmat[14][14] += mat_term;
		mat_term = pcl.N3 * pcl_m_f * pcl.N4 * coef_m;
		kmat[10][11] += mat_term;
		kmat[11][10] += mat_term;
		kmat[14][15] += mat_term;
		kmat[15][14] += mat_term;
		mat_term = pcl.N4 * pcl_m_f * pcl.N4 * coef_m;
		kmat[11][11] += mat_term;
		kmat[15][15] += mat_term;

		// damping matrix
		n2_miu_div_k_vol = pcl.n * pcl.n * pcl.miu / pcl.k * pcl.vol;
#define FILL_SEEP_MAT(x_id, y_id)           \
		kmat[x_id][y_id]       += mat_term; \
		kmat[x_id][y_id+8]     -= mat_term; \
		kmat[x_id+4][y_id+4]   += mat_term; \
		kmat[x_id+4][y_id+12]  -= mat_term; \
		kmat[x_id+8][y_id]     -= mat_term; \
		kmat[x_id+8][y_id+8]   += mat_term; \
		kmat[x_id+12][y_id+4]  -= mat_term; \
		kmat[x_id+12][y_id+12] += mat_term
		mat_term = pcl.N1 * n2_miu_div_k_vol * pcl.N1 * coef_c;
		FILL_SEEP_MAT(0, 0);
		mat_term = pcl.N1 * n2_miu_div_k_vol * pcl.N2 * coef_c;
		FILL_SEEP_MAT(0, 1);
		FILL_SEEP_MAT(1, 0);
		mat_term = pcl.N1 * n2_miu_div_k_vol * pcl.N3 * coef_c;
		FILL_SEEP_MAT(0, 2);
		FILL_SEEP_MAT(2, 0);
		mat_term = pcl.N1 * n2_miu_div_k_vol * pcl.N4 * coef_c;
		FILL_SEEP_MAT(0, 3);
		FILL_SEEP_MAT(3, 0);
		mat_term = pcl.N2 * n2_miu_div_k_vol * pcl.N2 * coef_c;
		FILL_SEEP_MAT(1, 1);
		mat_term = pcl.N2 * n2_miu_div_k_vol * pcl.N3 * coef_c;
		FILL_SEEP_MAT(1, 2);
		FILL_SEEP_MAT(2, 1);
		mat_term = pcl.N2 * n2_miu_div_k_vol * pcl.N4 * coef_c;
		FILL_SEEP_MAT(1, 3);
		FILL_SEEP_MAT(3, 1);
		mat_term = pcl.N3 * n2_miu_div_k_vol * pcl.N3 * coef_c;
		FILL_SEEP_MAT(2, 2);
		mat_term = pcl.N3 * n2_miu_div_k_vol * pcl.N4 * coef_c;
		FILL_SEEP_MAT(2, 3);
		FILL_SEEP_MAT(3, 2);
		mat_term = pcl.N4 * n2_miu_div_k_vol * pcl.N4 * coef_c;
		FILL_SEEP_MAT(3, 3);

		// stiffness matrix
		// BDB[8][8]
		form_E_matrix(E_mat, pcl.E, pcl.niu);
		dN_dx[0][0] = pcl.dN1_dx;
		dN_dx[0][1] = pcl.dN2_dx;
		dN_dx[0][2] = pcl.dN3_dx;
		dN_dx[0][3] = pcl.dN4_dx;
		dN_dx[0][4] = 0.0;
		dN_dx[0][5] = 0.0;
		dN_dx[0][6] = 0.0;
		dN_dx[0][7] = 0.0;
		dN_dx[1][0] = 0.0;
		dN_dx[1][1] = 0.0;
		dN_dx[1][2] = 0.0;
		dN_dx[1][3] = 0.0;
		dN_dx[1][4] = pcl.dN1_dy;
		dN_dx[1][5] = pcl.dN2_dy;
		dN_dx[1][6] = pcl.dN3_dy;
		dN_dx[1][7] = pcl.dN4_dy;
		dN_dx[2][0] = pcl.dN1_dy;
		dN_dx[2][1] = pcl.dN2_dy;
		dN_dx[2][2] = pcl.dN3_dy;
		dN_dx[2][3] = pcl.dN4_dy;
		dN_dx[2][4] = pcl.dN1_dx;
		dN_dx[2][5] = pcl.dN2_dx;
		dN_dx[2][6] = pcl.dN3_dx;
		dN_dx[2][7] = pcl.dN4_dx;
		for (size_t i = 0; i < 3; ++i)
			for (size_t j = 0; j < 8; ++j)
				mid_mat[i][j] = E_mat[i][0] * dN_dx[0][j]
							  + E_mat[i][1] * dN_dx[1][j]
							  + E_mat[i][2] * dN_dx[2][j];
		for (size_t i = 0; i < 8; ++i)
			for (size_t j = 0; j < 8; ++j)
				kmat[i][j] += (dN_dx[0][i] * mid_mat[0][j]
							 + dN_dx[1][i] * mid_mat[1][j]
							 + dN_dx[2][i] * mid_mat[2][j]) * pcl.vol;
		// Gsx
		mat_term = -pcl.dN1_dx * (1.0 - pcl.n) * pcl.N1 * pcl.vol;
		kmat[0][16] += mat_term;
		kmat[16][0] += mat_term;
		mat_term = -pcl.dN1_dx * (1.0 - pcl.n) * pcl.N2 * pcl.vol;
		kmat[0][17] += mat_term;
		kmat[17][0] += mat_term;
		mat_term = -pcl.dN1_dx * (1.0 - pcl.n) * pcl.N3 * pcl.vol;
		kmat[0][18] += mat_term;
		kmat[18][0] += mat_term;
		mat_term = -pcl.dN1_dx * (1.0 - pcl.n) * pcl.N4 * pcl.vol;
		kmat[0][19] += mat_term;
		kmat[19][0] += mat_term;
		mat_term = -pcl.dN2_dx * (1.0 - pcl.n) * pcl.N1 * pcl.vol;
		kmat[1][16] += mat_term;
		kmat[16][1] += mat_term;
		mat_term = -pcl.dN2_dx * (1.0 - pcl.n) * pcl.N2 * pcl.vol;
		kmat[1][17] += mat_term;
		kmat[17][1] += mat_term;
		mat_term = -pcl.dN2_dx * (1.0 - pcl.n) * pcl.N3 * pcl.vol;
		kmat[1][18] += mat_term;
		kmat[18][1] += mat_term;
		mat_term = -pcl.dN2_dx * (1.0 - pcl.n) * pcl.N4 * pcl.vol;
		kmat[1][19] += mat_term;
		kmat[19][1] += mat_term;
		mat_term = -pcl.dN3_dx * (1.0 - pcl.n) * pcl.N1 * pcl.vol;
		kmat[2][16] += mat_term;
		kmat[16][2] += mat_term;
		mat_term = -pcl.dN3_dx * (1.0 - pcl.n) * pcl.N2 * pcl.vol;
		kmat[2][17] += mat_term;
		kmat[17][2] += mat_term;
		mat_term = -pcl.dN3_dx * (1.0 - pcl.n) * pcl.N3 * pcl.vol;
		kmat[2][18] += mat_term;
		kmat[18][2] += mat_term;
		mat_term = -pcl.dN3_dx * (1.0 - pcl.n) * pcl.N4 * pcl.vol;
		kmat[2][19] += mat_term;
		kmat[19][2] += mat_term;
		mat_term = -pcl.dN4_dx * (1.0 - pcl.n) * pcl.N1 * pcl.vol;
		kmat[3][16] += mat_term;
		kmat[16][3] += mat_term;
		mat_term = -pcl.dN4_dx * (1.0 - pcl.n) * pcl.N2 * pcl.vol;
		kmat[3][17] += mat_term;
		kmat[17][3] += mat_term;
		mat_term = -pcl.dN4_dx * (1.0 - pcl.n) * pcl.N3 * pcl.vol;
		kmat[3][18] += mat_term;
		kmat[18][3] += mat_term;
		mat_term = -pcl.dN4_dx * (1.0 - pcl.n) * pcl.N4 * pcl.vol;
		kmat[3][19] += mat_term;
		kmat[19][3] += mat_term;
		// Gsy
		mat_term = -pcl.dN1_dy * (1.0 - pcl.n) * pcl.N1 * pcl.vol;
		kmat[4][16] += mat_term;
		kmat[16][4] += mat_term;
		mat_term = -pcl.dN1_dy * (1.0 - pcl.n) * pcl.N2 * pcl.vol;
		kmat[4][17] += mat_term;
		kmat[17][4] += mat_term;
		mat_term = -pcl.dN1_dy * (1.0 - pcl.n) * pcl.N3 * pcl.vol;
		kmat[4][18] += mat_term;
		kmat[18][4] += mat_term;
		mat_term = -pcl.dN1_dy * (1.0 - pcl.n) * pcl.N4 * pcl.vol;
		kmat[4][19] += mat_term;
		kmat[19][4] += mat_term;
		mat_term = -pcl.dN2_dy * (1.0 - pcl.n) * pcl.N1 * pcl.vol;
		kmat[5][16] += mat_term;
		kmat[16][5] += mat_term;
		mat_term = -pcl.dN2_dy * (1.0 - pcl.n) * pcl.N2 * pcl.vol;
		kmat[5][17] += mat_term;
		kmat[17][5] += mat_term;
		mat_term = -pcl.dN2_dy * (1.0 - pcl.n) * pcl.N3 * pcl.vol;
		kmat[5][18] += mat_term;
		kmat[18][5] += mat_term;
		mat_term = -pcl.dN2_dy * (1.0 - pcl.n) * pcl.N4 * pcl.vol;
		kmat[5][19] += mat_term;
		kmat[19][5] += mat_term;
		mat_term = -pcl.dN3_dy * (1.0 - pcl.n) * pcl.N1 * pcl.vol;
		kmat[6][16] += mat_term;
		kmat[16][6] += mat_term;
		mat_term = -pcl.dN3_dy * (1.0 - pcl.n) * pcl.N2 * pcl.vol;
		kmat[6][17] += mat_term;
		kmat[17][6] += mat_term;
		mat_term = -pcl.dN3_dy * (1.0 - pcl.n) * pcl.N3 * pcl.vol;
		kmat[6][18] += mat_term;
		kmat[18][6] += mat_term;
		mat_term = -pcl.dN3_dy * (1.0 - pcl.n) * pcl.N4 * pcl.vol;
		kmat[6][19] += mat_term;
		kmat[19][6] += mat_term;
		mat_term = -pcl.dN4_dy * (1.0 - pcl.n) * pcl.N1 * pcl.vol;
		kmat[7][16] += mat_term;
		kmat[16][7] += mat_term;
		mat_term = -pcl.dN4_dy * (1.0 - pcl.n) * pcl.N2 * pcl.vol;
		kmat[7][17] += mat_term;
		kmat[17][7] += mat_term;
		mat_term = -pcl.dN4_dy * (1.0 - pcl.n) * pcl.N3 * pcl.vol;
		kmat[7][18] += mat_term;
		kmat[18][7] += mat_term;
		mat_term = -pcl.dN4_dy * (1.0 - pcl.n) * pcl.N4 * pcl.vol;
		kmat[7][19] += mat_term;
		kmat[19][7] += mat_term;
		// Gfx
		mat_term = pcl.dN1_dx * pcl.n * pcl.N1 * -pcl.vol;
		kmat[8][16] += mat_term;
		kmat[16][8] += mat_term;
		mat_term = pcl.dN1_dx * pcl.n * pcl.N2 * -pcl.vol;
		kmat[8][17] += mat_term;
		kmat[17][8] += mat_term;
		mat_term = pcl.dN1_dx * pcl.n * pcl.N3 * -pcl.vol;
		kmat[8][18] += mat_term;
		kmat[18][8] += mat_term;
		mat_term = pcl.dN1_dx * pcl.n * pcl.N4 * -pcl.vol;
		kmat[8][19] += mat_term;
		kmat[19][8] += mat_term;
		mat_term = pcl.dN2_dx * pcl.n * pcl.N1 * -pcl.vol;
		kmat[9][16] += mat_term;
		kmat[16][9] += mat_term;
		mat_term = pcl.dN2_dx * pcl.n * pcl.N2 * -pcl.vol;
		kmat[9][17] += mat_term;
		kmat[17][9] += mat_term;
		mat_term = pcl.dN2_dx * pcl.n * pcl.N3 * -pcl.vol;
		kmat[9][18] += mat_term;
		kmat[18][9] += mat_term;
		mat_term = pcl.dN2_dx * pcl.n * pcl.N4 * -pcl.vol;
		kmat[9][19] += mat_term;
		kmat[19][9] += mat_term;
		mat_term = pcl.dN3_dx * pcl.n * pcl.N1 * -pcl.vol;
		kmat[10][16] += mat_term;
		kmat[16][10] += mat_term;
		mat_term = pcl.dN3_dx * pcl.n * pcl.N2 * -pcl.vol;
		kmat[10][17] += mat_term;
		kmat[17][10] += mat_term;
		mat_term = pcl.dN3_dx * pcl.n * pcl.N3 * -pcl.vol;
		kmat[10][18] += mat_term;
		kmat[18][10] += mat_term;
		mat_term = pcl.dN3_dx * pcl.n * pcl.N4 * -pcl.vol;
		kmat[10][19] += mat_term;
		kmat[19][10] += mat_term;
		mat_term = pcl.dN4_dx * pcl.n * pcl.N1 * -pcl.vol;
		kmat[11][16] += mat_term;
		kmat[16][11] += mat_term;
		mat_term = pcl.dN4_dx * pcl.n * pcl.N2 * -pcl.vol;
		kmat[11][17] += mat_term;
		kmat[17][11] += mat_term;
		mat_term = pcl.dN4_dx * pcl.n * pcl.N3 * -pcl.vol;
		kmat[11][18] += mat_term;
		kmat[18][11] += mat_term;
		mat_term = pcl.dN4_dx * pcl.n * pcl.N4 * -pcl.vol;
		kmat[11][19] += mat_term;
		kmat[19][11] += mat_term;
		// Gfy
		mat_term = pcl.dN1_dy * pcl.n * pcl.N1 * -pcl.vol;
		kmat[12][16] += mat_term;
		kmat[16][12] += mat_term;
		mat_term = pcl.dN1_dy * pcl.n * pcl.N2 * -pcl.vol;
		kmat[12][17] += mat_term;
		kmat[17][12] += mat_term;
		mat_term = pcl.dN1_dy * pcl.n * pcl.N3 * -pcl.vol;
		kmat[12][18] += mat_term;
		kmat[18][12] += mat_term;
		mat_term = pcl.dN1_dy * pcl.n * pcl.N4 * -pcl.vol;
		kmat[12][19] += mat_term;
		kmat[19][12] += mat_term;
		mat_term = pcl.dN2_dy * pcl.n * pcl.N1 * -pcl.vol;
		kmat[13][16] += mat_term;
		kmat[16][13] += mat_term;
		mat_term = pcl.dN2_dy * pcl.n * pcl.N2 * -pcl.vol;
		kmat[13][17] += mat_term;
		kmat[17][13] += mat_term;
		mat_term = pcl.dN2_dy * pcl.n * pcl.N3 * -pcl.vol;
		kmat[13][18] += mat_term;
		kmat[18][13] += mat_term;
		mat_term = pcl.dN2_dy * pcl.n * pcl.N4 * -pcl.vol;
		kmat[13][19] += mat_term;
		kmat[19][13] += mat_term;
		mat_term = pcl.dN3_dy * pcl.n * pcl.N1 * -pcl.vol;
		kmat[14][16] += mat_term;
		kmat[16][14] += mat_term;
		mat_term = pcl.dN3_dy * pcl.n * pcl.N2 * -pcl.vol;
		kmat[14][17] += mat_term;
		kmat[17][14] += mat_term;
		mat_term = pcl.dN3_dy * pcl.n * pcl.N3 * -pcl.vol;
		kmat[14][18] += mat_term;
		kmat[18][14] += mat_term;
		mat_term = pcl.dN3_dy * pcl.n * pcl.N4 * -pcl.vol;
		kmat[14][19] += mat_term;
		kmat[19][14] += mat_term;
		mat_term = pcl.dN4_dy * pcl.n * pcl.N1 * -pcl.vol;
		kmat[15][16] += mat_term;
		kmat[16][15] += mat_term;
		mat_term = pcl.dN4_dy * pcl.n * pcl.N2 * -pcl.vol;
		kmat[15][17] += mat_term;
		kmat[17][15] += mat_term;
		mat_term = pcl.dN4_dy * pcl.n * pcl.N3 * -pcl.vol;
		kmat[15][18] += mat_term;
		kmat[18][15] += mat_term;
		mat_term = pcl.dN4_dy * pcl.n * pcl.N4 * -pcl.vol;
		kmat[15][19] += mat_term;
		kmat[19][15] += mat_term;
		// P
		double n_div_Kf = pcl.n / pcl.Kf;
		mat_term = pcl.N1 * n_div_Kf * pcl.N1 * -pcl.vol;
		kmat[16][16] += mat_term;
		mat_term = pcl.N1 * n_div_Kf * pcl.N2 * -pcl.vol;
		kmat[16][17] += mat_term;
		kmat[17][16] += mat_term;
		mat_term = pcl.N1 * n_div_Kf * pcl.N3 * -pcl.vol;
		kmat[16][18] += mat_term;
		kmat[18][16] += mat_term;
		mat_term = pcl.N1 * n_div_Kf * pcl.N4 * -pcl.vol;
		kmat[16][19] += mat_term;
		kmat[19][16] += mat_term;
		mat_term = pcl.N2 * n_div_Kf * pcl.N2 * -pcl.vol;
		kmat[17][17] += mat_term;
		mat_term = pcl.N2 * n_div_Kf * pcl.N3 * -pcl.vol;
		kmat[17][18] += mat_term;
		kmat[18][17] += mat_term;
		mat_term = pcl.N2 * n_div_Kf * pcl.N4 * -pcl.vol;
		kmat[17][19] += mat_term;
		kmat[19][17] += mat_term;
		mat_term = pcl.N3 * n_div_Kf * pcl.N3 * -pcl.vol;
		kmat[18][18] += mat_term;
		mat_term = pcl.N3 * n_div_Kf * pcl.N4 * -pcl.vol;
		kmat[18][19] += mat_term;
		kmat[19][18] += mat_term;
		mat_term = pcl.N4 * n_div_Kf * pcl.N4 * -pcl.vol;
		kmat[19][19] += mat_term;
				
		// Elemental force vector
		// f_int
		add_internal_force_vec(fvec, pcl);
		
		// a
		pcl_ax_s = pcl.N1 * n1.ax_s + pcl.N2 * n2.ax_s + pcl.N3 * n3.ax_s + pcl.N4 * n4.ax_s;
		pcl_ay_s = pcl.N1 * n1.ay_s + pcl.N2 * n2.ay_s + pcl.N3 * n3.ay_s + pcl.N4 * n4.ay_s;
		pcl_ax_f = pcl.N1 * n1.ax_f + pcl.N2 * n2.ax_f + pcl.N3 * n3.ax_f + pcl.N4 * n4.ax_f;
		pcl_ay_f = pcl.N1 * n1.ay_f + pcl.N2 * n2.ay_f + pcl.N3 * n3.ay_f + pcl.N4 * n4.ay_f;
		pcl_max_s = pcl.m_s * pcl_ax_s;
		pcl_may_s = pcl.m_s * pcl_ay_s;
		pcl_max_f = pcl_m_f * pcl_ax_f;
		pcl_may_f = pcl_m_f * pcl_ay_f;
		pcl_n2_miu_div_k_vol_ax = n2_miu_div_k_vol * (pcl_ax_s - pcl_ax_f);
		pcl_n2_miu_div_k_vol_ay = n2_miu_div_k_vol * (pcl_ay_s - pcl_ay_f);
		// solid x
		vec_term1 = pcl.N1 * pcl_max_s;
		vec_term2 = pcl.N1 * pcl_n2_miu_div_k_vol_ax;
		fvec[0] += coef_ma * vec_term1 + coef_ca * vec_term2;
		vec_term1 = pcl.N2 * pcl_max_s;
		vec_term2 = pcl.N2 * pcl_n2_miu_div_k_vol_ax;
		fvec[1] += coef_ma * vec_term1 + coef_ca * vec_term2;
		vec_term1 = pcl.N3 * pcl_max_s;
		vec_term2 = pcl.N3 * pcl_n2_miu_div_k_vol_ax;
		fvec[2] += coef_ma * vec_term1 + coef_ca * vec_term2;
		vec_term1 = pcl.N4 * pcl_max_s;
		vec_term2 = pcl.N4 * pcl_n2_miu_div_k_vol_ax;
		fvec[3] += coef_ma * vec_term1 + coef_ca * vec_term2;
		// solid y
		vec_term1 = pcl.N1 * pcl_may_s;
		vec_term2 = pcl.N1 * pcl_n2_miu_div_k_vol_ay;
		fvec[4] += coef_ma * vec_term1 + coef_ca * vec_term2;
		vec_term1 = pcl.N2 * pcl_may_s;
		vec_term2 = pcl.N2 * pcl_n2_miu_div_k_vol_ay;
		fvec[5] += coef_ma * vec_term1 + coef_ca * vec_term2;
		vec_term1 = pcl.N3 * pcl_may_s;
		vec_term2 = pcl.N3 * pcl_n2_miu_div_k_vol_ay;
		fvec[6] += coef_ma * vec_term1 + coef_ca * vec_term2;
		vec_term1 = pcl.N4 * pcl_may_s;
		vec_term2 = pcl.N4 * pcl_n2_miu_div_k_vol_ay;
		fvec[7] += coef_ma * vec_term1 + coef_ca * vec_term2;
		// fluid x
		vec_term1 = pcl.N1 * pcl_max_f;
		vec_term2 = pcl.N1 * -pcl_n2_miu_div_k_vol_ax;
		fvec[8] += coef_ma * vec_term1 + coef_ca * vec_term2;
		vec_term1 = pcl.N2 * pcl_max_f;
		vec_term2 = pcl.N2 * -pcl_n2_miu_div_k_vol_ax;
		fvec[9] += coef_ma * vec_term1 + coef_ca * vec_term2;
		vec_term1 = pcl.N3 * pcl_max_f;
		vec_term2 = pcl.N3 * -pcl_n2_miu_div_k_vol_ax;
		fvec[10] += coef_ma * vec_term1 + coef_ca * vec_term2;
		vec_term1 = pcl.N4 * pcl_max_f;
		vec_term2 = pcl.N4 * -pcl_n2_miu_div_k_vol_ax;
		fvec[11] += coef_ma * vec_term1 + coef_ca * vec_term2;
		// fluid y
		vec_term1 = pcl.N1 * pcl_may_f;
		vec_term2 = pcl.N1 * -pcl_n2_miu_div_k_vol_ay;
		fvec[12] += coef_ma * vec_term1 + coef_ca * vec_term2;
		vec_term1 = pcl.N2 * pcl_may_f;
		vec_term2 = pcl.N2 * -pcl_n2_miu_div_k_vol_ay;
		fvec[13] += coef_ma * vec_term1 + coef_ca * vec_term2;
		vec_term1 = pcl.N3 * pcl_may_f;
		vec_term2 = pcl.N3 * -pcl_n2_miu_div_k_vol_ay; 
		fvec[14] += coef_ma * vec_term1 + coef_ca * vec_term2;
		vec_term1 = pcl.N4 * pcl_may_f; 
		vec_term2 = pcl.N4 * -pcl_n2_miu_div_k_vol_ay;
		fvec[15] += coef_ma * vec_term1 + coef_ca * vec_term2;

		// v
		pcl_vx_s = pcl.N1 * n1.vx_s + pcl.N2 * n2.vx_s + pcl.N3 * n3.vx_s + pcl.N4 * n4.vx_s;
		pcl_vy_s = pcl.N1 * n1.vy_s + pcl.N2 * n2.vy_s + pcl.N3 * n3.vy_s + pcl.N4 * n4.vy_s;
		pcl_vx_f = pcl.N1 * n1.vx_f + pcl.N2 * n2.vx_f + pcl.N3 * n3.vx_f + pcl.N4 * n4.vx_f;
		pcl_vy_f = pcl.N1 * n1.vy_f + pcl.N2 * n2.vy_f + pcl.N3 * n3.vy_f + pcl.N4 * n4.vy_f;
		pcl_mvx_s = pcl.m_s * pcl_vx_s;
		pcl_mvy_s = pcl.m_s * pcl_vy_s;
		pcl_mvx_f = pcl_m_f * pcl_vx_f;
		pcl_mvy_f = pcl_m_f * pcl_vy_f;
		pcl_n2_miu_div_k_vol_vx = n2_miu_div_k_vol * (pcl_vx_s - pcl_vx_f);
		pcl_n2_miu_div_k_vol_vy = n2_miu_div_k_vol * (pcl_vy_s - pcl_vy_f);
		// solid x
		vec_term1 = pcl.N1 * pcl_mvx_s;
		vec_term2 = pcl.N1 * pcl_n2_miu_div_k_vol_vx;
		fvec[0] += coef_mv * vec_term1 + coef_cv * vec_term2;
		vec_term1 = pcl.N2 * pcl_mvx_s;
		vec_term2 = pcl.N2 * pcl_n2_miu_div_k_vol_vx;
		fvec[1] += coef_mv * vec_term1 + coef_cv * vec_term2;
		vec_term1 = pcl.N3 * pcl_mvx_s;
		vec_term2 = pcl.N3 * pcl_n2_miu_div_k_vol_vx;
		fvec[2] += coef_mv * vec_term1 + coef_cv * vec_term2;
		vec_term1 = pcl.N4 * pcl_mvx_s;
		vec_term2 = pcl.N4 * pcl_n2_miu_div_k_vol_vx;
		fvec[3] += coef_mv * vec_term1 + coef_cv * vec_term2;
		// solid y
		vec_term1 = pcl.N1 * pcl_mvy_s;
		vec_term2 = pcl.N1 * pcl_n2_miu_div_k_vol_vy;
		fvec[4] += coef_mv * vec_term1 + coef_cv * vec_term2;
		vec_term1 = pcl.N2 * pcl_mvy_s;
		vec_term2 = pcl.N2 * pcl_n2_miu_div_k_vol_vy;
		fvec[5] += coef_mv * vec_term1 + coef_cv * vec_term2;
		vec_term1 = pcl.N3 * pcl_mvy_s;
		vec_term2 = pcl.N3 * pcl_n2_miu_div_k_vol_vy;
		fvec[6] += coef_mv * vec_term1 + coef_cv * vec_term2;
		vec_term1 = pcl.N4 * pcl_mvy_s;
		vec_term2 = pcl.N4 * pcl_n2_miu_div_k_vol_vy;
		fvec[7] += coef_mv * vec_term1 + coef_cv * vec_term2;
		// fluid x
		vec_term1 = pcl.N1 * pcl_mvx_f;
		vec_term2 = pcl.N1 * -pcl_n2_miu_div_k_vol_vx;
		fvec[8] += coef_mv * vec_term1 + coef_cv * vec_term2;
		vec_term1 = pcl.N2 * pcl_mvx_f;
		vec_term2 = pcl.N2 * -pcl_n2_miu_div_k_vol_vx;
		fvec[9] += coef_mv * vec_term1 + coef_cv * vec_term2;
		vec_term1 = pcl.N3 * pcl_mvx_f;
		vec_term2 = pcl.N3 * -pcl_n2_miu_div_k_vol_vx;
		fvec[10] += coef_mv * vec_term1 + coef_cv * vec_term2;
		vec_term1 = pcl.N4 * pcl_mvx_f;
		vec_term2 = pcl.N4 * -pcl_n2_miu_div_k_vol_vx;
		fvec[11] += coef_mv * vec_term1 + coef_cv * vec_term2;
		// fluid y
		vec_term1 = pcl.N1 * pcl_mvy_f;
		vec_term2 = pcl.N1 * -pcl_n2_miu_div_k_vol_vy;
		fvec[12] += coef_mv * vec_term1 + coef_cv * vec_term2;
		vec_term1 = pcl.N2 * pcl_mvy_f;
		vec_term2 = pcl.N2 * -pcl_n2_miu_div_k_vol_vy;
		fvec[13] += coef_mv * vec_term1 + coef_cv * vec_term2;
		vec_term1 = pcl.N3 * pcl_mvy_f;
		vec_term2 = pcl.N3 * -pcl_n2_miu_div_k_vol_vy;
		fvec[14] += coef_mv * vec_term1 + coef_cv * vec_term2;
		vec_term1 = pcl.N4 * pcl_mvy_f;
		vec_term2 = pcl.N4 * -pcl_n2_miu_div_k_vol_vy;
		fvec[15] += coef_mv * vec_term1 + coef_cv * vec_term2;
	}

	//print_mat(kmat, out_file);
}


// =============================================================================
// =============================== for debugging ===============================
// =============================================================================
void Step_S2D_CHM_s_uUp::elem_residual_force_vec(Element &e, double fvec[20])
{
	double pcl_m_f, pcl_max_s, pcl_may_s, pcl_max_f, pcl_may_f;
	double n2_miu_div_k_vol, spfx, spfy;
	memset(fvec, 0, sizeof(double) * 20);
	for (Particle *pcl_iter = e.pcls; pcl_iter; pcl_iter = pcl_iter->next)
	{
		Particle &pcl = *pcl_iter;
		Node &n1 = model->nodes[e.n1_id];
		Node &n2 = model->nodes[e.n2_id];
		Node &n3 = model->nodes[e.n3_id];
		Node &n4 = model->nodes[e.n4_id];
		
		// f_int
		double dux_s_dx = pcl.dN1_dx * n1.dux_s + pcl.dN2_dx * n2.dux_s
						+ pcl.dN3_dx * n3.dux_s + pcl.dN4_dx * n4.dux_s;
		double duy_s_dy = pcl.dN1_dy * n1.duy_s + pcl.dN2_dy * n2.duy_s
						+ pcl.dN3_dy * n3.duy_s + pcl.dN4_dy * n4.duy_s;
		double dux_f_dx = pcl.dN1_dx * n1.dux_f + pcl.dN2_dx * n2.dux_f
						+ pcl.dN3_dx * n3.dux_f + pcl.dN4_dx * n4.dux_f;
		double duy_f_dy = pcl.dN1_dy * n1.duy_f + pcl.dN2_dy * n2.duy_f
						+ pcl.dN3_dy * n3.duy_f + pcl.dN4_dy * n4.duy_f;
		double dp = pcl.N1 * n1.dp + pcl.N2 * n2.dp + pcl.N3 * n3.dp + pcl.N4 * n4.dp;
		//std::cout << dux_s_dx << " " << duy_s_dy << " " << dux_f_dx << " "
		//		  << duy_f_dy << " " << dp << "\n";
		//std::cout << pcl.N1 << " " << pcl.N2 << " " << pcl.N3 << " " << pcl.N4 << " " << pcl.n << " " << pcl.vol << "\n";
		add_internal_force_vec(fvec, pcl);
		fvec[16] += pcl.N1 * ((1.0 - pcl.n) * (dux_s_dx + duy_s_dy) + pcl.n * (dux_f_dx + duy_f_dy) + pcl.n/pcl.Kf * dp) * pcl.vol;
		fvec[17] += pcl.N2 * ((1.0 - pcl.n) * (dux_s_dx + duy_s_dy) + pcl.n * (dux_f_dx + duy_f_dy) + pcl.n/pcl.Kf * dp) * pcl.vol;
		fvec[18] += pcl.N3 * ((1.0 - pcl.n) * (dux_s_dx + duy_s_dy) + pcl.n * (dux_f_dx + duy_f_dy) + pcl.n/pcl.Kf * dp) * pcl.vol;
		fvec[19] += pcl.N4 * ((1.0 - pcl.n) * (dux_s_dx + duy_s_dy) + pcl.n * (dux_f_dx + duy_f_dy) + pcl.n/pcl.Kf * dp) * pcl.vol;

		// Ma
		pcl_max_s = (n1.ax_s * pcl.N1 + n2.ax_s * pcl.N2 + n3.ax_s * pcl.N3 + n4.ax_s * pcl.N4) * pcl.m_s;
		pcl_may_s = (n1.ay_s * pcl.N1 + n2.ay_s * pcl.N2 + n3.ay_s * pcl.N3 + n4.ay_s * pcl.N4) * pcl.m_s;
		pcl_m_f = pcl.n * pcl.density_f * pcl.vol;
		pcl_max_f = (n1.ax_f * pcl.N1 + n2.ax_f * pcl.N2 + n3.ax_f * pcl.N3 + n4.ax_f * pcl.N4) * pcl_m_f;
		pcl_may_f = (n1.ay_f * pcl.N1 + n2.ay_f * pcl.N2 + n3.ay_f * pcl.N3 + n4.ay_f * pcl.N4) * pcl_m_f;
		fvec[0]  -= pcl.N1 * pcl_max_s;
		fvec[1]  -= pcl.N2 * pcl_max_s;
		fvec[2]  -= pcl.N3 * pcl_max_s;
		fvec[3]  -= pcl.N4 * pcl_max_s;
		fvec[4]  -= pcl.N1 * pcl_may_s;
		fvec[5]  -= pcl.N2 * pcl_may_s;
		fvec[6]  -= pcl.N3 * pcl_may_s;
		fvec[7]  -= pcl.N4 * pcl_may_s;
		fvec[8]  -= pcl.N1 * pcl_max_f;
		fvec[9]  -= pcl.N2 * pcl_max_f;
		fvec[10] -= pcl.N3 * pcl_max_f;
		fvec[11] -= pcl.N4 * pcl_max_f;
		fvec[12] -= pcl.N1 * pcl_may_f;
		fvec[13] -= pcl.N2 * pcl_may_f;
		fvec[14] -= pcl.N3 * pcl_may_f;
		fvec[15] -= pcl.N4 * pcl_may_f;

		//Cv
		n2_miu_div_k_vol = pcl.n * pcl.n * pcl.miu / pcl.k * pcl.vol;
		spfx = ((n1.vx_s - n1.vx_f) * pcl.N1 + (n2.vx_s - n2.vx_f) * pcl.N2
			  + (n3.vx_s - n3.vx_f) * pcl.N3 + (n4.vx_s - n4.vx_f) * pcl.N4)
			   * n2_miu_div_k_vol;
		spfy = ((n1.vy_s - n1.vy_f) * pcl.N1 + (n2.vy_s - n2.vy_f) * pcl.N2
			  + (n3.vy_s - n3.vy_f) * pcl.N3 + (n4.vy_s - n4.vy_f) * pcl.N4)
			   * n2_miu_div_k_vol;
		fvec[0]  -= pcl.N1 * spfx;
		fvec[1]  -= pcl.N2 * spfx;
		fvec[2]  -= pcl.N3 * spfx;
		fvec[3]  -= pcl.N4 * spfx;
		fvec[4]  -= pcl.N1 * spfy;
		fvec[5]  -= pcl.N2 * spfy;
		fvec[6]  -= pcl.N3 * spfy;
		fvec[7]  -= pcl.N4 * spfy;
		fvec[8]  += pcl.N1 * spfx;
		fvec[9]  += pcl.N2 * spfx;
		fvec[10] += pcl.N3 * spfx;
		fvec[11] += pcl.N4 * spfx;
		fvec[12] += pcl.N1 * spfy;
		fvec[13] += pcl.N2 * spfy;
		fvec[14] += pcl.N3 * spfy;
		fvec[15] += pcl.N4 * spfy;
	}
}

void Step_S2D_CHM_s_uUp::form_global_residual_force(void)
{
	static size_t call_id = 0;
	++call_id;

	Model_S2D_CHM_s_uUp &md = *model;
	Eigen::VectorXd g_fvec(dof_num);
	g_fvec.setZero();

	size_t node_g_id, l2g_id_map[20];
	double e_fvec[20];
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
		if (e.pcls)
		{
			// form elemental stiffness matrix and force vector
			elem_residual_force_vec(e, e_fvec);

			// form map from global id to local id
			node_g_id = md.nodes[e.n1_id].g_id;
			l2g_id_map[0]  = n_id_to_dof_id(node_g_id, DOF::ux_s);
			l2g_id_map[4]  = n_id_to_dof_id(node_g_id, DOF::uy_s);
			l2g_id_map[8]  = n_id_to_dof_id(node_g_id, DOF::ux_f);
			l2g_id_map[12] = n_id_to_dof_id(node_g_id, DOF::uy_f);
			l2g_id_map[16] = n_id_to_dof_id(node_g_id, DOF::p);
			node_g_id = md.nodes[e.n2_id].g_id;
			l2g_id_map[1]  = n_id_to_dof_id(node_g_id, DOF::ux_s);
			l2g_id_map[5]  = n_id_to_dof_id(node_g_id, DOF::uy_s);
			l2g_id_map[9]  = n_id_to_dof_id(node_g_id, DOF::ux_f);
			l2g_id_map[13] = n_id_to_dof_id(node_g_id, DOF::uy_f);
			l2g_id_map[17] = n_id_to_dof_id(node_g_id, DOF::p);
			node_g_id = md.nodes[e.n3_id].g_id;
			l2g_id_map[2]  = n_id_to_dof_id(node_g_id, DOF::ux_s);
			l2g_id_map[6]  = n_id_to_dof_id(node_g_id, DOF::uy_s);
			l2g_id_map[10] = n_id_to_dof_id(node_g_id, DOF::ux_f);
			l2g_id_map[14] = n_id_to_dof_id(node_g_id, DOF::uy_f);
			l2g_id_map[18] = n_id_to_dof_id(node_g_id, DOF::p);
			node_g_id = md.nodes[e.n4_id].g_id;
			l2g_id_map[3]  = n_id_to_dof_id(node_g_id, DOF::ux_s);
			l2g_id_map[7]  = n_id_to_dof_id(node_g_id, DOF::uy_s);
			l2g_id_map[11] = n_id_to_dof_id(node_g_id, DOF::ux_f);
			l2g_id_map[15] = n_id_to_dof_id(node_g_id, DOF::uy_f);
			l2g_id_map[19] = n_id_to_dof_id(node_g_id, DOF::p);

			// add to global matrix and vector
			for (size_t l_id1 = 0; l_id1 < 20; ++l_id1)
			{
				size_t g_id1 = l2g_id_map[l_id1];
				g_fvec[g_id1] += e_fvec[l_id1];
			}
		}
	}

	// apply external force
	// traction
	size_t dof_g_id;
	for (size_t t_id = 0; t_id < md.tx_num; ++t_id)
	{
		TractionBC_MPM &tx = md.txs[t_id];
		Particle &pcl = md.pcls[tx.pcl_id];
		if (pcl.pe)
		{
			Element &elem = *pcl.pe;
			// node 1
			node_g_id = md.nodes[elem.n1_id].g_id;
			dof_g_id = n_id_to_dof_id(node_g_id, DOF::ux_s);
			g_fvec[dof_g_id] += pcl.N1 * tx.t;
			// node 2
			node_g_id = md.nodes[elem.n2_id].g_id;
			dof_g_id = n_id_to_dof_id(node_g_id, DOF::ux_s);
			g_fvec[dof_g_id] += pcl.N2 * tx.t;
			// node 3
			node_g_id = md.nodes[elem.n3_id].g_id;
			dof_g_id = n_id_to_dof_id(node_g_id, DOF::ux_s);
			g_fvec[dof_g_id] += pcl.N3 * tx.t;
			// node 4
			node_g_id = md.nodes[elem.n4_id].g_id;
			dof_g_id = n_id_to_dof_id(node_g_id, DOF::ux_s);
			g_fvec[dof_g_id] += pcl.N4 * tx.t;
		}
	}
	for (size_t t_id = 0; t_id < md.ty_num; ++t_id)
	{
		TractionBC_MPM &ty = md.tys[t_id];
		Particle &pcl = md.pcls[ty.pcl_id];
		if (pcl.pe)
		{
			Element &elem = *pcl.pe;
			// node 1
			node_g_id = md.nodes[elem.n1_id].g_id;
			dof_g_id = n_id_to_dof_id(node_g_id, DOF::uy_s);
			g_fvec[dof_g_id] += pcl.N1 * ty.t;
			// node 2
			node_g_id = md.nodes[elem.n2_id].g_id;
			dof_g_id = n_id_to_dof_id(node_g_id, DOF::uy_s);
			g_fvec[dof_g_id] += pcl.N2 * ty.t;
			// node 3
			node_g_id = md.nodes[elem.n3_id].g_id;
			dof_g_id = n_id_to_dof_id(node_g_id, DOF::uy_s);
			g_fvec[dof_g_id] += pcl.N3 * ty.t;
			// node 4
			node_g_id = md.nodes[elem.n4_id].g_id;
			dof_g_id = n_id_to_dof_id(node_g_id, DOF::uy_s);
			g_fvec[dof_g_id] += pcl.N4 * ty.t;
		}
	}
	// body force to be finished...

	//size_t node_num = md.node_num;
	//// Apply displacement boundary condition
	//for (size_t bc_id = 0; bc_id < md.usx_num; ++bc_id)
	//{
	//	DisplacementBC &dbc = md.usxs[bc_id];
	//	node_g_id = md.nodes[dbc.node_id].g_id;
	//	if (node_g_id != node_num)
	//	{
	//		dof_g_id = n_id_to_dof_id(node_g_id, DOF::ux_s);
	//		g_fvec[dof_g_id] = dbc.u;
	//	}
	//}
	//for (size_t bc_id = 0; bc_id < md.usy_num; ++bc_id)
	//{
	//	DisplacementBC &dbc = md.usys[bc_id];
	//	node_g_id = md.nodes[dbc.node_id].g_id;
	//	if (node_g_id != node_num)
	//	{
	//		dof_g_id = n_id_to_dof_id(node_g_id, DOF::uy_s);
	//		g_fvec[dof_g_id] = dbc.u;
	//	}
	//}
	//for (size_t bc_id = 0; bc_id < md.ufx_num; ++bc_id)
	//{
	//	DisplacementBC &dbc = md.ufxs[bc_id];
	//	node_g_id = md.nodes[dbc.node_id].g_id;
	//	if (node_g_id != node_num)
	//	{
	//		dof_g_id = n_id_to_dof_id(node_g_id, DOF::ux_f);
	//		g_fvec[dof_g_id] = dbc.u;
	//	}
	//}
	//for (size_t bc_id = 0; bc_id < md.ufy_num; ++bc_id)
	//{
	//	DisplacementBC &dbc = md.usys[bc_id];
	//	node_g_id = md.nodes[dbc.node_id].g_id;
	//	if (node_g_id != node_num)
	//	{
	//		dof_g_id = n_id_to_dof_id(node_g_id, DOF::uy_f);
	//		g_fvec[dof_g_id] = dbc.u;
	//	}
	//}

	std::cout << "Rf " << call_id << "\n" << g_fvec << "\n";

	std::cout << "ax_s: ";
	for (size_t i = 0; i < md.node_num; i++)
	{
		Node &n = md.nodes[i];
		std::cout << n.ax_s << ", ";
	}
	std::cout << "\n";

	std::cout << "ay_s: ";
	for (size_t i = 0; i < md.node_num; i++)
	{
		Node &n = md.nodes[i];
		std::cout << n.ay_s << ", ";
	}
	std::cout << "\n";

	std::cout << "ax_f: ";
	for (size_t i = 0; i < md.node_num; i++)
	{
		Node &n = md.nodes[i];
		std::cout << n.ax_f << ", ";
	}
	std::cout << "\n";

	std::cout << "ay_f: ";
	for (size_t i = 0; i < md.node_num; i++)
	{
		Node &n = md.nodes[i];
		std::cout << n.ay_f << ", ";
	}
	std::cout << "\n";

	std::cout << "vx_s: ";
	for (size_t i = 0; i < md.node_num; i++)
	{
		Node &n = md.nodes[i];
		std::cout << n.vx_s << ", ";
	}
	std::cout << "\n";

	std::cout << "vy_s: ";
	for (size_t i = 0; i < md.node_num; i++)
	{
		Node &n = md.nodes[i];
		std::cout << n.vy_s << ", ";
	}
	std::cout << "\n";

	std::cout << "vx_f: ";
	for (size_t i = 0; i < md.node_num; i++)
	{
		Node &n = md.nodes[i];
		std::cout << n.vx_f << ", ";
	}
	std::cout << "\n";

	std::cout << "vy_f: ";
	for (size_t i = 0; i < md.node_num; i++)
	{
		Node &n = md.nodes[i];
		std::cout << n.vy_f << ", ";
	}
	std::cout << "\n";

	std::cout << "p: ";
	for (size_t i = 0; i < md.node_num; i++)
	{
		Node &n = md.nodes[i];
		std::cout << n.p << ", ";
	}
	std::cout << "\n";
}