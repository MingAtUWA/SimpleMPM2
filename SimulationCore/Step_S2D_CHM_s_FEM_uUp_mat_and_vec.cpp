#include "SimulationCore_pcp.h"

#define Keep_Newmark_Coefficients
#include "Step_S2D_CHM_s_FEM_uUp.h"

namespace
{
	typedef Model_S2D_CHM_s_FEM_uUp::ShapeFuncValue ShapeFuncValue;
	typedef Model_S2D_CHM_s_FEM_uUp::GaussPoint GaussPoint_fem;
	typedef Model_S2D_CHM_s_FEM_uUp::Element Element_fem;
	typedef Model_S2D_CHM_s_FEM_uUp::Node Node_fem;
	typedef Model_S2D_CHM_s_FEM_uUp::DOF DOF;

	void print_mat(double mat[3][3])
	{
		for (size_t i = 0; i < 3; i++)
		{
			for (size_t j = 0; j < 3; j++)
				std::cout << mat[i][j] << ", ";
			std::cout << "\n";
		}
	}

	void print_mat(double mat[20][20],
				   std::fstream &out_file,
				   const char *mat_name = nullptr)
	{
		if (mat_name)
			out_file << mat_name << "\n";
		for (size_t i = 0; i < 20; i++)
		{
			for (size_t j = 0; j < 20; j++)
				out_file << mat[i][j] << ", ";
			out_file << "\n";
		}
	}
	
	void form_E_mat(GaussPoint_fem &gp, double E_mat[3][3])
	{
		double E_tmp = gp.E / (1.0 + gp.niu) / (1.0 - 2.0 * gp.niu);
		E_mat[0][0] = E_tmp * (1.0 - gp.niu);
		E_mat[0][1] = E_tmp * gp.niu;
		E_mat[0][2] = 0.0;
		E_mat[1][0] = E_tmp * gp.niu;
		E_mat[1][1] = E_tmp * (1.0 - gp.niu);
		E_mat[1][2] = 0.0;
		E_mat[2][0] = 0.0;
		E_mat[2][1] = 0.0;
		E_mat[2][2] = gp.E / (2.0 * (1.0 + gp.niu));
	}

	void cal_stiffness_mat(double kmat[20][20], 
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
				kmat[i][j] += (dN_dx[0][i] * mid_mat[0][j]
							 + dN_dx[1][i] * mid_mat[1][j]
							 + dN_dx[2][i] * mid_mat[2][j]) * weight;
	}
};

void Step_S2D_CHM_s_FEM_uUp::form_elem_stiffness_mat_and_force_vec(
	Element_fem &e, double kmat[20][20], double fvec[20])
{
	ShapeFuncValue &sf1 = model->gp1_sf;
	ShapeFuncValue &sf2 = model->gp2_sf;
	ShapeFuncValue &sf3 = model->gp3_sf;
	ShapeFuncValue &sf4 = model->gp4_sf;
	GaussPoint_fem &gp1 = e.gp1;
	GaussPoint_fem &gp2 = e.gp2;
	GaussPoint_fem &gp3 = e.gp3;
	GaussPoint_fem &gp4 = e.gp4;
	double gp_w = model->gp_w;
	double mat_term;
	double coef;

	memset(kmat, 0, sizeof(double) * 20 * 20);

	// mass matrix
	coef = 1.0 / (beta * dtime * dtime);
	// solid
	mat_term = (sf1.N1 * (1.0 - gp1.n) * gp1.density_s * sf1.N1
			  + sf2.N1 * (1.0 - gp2.n) * gp2.density_s * sf2.N1
			  + sf3.N1 * (1.0 - gp3.n) * gp3.density_s * sf3.N1
			  + sf4.N1 * (1.0 - gp4.n) * gp4.density_s * sf4.N1) * gp_w * coef;
	kmat[0][0] = mat_term;
	kmat[4][4] = mat_term;
	mat_term = (sf1.N1 * (1.0 - gp1.n) * gp1.density_s * sf1.N2
			  + sf2.N1 * (1.0 - gp2.n) * gp2.density_s * sf2.N2
			  + sf3.N1 * (1.0 - gp3.n) * gp3.density_s * sf3.N2
			  + sf4.N1 * (1.0 - gp4.n) * gp4.density_s * sf4.N2) * gp_w * coef;
	kmat[0][1] = mat_term;
	kmat[1][0] = mat_term;
	kmat[4][5] = mat_term;
	kmat[5][4] = mat_term;
	mat_term = (sf1.N1 * (1.0 - gp1.n) * gp1.density_s * sf1.N3
			  + sf2.N1 * (1.0 - gp2.n) * gp2.density_s * sf2.N3
			  + sf3.N1 * (1.0 - gp3.n) * gp3.density_s * sf3.N3
			  + sf4.N1 * (1.0 - gp4.n) * gp4.density_s * sf4.N3) * gp_w * coef;
	kmat[0][2] = mat_term;
	kmat[2][0] = mat_term;
	kmat[4][6] = mat_term;
	kmat[6][4] = mat_term;
	mat_term = (sf1.N1 * (1.0 - gp1.n) * gp1.density_s * sf1.N4
			  + sf2.N1 * (1.0 - gp2.n) * gp2.density_s * sf2.N4
			  + sf3.N1 * (1.0 - gp3.n) * gp3.density_s * sf3.N4
			  + sf4.N1 * (1.0 - gp4.n) * gp4.density_s * sf4.N4) * gp_w * coef;
	kmat[0][3] = mat_term;
	kmat[3][0] = mat_term;
	kmat[4][7] = mat_term;
	kmat[7][4] = mat_term;
	mat_term = (sf1.N2 * (1.0 - gp1.n) * gp1.density_s * sf1.N2
			  + sf2.N2 * (1.0 - gp2.n) * gp2.density_s * sf2.N2
			  + sf3.N2 * (1.0 - gp3.n) * gp3.density_s * sf3.N2
			  + sf4.N2 * (1.0 - gp4.n) * gp4.density_s * sf4.N2) * gp_w * coef;
	kmat[1][1] = mat_term;
	kmat[5][5] = mat_term;
	mat_term = (sf1.N2 * (1.0 - gp1.n) * gp1.density_s * sf1.N3
			  + sf2.N2 * (1.0 - gp2.n) * gp2.density_s * sf2.N3
			  + sf3.N2 * (1.0 - gp3.n) * gp3.density_s * sf3.N3
			  + sf4.N2 * (1.0 - gp4.n) * gp4.density_s * sf4.N3) * gp_w * coef;
	kmat[1][2] = mat_term;
	kmat[2][1] = mat_term;
	kmat[5][6] = mat_term;
	kmat[6][5] = mat_term;
	mat_term = (sf1.N2 * (1.0 - gp1.n) * gp1.density_s * sf1.N4
			  + sf2.N2 * (1.0 - gp2.n) * gp2.density_s * sf2.N4
			  + sf3.N2 * (1.0 - gp3.n) * gp3.density_s * sf3.N4
			  + sf4.N2 * (1.0 - gp4.n) * gp4.density_s * sf4.N4) * gp_w * coef;
	kmat[1][3] = mat_term;
	kmat[3][1] = mat_term;
	kmat[5][7] = mat_term;
	kmat[7][5] = mat_term;
	mat_term = (sf1.N3 * (1.0 - gp1.n) * gp1.density_s * sf1.N3
			  + sf2.N3 * (1.0 - gp2.n) * gp2.density_s * sf2.N3
			  + sf3.N3 * (1.0 - gp3.n) * gp3.density_s * sf3.N3
			  + sf4.N3 * (1.0 - gp4.n) * gp4.density_s * sf4.N3) * gp_w * coef;
	kmat[2][2] = mat_term;
	kmat[6][6] = mat_term;
	mat_term = (sf1.N3 * (1.0 - gp1.n) * gp1.density_s * sf1.N4
			  + sf2.N3 * (1.0 - gp2.n) * gp2.density_s * sf2.N4
			  + sf3.N3 * (1.0 - gp3.n) * gp3.density_s * sf3.N4
			  + sf4.N3 * (1.0 - gp4.n) * gp4.density_s * sf4.N4) * gp_w * coef;
	kmat[2][3] = mat_term;
	kmat[3][2] = mat_term;
	kmat[6][7] = mat_term;
	kmat[7][6] = mat_term;
	mat_term = (sf1.N4 * (1.0 - gp1.n) * gp1.density_s * sf1.N4
			  + sf2.N4 * (1.0 - gp2.n) * gp2.density_s * sf2.N4
			  + sf3.N4 * (1.0 - gp3.n) * gp3.density_s * sf3.N4
			  + sf4.N4 * (1.0 - gp4.n) * gp4.density_s * sf4.N4) * gp_w * coef;
	kmat[3][3] = mat_term;
	kmat[7][7] = mat_term;
	// fluid
	mat_term = (sf1.N1 * gp1.n * gp1.density_f * sf1.N1
			  + sf2.N1 * gp2.n * gp2.density_f * sf2.N1
			  + sf3.N1 * gp3.n * gp3.density_f * sf3.N1
			  + sf4.N1 * gp4.n * gp4.density_f * sf4.N1) * gp_w * coef;
	kmat[8][8] = mat_term;
	kmat[12][12] = mat_term;
	mat_term = (sf1.N1 * gp1.n * gp1.density_f * sf1.N2
			  + sf2.N1 * gp2.n * gp2.density_f * sf2.N2
			  + sf3.N1 * gp3.n * gp3.density_f * sf3.N2
			  + sf4.N1 * gp4.n * gp4.density_f * sf4.N2) * gp_w * coef;
	kmat[8][9] = mat_term;
	kmat[9][8] = mat_term;
	kmat[12][13] = mat_term;
	kmat[13][12] = mat_term;
	mat_term = (sf1.N1 * gp1.n * gp1.density_f * sf1.N3
			  + sf2.N1 * gp2.n * gp2.density_f * sf2.N3
			  + sf3.N1 * gp3.n * gp3.density_f * sf3.N3
			  + sf4.N1 * gp4.n * gp4.density_f * sf4.N3) * gp_w * coef;
	kmat[8][10] = mat_term;
	kmat[10][8] = mat_term;
	kmat[12][14] = mat_term;
	kmat[14][12] = mat_term;
	mat_term = (sf1.N1 * gp1.n * gp1.density_f * sf1.N4
			  + sf2.N1 * gp2.n * gp2.density_f * sf2.N4
			  + sf3.N1 * gp3.n * gp3.density_f * sf3.N4
			  + sf4.N1 * gp4.n * gp4.density_f * sf4.N4) * gp_w * coef;
	kmat[8][11] = mat_term;
	kmat[11][8] = mat_term;
	kmat[12][15] = mat_term;
	kmat[15][12] = mat_term;
	mat_term = (sf1.N2 * gp1.n * gp1.density_f * sf1.N2
			  + sf2.N2 * gp2.n * gp2.density_f * sf2.N2
			  + sf3.N2 * gp3.n * gp3.density_f * sf3.N2
			  + sf4.N2 * gp4.n * gp4.density_f * sf4.N2) * gp_w * coef;
	kmat[9][9] = mat_term;
	kmat[13][13] = mat_term;
	mat_term = (sf1.N2 * gp1.n * gp1.density_f * sf1.N3
			  + sf2.N2 * gp2.n * gp2.density_f * sf2.N3
			  + sf3.N2 * gp3.n * gp3.density_f * sf3.N3
			  + sf4.N2 * gp4.n * gp4.density_f * sf4.N3) * gp_w * coef;
	kmat[9][10] = mat_term;
	kmat[10][9] = mat_term;
	kmat[13][14] = mat_term;
	kmat[14][13] = mat_term;
	mat_term = (sf1.N2 * gp1.n * gp1.density_f * sf1.N4
			  + sf2.N2 * gp2.n * gp2.density_f * sf2.N4
			  + sf3.N2 * gp3.n * gp3.density_f * sf3.N4
			  + sf4.N2 * gp4.n * gp4.density_f * sf4.N4) * gp_w * coef;
	kmat[9][11] = mat_term;
	kmat[11][9] = mat_term;
	kmat[13][15] = mat_term;
	kmat[15][13] = mat_term;
	mat_term = (sf1.N3 * gp1.n * gp1.density_f * sf1.N3
			  + sf2.N3 * gp2.n * gp2.density_f * sf2.N3
			  + sf3.N3 * gp3.n * gp3.density_f * sf3.N3
			  + sf4.N3 * gp4.n * gp4.density_f * sf4.N3) * gp_w * coef;
	kmat[10][10] = mat_term;
	kmat[14][14] = mat_term;
	mat_term = (sf1.N3 * gp1.n * gp1.density_f * sf1.N4
			  + sf2.N3 * gp2.n * gp2.density_f * sf2.N4
			  + sf3.N3 * gp3.n * gp3.density_f * sf3.N4
			  + sf4.N3 * gp4.n * gp4.density_f * sf4.N4) * gp_w * coef;
	kmat[10][11] = mat_term;
	kmat[11][10] = mat_term;
	kmat[14][15] = mat_term;
	kmat[15][14] = mat_term;
	mat_term = (sf1.N4 * gp1.n * gp1.density_f * sf1.N4
			  + sf2.N4 * gp2.n * gp2.density_f * sf2.N4
			  + sf3.N4 * gp3.n * gp3.density_f * sf3.N4
			  + sf4.N4 * gp4.n * gp4.density_f * sf4.N4) * gp_w * coef;
	kmat[11][11] = mat_term;
	kmat[15][15] = mat_term;

	// damping matrix
	coef = gamma / (beta * dtime);
	double sp_f1 = gp1.n * gp1.n * gp1.miu / gp1.k;
	double sp_f2 = gp2.n * gp2.n * gp2.miu / gp2.k;
	double sp_f3 = gp3.n * gp3.n * gp3.miu / gp3.k;
	double sp_f4 = gp4.n * gp4.n * gp4.miu / gp4.k;
#define FILL_SEEP_MAT(x_id, y_id) \
	kmat[x_id][y_id] += mat_term; \
	kmat[x_id][y_id+8] -= mat_term; \
	kmat[x_id+4][y_id+4] += mat_term; \
	kmat[x_id+4][y_id+12] -= mat_term; \
	kmat[x_id+8][y_id] -= mat_term; \
	kmat[x_id+8][y_id+8] += mat_term; \
	kmat[x_id+12][y_id+4] -= mat_term; \
	kmat[x_id+12][y_id+12] += mat_term
	mat_term = (sf1.N1 * sp_f1 * sf1.N1 + sf2.N1 * sp_f2 * sf2.N1
			  + sf3.N1 * sp_f3 * sf3.N1 + sf4.N1 * sp_f4 * sf4.N1) * gp_w * coef;
	FILL_SEEP_MAT(0, 0);
	mat_term = (sf1.N1 * sp_f1 * sf1.N2 + sf2.N1 * sp_f2 * sf2.N2
			   + sf3.N1 * sp_f3 * sf3.N2 + sf4.N1 * sp_f4 * sf4.N2) * gp_w * coef;
	FILL_SEEP_MAT(0, 1);
	mat_term = (sf1.N1 * sp_f1 * sf1.N3 + sf2.N1 * sp_f2 * sf2.N3
			  + sf3.N1 * sp_f3 * sf3.N3 + sf4.N1 * sp_f4 * sf4.N3) * gp_w * coef;
	FILL_SEEP_MAT(0, 2);
	mat_term = (sf1.N1 * sp_f1 * sf1.N4 + sf2.N1 * sp_f2 * sf2.N4
			  + sf3.N1 * sp_f3 * sf3.N4 + sf4.N1 * sp_f4 * sf4.N4) * gp_w * coef;
	FILL_SEEP_MAT(0, 3);
	mat_term = (sf1.N2 * sp_f1 * sf1.N1 + sf2.N2 * sp_f2 * sf2.N1
			  + sf3.N2 * sp_f3 * sf3.N1 + sf4.N2 * sp_f4 * sf4.N1) * gp_w * coef;
	FILL_SEEP_MAT(1, 0);
	mat_term = (sf1.N2 * sp_f1 * sf1.N2 + sf2.N2 * sp_f2 * sf2.N2
			  + sf3.N2 * sp_f3 * sf3.N2 + sf4.N2 * sp_f4 * sf4.N2) * gp_w * coef;
	FILL_SEEP_MAT(1, 1);
	mat_term = (sf1.N2 * sp_f1 * sf1.N3 + sf2.N2 * sp_f2 * sf2.N3
			  + sf3.N2 * sp_f3 * sf3.N3 + sf4.N2 * sp_f4 * sf4.N3) * gp_w * coef;
	FILL_SEEP_MAT(1, 2);
	mat_term = (sf1.N2 * sp_f1 * sf1.N4 + sf2.N2 * sp_f2 * sf2.N4
			  + sf3.N2 * sp_f3 * sf3.N4 + sf4.N2 * sp_f4 * sf4.N4) * gp_w * coef;
	FILL_SEEP_MAT(1, 3);
	mat_term = (sf1.N3 * sp_f1 * sf1.N1 + sf2.N3 * sp_f2 * sf2.N1
			  + sf3.N3 * sp_f3 * sf3.N1 + sf4.N3 * sp_f4 * sf4.N1) * gp_w * coef;
	FILL_SEEP_MAT(2, 0);
	mat_term = (sf1.N3 * sp_f1 * sf1.N2 + sf2.N3 * sp_f2 * sf2.N2
			  + sf3.N3 * sp_f3 * sf3.N2 + sf4.N3 * sp_f4 * sf4.N2) * gp_w * coef;
	FILL_SEEP_MAT(2, 1);
	mat_term = (sf1.N3 * sp_f1 * sf1.N3 + sf2.N3 * sp_f2 * sf2.N3
			  + sf3.N3 * sp_f3 * sf3.N3 + sf4.N3 * sp_f4 * sf4.N3) * gp_w * coef;
	FILL_SEEP_MAT(2, 2);
	mat_term = (sf1.N3 * sp_f1 * sf1.N4 + sf2.N3 * sp_f2 * sf2.N4
			 + sf3.N3 * sp_f3 * sf3.N4 + sf4.N3 * sp_f4 * sf4.N4) * gp_w * coef;
	FILL_SEEP_MAT(2, 3);
	mat_term = (sf1.N4 * sp_f1 * sf1.N1 + sf2.N4 * sp_f2 * sf2.N1
			  + sf3.N4 * sp_f3 * sf3.N1 + sf4.N4 * sp_f4 * sf4.N1) * gp_w * coef;
	FILL_SEEP_MAT(3, 0);
	mat_term = (sf1.N4 * sp_f1 * sf1.N2 + sf2.N4 * sp_f2 * sf2.N2
			  + sf3.N4 * sp_f3 * sf3.N2 + sf4.N4 * sp_f4 * sf4.N2) * gp_w * coef;
	FILL_SEEP_MAT(3, 1);
	mat_term = (sf1.N4 * sp_f1 * sf1.N3 + sf2.N4 * sp_f2 * sf2.N3
			  + sf3.N4 * sp_f3 * sf3.N3 + sf4.N4 * sp_f4 * sf4.N3) * gp_w * coef;
	FILL_SEEP_MAT(3, 2);
	mat_term = (sf1.N4 * sp_f1 * sf1.N4 + sf2.N4 * sp_f2 * sf2.N4
			  + sf3.N4 * sp_f3 * sf3.N4 + sf4.N4 * sp_f4 * sf4.N4) * gp_w * coef;
	FILL_SEEP_MAT(3, 3);

	// stiffness matrix
	double E_mat[3][3];
	// K
	form_E_mat(gp1, E_mat);
	cal_stiffness_mat(kmat, E_mat, model->dN_dx_mat1, gp_w);
	form_E_mat(gp2, E_mat);
	cal_stiffness_mat(kmat, E_mat, model->dN_dx_mat2, gp_w);
	form_E_mat(gp3, E_mat);
	cal_stiffness_mat(kmat, E_mat, model->dN_dx_mat3, gp_w);
	form_E_mat(gp4, E_mat);
	cal_stiffness_mat(kmat, E_mat, model->dN_dx_mat4, gp_w);
	// Gsx
	mat_term = (sf1.dN1_dx * (1.0 - gp1.n) * sf1.N1 + sf2.dN1_dx * (1.0 - gp2.n) * sf2.N1
			  + sf3.dN1_dx * (1.0 - gp3.n) * sf3.N1 + sf4.dN1_dx * (1.0 - gp4.n) * sf4.N1) * -gp_w;
	kmat[0][16] += mat_term;
	kmat[16][0] += mat_term;
	mat_term = (sf1.dN1_dx * (1.0 - gp1.n) * sf1.N2 + sf2.dN1_dx * (1.0 - gp2.n) * sf2.N2
			  + sf3.dN1_dx * (1.0 - gp3.n) * sf3.N2 + sf4.dN1_dx * (1.0 - gp4.n) * sf4.N2) * -gp_w;
	kmat[0][17] += mat_term;
	kmat[17][0] += mat_term;
	mat_term = (sf1.dN1_dx * (1.0 - gp1.n) * sf1.N3 + sf2.dN1_dx * (1.0 - gp2.n) * sf2.N3
			  + sf3.dN1_dx * (1.0 - gp3.n) * sf3.N3 + sf4.dN1_dx * (1.0 - gp4.n) * sf4.N3) * -gp_w;
	kmat[0][18] += mat_term;
	kmat[18][0] += mat_term;
	mat_term = (sf1.dN1_dx * (1.0 - gp1.n) * sf1.N4 + sf2.dN1_dx * (1.0 - gp2.n) * sf2.N4
			  + sf3.dN1_dx * (1.0 - gp3.n) * sf3.N4 + sf4.dN1_dx * (1.0 - gp4.n) * sf4.N4) * -gp_w;
	kmat[0][19] += mat_term;
	kmat[19][0] += mat_term;
	mat_term = (sf1.dN2_dx * (1.0 - gp1.n) * sf1.N1 + sf2.dN2_dx * (1.0 - gp2.n) * sf2.N1
			  + sf3.dN2_dx * (1.0 - gp3.n) * sf3.N1 + sf4.dN2_dx * (1.0 - gp4.n) * sf4.N1) * -gp_w;
	kmat[1][16] += mat_term;
	kmat[16][1] += mat_term;
	mat_term = (sf1.dN2_dx * (1.0 - gp1.n) * sf1.N2 + sf2.dN2_dx * (1.0 - gp2.n) * sf2.N2
			  + sf3.dN2_dx * (1.0 - gp3.n) * sf3.N2 + sf4.dN2_dx * (1.0 - gp4.n) * sf4.N2) * -gp_w;
	kmat[1][17] += mat_term;
	kmat[17][1] += mat_term;
	mat_term = (sf1.dN2_dx * (1.0 - gp1.n) * sf1.N3 + sf2.dN2_dx * (1.0 - gp2.n) * sf2.N3
			  + sf3.dN2_dx * (1.0 - gp3.n) * sf3.N3 + sf4.dN2_dx * (1.0 - gp4.n) * sf4.N3) * -gp_w;
	kmat[1][18] += mat_term;
	kmat[18][1] += mat_term;
	mat_term = (sf1.dN2_dx * (1.0 - gp1.n) * sf1.N4 + sf2.dN2_dx * (1.0 - gp2.n) * sf2.N4
			  + sf3.dN2_dx * (1.0 - gp3.n) * sf3.N4 + sf4.dN2_dx * (1.0 - gp4.n) * sf4.N4) * -gp_w;
	kmat[1][19] += mat_term;
	kmat[19][1] += mat_term;
	mat_term = (sf1.dN3_dx * (1.0 - gp1.n) * sf1.N1 + sf2.dN3_dx * (1.0 - gp2.n) * sf2.N1
			  + sf3.dN3_dx * (1.0 - gp3.n) * sf3.N1 + sf4.dN3_dx * (1.0 - gp4.n) * sf4.N1) * -gp_w;
	kmat[2][16] += mat_term;
	kmat[16][2] += mat_term;
	mat_term = (sf1.dN3_dx * (1.0 - gp1.n) * sf1.N2 + sf2.dN3_dx * (1.0 - gp2.n) * sf2.N2
			  + sf3.dN3_dx * (1.0 - gp3.n) * sf3.N2 + sf4.dN3_dx * (1.0 - gp4.n) * sf4.N2) * -gp_w;
	kmat[2][17] += mat_term;
	kmat[17][2] += mat_term;
	mat_term = (sf1.dN3_dx * (1.0 - gp1.n) * sf1.N3 + sf2.dN3_dx * (1.0 - gp2.n) * sf2.N3
			  + sf3.dN3_dx * (1.0 - gp3.n) * sf3.N3 + sf4.dN3_dx * (1.0 - gp4.n) * sf4.N3) * -gp_w;
	kmat[2][18] += mat_term;
	kmat[18][2] += mat_term;
	mat_term = (sf1.dN3_dx * (1.0 - gp1.n) * sf1.N4 + sf2.dN3_dx * (1.0 - gp2.n) * sf2.N4
			  + sf3.dN3_dx * (1.0 - gp3.n) * sf3.N4 + sf4.dN3_dx * (1.0 - gp4.n) * sf4.N4) * -gp_w;
	kmat[2][19] += mat_term;
	kmat[19][2] += mat_term;
	mat_term = (sf1.dN4_dx * (1.0 - gp1.n) * sf1.N1 + sf2.dN4_dx * (1.0 - gp2.n) * sf2.N1
			  + sf3.dN4_dx * (1.0 - gp3.n) * sf3.N1 + sf4.dN4_dx * (1.0 - gp4.n) * sf4.N1) * -gp_w;
	kmat[3][16] += mat_term;
	kmat[16][3] += mat_term;
	mat_term = (sf1.dN4_dx * (1.0 - gp1.n) * sf1.N2 + sf2.dN4_dx * (1.0 - gp2.n) * sf2.N2
			  + sf3.dN4_dx * (1.0 - gp3.n) * sf3.N2 + sf4.dN4_dx * (1.0 - gp4.n) * sf4.N2) * -gp_w;
	kmat[3][17] += mat_term;
	kmat[17][3] += mat_term;
	mat_term = (sf1.dN4_dx * (1.0 - gp1.n) * sf1.N3 + sf2.dN4_dx * (1.0 - gp2.n) * sf2.N3
			  + sf3.dN4_dx * (1.0 - gp3.n) * sf3.N3 + sf4.dN4_dx * (1.0 - gp4.n) * sf4.N3) * -gp_w;
	kmat[3][18] += mat_term;
	kmat[18][3] += mat_term;
	mat_term = (sf1.dN4_dx * (1.0 - gp1.n) * sf1.N4 + sf2.dN4_dx * (1.0 - gp2.n) * sf2.N4
			  + sf3.dN4_dx * (1.0 - gp3.n) * sf3.N4 + sf4.dN4_dx * (1.0 - gp4.n) * sf4.N4) * -gp_w;
	kmat[3][19] += mat_term;
	kmat[19][3] += mat_term;
	// Gsy
	mat_term = (sf1.dN1_dy * (1.0 - gp1.n) * sf1.N1 + sf2.dN1_dy * (1.0 - gp2.n) * sf2.N1
			  + sf3.dN1_dy * (1.0 - gp3.n) * sf3.N1 + sf4.dN1_dy * (1.0 - gp4.n) * sf4.N1) * -gp_w;
	kmat[4][16] += mat_term;
	kmat[16][4] += mat_term;
	mat_term = (sf1.dN1_dy * (1.0 - gp1.n) * sf1.N2 + sf2.dN1_dy * (1.0 - gp2.n) * sf2.N2
			  + sf3.dN1_dy * (1.0 - gp3.n) * sf3.N2 + sf4.dN1_dy * (1.0 - gp4.n) * sf4.N2) * -gp_w;
	kmat[4][17] += mat_term;
	kmat[17][4] += mat_term;
	mat_term = (sf1.dN1_dy * (1.0 - gp1.n) * sf1.N3 + sf2.dN1_dy * (1.0 - gp2.n) * sf2.N3
			  + sf3.dN1_dy * (1.0 - gp3.n) * sf3.N3 + sf4.dN1_dy * (1.0 - gp4.n) * sf4.N3) * -gp_w;
	kmat[4][18] += mat_term;
	kmat[18][4] += mat_term;
	mat_term = (sf1.dN1_dy * (1.0 - gp1.n) * sf1.N4 + sf2.dN1_dy * (1.0 - gp2.n) * sf2.N4
			  + sf3.dN1_dy * (1.0 - gp3.n) * sf3.N4 + sf4.dN1_dy * (1.0 - gp4.n) * sf4.N4) * -gp_w;
	kmat[4][19] += mat_term;
	kmat[19][4] += mat_term;
	mat_term = (sf1.dN2_dy * (1.0 - gp1.n) * sf1.N1 + sf2.dN2_dy * (1.0 - gp2.n) * sf2.N1
			  + sf3.dN2_dy * (1.0 - gp3.n) * sf3.N1 + sf4.dN2_dy * (1.0 - gp4.n) * sf4.N1) * -gp_w;
	kmat[5][16] += mat_term;
	kmat[16][5] += mat_term;
	mat_term = (sf1.dN2_dy * (1.0 - gp1.n) * sf1.N2 + sf2.dN2_dy * (1.0 - gp2.n) * sf2.N2
			  + sf3.dN2_dy * (1.0 - gp3.n) * sf3.N2 + sf4.dN2_dy * (1.0 - gp4.n) * sf4.N2) * -gp_w;
	kmat[5][17] += mat_term;
	kmat[17][5] += mat_term;
	mat_term = (sf1.dN2_dy * (1.0 - gp1.n) * sf1.N3 + sf2.dN2_dy * (1.0 - gp2.n) * sf2.N3
			  + sf3.dN2_dy * (1.0 - gp3.n) * sf3.N3 + sf4.dN2_dy * (1.0 - gp4.n) * sf4.N3) * -gp_w;
	kmat[5][18] += mat_term;
	kmat[18][5] += mat_term;
	mat_term = (sf1.dN2_dy * (1.0 - gp1.n) * sf1.N4 + sf2.dN2_dy * (1.0 - gp2.n) * sf2.N4
			  + sf3.dN2_dy * (1.0 - gp3.n) * sf3.N4 + sf4.dN2_dy * (1.0 - gp4.n) * sf4.N4) * -gp_w;
	kmat[5][19] += mat_term;
	kmat[19][5] += mat_term;
	mat_term = (sf1.dN3_dy * (1.0 - gp1.n) * sf1.N1 + sf2.dN3_dy * (1.0 - gp2.n) * sf2.N1
			  + sf3.dN3_dy * (1.0 - gp3.n) * sf3.N1 + sf4.dN3_dy * (1.0 - gp4.n) * sf4.N1) * -gp_w;
	kmat[6][16] += mat_term;
	kmat[16][6] += mat_term;
	mat_term = (sf1.dN3_dy * (1.0 - gp1.n) * sf1.N2 + sf2.dN3_dy * (1.0 - gp2.n) * sf2.N2
			  + sf3.dN3_dy * (1.0 - gp3.n) * sf3.N2 + sf4.dN3_dy * (1.0 - gp4.n) * sf4.N2) * -gp_w;
	kmat[6][17] += mat_term;
	kmat[17][6] += mat_term;
	mat_term = (sf1.dN3_dy * (1.0 - gp1.n) * sf1.N3 + sf2.dN3_dy * (1.0 - gp2.n) * sf2.N3
			  + sf3.dN3_dy * (1.0 - gp3.n) * sf3.N3 + sf4.dN3_dy * (1.0 - gp4.n) * sf4.N3) * -gp_w;
	kmat[6][18] += mat_term;
	kmat[18][6] += mat_term;
	mat_term = (sf1.dN3_dy * (1.0 - gp1.n) * sf1.N4 + sf2.dN3_dy * (1.0 - gp2.n) * sf2.N4
			  + sf3.dN3_dy * (1.0 - gp3.n) * sf3.N4 + sf4.dN3_dy * (1.0 - gp4.n) * sf4.N4) * -gp_w;
	kmat[6][19] += mat_term;
	kmat[19][6] += mat_term;
	mat_term = (sf1.dN4_dy * (1.0 - gp1.n) * sf1.N1 + sf2.dN4_dy * (1.0 - gp2.n) * sf2.N1
			  + sf3.dN4_dy * (1.0 - gp3.n) * sf3.N1 + sf4.dN4_dy * (1.0 - gp4.n) * sf4.N1) * -gp_w;
	kmat[7][16] += mat_term;
	kmat[16][7] += mat_term;
	mat_term = (sf1.dN4_dy * (1.0 - gp1.n) * sf1.N2 + sf2.dN4_dy * (1.0 - gp2.n) * sf2.N2
			  + sf3.dN4_dy * (1.0 - gp3.n) * sf3.N2 + sf4.dN4_dy * (1.0 - gp4.n) * sf4.N2) * -gp_w;
	kmat[7][17] += mat_term;
	kmat[17][7] += mat_term;
	mat_term = (sf1.dN4_dy * (1.0 - gp1.n) * sf1.N3 + sf2.dN4_dy * (1.0 - gp2.n) * sf2.N3
			  + sf3.dN4_dy * (1.0 - gp3.n) * sf3.N3 + sf4.dN4_dy * (1.0 - gp4.n) * sf4.N3) * -gp_w;
	kmat[7][18] += mat_term;
	kmat[18][7] += mat_term;
	mat_term = (sf1.dN4_dy * (1.0 - gp1.n) * sf1.N4 + sf2.dN4_dy * (1.0 - gp2.n) * sf2.N4
			  + sf3.dN4_dy * (1.0 - gp3.n) * sf3.N4 + sf4.dN4_dy * (1.0 - gp4.n) * sf4.N4) * -gp_w;
	kmat[7][19] += mat_term;
	kmat[19][7] += mat_term;
	// Gfx
	mat_term = (sf1.dN1_dx * gp1.n * sf1.N1 + sf2.dN1_dx * gp2.n * sf2.N1
			  + sf3.dN1_dx * gp3.n * sf3.N1 + sf4.dN1_dx * gp4.n * sf4.N1) * -gp_w;
	kmat[8][16] += mat_term;
	kmat[16][8] += mat_term;
	mat_term = (sf1.dN1_dx * gp1.n * sf1.N2 + sf2.dN1_dx * gp2.n * sf2.N2
			  + sf3.dN1_dx * gp3.n * sf3.N2 + sf4.dN1_dx * gp4.n * sf4.N2) * -gp_w;
	kmat[8][17] += mat_term;
	kmat[17][8] += mat_term;
	mat_term = (sf1.dN1_dx * gp1.n * sf1.N3 + sf2.dN1_dx * gp2.n * sf2.N3
			  + sf3.dN1_dx * gp3.n * sf3.N3 + sf4.dN1_dx * gp4.n * sf4.N3) * -gp_w;
	kmat[8][18] += mat_term;
	kmat[18][8] += mat_term;
	mat_term = (sf1.dN1_dx * gp1.n * sf1.N4 + sf2.dN1_dx * gp2.n * sf2.N4
			  + sf3.dN1_dx * gp3.n * sf3.N4 + sf4.dN1_dx * gp4.n * sf4.N4) * -gp_w;
	kmat[8][19] += mat_term;
	kmat[19][8] += mat_term;
	mat_term = (sf1.dN2_dx * gp1.n * sf1.N1 + sf2.dN2_dx * gp2.n * sf2.N1
			  + sf3.dN2_dx * gp3.n * sf3.N1 + sf4.dN2_dx * gp4.n * sf4.N1) * -gp_w;
	kmat[9][16] += mat_term;
	kmat[16][9] += mat_term;
	mat_term = (sf1.dN2_dx * gp1.n * sf1.N2 + sf2.dN2_dx * gp2.n * sf2.N2
			  + sf3.dN2_dx * gp3.n * sf3.N2 + sf4.dN2_dx * gp4.n * sf4.N2) * -gp_w;
	kmat[9][17] += mat_term;
	kmat[17][9] += mat_term;
	mat_term = (sf1.dN2_dx * gp1.n * sf1.N3 + sf2.dN2_dx * gp2.n * sf2.N3
			  + sf3.dN2_dx * gp3.n * sf3.N3 + sf4.dN2_dx * gp4.n * sf4.N3) * -gp_w;
	kmat[9][18] += mat_term;
	kmat[18][9] += mat_term;
	mat_term = (sf1.dN2_dx * gp1.n * sf1.N4 + sf2.dN2_dx * gp2.n * sf2.N4
			  + sf3.dN2_dx * gp3.n * sf3.N4 + sf4.dN2_dx * gp4.n * sf4.N4) * -gp_w;
	kmat[9][19] += mat_term;
	kmat[19][9] += mat_term;
	mat_term = (sf1.dN3_dx * gp1.n * sf1.N1 + sf2.dN3_dx * gp2.n * sf2.N1
			  + sf3.dN3_dx * gp3.n * sf3.N1 + sf4.dN3_dx * gp4.n * sf4.N1) * -gp_w;
	kmat[10][16] += mat_term;
	kmat[16][10] += mat_term;
	mat_term = (sf1.dN3_dx * gp1.n * sf1.N2 + sf2.dN3_dx * gp2.n * sf2.N2
			  + sf3.dN3_dx * gp3.n * sf3.N2 + sf4.dN3_dx * gp4.n * sf4.N2) * -gp_w;
	kmat[10][17] += mat_term;
	kmat[17][10] += mat_term;
	mat_term = (sf1.dN3_dx * gp1.n * sf1.N3 + sf2.dN3_dx * gp2.n * sf2.N3
			  + sf3.dN3_dx * gp3.n * sf3.N3 + sf4.dN3_dx * gp4.n * sf4.N3) * -gp_w;
	kmat[10][18] += mat_term;
	kmat[18][10] += mat_term;
	mat_term = (sf1.dN3_dx * gp1.n * sf1.N4 + sf2.dN3_dx * gp2.n * sf2.N4
			  + sf3.dN3_dx * gp3.n * sf3.N4 + sf4.dN3_dx * gp4.n * sf4.N4) * -gp_w;
	kmat[10][19] += mat_term;
	kmat[19][10] += mat_term;
	mat_term = (sf1.dN4_dx * gp1.n * sf1.N1 + sf2.dN4_dx * gp2.n * sf2.N1
			  + sf3.dN4_dx * gp3.n * sf3.N1 + sf4.dN4_dx * gp4.n * sf4.N1) * -gp_w;
	kmat[11][16] += mat_term;
	kmat[16][11] += mat_term;
	mat_term = (sf1.dN4_dx * gp1.n * sf1.N2 + sf2.dN4_dx * gp2.n * sf2.N2
			  + sf3.dN4_dx * gp3.n * sf3.N2 + sf4.dN4_dx * gp4.n * sf4.N2) * -gp_w;
	kmat[11][17] += mat_term;
	kmat[17][11] += mat_term;
	mat_term = (sf1.dN4_dx * gp1.n * sf1.N3 + sf2.dN4_dx * gp2.n * sf2.N3
			  + sf3.dN4_dx * gp3.n * sf3.N3 + sf4.dN4_dx * gp4.n * sf4.N3) * -gp_w;
	kmat[11][18] += mat_term;
	kmat[18][11] += mat_term;
	mat_term = (sf1.dN4_dx * gp1.n * sf1.N4 + sf2.dN4_dx * gp2.n * sf2.N4
			  + sf3.dN4_dx * gp3.n * sf3.N4 + sf4.dN4_dx * gp4.n * sf4.N4) * -gp_w;
	kmat[11][19] += mat_term;
	kmat[19][11] += mat_term;
	// Gfy
	mat_term = (sf1.dN1_dy * gp1.n * sf1.N1 + sf2.dN1_dy * gp2.n * sf2.N1
			  + sf3.dN1_dy * gp3.n * sf3.N1 + sf4.dN1_dy * gp4.n * sf4.N1) * -gp_w;
	kmat[12][16] += mat_term;
	kmat[16][12] += mat_term;
	mat_term = (sf1.dN1_dy * gp1.n * sf1.N2 + sf2.dN1_dy * gp2.n * sf2.N2
			  + sf3.dN1_dy * gp3.n * sf3.N2 + sf4.dN1_dy * gp4.n * sf4.N2) * -gp_w;
	kmat[12][17] += mat_term;
	kmat[17][12] += mat_term;
	mat_term = (sf1.dN1_dy * gp1.n * sf1.N3 + sf2.dN1_dy * gp2.n * sf2.N3
			  + sf3.dN1_dy * gp3.n * sf3.N3 + sf4.dN1_dy * gp4.n * sf4.N3) * -gp_w;
	kmat[12][18] += mat_term;
	kmat[18][12] += mat_term;
	mat_term = (sf1.dN1_dy * gp1.n * sf1.N4 + sf2.dN1_dy * gp2.n * sf2.N4
			  + sf3.dN1_dy * gp3.n * sf3.N4 + sf4.dN1_dy * gp4.n * sf4.N4) * -gp_w;
	kmat[12][19] += mat_term;
	kmat[19][12] += mat_term;
	mat_term = (sf1.dN2_dy * gp1.n * sf1.N1 + sf2.dN2_dy * gp2.n * sf2.N1
			  + sf3.dN2_dy * gp3.n * sf3.N1 + sf4.dN2_dy * gp4.n * sf4.N1) * -gp_w;
	kmat[13][16] += mat_term;
	kmat[16][13] += mat_term;
	mat_term = (sf1.dN2_dy * gp1.n * sf1.N2 + sf2.dN2_dy * gp2.n * sf2.N2
			  + sf3.dN2_dy * gp3.n * sf3.N2 + sf4.dN2_dy * gp4.n * sf4.N2) * -gp_w;
	kmat[13][17] += mat_term;
	kmat[17][13] += mat_term;
	mat_term = (sf1.dN2_dy * gp1.n * sf1.N3 + sf2.dN2_dy * gp2.n * sf2.N3
			  + sf3.dN2_dy * gp3.n * sf3.N3 + sf4.dN2_dy * gp4.n * sf4.N3) * -gp_w;
	kmat[13][18] += mat_term;
	kmat[18][13] += mat_term;
	mat_term = (sf1.dN2_dy * gp1.n * sf1.N4 + sf2.dN2_dy * gp2.n * sf2.N4
			  + sf3.dN2_dy * gp3.n * sf3.N4 + sf4.dN2_dy * gp4.n * sf4.N4) * -gp_w;
	kmat[13][19] += mat_term;
	kmat[19][13] += mat_term;
	mat_term = (sf1.dN3_dy * gp1.n * sf1.N1 + sf2.dN3_dy * gp2.n * sf2.N1
			  + sf3.dN3_dy * gp3.n * sf3.N1 + sf4.dN3_dy * gp4.n * sf4.N1) * -gp_w;
	kmat[14][16] += mat_term;
	kmat[16][14] += mat_term;
	mat_term = (sf1.dN3_dy * gp1.n * sf1.N2 + sf2.dN3_dy * gp2.n * sf2.N2
			  + sf3.dN3_dy * gp3.n * sf3.N2 + sf4.dN3_dy * gp4.n * sf4.N2) * -gp_w;
	kmat[14][17] += mat_term;
	kmat[17][14] += mat_term;
	mat_term = (sf1.dN3_dy * gp1.n * sf1.N3 + sf2.dN3_dy * gp2.n * sf2.N3
			  + sf3.dN3_dy * gp3.n * sf3.N3 + sf4.dN3_dy * gp4.n * sf4.N3) * -gp_w;
	kmat[14][18] += mat_term;
	kmat[18][14] += mat_term;
	mat_term = (sf1.dN3_dy * gp1.n * sf1.N4 + sf2.dN3_dy * gp2.n * sf2.N4
			  + sf3.dN3_dy * gp3.n * sf3.N4 + sf4.dN3_dy * gp4.n * sf4.N4) * -gp_w;
	kmat[14][19] += mat_term;
	kmat[19][14] += mat_term;
	mat_term = (sf1.dN4_dy * gp1.n * sf1.N1 + sf2.dN4_dy * gp2.n * sf2.N1
			  + sf3.dN4_dy * gp3.n * sf3.N1 + sf4.dN4_dy * gp4.n * sf4.N1) * -gp_w;
	kmat[15][16] += mat_term;
	kmat[16][15] += mat_term;
	mat_term = (sf1.dN4_dy * gp1.n * sf1.N2 + sf2.dN4_dy * gp2.n * sf2.N2
			  + sf3.dN4_dy * gp3.n * sf3.N2 + sf4.dN4_dy * gp4.n * sf4.N2) * -gp_w;
	kmat[15][17] += mat_term;
	kmat[17][15] += mat_term;
	mat_term = (sf1.dN4_dy * gp1.n * sf1.N3 + sf2.dN4_dy * gp2.n * sf2.N3
			  + sf3.dN4_dy * gp3.n * sf3.N3 + sf4.dN4_dy * gp4.n * sf4.N3) * -gp_w;
	kmat[15][18] += mat_term;
	kmat[18][15] += mat_term;
	mat_term = (sf1.dN4_dy * gp1.n * sf1.N4 + sf2.dN4_dy * gp2.n * sf2.N4
			  + sf3.dN4_dy * gp3.n * sf3.N4 + sf4.dN4_dy * gp4.n * sf4.N4) * -gp_w;
	kmat[15][19] += mat_term;
	kmat[19][15] += mat_term;
	// P
	mat_term = (sf1.N1 * gp1.n / gp1.Kf * sf1.N1 + sf2.N1 * gp2.n / gp2.Kf * sf2.N1
			  + sf3.N1 * gp3.n / gp3.Kf * sf3.N1 + sf4.N1 * gp4.n / gp4.Kf * sf4.N1) * -gp_w;
	kmat[16][16] += mat_term;
	mat_term = (sf1.N1 * gp1.n / gp1.Kf * sf1.N2 + sf2.N1 * gp2.n / gp2.Kf * sf2.N2
			  + sf3.N1 * gp3.n / gp3.Kf * sf3.N2 + sf4.N1 * gp4.n / gp4.Kf * sf4.N2) * -gp_w;
	kmat[16][17] += mat_term;
	kmat[17][16] += mat_term;
	mat_term = (sf1.N1 * gp1.n / gp1.Kf * sf1.N3 + sf2.N1 * gp2.n / gp2.Kf * sf2.N3
			  + sf3.N1 * gp3.n / gp3.Kf * sf3.N3 + sf4.N1 * gp4.n / gp4.Kf * sf4.N3) * -gp_w;
	kmat[16][18] += mat_term;
	kmat[18][16] += mat_term;
	mat_term = (sf1.N1 * gp1.n / gp1.Kf * sf1.N4 + sf2.N1 * gp2.n / gp2.Kf * sf2.N4
			  + sf3.N1 * gp3.n / gp3.Kf * sf3.N4 + sf4.N1 * gp4.n / gp4.Kf * sf4.N4) * -gp_w;
	kmat[16][19] += mat_term;
	kmat[19][16] += mat_term;
	mat_term = (sf1.N2 * gp1.n / gp1.Kf * sf1.N2 + sf2.N2 * gp2.n / gp2.Kf * sf2.N2
			  + sf3.N2 * gp3.n / gp3.Kf * sf3.N2 + sf4.N2 * gp4.n / gp4.Kf * sf4.N2) * -gp_w;
	kmat[17][17] += mat_term;
	mat_term = (sf1.N2 * gp1.n / gp1.Kf * sf1.N3 + sf2.N2 * gp2.n / gp2.Kf * sf2.N3
			  + sf3.N2 * gp3.n / gp3.Kf * sf3.N3 + sf4.N2 * gp4.n / gp4.Kf * sf4.N3) * -gp_w;
	kmat[17][18] += mat_term;
	kmat[18][17] += mat_term;
	mat_term = (sf1.N2 * gp1.n / gp1.Kf * sf1.N4 + sf2.N2 * gp2.n / gp2.Kf * sf2.N4
			  + sf3.N2 * gp3.n / gp3.Kf * sf3.N4 + sf4.N2 * gp4.n / gp4.Kf * sf4.N4) * -gp_w;
	kmat[17][19] += mat_term;
	kmat[19][17] += mat_term;
	mat_term = (sf1.N3 * gp1.n / gp1.Kf * sf1.N3 + sf2.N3 * gp2.n / gp2.Kf * sf2.N3
			  + sf3.N3 * gp3.n / gp3.Kf * sf3.N3 + sf4.N3 * gp4.n / gp4.Kf * sf4.N3) * -gp_w;
	kmat[18][18] += mat_term;
	mat_term = (sf1.N3 * gp1.n / gp1.Kf * sf1.N4 + sf2.N3 * gp2.n / gp2.Kf * sf2.N4
			  + sf3.N3 * gp3.n / gp3.Kf * sf3.N4 + sf4.N3 * gp4.n / gp4.Kf * sf4.N4) * -gp_w;
	kmat[18][19] += mat_term;
	kmat[19][18] += mat_term;
	mat_term = (sf1.N4 * gp1.n / gp1.Kf * sf1.N4 + sf2.N4 * gp2.n / gp2.Kf * sf2.N4
			  + sf3.N4 * gp3.n / gp3.Kf * sf3.N4 + sf4.N4 * gp4.n / gp4.Kf * sf4.N4) * -gp_w;
	kmat[19][19] += mat_term;

	//print_mat(kmat, out_file);

	// Form elem force vector
	// f_int
	// solid x
	fvec[0] = (sf1.dN1_dx * (gp1.s11 - (1.0 - gp1.n) * gp1.p) + sf1.dN1_dy * gp1.s12
			 + sf2.dN1_dx * (gp2.s11 - (1.0 - gp2.n) * gp2.p) + sf2.dN1_dy * gp2.s12
			 + sf3.dN1_dx * (gp3.s11 - (1.0 - gp3.n) * gp3.p) + sf3.dN1_dy * gp3.s12
			 + sf4.dN1_dx * (gp4.s11 - (1.0 - gp4.n) * gp4.p) + sf4.dN1_dy * gp4.s12) * -gp_w;
	fvec[1] = (sf1.dN2_dx * (gp1.s11 - (1.0 - gp1.n) * gp1.p) + sf1.dN2_dy * gp1.s12
			 + sf2.dN2_dx * (gp2.s11 - (1.0 - gp2.n) * gp2.p) + sf2.dN2_dy * gp2.s12
			 + sf3.dN2_dx * (gp3.s11 - (1.0 - gp3.n) * gp3.p) + sf3.dN2_dy * gp3.s12
			 + sf4.dN2_dx * (gp4.s11 - (1.0 - gp4.n) * gp4.p) + sf4.dN2_dy * gp4.s12) * -gp_w;
	fvec[2] = (sf1.dN3_dx * (gp1.s11 - (1.0 - gp1.n) * gp1.p) + sf1.dN3_dy * gp1.s12
			 + sf2.dN3_dx * (gp2.s11 - (1.0 - gp2.n) * gp2.p) + sf2.dN3_dy * gp2.s12
			 + sf3.dN3_dx * (gp3.s11 - (1.0 - gp3.n) * gp3.p) + sf3.dN3_dy * gp3.s12
			 + sf4.dN3_dx * (gp4.s11 - (1.0 - gp4.n) * gp4.p) + sf4.dN3_dy * gp4.s12) * -gp_w;
	fvec[3] = (sf1.dN4_dx * (gp1.s11 - (1.0 - gp1.n) * gp1.p) + sf1.dN4_dy * gp1.s12
			 + sf2.dN4_dx * (gp2.s11 - (1.0 - gp2.n) * gp2.p) + sf2.dN4_dy * gp2.s12
			 + sf3.dN4_dx * (gp3.s11 - (1.0 - gp3.n) * gp3.p) + sf3.dN4_dy * gp3.s12
			 + sf4.dN4_dx * (gp4.s11 - (1.0 - gp4.n) * gp4.p) + sf4.dN4_dy * gp4.s12) * -gp_w;
	// solid y
	fvec[4] = (sf1.dN1_dx * gp1.s12 + sf1.dN1_dy * (gp1.s22 - (1.0 - gp1.n) * gp1.p)
			 + sf2.dN1_dx * gp2.s12 + sf2.dN1_dy * (gp2.s22 - (1.0 - gp2.n) * gp2.p)
			 + sf3.dN1_dx * gp3.s12 + sf3.dN1_dy * (gp3.s22 - (1.0 - gp3.n) * gp3.p)
			 + sf4.dN1_dx * gp4.s12 + sf4.dN1_dy * (gp4.s22 - (1.0 - gp4.n) * gp4.p)) * -gp_w;
	fvec[5] = (sf1.dN2_dx * gp1.s12 + sf1.dN2_dy * (gp1.s22 - (1.0 - gp1.n) * gp1.p)
			 + sf2.dN2_dx * gp2.s12 + sf2.dN2_dy * (gp2.s22 - (1.0 - gp2.n) * gp2.p)
			 + sf3.dN2_dx * gp3.s12 + sf3.dN2_dy * (gp3.s22 - (1.0 - gp3.n) * gp3.p)
			 + sf4.dN2_dx * gp4.s12 + sf4.dN2_dy * (gp4.s22 - (1.0 - gp4.n) * gp4.p)) * -gp_w;
	fvec[6] = (sf1.dN3_dx * gp1.s12 + sf1.dN3_dy * (gp1.s22 - (1.0 - gp1.n) * gp1.p)
			 + sf2.dN3_dx * gp2.s12 + sf2.dN3_dy * (gp2.s22 - (1.0 - gp2.n) * gp2.p)
			 + sf3.dN3_dx * gp3.s12 + sf3.dN3_dy * (gp3.s22 - (1.0 - gp3.n) * gp3.p)
			 + sf4.dN3_dx * gp4.s12 + sf4.dN3_dy * (gp4.s22 - (1.0 - gp4.n) * gp4.p)) * -gp_w;
	fvec[7] = (sf1.dN4_dx * gp1.s12 + sf1.dN4_dy * (gp1.s22 - (1.0 - gp1.n) * gp1.p)
			 + sf2.dN4_dx * gp2.s12 + sf2.dN4_dy * (gp2.s22 - (1.0 - gp2.n) * gp2.p)
			 + sf3.dN4_dx * gp3.s12 + sf3.dN4_dy * (gp3.s22 - (1.0 - gp3.n) * gp3.p)
			 + sf4.dN4_dx * gp4.s12 + sf4.dN4_dy * (gp4.s22 - (1.0 - gp4.n) * gp4.p)) * -gp_w;
	// fluid x
	fvec[8] = (sf1.dN1_dx * gp1.n * -gp1.p
			 + sf2.dN1_dx * gp2.n * -gp2.p
			 + sf3.dN1_dx * gp3.n * -gp3.p
			 + sf4.dN1_dx * gp4.n * -gp4.p) * -gp_w;
	fvec[9] = (sf1.dN2_dx * gp1.n * -gp1.p
			 + sf2.dN2_dx * gp2.n * -gp2.p
			 + sf3.dN2_dx * gp3.n * -gp3.p
			 + sf4.dN2_dx * gp4.n * -gp4.p) * -gp_w;
	fvec[10] = (sf1.dN3_dx * gp1.n * -gp1.p
			  + sf2.dN3_dx * gp2.n * -gp2.p
			  + sf3.dN3_dx * gp3.n * -gp3.p
			  + sf4.dN3_dx * gp4.n * -gp4.p) * -gp_w;
	fvec[11] = (sf1.dN4_dx * gp1.n * -gp1.p
			  + sf2.dN4_dx * gp2.n * -gp2.p
			  + sf3.dN4_dx * gp3.n * -gp3.p
			  + sf4.dN4_dx * gp4.n * -gp4.p) * -gp_w;
	// fluid y
	fvec[12] = (sf1.dN1_dy * gp1.n * -gp1.p
			  + sf2.dN1_dy * gp2.n * -gp2.p
			  + sf3.dN1_dy * gp3.n * -gp3.p
			  + sf4.dN1_dy * gp4.n * -gp4.p) * -gp_w;
	fvec[13] = (sf1.dN2_dy * gp1.n * -gp1.p
			  + sf2.dN2_dy * gp2.n * -gp2.p
			  + sf3.dN2_dy * gp3.n * -gp3.p
			  + sf4.dN2_dy * gp4.n * -gp4.p) * -gp_w;
	fvec[14] = (sf1.dN3_dy * gp1.n * -gp1.p
			  + sf2.dN3_dy * gp2.n * -gp2.p
			  + sf3.dN3_dy * gp3.n * -gp3.p
			  + sf4.dN3_dy * gp4.n * -gp4.p) * -gp_w;
	fvec[15] = (sf1.dN4_dy * gp1.n * -gp1.p
			  + sf2.dN4_dy * gp2.n * -gp2.p
			  + sf3.dN4_dy * gp3.n * -gp3.p
			  + sf4.dN4_dy * gp4.n * -gp4.p) * -gp_w;
	// p
	fvec[16] = 0.0;
	fvec[17] = 0.0;
	fvec[18] = 0.0;
	fvec[19] = 0.0;

	double coef1, coef2, tmp1, tmp2;
	// a
	coef1 = 1.0 / (2.0 * beta) - 1.0;
	coef2 = - dtime * (1.0 - gamma / (2.0 * beta));
	// solid x
	tmp1 = (sf1.N1 * (1.0 - gp1.n) * gp1.density_s * gp1.ax_s
		  + sf2.N1 * (1.0 - gp2.n) * gp2.density_s * gp2.ax_s
		  + sf3.N1 * (1.0 - gp3.n) * gp3.density_s * gp3.ax_s
		  + sf4.N1 * (1.0 - gp4.n) * gp4.density_s * gp4.ax_s) * gp_w;
	tmp2 = (sf1.N1 * gp1.n * gp1.n * gp1.miu / gp1.k * (gp1.ax_s - gp1.ax_f)
		  + sf2.N1 * gp2.n * gp2.n * gp2.miu / gp2.k * (gp2.ax_s - gp2.ax_f)
		  + sf3.N1 * gp3.n * gp3.n * gp3.miu / gp3.k * (gp3.ax_s - gp3.ax_f)
		  + sf4.N1 * gp4.n * gp4.n * gp4.miu / gp4.k * (gp4.ax_s - gp4.ax_f)) * gp_w;
	fvec[0] += coef1 * tmp1 + coef2 * tmp2;
	tmp1 = (sf1.N2 * (1.0 - gp1.n) * gp1.density_s * gp1.ax_s
		  + sf2.N2 * (1.0 - gp2.n) * gp2.density_s * gp2.ax_s
		  + sf3.N2 * (1.0 - gp3.n) * gp3.density_s * gp3.ax_s
		  + sf4.N2 * (1.0 - gp4.n) * gp4.density_s * gp4.ax_s) * gp_w;
	tmp2 = (sf1.N2 * gp1.n * gp1.n * gp1.miu / gp1.k * (gp1.ax_s - gp1.ax_f)
		  + sf2.N2 * gp2.n * gp2.n * gp2.miu / gp2.k * (gp2.ax_s - gp2.ax_f)
		  + sf3.N2 * gp3.n * gp3.n * gp3.miu / gp3.k * (gp3.ax_s - gp3.ax_f)
		  + sf4.N2 * gp4.n * gp4.n * gp4.miu / gp4.k * (gp4.ax_s - gp4.ax_f)) * gp_w;
	fvec[1] += coef1 * tmp1 + coef2 * tmp2;
	tmp1 = (sf1.N3 * (1.0 - gp1.n) * gp1.density_s * gp1.ax_s
		  + sf2.N3 * (1.0 - gp2.n) * gp2.density_s * gp2.ax_s
		  + sf3.N3 * (1.0 - gp3.n) * gp3.density_s * gp3.ax_s
		  + sf4.N3 * (1.0 - gp4.n) * gp4.density_s * gp4.ax_s) * gp_w;
	tmp2 = (sf1.N3 * gp1.n * gp1.n * gp1.miu / gp1.k * (gp1.ax_s - gp1.ax_f)
		  + sf2.N3 * gp2.n * gp2.n * gp2.miu / gp2.k * (gp2.ax_s - gp2.ax_f)
		  + sf3.N3 * gp3.n * gp3.n * gp3.miu / gp3.k * (gp3.ax_s - gp3.ax_f)
		  + sf4.N3 * gp4.n * gp4.n * gp4.miu / gp4.k * (gp4.ax_s - gp4.ax_f)) * gp_w;
	fvec[2] += coef1 * tmp1 + coef2 * tmp2;
	tmp1 = (sf1.N4 * (1.0 - gp1.n) * gp1.density_s * gp1.ax_s
		  + sf2.N4 * (1.0 - gp2.n) * gp2.density_s * gp2.ax_s
		  + sf3.N4 * (1.0 - gp3.n) * gp3.density_s * gp3.ax_s
		  + sf4.N4 * (1.0 - gp4.n) * gp4.density_s * gp4.ax_s) * gp_w;
	tmp2 = (sf1.N4 * gp1.n * gp1.n * gp1.miu / gp1.k * (gp1.ax_s - gp1.ax_f)
		  + sf2.N4 * gp2.n * gp2.n * gp2.miu / gp2.k * (gp2.ax_s - gp2.ax_f)
		  + sf3.N4 * gp3.n * gp3.n * gp3.miu / gp3.k * (gp3.ax_s - gp3.ax_f)
		  + sf4.N4 * gp4.n * gp4.n * gp4.miu / gp4.k * (gp4.ax_s - gp4.ax_f)) * gp_w;
	fvec[3] += coef1 * tmp1 + coef2 * tmp2;
	// solid y
	tmp1 = (sf1.N1 * (1.0 - gp1.n) * gp1.density_s * gp1.ay_s
		  + sf2.N1 * (1.0 - gp2.n) * gp2.density_s * gp2.ay_s
		  + sf3.N1 * (1.0 - gp3.n) * gp3.density_s * gp3.ay_s
		  + sf4.N1 * (1.0 - gp4.n) * gp4.density_s * gp4.ay_s) * gp_w;
	tmp2 = (sf1.N1 * gp1.n * gp1.n * gp1.miu / gp1.k * (gp1.ay_s - gp1.ay_f)
		  + sf2.N1 * gp2.n * gp2.n * gp2.miu / gp2.k * (gp2.ay_s - gp2.ay_f)
		  + sf3.N1 * gp3.n * gp3.n * gp3.miu / gp3.k * (gp3.ay_s - gp3.ay_f)
		  + sf4.N1 * gp4.n * gp4.n * gp4.miu / gp4.k * (gp4.ay_s - gp4.ay_f)) * gp_w;
	fvec[4] += coef1 * tmp1 + coef2 * tmp2;
	tmp1 = (sf1.N2 * (1.0 - gp1.n) * gp1.density_s * gp1.ay_s
		  + sf2.N2 * (1.0 - gp2.n) * gp2.density_s * gp2.ay_s
		  + sf3.N2 * (1.0 - gp3.n) * gp3.density_s * gp3.ay_s
		  + sf4.N2 * (1.0 - gp4.n) * gp4.density_s * gp4.ay_s) * gp_w;
	tmp2 = (sf1.N2 * gp1.n * gp1.n * gp1.miu / gp1.k * (gp1.ay_s - gp1.ay_f)
		  + sf2.N2 * gp2.n * gp2.n * gp2.miu / gp2.k * (gp2.ay_s - gp2.ay_f)
		  + sf3.N2 * gp3.n * gp3.n * gp3.miu / gp3.k * (gp3.ay_s - gp3.ay_f)
		  + sf4.N2 * gp4.n * gp4.n * gp4.miu / gp4.k * (gp4.ay_s - gp4.ay_f)) * gp_w;
	fvec[5] += coef1 * tmp1 + coef2 * tmp2;
	tmp1 = (sf1.N3 * (1.0 - gp1.n) * gp1.density_s * gp1.ay_s
		  + sf2.N3 * (1.0 - gp2.n) * gp2.density_s * gp2.ay_s
		  + sf3.N3 * (1.0 - gp3.n) * gp3.density_s * gp3.ay_s
		  + sf4.N3 * (1.0 - gp4.n) * gp4.density_s * gp4.ay_s) * gp_w;
	tmp2 = (sf1.N3 * gp1.n * gp1.n * gp1.miu / gp1.k * (gp1.ay_s - gp1.ay_f)
		  + sf2.N3 * gp2.n * gp2.n * gp2.miu / gp2.k * (gp2.ay_s - gp2.ay_f)
		  + sf3.N3 * gp3.n * gp3.n * gp3.miu / gp3.k * (gp3.ay_s - gp3.ay_f)
		  + sf4.N3 * gp4.n * gp4.n * gp4.miu / gp4.k * (gp4.ay_s - gp4.ay_f)) * gp_w;
	fvec[6] += coef1 * tmp1 + coef2 * tmp2;
	tmp1 = (sf1.N4 * (1.0 - gp1.n) * gp1.density_s * gp1.ay_s
		  + sf2.N4 * (1.0 - gp2.n) * gp2.density_s * gp2.ay_s
		  + sf3.N4 * (1.0 - gp3.n) * gp3.density_s * gp3.ay_s
		  + sf4.N4 * (1.0 - gp4.n) * gp4.density_s * gp4.ay_s) * gp_w;
	tmp2 = (sf1.N4 * gp1.n * gp1.n * gp1.miu / gp1.k * (gp1.ay_s - gp1.ay_f)
		  + sf2.N4 * gp2.n * gp2.n * gp2.miu / gp2.k * (gp2.ay_s - gp2.ay_f)
		  + sf3.N4 * gp3.n * gp3.n * gp3.miu / gp3.k * (gp3.ay_s - gp3.ay_f)
		  + sf4.N4 * gp4.n * gp4.n * gp4.miu / gp4.k * (gp4.ay_s - gp4.ay_f)) * gp_w;
	fvec[7] += coef1 * tmp1 + coef2 * tmp2;
	// fluid x
	tmp1 = (sf1.N1 * gp1.n * gp1.density_f * gp1.ax_f
		  + sf2.N1 * gp2.n * gp2.density_f * gp2.ax_f
		  + sf3.N1 * gp3.n * gp3.density_f * gp3.ax_f
		  + sf4.N1 * gp4.n * gp4.density_f * gp4.ax_f) * gp_w;
	tmp2 = (sf1.N1 * gp1.n * gp1.n * gp1.miu / gp1.k * (gp1.ax_f - gp1.ax_s)
		  + sf2.N1 * gp2.n * gp2.n * gp2.miu / gp2.k * (gp2.ax_f - gp2.ax_s)
		  + sf3.N1 * gp3.n * gp3.n * gp3.miu / gp3.k * (gp3.ax_f - gp3.ax_s)
		  + sf4.N1 * gp4.n * gp4.n * gp4.miu / gp4.k * (gp4.ax_f - gp4.ax_s)) * gp_w;
	fvec[8] += coef1 * tmp1 + coef2 * tmp2;
	tmp1 = (sf1.N2 * gp1.n * gp1.density_f * gp1.ax_f
		  + sf2.N2 * gp2.n * gp2.density_f * gp2.ax_f
		  + sf3.N2 * gp3.n * gp3.density_f * gp3.ax_f
		  + sf4.N2 * gp4.n * gp4.density_f * gp4.ax_f) * gp_w;
	tmp2 = (sf1.N2 * gp1.n * gp1.n * gp1.miu / gp1.k * (gp1.ax_f - gp1.ax_s)
		  + sf2.N2 * gp2.n * gp2.n * gp2.miu / gp2.k * (gp2.ax_f - gp2.ax_s)
		  + sf3.N2 * gp3.n * gp3.n * gp3.miu / gp3.k * (gp3.ax_f - gp3.ax_s)
		  + sf4.N2 * gp4.n * gp4.n * gp4.miu / gp4.k * (gp4.ax_f - gp4.ax_s)) * gp_w;
	fvec[9] += coef1 * tmp1 + coef2 * tmp2;
	tmp1 = (sf1.N3 * gp1.n * gp1.density_f * gp1.ax_f
		  + sf2.N3 * gp2.n * gp2.density_f * gp2.ax_f
		  + sf3.N3 * gp3.n * gp3.density_f * gp3.ax_f
		  + sf4.N3 * gp4.n * gp4.density_f * gp4.ax_f) * gp_w;
	tmp2 = (sf1.N3 * gp1.n * gp1.n * gp1.miu / gp1.k * (gp1.ax_f - gp1.ax_s)
		  + sf2.N3 * gp2.n * gp2.n * gp2.miu / gp2.k * (gp2.ax_f - gp2.ax_s)
		  + sf3.N3 * gp3.n * gp3.n * gp3.miu / gp3.k * (gp3.ax_f - gp3.ax_s)
		  + sf4.N3 * gp4.n * gp4.n * gp4.miu / gp4.k * (gp4.ax_f - gp4.ax_s)) * gp_w;
	fvec[10] += coef1 * tmp1 + coef2 * tmp2;
	tmp1 = (sf1.N4 * gp1.n * gp1.density_f * gp1.ax_f
		  + sf2.N4 * gp2.n * gp2.density_f * gp2.ax_f
		  + sf3.N4 * gp3.n * gp3.density_f * gp3.ax_f
		  + sf4.N4 * gp4.n * gp4.density_f * gp4.ax_f) * gp_w;
	tmp2 = (sf1.N4 * gp1.n * gp1.n * gp1.miu / gp1.k * (gp1.ax_f - gp1.ax_s)
		  + sf2.N4 * gp2.n * gp2.n * gp2.miu / gp2.k * (gp2.ax_f - gp2.ax_s)
		  + sf3.N4 * gp3.n * gp3.n * gp3.miu / gp3.k * (gp3.ax_f - gp3.ax_s)
		  + sf4.N4 * gp4.n * gp4.n * gp4.miu / gp4.k * (gp4.ax_f - gp4.ax_s)) * gp_w;
	fvec[11] += coef1 * tmp1 + coef2 * tmp2;
	// fluid y
	tmp1 = (sf1.N1 * gp1.n * gp1.density_f * gp1.ay_f
		  + sf2.N1 * gp2.n * gp2.density_f * gp2.ay_f
		  + sf3.N1 * gp3.n * gp3.density_f * gp3.ay_f
		  + sf4.N1 * gp4.n * gp4.density_f * gp4.ay_f) * gp_w;
	tmp2 = (sf1.N1 * gp1.n * gp1.n * gp1.miu / gp1.k * (gp1.ay_f - gp1.ay_s)
		  + sf2.N1 * gp2.n * gp2.n * gp2.miu / gp2.k * (gp2.ay_f - gp2.ay_s)
		  + sf3.N1 * gp3.n * gp3.n * gp3.miu / gp3.k * (gp3.ay_f - gp3.ay_s)
		  + sf4.N1 * gp4.n * gp4.n * gp4.miu / gp4.k * (gp4.ay_f - gp4.ay_s)) * gp_w;
	fvec[12] += coef1 * tmp1 + coef2 * tmp2;
	tmp1 = (sf1.N2 * gp1.n * gp1.density_f * gp1.ay_f
		  + sf2.N2 * gp2.n * gp2.density_f * gp2.ay_f
		  + sf3.N2 * gp3.n * gp3.density_f * gp3.ay_f
		  + sf4.N2 * gp4.n * gp4.density_f * gp4.ay_f) * gp_w;
	tmp2 = (sf1.N2 * gp1.n * gp1.n * gp1.miu / gp1.k * (gp1.ay_f - gp1.ay_s)
		  + sf2.N2 * gp2.n * gp2.n * gp2.miu / gp2.k * (gp2.ay_f - gp2.ay_s)
		  + sf3.N2 * gp3.n * gp3.n * gp3.miu / gp3.k * (gp3.ay_f - gp3.ay_s)
		  + sf4.N2 * gp4.n * gp4.n * gp4.miu / gp4.k * (gp4.ay_f - gp4.ay_s)) * gp_w;
	fvec[13] += coef1 * tmp1 + coef2 * tmp2;
	tmp1 = (sf1.N3 * gp1.n * gp1.density_f * gp1.ay_f
		  + sf2.N3 * gp2.n * gp2.density_f * gp2.ay_f
		  + sf3.N3 * gp3.n * gp3.density_f * gp3.ay_f
		  + sf4.N3 * gp4.n * gp4.density_f * gp4.ay_f) * gp_w;
	tmp2 = (sf1.N3 * gp1.n * gp1.n * gp1.miu / gp1.k * (gp1.ay_f - gp1.ay_s)
		  + sf2.N3 * gp2.n * gp2.n * gp2.miu / gp2.k * (gp2.ay_f - gp2.ay_s)
		  + sf3.N3 * gp3.n * gp3.n * gp3.miu / gp3.k * (gp3.ay_f - gp3.ay_s)
		  + sf4.N3 * gp4.n * gp4.n * gp4.miu / gp4.k * (gp4.ay_f - gp4.ay_s)) * gp_w;
	fvec[14] += coef1 * tmp1 + coef2 * tmp2;
	tmp1 = (sf1.N4 * gp1.n * gp1.density_f * gp1.ay_f
		  + sf2.N4 * gp2.n * gp2.density_f * gp2.ay_f
		  + sf3.N4 * gp3.n * gp3.density_f * gp3.ay_f
		  + sf4.N4 * gp4.n * gp4.density_f * gp4.ay_f) * gp_w;
	tmp2 = (sf1.N4 * gp1.n * gp1.n * gp1.miu / gp1.k * (gp1.ay_f - gp1.ay_s)
		  + sf2.N4 * gp2.n * gp2.n * gp2.miu / gp2.k * (gp2.ay_f - gp2.ay_s)
		  + sf3.N4 * gp3.n * gp3.n * gp3.miu / gp3.k * (gp3.ay_f - gp3.ay_s)
		  + sf4.N4 * gp4.n * gp4.n * gp4.miu / gp4.k * (gp4.ay_f - gp4.ay_s)) * gp_w;
	fvec[15] += coef1 * tmp1 + coef2 * tmp2;

	// v
	coef1 = 1.0 / (beta * dtime);
	coef2 = gamma / beta - 1.0;
	// solid x
	tmp1 = (sf1.N1 * (1.0 - gp1.n) * gp1.density_s * gp1.vx_s
		  + sf2.N1 * (1.0 - gp2.n) * gp2.density_s * gp2.vx_s
		  + sf3.N1 * (1.0 - gp3.n) * gp3.density_s * gp3.vx_s
		  + sf4.N1 * (1.0 - gp4.n) * gp4.density_s * gp4.vx_s) * gp_w;
	tmp2 = (sf1.N1 * gp1.n * gp1.n * gp1.miu / gp1.k * (gp1.vx_s - gp1.vx_f)
		  + sf2.N1 * gp2.n * gp2.n * gp2.miu / gp2.k * (gp2.vx_s - gp2.vx_f)
		  + sf3.N1 * gp3.n * gp3.n * gp3.miu / gp3.k * (gp3.vx_s - gp3.vx_f)
		  + sf4.N1 * gp4.n * gp4.n * gp4.miu / gp4.k * (gp4.vx_s - gp4.vx_f)) * gp_w;
	fvec[0] += coef1 * tmp1 + coef2 * tmp2;
	tmp1 = (sf1.N2 * (1.0 - gp1.n) * gp1.density_s * gp1.vx_s
		  + sf2.N2 * (1.0 - gp2.n) * gp2.density_s * gp2.vx_s
		  + sf3.N2 * (1.0 - gp3.n) * gp3.density_s * gp3.vx_s
		  + sf4.N2 * (1.0 - gp4.n) * gp4.density_s * gp4.vx_s) * gp_w;
	tmp2 = (sf1.N2 * gp1.n * gp1.n * gp1.miu / gp1.k * (gp1.vx_s - gp1.vx_f)
		  + sf2.N2 * gp2.n * gp2.n * gp2.miu / gp2.k * (gp2.vx_s - gp2.vx_f)
		  + sf3.N2 * gp3.n * gp3.n * gp3.miu / gp3.k * (gp3.vx_s - gp3.vx_f)
		  + sf4.N2 * gp4.n * gp4.n * gp4.miu / gp4.k * (gp4.vx_s - gp4.vx_f)) * gp_w;
	fvec[1] += coef1 * tmp1 + coef2 * tmp2;
	tmp1 = (sf1.N3 * (1.0 - gp1.n) * gp1.density_s * gp1.vx_s
		  + sf2.N3 * (1.0 - gp2.n) * gp2.density_s * gp2.vx_s
		  + sf3.N3 * (1.0 - gp3.n) * gp3.density_s * gp3.vx_s
		  + sf4.N3 * (1.0 - gp4.n) * gp4.density_s * gp4.vx_s) * gp_w;
	tmp2 = (sf1.N3 * gp1.n * gp1.n * gp1.miu / gp1.k * (gp1.vx_s - gp1.vx_f)
		  + sf2.N3 * gp2.n * gp2.n * gp2.miu / gp2.k * (gp2.vx_s - gp2.vx_f)
		  + sf3.N3 * gp3.n * gp3.n * gp3.miu / gp3.k * (gp3.vx_s - gp3.vx_f)
		  + sf4.N3 * gp4.n * gp4.n * gp4.miu / gp4.k * (gp4.vx_s - gp4.vx_f)) * gp_w;
	fvec[2] += coef1 * tmp1 + coef2 * tmp2;
	tmp1 = (sf1.N4 * (1.0 - gp1.n) * gp1.density_s * gp1.vx_s
		  + sf2.N4 * (1.0 - gp2.n) * gp2.density_s * gp2.vx_s
		  + sf3.N4 * (1.0 - gp3.n) * gp3.density_s * gp3.vx_s
		  + sf4.N4 * (1.0 - gp4.n) * gp4.density_s * gp4.vx_s) * gp_w;
	tmp2 = (sf1.N4 * gp1.n * gp1.n * gp1.miu / gp1.k * (gp1.vx_s - gp1.vx_f)
		  + sf2.N4 * gp2.n * gp2.n * gp2.miu / gp2.k * (gp2.vx_s - gp2.vx_f)
		  + sf3.N4 * gp3.n * gp3.n * gp3.miu / gp3.k * (gp3.vx_s - gp3.vx_f)
		  + sf4.N4 * gp4.n * gp4.n * gp4.miu / gp4.k * (gp4.vx_s - gp4.vx_f)) * gp_w;
	fvec[3] += coef1 * tmp1 + coef2 * tmp2;
	// solid y
	tmp1 = (sf1.N1 * (1.0 - gp1.n) * gp1.density_s * gp1.vy_s
		  + sf2.N1 * (1.0 - gp2.n) * gp2.density_s * gp2.vy_s
		  + sf3.N1 * (1.0 - gp3.n) * gp3.density_s * gp3.vy_s
		  + sf4.N1 * (1.0 - gp4.n) * gp4.density_s * gp4.vy_s) * gp_w;
	tmp2 = (sf1.N1 * gp1.n * gp1.n * gp1.miu / gp1.k * (gp1.vy_s - gp1.vy_f)
		  + sf2.N1 * gp2.n * gp2.n * gp2.miu / gp2.k * (gp2.vy_s - gp2.vy_f)
		  + sf3.N1 * gp3.n * gp3.n * gp3.miu / gp3.k * (gp3.vy_s - gp3.vy_f)
		  + sf4.N1 * gp4.n * gp4.n * gp4.miu / gp4.k * (gp4.vy_s - gp4.vy_f)) * gp_w;
	fvec[4] += coef1 * tmp1 + coef2 * tmp2;
	tmp1 = (sf1.N2 * (1.0 - gp1.n) * gp1.density_s * gp1.vy_s
		  + sf2.N2 * (1.0 - gp2.n) * gp2.density_s * gp2.vy_s
		  + sf3.N2 * (1.0 - gp3.n) * gp3.density_s * gp3.vy_s
		  + sf4.N2 * (1.0 - gp4.n) * gp4.density_s * gp4.vy_s) * gp_w;
	tmp2 = (sf1.N2 * gp1.n * gp1.n * gp1.miu / gp1.k * (gp1.vy_s - gp1.vy_f)
		  + sf2.N2 * gp2.n * gp2.n * gp2.miu / gp2.k * (gp2.vy_s - gp2.vy_f)
		  + sf3.N2 * gp3.n * gp3.n * gp3.miu / gp3.k * (gp3.vy_s - gp3.vy_f)
		  + sf4.N2 * gp4.n * gp4.n * gp4.miu / gp4.k * (gp4.vy_s - gp4.vy_f)) * gp_w;
	fvec[5] += coef1 * tmp1 + coef2 * tmp2;
	tmp1 = (sf1.N3 * (1.0 - gp1.n) * gp1.density_s * gp1.vy_s
		  + sf2.N3 * (1.0 - gp2.n) * gp2.density_s * gp2.vy_s
		  + sf3.N3 * (1.0 - gp3.n) * gp3.density_s * gp3.vy_s
		  + sf4.N3 * (1.0 - gp4.n) * gp4.density_s * gp4.vy_s) * gp_w;
	tmp2 = (sf1.N3 * gp1.n * gp1.n * gp1.miu / gp1.k * (gp1.vy_s - gp1.vy_f)
		  + sf2.N3 * gp2.n * gp2.n * gp2.miu / gp2.k * (gp2.vy_s - gp2.vy_f)
		  + sf3.N3 * gp3.n * gp3.n * gp3.miu / gp3.k * (gp3.vy_s - gp3.vy_f)
		  + sf4.N3 * gp4.n * gp4.n * gp4.miu / gp4.k * (gp4.vy_s - gp4.vy_f)) * gp_w;
	fvec[6] += coef1 * tmp1 + coef2 * tmp2;
	tmp1 = (sf1.N4 * (1.0 - gp1.n) * gp1.density_s * gp1.vy_s
		  + sf2.N4 * (1.0 - gp2.n) * gp2.density_s * gp2.vy_s
		  + sf3.N4 * (1.0 - gp3.n) * gp3.density_s * gp3.vy_s
		  + sf4.N4 * (1.0 - gp4.n) * gp4.density_s * gp4.vy_s) * gp_w;
	tmp2 = (sf1.N4 * gp1.n * gp1.n * gp1.miu / gp1.k * (gp1.vy_s - gp1.vy_f)
		  + sf2.N4 * gp2.n * gp2.n * gp2.miu / gp2.k * (gp2.vy_s - gp2.vy_f)
		  + sf3.N4 * gp3.n * gp3.n * gp3.miu / gp3.k * (gp3.vy_s - gp3.vy_f)
		  + sf4.N4 * gp4.n * gp4.n * gp4.miu / gp4.k * (gp4.vy_s - gp4.vy_f)) * gp_w;
	fvec[7] += coef1 * tmp1 + coef2 * tmp2;
	// fluid x
	tmp1 = (sf1.N1 * gp1.n * gp1.density_f * gp1.vx_f
		  + sf2.N1 * gp2.n * gp2.density_f * gp2.vx_f
		  + sf3.N1 * gp3.n * gp3.density_f * gp3.vx_f
		  + sf4.N1 * gp4.n * gp4.density_f * gp4.vx_f) * gp_w;
	tmp2 = (sf1.N1 * gp1.n * gp1.n * gp1.miu / gp1.k * (gp1.vx_f - gp1.vx_s)
		  + sf2.N1 * gp2.n * gp2.n * gp2.miu / gp2.k * (gp2.vx_f - gp2.vx_s)
		  + sf3.N1 * gp3.n * gp3.n * gp3.miu / gp3.k * (gp3.vx_f - gp3.vx_s)
		  + sf4.N1 * gp4.n * gp4.n * gp4.miu / gp4.k * (gp4.vx_f - gp4.vx_s)) * gp_w;
	fvec[8] += coef1 * tmp1 + coef2 * tmp2;
	tmp1 = (sf1.N2 * gp1.n * gp1.density_f * gp1.vx_f
		  + sf2.N2 * gp2.n * gp2.density_f * gp2.vx_f
		  + sf3.N2 * gp3.n * gp3.density_f * gp3.vx_f
		  + sf4.N2 * gp4.n * gp4.density_f * gp4.vx_f) * gp_w;
	tmp2 = (sf1.N2 * gp1.n * gp1.n * gp1.miu / gp1.k * (gp1.vx_f - gp1.vx_s)
		  + sf2.N2 * gp2.n * gp2.n * gp2.miu / gp2.k * (gp2.vx_f - gp2.vx_s)
		  + sf3.N2 * gp3.n * gp3.n * gp3.miu / gp3.k * (gp3.vx_f - gp3.vx_s)
		  + sf4.N2 * gp4.n * gp4.n * gp4.miu / gp4.k * (gp4.vx_f - gp4.vx_s)) * gp_w;
	fvec[9] += coef1 * tmp1 + coef2 * tmp2;
	tmp1 = (sf1.N3 * gp1.n * gp1.density_f * gp1.vx_f
		  + sf2.N3 * gp2.n * gp2.density_f * gp2.vx_f
		  + sf3.N3 * gp3.n * gp3.density_f * gp3.vx_f
		  + sf4.N3 * gp4.n * gp4.density_f * gp4.vx_f) * gp_w;
	tmp2 = (sf1.N3 * gp1.n * gp1.n * gp1.miu / gp1.k * (gp1.vx_f - gp1.vx_s)
		  + sf2.N3 * gp2.n * gp2.n * gp2.miu / gp2.k * (gp2.vx_f - gp2.vx_s)
		  + sf3.N3 * gp3.n * gp3.n * gp3.miu / gp3.k * (gp3.vx_f - gp3.vx_s)
		  + sf4.N3 * gp4.n * gp4.n * gp4.miu / gp4.k * (gp4.vx_f - gp4.vx_s)) * gp_w;
	fvec[10] += coef1 * tmp1 + coef2 * tmp2;
	tmp1 = (sf1.N4 * gp1.n * gp1.density_f * gp1.vx_f
		  + sf2.N4 * gp2.n * gp2.density_f * gp2.vx_f
		  + sf3.N4 * gp3.n * gp3.density_f * gp3.vx_f
		  + sf4.N4 * gp4.n * gp4.density_f * gp4.vx_f) * gp_w;
	tmp2 = (sf1.N4 * gp1.n * gp1.n * gp1.miu / gp1.k * (gp1.vx_f - gp1.vx_s)
		  + sf2.N4 * gp2.n * gp2.n * gp2.miu / gp2.k * (gp2.vx_f - gp2.vx_s)
		  + sf3.N4 * gp3.n * gp3.n * gp3.miu / gp3.k * (gp3.vx_f - gp3.vx_s)
		  + sf4.N4 * gp4.n * gp4.n * gp4.miu / gp4.k * (gp4.vx_f - gp4.vx_s)) * gp_w;
	fvec[11] += coef1 * tmp1 + coef2 * tmp2;
	// fluid y
	tmp1 = (sf1.N1 * gp1.n * gp1.density_f * gp1.vy_f
		  + sf2.N1 * gp2.n * gp2.density_f * gp2.vy_f
		  + sf3.N1 * gp3.n * gp3.density_f * gp3.vy_f
		  + sf4.N1 * gp4.n * gp4.density_f * gp4.vy_f) * gp_w;
	tmp2 = (sf1.N1 * gp1.n * gp1.n * gp1.miu / gp1.k * (gp1.vy_f - gp1.vy_s)
		  + sf2.N1 * gp2.n * gp2.n * gp2.miu / gp2.k * (gp2.vy_f - gp2.vy_s)
		  + sf3.N1 * gp3.n * gp3.n * gp3.miu / gp3.k * (gp3.vy_f - gp3.vy_s)
		  + sf4.N1 * gp4.n * gp4.n * gp4.miu / gp4.k * (gp4.vy_f - gp4.vy_s)) * gp_w;
	fvec[12] += coef1 * tmp1 + coef2 * tmp2;
	tmp1 = (sf1.N2 * gp1.n * gp1.density_f * gp1.vy_f
		  + sf2.N2 * gp2.n * gp2.density_f * gp2.vy_f
		  + sf3.N2 * gp3.n * gp3.density_f * gp3.vy_f
		  + sf4.N2 * gp4.n * gp4.density_f * gp4.vy_f) * gp_w;
	tmp2 = (sf1.N2 * gp1.n * gp1.n * gp1.miu / gp1.k * (gp1.vy_f - gp1.vy_s)
		  + sf2.N2 * gp2.n * gp2.n * gp2.miu / gp2.k * (gp2.vy_f - gp2.vy_s)
		  + sf3.N2 * gp3.n * gp3.n * gp3.miu / gp3.k * (gp3.vy_f - gp3.vy_s)
		  + sf4.N2 * gp4.n * gp4.n * gp4.miu / gp4.k * (gp4.vy_f - gp4.vy_s)) * gp_w;
	fvec[13] += coef1 * tmp1 + coef2 * tmp2;
	tmp1 = (sf1.N3 * gp1.n * gp1.density_f * gp1.vy_f
		  + sf2.N3 * gp2.n * gp2.density_f * gp2.vy_f
		  + sf3.N3 * gp3.n * gp3.density_f * gp3.vy_f
		  + sf4.N3 * gp4.n * gp4.density_f * gp4.vy_f) * gp_w;
	tmp2 = (sf1.N3 * gp1.n * gp1.n * gp1.miu / gp1.k * (gp1.vy_f - gp1.vy_s)
		  + sf2.N3 * gp2.n * gp2.n * gp2.miu / gp2.k * (gp2.vy_f - gp2.vy_s)
		  + sf3.N3 * gp3.n * gp3.n * gp3.miu / gp3.k * (gp3.vy_f - gp3.vy_s)
		  + sf4.N3 * gp4.n * gp4.n * gp4.miu / gp4.k * (gp4.vy_f - gp4.vy_s)) * gp_w;
	fvec[14] += coef1 * tmp1 + coef2 * tmp2;
	tmp1 = (sf1.N4 * gp1.n * gp1.density_f * gp1.vy_f
		  + sf2.N4 * gp2.n * gp2.density_f * gp2.vy_f
		  + sf3.N4 * gp3.n * gp3.density_f * gp3.vy_f
		  + sf4.N4 * gp4.n * gp4.density_f * gp4.vy_f) * gp_w;
	tmp2 = (sf1.N4 * gp1.n * gp1.n * gp1.miu / gp1.k * (gp1.vy_f - gp1.vy_s)
		  + sf2.N4 * gp2.n * gp2.n * gp2.miu / gp2.k * (gp2.vy_f - gp2.vy_s)
		  + sf3.N4 * gp3.n * gp3.n * gp3.miu / gp3.k * (gp3.vy_f - gp3.vy_s)
		  + sf4.N4 * gp4.n * gp4.n * gp4.miu / gp4.k * (gp4.vy_f - gp4.vy_s)) * gp_w;
	fvec[15] += coef1 * tmp1 + coef2 * tmp2;
}