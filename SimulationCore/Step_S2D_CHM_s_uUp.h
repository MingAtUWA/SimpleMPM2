#ifndef __STEP_S2D_CHM_s_uUp_H__
#define __STEP_S2D_CHM_s_uUp_H__

#include <fstream>
#include <Eigen/Sparse>

#include "ItemArrayFast.hpp"
#include "MatrixCoefficientSet.hpp"

#include "Step.h"
#include "Model_S2D_CHM_s_uUp.h"

int solve_substep_S2D_CHM_s_uUp(void *_self);

// for single object only
class Step_S2D_CHM_s_uUp : public Step
{
public:
	typedef Model_S2D_CHM_s_uUp::ShapeFuncValue ShapeFuncValue;
	typedef Model_S2D_CHM_s_uUp::Particle Particle;
	typedef Model_S2D_CHM_s_uUp::Element Element;
	typedef Model_S2D_CHM_s_uUp::Node Node;
	typedef Model_S2D_CHM_s_uUp::DOF DOF;

public:
	Step_S2D_CHM_s_uUp();
	~Step_S2D_CHM_s_uUp();

	inline void set_model(Model_S2D_CHM_s_uUp &md)
	{
		Step::set_model(md);
		model = &md;
	}

	// Restart from previous step
	inline void set_prev_step(Step_S2D_CHM_s_uUp &prev_step)
	{
		Step::set_prev_step(prev_step);
		model = prev_step.model;
	}

	inline void set_dt(double _dt,
		double dt_max_min_raio = 0.1 /*ad hoc number*/,
		double t_tol_r = 0.01)
	{
		max_dt = _dt;
		min_dt = max_dt * dt_max_min_raio;
		time_tol_ratio = t_tol_r;
		time_tol = min_dt * t_tol_r;
		dtime = max_dt;
	}

protected:
	Model_S2D_CHM_s_uUp *model;
	double min_dt, max_dt;

	MatrixCoefficientSet<> g_kmat_coefs;
	MemoryUtilities::ItemArrayFast<size_t> node_g_id_map_mem;
	// for Dirichlet boundary conditions
	MemoryUtilities::ItemArrayFast<double> kmat_col_mem;

protected:
	int init_calculation(void) override;
	friend int solve_substep_S2D_CHM_s_uUp(void *_self);
	int finalize_calculation(void) override;
	
	void map_from_pcl_to_node(void);
	
	void add_internal_force_vec(double fvec[20], Particle &pcl);
	void form_elem_stiffness_mat_and_force_vec(Element &e, double kmat[20][20], double fvec[20]);
	void form_and_solve_problem(Eigen::VectorXd &g_du_vec);

	int requilibration(void);
	int time_marching(void);

public:
	size_t cal_node_num, dof_num;
	inline size_t n_id_to_dof_id(size_t n_id, DOF dof_type) const { return size_t(dof_type) * cal_node_num + n_id; }

public: // for debugging
	void elem_residual_force_vec(Element &e, double fvec[20]);
	void form_global_residual_force(void);
	std::fstream out_file;
};

#endif
#ifdef Keep_Newmark_Coefficients

// constant of Newmark-beta method
#define beta 0.3025
#define gamma 0.6

#endif
#ifdef Keep_Debug_Mat_And_Vec_Output

void print_mat(double mat[3][3])
{
	for (size_t i = 0; i < 3; i++)
	{
		for (size_t j = 0; j < 3; j++)
			std::cout << mat[i][j] << ", ";
		std::cout << "\n";
	}
}

void print_mat(double mat[3][8])
{
	for (size_t i = 0; i < 3; i++)
	{
		for (size_t j = 0; j < 8; j++)
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

void print_sparse_mat(Eigen::SparseMatrix<double> &g_kmat, std::fstream &out_file)
{
	size_t row_num = g_kmat.rows();
	size_t col_num = g_kmat.cols();
	for (size_t row_id = 0; row_id < row_num; ++row_id)
	{
		for (size_t col_id = 0; col_id < col_num; ++col_id)
			out_file << g_kmat.coeff(row_id, col_id) << ", ";
		out_file << "\n";
	}
}

#endif