#ifndef __STEP_S2D_ME_s_up_H__
#define __STEP_S2D_ME_s_up_H__

#include <fstream>
#include <Eigen/Sparse>

#include "ItemArrayFast.hpp"
#include "MatrixCoefficientSet.hpp"

#include "Step.h"
#include "Model_S2D_ME_s_up.h"

int substep_requilibration_S2D_ME_s_up(void *_self, Eigen::VectorXd& g_fvec);
int solve_substep_S2D_ME_s_up(void *_self);

// for single object only
class Step_S2D_ME_s_up : public Step
{
protected:
	int init_calculation(void) override;
	friend int substep_requilibration_S2D_ME_s_up(void *_self, Eigen::VectorXd& g_fvec);
	friend int solve_substep_S2D_ME_s_up(void *_self);
	int finalize_calculation(void) override;

public:
	Step_S2D_ME_s_up();
	~Step_S2D_ME_s_up();

	inline void set_model(Model_S2D_ME_s_up &md)
	{
		Step::set_model(md);
		model = &md;
	}

	// Restart from previous step
	inline void set_prev_step(Step_S2D_ME_s_up &prev_step)
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
		time_tol = max_dt * t_tol_r;
		dtime = max_dt;
	}

protected:
	Model_S2D_ME_s_up *model;
	double min_dt, max_dt;

	MatrixCoefficientSet<> g_kmat_coefs;
	MemoryUtilities::ItemArrayFast<size_t> node_g_id_map_mem;
	// for Dirichlet boundary conditions
	MemoryUtilities::ItemArrayFast<double> kmat_col_mem;
	
	bool is_first_substep;
	void reequilibration_elem_stiffness_mat_and_force_vec(Model_S2D_ME_s_up::Element &e, double kmat[12][12], double fvec[12]);
	void form_elem_stiffness_mat_and_force_vec(Model_S2D_ME_s_up::Element &e, double kmat[12][12], double fvec[12]);
	
public:
	typedef Model_S2D_ME_s_up::DOF DOF;
	size_t cal_node_num;
	inline size_t n_id_to_dof_id(size_t n_id, DOF dof_type) const { return size_t(dof_type) * cal_node_num + n_id; }

public: // for debugging
	std::fstream out_file;
};

#endif
#ifdef KEEP_NEWMARK_BETA_COEFFICIENT
// constant of Newmark-beta method
#define beta 0.3025
#define gamma 0.6
#endif