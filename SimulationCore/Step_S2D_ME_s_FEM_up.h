#ifndef __STEP_S2D_ME_s_FEM_up_H__
#define __STEP_S2D_ME_s_FEM_up_H__

#include "MatrixCoefficientSet.hpp"

#include "Step.h"
#include "Model_S2D_ME_s_FEM_up.h"

// standard MPM
int solve_substep_S2D_ME_s_FEM_up(void *_self);

// for single object only
class Step_S2D_ME_s_FEM_up : public Step
{
protected:
	int init_calculation(void) override;
	friend int solve_substep_S2D_ME_s_FEM_up(void *_self);
	int finalize_calculation(void) override;

public:
	Step_S2D_ME_s_FEM_up();
	~Step_S2D_ME_s_FEM_up();

	inline void set_model(Model_S2D_ME_s_FEM_up &md)
	{
		Step::set_model(md);
		model = &md;
	}

	// Restart from previous step
	inline void set_prev_step(Step_S2D_ME_s_FEM_up &prev_step)
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
		time_tol = min_dt * time_tol_ratio;
		dtime = max_dt;
	}

protected:
	Model_S2D_ME_s_FEM_up *model;
	double min_dt, max_dt;
	
	double *kmat_col;
	MatrixCoefficientSet<> g_kmat_coefs;
	bool is_first_substep;

protected:
	void form_elem_stiffness_mat_and_force_vec(Model_S2D_ME_s_FEM_up::Element &e, double kmat[12][12], double fvec[12]);
};

#endif