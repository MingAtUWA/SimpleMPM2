#ifndef __STEP_S2D_CHM_S_FEM_UUP_H__
#define __STEP_S2D_CHM_S_FEM_UUP_H__

#include "Step.h"
#include "Model_S2D_CHM_s_FEM_uUp.h"

// standard MPM
int solve_substep_S2D_CHM_s_FEM_uUp(void *_self);

// for single object only
class Step_S2D_CHM_s_FEM_uUp : public Step
{
protected:
	int init_calculation(void) override;
	friend int solve_substep_S2D_CHM_s_FEM_uUp(void *_self);
	int finalize_calculation(void) override;

public:
	Step_S2D_CHM_s_FEM_uUp();
	~Step_S2D_CHM_s_FEM_uUp();

	inline void set_model(Model_S2D_CHM_s_FEM_uUp &md)
	{
		Step::set_model(md);
		model = &md;
	}

	// Restart from previous step
	inline void set_prev_step(Step_S2D_CHM_s_FEM_uUp &prev_step)
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
	Model_S2D_CHM_s_FEM_uUp *model;
	double min_dt, max_dt;

protected:
	enum class DOF : size_t
	{
		usx = 0,
		usy = 1,
		ufx = 2,
		ufy = 3,
		p = 4
	};
	inline size_t n_id_to_dof_id(size_t n_id, DOF dof_type) const
	{
		return size_t(dof_type) * model->node_num + n_id;
	}

	void form_elem_stiffness_mat(Model_S2D_CHM_s_FEM_uUp::Element &e, double kmat[20][20]);
	void form_elem_force_vec(Model_S2D_CHM_s_FEM_uUp::Element &e, double dt, double fvec[20]);
};

#endif