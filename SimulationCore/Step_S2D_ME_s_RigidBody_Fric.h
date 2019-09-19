#ifndef __STEP_S2D_ME_S_RIGIDBODY_FRIC_H__
#define __STEP_S2D_ME_S_RIGIDBODY_FRIC_H__

#include "Model_S2D_ME_s_RigidBody_Fric.h"
#include "Step.h"

int solve_substep_S2D_ME_s_RigidBody_Fric(void *_self);

// frictional contact between material point and RigidBody
class Step_S2D_ME_s_RigidBody_Fric : public Step
{
protected:
	Model_S2D_ME_s_RigidBody_Fric *model;
	double min_dt, max_dt;
	double h_elem_raio, h_pcl_ratio;
	// pcl_ratio part not yet finished

public:
	Step_S2D_ME_s_RigidBody_Fric();
	~Step_S2D_ME_s_RigidBody_Fric();

	inline void set_model(Model_S2D_ME_s_RigidBody_Fric &md)
	{
		Step::set_model(md);
		model = &md;
	}

	// Restart from previous step
	inline void set_prev_step(Step_S2D_ME_s_RigidBody_Fric &prev_step)
	{
		Step::set_prev_step(prev_step);
		model = prev_step.model;
	}

	inline void set_dtime(double _dt,
		double dt_max_min_raio = 1.0e-3 /*ad hoc number*/,
		double t_tol_r = 0.01)
	{
		max_dt = _dt;
		min_dt = max_dt * dt_max_min_raio;
		time_tol_ratio = t_tol_r;
		time_tol = min_dt * t_tol_r;
	}

	inline void set_dtime_ratio(double _h_elem_ratio, double _h_pcl_ratio)
	{
		h_elem_raio = _h_elem_ratio;
		h_pcl_ratio = _h_pcl_ratio;
	}

protected:
	int init_calculation(void) override;
	friend int solve_substep_S2D_ME_s_RigidBody_Fric(void *_self);
	int finalize_calculation(void) override;

protected: // keep tract of tangential contact
	ContactStateList contact_state_list;
};

#endif