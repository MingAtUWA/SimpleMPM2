#ifndef __STEP_S2D_CHM_S_H__
#define __STEP_S2D_CHM_S_H__

#include "Step.h"
#include "Model_S2D_CHM_s.h"

int solve_substep_S2D_CHM_s(void *_self);
int solve_substep_S2D_CHM_s_GIMP(void *_self);

// for single object only
class Step_S2D_CHM_s : public Step
{
protected:
	int init_calculation(void) override;
	friend int solve_substep_S2D_CHM_s(void *_self);
	friend int solve_substep_S2D_CHM_s_GIMP(void *_self);
	int finalize_calculation(void) override;

public:
	Step_S2D_CHM_s();
	~Step_S2D_CHM_s();

	inline void use_standard_mpm(void) noexcept {}
	inline void use_gimp(void) noexcept {}

	inline void set_model(Model_S2D_CHM_s &md)
	{
		Step::set_model(md);
		model = &md;
	}

	// Restart from previous step
	inline void set_prev_step(Step_S2D_CHM_s &prev_step)
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

	inline void set_dt_ratio(double _h_elem_ratio, double _h_pcl_ratio)
	{
		h_elem_raio = _h_elem_ratio;
		h_pcl_ratio = _h_pcl_ratio;
	}

protected:
	Model_S2D_CHM_s *model;
	double min_dt, max_dt;
	double h_elem_raio, h_pcl_ratio;
};

#endif