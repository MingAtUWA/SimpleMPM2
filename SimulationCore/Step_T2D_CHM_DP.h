#ifndef __Step_T2D_CHM_DP_H__
#define __Step_T2D_CHM_DP_H__

#include "Step.h"
#include "Model_T2D_CHM_DP.h"

int solve_substep_T2D_CHM_DP(void *_self);

// for single object only
class Step_T2D_CHM_DP : public Step
{
protected:
	typedef Model_T2D_CHM_DP::Node Node;
	typedef Model_T2D_CHM_DP::Element Element;
	typedef Model_T2D_CHM_DP::SolidParticle SolidParticle;
	typedef Model_T2D_CHM_DP::FluidParticle FluidParticle;

	int init_calculation(void) override;
	friend int solve_substep_T2D_CHM_DP(void *_self);
	int finalize_calculation(void) override;

public:
	Step_T2D_CHM_DP();
	~Step_T2D_CHM_DP();

	inline void set_model(Model_T2D_CHM_DP &md) { model = &md; Step::set_model(md); }
	inline void set_prev_step(Step_T2D_CHM_DP &prev_step) { model = prev_step.model; Step::set_prev_step(prev_step); }

	inline void set_dtime(
		double _dt,
		double dt_max_min_raio = 0.1,
		double t_tol_r = 0.01
		)
	{
		max_dt = _dt;
		min_dt = max_dt * dt_max_min_raio;
		time_tol_ratio = t_tol_r;
		time_tol = min_dt * t_tol_r;
		dtime = max_dt;
	}

protected:
	Model_T2D_CHM_DP *model;
	double min_dt, max_dt;
};

#endif