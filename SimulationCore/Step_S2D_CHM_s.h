#ifndef __STEP_S2D_CHM_S_H__
#define __STEP_S2D_CHM_S_H__

#include "Step.h"
#include "Model_S2D_CHM_s.h"

// standard MPM
int solve_substep_S2D_CHM_s(void *_self);

// average stress to the element centre:
//   avg_stress = sum(pcl.stress * pcl.vol) / sum(pcl.vol)
int solve_substep_S2D_CHM_s_avg_stress1(void *_self);
// average stress to the element centre:
//   avg_dN_dx_stress = sum(pcl.stress * pcl.dN_dx * pcl.vol) / sum(pcl.dN_dx * pcl.vol)
//   f_int = avg_dN_dx_stress * dN_dx_elem_centre * elem_pcl_vol
//   as accurate as avg_stress1, but less stable and efficient
int solve_substep_S2D_CHM_s_avg_stress2(void *_self);
// average all variables except traction and solve as FEM
// similar to ALE based on operator splitting
// not yet finished, may lost lots of precision but most robust
int solve_substep_S2D_CHM_s_avg_allvars(void *_self);

// "standard" GIMP, instable, no good
int solve_substep_S2D_CHM_s_GIMP(void *_self);
// average stress scheme combined with GIMP
int solve_substep_S2D_CHM_s_avg_stress_GIMP(void *_self);

// for single object only
class Step_S2D_CHM_s : public Step
{
protected:
	int init_calculation(void) override;
	friend int solve_substep_S2D_CHM_s(void *_self);
	friend int solve_substep_S2D_CHM_s_avg_stress1(void *_self);
	friend int solve_substep_S2D_CHM_s_avg_stress2(void *_self);
	friend int solve_substep_S2D_CHM_s_avg_allvars(void *_self);
	friend int solve_substep_S2D_CHM_s_GIMP(void *_self);
	friend int solve_substep_S2D_CHM_s_avg_stress_GIMP(void *_self);
	int finalize_calculation(void) override;

public:
	Step_S2D_CHM_s();
	~Step_S2D_CHM_s();

	inline void use_standard_mpm(void) noexcept { solve_substep = &solve_substep_S2D_CHM_s; }
	inline void use_avg_stress1(void) noexcept { solve_substep = &solve_substep_S2D_CHM_s_avg_stress1; }
	inline void use_avg_stress2(void) noexcept { solve_substep = &solve_substep_S2D_CHM_s_avg_stress2; }
	inline void use_avg_allvars(void) noexcept { solve_substep = &solve_substep_S2D_CHM_s_avg_allvars; }
	inline void use_gimp(void) noexcept { solve_substep = &solve_substep_S2D_CHM_s_GIMP; }
	inline void use_avg_stress_gimp(void) noexcept { solve_substep = &solve_substep_S2D_CHM_s_avg_stress_GIMP; }

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