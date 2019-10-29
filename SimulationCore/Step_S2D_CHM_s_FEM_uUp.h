#ifndef __Step_S2D_CHM_s_FEM_uUp_H__
#define __Step_S2D_CHM_s_FEM_uUp_H__

#include <fstream>

#include "MatrixCoefficientSet.hpp"

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
		time_tol = min_dt * time_tol_ratio;
		dtime = max_dt;
	}

protected:
	Model_S2D_CHM_s_FEM_uUp *model;
	double min_dt, max_dt;
	MatrixCoefficientSet<> g_kmat_coefs;
	double *kmat_col;
	
protected:
	void form_elem_stiffness_mat_and_force_vec(
		Model_S2D_CHM_s_FEM_uUp::Element &e,
		double kmat[20][20], double fvec[20]);
	
	void update_gauss_point(
		Model_S2D_CHM_s_FEM_uUp::GaussPoint &gp,
		Model_S2D_CHM_s_FEM_uUp::ShapeFuncValue &sf,
		double dt, double dt2,
		double dux_s1, double dux_s2, double dux_s3, double dux_s4,
		double duy_s1, double duy_s2, double duy_s3, double duy_s4,
		double dux_f1, double dux_f2, double dux_f3, double dux_f4,
		double duy_f1, double duy_f2, double duy_f3, double duy_f4,
		double dp1, double dp2, double dp3, double dp4);

public: // for debugging
	std::fstream out_file;
};

#endif

#ifdef Keep_Newmark_Coefficients
#define beta 0.3025
#define gamma 0.6
#endif