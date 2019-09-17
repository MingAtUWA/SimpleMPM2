#ifndef _STEP_H_
#define _STEP_H_

#include <string>
#include "Model.h"
#include "TimeHistory.h"

// Avoid using virtual function
typedef int(*SolveSubstepFunc)(void *_self);

/* ========================================================
Class Step:
	Functions needs to be rewritten in children classes:
	1. init()
	2. solve_substep()
	3. finalize()
 =========================================================== */
int solve_substep_base(void *_self);

class Step
{
protected:
	const char *type;
	std::string name;

	bool is_first_step;
	// substep index at start of this step
	size_t start_substep_index;
	// number of substep from the start of this step
	size_t substep_num;
	// time length of this step
	double step_time;
	// start time for this step
	double start_time;
	// time from the start of this step
	double current_time;
	// time increment
	double dtime; // time increment
	double time_tol_ratio;
	double time_tol; // = dt * time_tol_ratio
	// model
	Model *model;

public:
	Step(SolveSubstepFunc solve_substep_func = &solve_substep_base);
	~Step();
	// type
	inline const char *get_type(void) const { return type; }
	// step name
	inline const char *get_name(void) const noexcept { return name.c_str(); }
	inline void set_name(const char *_name) noexcept { name = _name; }
	// step time length
	inline void set_time(double _time) noexcept { step_time = _time; }
	inline void set_dtime(double _dtime, double t_tol_r = 0.01) noexcept
	{
		dtime = _dtime, time_tol_ratio = t_tol_r, time_tol = dtime * t_tol_r;
	}
	// time from the start of this step
	inline double get_current_time(void) { return current_time; }
	// total time from the start of the whole simulation
	inline double get_total_time(void) { return start_time + current_time; }
	// number of substep from teh start of this step
	inline size_t get_substep_num(void) { return substep_num; }
	// total number of substep from the start of the whole simulation
	inline size_t get_total_substep_num(void) { return start_substep_index + substep_num; }
	// time length of this step
	inline double get_step_time(void) { return step_time; }
	// size of time increment
	inline double get_dtime(void) { return dtime; }

	// set model
	inline void set_model(Model &md) noexcept { model = &md; }
	inline Model &get_model(void) const noexcept { return *model; }
	// continuate from prev step
	void set_prev_step(Step &prev_step)
	{
		model = prev_step.model;
		is_first_step = false;
		start_substep_index = prev_step.get_total_substep_num();
		start_time = prev_step.get_total_time();
	}

protected:
	// initialization before calculation
	virtual int init_calculation(void) { return 0; }
	// calculation of each substep
	SolveSubstepFunc solve_substep;
	// finalization after calculation
	virtual int finalize_calculation(void) { return 0; }

public: // main functions
	virtual int solve(void);
	
	// Time history Utilities
public:
	inline void add_time_history(TimeHistory &th) noexcept
	{
		th.step = this;
		th.model = model;
		th.next = time_history_top;
		time_history_top = &th;
	}
	inline void clear_time_history(void) noexcept
	{
		time_history_top = nullptr;
	}

protected:
	TimeHistory *time_history_top;
	void output_all_time_history(void);
	void output_time_history(void);
};

#endif