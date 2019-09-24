#include "SimulationCore_pcp.h"

#include "Step.h"

Step::Step(SolveSubstepFunc solve_substep_func) :
	solve_substep(solve_substep_func),
	name(20, '\0'), model(nullptr),
	is_first_step(true), start_substep_index(0), substep_num(0),
	step_time(0.0), start_time(0.0), current_time(0.0),
	dtime(0.0), time_tol_ratio(0.01), time_tol(0.0),
	time_history_top(nullptr) {}

Step::~Step() {}

int Step::solve(void)
{
	substep_num = 0;
	current_time = 0.0;
	
	// initialize calculation
	init_calculation();

	// init time history
	for (TimeHistory *pth = time_history_top; pth; pth = pth->next)
	{
		if (pth->need_output_init_state) pth->output();
		pth->interval_time = step_time / double(pth->interval_num);
		pth->next_time = pth->interval_time;
		pth->init_per_step();
	}

	double time_diff_tmp;
	do
	{
		// solve substep
		(*solve_substep)(this);

		++substep_num;
		current_time += dtime;
		time_diff_tmp = current_time - step_time;
		if (time_diff_tmp > time_tol)
		{
			dtime -= time_diff_tmp;
			current_time = step_time;
		}

		// output time history if needed
		output_time_history();

	} while (-time_diff_tmp > time_tol);

	// finalize time history
	for (TimeHistory *pth = time_history_top; pth; pth = pth->next)
		pth->finalize_per_step();

	// finalize calculation
	finalize_calculation();

	return 0;
}

int solve_substep_base(void *_self) { return 0; }

void Step::output_all_time_history(void)
{
	for (TimeHistory *pth = time_history_top; pth; pth = pth->next)
		pth->output();
}

void Step::output_time_history(void)
{
	for (TimeHistory *pth = time_history_top; pth; pth = pth->next)
	{
		if (pth->next_time <= current_time + time_tol)
		{
			pth->output();
			pth->next_time += pth->interval_time;
			if (pth->next_time < current_time)
				pth->next_time = current_time + pth->interval_time;
		}
	}
}
