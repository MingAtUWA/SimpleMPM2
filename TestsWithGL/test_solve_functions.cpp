#include "TestsWithGL_pcp.h"

#include "Step.h"
#include "ModelDataOutput.h"
#include "TimeHistoryOutput.h"
#include "TimeHistoryOutput_ConsoleProgressBar.h"

#include "test_sim_core.h"

#include <Windows.h>

int solve_substep_test(void *_self);
class Step_test : public Step
{
public:
	friend int solve_substep_test(void *_self);
	Step_test() : Step(&solve_substep_test, "Step_test") {}
};
int solve_substep_test(void *_self)
{
	Step_test &self = *reinterpret_cast<Step_test *>(_self);
	//std::cout << "Substep " << self.substep_num
	//		  << " time: " << self.current_time
	//		  << " dt: " << self.dtime << "\n";
	Sleep(500);
	return 0;
}

int model_data_output_func_test(ModelDataOutput &_self);
class ModelDataOutput_test : public ModelDataOutput
{
	static size_t cur_id;
	size_t id;
public:
	friend int model_data_output_func_test(ModelDataOutput &_self);
	ModelDataOutput_test(const char *_name) :
		ModelDataOutput(_name, "ModelDataOutput_test", &model_data_output_func_test), id(cur_id++) {}
};
size_t ModelDataOutput_test::cur_id = 0;
int model_data_output_func_test(ModelDataOutput &_self)
{
	ModelDataOutput_test &self = static_cast<ModelDataOutput_test &>(_self);
	size_t md_num = 1;
	for (ModelDataOutput *pmd = self.next; pmd; pmd = ((ModelDataOutput_test *)pmd)->next)
		++md_num;
	std::cout << "ModelData " << self.id << " time: " << self.current_time
			  << " md_num: " << md_num << "\n"; return 0;
}

int time_history_output_func_test(TimeHistoryOutput &_self);
class TimeHistoryOutput_test : public TimeHistoryOutput
{
	static size_t cur_id;
	size_t id;
public:
	friend int time_history_output_func_test(TimeHistoryOutput &_self);
	TimeHistoryOutput_test(const char *_name) :
		TimeHistoryOutput(_name, "TimeHistoryOutput_test", &time_history_output_func_test), id(cur_id++) {}
};
size_t TimeHistoryOutput_test::cur_id = 0;
int time_history_output_func_test(TimeHistoryOutput &_self)
{
	TimeHistoryOutput_test &self = static_cast<TimeHistoryOutput_test &>(_self);
	std::cout << "TimeHistory " << self.id << " time: " << self.get_step().get_current_time() << "\n"; return 0;
}

void test_solve_functions(void)
{
	ModelDataOutput_test md1("md1");
	md1.set_output_time(0.1);
	
	ModelDataOutput_test md2("md2");
	md2.set_output_time(1.0);

	TimeHistoryOutput_test th1("th1");
	th1.set_interval_num(5);
	th1.set_output_init_state();

	TimeHistoryOutput_test th2("th2");
	th2.set_interval_num(2);
	
	TimeHistoryOutput_ConsoleProgressBar cpb;

	Step_test step;
	step.set_time(1.0);
	step.set_dtime(0.1);
	//step.use_solve_th_and_md();
	//step.use_solve_th_only();

	//step.add_time_history(th1);
	//step.add_time_history(th2);
	//step.add_model_data(md1);
	//step.add_model_data(md2);
	
	step.add_time_history(cpb);

	step.solve();

	system("pause");
}
