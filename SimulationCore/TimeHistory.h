#ifndef _TIME_HISTORY_H_
#define _TIME_HISTORY_H_

#include <string>
class Model;
class Step;
#include "ResultFile.h"

/*=============================================================
Class TimeHistory
	Functions needs to be rewritten: output();
 ==============================================================*/
class TimeHistory
{
	friend Step;
protected:
	const char *type;
	std::string name;
	// number of output this step:
	size_t interval_num;
	double interval_time;
	// output schedule
	double next_time;
	bool need_output_init_state; // true if output the initial state
	
	ResultFile *res_file;
	OutputFunc output_func;

public:
	TimeHistory(const char *_type = "TimeHistory") :
		type(_type), name(20, '\0'),
		need_output_init_state(false), interval_num(1),
		step(nullptr), next(nullptr) {}
	~TimeHistory() {}
	
	inline const char *get_type(void) const { return type; }
	inline void set_name(const char *_name) { name = _name; }
	inline const char *get_name(void) const { return name.c_str(); }
	inline void set_interval_num(size_t num) { interval_num = num; }
	inline size_t get_interval_num(void) const { return interval_num; }
	inline void set_output_init_state(bool _need = true) noexcept { need_output_init_state = _need; }
	inline Model &get_model(void) const noexcept { return *model; }
	inline Step &get_step(void) const noexcept { return *step; }

	// Initialize each steps
	virtual int init_per_step(void) { return 0; }
	// Finalize each steps
	virtual void finalize_per_step(void) {}

	// Output per substep
	virtual int output(void) { return 0; }

protected: // set and used by Step class
	Step *step;
	Model *model;
	TimeHistory *next; // used by step
};

#endif