#ifndef __Model_Data_Output_H__
#define __Model_Data_Output_H__

#include <string>

#include "ResultFile.h"

class Model;
class Step;
class ModelDataOutput;

typedef int (*ModelDataOutputFunc)(ModelDataOutput &_self);
int model_data_output_func_null(ModelDataOutput &_self);

/*=============================================================
Class ModelDataOutput
==============================================================*/
class ModelDataOutput
{
	friend Step;

protected:
	const char *type;
	std::string name;
	double current_time;
	double total_time;

	Model *model;
	ResultFile *res_file;
	ModelDataOutputFunc output_func;

public:
	ModelDataOutput(const char *_name,
		const char *_type = "ModelDataOutput",
		ModelDataOutputFunc _output_func = &model_data_output_func_null) :
		name(_name), type(_type), current_time(0.0), total_time(0.0),
		model(nullptr), res_file(nullptr), next(nullptr),
		output_func(_output_func) {}
	~ModelDataOutput() {}

	// output time is respect to start of this step
	inline void set_output_time(double _time) noexcept { current_time = _time; }
	inline double get_output_time(void) const noexcept { return current_time; }
	inline void set_model(Model &_model) noexcept { model = &_model; }
	inline Model &get_model(void) noexcept { return *model; }
	inline void set_name(const char *_name) { name = _name; }
	inline const char *get_name(void) { return name.c_str(); }
	
	inline int output(void) { return (*output_func)(*this); }

protected:
	ModelDataOutput *next;
};

#endif