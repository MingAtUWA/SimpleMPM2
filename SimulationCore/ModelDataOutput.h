#ifndef __MODEL_DATA_OUTPUT_H__
#define __MODEL_DATA_OUTPUT_H__

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
	double time;

	Model *model;
	ResultFile *res_file;
	ModelDataOutputFunc output_func;

public:
	ModelDataOutput(const char *_type = "ModelDataOutput",
		ModelDataOutputFunc _output_func = &model_data_output_func_null) :
		type(_type), name(20, '\0'), time(0.0),
		model(nullptr), res_file(nullptr), output_func(_output_func),
		next(nullptr) {}
	~ModelDataOutput() {}
	inline void set_output_time(double _time) noexcept { time = _time; }
	inline double get_output_time(void) const noexcept { return time; }
	inline void set_model(Model &_model) noexcept { model = &_model; }
	inline Model &get_model(void) noexcept { return *model; }
	inline int output(void) { return (*output_func)(*this); }

protected:
	ModelDataOutput *next;
};

#endif