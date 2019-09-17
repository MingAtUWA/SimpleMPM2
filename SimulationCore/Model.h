#ifndef _MODEL_H_
#define _MODEL_H_

#include <string>

struct Model
{
protected:
	std::string name;
	const char *type;

public:
	Model(const char *type_name = "Model") :
		name(20, '\0'), type(type_name) {}
	~Model() {}
	inline void set_name(const char *na) { name = na; }
	inline const char *get_name(void) const { return name.c_str(); }
	inline const char *get_type(void) const { return type; }
};

#endif