#ifndef __SHADER_PROGRAM_H__
#define __SHADER_PROGRAM_H__

#include <string>
#include <map>

#include "OpenGL_headers.h"

#include "ItemArray.hpp"

class ShaderProgram
{
protected:
	GLuint vs_id, fs_id, program_id;
	// Uniform variables
	typedef std::pair<std::string, GLint> UniformMapPair;
	typedef std::map<std::string, GLint> UniformMap;
	UniformMap uniform_map;

public:
	ShaderProgram() : vs_id(0), fs_id(0), program_id(0) {}
	~ShaderProgram() { clear(); }
	void init_from_code(const char *vs_code, const char *fs_code);
	void init_from_file(const char *vs_file_name, const char *fs_file_name);
	inline void use(void) { glUseProgram(program_id); }
	inline void unuse(void) { glUseProgram(0); }
	void clear(void);

	// Uniform variables operations
	// Need to first call use() before calling function below
	int init_uniform(const char *uniform_name);
	int set_uniform_vec4f(GLint uniform_loc, const GLfloat *mat_value);
	int set_uniform_matrix4f(GLint uniform_loc, const GLfloat *mat_value);
	int set_uniform_matrix4f(const char *uniform_name, const GLfloat *mat_value);
	// need to first init_uniform() before get_uniform_loc()
	GLint get_uniform_loc(const char * uniform_name);

protected:
	void compile_vertex_shader(const char *shader_source);
	void compile_fragment_shader(const char *shader_source);
	void link_shader_program(void);
};

#endif