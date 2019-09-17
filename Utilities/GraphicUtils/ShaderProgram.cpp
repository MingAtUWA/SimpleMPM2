#include "Utilities_pcp.h"

#include <fstream>
#include "ShaderProgram.h"

void ShaderProgram::init_from_code(const char *vs_code, const char *fs_code)
{
	clear();
	compile_vertex_shader(vs_code);
	compile_fragment_shader(fs_code);
	link_shader_program();
}

void ShaderProgram::init_from_file(const char *vs_file_name, const char *fs_file_name)
{
	clear();
	std::ifstream vs_file(vs_file_name, std::ios::binary);
	std::ifstream fs_file(fs_file_name, std::ios::binary);
	if (!vs_file.is_open())
		throw std::exception("Error opening the vertex shader file.");
	if (!fs_file.is_open())
		throw std::exception("Error opening the fragment shader file.");

	size_t shader_file_len;
	char *shader_code;

	vs_file.seekg(0, SEEK_END);
	shader_file_len = size_t(vs_file.tellg());
	shader_code = new char[shader_file_len+1];
	vs_file.seekg(0, SEEK_SET);
	vs_file.read(shader_code, shader_file_len);
	shader_code[shader_file_len] = '\0';
	//std::cout << shader_code << "\n";
	compile_vertex_shader(shader_code);
	delete[] shader_code;

	fs_file.seekg(0, SEEK_END);
	shader_file_len = size_t(fs_file.tellg());
	shader_code = new char[shader_file_len+1];
	fs_file.seekg(0, SEEK_SET);
	fs_file.read(shader_code, shader_file_len);
	shader_code[shader_file_len] = '\0';
	//std::cout << shader_code << "\n";
	compile_fragment_shader(shader_code);
	delete[] shader_code;

	link_shader_program();	
}

void ShaderProgram::clear(void)
{
	if (program_id)
	{
		glDeleteProgram(program_id);
		program_id = 0;
	}
}

int ShaderProgram::init_uniform(const char *uniform_name)
{
	GLint uniform_loc = glGetUniformLocation(program_id, uniform_name);
	if (uniform_loc >= 0)
	{
		uniform_map.insert(UniformMapPair(uniform_name, uniform_loc));
		return uniform_loc;
	}
	return -1;
}

int ShaderProgram::set_uniform_vec4f(GLint uniform_loc, const GLfloat *mat_value)
{
	glUniform4fv(uniform_loc, 1, mat_value);
	return 0;
}

int ShaderProgram::set_uniform_matrix4f(GLint uniform_loc, const GLfloat *mat_value)
{
	glUniformMatrix4fv(uniform_loc, 1, GL_FALSE, mat_value);
	return 0;
}

int ShaderProgram::set_uniform_matrix4f(const char *uniform_name, const GLfloat *mat_value)
{
	GLint uniform_loc = get_uniform_loc(uniform_name);
	if (uniform_loc >= 0)
	{
		glUniformMatrix4fv(uniform_loc, 1, GL_FALSE, mat_value);
	}
	return -1;
}

GLint ShaderProgram::get_uniform_loc(const char *uniform_name)
{
	UniformMap::iterator res = uniform_map.find(std::string(uniform_name));
	if (res != uniform_map.end())
		return res->second;
	return -1;
}


void ShaderProgram::compile_vertex_shader(const char *shader_source)
{
	vs_id = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(vs_id, 1, &shader_source, nullptr);
	glCompileShader(vs_id);
	GLint res_state, res_info_log_len;
	char *res_info_log;
	glGetShaderiv(vs_id, GL_COMPILE_STATUS, &res_state);
	if (!res_state)
	{
		glGetShaderiv(vs_id, GL_INFO_LOG_LENGTH, &res_info_log_len);
		if (res_info_log_len)
		{
			res_info_log = (char *)alloca(res_info_log_len);
			glGetShaderInfoLog(vs_id, res_info_log_len, nullptr, res_info_log);
			std::string error_msg("Error compiling vertex shader: ");
			error_msg += res_info_log;
			error_msg += "\n";
			throw std::exception(error_msg.c_str());
		}
	}
}

void ShaderProgram::compile_fragment_shader(const char *shader_source)
{
	fs_id = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(fs_id, 1, &shader_source, nullptr);
	glCompileShader(fs_id);
	GLint res_state, res_info_log_len;
	char *res_info_log;
	glGetShaderiv(fs_id, GL_COMPILE_STATUS, &res_state);
	if (!res_state)
	{
		glGetShaderiv(fs_id, GL_INFO_LOG_LENGTH, &res_info_log_len);
		if (res_info_log_len)
		{
			res_info_log = (char *)alloca(res_info_log_len);
			glGetShaderInfoLog(fs_id, res_info_log_len, nullptr, res_info_log);
			std::string error_msg("Error compiling fragment shader: ");
			error_msg += res_info_log;
			error_msg += "\n";
			throw std::exception(error_msg.c_str());
		}
	}
}

void ShaderProgram::link_shader_program(void)
{
	program_id = glCreateProgram();
	glAttachShader(program_id, vs_id);
	glAttachShader(program_id, fs_id);
	glLinkProgram(program_id);
	GLint res_state, res_info_log_len;
	char *res_info_log;
	glGetProgramiv(program_id, GL_LINK_STATUS, &res_state);
	if (!res_state)
	{
		glGetProgramiv(program_id, GL_INFO_LOG_LENGTH, &res_info_log_len);
		if (res_info_log_len)
		{
			res_info_log = (char *)alloca(res_info_log_len);
			glGetProgramInfoLog(program_id, res_info_log_len, nullptr, res_info_log);
			std::string error_msg("Error linking program: ");
			error_msg += res_info_log;
			error_msg += "\n";
			throw std::exception(error_msg.c_str());
		}
	}
	glDeleteShader(vs_id);
	vs_id = 0;
	glDeleteShader(fs_id);
	fs_id = 0;
}
