#ifndef __Generate_Animation_H__
#define __Generate_Animation_H__

#include <fstream>

#include "ItemArray.hpp"

#include "SquareParticleSystem.h"
#include "ResultFile_PlainBin_DataStruct.h"
#include "ShaderProgram.h"

class GenerateAnimation
{
public:
	static GenerateAnimation *cur_generator;

protected:
	// shader program
	ShaderProgram shader;
	// color
	GLint color_id;
	// model/view matrix id in shader
	GLint mv_mat_id;
	// projection matrix id in shader
	GLint proj_mat_id;

	// plain buffer
	std::fstream res_file;

	struct BufferObject
	{
		GLuint vao, vbo, veo;
		BufferObject() : vao(0), vbo(0), veo(0) {}
		~BufferObject() { clear(); }
		void init_array_buffer(GLfloat *coords, size_t coords_len)
		{
			glGenVertexArrays(1, &vao);
			glBindVertexArray(vao);
			glGenBuffers(1, &vbo);
			glBindBuffer(GL_ARRAY_BUFFER, vbo);
			glBufferData(GL_ARRAY_BUFFER, coords_len * sizeof(GLfloat), coords, GL_STATIC_DRAW);
		}
		// after init_araray_buffer()
		void init_elem_array_buffer(GLuint *indices, size_t indices_num)
		{
			glBindVertexArray(vao);
			glGenBuffers(1, &veo);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, veo);
			glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices_num * sizeof(GLuint), indices, GL_STATIC_DRAW);
		}
		void clear(void)
		{
			if (vao) glDeleteVertexArrays(1, &vao);
			if (vbo) glDeleteBuffers(1, &vbo);
			if (veo) glDeleteBuffers(1, &veo);
		}
		// additional element array()
		void init_add_elem_array_buffer(GLuint &add_veo, GLuint *indices, size_t indices_num)
		{
			glBindVertexArray(vao);
			glGenBuffers(1, &add_veo);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, add_veo);
			glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices_num * sizeof(GLuint), indices, GL_STATIC_DRAW);
		}
		void clear_add_buffer(GLuint &add_vbo)
		{
			if (add_vbo) glDeleteBuffers(1, &add_vbo);
		}
		inline void use(void) const { glBindVertexArray(vao); }
	};
	BufferObject bg_grid_data;
	BufferObject rigid_body_data;
	BufferObject mp_data;

protected: // time history data
	typedef ResultFile_PlainBin_DataStruct::TimeHistoryHeader TimeHistoryHeader;
	TimeHistoryHeader *time_rcds;
	size_t time_rcd_num;

public: // animation generation
	double real_time;
	double animation_time;
	double ani_real_ratio;
	double min_delay_real; // minimum delay (in real time)

	// index of current frame
	size_t cur_time_rcd_id;
	// delay time between current frame and next frame
	double delay_ani;
	unsigned short int delay_ani_100th;

	// window
	GLsizei win_width, win_height;
	// viewport
	double vp_hw_ratio, vp_wh_ratio;

public:
	GenerateAnimation(GLsizei win_w = 600, GLsizei win_h = 600);
	virtual ~GenerateAnimation();
	inline void set_init_win_size(GLsizei _w, GLsizei _h) noexcept { win_width = _w, win_height = _h; }
	inline void make_current(void) noexcept { cur_generator = this; }
	// main function
	int generate(double ani_time, double xl, double xu, double yl, double yu,
				 const char *res_file_name, const char *gif_name = nullptr);

protected: // helper functions for generate()
	virtual int init(const char *res_file_name) { return 0; }
	virtual int render_frame(double xl, double xu, double yl, double yu) { return 0; }
	// return true if there is still next frame
	// return false if already reach the last frame
	bool find_next_frame(void);
};

#endif