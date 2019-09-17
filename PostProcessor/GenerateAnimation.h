#pragma once

#include <fstream>
#include "ItemArray.hpp"

#include "ShaderProgram.h"

class SquareParticleSystem
{
protected:
	MemoryUtilities::ItemArray<GLfloat> coords;
	MemoryUtilities::ItemArray<GLuint> indices;
	MemoryUtilities::ItemArray<GLuint> bl_indices;
	size_t index_num;

public:
	SquareParticleSystem() : index_num(0)
	{
		coords.set_page_size(30);
		indices.set_page_size(30);
		bl_indices.set_page_size(30);
	}
	~SquareParticleSystem() { clear(); }
	void add_pcl(double x, double y, double vol)
	{
		double half_len = sqrt(vol) / 2.0;
		GLfloat *new_coords = coords.alloc(12);
		new_coords[0] = GLfloat(x - half_len);
		new_coords[1] = GLfloat(y - half_len);
		new_coords[2] = 0.0f;
		new_coords[3] = GLfloat(x + half_len);
		new_coords[4] = GLfloat(y - half_len);
		new_coords[5] = 0.0f;
		new_coords[6] = GLfloat(x - half_len);
		new_coords[7] = GLfloat(y + half_len);
		new_coords[8] = 0.0f;
		new_coords[9] = GLfloat(x + half_len);
		new_coords[10] = GLfloat(y + half_len);
		new_coords[11] = 0.0f;
		GLuint *new_indices = indices.alloc(6);
		new_indices[0] = GLuint(index_num);
		new_indices[1] = GLuint(index_num + 1);
		new_indices[2] = GLuint(index_num + 2);
		new_indices[3] = GLuint(index_num + 1);
		new_indices[4] = GLuint(index_num + 3);
		new_indices[5] = GLuint(index_num + 2);
		GLuint *new_bl_indices = bl_indices.alloc(8);
		new_bl_indices[0] = GLuint(index_num);
		new_bl_indices[1] = GLuint(index_num + 1);
		new_bl_indices[2] = GLuint(index_num + 1);
		new_bl_indices[3] = GLuint(index_num + 3);
		new_bl_indices[4] = GLuint(index_num + 3);
		new_bl_indices[5] = GLuint(index_num + 2);
		new_bl_indices[6] = GLuint(index_num + 2);
		new_bl_indices[7] = GLuint(index_num);
		index_num += 4;
	}
	inline GLfloat *get_pcls(void) const { return coords.get_mem(); }
	inline size_t get_point_num(void) const { return coords.get_num(); }
	inline GLuint *get_indices(void) const { return indices.get_mem(); }
	inline GLsizei get_index_num(void) const { return GLsizei(indices.get_num()); }
	inline GLuint *get_bl_indices(void) const { return bl_indices.get_mem(); }
	inline GLsizei get_bl_index_num(void) const { return GLsizei(bl_indices.get_num()); }
	void reset(void)
	{
		coords.reset();
		indices.reset();
		bl_indices.reset();
		index_num = 0;
	}
	void clear(void)
	{
		coords.clear();
		indices.clear();
		bl_indices.clear();
		index_num = 0;
	}
};

class GenerateAnimation
{
public:
	// shader program
	ShaderProgram shader;
	// color
	GLint color_id;
	// model/view matrix id in shader
	GLint mv_mat_id;
	// projection matrix id in shader
	GLint proj_mat_id;

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

	// plain buffer
	std::fstream res_file;

public: // model data
	struct MeshHeader
	{
		double h, x0, xn, y0, yn;
		unsigned long long elem_x_num, elem_y_num;
	};
	struct RigidBodyHeader
	{
		unsigned long long node_num, elem_num;
		double x_mc, y_mc;
		// node coordinates, elem indices
	};
	struct MPObjectHeader
	{
		unsigned long long pcl_num;
		// x coord, y coord, vol, ...
	};

	MeshHeader mh;
	GLsizei grid_line_points_num;
	RigidBodyHeader rbh;
	GLsizei rb_elem_point_num;
	MPObjectHeader mph;

public: // time history data
	struct TimeHistoryHeader
	{
		size_t substep_num;
		size_t total_substep_num;
		double current_time;
		double total_time;
	};
	struct RigidBodyMotionHeader
	{
		double x, y, theta;
		double vx, vy, v_theta;
	};

	// particle data
	double *mp_x_data, *mp_y_data, *mp_vol_data;
	SquareParticleSystem pcls_mem;

	// start position of current time record
	size_t first_time_rcd_file_pos;
	// length of each time record
	size_t time_rcd_len;
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

public:
	GenerateAnimation();
	~GenerateAnimation();
	// main function
	int generate(double ani_time, double xl, double xu, double yl, double yu,
				 const char *res_file_name, const char *gif_name = nullptr);

protected: // helper functions of generate()
	// assume pcl_num is the same in each time history output
	int init(const char *res_file_name);
	int render_frame(double xl, double xu, double yl, double yu);
	// return true if there is still next frame
	// return false if already reach the last frame
	bool find_next_frame(void);
};