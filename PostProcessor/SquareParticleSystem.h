#ifndef __Square_Particle_System_H__
#define __Square_Particle_System_H__

#include "ItemArray.hpp"

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
	inline void set_pcl_num(size_t pcl_num)
	{
		coords.reserve(pcl_num * 12);
		indices.set_page_size(pcl_num * 6);
		bl_indices.set_page_size(pcl_num * 8);
	}
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
	inline void reset(void)
	{
		coords.reset();
		indices.reset();
		bl_indices.reset();
		index_num = 0;
	}
	inline void clear(void)
	{
		coords.clear();
		indices.clear();
		bl_indices.clear();
		index_num = 0;
	}
};

#endif