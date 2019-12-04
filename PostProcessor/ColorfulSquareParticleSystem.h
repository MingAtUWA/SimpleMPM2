#ifndef __Colorful_Square_Particle_System_H__
#define __Colorful_Square_Particle_System_H__

#include "ItemArray.hpp"

class ColorfulSquareParticleSystem
{
protected:
	MemoryUtilities::ItemArray<GLfloat> coord_and_color;
	MemoryUtilities::ItemArray<GLuint> indices;
	//MemoryUtilities::ItemArray<GLuint> bl_indices;
	size_t index_num;

public:
	ColorfulSquareParticleSystem() : index_num(0)
	{
		coord_and_color.set_page_size(60);
		indices.set_page_size(30);
		//bl_indices.set_page_size(30);
	}
	~ColorfulSquareParticleSystem() { clear(); }
	inline void set_pcl_num(size_t pcl_num)
	{
		coord_and_color.reserve(pcl_num * 24);
		indices.set_page_size(pcl_num * 6);
		//bl_indices.set_page_size(pcl_num * 8);
	}
	void add_pcl(double x, double y, double vol, float r, float g, float b)
	{
		double half_len = sqrt(vol) / 2.0;
		// coordinates
		GLfloat *pcl_coord_and_color = coord_and_color.alloc(24);
		// pcl 1
		pcl_coord_and_color[0] = GLfloat(x - half_len);
		pcl_coord_and_color[1] = GLfloat(y - half_len);
		pcl_coord_and_color[2] = 0.0f;
		pcl_coord_and_color[3] = GLfloat(r);
		pcl_coord_and_color[4] = GLfloat(g);
		pcl_coord_and_color[5] = GLfloat(b);
		// pcl 2
		pcl_coord_and_color[6] = GLfloat(x + half_len);
		pcl_coord_and_color[7] = GLfloat(y - half_len);
		pcl_coord_and_color[8] = 0.0f;
		pcl_coord_and_color[9] = GLfloat(r);
		pcl_coord_and_color[10] = GLfloat(g);
		pcl_coord_and_color[11] = GLfloat(b);
		// pcl 3
		pcl_coord_and_color[12] = GLfloat(x - half_len);
		pcl_coord_and_color[13] = GLfloat(y + half_len);
		pcl_coord_and_color[14] = 0.0f;
		pcl_coord_and_color[15] = GLfloat(r);
		pcl_coord_and_color[16] = GLfloat(g);
		pcl_coord_and_color[17] = GLfloat(b);
		// pcl 4
		pcl_coord_and_color[18] = GLfloat(x + half_len);
		pcl_coord_and_color[19] = GLfloat(y + half_len);
		pcl_coord_and_color[20] = 0.0f;
		pcl_coord_and_color[21] = GLfloat(r);
		pcl_coord_and_color[22] = GLfloat(g);
		pcl_coord_and_color[23] = GLfloat(b);
		// indices
		GLuint *new_indices = indices.alloc(6);
		new_indices[0] = GLuint(index_num);
		new_indices[1] = GLuint(index_num + 1);
		new_indices[2] = GLuint(index_num + 2);
		new_indices[3] = GLuint(index_num + 1);
		new_indices[4] = GLuint(index_num + 3);
		new_indices[5] = GLuint(index_num + 2);
		//GLuint *new_bl_indices = bl_indices.alloc(8);
		//new_bl_indices[0] = GLuint(index_num);
		//new_bl_indices[1] = GLuint(index_num + 1);
		//new_bl_indices[2] = GLuint(index_num + 1);
		//new_bl_indices[3] = GLuint(index_num + 3);
		//new_bl_indices[4] = GLuint(index_num + 3);
		//new_bl_indices[5] = GLuint(index_num + 2);
		//new_bl_indices[6] = GLuint(index_num + 2);
		//new_bl_indices[7] = GLuint(index_num);
		index_num += 4;
	}
	inline GLfloat *get_coord_and_color(void) const { return coord_and_color.get_mem(); }
	inline size_t get_coord_and_color_size(void) const { return coord_and_color.get_num(); }
	inline GLuint *get_indices(void) const { return indices.get_mem(); }
	inline GLsizei get_index_num(void) const { return GLsizei(indices.get_num()); }
	//inline GLuint *get_bl_indices(void) const { return bl_indices.get_mem(); }
	//inline GLsizei get_bl_index_num(void) const { return GLsizei(bl_indices.get_num()); }
	inline void reset(void)
	{
		coord_and_color.reset();
		indices.reset();
		//bl_indices.reset();
		index_num = 0;
	}
	inline void clear(void)
	{
		coord_and_color.clear();
		indices.clear();
		//bl_indices.clear();
		index_num = 0;
	}
};

#endif