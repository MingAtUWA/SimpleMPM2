#ifndef __DISPLAY_MODEL_H__
#define __DISPLAY_MODEL_H__

#include "ShaderProgram.h"
#include "Geometry.h"
#include "RigidBody.h"

// Usage:
// disp_model.init_win();
// disp_model.init_rigid_body(rigid_body);
// disp_model.init_bg_mesh(h, elem_x_num, elem_y_num, x0, y0);
// disp_model.init_particles(pcls, pcl_num);
// disp_model.display(left, right, bottom, top);
class DisplayModel
{
protected:
	ShaderProgram shader;
	// rigid body
	GLuint rb_vao;
	GLuint rb_vbo;
	GLuint rb_ibo;
	GLsizei rb_id_num;
	glm::mat4 rb_mv_mat;
	// background mesh
	GLuint bgm_vao;
	GLuint bgm_vbo;
	GLsizei grid_line_points_num;
	glm::mat4 bgm_mv_mat;
	// particles
	GLuint pcls_vao;
	GLuint pcls_vbo;
	GLsizei pcls_id_num;
	glm::mat4 pcls_mv_mat;
	// shader program
	GLint color_id;
	GLint mv_mat_id;
	GLint proj_mat_id;

	GLFWwindow *window;

public:
	DisplayModel();
	~DisplayModel();

	int init_win(void);
	int init_rigid_body(RigidBody &rb);
	int init_bg_mesh(double grid_size,
					 size_t _elem_x_num, size_t _elem_y_num,
					 double x_start = 0.0, double y_start = 0.0);

	template <typename ParticleType>
	int init_particles(ParticleType *pcls, size_t pcl_num)
	{
		pcls_mv_mat = glm::mat4(1.0f);

		pcls_id_num = pcl_num;
		glGenVertexArrays(1, &pcls_vao);
		glBindVertexArray(pcls_vao);
		glGenBuffers(1, &pcls_vbo);
		glBindBuffer(GL_ARRAY_BUFFER, pcls_vbo);
		GLfloat *pcl_coords_data = new GLfloat[pcl_num * 3];
		for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
		{
			pcl_coords_data[pcl_id * 3]     = GLfloat(pcls[pcl_id].x);
			pcl_coords_data[pcl_id * 3 + 1] = GLfloat(pcls[pcl_id].y);
			pcl_coords_data[pcl_id * 3 + 2] = 0.0f;
		}
		glBufferData(GL_ARRAY_BUFFER, pcls_id_num * 3 * sizeof(GLfloat), pcl_coords_data, GL_STATIC_DRAW);
		delete[] pcl_coords_data;
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 3, (GLvoid *)0);
		glEnableVertexAttribArray(0);

		return 0;
	}

	int display(double left, double right, double bottom, double top);

protected:
	void render(void);
};

#endif