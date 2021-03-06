#ifndef __Display_Model_S2D_H__
#define __Display_Model_S2D_H__

#include "ShaderProgram.h"
#include "Geometry.h"
#include "SquareParticleSystem.h"

#include "Model_S2D_ME_s.h"

// Usage:
// disp_model.init_win();
// disp_model.init_model(model);
// disp_model.add_point(coords, num);
// disp_model.display(left, right, bottom, top);
class DisplayModel_S2D
{
protected:
	ShaderProgram shader;
	GLFWwindow *window;

public:
	DisplayModel_S2D();
	~DisplayModel_S2D();

	int init_win(void);
	int display(double left, double right, double bottom, double top);

	int init_model(Model_S2D_ME_s &md);

	int init_points(GLfloat *n_coords, size_t num); // num = sizeof(n_coords) * 3

protected:
	// shader program
	GLint color_id;
	GLint mv_mat_id;
	GLint proj_mat_id;
	
	// background mesh
	GLuint bgm_vao;
	GLuint bgm_vbo, bgm_ibo;
	GLsizei bgm_elem_n_id_num;
	glm::mat4 bgm_mv_mat;

	// particles
	GLuint pcls_vao;
	GLuint pcls_vbo;
	GLuint pcls_veo;
	GLsizei pcls_id_num;
	glm::mat4 pcls_mv_mat;
	
	// rigid circle
	GLuint rc_vao;
	GLuint rc_vbo;
	GLsizei rc_pcls_id_num;
	glm::mat4 rc_mv_mat;

	// points
	GLuint pt_vao;
	GLuint pt_vbo;
	GLsizei pt_num;
	glm::mat4 pt_mv_mat;

	SquareParticleSystem pcls_mem;

protected:
	int init_bg_mesh(Model_S2D_ME_s &md);
	int init_pcls(Model_S2D_ME_s &md);

protected:
	void render(void);
};

#endif