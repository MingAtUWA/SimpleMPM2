#pragma once

#include "OpenGL_headers.h"

#include "ShaderProgram.h"
#include "TriangleMesh.h"

// Draw triangular mesh
class DrawTriangleMesh
{
protected:
	// For drawing elements
	GLuint va_elem;
	GLuint nb_elem; // nodal buffer: coordinates / colors
	GLuint ib_elem; // index buffer
	// For drawing boundary lines
	GLuint va_bl;
	GLuint nb_bl;   // nodal buffer: coordinates / colors
	GLuint ib_bl;   // index buffer
	// For drawing back-ground grid
	GLuint va_bg_grid;
	GLuint nb_bg_grid;
	GLuint va_grid_lines;
	GLuint nb_grid_line_coords;

	GLsizei node_num;
	GLsizei elem_num;
	GLsizei line_num;
	GLsizei grid_num;
	GLsizei grid_line_num;

protected:
	ShaderProgram unicolor_shader;
	GLint uc_mv_mat_id, uc_proj_mat_id, uc_color_vec_id;
	ShaderProgram multicolor_shader;
	GLint mc_mv_mat_id, mc_proj_mat_id;

	double rb_x_mc, rb_y_mc;

public:
	DrawTriangleMesh();
	~DrawTriangleMesh();
	void init(TriangleMesh &tri_mesh,
			  const char *vs_name, const char *fs_name,
			  const char *bg_grid_vs_name, const char *bg_grid_fs_name);
	void draw(bool disp_tri_mesh = true, bool dis_bl = true, bool disp_bg_grid = true);

protected:
	void init_tri_mesh(TriangleMesh &tri_mesh);
	void init_bl(TriangleMesh &tri_mesh);
	void init_bg_grid(TriangleMesh &tri_mesh);
	void init_grid_lines(TriangleMesh &tri_mesh);

	void draw_tri_mesh(void);
	void draw_bl(void);
	void draw_bg_grid(void);
	void draw_grid_lines(void);

protected:
	GLuint va_line_and_point;
	GLuint nb_line_and_point, ib_line;

public:
	void init_line_and_point(TriangleMesh &tri_mesh, TriangleMesh::Edge &edge, Point &p);
	void draw_line_and_point(void);
};
