#include "TestsWithGL_pcp.h"

#include "DrawTriangleMesh.h"

DrawTriangleMesh::DrawTriangleMesh() : 
	va_elem(0), nb_elem(0), ib_elem(0),
	va_bl(0), nb_bl(0), ib_bl(0),
	node_num(0), elem_num(0), line_num(0),
	va_bg_grid(0), nb_bg_grid(0),
	va_grid_lines(0), nb_grid_line_coords(0),
	grid_num(0), grid_line_num(0),
	va_line_and_point(0), nb_line_and_point(0), ib_line(0) {}

DrawTriangleMesh::~DrawTriangleMesh()
{
	// elements
	if (va_elem) { glDeleteVertexArrays(1, &va_elem); va_elem = 0; }
	if (nb_elem) { glDeleteBuffers(1, &nb_elem); nb_elem = 0; }
	if (ib_elem) { glDeleteBuffers(1, &ib_elem); ib_elem = 0; }
	// boundary lines
	if (va_bl) { glDeleteVertexArrays(1, &va_bl); va_bl = 0; }
	if (nb_bl) { glDeleteBuffers(1, &nb_bl); nb_bl = 0; }
	if (ib_bl) { glDeleteBuffers(1, &ib_bl); ib_bl = 0; }
	// background grid
	if (va_bg_grid) { glDeleteVertexArrays(1, &va_bg_grid); va_bg_grid = 0; }
	if (nb_bg_grid) { glDeleteBuffers(1, &nb_bg_grid); nb_bg_grid = 0; }
	if (va_grid_lines) { glDeleteVertexArrays(1, &va_grid_lines); va_grid_lines = 0; }
	if (nb_grid_line_coords) { glDeleteBuffers(1, &nb_grid_line_coords); nb_grid_line_coords = 0; }
	
	// point and line
	if (va_line_and_point) { glDeleteVertexArrays(1, &va_line_and_point); va_line_and_point = 0; }
	if (nb_line_and_point) { glDeleteBuffers(1, &nb_line_and_point); nb_line_and_point = 0; }
	if (ib_line) { glDeleteBuffers(1, &ib_line); ib_line = 0; }
}

void DrawTriangleMesh::init(TriangleMesh &tri_mesh,
	const char *vs_name, const char *fs_name,
	const char *bg_grid_vs_name, const char *bg_grid_fs_name)
{
	// uniform color shader
	unicolor_shader.init_from_file(vs_name, fs_name);
	// init shader uniform variables
	unicolor_shader.use();
	// model/view matrix
	glm::mat4 mv_mat(1.0f);
	uc_mv_mat_id = unicolor_shader.init_uniform("mv_mat");
	unicolor_shader.set_uniform_matrix4f(uc_mv_mat_id, glm::value_ptr(mv_mat));
	// projection matrix
	Rect disp_range = tri_mesh.get_display_range();
	GLfloat range = (GLfloat)(disp_range.xu > disp_range.yu ? disp_range.xu : disp_range.yu);
	GLfloat padding = range * 0.05f;
	glm::mat4 proj_mat = glm::ortho(-padding, range + padding, -padding, range + padding);
	uc_proj_mat_id = unicolor_shader.init_uniform("proj_mat");
	unicolor_shader.set_uniform_matrix4f(uc_proj_mat_id, glm::value_ptr(proj_mat));
	// color
	uc_color_vec_id = unicolor_shader.init_uniform("color_vec");

	// shader for grid of background mesh
	multicolor_shader.init_from_file(bg_grid_vs_name, bg_grid_fs_name);
	// init shader uniform variables
	multicolor_shader.use();
	// model/view matrix
	mc_mv_mat_id = multicolor_shader.init_uniform("mv_mat");
	unicolor_shader.set_uniform_matrix4f(mc_mv_mat_id, glm::value_ptr(mv_mat));
	// projection matrix
	mc_proj_mat_id = multicolor_shader.init_uniform("proj_mat");
	multicolor_shader.set_uniform_matrix4f(mc_proj_mat_id, glm::value_ptr(proj_mat));

	// init vao buffer
	node_num = GLsizei(tri_mesh.node_num);
	elem_num = GLsizei(tri_mesh.elem_num);
	line_num = GLsizei(tri_mesh.boundary_edge_num);

	init_tri_mesh(tri_mesh);
	init_bl(tri_mesh);
	init_bg_grid(tri_mesh);
	init_grid_lines(tri_mesh);
}

void DrawTriangleMesh::draw(bool disp_tri_mesh, bool dis_bl, bool disp_bg_grid)
{
	if (disp_bg_grid)
	{
		draw_bg_grid();
		draw_grid_lines();
	}

	if (disp_tri_mesh)
		draw_tri_mesh();
	if (dis_bl)
		draw_bl();
}

// draw line instead of triangle next time
void DrawTriangleMesh::init_tri_mesh(TriangleMesh &tri_mesh)
{
	rb_x_mc = tri_mesh.x_mc;
	rb_y_mc = tri_mesh.y_mc;
	
	glGenVertexArrays(1, &va_elem);
	glBindVertexArray(va_elem);

	struct NodeInfo { GLfloat x, y, z; };
	NodeInfo *node_info_array = new NodeInfo[node_num];
	for (size_t i = 0; i < node_num; ++i)
	{
		NodeInfo &node_info = node_info_array[i];
		node_info.x = GLfloat(tri_mesh.nodes[i].x);
		node_info.y = GLfloat(tri_mesh.nodes[i].y);
		node_info.z = 0.0f;
	}
	glGenBuffers(1, &nb_elem);
	glBindBuffer(GL_ARRAY_BUFFER, nb_elem);
	glBufferData(GL_ARRAY_BUFFER, sizeof(NodeInfo)*node_num, (GLfloat *)node_info_array, GL_STATIC_DRAW);
	delete[] node_info_array;
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid *)0);
	glEnableVertexAttribArray(0);

	struct ElementInfo { GLuint n1, n2, n3; };
	ElementInfo *elem_info_array = new ElementInfo[elem_num];
	for (size_t i = 0; i < elem_num; ++i)
	{
		ElementInfo &elem_info = elem_info_array[i];
		elem_info.n1 = GLuint(tri_mesh.elems[i].n1);
		elem_info.n2 = GLuint(tri_mesh.elems[i].n2);
		elem_info.n3 = GLuint(tri_mesh.elems[i].n3);
	}
	glGenBuffers(1, &ib_elem);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ib_elem);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(ElementInfo)*elem_num, (GLuint *)elem_info_array, GL_STATIC_DRAW);
	delete[] elem_info_array;

	glBindVertexArray(0);
}

void DrawTriangleMesh::draw_tri_mesh(void)
{
	unicolor_shader.use();
	// position (test)
	//glm::mat4 rb_mv_mat = glm::rotate(glm::mat4(1.0f), glm::radians(10.0f), glm::vec3(0.0f, 0.0f, 1.0f))
	//					* glm::translate(glm::mat4(1.0f), glm::vec3(-rb_x_mc, -rb_y_mc, 0.0f));
	//unicolor_shader.set_uniform_matrix4f(uc_mv_mat_id, glm::value_ptr(rb_mv_mat));
	// color
	glm::vec4 light_yellow(1.0f, 0.95f, 0.5f, 1.0f);
	unicolor_shader.set_uniform_vec4f(uc_color_vec_id, glm::value_ptr(light_yellow));
	glBindVertexArray(va_elem);
	glLineWidth(3);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glDrawElements(GL_TRIANGLES, elem_num * 3, GL_UNSIGNED_INT, 0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

// need optimization for better memory usage
void DrawTriangleMesh::init_bl(TriangleMesh &tri_mesh)
{
	glGenVertexArrays(1, &va_bl);
	glBindVertexArray(va_bl);

	struct NodeInfo	{ GLfloat x, y, z; };
	NodeInfo *node_info_array = new NodeInfo[node_num];
	for (size_t i = 0; i < node_num; ++i)
	{
		NodeInfo &node_info = node_info_array[i];
		node_info.x = GLfloat(tri_mesh.nodes[i].x);
		node_info.y = GLfloat(tri_mesh.nodes[i].y);
		node_info.z = 0.0;
	}
	glGenBuffers(1, &nb_bl);
	glBindBuffer(GL_ARRAY_BUFFER, nb_bl);
	glBufferData(GL_ARRAY_BUFFER, sizeof(NodeInfo) * node_num, (GLfloat *)node_info_array, GL_STATIC_DRAW);
	delete[] node_info_array;
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 3, (GLvoid *)0);
	glEnableVertexAttribArray(0);

	glGenBuffers(1, &ib_bl);
	struct LineInfo { GLuint n1, n2; };
	LineInfo *line_info_array = new LineInfo[line_num];
	for (size_t i = 0; i < line_num; ++i)
	{
		LineInfo &line_info = line_info_array[i];
		line_info.n1 = GLsizei(tri_mesh.boundary_edges[i].n1);
		line_info.n2 = GLsizei(tri_mesh.boundary_edges[i].n2);
	}
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ib_bl);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(LineInfo)* line_num, (LineInfo *)line_info_array, GL_STATIC_DRAW);
	delete[] line_info_array;

	glBindVertexArray(0);
}

void DrawTriangleMesh::draw_bl(void)
{
	unicolor_shader.use();
	glm::vec4 red(1.0f, 0.0f, 0.0f, 1.0f);
	unicolor_shader.set_uniform_vec4f(uc_color_vec_id, glm::value_ptr(red));
	glLineWidth(4);
	glBindVertexArray(va_bl);
	glDrawElements(GL_LINES, line_num * 2, GL_UNSIGNED_INT, 0);
}

void DrawTriangleMesh::init_bg_grid(TriangleMesh &tri_mesh)
{
	TriangleMesh::BgGrid &bg_grid = tri_mesh.bg_grid;

	glGenVertexArrays(1, &va_bg_grid);
	glBindVertexArray(va_bg_grid);

	grid_num = GLsizei(bg_grid.num);
	GLsizei data_num = GLsizei(bg_grid.num * 36);
	GLfloat *elem_info = new GLfloat[data_num];
	GLfloat *cur_elem_info = elem_info;
	struct { GLfloat r, g, b; } grid_color;
	size_t grid_id = 0;
	for (long long j = 0; j < bg_grid.y_num; ++j)
		for (long long i = 0; i < bg_grid.x_num; ++i)
		{
			TriangleMesh::BgGrid::Grid &grid = bg_grid.grids[grid_id];
			switch (grid.pos_type)
			{
			case TriangleMesh::BgGrid::PosType::Outside:
				// Grey
				grid_color.r = 0.9f;
				grid_color.g = 0.9f;
				grid_color.b = 0.9f;
				break;
			case TriangleMesh::BgGrid::PosType::AtBoundary:
				// Dark Orange
				grid_color.r = 1.0f;
				grid_color.g = 0.54902f;
				grid_color.b = 0.0f;
				break;
			case TriangleMesh::BgGrid::PosType::Inside:
				// Cornflower Blue
				grid_color.r = 0.39216f;
				grid_color.g = 0.58431f;
				grid_color.b = 0.92941f;
				break;
			default:
				// Grey
				grid_color.r = 0.5f;
				grid_color.g = 0.5f;
				grid_color.b = 0.5f;
				break;
			}

			cur_elem_info[0] = GLfloat(i * bg_grid.h);
			cur_elem_info[1] = GLfloat(j * bg_grid.h);
			cur_elem_info[2] = 0.0;
			cur_elem_info[3] = grid_color.r;
			cur_elem_info[4] = grid_color.g;
			cur_elem_info[5] = grid_color.b;

			cur_elem_info[6] = GLfloat((i + 1) * bg_grid.h);
			cur_elem_info[7] = GLfloat(j    * bg_grid.h);
			cur_elem_info[8] = 0.0;
			cur_elem_info[9] = grid_color.r;
			cur_elem_info[10] = grid_color.g;
			cur_elem_info[11] = grid_color.b;

			cur_elem_info[12] = GLfloat((i + 1) * bg_grid.h);
			cur_elem_info[13] = GLfloat((j + 1) * bg_grid.h);
			cur_elem_info[14] = 0.0;
			cur_elem_info[15] = grid_color.r;
			cur_elem_info[16] = grid_color.g;
			cur_elem_info[17] = grid_color.b;

			cur_elem_info[18] = GLfloat(i * bg_grid.h);
			cur_elem_info[19] = GLfloat(j * bg_grid.h);
			cur_elem_info[20] = 0.0;
			cur_elem_info[21] = grid_color.r;
			cur_elem_info[22] = grid_color.g;
			cur_elem_info[23] = grid_color.b;

			cur_elem_info[24] = GLfloat((i + 1) * bg_grid.h);
			cur_elem_info[25] = GLfloat((j + 1) * bg_grid.h);
			cur_elem_info[26] = 0.0;
			cur_elem_info[27] = grid_color.r;
			cur_elem_info[28] = grid_color.g;
			cur_elem_info[29] = grid_color.b;

			cur_elem_info[30] = GLfloat(i    * bg_grid.h);
			cur_elem_info[31] = GLfloat((j + 1) * bg_grid.h);
			cur_elem_info[32] = 0.0;
			cur_elem_info[33] = grid_color.r;
			cur_elem_info[34] = grid_color.g;
			cur_elem_info[35] = grid_color.b;

			cur_elem_info += 36;

			++grid_id;
		}
	glGenBuffers(1, &nb_bg_grid);
	glBindBuffer(GL_ARRAY_BUFFER, nb_bg_grid);
	glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)* data_num, elem_info, GL_STATIC_DRAW);
	delete[] elem_info;
	// coordinates
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GL_FLOAT), (GLvoid *)0);
	glEnableVertexAttribArray(0);
	// color
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GL_FLOAT), (GLvoid *)(sizeof(GL_FLOAT) * 3));
	glEnableVertexAttribArray(1);
}

void DrawTriangleMesh::draw_bg_grid(void)
{
	multicolor_shader.use();
	glBindVertexArray(va_bg_grid);
	glDrawArrays(GL_TRIANGLES, 0, grid_num * 6);
}

void DrawTriangleMesh::init_grid_lines(TriangleMesh &tri_mesh)
{
	TriangleMesh::BgGrid &bg_grid = tri_mesh.bg_grid;

	glGenVertexArrays(1, &va_grid_lines);
	glBindVertexArray(va_grid_lines);

	struct LineInfo	{ GLfloat x1, y1, z1, x2, y2, z2; };
	grid_line_num = GLuint(bg_grid.x_num + bg_grid.y_num + 2);
	LineInfo *line_infos = new LineInfo[grid_line_num];
	size_t min_line_num, max_line_num;
	// vertical lines
	min_line_num = 0;
	max_line_num = bg_grid.x_num + 1;
	for (size_t i = min_line_num; i < max_line_num; ++i)
	{
		LineInfo &line_info = line_infos[i];
		// node 1
		line_info.x1 = GLfloat(i * bg_grid.h);
		line_info.y1 = 0.0;
		line_info.z1 = 0.0;
		// node 2
		line_info.x2 = GLfloat(i * bg_grid.h);
		line_info.y2 = GLfloat(bg_grid.yn);
		line_info.z2 = 0.0;
	}
	// horizontal lines
	min_line_num = bg_grid.x_num + 1;
	max_line_num = grid_line_num;
	for (size_t i = min_line_num; i < max_line_num; ++i)
	{
		LineInfo &line_info = line_infos[i];
		// node 1
		line_info.x1 = 0.0;
		line_info.y1 = GLfloat((i - min_line_num) * bg_grid.h);
		line_info.z1 = 0.0;
		// node 2
		line_info.x2 = GLfloat(bg_grid.xn);
		line_info.y2 = GLfloat((i - min_line_num) * bg_grid.h);
		line_info.z2 = 0.0;
	}
	glGenBuffers(1, &nb_grid_line_coords);
	glBindBuffer(GL_ARRAY_BUFFER, nb_grid_line_coords);
	glBufferData(GL_ARRAY_BUFFER, sizeof(LineInfo) * grid_line_num, line_infos, GL_STATIC_DRAW);
	delete[] line_infos;

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (void *)0);
	glEnableVertexAttribArray(0);
}

void DrawTriangleMesh::draw_grid_lines(void)
{
	unicolor_shader.use();
	glm::vec4 grey(0.7f, 0.7f, 0.7f, 1.0f);
	unicolor_shader.set_uniform_vec4f(uc_color_vec_id, glm::value_ptr(grey));
	glBindVertexArray(va_grid_lines);
	glLineWidth(2);
	glDrawArrays(GL_LINES, 0, grid_line_num * 2);
}

void DrawTriangleMesh::init_line_and_point(TriangleMesh &tri_mesh, 
										   TriangleMesh::Edge &edge, Point &p)
{
	glGenVertexArrays(1, &va_line_and_point);
	glBindVertexArray(va_line_and_point);
	
	TriangleMesh::Node &n1 = tri_mesh.nodes[edge.n1];
	TriangleMesh::Node &n2 = tri_mesh.nodes[edge.n2];
	GLfloat data[] = {
		GLfloat(n1.x), GLfloat(n1.y), 0.0f,
		GLfloat(n2.x), GLfloat(n2.y), 0.0f,
		GLfloat( p.x), GLfloat( p.y), 0.0f,
	};
	GLuint index[] = { 0, 1 };
	glGenBuffers(1, &nb_line_and_point);
	glBindBuffer(GL_ARRAY_BUFFER, nb_line_and_point);
	glBufferData(GL_ARRAY_BUFFER, sizeof(data), data, GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 3, (GLvoid *)0);
	glEnableVertexAttribArray(0);

	glGenBuffers(1, &ib_line);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ib_line);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(index), index, GL_STATIC_DRAW);

	// adjust display range
	Rect disp_range = tri_mesh.get_display_range();
	if (disp_range.xl > p.x)
		disp_range.xl = p.x;
	if (disp_range.xu < p.x)
		disp_range.xu = p.x;
	if (disp_range.yl > p.y)
		disp_range.yl = p.y;
	if (disp_range.yu < p.y)
		disp_range.yu = p.y;
	double x_range = disp_range.xu - disp_range.xl;
	double y_range = disp_range.yu - disp_range.yl;
	GLfloat range = (GLfloat)(x_range > y_range ? x_range : y_range);
	GLfloat padding = range * 0.05f;
	GLfloat view_xl = GLfloat(disp_range.xl) - padding;
	GLfloat view_xu = GLfloat(disp_range.xl) + range + padding;
	GLfloat view_yl = GLfloat(disp_range.yl) - padding;
	GLfloat view_yu = GLfloat(disp_range.yl) + range + padding;
	glm::mat4 proj_mat = glm::ortho(view_xl, view_xu, view_yl, view_yu);
	unicolor_shader.use();
	unicolor_shader.set_uniform_matrix4f(uc_proj_mat_id, glm::value_ptr(proj_mat));
	multicolor_shader.use();
	multicolor_shader.set_uniform_matrix4f(mc_proj_mat_id, glm::value_ptr(proj_mat));
}

void DrawTriangleMesh::draw_line_and_point(void)
{
	unicolor_shader.use();
	//glm::vec4 purple(0.2941f, 0.0f, 0.5098f, 1.0f);
	glm::vec4 pink(1.0f, 0.752941f, 0.796078, 1.0f);
	unicolor_shader.set_uniform_vec4f(uc_color_vec_id, glm::value_ptr(pink));

	glBindVertexArray(va_line_and_point);

	glPointSize(10);
	glDrawArrays(GL_POINTS, 2, 1);
	glPointSize(1); // initial value

	glLineWidth(4);
	glDrawElements(GL_LINES, 2, GL_UNSIGNED_INT, 0);
}