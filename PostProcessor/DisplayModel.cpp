#include "PostProcessor_pcp.h"

#include "DisplayModel.h"

DisplayModel::DisplayModel() :
	rb_vao(0), rb_vbo(0), rb_ibo(0), rb_id_num(0),
	bgm_vao(0), bgm_vbo(0), grid_line_points_num(0),
	pcls_vao(0), pcls_vbo(0), pcls_id_num(0) {}
DisplayModel::~DisplayModel()
{
	if (rb_vao)
	{
		glDeleteVertexArrays(1, &rb_vao);
		glDeleteBuffers(1, &rb_vbo);
		glDeleteBuffers(1, &rb_ibo);
	}
	if (bgm_vao)
	{
		glDeleteVertexArrays(1, &bgm_vao);
		glDeleteBuffers(1, &bgm_vbo);
	}
	if (pcls_vao)
	{
		glDeleteVertexArrays(1, &pcls_vao);
		glDeleteBuffers(1, &pcls_vbo);
	}
}

static void processInput(GLFWwindow *window)
{
	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
		glfwSetWindowShouldClose(window, true);
}

static double h_div_w, w_div_h;
static void resize_window(GLFWwindow *window, GLsizei width, GLsizei height)
{
	GLsizei vp_width, vp_height, padding;
	vp_height = width * h_div_w;
	if (vp_height < height)
	{
		vp_width = width;
		padding = (height - vp_height) / 2;
		glViewport(0, padding, vp_width, vp_height);
	}
	else
	{
		vp_height = height;
		vp_width = height * w_div_h;
		padding = (width - vp_width) / 2;
		glViewport(padding, 0, vp_width, vp_height);
	}
}

static const size_t SCREEN_WIDTH = 600;
static const size_t SCREEN_HEIGHT = 600;

int DisplayModel::init_win(void)
{
	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	window = glfwCreateWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "model display", nullptr, nullptr);
	if (!window)
	{
		std::cout << "Failed to create GLFW window.\n";
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);
	glfwSetFramebufferSizeCallback(window, resize_window);

	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
	{
		std::cout << "Failed to initialize GLAD.\n";
		glfwTerminate();
		return -2;
	}

	shader.init_from_file("..\\..\\Asset\\unicolor_vshader.txt",
						  "..\\..\\Asset\\unicolor_fshader.txt");
	shader.use();
	color_id = shader.init_uniform("color_vec");
	if (color_id < 0) return -3;
	mv_mat_id = shader.init_uniform("mv_mat");
	if (mv_mat_id < 0) return -3;
	shader.set_uniform_matrix4f(mv_mat_id, glm::value_ptr(glm::mat4(1.0f)));
	proj_mat_id = shader.init_uniform("proj_mat");
	if (proj_mat_id < 0) return -3;

	return 0;
}

int DisplayModel::init_rigid_body(RigidBody &rb)
{
	TriangleMesh &mesh = rb.mesh;

	glm::mat4 identity_mat(1.0f);
	rb_mv_mat = glm::translate(identity_mat, glm::vec3(rb.x, rb.y, 0.0f))
		* glm::rotate(identity_mat, GLfloat(rb.theta), glm::vec3(0.0f, 0.0f, 1.0f))
		* glm::translate(identity_mat, glm::vec3(-mesh.get_x_mc(), -mesh.get_y_mc(), 0.0f));

	glGenVertexArrays(1, &rb_vao);
	glBindVertexArray(rb_vao);
	glGenBuffers(1, &rb_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, rb_vbo);
	size_t node_num = mesh.get_node_num();
	TriangleMesh::Node *nodes = mesh.get_nodes();
	GLfloat *mesh_node_coords = new GLfloat[node_num * 3];
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		mesh_node_coords[3 * n_id] = GLfloat(nodes[n_id].x);
		mesh_node_coords[3 * n_id + 1] = GLfloat(nodes[n_id].y);
		mesh_node_coords[3 * n_id + 2] = 0.0f;
	}
	glBufferData(GL_ARRAY_BUFFER, node_num * 3 * sizeof(GLfloat), mesh_node_coords, GL_STATIC_DRAW);
	delete[] mesh_node_coords;
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 3, (GLvoid *)0);
	glEnableVertexAttribArray(0);
	glGenBuffers(1, &rb_ibo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, rb_ibo);
	size_t elem_num = mesh.get_elem_num();
	rb_id_num = elem_num * 3;
	TriangleMesh::Element *elems = mesh.get_elems();
	GLuint *mesh_elem_indices = new GLuint[rb_id_num];
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		mesh_elem_indices[3 * e_id] = GLuint(elems[e_id].n1);
		mesh_elem_indices[3 * e_id + 1] = GLuint(elems[e_id].n2);
		mesh_elem_indices[3 * e_id + 2] = GLuint(elems[e_id].n3);
	}
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, rb_id_num * sizeof(GLuint), mesh_elem_indices, GL_STATIC_DRAW);
	delete[] mesh_elem_indices;
	return 0;
}

int DisplayModel::init_bg_mesh(double grid_size,
	size_t _elem_x_num, size_t _elem_y_num, double x_start, double y_start)
{
	bgm_mv_mat = glm::mat4(1.0f);

	size_t node_x_num = _elem_x_num + 1;
	size_t node_y_num = _elem_y_num + 1;
	grid_line_points_num = GLsizei((node_y_num + node_x_num) * 2);
	GLfloat *mesh_line_coords = new GLfloat[grid_line_points_num * 3];
	for (size_t line_id = 0; line_id < node_y_num; ++line_id)
	{
		mesh_line_coords[6 * line_id] = GLfloat(x_start);
		mesh_line_coords[6 * line_id + 1] = GLfloat(y_start + grid_size * line_id);
		mesh_line_coords[6 * line_id + 2] = 0.0f;
		mesh_line_coords[6 * line_id + 3] = GLfloat(x_start + grid_size * _elem_x_num);
		mesh_line_coords[6 * line_id + 4] = GLfloat(y_start + grid_size * line_id);
		mesh_line_coords[6 * line_id + 5] = 0.0f;
	}
	size_t line_coord_num = node_y_num * 6;
	for (size_t col_id = 0; col_id < node_x_num; ++col_id)
	{
		mesh_line_coords[line_coord_num + 6 * col_id] = GLfloat(x_start + grid_size * col_id);
		mesh_line_coords[line_coord_num + 6 * col_id + 1] = GLfloat(y_start);
		mesh_line_coords[line_coord_num + 6 * col_id + 2] = 0.0f;
		mesh_line_coords[line_coord_num + 6 * col_id + 3] = GLfloat(x_start + grid_size * col_id);
		mesh_line_coords[line_coord_num + 6 * col_id + 4] = GLfloat(y_start + grid_size * _elem_y_num);
		mesh_line_coords[line_coord_num + 6 * col_id + 5] = 0.0f;
	}
	glGenVertexArrays(1, &bgm_vao);
	glBindVertexArray(bgm_vao);
	glGenBuffers(1, &bgm_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, bgm_vbo);
	glBufferData(GL_ARRAY_BUFFER, grid_line_points_num * 3 * sizeof(GLfloat), mesh_line_coords, GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 3, (GLvoid *)0);
	glEnableVertexAttribArray(0);
	delete[] mesh_line_coords;
	return 0;
}


int DisplayModel::display(double left, double right, double bottom, double top)
{
	// adjust viewport
	double cam_width, cam_height;
	cam_width = right - left;
	cam_height = top - bottom;
	h_div_w = cam_height / cam_width;
	w_div_h = cam_width / cam_height;
	resize_window(window, SCREEN_WIDTH, SCREEN_HEIGHT);

	shader.use();
	glm::mat4 proj_mat = glm::ortho(left, right, bottom, top);
	shader.set_uniform_matrix4f(proj_mat_id, glm::value_ptr(proj_mat));
	
	while (!glfwWindowShouldClose(window))
	{
		processInput(window);
		render();
		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	glfwTerminate();
	return 0;
}


void DisplayModel::render(void)
{
	// background
	glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	shader.use();
	// background mesh
	if (bgm_vao)
	{
		shader.set_uniform_matrix4f(mv_mat_id, glm::value_ptr(bgm_mv_mat));
		glm::vec4 grey(0.6f, 0.6f, 0.6f, 1.0f);
		shader.set_uniform_vec4f(color_id, glm::value_ptr(grey));
		glBindVertexArray(bgm_vao);
		glDrawArrays(GL_LINES, 0, grid_line_points_num);
	}
	// particles
	if (pcls_vao)
	{
		shader.set_uniform_matrix4f(mv_mat_id, glm::value_ptr(pcls_mv_mat));
		glm::vec4 yellow(1.0f, 0.8f, 0.0f, 1.0f);
		shader.set_uniform_vec4f(color_id, glm::value_ptr(yellow));
		glBindVertexArray(pcls_vao);
		glPointSize(10.0f);
		glDrawArrays(GL_POINTS, 0, pcls_id_num);
	}
	// rigid body
	if (rb_vao)
	{
		shader.set_uniform_matrix4f(mv_mat_id, glm::value_ptr(rb_mv_mat));
		glm::vec4 orange(1.0f, 0.5f, 0.0f, 1.0f);
		shader.set_uniform_vec4f(color_id, glm::value_ptr(orange));
		glBindVertexArray(rb_vao);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glDrawElements(GL_TRIANGLES, rb_id_num, GL_UNSIGNED_INT, 0);
	}
}
