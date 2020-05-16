#include "PostProcessor_pcp.h"

#include "DisplayModel_S2D.h"

DisplayModel_S2D::DisplayModel_S2D() :
	bgm_vao(0), bgm_vbo(0), bgm_elem_n_id_num(0),
	pcls_vao(0), pcls_vbo(0), pcls_id_num(0),
	pt_vao(0), pt_vbo(0), pt_num(0) {}

DisplayModel_S2D::~DisplayModel_S2D()
{
	if (bgm_vao)
	{
		glDeleteBuffers(1, &bgm_vbo);
		glDeleteBuffers(1, &bgm_ibo);
		glDeleteVertexArrays(1, &bgm_vao);
	}
	bgm_vbo = 0;
	bgm_ibo = 0;
	bgm_vao = 0;
	if (pcls_vao)
	{
		glDeleteBuffers(1, &pcls_vbo);
		glDeleteBuffers(1, &pcls_veo);
		glDeleteVertexArrays(1, &pcls_vao);
	}
	pcls_vbo = 0;
	pcls_veo = 0;
	pcls_vao = 0;
	if (pt_vao)
	{
		glDeleteBuffers(1, &pt_vbo);
		glDeleteVertexArrays(1, &pt_vao);
	}
	pt_vbo = 0;
	pt_vao = 0;
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

int DisplayModel_S2D::init_win(void)
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

int DisplayModel_S2D::display(double left, double right, double bottom, double top)
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

int DisplayModel_S2D::init_model(Model_S2D_ME_s &md)
{
	int res;
	res = init_bg_mesh(md);
	if (res) return res;
	res = init_pcls(md);
	if (res) return res;

	return 0;
}

int DisplayModel_S2D::init_points(GLfloat *n_coords, size_t num)
{
	pt_mv_mat = glm::mat4(1.0f);

	pt_num = num;
	glGenVertexArrays(1, &pt_vao);
	glBindVertexArray(pt_vao);
	glGenBuffers(1, &pt_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, pt_vbo);
	glBufferData(GL_ARRAY_BUFFER, num * 3 * sizeof(GLfloat), n_coords, GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 3, (GLvoid *)0);
	glEnableVertexAttribArray(0);

	return 0;
}

int DisplayModel_S2D::init_bg_mesh(Model_S2D_ME_s &md)
{
	bgm_mv_mat = glm::mat4(1.0f);

	glGenVertexArrays(1, &bgm_vao);
	glBindVertexArray(bgm_vao);
	
	size_t node_x_num = md.get_elem_x_num() + 1;
	size_t node_y_num = md.get_elem_y_num() + 1;
	bgm_elem_n_id_num = (node_x_num + node_y_num) * 2;
	Model_S2D_ME_s::Node *nodes = md.get_nodes();
	GLfloat *node_coords = new GLfloat[bgm_elem_n_id_num * 3];
	GLuint *elem_indices = new GLuint[bgm_elem_n_id_num];
	GLfloat x0 = md.get_x0();
	GLfloat xn = md.get_xn();
	GLfloat y0 = md.get_y0();
	GLfloat yn = md.get_yn();
	GLfloat hx = md.get_hx();
	GLfloat hy = md.get_hy();
	for (size_t n_id = 0; n_id < node_x_num; ++n_id)
	{
		node_coords[6*n_id]   = x0 + GLfloat(n_id) * hx;
		node_coords[6*n_id+1] = y0;
		node_coords[6*n_id+2] = 0.0f;
		node_coords[6*n_id+3] = x0 + GLfloat(n_id) * hx;
		node_coords[6*n_id+4] = yn;
		node_coords[6*n_id+5] = 0.0f;
		elem_indices[2*n_id] = 2 * n_id;
		elem_indices[2*n_id+1] = 2 * n_id + 1;
	}
	std::cout << "\n\n";
	size_t node_num_tmp = node_x_num + node_y_num;
	for (size_t n_id = node_x_num; n_id < node_num_tmp; ++n_id)
	{
		node_coords[6*n_id]   = x0;
		node_coords[6*n_id+1] = y0 + GLfloat(n_id - node_x_num) * hy;
		node_coords[6*n_id+2] = 0.0f;
		node_coords[6*n_id+3] = xn;
		node_coords[6*n_id+4] = y0 + GLfloat(n_id - node_x_num) * hy;
		node_coords[6*n_id+5] = 0.0f;
		elem_indices[2*n_id] = 2 * n_id;
		elem_indices[2*n_id+1] = 2 * n_id + 1;
	}

	glGenBuffers(1, &bgm_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, bgm_vbo);
	glBufferData(GL_ARRAY_BUFFER, bgm_elem_n_id_num * 3 * sizeof(GLfloat), node_coords, GL_STATIC_DRAW);
	delete[] node_coords;
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 3, (GLvoid *)0);
	glEnableVertexAttribArray(0);

	glGenBuffers(1, &bgm_ibo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bgm_ibo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, bgm_elem_n_id_num * sizeof(GLuint), elem_indices, GL_STATIC_DRAW);
	delete[] elem_indices;

	return 0;
}

int DisplayModel_S2D::init_pcls(Model_S2D_ME_s &md)
{
	pcls_mv_mat = glm::mat4(1.0f);

	glGenVertexArrays(1, &pcls_vao);
	glBindVertexArray(pcls_vao);
	pcls_id_num = md.get_pcl_num();
	pcls_mem.reset();
	Model_S2D_ME_s::Particle *pcls = md.get_pcls();
	for (size_t pcl_id = 0; pcl_id < pcls_id_num; ++pcl_id)
	{
		Model_S2D_ME_s::Particle &pcl = pcls[pcl_id];
		pcls_mem.add_pcl(pcl.x, pcl.y, pcl.m / pcl.density * 0.25);
	}
	glGenBuffers(1, &pcls_vbo);
	glBindBuffer(GL_ARRAY_BUFFER, pcls_vbo);
	glBufferData(GL_ARRAY_BUFFER, pcls_mem.get_point_num() * sizeof(GLfloat), pcls_mem.get_pcls(), GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 3, (GLvoid *)0);
	glEnableVertexAttribArray(0);
	glGenBuffers(1, &pcls_veo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, pcls_veo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, pcls_mem.get_index_num() * sizeof(GLuint), pcls_mem.get_indices(), GL_STATIC_DRAW);

	return 0;
}

void DisplayModel_S2D::render(void)
{
	// background
	//glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	shader.use();

	// background mesh
	if (bgm_vao)
	{
		shader.set_uniform_matrix4f(mv_mat_id, glm::value_ptr(bgm_mv_mat));
		glm::vec4 grey(0.6f, 0.6f, 0.6f, 1.0f);
		shader.set_uniform_vec4f(color_id, glm::value_ptr(grey));
		glLineWidth(1.5f);
		glBindVertexArray(bgm_vao);
		glDrawElements(GL_LINES, bgm_elem_n_id_num, GL_UNSIGNED_INT, 0);
	}

	// particles
	if (pcls_vao)
	{
		shader.set_uniform_matrix4f(mv_mat_id, glm::value_ptr(pcls_mv_mat));
		//glm::vec4 yellow(1.0f, 0.8f, 0.0f, 1.0f);
		glm::vec4 yellow(0.28235f, 0.23921f, 0.54510f, 1.0f);
		shader.set_uniform_vec4f(color_id, glm::value_ptr(yellow));
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glBindVertexArray(pcls_vao);
		glDrawElements(GL_TRIANGLES, pcls_mem.get_index_num(), GL_UNSIGNED_INT, (GLvoid *)0);
	}

	// pt
	if (pt_vao)
	{
		shader.set_uniform_matrix4f(mv_mat_id, glm::value_ptr(pt_mv_mat));
		glm::vec4 red(1.0f, 0.0f, 0.0f, 1.0f);
		shader.set_uniform_vec4f(color_id, glm::value_ptr(red));
		glBindVertexArray(pt_vao);
		glPointSize(10.0f);
		glDrawArrays(GL_POINTS, 0, pt_num);
	}
}
