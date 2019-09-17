#include "PostProcessor_pcp.h"

#include "gif.h"

#include "GenerateAnimation.h"

GenerateAnimation::GenerateAnimation() :
	time_rcds(nullptr), mp_x_data(nullptr) {}

GenerateAnimation::~GenerateAnimation()
{
	if (time_rcds)
	{
		delete[] time_rcds;
		time_rcds = nullptr;
	}
	if (mp_x_data)
	{
		delete[] mp_x_data;
		mp_x_data = nullptr;
	}
}

static void processInput(GLFWwindow *window)
{
	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
		glfwSetWindowShouldClose(window, true);
}

static void framebuffer_size_callback(GLFWwindow *window, int width, int height)
{
	GLsizei win_size;
	GLsizei win_padding;
	if (width < height)
	{
		win_size = width;
		win_padding = (height - width) / 2;
		glViewport(0, win_padding, win_size, win_size);
	}
	else
	{
		win_size = height;
		win_padding = (width - height) / 2;
		glViewport(win_padding, 0, win_size, win_size);
	}
}

static const size_t SCREEN_WIDTH = 600;
static const size_t SCREEN_HEIGHT = 600;

void reorder_buffer(unsigned char *RGBA_data, int width, int height)
{
	// RGBA data are 4 bytes long
	long *data = reinterpret_cast<long *>(RGBA_data);
	long *line1 = data;
	long *line2 = data + (height - 1) * width;
	long data_tmp;
	while (line1 < line2)
	{
		for (size_t i = 0; i < width; i++)
		{
			data_tmp = line1[i];
			line1[i] = line2[i];
			line2[i] = data_tmp;
		}
		line1 += width;
		line2 -= width;
	}
}

int GenerateAnimation::generate(double ani_time, double xl, double xu, double yl, double yu,
								const char *res_file_name, const char *gif_name)
{
	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	GLFWwindow *window = glfwCreateWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "model display", nullptr, nullptr);
	if (window == nullptr)
	{
		std::cout << "Failed to create GLFW window.\n";
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
	{
		std::cout << "Failed to initialize GLAD.\n";
		glfwTerminate();
		return -2;
	}

	// init animation
	init(res_file_name);

	// animation_time
	real_time = time_rcds[time_rcd_num-1].total_time - time_rcds[0].total_time;
	animation_time = ani_time;
	ani_real_ratio = animation_time / real_time;
	min_delay_real = 0.01 / ani_real_ratio;

	// gif output
	GifWriter gif_file;
	unsigned char *pixels_data;
	if (gif_name)
	{
		GifBegin(&gif_file, gif_name, SCREEN_WIDTH, SCREEN_HEIGHT, 1);
		pixels_data = new unsigned char[SCREEN_WIDTH * SCREEN_HEIGHT * 4];
	}

	bool render_new_frame = true; // control whether it is time to swap and render new frame
	bool reach_last_frame = false; // whether run out of frame to render
	size_t draw_frame_id = 0;
	double prev_time;
	while (!glfwWindowShouldClose(window))
	{
		processInput(window);

		if (render_new_frame)
		{
			prev_time = glfwGetTime();
			glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
			glClear(GL_COLOR_BUFFER_BIT);
			render_frame(xl, xu, yl, yu);
			render_new_frame = false;
			// output to gif file
			if (!reach_last_frame)
			{
				std::cout << "frame " << draw_frame_id++
					<< " real time " << time_rcds[cur_time_rcd_id].total_time
					<< " ani time " << time_rcds[cur_time_rcd_id].total_substep_num << "\n";
				reach_last_frame = !find_next_frame();
				if (gif_name)
				{
					// read pixel from back buffer
					glReadPixels(0, 0, SCREEN_WIDTH, SCREEN_HEIGHT, GL_RGBA, GL_UNSIGNED_BYTE, pixels_data);
					reorder_buffer(pixels_data, SCREEN_WIDTH, SCREEN_HEIGHT);
					GifWriteFrame(&gif_file, pixels_data, SCREEN_WIDTH, SCREEN_HEIGHT, (delay_ani_100th ? delay_ani_100th : 1));
				}
			}
		}

		if (glfwGetTime() - prev_time >= delay_ani)
		{
			glfwSwapBuffers(window);
			render_new_frame = true;
		}

		glfwPollEvents();
	}

	if (gif_name)
	{
		GifEnd(&gif_file);
		delete[] pixels_data;
		pixels_data = nullptr;
	}

	glfwTerminate();
	return 0;
}


int GenerateAnimation::init(const char *res_file_name)
{
	res_file.open(res_file_name, std::ios_base::binary | std::ios_base::in);
	if (!res_file.is_open())
		return -1;

	// init shader
	shader.init_from_file("..\\..\\Asset\\unicolor_vshader.txt",
						  "..\\..\\Asset\\unicolor_fshader.txt");
	shader.use();
	color_id = shader.init_uniform("color_vec");
	if (color_id < 0)
		return -2;
	mv_mat_id = shader.init_uniform("mv_mat");
	if (mv_mat_id < 0)
		return -2;
	proj_mat_id = shader.init_uniform("proj_mat");
	if (proj_mat_id < 0)
		return -2;

	// mesh
	res_file.read(reinterpret_cast<char *>(&mh), sizeof(mh));
	// init bg grid data
	size_t node_x_num = mh.elem_x_num + 1;
	size_t node_y_num = mh.elem_y_num + 1;
	GLfloat *mesh_line_coords = new GLfloat[(node_y_num + node_x_num) * 6];
	grid_line_points_num = GLsizei((node_y_num + node_x_num) * 2);
	GLuint *mesh_line_indices = new GLuint[grid_line_points_num];
	for (size_t line_id = 0; line_id < node_y_num; ++line_id)
	{
		mesh_line_coords[6 * line_id] = GLfloat(mh.x0);
		mesh_line_coords[6 * line_id + 1] = GLfloat(mh.y0 + mh.h * line_id);
		mesh_line_coords[6 * line_id + 2] = 0.0f;
		mesh_line_coords[6 * line_id + 3] = GLfloat(mh.x0 + mh.h * mh.elem_x_num);
		mesh_line_coords[6 * line_id + 4] = GLfloat(mh.y0 + mh.h * line_id);
		mesh_line_coords[6 * line_id + 5] = 0.0f;
		mesh_line_indices[2 * line_id] = GLuint(2 * line_id);
		mesh_line_indices[2 * line_id + 1] = GLuint(2 * line_id + 1);
	}
	size_t line_coord_num = node_y_num * 6;
	size_t line_indices_num = node_y_num * 2;
	for (size_t col_id = 0; col_id < node_x_num; ++col_id)
	{
		mesh_line_coords[line_coord_num + 6 * col_id] = GLfloat(mh.x0 + mh.h * col_id);
		mesh_line_coords[line_coord_num + 6 * col_id + 1] = GLfloat(mh.y0);
		mesh_line_coords[line_coord_num + 6 * col_id + 2] = 0.0f;
		mesh_line_coords[line_coord_num + 6 * col_id + 3] = GLfloat(mh.x0 + mh.h * col_id);
		mesh_line_coords[line_coord_num + 6 * col_id + 4] = GLfloat(mh.y0 + mh.h * mh.elem_y_num);
		mesh_line_coords[line_coord_num + 6 * col_id + 5] = 0.0f;
		mesh_line_indices[line_indices_num + 2 * col_id] = GLuint(line_indices_num + 2 * col_id);
		mesh_line_indices[line_indices_num + 2 * col_id + 1] = GLuint(line_indices_num + 2 * col_id + 1);
	}
	bg_grid_data.init_array_buffer(mesh_line_coords, (node_y_num + node_x_num) * 6);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 3, (GLvoid *)0);
	glEnableVertexAttribArray(0);
	bg_grid_data.init_elem_array_buffer(mesh_line_indices, (node_y_num + node_x_num) * 2);
	//for (size_t i = 0; i < (node_y_num + node_x_num); i++)
	//{
	//	std::cout << "(" << mesh_line_coords[6 * i] << ", " << mesh_line_coords[6 * i + 1] << ", " << mesh_line_coords[6 * i + 2] << ") - ("
	//		<< mesh_line_coords[6 * i + 3] << ", " << mesh_line_coords[6 * i + 4] << ", " << mesh_line_coords[6 * i + 5] << ")  "
	//		<< mesh_line_indices[2 * i] << " - " << mesh_line_indices[2 * i + 1] << "\n";
	//}
	delete[] mesh_line_coords;
	delete[] mesh_line_indices;

	// rigid body
	res_file.read(reinterpret_cast<char *>(&rbh), sizeof(rbh));
	// node coordinates
	double *rb_node_coords = new double[rbh.node_num * 2];
	res_file.read(reinterpret_cast<char *>(rb_node_coords), rbh.node_num * 2 * sizeof(double));
	// element topology
	unsigned long long *rb_elem_indices = new unsigned long long[rbh.elem_num * 3];
	res_file.read(reinterpret_cast<char *>(rb_elem_indices), rbh.elem_num * 3 * sizeof(unsigned long long));
	rb_elem_point_num = GLsizei(rbh.elem_num * 3);
	// init rigid body data
	GLfloat *rb_coords = new GLfloat[rbh.node_num * 3];
	GLuint *rb_indices = new GLuint[rbh.elem_num * 3];
	for (size_t n_id = 0; n_id < rbh.node_num; ++n_id)
	{
		rb_coords[n_id * 3] = GLfloat(rb_node_coords[n_id * 2]);
		rb_coords[n_id * 3 + 1] = GLfloat(rb_node_coords[n_id * 2 + 1]);
		rb_coords[n_id * 3 + 2] = 0.0f;
		//std::cout << n_id << ", " << rb_coords[n_id * 3] << ", " << rb_coords[n_id * 3 + 1]
		//	<< ", " << rb_coords[n_id * 3 + 2] << "\n";
	}
	for (size_t e_id = 0; e_id < rbh.elem_num; ++e_id)
	{
		rb_indices[e_id * 3] = GLuint(rb_elem_indices[e_id * 3]);
		rb_indices[e_id * 3 + 1] = GLuint(rb_elem_indices[e_id * 3 + 1]);
		rb_indices[e_id * 3 + 2] = GLuint(rb_elem_indices[e_id * 3 + 2]);
		//std::cout << e_id << ", " << rb_indices[e_id * 3] << ", " << rb_indices[e_id * 3 + 1]
		//		<< ", " << rb_indices[e_id * 3 + 2] << "\n";
	}
	rigid_body_data.init_array_buffer(rb_coords, rbh.node_num * 3);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 3, (GLvoid *)0);
	glEnableVertexAttribArray(0);
	rigid_body_data.init_elem_array_buffer(rb_indices, rbh.elem_num * 3);
	delete[] rb_coords;
	delete[] rb_indices;
	delete[] rb_node_coords;
	delete[] rb_elem_indices;

	// material point object
	res_file.read(reinterpret_cast<char *>(&mph), sizeof(mph));
	if (!mp_x_data)
	{
		delete[] mp_x_data;
		mp_x_data = nullptr;
	}
	if (mph.pcl_num)
	{
		mp_x_data = new double[mph.pcl_num * 3];
		mp_y_data = mp_x_data + mph.pcl_num;
		mp_vol_data = mp_y_data + mph.pcl_num;
	}

	// get total number of time record
	size_t file_len;
	first_time_rcd_file_pos = res_file.tellg();
	res_file.seekg(0, SEEK_END);
	file_len = res_file.tellg();
	time_rcd_len = sizeof(TimeHistoryHeader) + sizeof(RigidBodyMotionHeader) + sizeof(MPObjectHeader) + mph.pcl_num * sizeof(double) * 3;
	time_rcd_num = (file_len - first_time_rcd_file_pos) / time_rcd_len;
	// init time record
	if (time_rcds)
	{
		delete[] time_rcds;
		time_rcds = nullptr;
	}
	if (time_rcd_num)
	{
		size_t cur_file_pos = first_time_rcd_file_pos;
		time_rcds = new TimeHistoryHeader[time_rcd_num];
		for (size_t tr_id = 0; tr_id < time_rcd_num; ++tr_id)
		{
			res_file.seekg(cur_file_pos, SEEK_SET);
			res_file.read(reinterpret_cast<char *>(time_rcds + tr_id), sizeof(TimeHistoryHeader));
			cur_file_pos += time_rcd_len;
		}
	}
	
	cur_time_rcd_id = 0;

	return 0;
}

int GenerateAnimation::render_frame(double xl, double xu, double yl, double yu)
{
	shader.use();
	// set display scope
	glm::mat4 proj_mat = glm::ortho(xl, xu, yl, yu); // to be improved
	shader.set_uniform_matrix4f(proj_mat_id, glm::value_ptr(proj_mat));

	// draw bg grid
	// model/view matrix
	glm::mat4 identity_mat = glm::mat4(1.0f);
	shader.set_uniform_matrix4f(mv_mat_id, glm::value_ptr(identity_mat));
	// color
	glm::vec4 white(1.0f, 1.0f, 1.0f, 1.0f);
	shader.set_uniform_vec4f(color_id, glm::value_ptr(white));
	bg_grid_data.use();
	// draw
	//glLineWidth(1);
	glDrawElements(GL_LINES, grid_line_points_num, GL_UNSIGNED_INT, (GLvoid *)0);
	
	// draw material points
	// model/view matrix
	shader.set_uniform_matrix4f(mv_mat_id, glm::value_ptr(identity_mat));
	// color
	glm::vec4 moccasin(1.0f, 0.894118f, 0.709804f, 1.0f); // yellow
	shader.set_uniform_vec4f(color_id, glm::value_ptr(moccasin));
	// init object buffer
	RigidBodyMotionHeader rbmh;
	MPObjectHeader mph;
	res_file.seekg(first_time_rcd_file_pos + cur_time_rcd_id * time_rcd_len + sizeof(TimeHistoryHeader), SEEK_SET);
	res_file.read(reinterpret_cast<char *>(&rbmh), sizeof(rbmh));
	res_file.read(reinterpret_cast<char *>(&mph), sizeof(mph));
	res_file.read(reinterpret_cast<char *>(mp_x_data), sizeof(double) * mph.pcl_num * 3);
	pcls_mem.reset();
	for (size_t pcl_id = 0; pcl_id < mph.pcl_num; ++pcl_id)
		pcls_mem.add_pcl(mp_x_data[pcl_id], mp_y_data[pcl_id], mp_vol_data[pcl_id]);
	mp_data.clear();
	mp_data.init_array_buffer(pcls_mem.get_pcls(), pcls_mem.get_point_num());
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 3, (GLvoid *)0);
	glEnableVertexAttribArray(0);
	mp_data.init_elem_array_buffer(pcls_mem.get_indices(), pcls_mem.get_index_num());
	// draw particles
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glDrawElements(GL_TRIANGLES, pcls_mem.get_index_num(), GL_UNSIGNED_INT, (GLvoid *)0);
	// draw particle boundary
	// color
	glm::vec4 red(1.0f, 0.0f, 0.0f, 1.0f);
	shader.set_uniform_vec4f(color_id, glm::value_ptr(red));
	GLuint pb_ibo = 0;
	mp_data.init_add_elem_array_buffer(pb_ibo, pcls_mem.get_bl_indices(), pcls_mem.get_bl_index_num());
	glLineWidth(2.0f);
	glDrawElements(GL_LINES, pcls_mem.get_bl_index_num(), GL_UNSIGNED_INT, (GLvoid *)0);
	mp_data.clear_add_buffer(pb_ibo);

	// draw rigid body
	// position
	glm::mat4 rb_mv_mat = glm::translate(identity_mat, glm::vec3(rbmh.x, rbmh.y, 0.0f))
		* glm::rotate(identity_mat, GLfloat(rbmh.theta), glm::vec3(0.0f, 0.0f, 1.0f))
		* glm::translate(identity_mat, glm::vec3(-rbh.x_mc, -rbh.y_mc, 0.0f));
	shader.set_uniform_matrix4f(mv_mat_id, glm::value_ptr(rb_mv_mat));
	// color
	glm::vec4 deepskyblue(0.0f, 0.74902f, 1.0f, 1.0f);
	shader.set_uniform_vec4f(color_id, glm::value_ptr(deepskyblue));
	// draw
	rigid_body_data.use();
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glDrawElements(GL_TRIANGLES, rb_elem_point_num, GL_UNSIGNED_INT, (GLvoid *)0);

	return 0;
}

bool GenerateAnimation::find_next_frame(void)
{
	double prev_time = time_rcds[cur_time_rcd_id].total_time;
	size_t frame_id;
	double diff_time;
	for (frame_id = cur_time_rcd_id + 1; frame_id < time_rcd_num; ++frame_id)
	{
		diff_time = time_rcds[frame_id].total_time - prev_time;
		if (diff_time >= min_delay_real)
		{
			delay_ani = diff_time * ani_real_ratio;
			delay_ani_100th = unsigned short int(delay_ani * 100);
			break;
		}
	}
	
	if (frame_id >= time_rcd_num)
	{
		delay_ani = 0.25; // update every 0.25s
		delay_ani_100th = 25;
		return false;
	}

	cur_time_rcd_id = frame_id;
	return true;
}
