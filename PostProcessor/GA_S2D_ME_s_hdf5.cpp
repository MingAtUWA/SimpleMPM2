#include "PostProcessor_pcp.h"

#include "ResultFile_hdf5_DataStruct.h"

#include "GA_S2D_ME_s_hdf5.h"

using namespace ResultFile_hdf5_DataStruct;
using namespace Model_S2D_ME_s_hdf5_io_utilities;

GA_S2D_ME_s_hdf5::GA_S2D_ME_s_hdf5(GLsizei win_w, GLsizei win_h) :
	GenerateAnimation(win_w, win_h),
	cbar_pt_data(nullptr),
	cbar_id_data(nullptr),
	th_id(-1),
	cbar_is_init(false) {}

GA_S2D_ME_s_hdf5::~GA_S2D_ME_s_hdf5()
{
	if (cbar_pt_data)
	{
		delete[] cbar_pt_data;
		cbar_pt_data = nullptr;
	}
	if (cbar_id_data)
	{
		delete[] cbar_id_data;
		cbar_id_data = nullptr;
	}
	if (th_id >= 0)
	{
		rf_hdf5.close_group(th_id);
		th_id = -1;
	}
}

int GA_S2D_ME_s_hdf5::_init(const char *file_name, const char *th_name)
{
	rf_hdf5.open(file_name);
	res_file_id = rf_hdf5.get_file_id();
	if (res_file_id < 0) return -1;

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

	shader_color.init_from_file("..\\..\\Asset\\multicolor_vshader.txt",
								"..\\..\\Asset\\multicolor_fshader.txt");
	shader_color.use();
	c_mv_mat_id = shader_color.init_uniform("mv_mat");
	if (c_mv_mat_id < 0)
		return -2;
	c_proj_mat_id = shader_color.init_uniform("proj_mat");
	if (c_proj_mat_id < 0)
		return -2;

	// model data
	model_data_id = rf_hdf5.get_model_data_grp_id();

	// bg mesh
	double x0, y0, hx, hy;
	size_t elem_x_num, elem_y_num;
	hid_t bg_mesh_id = rf_hdf5.open_group(model_data_id, "BackgroundMesh");
	rf_hdf5.read_attribute(bg_mesh_id, "x0", x0);
	rf_hdf5.read_attribute(bg_mesh_id, "y0", y0);
	rf_hdf5.read_attribute(bg_mesh_id, "elem_x_num", elem_x_num);
	rf_hdf5.read_attribute(bg_mesh_id, "elem_y_num", elem_y_num);
	rf_hdf5.read_attribute(bg_mesh_id, "hx", hx);
	rf_hdf5.read_attribute(bg_mesh_id, "hy", hy);
	rf_hdf5.close_group(bg_mesh_id);
	double xn, yn;
	xn = x0 + double(elem_x_num) * hx;
	yn = y0 + double(elem_y_num) * hy;

	size_t node_x_num = elem_x_num + 1;
	size_t node_y_num = elem_y_num + 1;
	elem_n_id_num = (node_x_num + node_y_num) * 2;
	GLfloat *node_coords = new GLfloat[elem_n_id_num * 3];
	GLuint *elem_indices = new GLuint[elem_n_id_num];
	for (size_t n_id = 0; n_id < node_x_num; ++n_id)
	{
		node_coords[6 * n_id] = x0 + GLfloat(n_id) * hx;
		node_coords[6 * n_id + 1] = y0;
		node_coords[6 * n_id + 2] = 0.0f;
		node_coords[6 * n_id + 3] = x0 + GLfloat(n_id) * hx;
		node_coords[6 * n_id + 4] = yn;
		node_coords[6 * n_id + 5] = 0.0f;
		elem_indices[2 * n_id] = 2 * n_id;
		elem_indices[2 * n_id + 1] = 2 * n_id + 1;
	}
	size_t node_num_tmp = node_x_num + node_y_num;
	for (size_t n_id = node_x_num; n_id < node_num_tmp; ++n_id)
	{
		node_coords[6 * n_id] = x0;
		node_coords[6 * n_id + 1] = y0 + GLfloat(n_id - node_x_num) * hy;
		node_coords[6 * n_id + 2] = 0.0f;
		node_coords[6 * n_id + 3] = xn;
		node_coords[6 * n_id + 4] = y0 + GLfloat(n_id - node_x_num) * hy;
		node_coords[6 * n_id + 5] = 0.0f;
		elem_indices[2 * n_id] = 2 * n_id;
		elem_indices[2 * n_id + 1] = 2 * n_id + 1;
	}
	bg_grid_data.init_array_buffer(node_coords, elem_n_id_num * 3);
	delete[] node_coords;
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 3, (GLvoid *)0);
	glEnableVertexAttribArray(0);
	bg_grid_data.init_elem_array_buffer(elem_indices, elem_n_id_num);
	delete[] elem_indices;

	// frame number
	time_history_id = rf_hdf5.get_time_history_grp_id();
	th_id = rf_hdf5.open_group(time_history_id, th_name);
	rf_hdf5.read_attribute(th_id, "output_num", time_rcd_num);
	time_rcds = new TimeHistoryHeader[time_rcd_num];
	hid_t frame_id;
	char frame_name[30];
	for (size_t tr_id = 0; tr_id < time_rcd_num; ++tr_id)
	{
		snprintf(frame_name, 30, "frame_%zu", tr_id);
		frame_id = rf_hdf5.open_group(th_id, frame_name);
		TimeHistoryHeader &thh = time_rcds[tr_id];
		rf_hdf5.read_attribute(frame_id, "current_time", thh.current_time);
		rf_hdf5.read_attribute(frame_id, "total_time", thh.total_time);
		rf_hdf5.read_attribute(frame_id, "substep_num", thh.substep_num);
		rf_hdf5.read_attribute(frame_id, "total_substep_num", thh.total_substep_num);
		rf_hdf5.close_group(frame_id);
	}

	// color bar
	if (cbar_is_init)
	{
		cbar_data.init_array_buffer(cbar_pt_data, color_num * 24);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 6, (GLvoid *)0);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 6, (GLvoid *)(3 * sizeof(GLfloat)));
		glEnableVertexAttribArray(1);
		cbar_data.init_elem_array_buffer(cbar_id_data, color_num * 6);
		delete[] cbar_pt_data;
		cbar_pt_data = nullptr;
		delete[] cbar_id_data;
		cbar_id_data = nullptr;
	}

	return 0;
}

int GA_S2D_ME_s_hdf5::init_color_graph(
	double bar_xpos, double bar_ypos, double bar_wid, double bar_len,
	double lower, double upper, ColorGraph::Colori *colors, size_t num,
	bool apply_out_of_bound_color)
{
	int res = color_graph.init(lower, upper, colors, num, apply_out_of_bound_color);

	color_num = num;

	GLfloat y_inv = bar_len / GLfloat(num);
	cbar_pt_data = new GLfloat[num * 24];
	cbar_id_data = new GLuint[num * 6];
	for (size_t i = 0; i < num; ++i)
	{
		ColorGraph::Colori &co = colors[i];
		cbar_pt_data[24 * i + 0] = bar_xpos;
		cbar_pt_data[24 * i + 1] = bar_ypos + GLfloat(i) * y_inv;
		cbar_pt_data[24 * i + 2] = 0.0f;
		cbar_pt_data[24 * i + 3] = co.r / 255.0;
		cbar_pt_data[24 * i + 4] = co.g / 255.0;
		cbar_pt_data[24 * i + 5] = co.b / 255.0;
		cbar_pt_data[24 * i + 6] = bar_xpos + bar_wid;
		cbar_pt_data[24 * i + 7] = bar_ypos + GLfloat(i) * y_inv;
		cbar_pt_data[24 * i + 8] = 0.0f;
		cbar_pt_data[24 * i + 9] = co.r / 255.0;
		cbar_pt_data[24 * i + 10] = co.g / 255.0;
		cbar_pt_data[24 * i + 11] = co.b / 255.0;
		cbar_pt_data[24 * i + 12] = bar_xpos + bar_wid;
		cbar_pt_data[24 * i + 13] = bar_ypos + GLfloat(i + 1) * y_inv;
		cbar_pt_data[24 * i + 14] = 0.0f;
		cbar_pt_data[24 * i + 15] = co.r / 255.0;
		cbar_pt_data[24 * i + 16] = co.g / 255.0;
		cbar_pt_data[24 * i + 17] = co.b / 255.0;
		cbar_pt_data[24 * i + 18] = bar_xpos;
		cbar_pt_data[24 * i + 19] = bar_ypos + GLfloat(i + 1) * y_inv;
		cbar_pt_data[24 * i + 20] = 0.0f;
		cbar_pt_data[24 * i + 21] = co.r / 255.0;
		cbar_pt_data[24 * i + 22] = co.g / 255.0;
		cbar_pt_data[24 * i + 23] = co.b / 255.0;
		cbar_id_data[6 * i + 0] = 4 * i;
		cbar_id_data[6 * i + 1] = 4 * i + 1;
		cbar_id_data[6 * i + 2] = 4 * i + 2;
		cbar_id_data[6 * i + 3] = 4 * i;
		cbar_id_data[6 * i + 4] = 4 * i + 2;
		cbar_id_data[6 * i + 5] = 4 * i + 3;
	}

	cbar_is_init = true;

	return res;
}

int GA_S2D_ME_s_hdf5::render_frame(double xl, double xu, double yl, double yu)
{
	//glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f); // white
	glClear(GL_COLOR_BUFFER_BIT);

	// set display scope
	glm::mat4 proj_mat = glm::ortho(xl, xu, yl, yu); // to be improved
	shader.use();
	shader.set_uniform_matrix4f(proj_mat_id, glm::value_ptr(proj_mat));
	shader_color.use();
	shader_color.set_uniform_matrix4f(c_proj_mat_id, glm::value_ptr(proj_mat));

	// draw bg grid
	// model/view matrix
	shader.use();
	glm::mat4 identity_mat = glm::mat4(1.0f);
	shader.set_uniform_matrix4f(mv_mat_id, glm::value_ptr(identity_mat));
	// color
	//glm::vec4 white(1.0f, 1.0f, 1.0f, 1.0f);
	glm::vec4 grey(0.3f, 0.3f, 0.3f, 1.0f);
	shader.set_uniform_vec4f(color_id, glm::value_ptr(grey));
	// draw
	bg_grid_data.use();
	glLineWidth(1.5f);
	glDrawElements(GL_LINES, elem_n_id_num, GL_UNSIGNED_INT, (GLvoid *)0);

	char frame_name[30];
	snprintf(frame_name, 30, "frame_%zu", cur_time_rcd_id);
	hid_t frame_id = rf_hdf5.open_group(th_id, frame_name);

	// draw material points
	// model/view matrix
	shader_color.use();
	shader_color.set_uniform_matrix4f(c_mv_mat_id, glm::value_ptr(identity_mat));
	// init object buffer
	size_t pcl_num;
	hid_t pcl_data_id = rf_hdf5.open_dataset(frame_id, "ParticleData");
	rf_hdf5.read_attribute(pcl_data_id, "pcl_num", pcl_num);
	rf_hdf5.close_dataset(pcl_data_id);
	pcls_data.resize(pcl_num);
	ParticleData *pds = pcls_data.get_mem();
	hid_t pcl_dt_id = get_pcl_dt_id();
	rf_hdf5.read_dataset(
		frame_id,
		"ParticleData",
		pcl_num,
		pds,
		pcl_dt_id
	);
	H5Tclose(pcl_dt_id);
	ColorGraph::Colorf pcl_color;
	pcls_mem.reset();
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		ParticleData &pd = pds[pcl_id];
		//pcl_color = color_graph.get_color(pd.p);
		pcl_color = color_graph.get_color(pd.s22);
		pcls_mem.add_pcl(
			pd.x,
			pd.y,
			pd.m / pd.density * 0.25, // vol
			pcl_color.r, // color
			pcl_color.g,
			pcl_color.b
		);
	}
	mp_data.clear();
	mp_data.init_array_buffer(pcls_mem.get_coord_and_color(), pcls_mem.get_coord_and_color_size());
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 6, (GLvoid *)0);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 6, (GLvoid *)(sizeof(GLfloat) * 3));
	glEnableVertexAttribArray(1);
	mp_data.init_elem_array_buffer(pcls_mem.get_indices(), pcls_mem.get_index_num());
	// draw particles
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glDrawElements(GL_TRIANGLES, pcls_mem.get_index_num(), GL_UNSIGNED_INT, (GLvoid *)0);

	rf_hdf5.close_group(frame_id);

	// color bar
	if (cbar_is_init)
	{
		// adjust viewport
		glViewport(0, 0, win_width, win_height);
		// mvp mat
		glm::mat4 win_proj_mat = glm::ortho(0.0f, GLfloat(win_width), 0.0f, GLfloat(win_height));
		shader_color.set_uniform_matrix4f(c_proj_mat_id, glm::value_ptr(win_proj_mat));
		// draw particles
		cbar_data.use();
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glDrawElements(GL_TRIANGLES, color_num * 6, GL_UNSIGNED_INT, (GLvoid *)0);
		// recover viewport
		glViewport(md_vp_x, md_vp_y, md_vp_width, md_vp_height);
	}

	return 0;
}

int GA_S2D_ME_s_hdf5::generate(double ani_time, double xl, double xu, double yl, double yu,
	const char *file_name, const char *th_name, const char *gif_name)
{
	th_name_str = th_name;
	GenerateAnimation::generate(ani_time, xl, xu, yl, yu, file_name, gif_name);
	return 0;
}
