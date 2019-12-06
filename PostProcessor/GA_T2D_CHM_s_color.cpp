#include "PostProcessor_pcp.h"

#include "GA_T2D_CHM_s_color.h"

GA_T2D_CHM_s_color::GA_T2D_CHM_s_color(GLsizei win_w, GLsizei win_h) :
	GenerateAnimation(win_w, win_h),
	rc_x_data(nullptr), mp_x_data(nullptr) {}
GA_T2D_CHM_s_color::~GA_T2D_CHM_s_color()
{
	if (rc_x_data)
	{
		delete[] rc_x_data;
		rc_x_data = nullptr;
	}
	if (mp_x_data)
	{
		delete[] mp_x_data;
		mp_x_data = nullptr;
	}
}

int GA_T2D_CHM_s_color::init(const char *res_file_name)
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
	res_file.read(reinterpret_cast<char *>(&mdh), sizeof(mdh));

	// bg mesh
	res_file.read(reinterpret_cast<char *>(&mh), sizeof(mh));
	// init bg grid data
	// nodes
	size_t node_num = mh.node_num;
	double *n_coords = new double[node_num * 2];
	res_file.read(reinterpret_cast<char *>(n_coords), node_num * 2 * sizeof(double));
	GLfloat *coords = new GLfloat[node_num * 3];
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		coords[n_id * 3]     = n_coords[n_id * 2];
		coords[n_id * 3 + 1] = n_coords[n_id * 2 + 1];
		coords[n_id * 3 + 2] = 0.0f;
		//std::cout << coords[n_id*3] << ", " << coords[n_id*3+1] << ", " << coords[n_id*3+2] << "\n";
	}
	bg_grid_data.init_array_buffer(coords, node_num * 3);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 3, (GLvoid *)0);
	glEnableVertexAttribArray(0);
	delete[] n_coords;
	delete[] coords;
	// elements
	size_t elem_num = mh.elem_num;
	elem_n_id_num = elem_num * 3;
	unsigned long long *e_indices = new unsigned long long[elem_n_id_num];
	res_file.read(reinterpret_cast<char *>(e_indices), elem_n_id_num * sizeof(unsigned long long));
	GLuint *indices = new GLuint[elem_n_id_num];
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		indices[e_id * 3]     = e_indices[e_id * 3];
		indices[e_id * 3 + 1] = e_indices[e_id * 3 + 1];
		indices[e_id * 3 + 2] = e_indices[e_id * 3 + 2];
		//std::cout << indices[e_id*3] << ", " << indices[e_id*3+1] << ", " << indices[e_id*3+2] << "\n";
	}
	bg_grid_data.init_elem_array_buffer(indices, elem_n_id_num);
	delete[] e_indices;
	delete[] indices;

	// rigid circle
	res_file.read(reinterpret_cast<char *>(&rch), sizeof(rch));
	rc_pcl_num = rch.pcl_num;
	if (rc_x_data)
	{
		delete[] rc_x_data;
		rc_x_data = nullptr;
	}
	if (rc_pcl_num)
	{
		rc_x_data = new double[rc_pcl_num * 3];
		rc_y_data = rc_x_data + rc_pcl_num;
		rc_vol_data = rc_y_data + rc_pcl_num;
	}

	// material point object
	res_file.read(reinterpret_cast<char *>(&mph), sizeof(mph));
	if (mp_x_data)
	{
		delete[] mp_x_data;
		mp_x_data = nullptr;
	}
	if (mph.pcl_num)
	{
		mp_x_data = new double[mph.pcl_num * 4];
		mp_y_data = mp_x_data + mph.pcl_num;
		mp_vol_data = mp_y_data + mph.pcl_num;
		mp_p_data = mp_vol_data + mph.pcl_num;
	}

	// get total number of time record
	size_t file_len;
	first_time_rcd_file_pos = res_file.tellg();
	res_file.seekg(0, SEEK_END);
	file_len = res_file.tellg();
	if (rch.pcl_num) // has rigid body
	{
		time_rcd_len = sizeof(TimeHistoryHeader)
					 + sizeof(DispConRigidCircleMotionHeader) + rc_pcl_num * 3 * sizeof(double)
					 + sizeof(MPObjectHeader) + mph.pcl_num * 4 * sizeof(double);
	}
	else
	{
		time_rcd_len = sizeof(TimeHistoryHeader)
					 + sizeof(MPObjectHeader) + mph.pcl_num * 4 * sizeof(double);
	}
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

	return 0;
}

int GA_T2D_CHM_s_color::render_frame(double xl, double xu, double yl, double yu)
{
	glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
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
	glm::vec4 white(1.0f, 1.0f, 1.0f, 1.0f);
	shader.set_uniform_vec4f(color_id, glm::value_ptr(white));
	// draw
	bg_grid_data.use();
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glDrawElements(GL_TRIANGLES, elem_n_id_num, GL_UNSIGNED_INT, (GLvoid *)0);

	// draw rigid circle
	if (rc_pcl_num)
	{
		// model/view matrix
		glm::mat4 identity_mat = glm::mat4(1.0f);
		shader.set_uniform_matrix4f(mv_mat_id, glm::value_ptr(identity_mat));
		// color
		glm::vec4 dodgerblue(0.11765f, 0.56471f, 1.0f, 1.0f);
		shader.set_uniform_vec4f(color_id, glm::value_ptr(dodgerblue));
		DispConRigidCircleMotionHeader rcmh;
		res_file.seekg(first_time_rcd_file_pos + cur_time_rcd_id * time_rcd_len + sizeof(TimeHistoryHeader), SEEK_SET);
		res_file.read(reinterpret_cast<char *>(&rcmh), sizeof(rcmh));
		res_file.read(reinterpret_cast<char *>(rc_x_data), sizeof(double) * rc_pcl_num * 3);
		rc_pcls_mem.reset();
		for (size_t pcl_id = 0; pcl_id < rc_pcl_num; ++pcl_id)
			rc_pcls_mem.add_pcl(rc_x_data[pcl_id], rc_y_data[pcl_id], rc_vol_data[pcl_id] * 0.25);
		rc_data.clear();
		rc_data.init_array_buffer(rc_pcls_mem.get_pcls(), rc_pcls_mem.get_point_num());
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 3, (GLvoid *)0);
		glEnableVertexAttribArray(0);
		rc_data.init_elem_array_buffer(rc_pcls_mem.get_indices(), rc_pcls_mem.get_index_num());
		// draw particles
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glDrawElements(GL_TRIANGLES, rc_pcls_mem.get_index_num(), GL_UNSIGNED_INT, (GLvoid *)0);
	}

	// draw material points
	// model/view matrix
	shader_color.use();
	shader_color.set_uniform_matrix4f(c_mv_mat_id, glm::value_ptr(identity_mat));
	// init object buffer
	MPObjectHeader mph;
	if (rc_pcl_num)
		res_file.seekg(first_time_rcd_file_pos + cur_time_rcd_id * time_rcd_len + sizeof(TimeHistoryHeader)
					 + sizeof(DispConRigidCircleMotionHeader) + rc_pcl_num * 3 * sizeof(double), SEEK_SET);
	else
		res_file.seekg(first_time_rcd_file_pos + cur_time_rcd_id * time_rcd_len + sizeof(TimeHistoryHeader), SEEK_SET);
	res_file.read(reinterpret_cast<char *>(&mph), sizeof(mph));
	res_file.read(reinterpret_cast<char *>(mp_x_data), sizeof(double) * mph.pcl_num * 4);
	ColorGraph::Colorf pcl_color;
	pcls_mem.reset();
	for (size_t pcl_id = 0; pcl_id < mph.pcl_num; ++pcl_id)
	{
		pcl_color = color_graph.get_color(mp_p_data[pcl_id]);
		pcls_mem.add_pcl(mp_x_data[pcl_id], mp_y_data[pcl_id], mp_vol_data[pcl_id] * 0.25, // pcl data
						 pcl_color.r, pcl_color.g, pcl_color.b); // pcl color
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
	// draw particle boundary
	// color
	//glm::vec4 red(1.0f, 0.0f, 0.0f, 1.0f);
	//shader.set_uniform_vec4f(color_id, glm::value_ptr(red));
	//GLuint pb_ibo = 0;
	//mp_data.init_add_elem_array_buffer(pb_ibo, pcls_mem.get_bl_indices(), pcls_mem.get_bl_index_num());
	//glLineWidth(2.0f);
	//glDrawElements(GL_LINES, pcls_mem.get_bl_index_num(), GL_UNSIGNED_INT, (GLvoid *)0);
	//mp_data.clear_add_buffer(pb_ibo);

	return 0;
}
