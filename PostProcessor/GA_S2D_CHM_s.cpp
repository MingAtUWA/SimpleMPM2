#include "PostProcessor_pcp.h"

#include "GA_S2D_CHM_s.h"

GA_S2D_CHM_s::GA_S2D_CHM_s() : mp_x_data(nullptr) {}

GA_S2D_CHM_s::~GA_S2D_CHM_s()
{
	if (mp_x_data)
	{
		delete[] mp_x_data;
		mp_x_data = nullptr;
	}
}

int GA_S2D_CHM_s::init(const char *res_file_name)
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

	// model data
	res_file.read(reinterpret_cast<char *>(&mdh), sizeof(mdh));

	// bg mesh
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
	
	// material point object
	res_file.read(reinterpret_cast<char *>(&mph), sizeof(mph));
	if (!mp_x_data)
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
	time_rcd_len = sizeof(TimeHistoryHeader) + sizeof(MPObjectHeader) + mph.pcl_num * 4 * sizeof(double);
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

int GA_S2D_CHM_s::render_frame(double xl, double xu, double yl, double yu)
{
	glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);

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
	MPObjectHeader mph;
	res_file.seekg(first_time_rcd_file_pos + cur_time_rcd_id * time_rcd_len + sizeof(TimeHistoryHeader), SEEK_SET);
	res_file.read(reinterpret_cast<char *>(&mph), sizeof(mph));
	res_file.read(reinterpret_cast<char *>(mp_x_data), sizeof(double) * mph.pcl_num * 4);
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

	return 0;
}
