#include "PostProcessor_pcp.h"

#include "GA_T2D_CHM_s_hdf5.h"

using namespace GA_T2D_CHM_s_hdf5_internal;

GA_T2D_CHM_s_hdf5::GA_T2D_CHM_s_hdf5(GLsizei win_w, GLsizei win_h) :
	GenerateAnimation(win_w, win_h), th_id(-1),
	cbar_pt_data(nullptr), cbar_id_data(nullptr), cbar_is_init(false) {}

GA_T2D_CHM_s_hdf5::~GA_T2D_CHM_s_hdf5()
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
	if (th_id >= 0) rf_hdf5.close_group(th_id);
}

int GA_T2D_CHM_s_hdf5::_init(const char *file_name, const char *th_name)
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
	hid_t bg_mesh_id = rf_hdf5.open_group(model_data_id, "BackgroundMesh");
	size_t node_num, elem_num;
	rf_hdf5.read_attribute(bg_mesh_id, "node_num", node_num);
	rf_hdf5.read_attribute(bg_mesh_id, "element_num", elem_num);
	// nodes
	struct NodeData
	{
		unsigned long long id;
		double x;
		double y;
	} *nodes_data;
	nodes_data = new NodeData[node_num];
	rf_hdf5.read_dataset(bg_mesh_id, "NodeCoordinate", node_num,
						 size_t(sizeof(NodeData)/sizeof(double)),
						 (double *)nodes_data);
	GLfloat *coords = new GLfloat[node_num * 3];
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		NodeData &nd = nodes_data[n_id];
		coords[n_id * 3 + 0] = GLfloat(nd.x);
		coords[n_id * 3 + 1] = GLfloat(nd.y);
		coords[n_id * 3 + 2] = 0.0f;
	}
	bg_grid_data.init_array_buffer(coords, node_num * 3);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 3, (GLvoid *)0);
	glEnableVertexAttribArray(0);
	delete[] coords;
	delete[] nodes_data;
	// elements
	struct ElemData
	{
		unsigned long long id;
		unsigned long long n1;
		unsigned long long n2;
		unsigned long long n3;
	} *elems_data;
	elems_data = new ElemData[elem_num];
	elem_n_id_num = elem_num * 3;
	rf_hdf5.read_dataset(bg_mesh_id, "ElementTopology", elem_num,
						 size_t(sizeof(ElemData)/sizeof(unsigned long long)),
						 (unsigned long long *)elems_data);
	GLuint *indices = new GLuint[elem_num * 3];
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		ElemData &ed = elems_data[e_id];
		indices[e_id * 3 + 0] = GLuint(ed.n1);
		indices[e_id * 3 + 1] = GLuint(ed.n2);
		indices[e_id * 3 + 2] = GLuint(ed.n3);
	}
	bg_grid_data.init_elem_array_buffer(indices, elem_num * 3);
	delete[] indices;
	delete[] elems_data;
	rf_hdf5.close_group(bg_mesh_id);

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

int GA_T2D_CHM_s_hdf5::init_color_graph(
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
		cbar_pt_data[24 * i + 13] = bar_ypos + GLfloat(i+1) * y_inv;
		cbar_pt_data[24 * i + 14] = 0.0f;
		cbar_pt_data[24 * i + 15] = co.r / 255.0;
		cbar_pt_data[24 * i + 16] = co.g / 255.0;
		cbar_pt_data[24 * i + 17] = co.b / 255.0;
		cbar_pt_data[24 * i + 18] = bar_xpos;
		cbar_pt_data[24 * i + 19] = bar_ypos + GLfloat(i+1) * y_inv;
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

int GA_T2D_CHM_s_hdf5::render_frame(double xl, double xu, double yl, double yu)
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
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glDrawElements(GL_TRIANGLES, elem_n_id_num, GL_UNSIGNED_INT, (GLvoid *)0);

	char frame_name[30];
	snprintf(frame_name, 30, "frame_%zu", cur_time_rcd_id);
	hid_t frame_id = rf_hdf5.open_group(th_id, frame_name);

	// draw rigid circle
	if (rf_hdf5.has_group(frame_id, "RigidBody"))
	{
		hid_t rb_id = rf_hdf5.open_group(frame_id, "RigidBody");

		// model/view matrix
		glm::mat4 identity_mat = glm::mat4(1.0f);
		shader.set_uniform_matrix4f(mv_mat_id, glm::value_ptr(identity_mat));
		// color
		glm::vec4 dodgerblue(0.11765f, 0.56471f, 1.0f, 1.0f);
		shader.set_uniform_vec4f(color_id, glm::value_ptr(dodgerblue));
		
		double cen_x, cen_y, theta;
		size_t rb_pcl_num;
		rf_hdf5.read_attribute(rb_id, "cen_x",  cen_x);
		rf_hdf5.read_attribute(rb_id, "cen_y", cen_y);
		rf_hdf5.read_attribute(rb_id, "theta", theta);
		rf_hdf5.read_attribute(rb_id, "pcl_num", rb_pcl_num);
		
		// read rigid body particle data
		rb_pcls_data.resize(rb_pcl_num);
		RigidBodyParticleData *rbpds = rb_pcls_data.get_mem();
		rf_hdf5.read_dataset(rb_id, "ParticleData", rb_pcl_num, 3, (double *)rbpds);

		double pcl_x, pcl_y;
		for (size_t pcl_id = 0; pcl_id < rb_pcl_num; ++pcl_id)
		{
			RigidBodyParticleData &rbpd = rbpds[pcl_id];
			pcl_x = cen_x + rbpd.xr * cos(theta) + rbpd.yr * -sin(theta);
			pcl_y = cen_y + rbpd.xr * sin(theta) + rbpd.yr *  cos(theta);
			rc_pcls_mem.add_pcl(pcl_x, pcl_y, rbpd.vol * 0.25);
		}

		rf_hdf5.close_group(rb_id);

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
	size_t pcl_num;
	hid_t pcl_data_id = rf_hdf5.open_dataset(frame_id, "ParticleData");
	rf_hdf5.read_attribute(pcl_data_id, "pcl_num", pcl_num);
	rf_hdf5.close_dataset(pcl_data_id);
	pcls_data.resize(pcl_num);
	ParticleData *pds = pcls_data.get_mem();
	rf_hdf5.read_dataset(frame_id, "ParticleData", pcl_num,
						 size_t(sizeof(ParticleData)/sizeof(double)),
						 (double *)pds);
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
			pd.m_s/((1.0-pd.n)*pd.density_s) * 0.25, // vol
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

int GA_T2D_CHM_s_hdf5::generate(double ani_time, double xl, double xu, double yl, double yu,
	const char *file_name, const char *th_name, const char *gif_name)
{
	th_name_str = th_name;
	GenerateAnimation::generate(ani_time, xl, xu, yl, yu, file_name, gif_name);
	return 0;
}
