#ifndef __GA_T2D_CHM_s_hdf5_H__
#define __GA_T2D_CHM_s_hdf5_H__

#include "ColorGraph.h"
#include "ColorfulSquareParticleSystem.h"

#include "ItemArrayFast.hpp"
#include "ResultFile_hdf5.h"

#include "GenerateAnimation.h"

namespace GA_T2D_CHM_s_hdf5_internal
{
	struct ParticleData
	{
		unsigned long long id;
		double m_s;
		double density_s;
		double density_f;
		double n;
		double x;
		double y;
		double vx_s;
		double vy_s;
		double vx_f;
		double vy_f;
		double s11;
		double s22;
		double s12;
		double p;
		double e11;
		double e22;
		double e12;
	};
	struct RigidBodyParticleData
	{
		double xr;
		double yr;
		double vol;
	};
};

class GA_T2D_CHM_s_hdf5 : public GenerateAnimation
{
protected:
	std::string th_name_str;
	ResultFile_hdf5 rf_hdf5;
	hid_t res_file_id;
	hid_t model_data_id;
	hid_t time_history_id;
	hid_t th_id;

	// model data
	size_t elem_n_id_num;

	// time history data
	// rigid circle data
	MemoryUtilities::ItemArrayFast<GA_T2D_CHM_s_hdf5_internal::RigidBodyParticleData> rb_pcls_data;
	SquareParticleSystem rc_pcls_mem;
	BufferObject rc_data;

	// particle data
	MemoryUtilities::ItemArrayFast<GA_T2D_CHM_s_hdf5_internal::ParticleData> pcls_data;
	ColorfulSquareParticleSystem pcls_mem;

	// shader for colorful point
	ShaderProgram shader_color;
	GLint c_mv_mat_id;
	GLint c_proj_mat_id;

	// color bar
	ColorGraph color_graph;
	size_t color_num;
	GLfloat *cbar_pt_data;
	GLuint *cbar_id_data;
	bool cbar_is_init;
	BufferObject cbar_data;

public:
	GA_T2D_CHM_s_hdf5(GLsizei win_w = 600, GLsizei win_h = 600);
	~GA_T2D_CHM_s_hdf5();
	inline int init_color_graph(double lower, double upper, ColorGraph::Colori *colors, size_t num,
								bool apply_out_of_bound_color = true)
	{
		return color_graph.init(lower, upper, colors, num, apply_out_of_bound_color);
	}
	int init_color_graph(double bar_xpos, double bar_ypos, double bar_wid,  double bar_len,
						 double lower, double upper, ColorGraph::Colori *colors, size_t num,
						 bool apply_out_of_bound_color = true);
	
	int generate(double ani_time, double xl, double xu, double yl, double yu,
		const char *file_name, const char *th_name, const char *gif_name);

protected: // helper functions of generate()
	// assume pcl_num is the same in each time history output
	int init(const char *file_name) override { _init(file_name, th_name_str.c_str()); return 0; }
	int _init(const char *file_name, const char *th_name);
	int render_frame(double xl, double xu, double yl, double yu) override;
};

#endif