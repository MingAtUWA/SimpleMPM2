#ifndef __GA_T2D_CHM_s_color_H__
#define __GA_T2D_CHM_s_color_H__

#include "ColorGraph.h"
#include "GenerateAnimation.h"

class GA_T2D_CHM_s_color : public GenerateAnimation
{
protected:
	typedef ResultFile_PlainBin_DataStruct::ModelDataHeader ModelDataHeader;
	typedef ResultFile_PlainBin_DataStruct::BackgroundMeshHeader_T2D BackgroundMeshHeader;
	typedef ResultFile_PlainBin_DataStruct::MPObjectHeader MPObjectHeader;

	// model data
	ModelDataHeader mdh;
	BackgroundMeshHeader mh;
	MPObjectHeader mph;
	size_t elem_n_id_num;

	// time history data
	// particle data
	double *mp_x_data, *mp_y_data, *mp_vol_data, *mp_p_data;
	SquareParticleSystem pcls_mem;

	// start position of current time record
	size_t first_time_rcd_file_pos;
	// length of each time record
	size_t time_rcd_len;

	ShaderProgram shader_color;
	GLint c_mv_mat_id;
	GLint c_proj_mat_id;

	ColorGraph color_graph;

public:
	GA_T2D_CHM_s_color();
	~GA_T2D_CHM_s_color();

protected: // helper functions of generate()
	// assume pcl_num is the same in each time history output
	int init(const char *res_file_name) override;
	int init_color_graph(ColorGraph::ValueColorPair *vcps, size_t num);
	int render_frame(double xl, double xu, double yl, double yu) override;
};

#endif