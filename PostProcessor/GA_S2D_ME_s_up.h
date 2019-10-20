#ifndef __GA_S2D_ME_s_up_H__
#define __GA_S2D_ME_s_up_H__

#include "GenerateAnimation.h"

class GA_S2D_ME_s_up : public GenerateAnimation
{
protected: // model data
	typedef ResultFile_PlainBin_DataStruct::ModelDataHeader ModelDataHeader;
	typedef ResultFile_PlainBin_DataStruct::BackgroundMeshHeader BackgroundMeshHeader;
	typedef ResultFile_PlainBin_DataStruct::MPObjectHeader MPObjectHeader;
	
	ModelDataHeader mdh;
	BackgroundMeshHeader mh;
	GLsizei grid_line_points_num;
	MPObjectHeader mph;

protected: // time history data
	// particle data
	double *mp_x_data, *mp_y_data, *mp_vol_data;
	double *mp_p_data, *mp_s11_data, *mp_s22_data, *mp_s12_data;
	SquareParticleSystem pcls_mem;

	// start position of current time record
	size_t first_time_rcd_file_pos;
	// length of each time record
	size_t time_rcd_len;

public:
	GA_S2D_ME_s_up();
	~GA_S2D_ME_s_up();

protected: // helper functions of generate()
	// assume pcl_num is the same in each time history output
	int init(const char *res_file_name) override;
	int render_frame(double xl, double xu, double yl, double yu) override;
};

#endif