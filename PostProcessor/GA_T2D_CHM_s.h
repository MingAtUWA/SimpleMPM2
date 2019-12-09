#ifndef __GA_T2D_CHM_s_H__
#define __GA_T2D_CHM_s_H__

#include "GenerateAnimation.h"

class GA_T2D_CHM_s : public GenerateAnimation
{
protected:
	typedef ResultFile_PlainBin_DataStruct::ModelDataHeader ModelDataHeader;
	typedef ResultFile_PlainBin_DataStruct::BackgroundMeshHeader_T2D BackgroundMeshHeader;
	typedef ResultFile_PlainBin_DataStruct::DispConRigidCircleHeader DispConRigidCircleHeader;
	typedef ResultFile_PlainBin_DataStruct::MPObjectHeader MPObjectHeader;
	typedef ResultFile_PlainBin_DataStruct::DispConRigidCircleMotionHeader DispConRigidCircleMotionHeader;

	// model data
	ModelDataHeader mdh;
	BackgroundMeshHeader mh;
	DispConRigidCircleHeader rch;
	MPObjectHeader mph;
	size_t elem_n_id_num;

	// time history data
	// rigid circle data
	size_t rc_pcl_num;
	double *rc_x_data, *rc_y_data, *rc_vol_data;
	SquareParticleSystem rc_pcls_mem;
	BufferObject rc_data;
	
	// particle data
	double *mp_x_data, *mp_y_data, *mp_vol_data, *mp_p_data;
	SquareParticleSystem pcls_mem;

	// start position of current time record
	size_t first_time_rcd_file_pos;
	// length of each time record
	size_t time_rcd_len;

public:
	GA_T2D_CHM_s(GLsizei win_w = 600, GLsizei win_h = 600);
	~GA_T2D_CHM_s();

protected: // helper functions of generate()
	// assume pcl_num is the same in each time history output
	int init(const char *res_file_name) override;
	int render_frame(double xl, double xu, double yl, double yu) override;
};

#endif