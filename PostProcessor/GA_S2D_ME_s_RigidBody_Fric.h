#ifndef __GA_S2D_ME_S_RIGIDBODY_FRIC_H__
#define __GA_S2D_ME_S_RIGIDBODY_FRIC_H__

#include "GenerateAnimation.h"

class GA_S2D_ME_s_RigidBody_Fric : public GenerateAnimation
{
protected: // model data
	typedef ResultFile_PlainBin_DataStruct::MeshHeader MeshHeader;
	typedef ResultFile_PlainBin_DataStruct::RigidBodyHeader RigidBodyHeader;
	typedef ResultFile_PlainBin_DataStruct::MPObjectHeader MPObjectHeader;
	
	MeshHeader mh;
	GLsizei grid_line_points_num;
	RigidBodyHeader rbh;
	GLsizei rb_elem_point_num;
	MPObjectHeader mph;

protected: // time history data
	typedef ResultFile_PlainBin_DataStruct::RigidBodyMotionHeader RigidBodyMotionHeader;
	typedef ResultFile_PlainBin_DataStruct::MPObjectHeader MPObjectHeader;

	// particle data
	double *mp_x_data, *mp_y_data, *mp_vol_data;
	SquareParticleSystem pcls_mem;

	// start position of current time record
	size_t first_time_rcd_file_pos;
	// length of each time record
	size_t time_rcd_len;

public:
	GA_S2D_ME_s_RigidBody_Fric();
	~GA_S2D_ME_s_RigidBody_Fric();

protected: // helper functions of generate()
	// assume pcl_num is the same in each time history output
	int init(const char *res_file_name) override;
	int render_frame(double xl, double xu, double yl, double yu) override;
};

#endif