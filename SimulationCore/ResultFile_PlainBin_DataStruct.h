#ifndef __RESULT_FILE_PLAIN_BIN_DATA_STRUCT_H__
#define __RESULT_FILE_PLAIN_BIN_DATA_STRUCT_H__

namespace ResultFile_PlainBin_DataStruct
{
	// ================= model data ================== 
	struct MeshHeader
	{
		double h, x0, xn, y0, yn;
		unsigned long long elem_x_num, elem_y_num;
	};
	struct RigidBodyHeader
	{
		unsigned long long node_num, elem_num;
		double x_mc, y_mc;
		// node coordinates, elem indices
	};
	struct MPObjectHeader
	{
		unsigned long long pcl_num;
		// x coord, y coord, vol, ...
	};
	
	// ================= time history ================= 
	struct TimeHistoryHeader
	{
		size_t substep_num;
		size_t total_substep_num;
		double current_time;
		double total_time;
	};
	struct RigidBodyMotionHeader
	{
		double x, y, theta;
		double vx, vy, v_theta;
		// resistence caused by contact
		double fx_con, fy_con, m_con;
	};
};

#endif