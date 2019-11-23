#ifndef __RESULT_FILE_PLAIN_BIN_DATA_STRUCT_H__
#define __RESULT_FILE_PLAIN_BIN_DATA_STRUCT_H__

namespace ResultFile_PlainBin_DataStruct
{
	// ================= model data ==================
	struct ModelDataHeader
	{
		char tag[32];
		double current_time;
		double total_time;
		inline void init(void)
		{
			memset(tag, 0, sizeof(tag) / sizeof(tag[0]));
			strcpy(tag, "MD");
		}
	};

	struct BackgroundMeshHeader
	{
		char tag[32];
		double h, x0, xn, y0, yn;
		unsigned long long elem_x_num, elem_y_num;
		inline void init(void)
		{
			memset(tag, 0, sizeof(tag) / sizeof(tag[0]));
			strcpy(tag, "BgMeshS2D");
		}
	};
	struct BackgroundMeshHeader_T2D
	{
		char tag[32];
		unsigned long long node_num, elem_num;
		inline void init(void)
		{
			memset(tag, 0, sizeof(tag) / sizeof(tag[0]));
			strcpy(tag, "BgMeshT2D");
		}
	};

	struct RigidBodyHeader
	{
		char tag[32];
		unsigned long long node_num, elem_num;
		double x_mc, y_mc;
		// node coordinates, elem indices ...
		inline void init(void)
		{
			memset(tag, 0, sizeof(tag) / sizeof(tag[0]));
			strcpy(tag, "RigidObj");
		}
	};
	
	// ================= time history ================= 
	struct TimeHistoryHeader
	{
		char tag[32];
		size_t substep_num;
		size_t total_substep_num;
		double current_time;
		double total_time;
		inline void init(void)
		{
			memset(tag, 0, sizeof(tag) / sizeof(tag[0]));
			strcpy(tag, "TH");
		}
	};

	struct RigidBodyMotionHeader
	{
		char tag[32];
		double x, y, theta;
		double vx, vy, v_theta;
		// resistence caused by contact
		double fx_con, fy_con, m_con;
		inline void init(void)
		{
			memset(tag, 0, sizeof(tag) / sizeof(tag[0]));
			strcpy(tag, "RigObjMotion");
		}
	};

	struct MPObjectHeader
	{
		char tag[32];
		unsigned long long pcl_num;
		unsigned long long fld_num;
		// x coord, y coord, vol, ...
		inline void init(void)
		{
			memset(tag, 0, sizeof(tag) / sizeof(tag[0]));
			strcpy(tag, "MPObj");
		}
	};

	struct MeshObjectHeader_2D4R
	{
		char tag[32];
		unsigned long long node_num;
		unsigned long long node_fld_num;
		unsigned long long elem_num;
		unsigned long long elem_fld_num;
		inline void init(void)
		{
			memset(tag, 0, sizeof(tag) / sizeof(tag[0]));
			strcpy(tag, "MeshObj");
		}
	};
};

#endif