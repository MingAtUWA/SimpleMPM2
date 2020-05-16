#ifndef __Model_S2D_ME_s_hdf5_io_utilities_h__
#define __Model_S2D_ME_s_hdf5_io_utilities_h__

#include "ResultFile_hdf5.h"

#include "Model_S2D_ME_s.h"
//#include "ModelContainer.h"

namespace Model_S2D_ME_s_hdf5_io_utilities
{
	struct ParticleData
	{
		unsigned long long id;
		double m;
		double density;
		double ar;
		double x;
		double y;
		double vx;
		double vy;
		double s11;
		double s22;
		double s12;
		double e11;
		double e22;
		double e12;
		void from_pcl(Model_S2D_ME_s::Particle &pcl)
		{
			m = pcl.m;
			density = pcl.density;
			ar = pcl.ar;
			x = pcl.x;
			y = pcl.y;
			vx = pcl.vx;
			vy = pcl.vy;
			s11 = pcl.s11;
			s22 = pcl.s22;
			s12 = pcl.s12;
			e11 = pcl.e11;
			e22 = pcl.e22;
			e12 = pcl.e12;
		}
		void to_pcl(Model_S2D_ME_s::Particle &pcl)
		{
			pcl.m = m;
			pcl.density = density;
			pcl.ar = ar;
			pcl.x = x;
			pcl.y = y;
			pcl.vx = vx;
			pcl.vy = vy;
			pcl.s11 = s11;
			pcl.s22 = s22;
			pcl.s12 = s12;
			pcl.e11 = e11;
			pcl.e22 = e22;
			pcl.e12 = e12;
		}
	};

	inline hid_t get_pcl_dt_id(void)
	{
		hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(ParticleData));
		H5Tinsert(res, "id", HOFFSET(ParticleData, id), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "m",  HOFFSET(ParticleData, m), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "density", HOFFSET(ParticleData, density), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "ar", HOFFSET(ParticleData, ar), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "x",  HOFFSET(ParticleData, x), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "y",  HOFFSET(ParticleData, y), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "vx", HOFFSET(ParticleData, vx), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "vy", HOFFSET(ParticleData, vy), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s11", HOFFSET(ParticleData, s11), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s22", HOFFSET(ParticleData, s22), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s12", HOFFSET(ParticleData, s12), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "e11", HOFFSET(ParticleData, e11), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "e22", HOFFSET(ParticleData, e22), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "e12", HOFFSET(ParticleData, e12), H5T_NATIVE_DOUBLE);
		return res;
	}

	int output_model_data_to_hdf5_file(Model_S2D_ME_s &md, ResultFile_hdf5 &rf);
	int load_model_data_from_hdf5_file(Model_S2D_ME_s &md, ResultFile_hdf5 &rf);

	//int output_model_container_to_hdf5_file(ModelContainer &mc, ResultFile_hdf5 &rf, hid_t frame_id /* frame id */);
	//int load_model_container_from_hdf5_file(Model_T2D_ME_s &md, ResultFile_hdf5 &rf, hid_t frame_id /* frame id */);

	int output_pcl_data_to_hdf5_file(Model_S2D_ME_s &md, ResultFile_hdf5 &rf, hid_t frame_id /* frame id */);
	int load_pcl_data_from_hdf5_file(Model_S2D_ME_s &md, ResultFile_hdf5 &rf, hid_t frame_id /* frame id */);

	//int load_me_s_model_from_hdf5_file(
	//	Model_S2D_ME_s &md,
	//	const char *hdf5_name,
	//	const char *th_name,
	//	size_t frame_id
	//	);
};

#endif