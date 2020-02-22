#ifndef __Model_T2D_fluid_hdf5_io_utilities_H__
#define __Model_T2D_fluid_hdf5_io_utilities_H__

#include "hdf5.h"
#include "ResultFile_hdf5.h"
#include "Model_T2D_fluid.h"

namespace Model_T2D_fluid_hdf5_io_utilities
{
struct ParticleData
{
	unsigned long long id;
	double m;
	double density;
	double x;
	double y;
	double vx;
	double vy;
	double t11;
	double t22;
	double t12;
	double p;
	void from_pcl(Model_T2D_fluid::Particle &pcl)
	{
		m = pcl.m;
		density = pcl.density;
		x = pcl.x;
		y = pcl.y;
		vx = pcl.vx;
		vy = pcl.vy;
		t11 = pcl.t11;
		t22 = pcl.t22;
		t12 = pcl.t12;
		p = pcl.p;
	}
	void to_pcl(Model_T2D_fluid::Particle &pcl)
	{
		pcl.m = m;
		pcl.density = density;
		pcl.x = x;
		pcl.y = y;
		pcl.vx = vx;
		pcl.vy = vy;
		pcl.t11 = t11;
		pcl.t22 = t22;
		pcl.t12 = t12;
		pcl.p = p;
	}
};

inline hid_t get_pcl_dt_id(void)
{
	hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(ParticleData));
	H5Tinsert(res, "id", HOFFSET(ParticleData, id), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "m", HOFFSET(ParticleData, m), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "density", HOFFSET(ParticleData, density), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "x", HOFFSET(ParticleData, x), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "y", HOFFSET(ParticleData, y), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "vx", HOFFSET(ParticleData, vx), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "vy", HOFFSET(ParticleData, vy), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "t11", HOFFSET(ParticleData, t11), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "t22", HOFFSET(ParticleData, t22), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "t12", HOFFSET(ParticleData, t12), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "p", HOFFSET(ParticleData, p), H5T_NATIVE_DOUBLE);
	return res;
}

int output_model_data_to_hdf5_file(Model_T2D_fluid &md, ResultFile_hdf5 &rf);
int load_model_data_from_hdf5_file(Model_T2D_fluid &md, ResultFile_hdf5 &rf);

int output_pcl_data_to_hdf5_file(Model_T2D_fluid &md, ResultFile_hdf5 &rf, hid_t frame_id /* frame id */);
int load_pcl_data_from_hdf5_file(Model_T2D_fluid &md, ResultFile_hdf5 &rf, hid_t frame_id /* frame id */);

//int load_fluid_model_from_hdf5_file(
//	Model_T2D_fluid &md,
//	const char *hdf5_name,
//	const char *th_name,
//	size_t frame_id
//	);
};

#endif