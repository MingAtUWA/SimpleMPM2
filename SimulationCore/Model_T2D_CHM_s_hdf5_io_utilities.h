#ifndef __Model_T2D_CHM_s_hdf5_io_utilities_H__
#define __Model_T2D_CHM_s_hdf5_io_utilities_H__

#include "hdf5.h"
#include "ResultFile_hdf5.h"
#include "ModelContainer.h"
#include "DispConRigidCircle.h"
#include "Model_T2D_CHM_s.h"

namespace Model_T2D_CHM_s_hdf5_io_utilities
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
	void from_pcl(Model_T2D_CHM_s::Particle &pcl)
	{
		m_s = pcl.m_s;
		density_s = pcl.density_s;
		density_f = pcl.density_f;
		n = pcl.n;
		x = pcl.x;
		y = pcl.y;
		vx_s = pcl.vx_s;
		vy_s = pcl.vy_s;
		vx_f = pcl.vx_f;
		vy_f = pcl.vy_f;
		s11 = pcl.s11;
		s22 = pcl.s22;
		s12 = pcl.s12;
		p = pcl.p;
		e11 = pcl.e11;
		e22 = pcl.e22;
		e12 = pcl.e12;
	}
	void to_pcl(Model_T2D_CHM_s::Particle &pcl)
	{
		pcl.m_s = m_s;
		pcl.density_s = density_s;
		pcl.density_f = density_f;
		pcl.n = n;
		pcl.x = x;
		pcl.y = y;
		pcl.vx_s = vx_s;
		pcl.vy_s = vy_s;
		pcl.vx_f = vx_f;
		pcl.vy_f = vy_f;
		pcl.s11 = s11;
		pcl.s22 = s22;
		pcl.s12 = s12;
		pcl.p = p;
		pcl.e11 = e11;
		pcl.e22 = e22;
		pcl.e12 = e12;
	}
};

inline hid_t get_pcl_dt_id(void)
{
	hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(ParticleData));
	H5Tinsert(res, "id", HOFFSET(ParticleData, id), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "m_s", HOFFSET(ParticleData, m_s), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "density_s", HOFFSET(ParticleData, density_s), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "density_f", HOFFSET(ParticleData, density_f), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "n", HOFFSET(ParticleData, n), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "x", HOFFSET(ParticleData, x), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "y", HOFFSET(ParticleData, y), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "vx_s", HOFFSET(ParticleData, vx_s), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "vy_s", HOFFSET(ParticleData, vy_s), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "vx_f", HOFFSET(ParticleData, vx_f), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "vy_f", HOFFSET(ParticleData, vy_f), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "s11", HOFFSET(ParticleData, s11), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "s22", HOFFSET(ParticleData, s22), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "s12", HOFFSET(ParticleData, s12), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "p", HOFFSET(ParticleData, p), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "e11", HOFFSET(ParticleData, e11), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "e22", HOFFSET(ParticleData, e22), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "e12", HOFFSET(ParticleData, e12), H5T_NATIVE_DOUBLE);
	return res;
}


int output_model_data_to_hdf5_file(Model_T2D_CHM_s &md, ResultFile_hdf5 &rf);
int load_model_data_from_hdf5_file(Model_T2D_CHM_s &md, ResultFile_hdf5 &rf);

int output_bcs_to_hdf5_file(Model_T2D_CHM_s& md, ResultFile_hdf5& rf);
int load_bcs_to_hdf5_file(Model_T2D_CHM_s& md, ResultFile_hdf5& rf);

int output_model_container_to_hdf5_file(ModelContainer &mc, ResultFile_hdf5 &rf, hid_t frame_id /* frame id */);
int load_model_container_from_hdf5_file(Model_T2D_CHM_s &md, ResultFile_hdf5 &rf, hid_t frame_id /* frame id */);

int output_rigid_ciricle_to_hdf5_file(DispConRigidCircle &rc, ResultFile_hdf5 &rf, hid_t frame_id /* frame id */);
int load_rigid_ciricle_to_hdf5_file(DispConRigidCircle &rc, ResultFile_hdf5 &rf, hid_t frame_id /* frame id */);

int output_pcl_data_to_hdf5_file(Model_T2D_CHM_s &md, ResultFile_hdf5 &rf, hid_t frame_id /* frame id */);
int load_pcl_data_from_hdf5_file(Model_T2D_CHM_s &md, ResultFile_hdf5 &rf, hid_t frame_id /* frame id */);

int load_chm_s_model_from_hdf5_file(Model_T2D_CHM_s &md, const char *hdf5_name, const char *th_name, size_t frame_id);

};

#endif