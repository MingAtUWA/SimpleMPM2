#ifndef __Model_T2D_CHM_s_hdf5_io_utilities_H__
#define __Model_T2D_CHM_s_hdf5_io_utilities_H__

#include "hdf5.h"
#include "ResultFile_hdf5.h"
#include "ModelContainer.h"
#include "DispConRigidCircle.h"
#include "Model_T2D_CHM_s.h"

int output_model_data_to_hdf5_file(Model_T2D_CHM_s &md, ResultFile_hdf5 &rf);
int load_model_data_from_hdf5_file(Model_T2D_CHM_s &md, ResultFile_hdf5 &rf);

int output_model_container_to_hdf5_file(ModelContainer &mc, ResultFile_hdf5 &rf, hid_t frame_id /* frame id */);
int load_model_container_from_hdf5_file(Model_T2D_CHM_s &md, ResultFile_hdf5 &rf, hid_t frame_id /* frame id */);

int output_rigid_ciricle_to_hdf5_file(DispConRigidCircle &rc, ResultFile_hdf5 &rf, hid_t frame_id /* frame id */);
int load_rigid_ciricle_to_hdf5_file(DispConRigidCircle &rc, ResultFile_hdf5 &rf, hid_t frame_id /* frame id */);

int ouput_pcl_data_to_hdf5_file(Model_T2D_CHM_s &md, ResultFile_hdf5 &rf, hid_t frame_id /* frame id */);
int load_pcl_data_from_hdf5_file(Model_T2D_CHM_s &md, ResultFile_hdf5 &rf, hid_t frame_id /* frame id */);

int load_chm_s_model_from_hdf5_file(Model_T2D_CHM_s &md, const char *hdf5_name, const char *th_name, size_t frame_id);

#endif