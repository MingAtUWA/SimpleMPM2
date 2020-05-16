#include "SimulationCore_pcp.h"

#include "Model_S2D_ME_s_hdf5_io_utilities.h"

namespace Model_S2D_ME_s_hdf5_io_utilities
{
	int output_model_data_to_hdf5_file(Model_S2D_ME_s &md, ResultFile_hdf5 &rf)
	{
		hid_t md_id = rf.get_model_data_grp_id();
		hid_t bg_mesh_id = rf.create_group(md_id, "BackgroundMesh");

		// bg mesh attributes
		const char *bg_mesh_type = "S2D";
		rf.write_attribute(bg_mesh_id, "type", strlen(bg_mesh_type), bg_mesh_type);
		rf.write_attribute(bg_mesh_id, "x0", md.get_x0());
		rf.write_attribute(bg_mesh_id, "y0", md.get_y0());
		rf.write_attribute(bg_mesh_id, "elem_x_num", md.get_elem_x_num());
		rf.write_attribute(bg_mesh_id, "elem_y_num", md.get_elem_y_num());
		rf.write_attribute(bg_mesh_id, "hx", md.get_hx());
		rf.write_attribute(bg_mesh_id, "hy", md.get_hy());
		
		rf.close_group(bg_mesh_id);

		return 0;
	}


	int load_model_data_from_hdf5_file(Model_S2D_ME_s &md, ResultFile_hdf5 &rf)
	{
		// model data output
		hid_t md_out_id = rf.get_model_data_grp_id();
		
		hid_t bg_mesh_id = rf.open_group(md_out_id, "BackgroundMesh");

		double x0, y0, hx, hy;
		size_t elem_x_num, elem_y_num;
		rf.read_attribute(bg_mesh_id, "x0", x0);
		rf.read_attribute(bg_mesh_id, "y0", y0);
		rf.read_attribute(bg_mesh_id, "elem_x_num", elem_x_num);
		rf.read_attribute(bg_mesh_id, "elem_y_num", elem_y_num);
		rf.read_attribute(bg_mesh_id, "hx", hx);
		rf.read_attribute(bg_mesh_id, "hy", hy);
		md.init_mesh(
			x0, y0, 
			x0 + hx * double(elem_x_num),
			y0 + hy * double(elem_y_num),
			elem_x_num, elem_y_num
			);

		rf.close_group(bg_mesh_id);

		return 0;
	}


	int output_pcl_data_to_hdf5_file(Model_S2D_ME_s &md, ResultFile_hdf5 &rf, hid_t frame_id)
	{
		size_t pcl_num = md.get_pcl_num();
		Model_S2D_ME_s::Particle *pcls = md.get_pcls();
		ParticleData *pcl_data = new ParticleData[pcl_num];
		for (size_t p_id = 0; p_id < pcl_num; ++p_id)
		{
			Model_S2D_ME_s::Particle &pcl = pcls[p_id];
			ParticleData &pd = pcl_data[p_id];
			pd.id = pcl.id;
			pd.from_pcl(pcl);
		}
		hid_t pcl_dt_id = get_pcl_dt_id();
		int res = rf.write_dataset(
			frame_id,
			"ParticleData",
			pcl_num,
			pcl_data,
			pcl_dt_id
		);
		H5Tclose(pcl_dt_id);
		delete[] pcl_data;
		// particle num
		hid_t pcl_data_id = rf.open_dataset(frame_id, "ParticleData");
		rf.write_attribute(pcl_data_id, "pcl_num", pcl_num);
		rf.close_dataset(pcl_data_id);
		return res;
	}

	int load_pcl_data_from_hdf5_file(Model_S2D_ME_s &md, ResultFile_hdf5 &rf, hid_t frame_id)
	{
		// particle num
		size_t pcl_num;
		hid_t pcl_data_id = rf.open_dataset(frame_id, "ParticleData");
		rf.read_attribute(pcl_data_id, "pcl_num", pcl_num);
		rf.close_dataset(pcl_data_id);

		ParticleData *pcls_data = new ParticleData[pcl_num];
		hid_t pcl_dt_id = get_pcl_dt_id();
		int res = rf.read_dataset(
			frame_id,
			"ParticleData",
			pcl_num,
			pcls_data,
			pcl_dt_id
		);
		H5Tclose(pcl_dt_id);
		if (res) return res;
		md.alloc_pcls(pcl_num);
		Model_S2D_ME_s::Particle *pcls = md.get_pcls();
		for (size_t p_id = 0; p_id < pcl_num; ++p_id)
		{
			ParticleData &pcl_data = pcls_data[p_id];
			Model_S2D_ME_s::Particle &pcl = pcls[p_id];
			pcl.id = pcl_data.id;
			pcl_data.to_pcl(pcl);
		}
		delete[] pcls_data;
		return 0;
	}

};