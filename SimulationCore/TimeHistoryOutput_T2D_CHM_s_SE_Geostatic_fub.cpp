#include "SimulationCore_pcp.h"

#include "Model_T2D_CHM_s.h"
#include "Step_T2D_CHM_s_SE_Geostatic.h"

#include "TimeHistoryOutput_T2D_CHM_s_SE_Geostatic_fub.h"

int TimeHistoryOutput_T2D_CHM_s_SE_Geostatic_fub::init(void)
{
	switch (res_file->get_type())
	{
	case ResultFileType::PlainBin:
		break;
	case ResultFileType::XML:
		break;
	case ResultFileType::Hdf5:
	{
		ResultFile_hdf5 &rf = *static_cast<ResultFile_hdf5 *>(res_file);
		hid_t th_grp_id = rf.get_time_history_grp_id();
		th_id = rf.create_group(th_grp_id, name.c_str());
	}
	break;
	default:
		break;
	}
	return 0;
}

int TimeHistoryOutput_T2D_CHM_s_SE_Geostatic_fub::close(void)
{
	switch (res_file->get_type())
	{
	case ResultFileType::PlainBin:
		break;
	case ResultFileType::XML:
		break;
	case ResultFileType::Hdf5:
		{
		ResultFile_hdf5 &rf = *static_cast<ResultFile_hdf5 *>(res_file);
		rf.write_attribute(th_id, "output_num", output_id);
		if (th_id > 0)
		{
			rf.close_group(th_id);
			th_id = -1;
		}
		}
		break;
	default:
		break;
	}
	return 0;
}

int time_history_output_func_t2d_chm_s_SE_geostatic_fub_to_plain_bin_res_file(TimeHistoryOutput &_self)
{
	return 0;
}

int time_history_output_func_t2d_chm_s_SE_geostatic_fub_to_xml_res_file(TimeHistoryOutput &_self)
{
	TimeHistoryOutput_T2D_CHM_s_SE_Geostatic_fub &th
		= static_cast<TimeHistoryOutput_T2D_CHM_s_SE_Geostatic_fub &>(_self);
	Step_T2D_CHM_s_SE_Geostatic &step
		= static_cast<Step_T2D_CHM_s_SE_Geostatic &>(th.get_step());
	ResultFile_XML &rf = static_cast<ResultFile_XML &>(*th.res_file);
	std::fstream &file = rf.get_file();
	
	char str_buffer[512];
#define str_buffer_len (sizeof(str_buffer) / sizeof(str_buffer[0]))

	// time history
	const char *time_history_info = ""
		"<TimeHistory>\n"
		"    <substep_num> %zu </substep_num>\n"
		"    <total_substep_num> %zu </total_substep_num>\n"
		"    <current_time> %16.10e </current_time>\n"
		"    <total_time> %16.10e </total_time>\n";
	snprintf(str_buffer, str_buffer_len, time_history_info,
		step.get_substep_num(),  step.get_total_substep_num(),
		step.get_current_time(), step.get_total_time());
	file.write(str_buffer, strlen(str_buffer));
	
	// output nodal force
	Model_T2D_CHM_s &model = static_cast<Model_T2D_CHM_s &>(th.get_model());
	const char *node_info = ""
		"    <Node>\n"
		"        <node_num> %zu </node_num>\n";
	snprintf(str_buffer, str_buffer_len, node_info, model.node_num);
	file.write(str_buffer, strlen(str_buffer));
	// nodal force: fx_ext, fy_ext, fx_int, fy_int, fx_ub, fy_ub
	file << "        <nodal_force>\n"
			"        <!-- fx_ext, fy_ext, fx_int, fy_int, fx_ub, fy_ub -->\n";
	for (size_t n_id = 0; n_id < model.node_num; ++n_id)
	{
		Model_T2D_CHM_s::Node &n = model.nodes[n_id];
		file << "        ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", n.fx_ext_s);
		file << str_buffer << ", ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", n.fy_ext_s);
		file << str_buffer << ", ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", n.fx_int_s);
		file << str_buffer << ", ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", n.fy_int_s);
		file << str_buffer << ", ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", n.fx_ext_s - n.fx_int_s);
		file << str_buffer << ", ";
		snprintf(str_buffer, str_buffer_len, "%16.10e", n.fy_ext_s - n.fy_int_s);
		file << str_buffer << "\n";
	}
	file << "        </nodal_force>\n"
			"    </Node>\n";

	// ending
	file << "</TimeHistory>\n";

	return 0;
}

namespace
{
struct NodalForce
{
	unsigned long long id;
	bool need_cal;
	bool has_vsx_bc;
	bool has_vsy_bc;
	bool has_vfx_bc;
	bool has_vfy_bc;
	double fx_ext;
	double fy_ext;
	double fx_int;
	double fy_int;
	double fx_ub;
	double fy_ub;
};
}

int time_history_output_func_t2d_chm_s_SE_geostatic_fub_to_hdf5_res_file(TimeHistoryOutput &_self)
{
	TimeHistoryOutput_T2D_CHM_s_SE_Geostatic_fub &th
		= static_cast<TimeHistoryOutput_T2D_CHM_s_SE_Geostatic_fub &>(_self);
	Step_T2D_CHM_s_SE_Geostatic &step
		= static_cast<Step_T2D_CHM_s_SE_Geostatic &>(th.get_step());
	Model_T2D_CHM_s &md = static_cast<Model_T2D_CHM_s &>(step.get_model());
	ResultFile_hdf5 &rf = static_cast<ResultFile_hdf5 &>(*th.res_file);

	char frame_name[30];
	snprintf(frame_name, 30, "frame_%zu", th.output_id);
	hid_t frame_grp_id = rf.create_group(th.th_id, frame_name);
	// step infos
	rf.write_attribute(frame_grp_id, "current_time", step.get_current_time());
	rf.write_attribute(frame_grp_id, "total_time", step.get_total_time());
	rf.write_attribute(frame_grp_id, "substep_num", step.get_substep_num());
	rf.write_attribute(frame_grp_id, "total_substep_num", step.get_total_substep_num());

	// output nodal force
	NodalForce *nfs = new NodalForce[md.node_num];
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		NodalForce &nf = nfs[n_id];
		Model_T2D_CHM_s::Node &n = md.nodes[n_id];
		nf.id = n.id;
		nf.need_cal = n.has_mp;
		nf.has_vsx_bc = n.has_vsx_bc;
		nf.has_vsy_bc = n.has_vsy_bc;
		nf.has_vfx_bc = n.has_vfx_bc;
		nf.has_vfy_bc = n.has_vfy_bc;
		nf.fx_ext = n.fx_ext_s;
		nf.fy_ext = n.fy_ext_s;
		nf.fx_int = n.fx_int_s;
		nf.fy_int = n.fy_int_s;
		nf.fx_ub = n.fx_ext_s - n.fx_int_s;
		nf.fy_ub = n.fy_ext_s - n.fy_int_s;
	}
	hid_t nf_type_id = H5Tcreate(H5T_COMPOUND, sizeof(NodalForce));
	H5Tinsert(nf_type_id, "id", HOFFSET(NodalForce, id), H5T_NATIVE_ULLONG);
	H5Tinsert(nf_type_id, "need_cal", HOFFSET(NodalForce, need_cal), H5T_NATIVE_HBOOL);
	H5Tinsert(nf_type_id, "has_vsx_bc", HOFFSET(NodalForce, has_vsx_bc), H5T_NATIVE_HBOOL);
	H5Tinsert(nf_type_id, "has_vsy_bc", HOFFSET(NodalForce, has_vsy_bc), H5T_NATIVE_HBOOL);
	H5Tinsert(nf_type_id, "has_vfx_bc", HOFFSET(NodalForce, has_vfx_bc), H5T_NATIVE_HBOOL);
	H5Tinsert(nf_type_id, "has_vfy_bc", HOFFSET(NodalForce, has_vfy_bc), H5T_NATIVE_HBOOL);
	H5Tinsert(nf_type_id, "fx_ext", HOFFSET(NodalForce, fx_ext), H5T_NATIVE_DOUBLE);
	H5Tinsert(nf_type_id, "fy_ext", HOFFSET(NodalForce, fy_ext), H5T_NATIVE_DOUBLE);
	H5Tinsert(nf_type_id, "fx_int", HOFFSET(NodalForce, fx_int), H5T_NATIVE_DOUBLE);
	H5Tinsert(nf_type_id, "fy_int", HOFFSET(NodalForce, fy_int), H5T_NATIVE_DOUBLE);
	H5Tinsert(nf_type_id, "fx_ub", HOFFSET(NodalForce, fx_ub), H5T_NATIVE_DOUBLE);
	H5Tinsert(nf_type_id, "fy_ub", HOFFSET(NodalForce, fy_ub), H5T_NATIVE_DOUBLE);
	hid_t dataspace_id = H5Screate_simple(1, &md.node_num, nullptr);
	hid_t nf_id = H5Dcreate(frame_grp_id, "NodalForce", nf_type_id,
						  dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Dwrite(nf_id, nf_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, nfs);
	H5Dclose(nf_id);
	H5Sclose(dataspace_id);
	H5Tclose(nf_type_id);
	delete[] nfs;

	rf.close_group(frame_grp_id);

	++th.output_id;
	return 0;
}
