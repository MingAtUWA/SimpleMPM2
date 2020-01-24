#include "SimulationCore_pcp.h"

#include "ResultFile_hdf5.h"

ResultFile_hdf5::ResultFile_hdf5() : file_id(-1)
{
}

ResultFile_hdf5::~ResultFile_hdf5()
{
	close();
}

int ResultFile_hdf5::init(const char *file_name, bool over_write)
{
	if (over_write)
		file_id = H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	else
	{
		file_id = H5Fopen(file_name, H5F_ACC_RDWR, H5P_DEFAULT);
		if (file_id < 0)
			file_id = H5Fcreate(file_name, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
	}
	return file_id >= 0 ? 0 : -1;
}

void ResultFile_hdf5::close(void) { if (file_id >= 0) H5Fclose(file_id); }

hid_t ResultFile_hdf5::create_group(hid_t parent_id, const char *name)
{
	return H5Gcreate(parent_id, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
}

hid_t ResultFile_hdf5::open_group(hid_t parent_id, const char *name)
{
	hid_t grp_id = H5Gopen(parent_id, name, H5P_DEFAULT);
	if (grp_id < 0) // if group doesn't exist, then create new group
		grp_id = H5Gcreate(parent_id, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	return grp_id;
}

void ResultFile_hdf5::close_group(hid_t id) { if (id < 0) H5Gclose(id); }

hid_t ResultFile_hdf5::create_dataset(hid_t grp_id, const char *name, size_t num)
{
	hid_t dataspace_id = H5Screate_simple(1, &num, NULL);
	hid_t dataset_id = H5Dcreate(grp_id, name, H5T_NATIVE_DOUBLE, dataspace_id,
								 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Sclose(dataspace_id);
	return dataset_id;
}

hid_t ResultFile_hdf5::create_dataset(hid_t grp_id, const char *name, size_t row_num, size_t col_num)
{
	hsize_t dims[2] = { row_num, col_num };
	hid_t dataspace_id = H5Screate_simple(2, dims, nullptr);
	hid_t dataset_id = H5Dcreate(grp_id, name, H5T_NATIVE_DOUBLE, dataspace_id,
								 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Sclose(dataspace_id);
	return dataset_id;
}
void ResultFile_hdf5::close_dataset(hid_t id) { if (id < 0) H5Dclose(id); }

int ResultFile_hdf5::write_dataset(hid_t set_id, size_t num, double *data)
{
	hid_t dataspace_id = H5Screate_simple(1, &num, nullptr);
	herr_t res = H5Dwrite(set_id, H5T_NATIVE_DOUBLE, dataspace_id, dataspace_id, H5P_DEFAULT, data);
	H5Sclose(dataspace_id);
	return res < 0 ? -1 : 0;
}

int ResultFile_hdf5::write_dataset(hid_t set_id, size_t row_num, size_t col_num, double *data)
{
	hsize_t dims[2] = { row_num, col_num };
	hid_t dataspace_id = H5Screate_simple(2, dims, nullptr);
	herr_t res = H5Dwrite(set_id, H5T_NATIVE_DOUBLE, dataspace_id, dataspace_id, H5P_DEFAULT, data);
	H5Sclose(dataspace_id);
	return res < 0 ? -1 : 0;
}

int ResultFile_hdf5::write_dataset(hid_t set_id, size_t row_num, size_t col_num, unsigned long long *data)
{
	hsize_t dims[2] = { row_num, col_num };
	hid_t dataspace_id = H5Screate_simple(2, dims, nullptr);
	herr_t res = H5Dwrite(set_id, H5T_NATIVE_ULLONG, dataspace_id, dataspace_id, H5P_DEFAULT, data);
	H5Sclose(dataspace_id);
	return res < 0 ? -1 : 0;
}

// ========================= create model group =========================
hid_t ResultFile_hdf5::create_model(const char *name)
{
	return create_group(file_id, name);
}
hid_t ResultFile_hdf5::open_model(const char *name)
{
	return open_group(file_id, name);
}
void ResultFile_hdf5::close_model(hid_t id) { close_group(id); }

hid_t ResultFile_hdf5::create_model_output(hid_t md_id, const char *name)
{
	return create_group(md_id, name);
}
void ResultFile_hdf5::close_model_output(hid_t id) { close_group(id); }

hid_t ResultFile_hdf5::create_time_history(hid_t md_id, const char *name)
{
	return create_group(md_id, name);
}
void ResultFile_hdf5::close_time_history(hid_t id) { close_group(id); }


int ResultFile_hdf5::add_attribute(hid_t grp_id, const char *name, double value)
{
	hid_t attr_id = H5Acreate(grp_id, name, H5T_NATIVE_DOUBLE,
							  H5S_SCALAR, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &value);
	H5Aclose(attr_id);
	return 0;
}

int ResultFile_hdf5::add_attribute(hid_t grp_id, const char *name, size_t value)
{
	unsigned long long _value = (unsigned long long)value;
	hid_t attr_id = H5Acreate(grp_id, name, H5T_NATIVE_ULLONG,
							  H5S_SCALAR, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attr_id, H5T_NATIVE_ULLONG, &_value);
	H5Aclose(attr_id);
	return 0;
}

int ResultFile_hdf5::add_attribute(hid_t grp_id, const char *name, const char *str)
{
	size_t str_len = strlen(str);
	hid_t dataspace_id = H5Screate_simple(1, &str_len, nullptr);
	hid_t attr_id = H5Acreate(grp_id, name, H5T_NATIVE_CHAR,
							  dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
	H5Sclose(dataspace_id);
	H5Awrite(attr_id, H5T_NATIVE_CHAR, str);
	H5Aclose(attr_id);
	return 0;
}
