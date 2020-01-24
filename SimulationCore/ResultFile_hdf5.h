#ifndef __Result_File_hdf5_H__
#define __Result_File_hdf5_H__

#include "hdf5.h"

#include "ResultFile.h"

class ResultFile_hdf5 : public ResultFile
{
protected:
	hid_t file_id;
	
public:
	ResultFile_hdf5();
	~ResultFile_hdf5();

	// open hdf5 file
	int init(const char *file_name, bool over_write = false);
	inline hid_t get_file_id(void) noexcept { return file_id; }
	void close(void);

	// group
	hid_t create_group(hid_t parent_id, const char *name);
	hid_t open_group(hid_t parent_id, const char *name);
	void close_group(hid_t id);
	// dataset
	// one dimension
	hid_t create_dataset(hid_t grp_id, const char *name, size_t num);
	// two dimension
	hid_t create_dataset(hid_t grp_id, const char *name, size_t row_num, size_t col_num);
	void close_dataset(hid_t id);
	int write_dataset(hid_t set_id, size_t num, double *data);
	int write_dataset(hid_t set_id, size_t row_num, size_t col_num, double *data);
	int write_dataset(hid_t set_id, size_t row_num, size_t col_num, unsigned long long *data);

	// create model
	hid_t create_model(const char *name);
	hid_t open_model(const char *name);
	void close_model(hid_t id);
	// create model output
	hid_t create_model_output(hid_t md_id, const char *name);
	void close_model_output(hid_t id);
	// create time history
	hid_t create_time_history(hid_t md_id, const char *name);
	void close_time_history(hid_t id);

	// add attribute
	int add_attribute(hid_t grp_id, const char *name, double value);
	int add_attribute(hid_t grp_id, const char *name, size_t value);
	int add_attribute(hid_t grp_id, const char *name, const char *str);
};

#endif