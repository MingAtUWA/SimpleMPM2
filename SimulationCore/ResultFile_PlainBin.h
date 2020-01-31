#ifndef __RESULT_FILE_PLAIN_BIN_H__
#define __RESULT_FILE_PLAIN_BIN_H__

#include "SimulationCore_pcp.h"

#include <fstream>
#include "ResultFile.h"

class ResultFile_PlainBin : public ResultFile
{
protected:
	std::fstream file;

public:
	ResultFile_PlainBin() : ResultFile(ResultFileType::PlainBin) {}
	~ResultFile_PlainBin() { finalize(); }
	int init(const char *file_name)
	{
		file.open(file_name, std::ios::binary | std::ios::out);
		if (file.is_open())
			return 0;
		return -1;
	}
	void finalize(void) { file.close(); }

	std::fstream &get_file(void) noexcept { return file; }
};

#endif